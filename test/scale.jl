# Tests for the scale Ops (`@upscale` and `@downscale`), including how biases survive the grid reallocation.

"A test bias which tracks its live states, to catch leaks and double-cleanups"
struct CountingBias <: MJ.AbstractMarkovBias end
const COUNTING_BIAS_LIVE_STATES = Set{Int}()
const COUNTING_BIAS_N_INITIALIZED = Ref(0)
function MJ.markov_bias_initialize(::CountingBias, algo, algo_state)
    COUNTING_BIAS_N_INITIALIZED[] += 1
    push!(COUNTING_BIAS_LIVE_STATES, COUNTING_BIAS_N_INITIALIZED[])
    return COUNTING_BIAS_N_INITIALIZED[]
end
function MJ.markov_bias_cleanup(::CountingBias, id::Int, algo, algo_state)
    @bp_check(id in COUNTING_BIAS_LIVE_STATES, "Bias state ", id, " was cleaned up twice")
    delete!(COUNTING_BIAS_LIVE_STATES, id)
    return nothing
end
MJ.markov_bias_calculate(::CountingBias, ::Int, algo, algo_state,
                         ::Union{MJ.CellLine{N}, MJ.CellRegion{N}}) where {N} = 1.0f0
MJ.dsl_format(::CountingBias) = "counting()"
MJ.parse_markovjunior_bias(::Val{:counting}, inputs::MJ.MacroParserInputs, loc, args) = CountingBias()


(() -> begin

c(ch::Char) = MJ.CELL_CODE_BY_CHAR[ch]

"Makes a 2D grid from rows of color chars, where each string is one row (so X is the column)"
grid_from_rows(rows::String...) = UInt8[ c(rows[y][x]) for x in 1:length(rows[1]), y in 1:length(rows) ]

"Checks a grid against rows of color chars, where `?` may be any color in `random_options`"
function check_rows(grid::Matrix{UInt8}, random_options::String, rows::String...)
    @bp_check(size(grid) == (length(rows[1]), length(rows)), "Size: ", size(grid))
    for x in 1:size(grid, 1), y in 1:size(grid, 2)
        expected = rows[y][x]
        actual = MJ.CELL_TYPES[grid[x, y] + 1].char
        @bp_check((expected == '?') ? (actual in random_options) : (actual == expected),
                  "At (", x, ", ", y, ") expected '", expected, "' but got '", actual, "'")
    end
end

function check_roundtrip(mj::MarkovAlgorithm)
    test_parsed_mj = markov_algo_parse(markov_algo_to_string(mj))
    @bp_check(mj == test_parsed_mj,
              "Parsing failed!\n\nSource:\n", mj,
              "\nDest:\n", test_parsed_mj,
              "\nSource as string: ", markov_algo_to_string(mj))
end

"Runs the given algorithm to completion on the given grid, skipping all tagged events"
function run_on(mj::MarkovAlgorithm, initial::Array{UInt8})
    check_roundtrip(mj)
    # Use a normal heap allocator so we can steal the grid after the algo is done.
    state_ch = markov_algo_run(mj, initial,
                               MarkovTickSettings(skip_tagged_events = MJ.SkipTaggedEvents.all_skippable),
                               seeds=0xaabbcc,
                               allocator=MJ.MarkovAllocatorHeap())
    return markov_algo_complete(copy, state_ch)
end

function check_parse_fails(src::String, expected_msg::String)
    try
        markov_algo_parse(src)
    catch e
        msg = sprint(showerror, e)
        @bp_check(occursin(expected_msg, msg), "Expected error containing '", expected_msg, "', got:\n", msg)
        return nothing
    end
    error("Should have failed to parse:\n", src)
end


# Plain stretching.
grid_3d = reshape(UInt8.(0:7), 2, 2, 2)
for (factors_str, inner) in (("(2, 3)", (2, 3, 1)), ("(2...)", (2, 2, 2)), ("(2)", (2, 1, 1)),
                             ("(3, 2...)", (3, 2, 2)), ("(1, 1, 2)", (1, 1, 2)))
    mj = markov_algo_parse("@markovjunior begin\n @upscale $factors_str\nend")
    grid = run_on(mj, grid_3d)
    @bp_check(grid == repeat(grid_3d, inner=inner), factors_str, " => ", grid)
end

# Rules: sets, wildcards, extrusion of size-1 axes, and the first matching rule wins.
grid = run_on(@markovjunior(begin
    @upscale (3, 2) begin
        R => [
            w b g
            M M M
        ]
        {GB} => [ _ {Yg} _ ]
        G => [ T T T ]
    end
end), grid_from_rows("RB", "Gw"))
check_rows(grid, "Yg",
    "wbgB?B",
    "MMMB?B",
    "G?Gwww",
    "G?Gwww"
)

# A single column.
grid = run_on(@markovjunior(begin
    @upscale (1, 2) Y => [ R ; G ]
end), grid_from_rows("Y"))
check_rows(grid, "", "R", "G")

# A 2D block extruded along Z.
grid = run_on(@markovjunior(begin
    @upscale (2, 2, 3) R => [
        w b
        g M
    ]
end), fill(c('R'), 1, 1, 2))
@bp_check(size(grid) == (2, 2, 6), size(grid))
for z in 1:6
    check_rows(grid[:, :, z], "", "wb", "gM")
end

# A YZ block extruded along X.
grid = run_on(@markovjunior(begin
    @upscale (1, 2, 2) R => [ w;b ;;; g;M ]
end), fill(c('R'), 2, 1, 1))
@bp_check(size(grid) == (2, 2, 2), size(grid))
check_rows(grid[:, :, 1], "", "ww", "bb")
check_rows(grid[:, :, 2], "", "gg", "MM")

# An upscale that does nothing should fail, so that 'repeat' can finish.
grid = run_on(@markovjunior(begin
    @sequence repeat begin
        @upscale (1, 1)
    end
end), grid_from_rows("RB"))
check_rows(grid, "", "RB")

# Parse errors.
check_parse_fails("@markovjunior begin\n @upscale (3, 2) R => [ w b ; g M ]\nend", "along the x axis")
check_parse_fails("@markovjunior begin\n @upscale (2) R => [ w b ; g M ]\nend", "along the y axis")
check_parse_fails("@markovjunior begin\n @upscale (0, 2)\nend", "positive integer factors")
check_parse_fails("@markovjunior begin\n @upscale (2.5)\nend", "positive integer factors")
check_parse_fails("@markovjunior begin\n @upscale ()\nend", "positive integer factors")
check_parse_fails("@markovjunior begin\n @upscale (2) R => G\nend", "Expected a block of cells")
check_parse_fails("@markovjunior begin\n @upscale (2) _ => [ w b ]\nend", "Unsupported color value '_'")
check_parse_fails("@markovjunior begin\n @upscale (2) R => [ [wb] b ]\nend", "List syntax")
check_parse_fails("@markovjunior begin\n @rewrite [ {RG} b ] => [ b b ]\nend", "try List syntax")

############
# Downscale

# The 'first' fallback, including partial blocks at the far edges.
grid_4x4 = reshape(UInt8.(0:15), 4, 4)
grid = run_on(@markovjunior(begin
    @downscale (2, 2) first
end), grid_4x4)
@bp_check(grid == grid_4x4[1:2:end, 1:2:end], grid)
grid_5x3 = reshape(UInt8.(0:14), 5, 3)
grid = run_on(@markovjunior(begin
    @downscale (2...) first
end), grid_5x3)
@bp_check(grid == grid_5x3[1:2:end, 1:2:end], grid)

# The default 'random' fallback picks from within each block.
grid = run_on(@markovjunior(begin
    @downscale (2, 2)
end), grid_4x4)
@bp_check(size(grid) == (2, 2), size(grid))
for x in 1:2, y in 1:2
    block = grid_4x4[(x*2 - 1):(x*2), (y*2 - 1):(y*2)]
    @bp_check(grid[x, y] in block, "At (", x, ", ", y, "): ", grid[x, y], " not in ", block)
end

# The 'common' and 'rare' fallbacks, with and without ties.
grid = run_on(@markovjunior(begin
    @downscale (2, 2) common
end), grid_from_rows("RRGBRG", "RGBBGR"))
check_rows(grid, "RG", "RB?")
grid = run_on(@markovjunior(begin
    @downscale (2, 2) rare
end), grid_from_rows("RRGBRG", "RGBBBB"))
check_rows(grid, "RG", "GG?")

# Rules: sets, wildcards, extrusion of size-1 axes, and the first matching rule wins.
grid = run_on(@markovjunior(begin
    @downscale (2, 2) begin
        R <= [
            w b
            b w
        ]
        {GB} <= [ _ M ; M _ ]
        T <= [ M M ; M M ]
        Y <= [ g ; g ]
        first
    end
end), grid_from_rows("wbMMgggw", "bwMMgggw"))
check_rows(grid, "GB", "R?Yg")

# Partial blocks never match block rules.
grid = run_on(@markovjunior(begin
    @downscale (2, 2) begin
        M <= [ _ _ ; _ _ ]
        first
    end
end), grid_from_rows("RGB", "wbg", "TYO"))
check_rows(grid, "", "MB", "TO")

# 'has' rules, interleaved with block rules, and the 'always' fallback.
# Partial blocks can match 'has' rules.
grid = run_on(@markovjunior(begin
    @downscale (2, 2) begin
        w <= has(w)
        N <= has(N)
        R <= [ G G ; G G ]
        {Gg} <= has(GT)
        always(gb)
    end
end), grid_from_rows("wNNGGGGBBB",
                     "GGGGGGBBbB",
                     "GwBBGGbNBB"))
check_rows(grid, "Ggb", "wNR??", "w??N?")
@bp_check(grid[4, 1] in (c('G'), c('g')) && grid[3, 2] in (c('G'), c('g')), grid)
@bp_check(grid[5, 1] in (c('g'), c('b')) && grid[2, 2] in (c('g'), c('b')) && grid[5, 2] in (c('g'), c('b')), grid)
grid = run_on(@markovjunior(begin
    @downscale (2, 2) always(g)
end), grid_4x4)
@bp_check(grid == fill(c('g'), 2, 2), grid)

# A 2D block extruded along Z.
grid_3d = cat(grid_from_rows("wb", "gM"), grid_from_rows("wb", "gM"), dims=3)
downscale_3d = @markovjunior begin
    @downscale (2, 2, 2) begin
        R <= [
            w b
            g M
        ]
        first
    end
end
@bp_check(run_on(downscale_3d, grid_3d) == fill(c('R'), 1, 1, 1))
grid_3d[2, 2, 2] = c('Y')
@bp_check(run_on(downscale_3d, grid_3d) == fill(c('w'), 1, 1, 1))

# A downscale that does nothing should fail, so that 'repeat' can finish.
grid = run_on(@markovjunior(begin
    @sequence repeat begin
        @downscale (1, 1)
    end
end), grid_from_rows("RB"))
check_rows(grid, "", "RB")

# Parse errors.
check_parse_fails("@markovjunior begin\n @downscale (2) _ <= [ w b ]\nend", "Unsupported color value '_'")
check_parse_fails("@markovjunior begin\n @downscale (2) R <= [ w b g ]\nend", "along the x axis")
check_parse_fails("@markovjunior begin\n @downscale (2) sometimes\nend", "Unknown fallback")
check_parse_fails("@markovjunior begin\n @downscale (2) R <= [ [wb] b ]\nend", "List syntax")
check_parse_fails("@markovjunior begin\n @downscale (2) R => [ w b ]\nend", "Expected a rule like `R <=")
check_parse_fails("@markovjunior begin\n @downscale (2) begin\n first\n R <= [ w b ]\n end\nend", "must be the last line")
check_parse_fails("@markovjunior begin\n @downscale (2) begin\n first\n common\n end\nend", "must be the last line")
check_parse_fails("@markovjunior begin\n @downscale (2) first R <= [ w b ]\nend", "got 3 arguments")
check_parse_fails("@markovjunior begin\n @downscale (2) R <= has({wb})\nend", "without braces")
check_parse_fails("@markovjunior begin\n @downscale (2) R <= has(w, b)\nend", "with one argument")
check_parse_fails("@markovjunior begin\n @downscale (2) always(_)\nend", "Unsupported color value '_'")


########
# Biases

# Biases should be re-created at the new grid size, including for nested sequences.
biased_mj = @markovjunior begin
    @sequence begin
        @sequence begin
            @upscale (2, 2)
        end field(R)
        @upscale (2, 1)
        @rewrite b => w
        @downscale (2, 1) first
    end counting()
end
COUNTING_BIAS_N_INITIALIZED[] = 0
grid = run_on(biased_mj, grid_from_rows("Rb", "bb"))
@bp_check(size(grid) == (4, 4), size(grid))
@bp_check(count(==(c('R')), grid) == 4, grid)
@bp_check(count(==(c('w')), grid) == 12, grid)
@bp_check(COUNTING_BIAS_N_INITIALIZED[] == 4, COUNTING_BIAS_N_INITIALIZED[])
@bp_check(isempty(COUNTING_BIAS_LIVE_STATES), COUNTING_BIAS_LIVE_STATES)

# Canceling right after any reallocation shouldn't leak or double-free bias states.
for cancel_at_grid in (2, 3, 4)
    ticking = MarkovTickSettings()
    ch = markov_algo_run(biased_mj, grid_from_rows("Rb", "bb"), ticking,
                         seeds=0xaabbcc, allocator=MJ.MarkovAllocatorHeap())
    n_grids = 0
    while true
        event = markov_algo_next(ch)
        if event isa MJ.CellGrid
            n_grids += 1
            if n_grids == cancel_at_grid
                ticking.cancel_algo = true
            end
        elseif event === MJ.TAG_ALGO_CANCELED
            break
        elseif event === MJ.TAG_ALGO_COMPLETED
            error("Algorithm should have been canceled at grid ", cancel_at_grid)
        end
    end
    markov_algo_cleanup(ch)
    @bp_check(isempty(COUNTING_BIAS_LIVE_STATES), "cancel_at_grid=", cancel_at_grid, ": ", COUNTING_BIAS_LIVE_STATES)
end

end)()
