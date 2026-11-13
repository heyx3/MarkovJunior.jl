# Tests for `@within_regions`, which runs Ops within each connected region of certain colors.
# Uses `CountingBias` from scale.jl, so it must be included after that file.

(() -> begin

c(ch::Char) = MJ.CELL_CODE_BY_CHAR[ch]
cells_of(grid, ch::Char) = Set(Tuple.(findall(==(c(ch)), grid)))

"Makes a 2D grid from rows of color chars, where each string is one row (so X is the column)"
grid_from_rows(rows::String...) = UInt8[ c(rows[y][x]) for x in 1:length(rows[1]), y in 1:length(rows) ]
"Makes a 1D grid from a string of color chars"
grid_1d(s::String) = UInt8[ c(ch) for ch in s ]
"Turns a grid back into rows of color chars"
rows_of(grid::Matrix{UInt8}) = [ String([ MJ.CELL_TYPES[grid[x, y] + 1].char for x in 1:size(grid, 1) ])
                                   for y in 1:size(grid, 2) ]
string_of(grid::Vector{UInt8}) = String([ MJ.CELL_TYPES[u + 1].char for u in grid ])

function check_roundtrip(mj::MarkovAlgorithm)
    test_parsed_mj = markov_algo_parse(markov_algo_to_string(mj))
    @bp_check(mj == test_parsed_mj,
              "Parsing failed!\n\nSource:\n", mj,
              "\nDest:\n", test_parsed_mj,
              "\nSource as string: ", markov_algo_to_string(mj))
end

"Runs the given algorithm to completion on the given grid, skipping all tagged events"
function run_on(mj::MarkovAlgorithm, initial::Array{UInt8}; seed = 0xaabbcc)
    check_roundtrip(mj)
    state_ch = markov_algo_run(mj, initial,
                               MarkovTickSettings(skip_tagged_events = MJ.SkipTaggedEvents.all_skippable),
                               seeds=seed,
                               allocator=MJ.MarkovAllocatorHeapReused())
    return markov_algo_complete(copy, state_ch)
end

"Checks that parsing fails (but not the error message, which may change)"
function check_parse_fails(src::String)
    try
        markov_algo_parse(src)
    catch e
        return nothing
    end
    error("Should have failed to parse:\n", src)
end


# Every region of each rule gets its own Ops; other cells are untouched
#    including the external cells near a non-rectangular region.
grid = run_on(@markovjunior(begin
    @within_regions border(0) begin
        R => @sequence begin
            @rewrite R => G
            @rewrite b => w
        end
        MO => @rewrite M => Y
    end
end), grid_from_rows("RRbbMO",
                     "RbbbOO",
                     "bbRbbb",
                     "MbRRbN"))
@bp_check(rows_of(grid) == [ "GGbbYO",
                             "GbbbOO",
                             "bbGbbb",
                             "YbGGbN" ], rows_of(grid))

# Regions don't mix rules.
grid = run_on(@markovjunior(begin
    @within_regions border(0) begin
        R => @rewrite R => G
        B => @rewrite B => Y
    end
end), grid_from_rows("RBBR"))
@bp_check(rows_of(grid) == [ "GYYG" ], rows_of(grid))

# `uv` boxes and thresholds are relative to each region.
grid = run_on(@markovjunior(begin
    @within_regions border(0) begin
        R => @sequence begin
            @fill 'Y' uv(center=0.5, size=texel)
            @rewrite (size_prod()/3) R => G
        end
    end
end), grid_from_rows("RRRbRRRRRR",
                     "RRRbRRRRRR",
                     "RRRbRRRRRR"))
@bp_check(cells_of(grid, 'Y') == Set([ (2, 2), (8, 2) ]), rows_of(grid))
@bp_check(count(==(c('G')), grid[1:3, :]) == 3, rows_of(grid))
@bp_check(count(==(c('G')), grid[5:10, :]) == 6, rows_of(grid))

# External cells within the region's bounding box are visible but have their edits reverted.
grid = run_on(@markovjunior(begin
    @within_regions border(0) begin
        R => @sequence begin
            @rewrite bR => bM
            @rewrite b => w
        end
    end
end), grid_from_rows("RRbN",
                     "RbbN",
                     "RRRN"))
@bp_check(rows_of(grid) == [ "RMbN",
                             "MbbN",
                             "RMMN" ], rows_of(grid))

# 3D regions.
grid_3d = fill(c('b'), 3, 3, 3)
grid_3d[1, 1, 1:3] .= c('R')
grid_3d[3, 3, 3] = c('R')
grid = run_on(@markovjunior(begin
    @within_regions border(0) begin
        R => @fill 'G' pixel(min=1, size=1)
    end
end), grid_3d)
@bp_check(cells_of(grid, 'G') == Set([ (1, 1, 1), (3, 3, 3) ]), grid)

# Nesting with `@within_box`, both ways.
grid = run_on(@markovjunior(begin
    @within_box pixel(min=1, size=(3, 2)) begin
        @within_regions border(0) begin
            R => @sequence begin
                @within_box pixel(min=1, size=1) begin
                    @rewrite R => G
                end
            end
        end
    end
end), grid_from_rows("RRbRR",
                     "bbRRR",
                     "RRRRR"))
@bp_check(rows_of(grid) == [ "GRbRR",
                             "bbGRR",
                             "RRRRR" ], rows_of(grid))

# Nothing to do.
grid = run_on(@markovjunior(begin
    @sequence repeat begin
        @within_regions border(0) begin
            R => @rewrite b => w
        end
    end
end), grid_from_rows("bbb"))
@bp_check(rows_of(grid) == [ "bbb" ], rows_of(grid))


#####################
# Threshold & biases

# A threshold limits how many regions run.
for seed in 1:5
    grid = run_on(@markovjunior(begin
        @within_regions border(0) 2 begin
            R => @rewrite R => G
        end
    end), grid_1d("RbRRbRbRRR"); seed=seed)
    s = string_of(grid)
    n_green_regions = count(r -> all(==('G'), r), split(s, 'b'))
    @bp_check(n_green_regions == 2, seed, ": ", s)
    @bp_check(count(r -> all(==('R'), r), split(s, 'b')) == 2, seed, ": ", s)
end

# Biases pick which regions run.
for seed in 1:5
    grid = run_on(@markovjunior(begin
        @within_regions border(0) 1 begin
            R => @rewrite R => G
        end field(Y)
    end), grid_1d("RRbbbRRbYbRRbbbbbRR"); seed=seed)
    @bp_check(string_of(grid) == "RRbbbGGbYbRRbbbbbRR" || string_of(grid) == "RRbbbRRbYbGGbbbbbRR",
              seed, ": ", string_of(grid))
end

# Live biases see each region's changes before the next one is picked.
for (live, expected) in ((false, "YbbbbRRbbbGGbbbYYYY"), (true, "YbbbbGGbbbRRbbbYYYY"))
    for seed in 1:5
        grid = run_on(markov_algo_parse("""@markovjunior 1 begin
            @within_regions border(0) 2 begin
                O => @fill 'Y'
                R => @rewrite R => G
            end field(-Y$(live ? ", live" : ""))
        end"""), grid_1d("YbbbbRRbbbRRbbbOOOO"); seed=seed)
        @bp_check(string_of(grid) == expected, "live=", live, " seed=", seed, ": ", string_of(grid))
    end
end

# Biases that forbid every region stop the Op.
grid = run_on(@markovjunior(begin
    @within_regions border(0) begin
        R => @rewrite R => G
    end field(Y->bR)
end), grid_1d("RRwYbbRR"))
@bp_check(string_of(grid) == "RRwYbbGG", string_of(grid))

# `inherit_biases` brings the Op's biases into each region (seeing only that region).
for seed in 1:5
    grid = run_on(@markovjunior(begin
        @within_regions border(0) inherit_biases begin
            RY => @rewrite 1 R => G
        end field(Y)
    end), grid_1d("RRRRYbYRRRR"); seed=seed)
    @bp_check(string_of(grid) == "RRRGYbYGRRR", seed, ": ", string_of(grid))
end

# Bias states shouldn't leak, whether or not they're inherited, even when canceled partway.
for inherit in (false, true)
    mj = markov_algo_parse("""@markovjunior 1 begin
        @sequence begin
            @within_regions border(0) $(inherit ? "inherit_biases" : "") begin
                R => @rewrite R => G counting()
            end counting()
        end counting()
    end""")
    COUNTING_BIAS_N_INITIALIZED[] = 0
    run_on(mj, grid_1d("RRbRRbRR"))
    @bp_check(isempty(COUNTING_BIAS_LIVE_STATES), "inherit=", inherit, ": ", COUNTING_BIAS_LIVE_STATES)
    # Without inheriting: 2 outer + 3 inner.
    # With inheriting: also a separate state of the 2 outer biases for each of the 3 regions.
    @bp_check(COUNTING_BIAS_N_INITIALIZED[] == (inherit ? (5 + 6) : 5),
              "inherit=", inherit, ": ", COUNTING_BIAS_N_INITIALIZED[])

    for cancel_at_tick in 1:4
        ticking = MarkovTickSettings()
        ch = markov_algo_run(mj, grid_1d("RRbRRbRR"), ticking, seeds=0xaabbcc,
                             allocator=MJ.MarkovAllocatorHeapReused())
        n_ticks = 0
        while true
            event = markov_algo_next(ch)
            if event isa Int
                n_ticks += 1
                (n_ticks == cancel_at_tick) && (ticking.cancel_algo = true)
            elseif event in (MJ.TAG_ALGO_CANCELED, MJ.TAG_ALGO_COMPLETED)
                break
            end
        end
        markov_algo_cleanup(ch)
        @bp_check(isempty(COUNTING_BIAS_LIVE_STATES),
                  "inherit=", inherit, ", cancel_at_tick=", cancel_at_tick, ": ", COUNTING_BIAS_LIVE_STATES)
    end
end


################
# Shared colors

# Rules can share colors: regions start from a rule's own colors, earlier rules win contested cells,
#    and leftover groups of shared colors go to a random rule which has them.
let str(grid) = String([ MJ.CELL_TYPES[v + 1].char for v in grid ]),
    initial = [ MJ.CELL_CODE_BY_CHAR[ch] for ch in "RYbYBbYbRYB" ],
    mj = markov_algo_parse("""@markovjunior 1 begin
        @within_regions border(0) begin
            RY => @rewrite Y => G
            BY => @rewrite Y => N
        end
    end""")
    lone_yellows = Set{Char}()
    for seed in 1:20
        s = str(run_on(mj, initial; seed=seed))
        @bp_check(occursin(r"^RGbNBb[GN]bRGB$", s), s)
        push!(lone_yellows, s[7])
    end
    @bp_check(lone_yellows == Set([ 'G', 'N' ]), lone_yellows)
end

# Inside a region, `area` counts the region's own cells rather than its bounding box.
let mj = markov_algo_parse("@markovjunior 2 begin
 @within_regions border(0) R => @rewrite (area - 1) R => G
end"),
    initial = fill(MJ.CELL_CODE_BY_CHAR['b'], 3, 3)
    initial[1, 1] = initial[2, 1] = initial[1, 2] = MJ.CELL_CODE_BY_CHAR['R'] # An L shape, in a 2x2 box
    grid = run_on(mj, initial)
    @bp_check(count(==(MJ.CELL_CODE_BY_CHAR['G']), grid) == 2, grid)
end

# By default the region's Ops can see one cell of real colors around it, but changes to those cells are undone.
let str(grid) = String([ MJ.CELL_TYPES[v + 1].char for v in grid ]),
    run_1D(src, row) = str(run_on(markov_algo_parse("@markovjunior 1 begin
 $src
end"),
                                   [ MJ.CELL_CODE_BY_CHAR[ch] for ch in row ]))
    @bp_check((s = run_1D("@within_regions R => @rewrite LR => LG", "LRRL")) == "LGGL", s) # Mirrored too
    @bp_check((s = run_1D("@within_regions border(0) R => @rewrite LR => LG", "LRRL")) == "LRRL", s)
    @bp_check((s = run_1D("@within_regions R => @rewrite RL => RY", "LRRL")) == "LRRL", s)
    # Farther cells are outside the view.
    @bp_check((s = run_1D("@within_regions R => @rewrite LR => LG", "LbRR")) == "LbRR", s)
    @bp_check((s = run_1D("@within_regions border(2) R => @rewrite LbR => LbG", "LbRR")) == "LbGR", s)
    # Other cells with the rule's colors are masked out with the 'C' color, even within the border.
    @bp_check((s = run_1D("@within_regions border(2) RG(w) => @rewrite Rbw => Gbw", "RbR")) == "GbG", s)
    @bp_check((s = run_1D("@within_regions border(2) RG(w) => @rewrite RbR => GbR", "RbR")) == "RbR", s)
end

# A concave region never sees another region of its colors inside its box.
let initial = [ MJ.CELL_CODE_BY_CHAR[ch] for ch in [ 'R' 'R' 'R' 'R' 'R'
                                                      'R' 'b' 'b' 'b' 'R'
                                                      'R' 'b' 'R' 'b' 'R'
                                                      'R' 'b' 'b' 'b' 'R' ] ],
    mj = markov_algo_parse("@markovjunior 2 begin\n @within_regions R => @rewrite 1 R => G\nend")
    for seed in 1:20
        grid = run_on(mj, initial; seed = seed)
        @bp_check(grid[3, 3] == MJ.CELL_CODE_BY_CHAR['G'], seed, ": ", grid)
        @bp_check(count(==(MJ.CELL_CODE_BY_CHAR['G']), grid) == 2, seed, ": ", grid)
    end
end

# Regions can be limited to spreading along certain axes.
let initial = fill(MJ.CELL_CODE_BY_CHAR['R'], 3, 2)
    for (setting, n_regions) in (("", 1), ("along(x)", 2), ("along(y)", 3), ("along(x, y)", 1))
        mj = markov_algo_parse("@markovjunior 2 begin
 @within_regions border(0) $setting R => @rewrite 1 R => G
end")
        grid = run_on(mj, initial)
        @bp_check(count(==(MJ.CELL_CODE_BY_CHAR['G']), grid) == n_regions, setting, ": ", grid)
    end
end

# A single rule can be written without a block.
let mj = markov_algo_parse("@markovjunior 1 begin\n @within_regions border(0) 1 R => @rewrite R => G\nend"),
    grid = run_on(mj, [ MJ.CELL_CODE_BY_CHAR[ch] for ch in "RRbRR" ])
    @bp_check(count(==(MJ.CELL_CODE_BY_CHAR['G']), grid) == 2, grid)
end


##############
# Parse errors

check_parse_fails("@markovjunior begin\n @within_regions border(0) begin\n R => @upscale (2)\n end\nend")
check_parse_fails("@markovjunior begin\n @within_regions border(0) begin\n R\n end\nend")
check_parse_fails("@markovjunior begin\n @within_regions border(0) begin\n {RG} => @rewrite R => G\n end\nend")
check_parse_fails("@markovjunior begin\n @within_regions border(0) begin\n RX => @rewrite R => G\n end\nend")
check_parse_fails("@markovjunior begin\n @within_regions border(0) begin\n end\nend")
check_parse_fails("@markovjunior begin\n @within_regions\nend")
check_parse_fails("@markovjunior begin\n @within_regions border(0) sometimes begin\n R => @rewrite R => G\n end\nend")
check_parse_fails("@markovjunior begin\n @within_regions border(0) 2 3 begin\n R => @rewrite R => G\n end\nend")
check_parse_fails("@markovjunior begin\n @within_regions border(0) begin\n R(X) => @rewrite R => G\n end\nend")

end)()
