# Tests for `@within_box`, which runs Ops on a view of part of the grid.
# Uses `CountingBias` from scale.jl, so it must be included after that file.

"A custom Op which doubles the grid size, but claims not to in `markov_op_reallocates_grid()`"
struct TestSneakyResizeOp <: MJ.AbstractMarkovOp end
MJ.markov_op_reallocates_grid(::TestSneakyResizeOp) = false
function MJ.markov_algo_run(::TestSneakyResizeOp, algo, algo_state,
                            inherited_biases, inherited_bias_states, grid = algo_state.grid)
    (grid, inherited_bias_states) = MJ.markov_opalgo_new_grid(
        algo, algo_state, size(grid) .* 2,
        inherited_biases, inherited_bias_states
    ) do old_grid, new_grid
        fill!(new_grid, zero(UInt8))
    end
    return true
end
MJ.dsl_format(::TestSneakyResizeOp) = "@test_sneaky_resize"
MJ.parse_markovjunior_op(::Val{Symbol("@test_sneaky_resize")}, inputs::MJ.MacroParserInputs,
                         loc, expr_args, expr_line) = TestSneakyResizeOp()

"A custom Op which doesn't implement `markov_op_reallocates_grid()`"
struct TestUnqueriedOp <: MJ.AbstractMarkovOp end
MJ.markov_algo_run(::TestUnqueriedOp, algo, algo_state, inherited_biases, inherited_bias_states) = nothing
MJ.dsl_format(::TestUnqueriedOp) = "@test_unqueried"
MJ.parse_markovjunior_op(::Val{Symbol("@test_unqueried")}, inputs::MJ.MacroParserInputs,
                         loc, expr_args, expr_line) = TestUnqueriedOp()

"A custom Op which records the current grid view's ID and depth into `TEST_VIEW_RECORDS`"
struct TestRecordViewOp <: MJ.AbstractMarkovOp end
const TEST_VIEW_RECORDS = Tuple{Int, Int}[ ]
MJ.markov_op_reallocates_grid(::TestRecordViewOp) = false
function MJ.markov_algo_run(::TestRecordViewOp, algo, algo_state, inherited_biases, inherited_bias_states)
    push!(TEST_VIEW_RECORDS, (algo_state.grid_view_id, algo_state.scope_depth))
    return nothing
end
MJ.dsl_format(::TestRecordViewOp) = "@test_record_view"
MJ.parse_markovjunior_op(::Val{Symbol("@test_record_view")}, inputs::MJ.MacroParserInputs,
                         loc, expr_args, expr_line) = TestRecordViewOp()

(() -> begin

c(ch::Char) = MJ.CELL_CODE_BY_CHAR[ch]
cells_of(grid, ch::Char) = Set(Tuple.(findall(==(c(ch)), grid)))
neighbors(x, y) = Set([ (x-1, y), (x+1, y), (x, y-1), (x, y+1) ])

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
    # The reusing allocator would choke if the grid was still a view at the end.
    state_ch = markov_algo_run(mj, initial,
                               MarkovTickSettings(skip_tagged_events = MJ.SkipTaggedEvents.all_skippable),
                               seeds=seed,
                               allocator=MJ.MarkovAllocatorHeapReused())
    return markov_algo_complete(copy, state_ch)
end

function check_parse_fails(src::String)
    try
        markov_algo_parse(src)
    catch e
        return nothing
    end
    error("Should have failed to parse:\n", src)
end


# Each grid view gets a unique ID, and leaving it restores the outer one's.
empty!(TEST_VIEW_RECORDS)
run_on(@markovjunior(begin
    @test_record_view
    @within_box pixel(min=1, size=2) begin
        @test_record_view
        @within_box pixel(min=1, size=1) begin
            @test_record_view
        end
        @test_record_view
    end
    @within_box pixel(min=1, size=2) begin # The same area as before, but a different view
        @test_record_view
    end
    @test_record_view
end), fill(c('b'), 4, 4))
@bp_check(TEST_VIEW_RECORDS == [ (0, 0), (1, 1), (2, 2), (1, 1), (3, 1), (0, 0) ], TEST_VIEW_RECORDS)

# Ops only touch the box.
grid = run_on(@markovjunior(begin
    @within_box pixel(min=(2, 3), size=(3, 2)) begin
        @rewrite b => w
    end
end), fill(c('b'), 6, 6))
@bp_check(cells_of(grid, 'w') == Set((x, y) for x in 2:4, y in 3:4), grid)

# Thresholds and `uv` boxes are relative to the box.
grid = run_on(@markovjunior(begin
    @within_box pixel(min=5, size=4) begin
        @fill 'Y' uv(center=0.5, size=0)
        @rewrite (area/2) b => w
    end
end), fill(c('b'), 8, 8))
@bp_check(count(==(c('w')), grid) == 8, grid)
@bp_check(all(I -> (I[1] >= 5) && (I[2] >= 5), findall(==(c('w')) , grid)), grid)
standalone_Y = only(cells_of(run_on(@markovjunior(begin
    @fill 'Y' uv(center=0.5, size=0)
end), fill(c('b'), 4, 4)), 'Y'))
@bp_check(cells_of(grid, 'Y') == Set([ standalone_Y .+ 4 ]), grid)

# `uv` boxes split the grid evenly, for both `@fill` and `@within_box`.
for (box, expected) in (("uv(min=0, size=0.5)", 1:4), ("uv(min=0.5, size=0.5)", 5:8),
                        ("uv(min=0.5, max=1)", 5:8), ("uv(min=0.25, size=0.25)", 3:4))
    filled = run_on(markov_algo_parse("@markovjunior 1 begin\n @fill 'w' $box\nend"), fill(c('b'), 8))
    @bp_check(findall(==(c('w')), filled) == collect(expected), box, ": ", filled)
    scoped = run_on(markov_algo_parse("@markovjunior 1 begin\n @within_box $box begin\n @rewrite b => w\n end\nend"),
                    fill(c('b'), 8))
    @bp_check(findall(==(c('w')), scoped) == collect(expected), box, ": ", scoped)
end

# A `uv` box covering the whole grid.
grid = run_on(@markovjunior(begin
    @within_box uv(min=0, max=1) begin
        @rewrite b => w
    end
end), fill(c('b'), 3, 5))
@bp_check(all(==(c('w')), grid), grid)

# A scalar box applies to every axis.
grid = run_on(@markovjunior(begin
    @within_box pixel(min=3, size=2) begin
        @rewrite b => w
    end
end), fill(c('b'), 4, 4, 4))
@bp_check(count(==(c('w')), grid) == 8, grid)
@bp_check(all(==(c('w')), grid[3:4, 3:4, 3:4]), grid)

# Nested sequences and nested scopes stay within the box.
grid = run_on(@markovjunior(begin
    @within_box pixel(min=3, size=4) begin
        @sequence begin
            @rewrite b => w
        end
        @within_box pixel(min=2, size=2) begin
            @rewrite w => R
        end
    end
end), fill(c('b'), 8, 8))
@bp_check(cells_of(grid, 'w') ∪ cells_of(grid, 'R') == Set((x, y) for x in 3:6, y in 3:6), grid)
@bp_check(cells_of(grid, 'R') == Set((x, y) for x in 4:5, y in 4:5), grid)

# An empty scope is cosmetic, so 'repeat' should look past it.
grid = run_on(@markovjunior(begin
    @sequence repeat begin
        @within_box pixel(min=1, size=1) begin end
        @rewrite 1 b => w
    end
end), fill(c('b'), 3, 3))
@bp_check(all(==(c('w')), grid), grid)


########
# Biases

grid_with_Y = let g = fill(c('b'), 8, 8)
    g[7, 7] = c('Y')
    g
end

# By default, inherited biases don't apply inside the scope;
#    the scope's own biases are local to it.
COUNTING_BIAS_N_INITIALIZED[] = 0
for seed in 1:5
    grid = run_on(@markovjunior(begin
        @sequence begin
            @within_box pixel(min=5, size=4) begin
                @rewrite 1 b => R field(Y)
                @rewrite 1 b => G
            end
            @rewrite 1 b => M
        end begin
            field(-Y)
            counting()
        end
    end), grid_with_Y; seed=seed)
    @bp_check(only(cells_of(grid, 'R')) in neighbors(7, 7), seed, ": ", grid)
    @bp_check(all(I -> (I[1] >= 5) && (I[2] >= 5), findall(==(c('G')), grid)), seed, ": ", grid)
    # The outer field pushes away from Y, which the scope doesn't see.
    @bp_check(only(cells_of(grid, 'M')) in Set([ (1, 1) ]), seed, ": ", grid)
end
@bp_check(COUNTING_BIAS_N_INITIALIZED[] == 5, COUNTING_BIAS_N_INITIALIZED[])
@bp_check(isempty(COUNTING_BIAS_LIVE_STATES), COUNTING_BIAS_LIVE_STATES)

# With `inherit_biases`, inherited biases get a separate state for the scope.
COUNTING_BIAS_N_INITIALIZED[] = 0
for seed in 1:5
    grid = run_on(@markovjunior(begin
        @sequence begin
            @within_box inherit_biases pixel(min=5, size=4) begin
                @rewrite 1 b => R
            end
            @rewrite 1 b => M
        end begin
            field(Y)
            counting()
        end
    end), grid_with_Y; seed=seed)
    @bp_check(only(cells_of(grid, 'R')) in neighbors(7, 7), seed, ": ", grid)
    @bp_check(only(cells_of(grid, 'M')) in neighbors(7, 7), seed, ": ", grid)
end
@bp_check(COUNTING_BIAS_N_INITIALIZED[] == 5 * 2, COUNTING_BIAS_N_INITIALIZED[])
@bp_check(isempty(COUNTING_BIAS_LIVE_STATES), COUNTING_BIAS_LIVE_STATES)

# Nested scopes which inherit biases don't disturb the states of the scope around them.
COUNTING_BIAS_N_INITIALIZED[] = 0
grid = run_on(@markovjunior(begin
    @sequence begin
        @within_box inherit_biases pixel(min=1, size=6) begin
            @within_box inherit_biases pixel(min=2, size=2) begin
                @rewrite 1 b => R
            end
            @rewrite 1 b => G
        end
        @rewrite 1 b => M
    end counting()
end), fill(c('b'), 8, 8))
@bp_check(count(==(c('R')), grid) == 1, grid)
@bp_check(count(==(c('G')), grid) == 1, grid)
@bp_check(count(==(c('M')), grid) == 1, grid)
@bp_check(COUNTING_BIAS_N_INITIALIZED[] == 3, COUNTING_BIAS_N_INITIALIZED[])
@bp_check(isempty(COUNTING_BIAS_LIVE_STATES), COUNTING_BIAS_LIVE_STATES)

# Without `inherit_biases`, live inherited biases still hear about the scope's changes.
for seed in 1:5
    grid = run_on(@markovjunior(begin
        @sequence begin
            @within_box pixel(min=1, size=8) begin
                @fill 'w' pixel(min=8, size=1)
                @fill 'Y' pixel(min=1, size=1)
            end
            @rewrite 1 b => M
        end field(Y, live)
    end), let g = fill(c('b'), 8, 8); g[8, 8] = c('Y'); g end; seed=seed)
    @bp_check(only(cells_of(grid, 'M')) in neighbors(1, 1), seed, ": ", grid)
end

# Canceling inside a scope shouldn't leak bias states or leave the grid as a view.
for inherit in (false, true)
    mj = markov_algo_parse("""@markovjunior begin
        @sequence begin
            @within_box $(inherit ? "inherit_biases" : "") pixel(min=2, size=2) begin
                @rewrite b => w counting()
            end
        end counting()
    end""")
    ticking = MarkovTickSettings()
    ch = markov_algo_run(mj, (4, 4), ticking, seeds=0xaabbcc, allocator=MJ.MarkovAllocatorHeapReused())
    n_ticks = 0
    while true
        event = markov_algo_next(ch)
        if event isa Int
            n_ticks += 1
            (n_ticks == 1) && (ticking.cancel_algo = true)
        elseif event === MJ.TAG_ALGO_CANCELED
            break
        elseif event === MJ.TAG_ALGO_COMPLETED
            error("Algorithm should have been canceled (inherit=", inherit, ")")
        end
    end
    markov_algo_cleanup(ch)
    @bp_check(isempty(COUNTING_BIAS_LIVE_STATES), "inherit=", inherit, ": ", COUNTING_BIAS_LIVE_STATES)
end


##########################
# Reallocation inside scopes

# A custom Op must say whether it reallocates before it can be scoped.
markov_algo_parse("@markovjunior begin\n @test_unqueried\nend")
check_parse_fails("@markovjunior begin\n @within_box pixel(min=1, size=2) begin\n @test_unqueried\n end\nend")

# A custom Op that resizes the grid without saying so at parse-time is caught at run-time.
mj = @markovjunior begin
    @within_box pixel(min=1, size=2) begin
        @test_sneaky_resize
    end
end
for allocator in (MJ.MarkovAllocatorHeap(), MJ.MarkovAllocatorHeapReused())
    failed = try
        ch = markov_algo_run(mj, (4, 4), MarkovTickSettings(skip_tagged_events = MJ.SkipTaggedEvents.all_skippable),
                             allocator=allocator)
        markov_algo_complete(identity, ch)
        false
    catch e
        true
    end
    @bp_check(failed, "Reallocating inside a scope should have failed with ", typeof(allocator))
end


##############
# Parse errors

check_parse_fails("@markovjunior begin\n @within_box pixel(min=1, size=2) begin\n @upscale (2)\n end\nend")
check_parse_fails("@markovjunior begin\n @within_box pixel(min=1, size=2) begin\n @sequence begin\n @downscale (2)\n end\n end\nend")
check_parse_fails("@markovjunior begin\n @within_box begin\n @rewrite b => w\n end\nend")
check_parse_fails("@markovjunior begin\n @within_box everywhere begin\n @rewrite b => w\n end\nend")
check_parse_fails("@markovjunior begin\n @within_box pixel(min=1, width=2) begin\n @rewrite b => w\n end\nend")
check_parse_fails("@markovjunior begin
 @within_box pixel(min=1, size=nope) begin
 @rewrite b => w
 end
end")
check_parse_fails("@markovjunior begin
 @within_box pixel(min=1, min=2) begin
 @rewrite b => w
 end
end")
check_parse_fails("@markovjunior 1 begin\n @within_box pixel(min=(1, 1), size=2) begin\n @rewrite b => w\n end\nend")

end)()
