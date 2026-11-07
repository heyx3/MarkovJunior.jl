(() -> begin

# Test for sequence handling of Ops that are cosmetic vs those which fail to change anything.

# Runs the given algorithm to completion on a black 3x3 grid, skipping all tagged events
function run_skipping_events(mj::MarkovAlgorithm)::Matrix{UInt8}
    test_parsed_mj = markov_algo_parse(markov_algo_to_string(mj))
    @bp_check(mj == test_parsed_mj,
                "Parsing failed!\n\nSource:\n", mj,
                "\nDest:\n", test_parsed_mj,
                "\nSource as string: ", markov_algo_to_string(mj))

    # Use a normal heap allocator so we can steal the grid after the algo is done.
    state_ch = markov_algo_run(mj, (3, 3),
                                MarkovTickSettings(skip_tagged_events = MJ.SkipTaggedEvents.all_skippable),
                                seeds=0xaabbcc,
                                allocator=MJ.MarkovAllocatorHeap())
    return markov_algo_complete(identity, state_ch)
end

const_grid(c::Char) = fill(MJ.CELL_CODE_BY_CHAR[c], 3, 3)

# Read-only events are cosmetic, so 'repeat' should look past them.
grid = run_skipping_events(@markovjunior 2 begin
    @sequence repeat begin
        @event A
        @rewrite 1 b=>w
    end
end)
@bp_check(grid == const_grid('w'), "Actual: ", grid)

# A skipped mutable event fails to change anything, so 'repeat' should exit immediately.
grid = run_skipping_events(@markovjunior 2 begin
    @sequence repeat begin
        @event !A
        @rewrite 1 b=>w
    end
end)
@bp_check(grid == const_grid('b'), "Actual: ", grid)

# Sequences that can never do anything should still finish.
grid = run_skipping_events(@markovjunior 2 begin
    @sequence repeat begin
        @event A
    end
    @sequence repeat begin end
    @sequence begin end
    @rewrite 1 b=>w
end)
@bp_check(count(==(MJ.CELL_CODE_BY_CHAR['w']), grid) == 1, "Actual: ", grid)

end)()
