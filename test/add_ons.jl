# Tests the ability to extend the algorithm with:
#   * A custom Bias reading from 'add_ons'
#   * A custom Op reading from 'pragmas'
#   * A custom Op reading/writing to 'data_store'

#TODO: Test custom events

using Setfield


"
Orders pixels by X, then Y, then Z, etc, where the first pixel is chosen first.
Ordering can be flipped by giving a named bool to the algo's `add_ons`.
"
struct CustomBias <: MJ.AbstractMarkovBias end
const ADDON_CUSTOM_BIAS_FLIPPED_KEY = :flip_custom_bias
function MJ.markov_bias_initialize(::CustomBias, algo::MJ.MarkovAlgorithm, algo_state::MJ.AlgoState)
    return haskey(algo.add_ons, ADDON_CUSTOM_BIAS_FLIPPED_KEY) &&
           convert(Bool, algo.add_ons[ADDON_CUSTOM_BIAS_FLIPPED_KEY])
end
function impl_custom_bias(pos::MJ.CellIdx{N}, size::MJ.CellIdx{N},
                          is_flipped::Bool
                          ;
                          axis::Int = 1,
                          slice_size::Int = 1) where {N}
    return ((is_flipped ? 1 : -1) * (pos[axis]-1) * slice_size) + (if axis == N
        is_flipped ? 0 : prod(size)
    else
        impl_custom_bias(pos, size, is_flipped,
                         axis = axis + 1,
                         slice_size = slice_size * size[axis])
    end)
end
MJ.markov_bias_calculate(::CustomBias, is_flipped::Bool,
                         algo::MJ.MarkovAlgorithm, algo_state::MJ.AlgoState,
                         line::MJ.CellLine{N}
                        ) where {N} = impl_custom_bias(line.start_cell, vsize(algo_state.grid), is_flipped)
MJ.markov_bias_calculate(::CustomBias, is_flipped::Bool,
                         algo::MJ.MarkovAlgorithm, algo_state::MJ.AlgoState,
                         area::MJ.CellRegion{N}
                        ) where {N} = impl_custom_bias(min_inclusive(area.b), vsize(algo_state.grid), is_flipped)
MJ.dsl_string(::CustomBias) = "custom()"
MJ.parse_markovjunior_bias(::Val{:custom}, inputs::MJ.MacroParserInputs,
                           loc::LineNumberNode, args
                          ) = CustomBias()

# Test the bias.
# Since pixel selection is normally random, we'll run the test multiple times.
(() -> begin
    test_mj = MJ.@markovjunior 'b' begin
        @rewrite b => w custom()
    end
    test_parsed_mj = markov_algo_parse(MJ.dsl_string(test_mj))
    @bp_check(test_mj == test_parsed_mj,
              "Custom Bias parsing failed!\n\nSource:\n", test_mj,
              "\nDest:\n", test_parsed_mj)
    test_tick_settings = MJ.MarkovTickSettings(
        # The tick for length-1 rewrite ops is normally compiled out.
        MJ.STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY - 1
    )

    for i in 1:20
        # Halfway through, flip the ordering.
        flipped::Bool = (i > 10)
        if flipped
            test_mj.add_ons[ADDON_CUSTOM_BIAS_FLIPPED_KEY] = true
        end
        log_info = (" i=", i, flipped ? " (flipped)" : "")

        testable_tick_output(x::Union{Int, Symbol}) = x
        testable_tick_output(::MJ.CellGrid{N}) where {N} = MJ.CellGrid{N}

        size_1 = (2, 3, 4, 3, 2)
        state_ch_1 = markov_algo_run(test_mj, size_1, test_tick_settings, seeds=(i, 0xaabbcc))
        first_ticks_1 = map(i->markov_algo_next(state_ch_1), 1:2)
        @bp_check(testable_tick_output.(first_ticks_1) == [
                      MJ.CellGrid{length(size_1)},
                      MJ.TAG_ALGO_STARTING
                  ], first_ticks_1)
        grid_1 = first_ticks_1[1]
        let tt = markov_algo_next(state_ch_1)
            @bp_check(tt == MJ.STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY-1,
                      typeof(tt), "(", tt, ")")
        end
        base_idx = flipped ? size_1 : (1, 1, 1, 1, 1)
        idx_desc = flipped ? "last" : "first"
        rest_range = flipped ? (1:(prod(size_1)-1)) : (2:(prod(size_1)))
        log_info_1 = (log_info..., "  rest_range=", rest_range, "  grid:\n", Int.(grid_1))
        @bp_check(grid_1[base_idx...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "1-step didn't write 'w' at the ", idx_desc, " pixel!", log_info_1...)
        @bp_check(all(u -> (u==MJ.CELL_CODE_BY_CHAR['b']),
                      grid_1[rest_range]),
                  "1-step didn't leave all other cells as 'b'/", MJ.CELL_CODE_BY_CHAR['b'],
                    "!", log_info_1...)
        close(state_ch_1)

        size_2 = (2, 3, 4, 3, 2)
        state_ch_2 = markov_algo_run(test_mj, size_2, test_tick_settings, seeds=(i, 0xafabbcfc))
        first_ticks_2 = map(i->markov_algo_next(state_ch_2), 1:2)
        @bp_check(testable_tick_output.(first_ticks_2) == [
                      MJ.CellGrid{length(size_2)},
                      MJ.TAG_ALGO_STARTING
                  ], first_ticks_2)
        grid_2 = first_ticks_2[1]
        for i in 1:2
            let tt = markov_algo_next(state_ch_2)
                @bp_check(tt == MJ.STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY-1,
                          i, ": ", typeof(tt), "(", tt, ")")
            end
        end
        base_idx = flipped ? size_2 : (1, 1, 1, 1, 1)
        idx_desc = flipped ? "last" : "first"
        rest_range = flipped ? (1:(prod(size_2)-2)) : (3:(prod(size_2)))
        idx_sign = flipped ? -1 : 1
        log_info_2 = (log_info..., "  grid:\n", grid_2)
        @bp_check(grid_2[base_idx...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "2-step didn't write 'w' at the ", idx_desc, " pixel!", log_info_2...)
        @bp_check(grid_2[@set(base_idx[1] += idx_sign)...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "2-step didn't write 'w' at the ", idx_desc, "-take-1 pixel!", log_info_2...)
        @bp_check(all(u -> (u==MJ.CELL_CODE_BY_CHAR['b']),
                      grid_2[rest_range]),
                  "2-step didn't leave all other cells as 'b'!", log_info_2...)
        close(state_ch_2)

        # This test will use a random seed; if this is the one that fails
        #    then perhaps the other seeds just got lucky.
        size_3 = (2, 2, 2)
        state_ch_3 = markov_algo_run(test_mj, size_3, test_tick_settings)
        first_ticks_3 = map(i->markov_algo_next(state_ch_3), 1:2)
        @bp_check(testable_tick_output.(first_ticks_3) == [
                      MJ.CellGrid{length(size_3)},
                      MJ.TAG_ALGO_STARTING
                  ], first_ticks_3)
        grid_3 = first_ticks_3[1]
        for i in 1:3
            let tt = markov_algo_next(state_ch_3)
                @bp_check(tt == MJ.STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY-1,
                          i, ": ", typeof(tt), "(", tt, ")")
            end
        end
        base_idx = flipped ? size_3 : (1, 1, 1)
        idx_desc = flipped ? "last" : "first"
        idx_sign = flipped ? -1 : 1
        rest_range = flipped ? (1:(prod(size_3)-3)) : (4:(prod(size_3)))
        log_info_3 = (log_info..., "  grid:\n", grid_3)
        @bp_check(grid_3[base_idx...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "3-step didn't write 'w' at the ", idx_desc, " pixel!", log_info_3...)
        @bp_check(grid_3[@set(base_idx[1] += idx_sign)...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "3-step didn't write 'w' at the ", idx_desc, "-take-1 pixel!", log_info_3...)
        @bp_check(grid_3[@set(base_idx[2] += idx_sign)...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "3-step didn't write 'w' at the ", idx_desc, "-take-2 pixel!", log_info_3...)
        @bp_check(all(u -> (u==MJ.CELL_CODE_BY_CHAR['b']),
                      grid_3[rest_range]),
                  "3-step didn't leave all other cells as 'b'!", log_info_3...)
        close(state_ch_3)
    end
end)()


"Outputs a checkerboard pattern"
struct CustomOp <: MJ.AbstractMarkovOp
    a::UInt8
    b::UInt8
end
struct CustomOpState{N}
    overridden_a::UInt8
    overridden_b::UInt8
    next_pos::MJ.CellIdx{N}
end
const PRAGMA_CUSTOM_OP_OVERRIDE_A = :cust_a_forced
const PRAGMA_CUSTOM_OP_OVERRIDE_B = :cust_b_forced
function MJ.markov_algo_run(co::CustomOp, algo::MJ.MarkovAlgorithm, algo_state::MJ.AlgoState,
                            inherited_biases::Tuple{Vararg{MJ.AbstractMarkovBias}},
                            inherited_bias_states::Tuple{Vararg{Any}},
                            # Turn some data into compile-time knowledge.
                            ::Val{N} = Val(ndims(algo_state.grid)),
                            ::Type{TGrid} = typeof(algo_state.grid)
                           )::Tuple{Bool, typeof(inherited_bias_states)} where {N, TGrid}
    # Look for pragma statements that force checkerboard values over the ones in the Op.
    a = co.a
    b = co.b
    if haskey(algo.pragmas_map, PRAGMA_CUSTOM_OP_OVERRIDE_A)
        (override, ) = algo.pragmas_map[PRAGMA_CUSTOM_OP_OVERRIDE_A][end]
        a = MJ.CELL_CODE_BY_SYMBOL[override]
    end
    if haskey(algo.pragmas_map, PRAGMA_CUSTOM_OP_OVERRIDE_B)
        (override, ) = algo.pragmas_map[PRAGMA_CUSTOM_OP_OVERRIDE_B][end]
        b = MJ.CELL_CODE_BY_SYMBOL[override]
    end

    grid::TGrid = algo_state.grid
    next_pos::MJ.CellIdx{N} = one(MJ.CellIdx{N})

    # Technically we are supposed to update biases after we're done changing pixels;
    #    even though we don't use biases in this test we may as well do it right.
    inherited_bias_states = MJ.markov_allocator_with_array(algo_state.allocator,
                                                           size(grid), eltype(grid)) do old_grid
    #begin
        old_grid .= grid

        for next_pos in one(MJ.CellIdx{N}) : convert(MJ.CellIdx{N}, vsize(grid))
            # Write the next checkerboard cell.
            grid[next_pos] = if iseven(sum(next_pos - 1))
                a
            else
                b
            end

            MJ.markov_algo_tick(algo_state, 2)
        end

        MJ.markov_bias_update.(
            inherited_biases, inherited_bias_states,
            Ref(algo), Ref(algo_state),
            Ref(BoxI{N}(min=one(MJ.CellIdx{N}), size=vsize(grid))),
            Ref(old_grid)
        )
    end

    MJ.markov_algo_tick(algo_state, MJ.STANDARD_END_OF_OP_TICK_PRIORITY)
    return (!isempty(grid), inherited_bias_states)
end
MJ.dsl_string(co::CustomOp) = "@cust $(MJ.dsl_string(co.a)) $(MJ.dsl_string(co.b))"
function MJ.parse_markovjunior_op(::Val{Symbol("@cust")},
                                  inputs::MJ.MacroParserInputs,
                                  loc, expr_args, full_line)
    (expr_a, expr_b) = expr_args
    return CustomOp(
        MJ.CELL_CODE_BY_SYMBOL[expr_a],
        MJ.CELL_CODE_BY_SYMBOL[expr_b]
    )
end

# Test the custom op.
(() -> begin
    test_mj1 = @markovjunior 'b' begin
        @cust R G
    end
    test_parsed_mj = markov_algo_parse(markov_algo_to_string(test_mj1))
    @bp_check(test_mj1 == test_parsed_mj,
              "Custom Op parsing failed!\n\nSource:\n", test_mj1,
                "\nDest:\n", test_parsed_mj,
                "\nSource as string: ", markov_algo_to_string(test_mj1))
    # Use a normal heap allocator so we can steal the grid after the algo is done.
    state_ch_1 = markov_algo_run(test_mj1, (3, 4, 2), seeds=0xaabbcc,
                                    allocator=MJ.MarkovAllocatorHeap())
    grid_1 = markov_algo_complete(identity, state_ch_1)
    @bp_check(grid_1 == getindex.(Ref(MJ.CELL_CODE_BY_SYMBOL), [
        :R :G :R :G
        :G :R :G :R
        :R :G :R :G ;;;
        :G :R :G :R
        :R :G :R :G
        :G :R :G :R
    ]), "Actual: ", grid_1)

    # This test uses a random seed, so if it fails then maybe the above tests are using lucky seeds!
    test_mj2 = @markovjunior 'b' begin
        @pragma cust_a_forced B
        @cust R G
        @fill 'g' uv(min=0, max=1) -B
    end
    test_parsed_mj = markov_algo_parse(markov_algo_to_string(test_mj2))
    @bp_check(test_mj2 == test_parsed_mj,
              "Custom Op parsing #2 failed!\n\nSource:\n", test_mj2,
                "\nDest:\n", test_parsed_mj,
                "\nSource as string: ", markov_algo_to_string(test_mj2))
    # Use a normal heap allocator so we can steal the grid after the algo is done.
    state_ch_2 = markov_algo_run(test_mj2, (10, ),
                                    allocator=MJ.MarkovAllocatorHeap())
    grid_2 = markov_algo_complete(identity, state_ch_2)
    @bp_check(grid_2 == getindex.(Ref(MJ.CELL_CODE_BY_SYMBOL), [
        :B, :g, :B, :g, :B, :g, :B, :g, :B, :g
    ]), "Actual: ", grid_2)
end)()


"
A custom Op that writes to a specific pixel.
The pixel location is parsed from the Op.
The output value is taken from the algorithm instance's `data_store`,
  starting at 0 and incrementing each time the Op is used.
"
struct CustomOp2{N} <: MJ.AbstractMarkovOp
    location::VecI{N}
end
function MJ.markov_algo_run(op::CustomOp2{N}, algo::MJ.MarkovAlgorithm, algo_state::MJ.AlgoState,
                            inherited_biases::NTuple{NB, MJ.AbstractMarkovBias},
                            inherited_bias_states::NTuple{NB, Any}
                           )::Tuple{Bool, typeof(inherited_bias_states)} where {N, NB}
    if any(op.location <= 0) || any(op.location > vsize(algo_state.grid))
        return (false, inherited_bias_states)
    end

    # Technically we are supposed to update biases after we're done changing pixels;
    #    even though we don't use biases in this test we may as well do it right.
    inherited_bias_states = MJ.markov_allocator_with_array(algo_state.allocator,
                                                           ntuple(i->1, ndims(algo_state.grid)),
                                                           eltype(algo_state.grid)) do old_grid
    #begin
        old_grid[1] = algo_state.grid[op.location]

        counter = get!(() -> Ref(zero(UInt8)), algo_state.data_store, :cust_op_2_counter)::Ref{UInt8}
        algo_state.grid[op.location] = counter[]
        counter[] += 1

        MJ.markov_bias_update(
            inherited_biases, inherited_bias_states,
            Ref(algo), Ref(algo_state),
            Ref(BoxI{N}(min=one(MJ.CellIdx{N}), size=vsize(algo_state.grid))),
            Ref(old_grid)
        )
    end

    return (true, inherited_bias_states)
end
MJ.dsl_string(co2::CustomOp2) = "@cust2 $(iter_join(co2.location, " ")...)"
function MJ.parse_markovjunior_op(::Val{Symbol("@cust2")},
                                  inputs::MJ.MacroParserInputs,
                                  loc, expr_args, full_line)
    @bp_check(all(x -> x isa Real, expr_args),
              "Coords should be numbers, but got ", expr_args)
    return CustomOp2(Vec{length(expr_args), Int32}(expr_args...))
end

# Test the custom Op.
(() -> begin
    test_mj = @markovjunior 2 'S' begin
        @cust2 1 2 # Writes 0
        @cust2 2 2 # Writes 1
        @cust2 1 1 # Writes 2
    end

    test_mj_str = markov_algo_to_string(test_mj)
    test_mj_parsed = markov_algo_parse(test_mj_str)
    @bp_check(test_mj == test_mj_parsed,
              "Failed to stringify/parse CustomOp2!",
                "\nSource:\n", test_mj,
                "\nSource to string: ", test_mj_str,
                "\nDest:\n", test_mj_parsed)

    # Test data_store multiple times, to be sure nothing is persistent.
    for i in 1:10
        # Use the basic Heap allocator so we can steal the grid at the end.
        state_ch = markov_algo_run(test_mj, (3, 4), seeds=0x0aa1bcde,
                                      allocator=MJ.MarkovAllocatorHeap())
        grid = markov_algo_complete(identity, state_ch)
        @bp_check(grid == UInt8[
            2 0 15 15
            15 1 15 15
            15 15 15 15
        ], "Iter ", i, ";  actual: ", grid)
    end
end)()