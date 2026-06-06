# Tests the ability to add onto the algorithm with pragmas and custom data.

using Setfield

"
Orders pixels by X, then Y, then Z, etc, where the first pixel is chosen first.
Ordering can be flipped by giving a named bool to the algo's `add_ons`.
"
struct CustomBias <: MJ.AbstractMarkovBias end
const ADDON_CUSTOM_BIAS_FLIPPED_KEY = :flip_custom_bias
function MJ.markov_bias_initialize(::CustomBias, grid::MJ.CellGrid, ::PRNG, ctx::MJ.MarkovBiasContext)
    return haskey(ctx.add_ons, ADDON_CUSTOM_BIAS_FLIPPED_KEY) &&
           convert(Bool, ctx.add_ons[ADDON_CUSTOM_BIAS_FLIPPED_KEY])
end
MJ.markov_bias_state_type(::Type{CustomBias}) = Bool
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
MJ.markov_bias_calculate(::CustomBias, is_flipped::Bool, grid::MJ.CellGrid{N},
                         line::MJ.CellLine{N}, ::PRNG) where {N} = impl_custom_bias(line.start_cell, vsize(grid), is_flipped)
MJ.markov_bias_calculate(::CustomBias, is_flipped::Bool, grid::MJ.CellGrid{N},
                         area::MJ.CellRegion{N}, ::PRNG) where {N} = impl_custom_bias(min_inclusive(area.b), vsize(grid), is_flipped)
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
    test_parsed_mj = parse_markovjunior(MJ.dsl_string(test_mj))
    @bp_check(test_mj == test_parsed_mj,
              "Custom Bias parsing failed!\n\nSource:\n", test_mj,
              "\nDest:\n", test_parsed_mj)

    for i in 1:20
        # Halfway through, flip the ordering.
        flipped::Bool = (i > 10)
        if flipped
            test_mj.add_ons[ADDON_CUSTOM_BIAS_FLIPPED_KEY] = true
        end
        log_info = (" i=", i, flipped ? " (flipped)" : "")

        size_1 = (2, 3, 4, 3, 2)
        state_1 = MJ.markov_algo_start(test_mj, size_1, (i, 0xaabbcc))
        MJ.markov_algo_step(test_mj, state_1, 1)
        base_idx = flipped ? size_1 : (1, 1, 1, 1, 1)
        idx_desc = flipped ? "last" : "first"
        rest_range = flipped ? (1:(prod(size_1)-1)) : (2:(prod(size_1)))
        log_info_1 = (log_info..., "  grid:\n", MJ.markov_algo_grid(state_1))
        @bp_check(MJ.markov_algo_grid(state_1)[base_idx...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "1-step didn't write 'w' at the ", idx_desc, " pixel!", log_info_1...)
        @bp_check(all(u -> (u==MJ.CELL_CODE_BY_CHAR['b']),
                      MJ.markov_algo_grid(state_1)[rest_range]),
                  "1-step didn't leave all other cells as 'b'!", log_info_1...)
        close(state_1, test_mj)

        size_2 = (2, 3, 4, 3, 2)
        state_2 = MJ.markov_algo_start(test_mj, size_2, (i, 0xafabbcfc))
        MJ.markov_algo_step(test_mj, state_2, 2)
        base_idx = flipped ? size_2 : (1, 1, 1, 1, 1)
        idx_desc = flipped ? "last" : "first"
        rest_range = flipped ? (1:(prod(size_2)-2)) : (3:(prod(size_2)))
        idx_sign = flipped ? -1 : 1
        log_info_2 = (log_info..., "  grid:\n", MJ.markov_algo_grid(state_2))
        @bp_check(MJ.markov_algo_grid(state_2)[base_idx...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "2-step didn't write 'w' at the ", idx_desc, " pixel!", log_info_2...)
        @bp_check(MJ.markov_algo_grid(state_2)[@set(base_idx[1] += idx_sign)...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "2-step didn't write 'w' at the ", idx_desc, "-take-1 pixel!", log_info_2...)
        @bp_check(all(u -> (u==MJ.CELL_CODE_BY_CHAR['b']),
                      MJ.markov_algo_grid(state_2)[rest_range]),
                  "2-step didn't leave all other cells as 'b'!", log_info_2...)
        close(state_2, test_mj)

        size_3 = (2, 2, 2)
        state_3 = MJ.markov_algo_start(test_mj, size_3, (i, 0xaf67bcfc))
        MJ.markov_algo_step(test_mj, state_3, 3)
        base_idx = flipped ? size_3 : (1, 1, 1)
        idx_desc = flipped ? "last" : "first"
        idx_sign = flipped ? -1 : 1
        rest_range = flipped ? (1:(prod(size_3)-3)) : (4:(prod(size_3)))
        log_info_3 = (log_info..., "  grid:\n", MJ.markov_algo_grid(state_3))
        @bp_check(MJ.markov_algo_grid(state_3)[base_idx...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "3-step didn't write 'w' at the ", idx_desc, " pixel!", log_info_3...)
        @bp_check(MJ.markov_algo_grid(state_3)[@set(base_idx[1] += idx_sign)...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "3-step didn't write 'w' at the ", idx_desc, "-take-1 pixel!", log_info_3...)
        @bp_check(MJ.markov_algo_grid(state_3)[@set(base_idx[2] += idx_sign)...] == MJ.CELL_CODE_BY_CHAR['w'],
                  "3-step didn't write 'w' at the ", idx_desc, "-take-2 pixel!", log_info_3...)
        @bp_check(all(u -> (u==MJ.CELL_CODE_BY_CHAR['b']),
                      MJ.markov_algo_grid(state_3)[rest_range]),
                  "3-step didn't leave all other cells as 'b'!", log_info_3...)
        close(state_3, test_mj)
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
function MJ.markov_op_initialize(co::CustomOp, ::MJ.CellGrid{N}, ::PRNG, ctx::MJ.MarkovOpContext) where {N}
    # Look for pragma statements that force checkerboard values for the whole grid.
    a = co.a
    b = co.b
    for (name::Symbol, values::Vector{Any}) in ctx.pragmas
        if name == PRAGMA_CUSTOM_OP_OVERRIDE_A
            (override, ) = values
            a = MJ.CELL_CODE_BY_SYMBOL[override]
        elseif name == PRAGMA_CUSTOM_OP_OVERRIDE_B
            (override, ) = values
            b = MJ.CELL_CODE_BY_SYMBOL[override]
        end
    end

    return CustomOpState{N}(a, b, one(MJ.CellIdx{N}))
end
function MJ.markov_op_iterate(op::CustomOp, state::CustomOpState{N}, grid::MJ.CellGrid{N}, ::PRNG, ::MJ.MarkovOpContext) where {N}
    # Write the next checkerboard cell.
    grid[state.next_pos] = if iseven(sum(state.next_pos - 1))
        state.overridden_a
    else
        state.overridden_b
    end

    # Move on.
    next_pos = state.next_pos
    axis = 1
    while (axis <= N) && (next_pos[axis] == size(grid, axis))
        @set! next_pos[axis] = 1
        axis += 1
    end
    if axis > N
        return nothing
    else
        @set! next_pos[axis] += 1
        return CustomOpState{N}(state.overridden_a, state.overridden_b, next_pos)
    end
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
    test_parsed_mj = parse_markovjunior(markov_algo_to_string(test_mj1))
    @bp_check(test_mj1 == test_parsed_mj,
              "Custom Op parsing failed!\n\nSource:\n", test_mj1,
                "\nDest:\n", test_parsed_mj)
    state_1 = MJ.markov_algo_start(test_mj1, (3, 4, 2), 0xaabbcc)
    MJ.markov_algo_finish(test_mj1, state_1)
    grid_1 = MJ.markov_algo_grid(state_1)
    @bp_check(grid_1 == getindex.(Ref(MJ.CELL_CODE_BY_SYMBOL), [
        :R :G :R :G
        :G :R :G :R
        :R :G :R :G ;;;
        :G :R :G :R
        :R :G :R :G
        :G :R :G :R
    ]), "Actual: ", grid_1)

    test_mj2 = @markovjunior 'b' begin
        @pragma cust_a_forced B
        @cust R G
        @fill 'g' uv(min=0, max=1) -B
    end
    test_parsed_mj = parse_markovjunior(markov_algo_to_string(test_mj2))
    @bp_check(test_mj2 == test_parsed_mj,
              "Custom Op parsing #2 failed!\n\nSource:\n", test_mj2,
                "\nDest:\n", test_parsed_mj)
    state_2 = MJ.markov_algo_start(test_mj2, (10, ), 0x80110110)
    MJ.markov_algo_finish(test_mj2, state_2)
    grid_2 = MJ.markov_algo_grid(state_2)
    @bp_check(grid_2 == getindex.(Ref(MJ.CELL_CODE_BY_SYMBOL), [
        :B, :g, :B, :g, :B, :g, :B, :g, :B, :g
    ]), "Actual: ", grid_2)
end)()