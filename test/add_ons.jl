# Tests the ability to add onto the algorithm with pragmas and custom data.

using Setfield

"Orders pixels by X, then Y, then Z, etc, where the first pixel is chosen first"
struct CustomBias <: MJ.AbstractMarkovBias end
const ADDON_CUSTOM_BIAS_FLIPPED_KEY = :flip_custom_bias
function MJ.markov_bias_initialize(::CustomBias, grid::MJ.CellGrid, ::PRNG, ctx::MJ.MarkovBiasContext)
    return haskey(ctx.add_ons, ADDON_CUSTOM_BIAS_FLIPPED_KEY) &&
           convert(Bool, ctx.add_ons[ADDON_CUSTOM_BIAS_FLIPPED_KEY])
end
function impl_custom_bias(pos::MJ.CellIdx{N}, size::MJ.CellIdx{N},
                          is_flipped::Bool
                          ;
                          axis::Int = 1, slice_size::Int = 1) where {N}
    return ((is_flipped ? 1 : -1) * pos[axis] * slice_size) + (if axis == N
        is_flipped ? 0 : prod(size)
    else
        impl_custom_bias(pos, size, is_flipped, axis + 1, slice_size * size[axis])
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
#


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
            a = MJ.CELL_CODE_BY_CHAR[string(override)[1]]
        elseif name == PRAGMA_CUSTOM_OP_OVERRIDE_B
            (override, ) = values
            b = MJ.CELL_CODE_BY_CHAR[string(override)[1]]
        end
    end

    return CustomOpState{N}(a, b, one(MJ.CellIdx{N}))
end
function MJ.markov_op_iterate(op::CustomOp, state::CustomOpState{N}, grid::MJ.CellGrid{N}, ::PRNG, ::MJ.MarkovOpContext) where {N}
    # Write the next checkerboard cell.
    grid[state.next_pos] = if isodd(sum(state.next_pos - 1))
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
MJ.dsl_string(co::CustomOp) = "@cust $(dsl_string(co.a)) $(dsl_string(co.b))"
function MJ.parse_markovjunior_op(::Val{Symbol("@cust")},
                                  inputs::MJ.MacroParserInputs,
                                  loc, expr_args, full_line)
    (expr_a, expr_b) = expr_args
    return CustomOp(
        MJ.CELL_CODE_BY_SYMBOL[expr_a],
        MJ.CELL_CODE_BY_SYMBOL[expr_b]
    )
end


#TODO: Run algorithms using the above and test that they work