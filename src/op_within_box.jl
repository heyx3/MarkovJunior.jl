#NOTE: this file takes some things from op_sequence.jl, and should be included after it.

"Runs a sequence of Ops within a box-shaped part of the grid, as if that part were the entire grid"
struct MarkovOpWithinBox <: AbstractMarkovOp
    inherit_biases::Bool
    box::CalcBox # Same as `@fill`
    ops::Vector{AbstractMarkovOp}
end

function markov_algo_run(op::MarkovOpWithinBox,
                         algo::MarkovAlgorithm, algo_state::AlgoState,
                         inherited_biases::Tuple{Vararg{AbstractMarkovBias}},
                         inherited_bias_states::Tuple,
                         grid::CellGrid{NGrid} = algo_state.grid
                        )::OpChangedGrid where {NGrid}
    if isempty(op.ops)
        return nothing
    end

    area = eval_calc_box(op.box, grid, algo_state.rng, algo.calc_macros)
    if isnothing(area)
        @logic_logln("@within_box's area is empty")
        markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY + 1)
        return false
    end
    (area_min, area_max) = (min_inclusive(area), max_inclusive(area))
    # Views of views get flattened by Julia into one view of the original array,
    #    so nested scopes don't create new grid types.
    sub_grid = view(grid, ntuple(i -> area_min[i]:area_max[i], Val(NGrid))...)

    @logic_logln("Running ", length(op.ops), " ops within ", area_min, " to ", area_max)
    result = markov_opalgo_on_subgrid(algo, algo_state, sub_grid, area,
                                      op.inherit_biases, inherited_biases, inherited_bias_states
                                     ) do inner_biases, inner_bias_states
        # The grid can't be reallocated in here, so the bias states never get replaced.
        inner_result::OpChangedGrid = nothing
        for inner_op in op.ops
            inner_result = combine_op_results(inner_result, markov_algo_run(
                inner_op, algo, algo_state,
                inner_biases, inner_bias_states
            ))
        end
        return inner_result
    end

    markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY + 1)
    return result
end

# The inner Ops are rejected at parse-time if they reallocate.
markov_op_reallocates_grid(op::MarkovOpWithinBox) = any(markov_op_reallocates_grid, op.ops)


dsl_format(op::MarkovOpWithinBox) = string(
    "@within_box ",
    op.inherit_biases ? "inherit_biases " : "",
    dsl_format(op.box),
    " begin\n    ",
    iter_join(Iterators.map(dsl_format, op.ops), "\n    ")...,
    "\nend"
)

function parse_markovjunior_op(::Val{Symbol("@within_box")}, inputs::MacroParserInputs,
                               loc, expr_args, expr_line)
    args = collect(Any, expr_args)
    inherit_biases = !isempty(args) && (args[1] == :inherit_biases)
    inherit_biases && popfirst!(args)
    if length(args) != 2
        raise_parse_error(loc, inputs,
                          "Expected `@within_box [inherit_biases] S(A=N, B=M) begin ... end`; got ",
                            length(expr_args), " arguments")
    end
    (box_expr, ops_expr) = args

    box = with_parser_stacktrace(inputs, "Box") do
        parse_markovjunior_calc_box(inputs, loc, box_expr, CALC_BOX_ENV)
    end

    # Parse the inner Ops.
    # Without inherited biases, they shouldn't be validated against the outer biases either.
    outer_bias_stack = inputs.bias_stack
    inherit_biases || (inputs.bias_stack = Stack{Vector{AbstractMarkovBias}}(16))
    ops = try
        with_parser_stacktrace(inputs, "Inner Ops") do
            parse_markovjunior_sequence_main_block(ops_expr, loc, inputs)
        end
    finally
        inputs.bias_stack = outer_bias_stack
    end
    check_markovjunior_scoped_ops(inputs, loc, "@within_box", ops)

    return MarkovOpWithinBox(inherit_biases, box, ops)
end
