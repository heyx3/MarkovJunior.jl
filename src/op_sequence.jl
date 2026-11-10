const SequenceRepeatModeTag = Val{:repeat}

struct MarkovOpSequence{NBiases, TBiases<:NTuple{NBiases, AbstractMarkovBias}} <: AbstractMarkovOp
    ops::Vector{AbstractMarkovOp}
    threshold::Union{Nothing, SequenceRepeatModeTag, Threshold}
    biases::TBiases
end

function markov_algo_run(sequence::MarkovOpSequence{NSelfBiases},
                         algo::MarkovAlgorithm, algo_state::AlgoState,
                         inherited_biases::NTuple{NInheritedBiases, AbstractMarkovBias},
                         inherited_bias_states::NTuple{NInheritedBiases, Any},
                         grid::TGrid = algo_state.grid
                        )::OpChangedGrid where {
                            TGrid<:CellGrid,
                            NInheritedBiases, NSelfBiases
                        }
    # Set up the repetition counter.
    _repetitions_left = if isnothing(sequence.threshold)
        1
    elseif sequence.threshold isa SequenceRepeatModeTag
        sequence.threshold
    elseif sequence.threshold isa Threshold
        get_threshold(sequence.threshold, grid, algo_state.rng, algo.calc_macros)
    else
        error("Unhandled: ", typeof(sequence.threshold))
    end
    if isempty(sequence.ops) ||
       ((_repetitions_left isa Int) && (_repetitions_left <= 0))
        return isempty(sequence.ops) ? nothing : false
    end

    return markov_opbias_with_new_biases(sequence.biases, algo, algo_state,
                                         inherited_biases, inherited_bias_states
                                        ) do all_biases, all_bias_states
        repetitions_left = _repetitions_left
        known_grid_generation::Int = algo_state.grid_generation

        # Run the loop.
        @logic_logln("Starting sequence with threshold `", repetitions_left, "`")
        result::OpChangedGrid = nothing
        rep_i::Int = 0
        while (repetitions_left isa SequenceRepeatModeTag) || (repetitions_left > 0)
            if repetitions_left isa Int
                repetitions_left -= 1
            end
            rep_i += 1
            @logic_logln("REPETITION ", rep_i)

            # In 'repeat' mode, the first non-cosmetic op decides whether to keep going.
            pass_result::OpChangedGrid = nothing
            for op_i in 1:length(sequence.ops)
                @logic_logln("Inner op ", op_i, ":")
                @logic_tab_in()

                op_result = markov_algo_run(
                    sequence.ops[op_i], algo, algo_state,
                    all_biases, all_bias_states
                )
                # If the grid was reallocated, bias states need to be re-cached.
                if algo_state.grid_generation != known_grid_generation
                    known_grid_generation = algo_state.grid_generation
                    all_bias_states = markov_opbias_current_states(algo_state)
                end

                stop_early::Bool = false
                if (repetitions_left isa SequenceRepeatModeTag) && isnothing(pass_result) && (op_result === false)
                    @logic_logln("Inner op did nothing! The outer sequence will end now")
                    stop_early = true
                    repetitions_left = 0
                end
                pass_result = combine_op_results(pass_result, op_result)

                @logic_tab_out()
                stop_early && break
            end
            result = combine_op_results(result, pass_result)

            # If every op was cosmetic, repeating would never end.
            if (repetitions_left isa SequenceRepeatModeTag) && isnothing(pass_result)
                @logic_logln("Every inner op was cosmetic! The outer sequence will end now")
                repetitions_left = 0
            end
        end

        markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY + 1)
        return result
    end
end

markov_op_reallocates_grid(s::MarkovOpSequence) = any(markov_op_reallocates_grid, s.ops)

dsl_format(s::MarkovOpSequence) = string(
    "@sequence ",
    if isnothing(s.threshold)
        ""
    elseif s.threshold isa SequenceRepeatModeTag
        "repeat"
    elseif s.threshold isa Threshold
        dsl_format(s.threshold)
    else
        error("Unhandled: ", typeof(s.threshold))
    end,
     " ",
    "begin\n    ",
    iter_join(
        Iterators.map(dsl_format, s.ops),
        "\n    "
    )...,
    "\nend begin\n    ",
    iter_join(
        Iterators.map(dsl_format, s.biases),
        "\n    "
    )...,
    "\nend"
)

"
Parses the block of Ops for `@sequence` (and other Op which runs nested Ops).
`@calc` statements may be mixed in, but an error is thrown if they're the only kind of statement.
"
function parse_markovjunior_sequence_main_block(expr, src,
                                                inputs::MacroParserInputs,
                                               )::Vector{AbstractMarkovOp}
    function process_op(o_expr, i::Optional{Int})
        if !Base.isexpr(o_expr, :macrocall)
            raise_parse_error(src, inputs, "Operation isn't a macro-call: `", o_expr, "`")
        end
        (o_name, o_src, o_args...) = o_expr.args

        return with_parser_stacktrace(inputs, "Op $((exists(i) ? i : "")) $o_name") do
            return parse_markovjunior_op(Val(o_name), inputs, o_src, o_args, o_expr)
        end
    end

    if Base.isexpr(expr, :macrocall)
        if try_parse_markovjunior_calc!(inputs, src, expr)
            raise_parse_error(src, inputs,
                              "This block only has `@calc` statements, but it needs at least one Op")
        else
            return AbstractMarkovOp[  process_op(expr, nothing)  ]
        end
    elseif Base.isexpr(expr, :block)
        output = AbstractMarkovOp[ ]
        parse_markovjunior_block(expr.args) do inner_src, line
            if !try_parse_markovjunior_calc!(inputs, something(inner_src, src), line)
                push!(output, process_op(line, length(output)))
            end
        end

        if isempty(output)
            raise_parse_error(src, inputs, "This block has no Ops, e.g. only `@calc` statements??")
        else
            return output
        end
    else
        raise_parse_error(
            src, inputs,
            "Expected a block of operations; got `",
            if expr isa Expr
                tuple("`", expr.head, "` expression")
            else
                tuple(typeof(expr))
            end...
        )
    end
end

function parse_markovjunior_op(::Val{Symbol("@sequence")}, inputs::MacroParserInputs,
                               loc, expr_args, expr_line)
    (i_threshold, i_ops, i_bias) = if length(expr_args) == 3
        (1, 2, 3)
    elseif length(expr_args) == 2
        # Thresholds can be any calculation, so check whether the first argument is clearly the Ops instead:
        #   a single Op or a block of them.
        # A threshold wrapped in `threshold(...)` skips the guessing.
        first_is_ops = isnothing(unwrap_markovjunior_threshold(expr_args[1])) &&
                       Base.isexpr(expr_args[1], (:macrocall, :block))
        first_is_ops ? (nothing, 1, 2) : (1, 2, nothing)
    elseif length(expr_args) == 1
        (nothing, 1, nothing)
    else
        error("Expected 1-3 arguments; got ", length(expr_args))
    end

    generate_output(biases = AbstractMarkovBias[ ]) = MarkovOpSequence(
        parse_markovjunior_sequence_main_block(expr_args[i_ops], loc, inputs),
        with_parser_stacktrace(inputs, "Threshold") do
            if exists(i_threshold)
                if expr_args[i_threshold] == :repeat
                    SequenceRepeatModeTag()
                else
                    parse_markovjunior_threshold(inputs, loc, expr_args[i_threshold])
                end
            else
                nothing
            end
        end,
        Tuple(biases)
    )
    return if exists(i_bias)
        with_parsed_markovjunior_bias_statement(generate_output, inputs, loc, expr_args[i_bias])
    else
        generate_output()
    end
end