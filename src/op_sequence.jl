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
                         ::Val{NGrid} = Val(ndims(algo_state.grid))
                        )::Tuple{Bool, typeof(inherited_bias_states)} where {NGrid, NInheritedBiases, NSelfBiases}
    # Set up the repetition counter.
    repetitions_left = if isnothing(s.threshold)
        typemax(Int)
    elseif s.threshold isa SequenceRepeatModeTag
        s.threshold
    elseif s.threshold isa Threshold
        get_threshold(s.threshold, ThresholdInputs(
            convert(Float32, prod(size(grid))),
            convert(Float32, sum(size(grid), init=0) / N),
            rng
        )) - 1
    else
        error("Unhandled: ", typeof(s.threshold))
    end
    if (repetitions_left isa Int) && (repetitions_left < 0)
        return false
    end

    # Initialize our own biases.
    all_biases = tuple(
        inherited_biases...,
        biases...
    )
    all_bias_states = tuple(
        inherited_bias_states...,
        markov_bias_initialize.(biases, Ref(algo), Ref(algo_state))...
    )

    # Run the loop.
    @logic_logln("Starting sequence with threshold `", repetitions_left, "`")
    made_any_changes::Bool = false
    rep_i::Int = 0
    while (repetitions_left isa SequenceRepeatModeTag) || (repetitions_left > 0)
        if repetitions_left isa Int
            repetitions_left -= 1
        end
        rep_i += 1
        @logic_logln("REPETITION ", rep_i)

        for op_i in 1:length(sequence.ops)
            @logic_logln("Inner op ", op_i, ":")
            @logic_tab_in()

            (op_made_changes, all_bias_states) = markov_algo_run(
                sequence.ops[op_i], algo, algo_state,
                all_biases, all_bias_states
            )
            made_any_changes |= op_made_changes
            if !op_made_changes && (op_i == 1) && (repetitions_left isa SequenceRepeatModeTag)
                @logic_logln("Inner op did nothing! The outer sequence will end now")
                repetitions_left = 0
            end

            @logic_tab_out()
            (repetitions_left == 0) && break
        end
    end

    foreach(ntuple(identity, Val(NSelfBiases))) do i
        markov_bias_cleanup(all_biases[i + NInheritedBiases],
                            all_bias_states[i + NInheritedBiases],
                            algo, algo_state)
    end

    markov_algo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY + 1)
    return (made_any_changes, all_bias_states[1:NInheritedBiases])
end

dsl_string(s::MarkovOpSequence) = string(
    "@sequence ",
    if isnothing(s.threshold)
        ""
    elseif s.threshold isa SequenceRepeatModeTag
        "repeat"
    elseif s.threshold isa Threshold
        dsl_string(s.threshold)
    else
        error("Unhandled: ", typeof(s.threshold))
    end,
     " ",
    "begin\n    ",
    iter_join(
        Iterators.map(dsl_string, s.ops),
        "\n    "
    )...,
    "\nend begin\n    ",
    iter_join(
        Iterators.map(dsl_string, s.biases),
        "\n    "
    )...,
    "\nend"
)

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
        return AbstractMarkovOp[ process_op(expr, nothing) ]
    elseif Base.isexpr(expr, :block)
        line_idx = Ref(1)
        output = AbstractMarkovOp[ ]
        parse_markovjunior_block(expr.args) do inner_src, line
            push!(output, process_op(line, line_idx[]))
            line_idx[] += 1
        end
        return output
    else
        raise_parse_error(src, inputs,
                       "Expected a block of operations; got `",
                       if expr isa Expr
                           tuple("`", expr.head, "` expression")
                       else
                           tuple(typeof(expr))
                       end...)
    end
end

function parse_markovjunior_op(::Val{Symbol("@sequence")}, inputs::MacroParserInputs,
                               loc, expr_args, expr_line)
    (i_threshold, i_ops, i_bias) = if length(expr_args) == 3
        (1, 2, 3)
    elseif length(expr_args) == 2
        if check_markovjunior_threshold_appearance(expr_args[1]) || (expr_args[1] == :repeat)
            ( 1, 2, nothing)
        else
            (nothing, 1, 2)
        end
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