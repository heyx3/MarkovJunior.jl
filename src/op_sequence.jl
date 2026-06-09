const SequenceRepeatModeTag = Val{:repeat}

struct MarkovOpSequence <: AbstractMarkovOp
    ops::Vector{AbstractMarkovOp}
    threshold::Union{Nothing, SequenceRepeatModeTag, Threshold}
    biases::Vector{AbstractMarkovBias}
end


struct MarkovOpSequence_State
    repetitions_left::Union{SequenceRepeatModeTag, Int}
    current_op_idx::Int
    current_op_state::Any #TODO: Take advantage of markov_op_state_type() to make this type-stable
end

markov_op_state_type(::Type{MarkovOpSequence}, ::Val{NDims}) where {NDims} = MarkovOpSequence_State
function markov_op_initialize(s::MarkovOpSequence, expected_state_type::Type{MarkovOpSequence_State},
                              grid_ref::Ref{<:CellGrid{N}},
                              rng::PRNG, context::MarkovOpContext) where {N}
    grid = grid_ref[]

    # Set up the repetition counter.
    repetitions_left = if isnothing(s.threshold)
        0
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
        return nothing
    end

    # Append our biases to the stack.
    append!(context.all_biases, s.biases)
    n_biases = length(context.all_biases)
    n_self_biases = length(s.biases)
    clear_biases() = deleteat!(context.all_biases, (n_biases - n_self_biases + 1):n_biases)

    # Start up the first op.
    if isempty(s.ops)
        clear_biases()
        return nothing
    end
    op_idx = 1
    op_initial_state = markov_op_initialize(
        s.ops[1],
        markov_op_state_type(s.ops[1], typeof(grid), rng, context),
        grid, rng, context
    )
    # Check if the first op failed to start (e.g. because it had no matches).
    # As long as we're not in repeat mode, we should keep moving forward until an op succeeds.
    if !isa(repetitions_left, SequenceRepeatModeTag)
        while isnothing(op_initial_state) && (op_idx < length(s.ops))
            op_idx += 1
            op_initial_state = markov_op_initialize(
                s.ops[op_idx],
                markov_op_state_type(s.ops[op_idx], typeof(grid), rng, context),
                grid, rng, context
            )
        end
    end
    if isnothing(op_initial_state)
        clear_biases()
        return nothing
    end

    return MarkovOpSequence_State(repetitions_left, op_idx, op_initial_state)
end

function markov_op_cancel(s::MarkovOpSequence, state::MarkovOpSequence_State,
                          context::MarkovOpContext)
    markov_op_cancel(s.ops[state.current_op_idx], state.current_op_state, context)

    n_biases = length(context.all_biases)
    n_self_biases = length(s.biases)
    deleteat!(context.all_biases, (n_biases - n_self_biases + 1):n_biases)
end

function markov_op_iterate(s::MarkovOpSequence, state::MarkovOpSequence_State,
                           grid::Ref{<:CellGrid{N}}, rng::PRNG,
                           context::MarkovOpContext,
                           n_ticks_left::Ref{Optional{Int}}) where {N}
    # Run each sub-op within a type-stable context.
    function run_sub_op(sub_op::TOp, current_state::TState, idx::Int,
                        ::Val{BInfiniteTicks}
                       ) where {TOp, TState, BInfiniteTicks}
        # Note that we may need to handle cases where the initial state is `nothing`.
        @logic_logln("Running sub-op ", idx, BInfiniteTicks ? " infinitely" : "",
                     ": ", typeof(sub_op))
        next_state = current_state
        while exists(next_state) && (BInfiniteTicks || (n_ticks_left[] > 0))
            next_state = markov_op_iterate(sub_op, next_state, grid, rng, context, n_ticks_left)
        end
        if exists(next_state)
            @logic_logln("    Couldn't finish the sub-op")
        else
            @logic_logln("    Finished the sub-op!")
        end
        return next_state
    end
    b_infinite_ticks = Val(isnothing(n_ticks_left[]))

    @markovjunior_assert exists(state.current_op_state) "Iterating on a non-running Sequence!?"

    # Nest the sequence logic in a loop to handle repetitions and insta-quitting ops.
    while get_something(n_ticks_left[], 1) > 0
        # Advance the current sub-op.
        if exists(state.current_op_state)
            @logic_logln("Running sub-op ", state.current_op_idx, "...")
            @logic_tab_in()
            @set! state.current_op_state = run_sub_op(
                s.ops[state.current_op_idx], state.current_op_state,
                state.current_op_idx, b_infinite_ticks
            )
            @logic_tab_out()
        end

        # If it ended, advance to the next one (until we find one that can actually run).
        while isnothing(state.current_op_state) && (state.current_op_idx < length(s.ops))
            @logic_logln("Advancing to the next sub-op...")
            @logic_tab_in()
            state = MarkovOpSequence_State(
                state.repetitions_left,
                state.current_op_idx + 1,
                markov_op_initialize(
                    s.ops[state.current_op_idx + 1],
                    markov_op_state_type(s.ops[state.current_op_idx + 1],
                                         typeof(grid[]),
                                         rng, context),
                    grid, rng, context
                )
            )
            @logic_tab_out()
            if isnothing(state.current_op_state)
                @logic_logln("...and that sub-op failed to do anything.")
            end
        end

        # If we advanced through all ops, repeat or quit.
        if isnothing(state.current_op_state)
            @logic_logln("Finished one repetition of this sequence!")
            if state.repetitions_left == 0
                @logic_logln("   Now this sequence is complete.")
                markov_op_cancel(s, state, context)
                return nothing
            else
                @logic_logln("Starting first op again with ",
                            if state.repetitions_left isa Int
                                state.repetitions_left
                            elseif state.repetitions_left isa SequenceRepeatModeTag
                                "indefinite"
                            else
                                error("Unhandled: ", typeof(state.repetitions_left))
                            end,
                            " repetitions remaining")
                @logic_tab_in()
                state = MarkovOpSequence_State(
                    if state.repetitions_left isa Int
                        state.repetitions_left - 1
                    elseif state.repetitions_left isa SequenceRepeatModeTag
                        state.repetitions_left
                    else
                        error("Unhandled: ", typeof(state.repetitions_left))
                    end,
                    1,
                    markov_op_initialize(
                        s.ops[1],
                        markov_op_state_type(s.ops[1], typeof(grid[]), rng, context),
                        grid, rng, context
                    )
                )
                @logic_tab_out()

                # If the first op failed to match and we're using the 'repeat' threshold, we're done.
                if (state.repetitions_left isa SequenceRepeatModeTag) && isnothing(state.current_op_state)
                    @logic_logln("The first op failed to start and we're in indefinite 'repeat' mode, so we're done!")
                    markov_op_cancel(s, state, context)
                    return nothing
                end
            end
        end
    end

    @markovjunior_assert(
        exists(state.current_op_state),
        "Finished an iteration of @sequence without ending, but also without a running Op!"
    )
    return state
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
        biases
    )
    return if exists(i_bias)
        with_parsed_markovjunior_bias_statement(generate_output, inputs, loc, expr_args[i_bias])
    else
        generate_output()
    end
end