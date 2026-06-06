##################
#  dsl_string()

dsl_string(c::Char) = c
dsl_string(u::UInt8) = (u == CELL_CODE_INVALID) ? CELL_CHAR_INVALID : CELL_TYPES[u+1].char
dsl_string(s::CellTypeSet) = string(dsl_string.(s)...)

dsl_string(i::Int) = i
dsl_string(th::ThresholdByArea) = "(area*$(th.scale))"
dsl_string(th::ThresholdByLength) = "(length*$(th.scale))"
dsl_string(th::ThresholdRange) = "(($(dsl_string(th.a))):($(dsl_string(th.b))))"

dsl_string(ma::MarkovAlgorithm) = string(
    "@markovjunior '", dsl_string(ma.initial_fill), "' ",
    exists(ma.fixed_dimension) ? "$(ma.fixed_dimension) " : "",
    "begin
    ",
    iter_join(
        ("@pragma $name $(iter_join(args, " ")...)"
                for (name, args) in ma.pragmas),
        "\n    "
    )...,
    "
    ",
    iter_join(dsl_string.(ma.sequence), "\n    ")...,
    "
    end"
)


####################
#  DSL interface

mutable struct MacroParserInputs
    # Data given by the macro at the top level:
    initial_fill::Char
    fixed_dims::Optional{Int}

    # Data that operations can write into during parsing:
    min_dims::Int
    op_stack_trace::Stack{Any}
    bias_stack::Stack{Vector{AbstractMarkovBias}}
end

"
Parses a specific op (phrased in the DSL as a Julia macro expression),
  dispatching based on the op's name.
"
function parse_markovjunior_op(name_symbol_val,
                               macro_parser_inputs,
                               code_location,
                               expr_args,
                               original_expr
                              )::AbstractMarkovOp
    error("Unimplemented: ", typeof.((name_symbol_val, macro_parser_inputs, code_location, expr_args, original_expr)))
end

"
Parses a specific bias (phrased in the DSL as a function call),
  dispatching based on its name.
"
function parse_markovjunior_bias(name_symbol_val,
                                 macro_parser_inputs,
                                 code_location,
                                 expr_args
                                )::AbstractMarkovBias
    error("Unimplemented: ", typeof.((name_symbol_val, macro_parser_inputs, code_location, expr_args)))
end

"
Called once for every group of biases, for every type of bias in that group.
The intent is to allow new biases to add constraints on how they are used
  (e.g. throw error if more than one of themselves).

The group is understood to be an in-order accumulation of every group within `inputs.bias_stack`,
  where each sub-group is an inherited set of biases (e.g. nested sequences each with a bias section).
Note that you cannot modify the bias groups.
"
check_markovjunior_biases(type::Type, inputs::MacroParserInputs) = nothing


#####################
#  Main Logic

"Generates a markov algorithm using our DSL"
macro markovjunior(args...)
    return parse_markovjunior(args)
end

"Parses the arguments of a `@markovjunior` macro"
function parse_markovjunior(_macro_args)::MarkovAlgorithm
    macro_args = collect(Any, _macro_args) # Will delete them as they're processed

    # Decide on the initial fill value.
    initial_fill_char::Char = 'b'
    for (i, a) in enumerate(macro_args)
        if a isa Char
            initial_fill_char = a
            deleteat!(macro_args, i)
            break
        end
    end
    initial_fill = CELL_CODE_BY_CHAR[initial_fill_char]

    # Decide on the fixed-dimension.
    dims::Optional{Int} = nothing
    for (i, a) in enumerate(macro_args)
        if a isa Int
            @bp_check a > 0 "A Markov algorithm must be at least 1D; got $a"
            dims = a

            deleteat!(macro_args, i)
            break
        end
    end

    # Grab the main sequence.
    inputs = MacroParserInputs(initial_fill, dims, get_something(dims, 1), Stack{Any}(16), Stack{Vector{AbstractMarkovBias}}(16))
    push!(inputs.op_stack_trace, "Main Sequence")
    main_sequence = Vector{AbstractMarkovOp}()
    pragmas = Vector{Pair{Symbol, Vector{Any}}}()
    for (i, a) in enumerate(macro_args)
        if a isa Expr && a.head == :block
            main_sequence = parse_markovjunior_sequence(inputs, a.args) do location, line
                if @capture(line, @pragma prg__)
                    if (length(prg) < 1) || !isa(prg[1], Symbol)
                        raise_parse_error(line, inputs,
                                       "@pragma statement must have a name, e.g. @pragma Viz")
                    else
                        push!(pragmas, prg[1]=>collect(Any, prg[2:end]))
                    end
                    return true
                else
                    return false
                end
            end
            deleteat!(macro_args, i)
            break
        end
    end

    # Finish up.
    @bp_check isempty(macro_args) "Unexpected arguments: $macro_args"
    if exists(inputs.fixed_dims)
        if inputs.fixed_dims < inputs.min_dims
            error("The dimensionality is fixed at ", inputs.fixed_dims, "D, ",
                    "but the operations require at least ", inputs.min_dims, "D")
        else
            inputs.min_dims = inputs.fixed_dims
        end
    end
    return MarkovAlgorithm(initial_fill, inputs.fixed_dims, inputs.min_dims,
                           main_sequence, pragmas, Dict{Symbol, Any}())
end
"Tries to evaluate a `@markovjunior` macro, throwing an error if that's not what was parsed"
function parse_markovjunior(syntax::Union{String, Expr})::MarkovAlgorithm
    if syntax isa String
        return parse_markovjunior(Meta.parse(syntax))
    elseif !Base.isexpr(syntax, :macrocall) || (syntax.args[1] != Symbol("@markovjunior"))
        error("Expression wasn't a `@markovjunior` macro")
    else
        return parse_markovjunior(Tuple(syntax.args[3:end]))
    end
end

#


####################################
##   Parser Utilities

"Raises an error using the given LineNumberNode to point to user source"
function raise_parse_error(src::Optional{LineNumberNode}, state::MacroParserInputs, msg...)
    stringify_ast(x) = if x isa Union{Symbol, Expr}
        string(x)
    else
        x
    end
    error_args = Any[
        "ERROR at ",
        stringify_ast.(iter_join(
            reverse!(collect(state.op_stack_trace)),
            " \\\\ "
        ))...,
        ":\n\t",
        stringify_ast.(msg)...
    ]

    if isnothing(src)
        error(error_args...)
    else
        eval(Expr(:block, src, :( error($(error_args...)) )))
    end
end

function with_parser_stacktrace(to_do, inputs::MacroParserInputs, new_stack_entry)
    push!(inputs.op_stack_trace, new_stack_entry)
    try
        return to_do()
    finally
        pop!(inputs.op_stack_trace)
    end
end

"
Processes each of the given lines, tracking the most recent LineNumberNode.

First it lets your lambda handle the line, if you gave one.
You should return whether or not you handled that line.

If not, then this function looks for a macro call and assumes it's an op.

Any other unhandled lines turn into an error.
"
function parse_markovjunior_sequence(try_handle_line, inputs::MacroParserInputs, block_args)::Vector{AbstractMarkovOp}
    output = Vector{AbstractMarkovOp}()
    i = Ref(0)
    parse_markovjunior_block(block_args) do location, line
        i[] += 1
        if try_handle_line(location, line)
            # Do nothing; the line was handled.
        elseif (line isa Expr) && (line.head == :macrocall)
            push!(inputs.op_stack_trace, "Op $(i[]) `$(line.args[1])`")
            push!(output, parse_markovjunior_op(
                Val(line.args[1]::Symbol), inputs,
                line.args[2]::LineNumberNode,
                line.args[3:end], line
            ))
            pop!(inputs.op_stack_trace)
        else
            raise_parse_error(location, inputs,
                           "Unexpected sequence expression: '", line, "'")
        end
    end
    return output
end
function parse_markovjunior_sequence(inputs::MacroParserInputs, block_args)::Vector{AbstractMarkovOp}
    return parse_markovjunior_sequence((loc, line) -> false, inputs, block_args)
end

"
Processes a bias statement/block-of-statements,
  pushing them onto the end of `inputs.bias_stack` and validating the result
  before returning it.

Make sure to pop this off the stack once you're done parsing your op!

It's recommended to use `with_parsed_markovjunior_bias_statement` instead.
"
function push_parsed_markovjunior_bias_statement(inputs::MacroParserInputs, location, biases)
    # Define how to process each statement.
    output = Vector{AbstractMarkovBias}()
    function process_line(location, line)
        push!(inputs.op_stack_trace, "Bias \"$line\"")
        try
            if @capture line f_Symbol(args__)
                push!(output, parse_markovjunior_bias(Val(f), inputs, location, args))
            else
                raise_parse_error(location, inputs,
                               "Invalid bias syntax! Expected a function call, got:", line)
            end
        finally
            pop!(inputs.op_stack_trace)
        end
    end

    # Process each statement.
    if @capture biases f_Symbol(args__)
        process_line(location, biases)
    elseif Base.isexpr(biases, :block)
        parse_markovjunior_block(process_line, biases.args)
    else
        raise_parse_error(location, inputs,
                       "Unexpected expression where bias was expected: ", biases)
    end

    # Validate the result.
    if !isempty(output)
        push!(inputs.bias_stack, output)
        for T in unique(typeof.(Iterators.flatten(inputs.bias_stack)))
            check_markovjunior_biases(T, inputs)
        end
    end

    return output
end

"
Parses the given group of biases, pushes them onto the stack,
  passes the new group into your lambda, then pops them off the stack again.

You should use this when parsing ops that define a `bias` section.
"
function with_parsed_markovjunior_bias_statement(to_do, inputs::MacroParserInputs,
                                                 location, expr)
    biases = push_parsed_markovjunior_bias_statement(inputs, location, expr)
    try
        return to_do(biases)
    finally
        if !isempty(biases)
            pop!(inputs.bias_stack)
        end
    end
end

"
Checks whether an expression looks to be a Threshold value.

This doesn't guarantee the threshold is well-formed!
It's meant to disambiguate statements containing optional thresholds and other optional things.
You need to make sure those other things can't look like thresholds.
"
check_markovjunior_threshold_appearance(expr)::Bool =
    #NOTE: The '|' operator is broken in @capture sadly.
    (expr isa Integer) || @capture(expr, a_:b_) ||
    @capture(expr, length/x_) || @capture(expr, length*x_) || @capture(expr, x_*length) ||
    @capture(expr, area/x_) || @capture(expr, area*x_) || @capture(expr, x_*area)
parse_markovjunior_threshold(inputs::MacroParserInputs, location, threshold_expr)::Threshold = parse_markovjunior_threshold(
    expr -> nothing,
    inputs, location, threshold_expr
)
function parse_markovjunior_threshold(try_handle,
                                      inputs::MacroParserInputs, location, threshold_expr)
    push!(inputs.op_stack_trace, "Threshold statement")
    try # Ensure stack trace is popped at end
        user_attempt = try_handle(threshold_expr)
        if exists(user_attempt)
            user_attempt
        elseif threshold_expr isa Integer
            return convert(Int, threshold_expr)
        elseif @capture threshold_expr (area/x_Real)
            return ThresholdByArea(convert(Float32, 1/x))
        elseif @capture(threshold_expr, (area*x_Real)) || @capture(threshold_expr, (x_Real*area))
            return ThresholdByArea(convert(Float32, x))
        elseif @capture threshold_expr (length/x_Real)
            return ThresholdByLength(convert(Float32, 1/x))
        elseif @capture(threshold_expr, (length*x_Real)) || @capture(threshold_expr, (x_Real*length))
            return ThresholdByLength(convert(Float32, x))
        elseif @capture threshold_expr (a_:b_)
            aa = begin
                push!(inputs.op_stack_trace, "Range start `$a`")
                result = parse_markovjunior_threshold(inputs, location, a)
                if !isa(result, ThresholdScalar)
                    raise_parse_error(location, inputs, "Value not a scalar")
                end
                pop!(inputs.op_stack_trace)
                result
            end
            bb = begin
                push!(inputs.op_stack_trace, "Range end `$b`")
                result = parse_markovjunior_threshold(inputs, location, b)
                if !isa(result, ThresholdScalar)
                    raise_parse_error(location, inputs, "Value not a scalar")
                end
                pop!(inputs.op_stack_trace)
                result
            end
            return ThresholdRange(aa, bb)
        else
            raise_parse_error(location, inputs, "Unexpected format for threshold; expected")
        end
    # Handle stack trace no matter what.
    finally
        pop!(inputs.op_stack_trace)
    end
end

"
Runs your lambda on every top-level expression that isn't a LineNumberNode,
  and includes the last LineNumberNode in each call.

Before the first LineNumberNode is hit, `nothing` is passed in place of it.
"
function parse_markovjunior_block(to_do, block_lines)
    last_src_line::Optional{LineNumberNode} = nothing
    for block_line in block_lines
        if block_line isa LineNumberNode
            last_src_line = block_line
        else
            to_do(last_src_line, block_line)
        end
    end
end

"Grabs, returns (or `nothing` if not found), and deletes the assignment expession with the given lhs"
function peel_markovjunior_block_assignment(inout_block_args, name::Symbol)::Optional
    for (i, a) in enumerate(inout_block_args)
        if Base.isexpr(a, :(=)) && a.args[1] == name
            value = a.args[2]
            deleteat!(inout_block_args, i)
            return value
        end
    end
    return nothing
end