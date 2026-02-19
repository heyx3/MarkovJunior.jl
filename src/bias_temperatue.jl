
struct MarkovBiasTemperature <: AbstractMarkovBias
    amount::Float32
end

markov_bias_calculate(t::MarkovBiasTemperature, ::Nothing,
                      ::CellGrid{N}, ::Any, ::CellLine{N},
                      rng::PRNG) where {N} = rand(rng, Float32) * t.amount

dsl_string(b::MarkovBiasTemperature) = "temperature($(b.amount))"

function parse_markovjunior_bias(::Val{:temperature}, inputs::MacroParserInputs,
                                 loc::LineNumberNode, args
                                )::MarkovBiasTemperature
    push!(inputs.op_stack_trace, "temperature(...)")
    if length(args) != 1
        raise_error_at(loc, inputs, "Expected one parameter, got ", length(args))
    elseif !isa(args[1], Real)
        raise_error_at(loc, inputs, "Expected a number, got ", type(args[1]))
    else
        pop!(inputs.op_stack_trace)
        return MarkovBiasTemperature(args[1])
    end
end

function check_markovjunior_biases(::Type{<:MarkovBiasTemperature},
                                   inputs::MacroParserInputs)
    # Within a single group, there should only be one temperature parameter.
    # Across inherited groups, it's fine to have multiple.
    if count(t -> t isa MarkovBiasTemperature, top(inputs.bias_stack)) > 1
        raise_error_at(nothing, inputs, "Found more than one temperature() bias!")
    end
end