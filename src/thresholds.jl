struct ThresholdByArea
    scale::Float32
end
struct ThresholdByLength
    scale::Float32
end
const ThresholdScalar = Union{Int, ThresholdByArea, ThresholdByLength}

"A randomly chosen threshold between two (unordered) values"
struct ThresholdRange
    a::ThresholdScalar
    b::ThresholdScalar
end

const Threshold = Union{ThresholdScalar, ThresholdRange}


struct ThresholdInputs
    area::Float32
    length::Float32
    rng::PRNG
end

get_threshold(th::Threshold, in::ThresholdInputs)::Int = max(1, round(Int, get_raw_threshold(th, in)))
get_threshold(th::Threshold, grid::CellGrid, rng::PRNG)::Int = get_threshold(th, ThresholdInputs(
    convert(Float32, prod(vsize(grid))),
    convert(Float32, sum(vsize(grid)) / ndims(grid)),
    rng
))

get_raw_threshold(th::Int, in::ThresholdInputs) = convert(Float32, th)
get_raw_threshold(th::ThresholdByArea, in::ThresholdInputs) = convert(Float32, in.area) * th.scale
get_raw_threshold(th::ThresholdByLength, in::ThresholdInputs) = convert(Float32, in.length) * th.scale
function get_raw_threshold(th::ThresholdRange, in::ThresholdInputs)
    a::Float32 = get_raw_threshold(th.a, in)
    b::Float32 = get_raw_threshold(th.b, in)
    return lerp(a, b, rand(in.rng, Float32))
end


dsl_string(th::ThresholdByArea) = "(area*$(th.scale))"
dsl_string(th::ThresholdByLength) = "(length*$(th.scale))"
dsl_string(th::ThresholdRange) = "(($(dsl_string(th.a))):($(dsl_string(th.b))))"

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
