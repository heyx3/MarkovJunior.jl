# Thresholds are an expression indicating how many times some Op should run.

"
What a threshold calculation can reference:
  * `size`: the size of the current grid, as an array,
      so `size_prod()` is its number of cells and `size_average()` its average side length.
  * `area`: how many cells are being worked on. Usually that's `size_prod()`,
      but sometimes (e.g. `@within_regions`) we are working on a non-rectangular subset of the grid.
  * `ndims`: the number of dimensions, i.e. `size_len()`.
  * The algorithm's user functions (from `@calc`), looked up whenever they're called.
"
const THRESHOLD_CALC_ENV = CalcEnvironment(
    inputs = (:area, :ndims),
    arrays = (:size, ),
    macros = true
)

"
How many times something should happen (e.g. how many rewrites, or how many sequence repetitions).
It's a calculation (see `CalcInterpreter`), evaluated whenever the Op starts,
  then rounded to the nearest `Int`.
A threshold of 0 or less means the Op does nothing.
"
struct Threshold
    calc::CalcProgram
end
Base.:(==)(a::Threshold, b::Threshold) = (a.calc == b.calc)
Base.hash(t::Threshold, h::UInt) = hash(t.calc, hash(Threshold, h))
Base.show(io::IO, t::Threshold) = print(io, "Threshold(`", t.calc, "`)")

# Thresholds are always written wrapped in `threshold(...)`, so their position never needs guessing.
dsl_format(th::Threshold) = string("threshold(", th.calc, ")")

"Logs a problem with a threshold during a run, without stopping the algorithm"
threshold_log_error(msg...) = println(stderr, "Threshold error: ", msg...)

"
Evaluates the threshold for the given grid, where `area` is how many of its cells are being worked on
  (see `AlgoState.grid_area`).
Errors are logged, and the threshold falls back to 1.
"
function get_threshold(th::Threshold, grid::CellGrid, area::Int, rng::PRNG,
                       calc_macros::CalcMacros)::Int
    raw = calc_evaluate(
        th.calc, rng,
        e -> begin
            threshold_log_error(e.msg, "; using 1")
            1
        end,
        calc_macros
        ;
        size = size(grid),
        area = area,
        ndims = ndims(grid)
    )
    if (raw isa Bool) || !isfinite(raw)
        threshold_log_error("`", th.calc, "` should give a finite number; got `", raw, "`. Using 1")
        return 1
    end
    # Clamp before converting, so huge values don't overflow.
    return round(Int, clamp(raw, -1, 1e15))
end


"
If the expression is wrapped in `threshold(...)`, returns what's inside it; otherwise returns `nothing`.
That syntax is used to disambiguate thresholds from other arguments when necessary.
"
unwrap_markovjunior_threshold(expr) = if Base.isexpr(expr, :call, 2) && (expr.args[1] == :threshold)
    expr.args[2]
else
    nothing
end

"
Compiles a threshold (optionally wrapped in `threshold(...)`), turning calculation errors into parse errors.
Like any calculation, user functions it calls aren't checked until it runs.
"
function parse_markovjunior_threshold(inputs::MacroParserInputs, location, threshold_expr)::Threshold
    inner = unwrap_markovjunior_threshold(threshold_expr)
    calc_expr = isnothing(inner) ? threshold_expr : inner
    with_parser_stacktrace(inputs, "Threshold") do
        Threshold(parse_markovjunior_calc(
            inputs, location, calc_expr,
            THRESHOLD_CALC_ENV, "Invalid syntax"
        ))
    end
end
