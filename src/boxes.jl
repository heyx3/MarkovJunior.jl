@bp_enum(DrawBoxSpace,
    uv,
    pixel
)

"
A box within the grid, used by ops like `@fill` and `@within_box`.
It's written like `uv(min=0.25, size=0.5)`, or `pixel(center=center, size=[3, 4])`.

The choice of uv or pixel is the 'space',
  and the parameters must be two of `min`, `max`, `size`, and `center`.
UV space ranges from 0-1, while pixel space starts at 1.

Each param value is a calculation (via `CalcInterpreter`), either a list with one per axis
  or a single expression for every axis.
The single-expression form is evaluated once and copied across every axis,
  EXCEPT when it's wrapped in `each(...)` or uses one of the per-axis inputs like `len`.

If the box covers fewer dimensions than the grid, it stretches all the way across the extra axes
  (e.g. a 2D box implicitly covers the entire Z axis).
"
struct CalcBox
    space::E_DrawBoxSpace
    # Two of `CALC_BOX_ARG_NAMES`, always in the same order as they appear there.
    sorted_param_names::NTuple{2, Symbol}
    # For each parameter, the per-axis or singleton expression.
    param_values::NTuple{2, Vector{CalcProgram}}
    param_is_list::NTuple{2, Bool}
    # For single values, whether they were written as `each(...)`, to be evaluated separately per axis.
    param_each::NTuple{2, Bool}
end
Base.:(==)(a::CalcBox, b::CalcBox) = all(f -> getfield(a, f) == getfield(b, f), fieldnames(CalcBox))
Base.hash(b::CalcBox, h::UInt) = hash(map(f -> getfield(b, f), fieldnames(CalcBox)), hash(CalcBox, h))

"
The parameters that are needed to define a box (take two).
NOTE: other code in this file assumes this order!
"
const CALC_BOX_ARG_NAMES = (:min, :max, :size, :center)

"
The inputs every box calculation gets, describing the axis it's being evaluated for:
  * `len`: the area's length along the axis, in pixels.
  * `center`: the middle of the axis, in the box's space (0.5 for UV, `len/2` for pixels).
  * `texel`: the size of one pixel in UV space (`1/len`).

Use of any of these in a scalar parameter causes it to get evaluated once per-axis,
  rather than evaluated once and *copied* across axes.
"
const CALC_BOX_AXIS_INPUTS = (:len, :center, :texel)

"
What `@fill` and `@within_box` box calculations can reference:
  * The axis inputs (see `CALC_BOX_AXIS_INPUTS`).
  * `size`: the resolution of the current grid -- an array of integers.
  * `ndims`: the number of dimensions, i.e. `size_len()`.
  * The algorithm's user functions (from `@calc`).
"
const CALC_BOX_ENV = CalcEnvironment(
    inputs = (:ndims, CALC_BOX_AXIS_INPUTS...),
    arrays = (:size, ),
    macros = true
)

"Logs a problem with box calculation during a run, without stopping the algorithm"
box_log_error(msg...) = println(stderr, "Box error: ", msg...)

"
Whether the box's `i`th parameter should be evaluated separately for each axis.
Returns true under the following conditions:

* The expression is actually a per-axis list of expressions
* The expression was wrapped in `each(...)`.
* The expression mentions a per-axis input, implicitly opting into per-axis evaluation.
"
calc_box_param_is_per_axis(b::CalcBox, i::Int, macros::CalcMacros)::Bool =
    b.param_is_list[i] || b.param_each[i] ||
    calc_mentions(b.param_values[i][1], CALC_BOX_AXIS_INPUTS, macros)


"
Calculates the pixel area covered by the given box within an area of the given size.
If the box fails to cover any pixels, `nothing` is returned.

The box's edges are rounded to the nearest pixel edge (rounding halfway up), then clamped to the area.

`evaluate(calc::CalcProgram, axis_inputs::NamedTuple) -> CalcValue`
  crunches one of the numbers defining the box.
This gives you the opportunity to inject your own inputs to the calculator's context.
`axis_inputs` holds the inputs for the axis being evaluated (see `CALC_BOX_AXIS_INPUTS`)
  IF evaluating for a single axis; otherwise it's empty.
The caller should put them into the calc context too.

Any time a calculation fails to produce a finite number,
  0 is used in its place and `log_error(msg...)` is called.

`macros` are user functions the calculations may call.
"
function eval_calc_box(evaluate, log_error, b::CalcBox, area_size::Vec{N, <:Integer},
                       macros::CalcMacros)::Optional{BoxI{N}} where {N}
    n_axes = if b.param_is_list[1]
        length(b.param_values[1])
    elseif b.param_is_list[2]
        length(b.param_values[2])
    else
        N
    end
    is_uv = (b.space == DrawBoxSpace.uv)

    # Computes the given parameter expression and interprets it as a finite number.
    function to_number(calc::CalcProgram, axis_inputs::NamedTuple)::Float64
        raw = evaluate(calc, axis_inputs)
        if !isa(raw, Real) || isa(raw, Bool) || !isfinite(raw)
            log_error("`", calc, "` should give a finite number; got `", raw, "`. Using 0")
            0.0
        else
            convert(Float64, raw)
        end
    end

    # Values shared by every axis are evaluated once, up front.
    param_shared_values::NTuple{2, Float64} = ntuple(Val(2)) do i
        if calc_box_param_is_per_axis(b, i, macros)
            NaN
        else
            to_number(b.param_values[i][1], NamedTuple())
        end
    end

    # Get each parameter's value along each axis,
    #   as continuous coordinates where the area goes from 0 to `len`.
    # *Not* `len-1` because it goes out to the max edge.
    # Extra grid axes are handled afterwards.
    axis_len(axis::Int)::Int = convert(Int, area_size[axis])
    param_coords(i::Int) = Vec{N, Float64}() do axis
        if axis > n_axes
            return 0.0
        end

        a_len = axis_len(axis)
        v = if !isnan(param_shared_values[i])
            param_shared_values[i]
        else
            program::CalcProgram = b.param_values[i][b.param_is_list[i] ? axis : 1]
            to_number(program, (
                len = a_len,
                center = is_uv ? 0.5 : (a_len / 2),
                texel = 1 / a_len
            ))
        end

        name = b.sorted_param_names[i]
        return if is_uv
            v * a_len
        elseif name == :min
            v - 1 # The low edge of the pixel
        elseif name == :center
            v - 0.5 # The middle of the pixel
        else # `max` is the high edge of the pixel, and sizes need no change
            @markovjunior_assert(name in (:max, :size))
            v
        end
    end
    continuous_box = Box{N, Float64}(NamedTuple{b.sorted_param_names}((param_coords(1), param_coords(2))))
    continuous_min = min_inclusive(continuous_box)
    continuous_max = max_exclusive(continuous_box)

    # Round to pixel edges.
    edges_per_side_per_axis::NTuple{N, NTuple{2, Int}} = ntuple(Val(N)) do axis
        a_len = axis_len(axis)
        if axis > n_axes
            return (0, a_len)
        end

        to_pixel_edge(x::Float64) = floor(Int, clamp(x, 0, a_len) + 0.5)
        return to_pixel_edge.((continuous_min[axis], continuous_max[axis]))
    end
    if any(((start, stop), ) -> (stop <= start), edges_per_side_per_axis)
        return nothing
    else
        return BoxI{N}(
            min = Vec{N, Int32}(a -> convert(Int32, edges_per_side_per_axis[a][1] + 1)),
            max = Vec{N, Int32}(a -> convert(Int32, edges_per_side_per_axis[a][2]))
        )
    end
end

"Calculates a box for the given grid, returning the pixels it covers or `nothing` if it's empty"
function eval_calc_box(b::CalcBox, grid::CellGrid{N}, rng::PRNG, calc_macros::CalcMacros)::Optional{BoxI{N}} where {N}
    grid_size = size(grid)
    return eval_calc_box(box_log_error, b, Vec{N, Int}(grid_size...), calc_macros) do calc, axis_inputs
        calc_evaluate(
            calc, rng,
            e -> begin
                box_log_error(e.msg, "; using 0")
                0
            end,
            calc_macros
            ;
            size = grid_size,
            ndims = N,
            axis_inputs...
        )
    end
end


function dsl_format_calc_box_value(b::CalcBox, i::Int)::String
    values = b.param_values[i]
    if b.param_is_list[i]
        return string("[", iter_join(Iterators.map(v -> "($v)", values), ", ")..., "]")
    elseif b.param_each[i]
        return "each(($(values[1])))"
    else
        return "($(values[1]))"
    end
end
"Formats a box like `uv(min=0.25, size=[0.5, 1])`, always writing per-axis lists in brackets"
dsl_format(b::CalcBox) = string(
    b.space, "(",
    b.sorted_param_names[1], "=", dsl_format_calc_box_value(b, 1), ", ",
    b.sorted_param_names[2], "=", dsl_format_calc_box_value(b, 2),
    ")"
)

"
Parses a box expression such as `uv(center=0.5, size=[0.3, 1])`.
Per-axis lists may be written with brackets or parentheses,
  and a single value may be wrapped in `each(...)` to evaluate it separately for each axis.

All per-axis lists must have the same length within one box,
  and the length also puts a floor on the algorithm's dimensionality.
"
function parse_markovjunior_calc_box(inputs::MacroParserInputs, loc, expr, env::CalcEnvironment)::CalcBox
    box_error() = raise_parse_error(
        loc, inputs,
        "Expected a box like `uv(min=0, size=0.5)` or ",
          "`pixel(center=center, max=[3, 4])`; got `", expr, "`"
    )
    if !Base.isexpr(expr, :call) || !isa(expr.args[1], Symbol)
        box_error()
    end

    # Extract the coordinate space.
    space = try
        DrawBoxSpace.from(expr.args[1])
    catch
        raise_parse_error(loc, inputs, "A box's space should be one of ",
                            iter_join(Iterators.map(s -> "`$s`", DrawBoxSpace.instances()), ", ")...,
                            "; got `", expr.args[1], "`")
    end

    # Extract the parameters.
    named_values = Dict{Symbol, Any}()
    add_value(kw) = begin
        if !Base.isexpr(kw, :kw) && !Base.isexpr(kw, :(=))
            box_error()
        end

        (name, value) = kw.args
        if !in(name, CALC_BOX_ARG_NAMES)
            raise_parse_error(
                loc, inputs,
                "A box's values are two of ",
                  iter_join(Iterators.map(n -> "`$n`", CALC_BOX_ARG_NAMES), ", ")...,
                  "; got `", name, "`"
            )
        elseif haskey(named_values, name)
            raise_parse_error(loc, inputs, "A box's `", name, "` is given more than once")
        end

        named_values[name] = value
    end
    for arg in expr.args[2:end]
        if Base.isexpr(arg, :parameters)
            foreach(add_value, arg.args)
        else
            add_value(arg)
        end
    end
    if length(named_values) != 2
        raise_parse_error(loc, inputs, "A box needs exactly two of ",
                            iter_join(Iterators.map(n -> "`$n`", CALC_BOX_ARG_NAMES), ", ")...,
                            "; got `", expr, "`")
    end
    sorted_names = Tuple(filter(n -> haskey(named_values, n), CALC_BOX_ARG_NAMES))

    # Compile the parameter values into calculator expressions.
    function compile_value(name::Symbol)::Tuple{Vector{CalcProgram}, Bool, Bool}
        val = named_values[name]
        is_list = Base.isexpr(val, :vect) || Base.isexpr(val, :tuple)
        is_each = Base.isexpr(val, :call) && (val.args[1] == :each)
        if is_each && (length(val.args) != 2)
            raise_parse_error(loc, inputs, "Expected `each(value)` in the box's `", name, "`; got `", val, "`")
        end

        param_exprs = if is_list
            val.args
        elseif is_each
            Any[ val.args[2] ]
        else
            Any[ val ]
        end
        if isempty(param_exprs)
            raise_parse_error(loc, inputs, "A box's `", name, "` list can't be empty")
        end
        param_calcs = CalcProgram[
            parse_markovjunior_calc(inputs, loc, pe, env, "Box `$name`")
              for pe in param_exprs
        ]

        if is_list
            inputs.min_dims = max(inputs.min_dims, length(param_calcs))
        end
        return (param_calcs, is_list, is_each)
    end

    ((param1, param1_is_list, param1_each), (param2, param2_is_list, param2_each)) = map(compile_value, sorted_names)
    if param1_is_list && param2_is_list && (length(param1) != length(param2))
        raise_parse_error(
            loc, inputs,
            "The box's `", sorted_names[1], "` and `", sorted_names[2],
              "` parameters are both per-axis yet different lengths!"
        )
    end

    return CalcBox(space, sorted_names, (param1, param2),
                   (param1_is_list, param2_is_list), (param1_each, param2_each))
end
