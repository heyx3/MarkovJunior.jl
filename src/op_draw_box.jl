#NOTE: this file takes some things from op_rewrite.jl and should be included after it

#############
#  Data

@bp_enum(DrawBoxSpace,
    uv,
    pixel
)

"Gets the clamped grid area covered by the given draw-box settings"
function get_draw_box_pixels(space::E_DrawBoxSpace,
                             box::Box{N2, R},
                             resolution::VecI{N},
                             box_is_1D_scalar::Bool
                            )::Box{N, Int32} where {N, N2, R<:Union{Float32, Int32}}
    # The box size becomes infinity on any unmentioned axis.
    if N2 != N
        mi = min_inclusive(box)
        sz = size(box)
        return get_draw_box_pixels(
            space,
            Box(
                min=Vec{N, R}() do i
                    if i > N2
                        box_is_1D_scalar ? mi[1] : 0
                    else
                        mi[i]
                    end
                end,
                size=Vec{N, R}() do i
                    if i > N2
                        if box_is_1D_scalar
                            sz[1]
                        elseif space == DrawBoxSpace.uv
                            1
                        elseif space == DrawBoxSpace.pixel
                            resolution[i] + 1
                        else
                            error("Unhandled: ", space)
                        end
                    else
                        sz[i]
                    end
                end
            ),
            resolution,
            false
        )
    else
        if space == DrawBoxSpace.uv
            return get_draw_box_pixels(DrawBoxSpace.pixel, Box(
                min = min_inclusive(box) * resolution,
                size = size(box) * resolution
            ), resolution, box_is_1D_scalar)
        elseif space == DrawBoxSpace.pixel
            mi = min_inclusive(box)
            ma = max_inclusive(box)
            V = Vec{N, Int32}
            if R == Int32
                b_min = convert(V, clamp(mi, one(V), resolution))
                b_max = convert(V, max(b_min, clamp(ma, one(V), resolution)))
                return Box(min=b_min, max=b_max)
            elseif R == Float32
                b_min = convert(V, clamp(floor(V, mi), one(V), resolution))
                b_max = convert(V, max(b_min, clamp(floor(V, ma), one(V), resolution)))
                return Box(min=b_min, max=b_max)
            else
                error("Unhandled: ", R)
            end
        else
            error("Unhandled: ", space)
        end
    end
end

const DrawBoxRule = Union{Tuple{Val{:whitelist}, CellTypeSet},
                          Tuple{Val{:blacklist}, CellTypeSet},
                          Nothing}
check_draw_box_rule(rule::Nothing, u::UInt8) = true
check_draw_box_rule(rule::Tuple{Val{:whitelist}, CellTypeSet}, u::UInt8) = u in rule[2]
check_draw_box_rule(rule::Tuple{Val{:blacklist}, CellTypeSet}, u::UInt8) = !(u in rule[2])


#######
#  Op

"Conditionally fills a color in a rectangular area of pixels"
struct MarkovOpDrawBox{N, TRule<:DrawBoxRule} <: AbstractMarkovOp
    value::UInt8

    space::E_DrawBoxSpace
    box_is_1D_scalar::Bool # Only important if N=1
    box::Union{Box{N, Float32}, Box{N, Int32}}

    rule::TRule
    mask::Union{Nothing, Float32, NTuple{2, Float32}}
end

struct MarkovOpDrawBox_State{N, TMask<:Optional{MaskGrid{N}}}
    pixels::Bplus.Math.VecRange{N, Int32}
    next_pixel::Optional{Vec{N, Int32}}
    mask_grid::TMask
    mask_level::Float32
end

function markov_op_initialize(b::MarkovOpDrawBox{NBox, TRule},
                              grid::CellGrid{NGrid},
                              rng::PRNG, context::MarkovOpContext
                             ) where {NBox, NGrid, TRule}
    box = get_draw_box_pixels(
        b.space, b.box,
        convert(Vec{NGrid, Int32}, vsize(grid)),
        b.box_is_1D_scalar
    )
    pixel_range = min_inclusive(box):max_inclusive(box)

    iter_start = iterate(pixel_range)
    if isnothing(iter_start)
        return nothing
    else
        (next_pos, next_state) = iter_start
        @markovjunior_assert(next_pos == next_state, "Iterator works differently than I thought")

        mask_grid = if isnothing(b.mask)
            nothing
        else
            a = markov_allocator_acquire_array(context.allocator, size(grid), Float32)
            rand!(rng, a)
            a
        end

        mask_level = if isnothing(b.mask)
            # Value doesn't matter
            1.0f0
        elseif b.mask isa Float32
            b.mask
        elseif b.mask isa NTuple{2, Float32}
            lerp(b.mask..., rand(rng, Float32))
        else
            error("Unhandled ", typeof(b.mask))
        end

        return MarkovOpDrawBox_State(pixel_range, next_pos, mask_grid, mask_level)
    end
end
function markov_op_iterate(b::MarkovOpDrawBox{NBox, TRule},
                           state::MarkovOpDrawBox_State{NGrid, TMaskGrid},
                           grid::CellGrid{NGrid},
                           rng::PRNG, context::MarkovOpContext,
                           ticks_left::Ref{Optional{Int}}
                          ) where {NBox, NGrid, TRule, TMaskGrid}
    function apply_at(cell_idx)
        if check_draw_box_rule(b.rule, grid[cell_idx]) &&
           (TMaskGrid == Nothing || state.mask_grid[cell_idx] < state.mask_level)
        #begin
            grid[cell_idx] = b.value
        end
    end

    # If the user wants to finish asap, and we made it to the beginning of a row,
    #    do a quick fill of the rest.
    if isnothing(ticks_left[]) && (state.next_pixel.x == first(state.pixels).x)
        foreach(apply_at, first(state.pixels):last(state.pixels))
        markov_op_cancel(b, state, context)
        return nothing
    # Otherwise apply to the next pixel and try to advance.
    else
        apply_at(state.next_pixel)
        if exists(ticks_left[])
            ticks_left[] -= 1
        end

        try_iterate = iterate(state.pixels, state.next_pixel)
        return if isnothing(try_iterate)
            markov_op_cancel(b, state, context)
            nothing
        else
            @markovjunior_assert(try_iterate[1] == try_iterate[2],
                                 "Iterator works differently than I thought")
            typeof(state)(
                state.pixels,
                try_iterate[1],
                state.mask_grid, state.mask_level
            )
        end
    end
end
function markov_op_cancel(b::MarkovOpDrawBox{NBox, TRule},
                          state::MarkovOpDrawBox_State{NGrid, TMaskGrid},
                          context::MarkovOpContext
                         ) where {NBox, NGrid, TRule, TMaskGrid}
    if TMaskGrid != Nothing
        markov_allocator_release_array(context.allocator, state.mask_grid)
    end
    return nothing
end
markov_op_min_dimension(b::MarkovOpDrawBox{NBox}) where {NBox} = if b.box_is_1D_scalar
    1
else
    NBox
end


#############
#  DSL

dsl_string_box_rule(r::Nothing) = ""
dsl_string_box_rule(r::Tuple{Val{:whitelist}, CellTypeSet}) = "+$(dsl_string(r[2]))"
dsl_string_box_rule(r::Tuple{Val{:blacklist}, CellTypeSet}) = "-$(dsl_string(r[2]))"

dsl_string(b::MarkovOpDrawBox) = string(
    "@fill '", dsl_string(b.value), "'",
    " ", b.space, "(",
    "min=",
        if b.box_is_1D_scalar
            min_inclusive(b.box).x
        else
            min_inclusive(b.box).data
        end,
    ", size=",
        if b.box_is_1D_scalar
            size(b.box).x
        else
            size(b.box).data
        end,
    ")",
    " ", dsl_string_box_rule(b.rule),
    " ", dsl_string_rewrite_mask(b.mask)
)

function parse_markovjunior_op(::Val{Symbol("@fill")},
                               inputs::MacroParserInputs,
                               loc, expr_args, full_line)
    # The @capture macro has a bug with the '|' operator so we don't use it here.
    exCol = nothing
    exSpace = nothing
    exArgNameA = nothing
    exArgValA = nothing
    exArgNameB = nothing
    exArgValB = nothing
    exRuleAdd = nothing
    exRuleSub = nothing
    exMask = nothing
    if !@capture(full_line,  @fill exCol_Char exSpace_Symbol(exArgNameA_=exArgValA_, exArgNameB_=exArgValB_)                       ) &&
        !@capture(full_line, @fill exCol_Char exSpace_Symbol(exArgNameA_=exArgValA_, exArgNameB_=exArgValB_) +exRuleAdd_           ) &&
        !@capture(full_line, @fill exCol_Char exSpace_Symbol(exArgNameA_=exArgValA_, exArgNameB_=exArgValB_) -exRuleSub_           ) &&
        !@capture(full_line, @fill exCol_Char exSpace_Symbol(exArgNameA_=exArgValA_, exArgNameB_=exArgValB_)               %exMask_) &&
        !@capture(full_line, @fill exCol_Char exSpace_Symbol(exArgNameA_=exArgValA_, exArgNameB_=exArgValB_) +exRuleAdd_   %exMask_) &&
        !@capture(full_line, @fill exCol_Char exSpace_Symbol(exArgNameA_=exArgValA_, exArgNameB_=exArgValB_) -exRuleSub_   %exMask_) &&
       !@capture(full_line,  @fill exCol_Char                                                                                      ) &&
        !@capture(full_line, @fill exCol_Char                                                                +exRuleAdd_           ) &&
        !@capture(full_line, @fill exCol_Char                                                                -exRuleSub_           ) &&
        !@capture(full_line, @fill exCol_Char                                                                              %exMask_) &&
        !@capture(full_line, @fill exCol_Char                                                                +exRuleAdd_   %exMask_) &&
        !@capture(full_line, @fill exCol_Char                                                                -exRuleSub_   %exMask_)
    #begin
        raise_parse_error(loc, inputs,
                          "Invalid format! Expected `@fill 'C' [S(A=N, B=M)] [rule] [mask]`")
    end

    # Fill in defaults.
    if isnothing(exSpace)
        exSpace = :uv
        exArgNameA = :min
        exArgValA = 0
        exArgNameB = :max
        exArgValB = 1
    end

    # Do some quick error-checking of expression types.
    !in(exArgNameA, (:min, :max, :size, :center)) && raise_parse_error(
        loc, inputs,
        "Box space arguments should be two of `min`, `max`, `size`, and `center`! Got `",
          exArgNameA, "`"
    )
    !in(exArgNameB, (:min, :max, :size, :center)) && raise_parse_error(
        loc, inputs,
        "Box space arguments should be two of `min`, `max`, `size`, and `center`! Got `",
          exArgNameB, "`"
    )
    !(exRuleAdd isa Optional{Symbol}) && raise_parse_error(
        loc, inputs,
        "Whitelist rule should be formatted like `+abc`, got `+", exRuleAdd, "`"
    )
    !(exRuleSub isa Optional{Symbol}) && raise_parse_error(
        loc, inputs,
        "Blacklist rule should be formatted like `-abc`, got `-", exRuleSub, "`"
    )
    (exArgNameA == exArgNameB) && raise_parse_error(
        loc, inputs,
        "Used the same space parameter twice: `", exArgNameA, "`"
    )

    # Parse the mask.
    mask = if isnothing(exMask)
        nothing
    elseif @capture(exMask, maskAExpr_Real:maskBExpr_Real)
        convert.(Ref(Float32), (maskAExpr, maskBExpr))
    elseif exMask isa Real
        convert(Float32, exMask)
    else
        raise_parse_error(loc, inputs, "Expected mask to be `%x` or `%(x:y)`; got `%$exMask`")
    end

    # Parse the rule.
    rule = if exists(exRuleAdd)
        (Val(:whitelist), CellTypeSet(string(exRuleAdd)))
    elseif exists(exRuleSub)
        (Val(:blacklist), CellTypeSet(string(exRuleSub)))
    else
        nothing
    end

    # Parse the space type.
    space = try
        DrawBoxSpace.from(Val(exSpace))
    catch e
        raise_parse_error(loc, inputs,
                       "Invald space! Got `", space, "`; expected one of `",
                          string.(DrawBoxSpace.instances()), "`")
    end

    # Parse the space arguments.
    parse_space_arg(arg, name)::Union{Int32, Float32, Tuple{Vararg{Int32}}, Tuple{Vararg{Float32}}} =
        if arg isa Integer
            convert(Int32, arg)
        elseif arg isa Real
            convert(Float32, arg)
        elseif Base.isexpr(arg, :tuple) && all(a -> a isa Real, arg.args)
            if all(a -> a isa Integer, arg.args)
                Tuple(convert(Int32, a) for a in arg.args)
            else
                Tuple(convert(Float32, a) for a in arg.args)
            end
        else
            raise_parse_error(loc, inputs,
                           "Invalid value for argument `", name, "`! ",
                             "Expected scalar or vector `(x, y, ...)`; got `",
                             arg, "`")
        end
    argValA = parse_space_arg(exArgValA, exArgNameA)
    argValB = parse_space_arg(exArgValB, exArgNameB)

    # Reconcile the space arguments:
    #   * If one is scalar and one is vector, promote the scalar
    if (argValA isa Tuple) && !(argValB isa Tuple)
        argValB = ntuple(i -> argValB, length(argValA))
    elseif !(argValA isa Tuple) && (argValB isa Tuple)
        argValA = ntuple(i -> argValA, length(argValB))
    end
    #   * If at least one is float, promote the other to float
    if (eltype(argValA) == Float32) || (eltype(argValB) == Float32)
        argValA = convert.(Ref(Float32), argValA)
        argValB = convert.(Ref(Float32), argValB)
    end
    # If the rect is in UV space, make sure the values are float and not integer.
    # Integer boxes with the same min and max are size 1,
    #    so thin lines would get casted to cover the entire space!
    if (space == DrawBoxSpace.uv) && all(f -> f isa Integer, tuple(argValA..., argValB...))
        argValA = convert.(Ref(Float32), argValA)
        argValB = convert.(Ref(Float32), argValB)
    end

    # Finally, put the space arguments into an actual value.
    box_is_1D_scalar = argValA isa Real
    space_box = if argValA isa Tuple
        Box(namedtuple(
            (exArgNameA, exArgNameB),
            (Vec(argValA...), Vec(argValB...))
        ))
    elseif argValA isa Real
        Box(namedtuple(
            (exArgNameA, exArgNameB),
            (Vec(argValA), Vec(argValB))
        ))
    else
        error("Unhandled: ", typeof(argValA))
    end

    # Get the value of the fill color.
    if !haskey(CELL_CODE_BY_CHAR, exCol)
        raise_parse_error(loc, inputs,
                       "Unsupported color value '", exCol, "'! ",
                         "Supported are [ ", iter_join(keys(CELL_CODE_BY_CHAR), ", ")..., "]")
    end
    col = CELL_CODE_BY_CHAR[exCol]

    return MarkovOpDrawBox(col, space, box_is_1D_scalar, space_box, rule, mask)
end