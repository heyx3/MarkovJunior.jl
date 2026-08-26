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

function markov_algo_run(op::MarkovOpDrawBox{NBox, TRule},
                         algo::MarkovAlgorithm, algo_state::AlgoState,
                         inherited_bias_tuple::Tuple{Vararg{AbstractMarkovBias}},
                         inherited_bias_state_tuple::Tuple,
                         ::Val{NGrid} = Val(ndims(algo_state.grid))
                        )::Tuple{Bool, typeof(inherited_bias_state_tuple)} where {NBox, TRule, NGrid}
    box = get_draw_box_pixels(
        op.space, op.box,
        convert(Vec{NGrid, Int32}, vsize(algo_state.grid)),
        op.box_is_1D_scalar
    )

    mask_grid = if isnothing(op.mask)
        nothing
    else
        a::Array{Float32, NGrid} = markov_allocator_acquire_array(algo_state.allocator,
                                                                  size(algo_state.grid), Float32)
        rand!(algo_state.rng, a)
        a
    end
    previous_value_grid = if isempty(inherited_bias_tuple)
        nothing
    else
        a = markov_allocator_acquire_array(algo_state.allocator, size(algo_state.grid), UInt8)
        box_a = min_inclusive(box)
        box_b = max_inclusive(box)
        grid_slice = ntuple(i -> box_a[i]:box_b[i], Val(NGrid))
        previous_value_grid .= algo_state.grid[grid_slice...]
        a
    end

    # Make sure to de-allocate the arrays no matter what.
    made_changes::Bool = try

        # Pick a mask.
        mask_level = if isnothing(op.mask)
            # Value doesn't matter
            1.0f0
        elseif op.mask isa Float32
            op.mask
        elseif op.mask isa NTuple{2, Float32}
            lerp(op.mask..., rand(algo_state.rng, Float32))
        else
            error("Unhandled ", typeof(op.mask))
        end

        # Make the mask data type-stable before entering the loop.
        ((mask_grid, mask_level) -> begin

            mc::Bool = false
            for pixel in min_inclusive(box):max_inclusive(box)
                if check_draw_box_rule(op.rule, algo_state.grid[pixel]) &&
                    (isnothing(mask_grid) || (mask_grid[pixel] < mask_level))
                #begin
                    mc = true
                    algo_state.grid[pixel] = op.value
                    markov_algo_tick(algo_state, 1)
                end
            end
            inherited_bias_state_tuple = markov_bias_update.(
                inherited_bias_tuple, inherited_bias_state_tuple,
                Ref(algo), Ref(algo_state),
                Ref(box), Ref(previous_value_grid)
            )
            return mc

        end)(mask_grid, mask_level)

    finally
        exists(previous_value_grid) && markov_allocator_release_array(algo_state.allocator, previous_value_grid)
        exists(mask_grid) && markov_allocator_release_array(algo_state.allocator, mask_grid)
    end

    markov_algo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY - 1)
    return (made_changes, inherited_bias_state_tuple)
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
        DrawBoxSpace.from(exSpace)
    catch e
        raise_parse_error(loc, inputs,
                          "Invald space! Got `", exSpace, "`; expected one of `",
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

    inputs.min_dims = max(inputs.min_dims,
        if box_is_1D_scalar
            1
        else
            length(min_inclusive(space_box))
        end
    )
    return MarkovOpDrawBox(col, space, box_is_1D_scalar, space_box, rule, mask)
end