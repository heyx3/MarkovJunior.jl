#NOTE: this file takes some things from op_rewrite.jl, and should be included after it

#############
#  Data

const DrawBoxRule = Union{Tuple{Val{:whitelist}, CellTypeSet},
                          Tuple{Val{:blacklist}, CellTypeSet},
                          Nothing}
check_draw_box_rule(rule::Nothing, u::UInt8) = true
check_draw_box_rule(rule::Tuple{Val{:whitelist}, CellTypeSet}, u::UInt8) = u in rule[2]
check_draw_box_rule(rule::Tuple{Val{:blacklist}, CellTypeSet}, u::UInt8) = !(u in rule[2])


#######
#  Op

"Conditionally fills a color in a rectangular area of pixels"
struct MarkovOpDrawBox{TRule<:DrawBoxRule} <: AbstractMarkovOp
    value::UInt8
    box::CalcBox
    rule::TRule
    mask::Union{Nothing, Float32, NTuple{2, Float32}}
end

function markov_algo_run(op::MarkovOpDrawBox{TRule},
                         algo::MarkovAlgorithm, algo_state::AlgoState,
                         inherited_bias_tuple::Tuple{Vararg{AbstractMarkovBias}},
                         inherited_bias_state_tuple::Tuple,
                         grid::CellGrid{NGrid} = algo_state.grid
                        )::Bool where {TRule, NGrid}
    box = eval_calc_box(op.box, grid, algo_state.rng, algo.calc_macros)
    if isnothing(box)
        markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY - 1)
        return false
    end
    (box_a, box_b) = (min_inclusive(box), max_inclusive(box))
    grid_slice_idcs = ntuple(i -> box_a[i]:box_b[i], Val(NGrid))

    mask_grid = if isnothing(op.mask)
        nothing
    else
        a::Array{Float32, NGrid} = markov_allocator_acquire_array(algo_state.allocator,
                                                                  size(box).data, Float32)
        rand!(algo_state.rng, a)
        a
    end
    previous_value_grid = if isempty(inherited_bias_tuple)
        nothing
    else
        a2 = markov_allocator_acquire_array(algo_state.allocator, size(box).data, UInt8)
        a2 .= @view grid[grid_slice_idcs...]
        a2
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
            for pixel in box_a:box_b
                local_pixel = pixel - box_a + one(Vec{NGrid, Int32})
                if check_draw_box_rule(op.rule, grid[pixel]) &&
                    (isnothing(mask_grid) || (mask_grid[local_pixel] < mask_level))
                #begin
                    mc = true
                    grid[pixel] = op.value
                    markov_opalgo_tick(algo_state, 1)
                end
            end
            markov_opbias_update_all!(inherited_bias_tuple, inherited_bias_state_tuple,
                                      algo, algo_state, box, previous_value_grid)
            return mc

        end)(mask_grid, mask_level)

    finally
        exists(previous_value_grid) && markov_allocator_release_array(algo_state.allocator, previous_value_grid)
        exists(mask_grid) && markov_allocator_release_array(algo_state.allocator, mask_grid)
    end

    markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY - 1)
    return made_changes
end

markov_op_reallocates_grid(::MarkovOpDrawBox) = false


#############
#  DSL

dsl_format_box_rule(r::Nothing) = ""
dsl_format_box_rule(r::Tuple{Val{:whitelist}, CellTypeSet}) = "+$(dsl_format(r[2]))"
dsl_format_box_rule(r::Tuple{Val{:blacklist}, CellTypeSet}) = "-$(dsl_format(r[2]))"

dsl_format(b::MarkovOpDrawBox) = string(
    "@fill '", dsl_format(b.value), "'",
    " ", dsl_format(b.box),
    " ", dsl_format_box_rule(b.rule),
    " ", dsl_format_rewrite_mask(b.mask)
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
    !(exRuleAdd isa Optional{Symbol}) && raise_parse_error(
        loc, inputs,
        "Whitelist rule should be formatted like `+abc`, got `+", exRuleAdd, "`"
    )
    !(exRuleSub isa Optional{Symbol}) && raise_parse_error(
        loc, inputs,
        "Blacklist rule should be formatted like `-abc`, got `-", exRuleSub, "`"
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
        rule_set = try
            CellTypeSet(string(exRuleAdd))
        catch e
            raise_parse_error(nothing, inputs, "Invalid colors in rule `+", exRuleAdd, "`")
        end
        (Val(:whitelist), rule_set)
    elseif exists(exRuleSub)
        rule_set = try
            CellTypeSet(string(exRuleSub))
        catch e
            raise_parse_error(nothing, inputs, "Invalid colors in rule `-", exRuleSub, "`")
        end
        (Val(:blacklist), rule_set)
    else
        nothing
    end

    # Parse the box.
    box = with_parser_stacktrace(inputs, "Box") do
        parse_markovjunior_calc_box(
            inputs, loc,
            Expr(:call, exSpace, Expr(:kw, exArgNameA, exArgValA), Expr(:kw, exArgNameB, exArgValB)),
            CALC_BOX_ENV
        )
    end

    # Get the value of the fill color.
    if !haskey(CELL_CODE_BY_CHAR, exCol)
        raise_parse_error(
            loc, inputs,
            "Unsupported color value '", exCol, "'! ",
              "Supported are [ ", iter_join(keys(CELL_CODE_BY_CHAR), ", ")..., "]"
        )
    end
    col = CELL_CODE_BY_CHAR[exCol]

    return MarkovOpDrawBox(col, box, rule, mask)
end