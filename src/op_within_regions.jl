#NOTE: this file takes some things from op_rewrite.jl and op_sequence.jl, and should be included after them.

"
An array wrapper which tracks whether each element has a particular Int32 value.
Used to mask the cell grid by 'region ID', for `@within_regions`.
"
struct RegionLabelMask{N, TLabels<:AbstractArray{Int32, N}} <: AbstractArray{Bool, N}
    labels::TLabels
    id::Int32
end
Base.size(m::RegionLabelMask) = size(m.labels)
Base.IndexStyle(::Type{<:RegionLabelMask}) = IndexCartesian()
Base.getindex(m::RegionLabelMask{N}, I::Vararg{Int, N}) where {N} = (m.labels[I...] == m.id)

"Data for one rule of a `@within_regions` Op"
struct WithinRegionsRule
    # The color set for regions this rule can apply to.
    colors::CellTypeSet
    # The "dummy color" that's used temporarily for pixels which are outside our regions
    #   but still nearby and with matching colors.
    mask_color::UInt8
    # The operation(s) to run on our regions.
    op::AbstractMarkovOp
end
Base.:(==)(a::WithinRegionsRule, b::WithinRegionsRule) = (a.colors == b.colors) &&
                                                         (a.mask_color == b.mask_color) &&
                                                         (a.op == b.op)
"The color which hides other regions' cells, if a rule doesn't specify one"
const WITHIN_REGIONS_DEFAULT_MASK = CELL_CODE_BY_CHAR['b']

"One region found by a `@within_regions` Op"
struct WithinRegionsElement{N}
    rule_idx::Int
    id::Int32 # A.k.a. its "label"
    n_cells::Int
    min::CellIdx{N}
    max::CellIdx{N}
end
"For each axis in the given region, returns the range of its bounding box"
(within_regions_ranges(r::WithinRegionsElement{N})::NTuple{N}) where {N} = ntuple(i -> r.min[i]:r.max[i], Val(N))


"Splits the grid into distinct regions (not necessarily rectangular), then runs Ops on each one"
struct MarkovOpWithinRegions{TBiases<:Tuple{Vararg{AbstractMarkovBias}}} <: AbstractMarkovOp
    threshold::Optional{Threshold}
    inherit_biases::Bool # Whether our biases get passed into each region's Ops
    border::Int # How much extra space to provide around each region, for border queries
    spread_axes::Optional{Vector{Int}} # Whitelist of axes the region can span
    rules::Vector{WithinRegionsRule}
    biases::TBiases
end
const WITHIN_REGIONS_DEFAULT_BORDER = 1

"
Finds the distinct regions for a `@within_regions` op,
  writing their data into `out_regions` and a segmentation into `out_labels`.
"
function find_within_regions!(out_labels::Array{Int32, N},
                              out_regions::Vector{WithinRegionsElement{N}},
                              queue::Vector{CellIdx{N}},
                              grid::CellGrid{N},
                              rules::Vector{WithinRegionsRule},
                              spread_axes::Optional{Vector{Int}},
                              rng::PRNG) where {N}
    VI = CellIdx{N}

    #=
        Rules may share colors, though each must have at least one unique color,
          so regions are found in two passes:

        1. Each rule in order claims the regions which start from one of its unique colors.
        2. Extra regions are built from shared colors in the leftover cells,
            and assigned to a random applicable rule.
    =#

    grid_size = convert(VI, vsize(grid))
    fill!(out_labels, zero(Int32))

    # Find the colors which are in more than one rule.
    shared_colors = CellTypeSet()
    seen_colors = CellTypeSet()
    for rule in rules
        shared_colors = union(shared_colors, intersect(seen_colors, rule.colors))
        seen_colors = union(seen_colors, rule.colors)
    end

    function claim_new_region!(start::CellIdx{N}, rule_i::Int)
        new_id = convert(Int32, length(out_regions) + 1)
        out_labels[start] = new_id

        (region_min, region_max) = (start, start)
        push!(queue, start)
        n_cells::Int = 0
        while !isempty(queue)
            next_cell = pop!(queue)
            n_cells += 1

            region_min = min(region_min, next_cell)
            region_max = max(region_max, next_cell)

            for axis in 1:N, dir in (-one(Int32), one(Int32))
                if exists(spread_axes) && !in(axis, spread_axes)
                    continue
                end

                neighbor_cell = next_cell + VI(i -> (i == axis) ? dir : zero(Int32))
                if (neighbor_cell[axis] in 1:grid_size[axis]) &&
                    iszero(out_labels[neighbor_cell]) && (grid[neighbor_cell] in rules[rule_i].colors)
                #begin
                    out_labels[neighbor_cell] = new_id
                    push!(queue, neighbor_cell)
                end
            end
        end

        push!(out_regions, WithinRegionsElement{N}(
            rule_i, new_id, n_cells,
            region_min, region_max
        ))
    end

    # 1. Each rule in order claims the regions which start from one of its unique colors.
    for (rule_i, rule) in enumerate(rules), start_cell::VI in one(VI):grid_size
        start_color = grid[start_cell]
        if (start_color in rule.colors) && !(start_color in shared_colors) && iszero(out_labels[start_cell])
            claim_new_region!(start_cell, rule_i)
        end
    end

    # 2. Extra regions are built from shared colors in the leftover cells,
    #      and assigned to a random applicable rule.
    if !isempty(shared_colors)
        for start_cell::VI in one(VI):grid_size
            start_color = grid[start_cell]

            if (start_color in shared_colors) && iszero(out_labels[start_cell])
                is_candidate(rule::WithinRegionsRule) = start_color in rule.colors
                pick::Int = rand(rng, 1:count(is_candidate, rules))
                rule_i = findfirst(rules) do r
                    if !is_candidate(r)
                        return false
                    else
                        pick -= 1
                        return pick == 0
                    end
                end

                claim_new_region!(start_cell, rule_i)
            end
        end
    end

    return nothing
end

"
Picks the next region to run, and returns its index in `remaining`.
If no region is available, returns `nothing`.

Biases are applied to influence our choice of region.
"
function pick_within_region(options::Vector{Int},
                            region_data::Vector{WithinRegionsElement{N}},
                            labels::Array{Int32, N},
                            biases::Tuple, bias_states::Tuple,
                            algo::MarkovAlgorithm, algo_state::AlgoState
                           )::Optional{Int} where {N}
    # With no choices or no biases, this is trivial.
    if isempty(options)
        return nothing
    elseif isempty(biases)
        return rand(algo_state.rng, 1:length(options))
    end

    # Otherwise we have to find the option with the strongest bias.
    best_i::Optional{Int} = nothing
    best_desirability = -Inf32
    n_ties = 0
    for (option_i, region_i) in enumerate(options)
        region = region_data[region_i]
        area = CellMask(
            RegionLabelMask(view(labels, within_regions_ranges(region)...), region.id),
            region.min - one(Int32)
        )

        desirability = sum_bias_desirabilities(biases, bias_states, algo, algo_state, area)
        if isnothing(desirability)
            continue
        elseif desirability > best_desirability
            (best_i, best_desirability, n_ties) = (option_i, desirability, 1)
        elseif desirability == best_desirability
            # Reservoir sampling, so every tied region is equally likely.
            n_ties += 1
            if rand(algo_state.rng, 1:n_ties) == 1
                best_i = option_i
            end
        end
    end
    return best_i
end

function markov_algo_run(op::MarkovOpWithinRegions,
                         algo::MarkovAlgorithm, algo_state::AlgoState,
                         inherited_biases::Tuple{Vararg{AbstractMarkovBias}},
                         inherited_bias_states::Tuple,
                         grid::CellGrid{N} = algo_state.grid
                        )::OpChangedGrid where {N}
    # Dispatch on the allocator type before making allocations.
    (alloc -> begin

    labels::Array{Int32, N} = markov_allocator_acquire_array(alloc, size(grid), Int32)
    regions::Vector{WithinRegionsElement{N}} = markov_allocator_acquire_array(alloc, (), WithinRegionsElement{N})
    queue::Vector{CellIdx{N}} = markov_allocator_acquire_array(alloc, (), CellIdx{N})
    remaining::Vector{Int} = markov_allocator_acquire_array(alloc, (), Int)
    try
        # Find the regions.
        find_within_regions!(labels, regions, queue, grid, op.rules, op.spread_axes, algo_state.rng)
        @logic_logln("Found ", length(regions), " regions")
        if isempty(regions)
            return false
        end
        append!(remaining, 1:length(regions))

        # Get a threshold.
        n_to_run = if isnothing(op.threshold)
            length(regions)
        else
            get_threshold(op.threshold, grid, algo_state.grid_area, algo_state.rng, algo.calc_macros)
        end

        # Execute.
        return markov_opbias_with_new_biases(op.biases, algo, algo_state,
                                             inherited_biases, inherited_bias_states
                                            ) do all_biases, all_bias_states
            total_result::OpChangedGrid = false
            n_ran = 0
            while n_ran < n_to_run
                # Pick the next region and remove it from the list.
                pick_i = pick_within_region(remaining, regions, labels, all_biases, all_bias_states, algo, algo_state)
                if isnothing(pick_i)
                    @logic_logln("No more regions are left (we did them all, or biases forbid the rest); stopping")
                    break
                end
                region = regions[remaining[pick_i]]
                deleteat!(remaining, pick_i)

                # Execute that region.
                region_result::OpChangedGrid = begin
                    border = convert(Int32, op.border)
                    view_min = max(region.min - border, one(CellIdx{N}))
                    view_max = min(region.max + border, convert(CellIdx{N}, vsize(grid)))
                    @logic_logln("Running region ", region.id, " (rule ", region.rule_idx,
                                 ") within ", view_min, " to ", view_max)

                    ranges = ntuple(i -> view_min[i]:view_max[i], Val(N))
                    sub_grid = view(grid, ranges...)
                    sub_labels = view(labels, ranges...)

                    markov_allocator_with_array(alloc, size(sub_grid), UInt8) do original_pixels::Array{UInt8, N}
                        copyto!(original_pixels, sub_grid)

                        rule = op.rules[region.rule_idx]
                        markov_opalgo_on_subgrid(
                            (inner_biases, inner_bias_states) -> markov_algo_run(
                                rule.op, algo, algo_state, inner_biases, inner_bias_states
                            ),
                            algo, algo_state,
                            sub_grid, BoxI{N}(min=view_min, max=view_max),
                            op.inherit_biases, all_biases, all_bias_states,

                            # Hide other cells with the rule's colors, then put back everything outside the region afterwards.
                            () -> for i in CartesianIndices(sub_grid)
                                if (sub_labels[i] != region.id) && (sub_grid[i] in rule.colors)
                                    sub_grid[i] = rule.mask_color
                                end
                            end,
                            () -> for i in CartesianIndices(sub_grid)
                                if sub_labels[i] != region.id
                                    sub_grid[i] = original_pixels[i]
                                end
                            end,

                            ;
                            n_cells = region.n_cells
                        )
                    end
                end

                # Update bookkeeping.
                total_result = (n_ran == 0) ?
                    region_result :
                    combine_op_results(total_result, region_result)
                n_ran += 1
                markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY)
            end

            markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY + 1)
            return total_result
        end

    finally
        # Release allocations.
        markov_allocator_release_array(alloc, labels)
        markov_allocator_release_array(alloc, regions)
        markov_allocator_release_array(alloc, queue)
        markov_allocator_release_array(alloc, remaining)
    end

    end)(algo_state.allocator) # (dispatch for type-stability)
end

# The inner Ops are rejected at parse-time if they reallocate.
markov_op_reallocates_grid(op::MarkovOpWithinRegions) = any(r -> markov_op_reallocates_grid(r.op), op.rules)


dsl_format(r::WithinRegionsRule) = string(
    dsl_format(r.colors),
    (r.mask_color == WITHIN_REGIONS_DEFAULT_MASK) ? "" : "($(dsl_format(r.mask_color)))",
    " => ", dsl_format(r.op)
)
dsl_format(op::MarkovOpWithinRegions) = string(
    "@within_regions ",

    # Args:
    isnothing(op.threshold) ? "" : "$(dsl_format(op.threshold)) ",
    op.inherit_biases ? "inherit_biases " : "",
    (op.border == WITHIN_REGIONS_DEFAULT_BORDER) ? "" : "border($(op.border)) ",
    isnothing(op.spread_axes) ? "" : string("along(", iter_join(op.spread_axes, ", ")..., ") "),

    # Rules:
    "begin\n    ",
    iter_join(Iterators.map(dsl_format, op.rules), "\n    ")...,
    "\nend",

    # Biases:
    if isempty(op.biases)
        ""
    else
        string(" begin\n    ", iter_join(Iterators.map(dsl_format, op.biases), "\n    ")..., "\nend")
    end
)

function parse_markovjunior_op(::Val{Symbol("@within_regions")}, inputs::MacroParserInputs,
                               loc, expr_args, expr_line)
    args = collect(Any, expr_args)

    # The beginning of the rules marks the end of the parameters.
    # Each rule is `colors => op` or `colors(mask) => op`,
    #    and they come as a single instance or a block of them.
    # None of the other args, including the threshold, can contain a `=>` expression or block.
    is_inline_rule(a) = Base.isexpr(a, :call) && (a.args[1] == :(=>))
    rule_expr_i = findfirst(a -> Base.isexpr(a, :block) || is_inline_rule(a), args)
    if isnothing(rule_expr_i)
        raise_parse_error(loc, inputs,
                          "Expected `@within_regions [threshold] [inherit_biases] [border(N)] [along(axes...)] ",
                            "begin ... end [biases]`, or a single rule in place of the block")
    elseif length(args) > rule_expr_i + 1
        raise_parse_error(loc, inputs, "Unexpected arguments after the biases: ", args[(rule_expr_i + 2):end])
    end

    # Parse the args.
    threshold::Optional{Threshold} = nothing
    inherit_biases = false
    border::Optional{Int} = nothing
    spread_axes::Optional{Vector{Int}} = nothing
    for setting in args[1:(rule_expr_i - 1)]
        if setting == :inherit_biases
            inherit_biases = true
        elseif @capture(setting, border(borderArgs__))
            if exists(border)
                raise_parse_error(loc, inputs, "`border()` is given more than once")
            elseif (length(borderArgs) != 1) || !isa(borderArgs[1], Integer) || (borderArgs[1] < 0)
                raise_parse_error(loc, inputs, "Expected `border(N)` with a whole number N >= 0; got `", setting, "`")
            end
            border = convert(Int, borderArgs[1])
        elseif @capture(setting, along(alongArgs__))
            if exists(spread_axes)
                raise_parse_error(loc, inputs, "`along()` is given more than once")
            elseif length(alongArgs) < 1
                raise_parse_error(loc, inputs, "`along()` needs at least one axis, like `along(x, z)`")
            end
            spread_axes = map(alongArgs) do a
                axis = if a isa Integer
                    convert(Int, a)
                elseif a isa Symbol
                    something(findfirst(==(a), (:x, :y, :z, :w)), 0)
                else
                    0
                end
                if axis < 1
                    raise_parse_error(loc, inputs, "Expected an axis like `x`, `y`, `z`, `w`, or a number; got `", a, "`")
                end
                inputs.min_dims = max(inputs.min_dims, axis)
                axis
            end
            if !allunique(spread_axes)
                raise_parse_error(loc, inputs, "`", setting, "` lists an axis more than once")
            end
        # Whatever is left must be a threshold expression, by elimination.
        elseif exists(threshold)
            raise_parse_error(loc, inputs,
                              "Unexpected argument (are there two thresholds?): `",
                                setting, "`")
        else
            threshold = parse_markovjunior_threshold(inputs, loc, setting)
        end
    end

    # Rule-parsing comes last, as it needs to validate the regional Ops' biases
    #   and that may include the outer ones.
    function parse_rules(outer_biases)
        rules = WithinRegionsRule[ ]

        # Without inherited biases, the inner Ops shouldn't be validated against the outer biases.
        outer_bias_stack = inputs.bias_stack
        if !inherit_biases
            inputs.bias_stack = Stack{Vector{AbstractMarkovBias}}(16)
        end
        try
            rule_lines = is_inline_rule(args[rule_expr_i]) ?
                Any[ args[rule_expr_i] ] :
                args[rule_expr_i].args
            parse_markovjunior_block(rule_lines) do rule_loc, line
                rule_loc = something(rule_loc, loc)
                with_parser_stacktrace(inputs, "Rule #$(length(rules) + 1)") do
                    colorsExpr = nothing
                    maskExpr = nothing
                    opExpr = nothing
                    if !@capture(line, colorsExpr_(maskExpr_) => opExpr_) &&
                       !@capture(line, colorsExpr_ => opExpr_)
                    #begin
                        raise_parse_error(rule_loc, inputs,
                                          "Expected syntax like `RGB => @rewrite ...`, ",
                                            "or `RGB(w) => @rewrite ...`, but got `",
                                            line, "`")
                    end

                    # Interpret the mask.
                    mask_color = WITHIN_REGIONS_DEFAULT_MASK
                    if !isnothing(maskExpr)
                        mask_str = string(maskExpr)
                        if !isa(maskExpr, Symbol) || (length(mask_str) != 1) ||
                           !haskey(CELL_CODE_BY_CHAR, mask_str[1])
                        #begin
                            raise_parse_error(rule_loc, inputs,
                                              "Expected one mask color like `w`; got `", maskExpr, "`")
                        end
                        mask_color = CELL_CODE_BY_CHAR[mask_str[1]]
                    end

                    # Interpret the color key.
                    if !isa(colorsExpr, Symbol)
                        raise_parse_error(rule_loc, inputs,
                                          "Expected a rule like `MO => @rewrite M => O`; got `", line, "`")
                    end

                    # Interpret the op.
                    if !Base.isexpr(opExpr, :macrocall)
                        raise_parse_error(rule_loc, inputs,
                                          "A rule runs one Op; for several, use `@sequence begin ... end`. Got `",
                                            opExpr, "`")
                    end
                    ops = parse_markovjunior_sequence_main_block(opExpr, rule_loc, inputs)

                    # Put them all together.
                    color_set = try
                        CellTypeSet(string(colorsExpr))
                    catch e
                        raise_parse_error(rule_loc, inputs, "Invalid colors `", colorsExpr, "`")
                    end
                    check_markovjunior_scoped_ops(inputs, rule_loc, "@within_regions", ops)
                    push!(rules, WithinRegionsRule(color_set, mask_color, only(ops)))
                end
            end
        finally
            inputs.bias_stack = outer_bias_stack
        end
        if isempty(rules)
            raise_parse_error(loc, inputs, "`@within_regions` needs at least one rule")
        end

        # Rules may share colors, but each one needs a color of its own to start its regions from.
        for (rule_i, rule) in enumerate(rules)
            other_colors = foldl(union, (r.colors for (i, r) in enumerate(rules) if i != rule_i);
                                 init = CellTypeSet())
            if isempty(setdiff(rule.colors, other_colors))
                raise_parse_error(loc, inputs,
                                  "Rule ", rule_i, " (`", dsl_format(rule.colors), "`) ",
                                    "only has colors which other rules also have; ",
                                    "each rule needs at least one color of its own")
            end
        end

        return MarkovOpWithinRegions(threshold, inherit_biases,
                                     something(border, WITHIN_REGIONS_DEFAULT_BORDER), spread_axes,
                                     rules, Tuple(outer_biases))
    end

    return if length(args) > rule_expr_i
        with_parsed_markovjunior_bias_statement(parse_rules, inputs, loc, args[rule_expr_i + 1])
    else
        parse_rules(AbstractMarkovBias[ ])
    end
end
