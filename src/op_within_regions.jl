#NOTE: this file takes some things from op_rewrite.jl and op_sequence.jl, and should be included after them.

"
A mask of the cells with a specific label, out of a grid of labels.
Lets a `CellMask` describe one region without allocating anything.
"
struct RegionLabelMask{N, TLabels<:AbstractArray{Int32, N}} <: AbstractArray{Bool, N}
    labels::TLabels
    id::Int32
end
Base.size(m::RegionLabelMask) = size(m.labels)
Base.IndexStyle(::Type{<:RegionLabelMask}) = IndexCartesian()
Base.getindex(m::RegionLabelMask{N}, I::Vararg{Int, N}) where {N} = (m.labels[I...] == m.id)

"One rule of a `@within_regions` Op: the colors which make up its regions, and the Ops to run on each region"
struct WithinRegionsRule
    colors::CellTypeSet
    ops::Vector{AbstractMarkovOp}
end
Base.:(==)(a::WithinRegionsRule, b::WithinRegionsRule) = (a.colors == b.colors) && (a.ops == b.ops)

"One region found by a `@within_regions` Op"
struct WithinRegionsRegion{N}
    rule_idx::Int
    id::Int32 # Its label in the grid of labels
    min::CellIdx{N}
    max::CellIdx{N}
end
within_regions_ranges(r::WithinRegionsRegion{N}) where {N} = ntuple(i -> r.min[i]:r.max[i], Val(N))


"
Runs a sequence of Ops within each connected region of certain colors, one region at a time,
  as if that region's bounding box were the entire grid.
Cells in the bounding box that aren't part of the region are temporarily replaced with `clear_color`.
"
struct MarkovOpWithinRegions{TBiases<:Tuple{Vararg{AbstractMarkovBias}}} <: AbstractMarkovOp
    threshold::Optional{Threshold}
    inherit_biases::Bool
    clear_color::UInt8
    rules::Vector{WithinRegionsRule}
    biases::TBiases

    # For each cell type, the index of the rule that includes it (or 0 for none).
    # Each color may only be in one rule.
    rule_by_cell::NTuple{N_CELL_TYPES, Int}
end
MarkovOpWithinRegions(threshold, inherit_biases, clear_color, rules, biases) = MarkovOpWithinRegions(
    threshold, inherit_biases, clear_color, rules, biases,
    ntuple(Val(N_CELL_TYPES)) do i
        something(findfirst(r -> (UInt8(i - 1) in r.colors), rules), 0)
    end
)

"
Labels each connected region of cells whose colors all map to the same rule,
  pushing info about each one into `out_regions`.
Cells are connected along each axis (not diagonally).
"
function find_within_regions!(out_labels::Array{Int32, N},
                              out_regions::Vector{WithinRegionsRegion{N}},
                              queue::Vector{CellIdx{N}},
                              grid::CellGrid{N},
                              rule_by_cell::NTuple{N_CELL_TYPES, Int}) where {N}
    VI = CellIdx{N}
    grid_size = convert(VI, vsize(grid))
    fill!(out_labels, zero(Int32))
    for start::VI in one(VI):grid_size
        rule_i = rule_by_cell[grid[start] + 1]
        if iszero(rule_i) || !iszero(out_labels[start])
            continue
        end

        id = convert(Int32, length(out_regions) + 1)
        (region_min, region_max) = (start, start)
        out_labels[start] = id
        push!(queue, start)
        while !isempty(queue)
            cell = pop!(queue)
            region_min = min(region_min, cell)
            region_max = max(region_max, cell)
            for axis in 1:N, dir in (-one(Int32), one(Int32))
                neighbor = cell + VI(i -> (i == axis) ? dir : zero(Int32))
                if (neighbor[axis] >= 1) && (neighbor[axis] <= grid_size[axis]) &&
                   iszero(out_labels[neighbor]) &&
                   (rule_by_cell[grid[neighbor] + 1] == rule_i)
                #begin
                    out_labels[neighbor] = id
                    push!(queue, neighbor)
                end
            end
        end
        push!(out_regions, WithinRegionsRegion{N}(rule_i, id, region_min, region_max))
    end
    return nothing
end

"
Picks the next region to run, as an index into `remaining`, or `nothing` if no region is allowed.
Like `@rewrite`, this takes the most desirable region according to the biases, breaking ties randomly.
"
function pick_within_region(remaining::Vector{Int}, regions::Vector{WithinRegionsRegion{N}},
                            labels::Array{Int32, N},
                            biases::Tuple, bias_states::Tuple,
                            algo::MarkovAlgorithm, algo_state::AlgoState
                           )::Optional{Int} where {N}
    if isempty(biases)
        return rand(algo_state.rng, 1:length(remaining))
    end

    best_i::Optional{Int} = nothing
    best_desirability = -Inf32
    n_tied = 0
    for (i, region_i) in enumerate(remaining)
        region = regions[region_i]
        area = CellMask(
            RegionLabelMask(view(labels, within_regions_ranges(region)...), region.id),
            region.min - one(Int32)
        )
        desirability = sum_bias_desirabilities(biases, bias_states, algo, algo_state, area)
        if isnothing(desirability)
            continue
        elseif desirability > best_desirability
            (best_i, best_desirability, n_tied) = (i, desirability, 1)
        elseif desirability == best_desirability
            # Reservoir sampling, so every tied region is equally likely.
            n_tied += 1
            if rand(algo_state.rng, 1:n_tied) == 1
                best_i = i
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
    alloc = algo_state.allocator
    labels::Array{Int32, N} = markov_allocator_acquire_array(alloc, size(grid), Int32)
    regions::Vector{WithinRegionsRegion{N}} = markov_allocator_acquire_array(alloc, (), WithinRegionsRegion{N})
    queue::Vector{CellIdx{N}} = markov_allocator_acquire_array(alloc, (), CellIdx{N})
    remaining::Vector{Int} = markov_allocator_acquire_array(alloc, (), Int)
    try
        # Find the regions up-front; each one only changes itself, so they never need to be found again.
        find_within_regions!(labels, regions, queue, grid, op.rule_by_cell)
        @logic_logln("Found ", length(regions), " regions")
        if isempty(regions)
            return false
        end
        append!(remaining, 1:length(regions))
        n_to_run = if isnothing(op.threshold)
            length(regions)
        else
            get_threshold(op.threshold, grid, algo_state.rng)
        end

        return markov_opbias_with_new_biases(op.biases, algo, algo_state,
                                             inherited_biases, inherited_bias_states
                                            ) do all_biases, all_bias_states
            run_within_regions(op, algo, algo_state, grid, labels, regions, remaining, n_to_run,
                               all_biases, all_bias_states)
        end
    finally
        markov_allocator_release_array(alloc, labels)
        markov_allocator_release_array(alloc, regions)
        markov_allocator_release_array(alloc, queue)
        markov_allocator_release_array(alloc, remaining)
    end
end

"Picks and runs up to `n_to_run` of the remaining regions, one at a time"
function run_within_regions(op::MarkovOpWithinRegions, algo::MarkovAlgorithm, algo_state::AlgoState,
                            grid::CellGrid{N}, labels::Array{Int32, N},
                            regions::Vector{WithinRegionsRegion{N}}, remaining::Vector{Int}, n_to_run::Int,
                            biases::Tuple, bias_states::Tuple
                           )::OpChangedGrid where {N}
    result::OpChangedGrid = false
    n_ran = 0
    while (n_ran < n_to_run) && !isempty(remaining)
        pick_i = pick_within_region(remaining, regions, labels, biases, bias_states, algo, algo_state)
        if isnothing(pick_i)
            @logic_logln("No regions are allowed by the biases; stopping")
            break
        end
        region = regions[remaining[pick_i]]
        deleteat!(remaining, pick_i)

        region_result = run_within_region(op, algo, algo_state, grid, labels, region, biases, bias_states)
        result = (n_ran == 0) ? region_result : combine_op_results(result, region_result)
        n_ran += 1
        markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY)
    end

    markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY + 1)
    return result
end

"Runs the Ops for one region, on a view of its bounding box"
function run_within_region(op::MarkovOpWithinRegions, algo::MarkovAlgorithm, algo_state::AlgoState,
                           grid::CellGrid{N}, labels::Array{Int32, N},
                           region::WithinRegionsRegion{N},
                           biases::Tuple, bias_states::Tuple
                          )::OpChangedGrid where {N}
    @logic_logln("Running region ", region.id, " (rule ", region.rule_idx, ") within ",
                   region.min, " to ", region.max)
    alloc = algo_state.allocator
    ranges = within_regions_ranges(region)
    sub_grid = view(grid, ranges...)
    sub_labels = view(labels, ranges...)

    original_pixels::Array{UInt8, N} = markov_allocator_acquire_array(alloc, size(sub_grid), UInt8)
    copyto!(original_pixels, sub_grid)
    try
        ops = op.rules[region.rule_idx].ops
        return markov_opalgo_on_subgrid(
            algo, algo_state,
            sub_grid, BoxI{N}(min=region.min, max=region.max),
            op.inherit_biases, biases, bias_states,
            # Blank out the bounding box's cells which aren't part of this region,
            #    then put them back afterwards.
            () -> replace_outside_region!(sub_grid, sub_labels, region.id, op.clear_color),
            () -> replace_outside_region!(sub_grid, sub_labels, region.id, original_pixels)
        ) do inner_biases, inner_bias_states
            # The grid can't be reallocated in here, so the bias states never get replaced.
            inner_result::OpChangedGrid = nothing
            for inner_op in ops
                inner_result = combine_op_results(inner_result, markov_algo_run(
                    inner_op, algo, algo_state,
                    inner_biases, inner_bias_states
                ))
            end
            return inner_result
        end
    finally
        markov_allocator_release_array(alloc, original_pixels)
    end
end

"Sets every cell of the view which isn't in the given region, either to one color or from a grid of colors"
function replace_outside_region!(sub_grid::CellGrid{N}, sub_labels::AbstractArray{Int32, N}, id::Int32,
                                 values::Union{UInt8, AbstractArray{UInt8, N}}) where {N}
    for i in CartesianIndices(sub_grid)
        if sub_labels[i] != id
            sub_grid[i] = (values isa UInt8) ? values : values[i]
        end
    end
    return nothing
end

# The inner Ops are rejected at parse-time if they reallocate.
markov_op_reallocates_grid(op::MarkovOpWithinRegions) = any(r -> any(markov_op_reallocates_grid, r.ops), op.rules)


dsl_format(r::WithinRegionsRule) = string(
    dsl_format(r.colors), " => begin\n        ",
    iter_join(Iterators.map(dsl_format, r.ops), "\n        ")...,
    "\n    end"
)
dsl_format(op::MarkovOpWithinRegions) = string(
    "@within_regions ",
    isnothing(op.threshold) ? "" : "$(dsl_format(op.threshold)) ",
    op.inherit_biases ? "inherit_biases " : "",
    (op.clear_color == CELL_CODE_BY_CHAR['b']) ? "" : "'$(dsl_format(op.clear_color))' ",
    "begin\n    ",
    iter_join(Iterators.map(dsl_format, op.rules), "\n    ")...,
    "\nend",
    if isempty(op.biases)
        ""
    else
        string(" begin\n    ", iter_join(Iterators.map(dsl_format, op.biases), "\n    ")..., "\nend")
    end
)

function parse_markovjunior_op(::Val{Symbol("@within_regions")}, inputs::MacroParserInputs,
                               loc, expr_args, expr_line)
    args = collect(Any, expr_args)
    i_rules = findfirst(a -> Base.isexpr(a, :block), args)
    if isnothing(i_rules)
        raise_parse_error(loc, inputs,
                          "Expected `@within_regions [threshold] [inherit_biases] ['C'] begin ... end [biases]`")
    elseif length(args) > i_rules + 1
        raise_parse_error(loc, inputs, "Unexpected arguments after the biases: ", args[(i_rules + 2):end])
    end

    # Parse the settings before the rules.
    threshold::Optional{Threshold} = nothing
    inherit_biases = false
    clear_color = CELL_CODE_BY_CHAR['b']
    for setting in args[1:(i_rules - 1)]
        if setting == :inherit_biases
            inherit_biases = true
        elseif setting isa Char
            if !haskey(CELL_CODE_BY_CHAR, setting)
                raise_parse_error(loc, inputs, "Unsupported color value '", setting, "'")
            end
            clear_color = CELL_CODE_BY_CHAR[setting]
        elseif check_markovjunior_threshold_appearance(setting)
            threshold = parse_markovjunior_threshold(inputs, loc, setting)
        else
            raise_parse_error(loc, inputs,
                              "Expected a threshold, `inherit_biases`, or a 'C' color to clear with; got `",
                                setting, "`")
        end
    end

    function parse_rules(biases)
        # Without inherited biases, the inner Ops shouldn't be validated against the outer biases either.
        outer_bias_stack = inputs.bias_stack
        inherit_biases || (inputs.bias_stack = Stack{Vector{AbstractMarkovBias}}(16))
        rules = WithinRegionsRule[ ]
        try
            parse_markovjunior_block(args[i_rules].args) do rule_loc, line
                rule_loc = something(rule_loc, loc)
                with_parser_stacktrace(inputs, "Rule $(length(rules) + 1)") do
                    if !@capture(line, colors_ => ops_expr_) || !isa(colors, Symbol)
                        raise_parse_error(rule_loc, inputs,
                                          "Expected a rule like `MO => begin ... end`; got `", line, "`")
                    end
                    color_set = try
                        CellTypeSet(string(colors))
                    catch e
                        raise_parse_error(rule_loc, inputs, "Invalid colors `", colors, "`")
                    end
                    for (other_i, other_rule) in enumerate(rules)
                        shared = intersect(color_set, other_rule.colors)
                        if !isempty(shared)
                            raise_parse_error(rule_loc, inputs,
                                              "Colors `", dsl_format(shared), "` are already used by rule ", other_i,
                                                "; each color can only be in one rule")
                        end
                    end
                    ops = parse_markovjunior_sequence_main_block(ops_expr, rule_loc, inputs)
                    check_markovjunior_scoped_ops(inputs, rule_loc, "@within_regions", ops)
                    push!(rules, WithinRegionsRule(color_set, ops))
                end
            end
        finally
            inputs.bias_stack = outer_bias_stack
        end
        if isempty(rules)
            raise_parse_error(loc, inputs, "`@within_regions` needs at least one rule")
        end

        return MarkovOpWithinRegions(threshold, inherit_biases, clear_color, rules, Tuple(biases))
    end

    return if length(args) > i_rules
        with_parsed_markovjunior_bias_statement(parse_rules, inputs, loc, args[i_rules + 1])
    else
        parse_rules(AbstractMarkovBias[ ])
    end
end
