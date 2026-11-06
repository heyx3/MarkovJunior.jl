##########################################
#  Shared by @upscale and @downscale

"
Per-axis integer factors for the scale Ops, like `(2, 3)`.
If `extend_last` is set (written as a trailing ellipsis, like `(2, 3...)`),
  the last factor applies to every further axis; otherwise further axes use 1.
"
struct ScaleOpsFactors
    factors::Vector{Int}
    extend_last::Bool
end
Base.:(==)(a::ScaleOpsFactors, b::ScaleOpsFactors) = (a.factors == b.factors) && (a.extend_last == b.extend_last)

"Gets the scale factor along the given axis"
scale_ops_factor(f::ScaleOpsFactors, axis::Int)::Int = if axis <= length(f.factors)
    f.factors[axis]
elseif f.extend_last
    f.factors[end]
else
    1
end

"Gets the scale factors along every axis of an N-dimensional grid"
scale_ops_factors(f::ScaleOpsFactors, ::Val{N}) where {N} = CellIdx{N}(a -> scale_ops_factor(f, a))

"
Gets the cell of a rule's block which covers the given pixel,
  where `local_pos` is the pixel's 0-based offset within its scaled block.
Size-1 block axes are extruded.
"
@inline function scale_ops_block_idx(block::Array{<:Any, NBlockDims},
                                     local_pos::CellIdx{NGrid}
                                    )::Vec{NBlockDims, Int32} where {NBlockDims, NGrid}
    return Vec{NBlockDims, Int32}() do a
        if (a > NGrid) || (size(block, a) == 1)
            one(Int32)
        else
            local_pos[a] + one(Int32)
        end
    end
end

dsl_format(f::ScaleOpsFactors) = string(
    "(", iter_join(f.factors, ", ")..., f.extend_last ? "..." : "", ")"
)
"Formats a set of colors, as a single color if it only has one"
dsl_format_scale_ops_set(s::RewriteRuleCell_Set) = dsl_format_rewrite_dest((length(s) == 1) ? first(s) : s)

function parse_markovjunior_scale_ops_factors(inputs::MacroParserInputs, loc, expr)::ScaleOpsFactors
    # `(2, 3)` is a tuple, while `(2...)` and `(2)` are just their single element.
    elements = Base.isexpr(expr, :tuple) ? collect(Any, expr.args) : Any[ expr ]
    extend_last = !isempty(elements) && Base.isexpr(elements[end], :...)
    if extend_last
        elements[end] = elements[end].args[1]
    end
    if isempty(elements) || !all(e -> (e isa Integer) && (e > 0), elements)
        raise_parse_error(loc, inputs,
                          "Expected positive integer factors like `(2, 3)` or `(2...)`; got `", expr, "`")
    end

    # Each named axis puts a floor on the grid's dimensionality.
    inputs.min_dims = max(inputs.min_dims, length(elements))

    return ScaleOpsFactors(convert(Vector{Int}, elements), extend_last)
end

"
Parses either a single rule, or a `begin ... end` block of them,
  using your lambda `(loc, rule_expr) -> rule`.
"
function parse_markovjunior_scale_ops_rules(parse_rule, loc, expr)::Vector{Any}
    rules = Any[ ]
    if Base.isexpr(expr, :block)
        parse_markovjunior_block(expr.args) do inner_loc, line
            push!(rules, parse_rule(something(inner_loc, loc), line))
        end
    else
        push!(rules, parse_rule(loc, expr))
    end
    return rules
end

"Checks that each axis of a rule's block either matches the scale factor or is extruded"
function check_scale_ops_block_size(inputs::MacroParserInputs, loc,
                                    block::Array, factors::ScaleOpsFactors)
    for axis in 1:ndims(block)
        factor = scale_ops_factor(factors, axis)
        if !in(size(block, axis), (1, factor))
            raise_parse_error(loc, inputs,
                              "Block is ", iter_join(size(block), "x")...,
                                ", but along the ", dsl_format_rewrite_axis(axis), " axis",
                                " it should be either 1 (to extrude) or ", factor, " (the scale factor)")
        end
    end
end

"
Picks a shared dimensionality for a scale Op's rule blocks, so the Op stays type-stable.
Each block's dimensionality ends at its last non-extruded axis,
  and puts a floor on the dimensionality of the grid itself.
"
function scale_ops_block_dims!(inputs::MacroParserInputs, blocks)::Int
    block_dims(block) = something(findlast(>(1), size(block)), 1)
    n_dims = maximum(block_dims, blocks; init=1)
    inputs.min_dims = max(inputs.min_dims, n_dims)
    return n_dims
end

"Converts a parsed block to the given element type and dimensionality (adding/dropping trailing size-1 axes)"
scale_ops_normalize_block(::Type{T}, block::Array, n_dims::Int) where {T} = reshape(
    collect(T, block),
    ntuple(a -> size(block, a), n_dims)
)


##########################################
#  @upscale

"A cell in an `@upscale` rule's output block"
const UpscaleCellDest = Union{UInt8, RewriteRuleCell_Set, RewriteRuleCell_Wildcard}

pick_upscale_value(cell::UInt8, source::UInt8, rng::PRNG) = cell
pick_upscale_value(cell::RewriteRuleCell_Set, source::UInt8, rng::PRNG) = rand(rng, cell)
pick_upscale_value(::RewriteRuleCell_Wildcard, source::UInt8, rng::PRNG) = source

"Upscales any of the source colors into a specific block of pixels"
struct UpscaleRule{NBlockDims}
    source::RewriteRuleCell_Set
    block::Array{UpscaleCellDest, NBlockDims}
end
Base.:(==)(a::UpscaleRule, b::UpscaleRule) = (a.source == b.source) && (a.block == b.block)

"Stretches each pixel into a block, optionally replacing specific colors with specific patterns"
struct MarkovOpUpscale{NBlockDims} <: AbstractMarkovOp
    factors::ScaleOpsFactors
    rules::Vector{UpscaleRule{NBlockDims}}
    # For each cell type, the index of the first rule it matches (or 0 for none).
    rule_by_cell::NTuple{N_CELL_TYPES, Int}
end
MarkovOpUpscale(factors::ScaleOpsFactors, rules::Vector{UpscaleRule{NBlockDims}}) where {NBlockDims} = MarkovOpUpscale{NBlockDims}(
    factors, rules,
    ntuple(Val(N_CELL_TYPES)) do i
        something(findfirst(r -> (UInt8(i - 1) in r.source), rules), 0)
    end
)

function markov_algo_run(op::MarkovOpUpscale{NBlockDims},
                         algo::MarkovAlgorithm, algo_state::AlgoState,
                         inherited_biases::Tuple{Vararg{AbstractMarkovBias}},
                         inherited_bias_states::Tuple,
                         grid::CellGrid{NGrid} = algo_state.grid
                        )::Bool where {NBlockDims, NGrid}
    VI = CellIdx{NGrid}
    scale_factor = scale_ops_factors(op.factors, Val(NGrid))
    new_size = vsize(grid) * scale_factor
    any_rule_applied = Ref(false)

    # No per-pixel ticks, because the new grid isn't visible externally until it's swapped in.
    markov_algo_new_grid(algo, algo_state, new_size,
                         inherited_biases, inherited_bias_states) do old_grid, new_grid
        for old_pos::VI in one(VI):vsize(old_grid)
            old_value::UInt8 = old_grid[old_pos]
            new_min = ((old_pos - one(Int32)) * scale_factor) + one(Int32)
            new_max = new_min + scale_factor - one(Int32)

            rule_i = op.rule_by_cell[old_value + 1]
            if rule_i == 0
                for new_pos::VI in new_min:new_max
                    new_grid[new_pos] = old_value
                end
            else
                any_rule_applied[] = true
                block = op.rules[rule_i].block
                for new_pos::VI in new_min:new_max
                    local_block_idx = scale_ops_block_idx(block, new_pos - new_min)
                    new_grid[new_pos] = pick_upscale_value(block[local_block_idx], old_value, algo_state.rng)
                end
            end
        end
    end

    markov_algo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY)
    return (scale_factor != one(VI)) || any_rule_applied[]
end


dsl_format(r::UpscaleRule) = string(
    dsl_format_scale_ops_set(r.source),
    " => ",
    dsl_format_rewrite_md_array(dsl_format_rewrite_dest, r.block)
)
dsl_format(op::MarkovOpUpscale) = string(
    "@upscale ", dsl_format(op.factors),
    if isempty(op.rules)
        ""
    else
        string(" begin\n    ", iter_join(dsl_format.(op.rules), "\n    ")..., "\nend")
    end
)

"Parses an `@upscale` rule into its source and (not yet normalized) block"
function parse_markovjunior_upscale_rule(inputs::MacroParserInputs, loc, expr,
                                         factors::ScaleOpsFactors
                                        )::Tuple{RewriteRuleCell_Set, Array{Any}}
    return with_parser_stacktrace(inputs, "Upscale rule") do
        if !@capture(expr, src_ => dest_)
            raise_parse_error(loc, inputs, "Expected a rule like `R => [ G B ]`; got `", expr, "`")
        end

        source_cell = with_parser_stacktrace(inputs, "Source") do
            parse_markovjunior_rewrite_cell(inputs, loc, src,
                                            allow_wildcard=false, allow_set=true,
                                            allow_list=false, allow_lookup=false)
        end
        source = (source_cell isa UInt8) ? RewriteRuleCell_Set(source_cell) : source_cell

        block = with_parser_stacktrace(inputs, "Destination") do
            parse_markovjunior_rewrite_block(inputs, loc, dest,
                                             allow_wildcard=true, allow_set=true,
                                             allow_list=false, allow_lookup=false)
        end
        check_scale_ops_block_size(inputs, loc, block, factors)

        return (source, block)
    end
end

function parse_markovjunior_op(::Val{Symbol("@upscale")}, inputs::MacroParserInputs,
                               loc, expr_args, expr_line)
    args = Vector{Any}()
    parse_markovjunior_block(expr_args) do inner_loc, line
        push!(args, line)
    end
    if !in(length(args), 1:2)
        raise_parse_error(loc, inputs, "Expected factors and then optionally some rules; got ", length(args), " arguments")
    end

    factors = with_parser_stacktrace(inputs, "Factors") do
        parse_markovjunior_scale_ops_factors(inputs, loc, args[1])
    end

    rule_data = if length(args) == 2
        parse_markovjunior_scale_ops_rules(loc, args[2]) do rule_loc, rule_expr
            parse_markovjunior_upscale_rule(inputs, rule_loc, rule_expr, factors)
        end
    else
        [ ]
    end

    n_block_dims = scale_ops_block_dims!(inputs, (block for (_, block) in rule_data))
    rules = UpscaleRule{n_block_dims}[
        UpscaleRule{n_block_dims}(source, scale_ops_normalize_block(UpscaleCellDest, block, n_block_dims))
          for (source, block) in rule_data
    ]
    return MarkovOpUpscale(factors, rules)
end


##########################################
#  @downscale

"""
How `@downscale` picks a color for blocks that no rule matches (including partial blocks at the grid's edges).
Ties are broken randomly.

* `random`: a random pixel from the block.
* `common`: the block's most common color.
* `rare`: the block's least common color (among the colors it has).
* `first`: the block's first pixel (its min corner).
"""
@bp_enum(DownscaleFallback,
    random,
    common,
    rare,
    first
)

"Picks the most (or least) common color in the given area of the grid, breaking ties randomly"
function pick_downscale_by_frequency(grid::CellGrid{N}, area_min::CellIdx{N}, area_max::CellIdx{N},
                                     pick_rarest::Bool, rng::PRNG)::UInt8 where {N}
    counts = zeros(MVector{N_CELL_TYPES, Int32})
    for pos::CellIdx{N} in area_min:area_max
        counts[grid[pos] + 1] += one(Int32)
    end

    best_value = zero(UInt8)
    best_count = zero(Int32)
    n_tied = 0
    for i in 1:N_CELL_TYPES
        n = counts[i]
        if iszero(n)
            continue
        elseif iszero(best_count) || (pick_rarest ? (n < best_count) : (n > best_count))
            (best_value, best_count, n_tied) = (UInt8(i - 1), n, 1)
        elseif n == best_count
            # Reservoir sampling, so every tied color is equally likely.
            n_tied += 1
            if rand(rng, 1:n_tied) == 1
                best_value = UInt8(i - 1)
            end
        end
    end
    return best_value
end

"Picks a color for an area of the grid that no rule matched"
function pick_downscale_fallback(fallback::E_DownscaleFallback,
                                 grid::CellGrid{N}, area_min::CellIdx{N}, area_max::CellIdx{N},
                                 rng::PRNG)::UInt8 where {N}
    if fallback == DownscaleFallback.first
        return grid[area_min]
    elseif fallback == DownscaleFallback.random
        return grid[area_min + CellIdx{N}(a -> rand(rng, zero(Int32):(area_max[a] - area_min[a])))]
    elseif fallback == DownscaleFallback.common
        return pick_downscale_by_frequency(grid, area_min, area_max, false, rng)
    elseif fallback == DownscaleFallback.rare
        return pick_downscale_by_frequency(grid, area_min, area_max, true, rng)
    else
        error("Unhandled: ", fallback)
    end
end
pick_downscale_fallback(always_colors::RewriteRuleCell_Set,
                        grid::CellGrid{N}, area_min::CellIdx{N}, area_max::CellIdx{N},
                        rng::PRNG) where {N} = rand(rng, always_colors)

"
Shrinks a specific block of pixels into one of the destination colors.
Each block cell is a set of allowed colors (a single color or wildcard are stored as sets too).
"
struct DownscaleRule_Block{NBlockDims}
    dest::RewriteRuleCell_Set
    block::Array{RewriteRuleCell_Set, NBlockDims}
end
Base.:(==)(a::DownscaleRule_Block, b::DownscaleRule_Block) = (a.dest == b.dest) && (a.block == b.block)

"Shrinks any block containing one of the given colors into one of the destination colors"
struct DownscaleRule_Has
    dest::RewriteRuleCell_Set
    colors::RewriteRuleCell_Set
end

const DownscaleRule{NBlockDims} = Union{DownscaleRule_Block{NBlockDims}, DownscaleRule_Has}

"The set of every color, which is what a wildcard `_` matches"
const DOWNSCALE_WILDCARD_SET = RewriteRuleCell_Set(collect(UInt8, 0:(N_CELL_TYPES - 1)))

downscale_source_set(cell::UInt8) = RewriteRuleCell_Set(cell)
downscale_source_set(cell::RewriteRuleCell_Set) = cell
downscale_source_set(::RewriteRuleCell_Wildcard) = DOWNSCALE_WILDCARD_SET

"Gets the set of colors in the given area of the grid"
function downscale_present_colors(grid::CellGrid{N}, area_min::CellIdx{N}, area_max::CellIdx{N}
                                 )::RewriteRuleCell_Set where {N}
    colors = RewriteRuleCell_Set()
    for pos::CellIdx{N} in area_min:area_max
        colors = push(colors, grid[pos])
    end
    return colors
end

"
Checks whether a rule matches the block of the grid from `block_min` to `block_max`,
  given the set of colors present in that block.
Partial blocks at the far edges of the grid can only match `has()` rules.
"
function downscale_rule_matches(rule::DownscaleRule_Block, grid::CellGrid{N},
                                block_min::CellIdx{N}, block_max::CellIdx{N},
                                scale_factor::CellIdx{N}, present_colors::RewriteRuleCell_Set
                               )::Bool where {N}
    if (block_max - block_min + one(Int32)) != scale_factor
        return false
    end
    for pos::CellIdx{N} in block_min:block_max
        local_block_idx = scale_ops_block_idx(rule.block, pos - block_min)
        if !in(grid[pos], rule.block[local_block_idx])
            return false
        end
    end
    return true
end
downscale_rule_matches(rule::DownscaleRule_Has, grid::CellGrid{N},
                       block_min::CellIdx{N}, block_max::CellIdx{N},
                       scale_factor::CellIdx{N}, present_colors::RewriteRuleCell_Set
                      ) where {N} = !iszero(rule.colors.bitfield & present_colors.bitfield)

"Shrinks each block of pixels into one pixel, optionally mapping specific patterns to specific colors"
struct MarkovOpDownscale{NBlockDims} <: AbstractMarkovOp
    factors::ScaleOpsFactors
    rules::Vector{DownscaleRule{NBlockDims}}
    # Either a fallback mode, or the colors from `always(...)`.
    fallback::Union{E_DownscaleFallback, RewriteRuleCell_Set}
end

function markov_algo_run(op::MarkovOpDownscale{NBlockDims},
                         algo::MarkovAlgorithm, algo_state::AlgoState,
                         inherited_biases::Tuple{Vararg{AbstractMarkovBias}},
                         inherited_bias_states::Tuple,
                         grid::CellGrid{NGrid} = algo_state.grid
                        )::Bool where {NBlockDims, NGrid}
    VI = CellIdx{NGrid}
    scale_factor = scale_ops_factors(op.factors, Val(NGrid))
    old_size = vsize(grid)
    # Round up, so leftover pixels at the far edges become partial blocks.
    new_size = VI(a -> cld(old_size[a], scale_factor[a]))
    any_rule_applied = Ref(false)
    needs_present_colors = any(r -> r isa DownscaleRule_Has, op.rules)

    # No per-pixel ticks, because the new grid isn't visible externally until it's swapped in.
    markov_algo_new_grid(algo, algo_state, new_size,
                         inherited_biases, inherited_bias_states) do old_grid, new_grid
        for new_pos::VI in one(VI):new_size
            old_min = ((new_pos - one(Int32)) * scale_factor) + one(Int32)
            old_max = min(old_min + scale_factor - one(Int32), old_size)

            present_colors = if needs_present_colors
                downscale_present_colors(old_grid, old_min, old_max)
            else
                RewriteRuleCell_Set()
            end

            matched_rule_i = findfirst(op.rules) do rule
                downscale_rule_matches(rule, old_grid, old_min, old_max, scale_factor, present_colors)
            end
            new_grid[new_pos] = if exists(matched_rule_i)
                any_rule_applied[] = true
                rand(algo_state.rng, op.rules[matched_rule_i].dest)
            else
                pick_downscale_fallback(op.fallback, old_grid, old_min, old_max, algo_state.rng)
            end
        end
    end

    markov_algo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY)
    return (new_size != old_size) || any_rule_applied[]
end


dsl_format_downscale_source(cell::RewriteRuleCell_Set) = if cell == DOWNSCALE_WILDCARD_SET
    "_"
else
    dsl_format_scale_ops_set(cell)
end
dsl_format(r::DownscaleRule_Block) = string(
    dsl_format_scale_ops_set(r.dest),
    " <= ",
    dsl_format_rewrite_md_array(dsl_format_downscale_source, r.block)
)
dsl_format(r::DownscaleRule_Has) = string(
    dsl_format_scale_ops_set(r.dest),
    " <= has(", dsl_format(r.colors), ")"
)
dsl_format_downscale_fallback(f::E_DownscaleFallback) = string(f)
dsl_format_downscale_fallback(f::RewriteRuleCell_Set) = "always($(dsl_format(f)))"
function dsl_format(op::MarkovOpDownscale)
    has_default_fallback = (op.fallback isa E_DownscaleFallback) && (op.fallback == DownscaleFallback.random)
    lines = String[ dsl_format.(op.rules)... ]
    if !has_default_fallback
        push!(lines, dsl_format_downscale_fallback(op.fallback))
    end

    return string(
        "@downscale ", dsl_format(op.factors),
        if isempty(lines)
            ""
        elseif isempty(op.rules)
            " $(lines[1])"
        else
            string(" begin\n    ", iter_join(lines, "\n    ")..., "\nend")
        end
    )
end

"Parses the colors inside a call like `has(wb)`, which are written without braces"
function parse_markovjunior_downscale_colors(inputs::MacroParserInputs, loc, expr,
                                             call_name::Symbol)::RewriteRuleCell_Set
    if !isa(expr, Symbol)
        raise_parse_error(loc, inputs,
                          "`", call_name, "()` takes colors without braces, like `", call_name, "(wb)`; got `",
                            expr, "`")
    end
    return parse_markovjunior_rewrite_cell(inputs, loc, Expr(:braces, expr),
                                           allow_wildcard=false, allow_set=true,
                                           allow_list=false, allow_lookup=false)
end

"Checks whether a line of a `@downscale` rule list is its fallback rather than a rule"
is_downscale_fallback_expr(expr) = (expr isa Symbol) || (Base.isexpr(expr, :call) && (expr.args[1] == :always))

function parse_markovjunior_downscale_fallback(inputs::MacroParserInputs, loc, expr
                                              )::Union{E_DownscaleFallback, RewriteRuleCell_Set}
    return with_parser_stacktrace(inputs, "Fallback") do
        if expr isa Symbol
            fallback = tryparse(E_DownscaleFallback, expr)
            if isnothing(fallback)
                raise_parse_error(loc, inputs,
                                  "Unknown fallback `", expr, "`; expected one of ",
                                    iter_join(DownscaleFallback.instances(), ", ")..., ", or `always(colors)`")
            end
            fallback
        elseif @capture(expr, always(colors_))
            parse_markovjunior_downscale_colors(inputs, loc, colors, :always)
        else
            raise_parse_error(loc, inputs, "Expected a fallback like `always(wb)`; got `", expr, "`")
        end
    end
end

"
Parses a `@downscale` rule.
`has()` rules are returned as-is, while block rules are returned as their destination and (not yet normalized) block.
"
function parse_markovjunior_downscale_rule(inputs::MacroParserInputs, loc, expr,
                                           factors::ScaleOpsFactors
                                          )::Union{DownscaleRule_Has, Tuple{RewriteRuleCell_Set, Array}}
    return with_parser_stacktrace(inputs, "Downscale rule") do
        if !@capture(expr, dest_ <= src_)
            raise_parse_error(loc, inputs, "Expected a rule like `R <= [ G B ]` or `R <= has(GB)`; got `", expr, "`")
        end

        dest_cell = with_parser_stacktrace(inputs, "Destination") do
            parse_markovjunior_rewrite_cell(inputs, loc, dest,
                                            allow_wildcard=false, allow_set=true,
                                            allow_list=false, allow_lookup=false)
        end
        dest_set = (dest_cell isa UInt8) ? RewriteRuleCell_Set(dest_cell) : dest_cell

        return with_parser_stacktrace(inputs, "Source") do
            if Base.isexpr(src, :call) && (src.args[1] == :has)
                if !@capture(src, has(colors_))
                    raise_parse_error(loc, inputs, "Expected `has(colors)` with one argument, like `has(wb)`; got `", src, "`")
                end
                DownscaleRule_Has(dest_set, parse_markovjunior_downscale_colors(inputs, loc, colors, :has))
            else
                block = parse_markovjunior_rewrite_block(inputs, loc, src,
                                                         allow_wildcard=true, allow_set=true,
                                                         allow_list=false, allow_lookup=false)
                check_scale_ops_block_size(inputs, loc, block, factors)
                (dest_set, map(downscale_source_set, block))
            end
        end
    end
end

function parse_markovjunior_op(::Val{Symbol("@downscale")}, inputs::MacroParserInputs,
                               loc, expr_args, expr_line)
    args = Vector{Any}()
    parse_markovjunior_block(expr_args) do inner_loc, line
        push!(args, line)
    end
    if !in(length(args), 1:2)
        raise_parse_error(loc, inputs,
                          "Expected factors and then optionally some rules (ending with an optional fallback); got ",
                            length(args), " arguments")
    end

    factors = with_parser_stacktrace(inputs, "Factors") do
        parse_markovjunior_scale_ops_factors(inputs, loc, args[1])
    end

    # The fallback is optional, and goes after all the rules.
    rule_data = [ ]
    fallback = DownscaleFallback.random
    if length(args) == 2
        lines = parse_markovjunior_scale_ops_rules(tuple, loc, args[2])
        for (i, (line_loc, line)) in enumerate(lines)
            if is_downscale_fallback_expr(line)
                if i != length(lines)
                    raise_parse_error(line_loc, inputs,
                                      "The fallback `", line, "` must be the last line, after all the rules")
                end
                fallback = parse_markovjunior_downscale_fallback(inputs, line_loc, line)
            else
                push!(rule_data, parse_markovjunior_downscale_rule(inputs, line_loc, line, factors))
            end
        end
    end

    block_rules_data = filter(r -> r isa Tuple, rule_data)
    n_block_dims = scale_ops_block_dims!(inputs, (block for (_, block) in block_rules_data))
    rules = DownscaleRule{n_block_dims}[
        if r isa DownscaleRule_Has
            r
        else
            (dest, block) = r
            DownscaleRule_Block{n_block_dims}(dest, scale_ops_normalize_block(RewriteRuleCell_Set, block, n_block_dims))
        end
          for r in rule_data
    ]
    return MarkovOpDownscale{n_block_dims}(factors, rules, fallback)
end
