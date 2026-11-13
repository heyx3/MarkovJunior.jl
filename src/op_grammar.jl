#NOTE: this file takes some things from op_rewrite.jl, op_scale.jl, and op_sequence.jl,
#   and should be included after them.


##################
#  Calculations

"
What a `@grammar` calculation can reference.

* `size`: the current region's resolution, as an array (so `size(y)`, `size_max()`, `size_max_axis()`, etc).
* `ndims`: the number of dimensions, i.e. `size_len()`.
* `grid_size`: the whole grid's size, as an array (so `grid_size(y)`, `grid_size_max()`, `grid_size_min_axis()`, etc).
* `depth`: how many expansions led to the current region (0 for the start symbol).
* The axis names `x`, `y`, `z`, and `w` (see `CALC_AXIS_CONSTANTS`).
* The algorithm's user functions (from `@calc`), looked up whenever they're called.
"
const GRAMMAR_CALC_ENV = CalcEnvironment(
    constants = CALC_AXIS_CONSTANTS,
    inputs = (:depth, :ndims),
    arrays = (:size, :grid_size),
    macros = true
)
"The calculations for a child's size also get `len`, the length of the parent region along the axis"
const GRAMMAR_CALC_SIZE_ENV = CalcEnvironment(
    constants = CALC_AXIS_CONSTANTS,
    inputs = (:depth, :ndims, :len),
    arrays = (:size, :grid_size),
    macros = true
)

"`place()` boxes get the same inputs as other `@grammar` calculations, plus the box inputs for each axis"
const GRAMMAR_CALC_BOX_ENV = CalcEnvironment(
    constants = CALC_AXIS_CONSTANTS,
    inputs = (:depth, :ndims, CALC_BOX_AXIS_INPUTS...),
    arrays = (:size, :grid_size),
    macros = true
)

"
Regions nested deeper than this become errors, painted with `GRAMMAR_ERROR_COLOR`,
  so rules that recurse forever (e.g. `A => B`, `B => A`) can't hang the run.
Regions at the same depth never overlap, so the expansion does at most `grid area * max depth` steps.
"
const GRAMMAR_MAX_DEPTH = 1000
const GRAMMAR_ERROR_COLOR = CELL_CODE_BY_CHAR['M']
"A leaf's rule index meaning the region went past `GRAMMAR_MAX_DEPTH`"
const GRAMMAR_ERROR_RULE = 0

"Logs a problem with a `@grammar` during a run, without stopping the algorithm"
grammar_log_error(msg...) = println(stderr, "@grammar error: ", msg...)


##################
#  Rules

"
A child region of a `split()`, `repeat()`, or `inset()`:
  the symbol it expands into, and (except for `inset()`) its size along the axis.
"
struct GrammarChild
    symbol::Int # Index into the grammar's symbols, or 0 for `_` (left untouched)
    size::CalcProgram
    is_weight::Bool # If false, the size is in pixels
end

"Splits the region along an axis into a sequence of children"
struct GrammarSplit
    axis::CalcProgram # Negative to lay the children out in reverse
    children::Vector{GrammarChild}
end
"Where the space left over by a `repeat()` goes"
struct GrammarLeftover
    # `:stretch` (the copies grow), `:spread` (the separators grow),
    #   or `:before`, `:after`, or `:around` (an extra child at the start, the end, or both).
    mode::Symbol
    symbol::Int # The extra child's symbol (or 0 for `_`), for the modes that add one
end
const GRAMMAR_LEFTOVER_MODES = (:stretch, :spread, :before, :after, :around)
const GRAMMAR_LEFTOVER_CHILD_MODES = (:before, :after, :around)

"Fills the region along an axis with as many copies of a child as will fit, optionally with separators"
struct GrammarRepeat
    axis::CalcProgram
    child::GrammarChild
    child_per_copy::Bool # If true, the size is evaluated again for every copy
    separator::Optional{GrammarChild} # Goes between each pair of copies
    separator_per_copy::Bool
    leftover::GrammarLeftover
end

"
One face of a box (an axis, negated for the face at the low end of it),
  and the symbol that should cover it instead of the default
"
struct GrammarFace
    face::CalcProgram
    symbol::Int
end

"An `inset()` amount along one axis, which may be different on the low and high sides"
struct GrammarInsetAmount
    min::CalcProgram
    max::Optional{CalcProgram} # If `nothing`, the `min` amount is used for both sides
end
"Splits the region into a shell (one box per face) and a core"
struct GrammarInset
    # Either one amount for every axis, or one amount per axis (missing axes use 0).
    amounts::Vector{GrammarInsetAmount}
    per_axis::Bool
    shell::Int
    core::Int
    faces::Vector{GrammarFace} # Faces covered by something other than `shell`
end

"Places a box inside the region; the rest is split into one box per face, like `inset()`"
struct GrammarPlace
    box::CalcBox # Relative to the region, in `GRAMMAR_CALC_BOX_ENV`
    inner::Int
    outer::Int
    faces::Vector{GrammarFace}
end

"
Expands one part of the region (the low half along each axis, or the high half for negative axes),
  then mirrors it onto the rest of the region.
"
struct GrammarMirror
    axes::Vector{CalcProgram}
    axes_is_list::Bool
    symbol::Int
    # If true, only the layout is mirrored, and the mirrored leaves run (and pick colors) independently.
    # Otherwise the source's final pixels are copied.
    layout::Bool
end
"Expands into a different symbol, covering the same region"
struct GrammarRename
    symbol::Int
end
"Fills the region with one color"
struct GrammarColor
    color::UInt8
end
"Leaves the region untouched"
struct GrammarKeep end
"Stretches a block of colors over the region (using nearest-neighbor sampling)"
struct GrammarBlock
    block::Array{UpscaleCellDest}
end
"Runs Ops within the region, like `@within_box`"
struct GrammarOps
    ops::Vector{AbstractMarkovOp}
end

const GrammarProduction = Union{GrammarSplit, GrammarRepeat, GrammarInset, GrammarPlace, GrammarMirror,
                                GrammarRename, GrammarColor, GrammarKeep, GrammarBlock, GrammarOps}

"One alternative for expanding a symbol"
struct GrammarRule
    guard::Optional{CalcProgram}
    weight::CalcProgram
    production::GrammarProduction
end

const GrammarStruct = Union{GrammarChild, GrammarLeftover, GrammarFace, GrammarInsetAmount,
                            GrammarProduction, GrammarRule}
Base.:(==)(a::GrammarStruct, b::GrammarStruct) =
    (typeof(a) == typeof(b)) && all(
        f -> getfield(a, f) == getfield(b, f),
        fieldnames(typeof(a))
    )


"
A split grammar: starting from the whole grid, recursively divides it into labeled boxes,
  then paints the leaf boxes and runs Ops within them.
"
struct MarkovOpGrammar <: AbstractMarkovOp
    inherit_biases::Bool
    symbols::Vector{Symbol} # The first one is always `GRAMMAR_ROOT_SYMBOL`
    rules::Vector{Vector{GrammarRule}} # The alternatives for each symbol
end
"The symbol whose box covers the whole grid"
const GRAMMAR_ROOT_SYMBOL = :ROOT


##################
#  Running

"
A region of the grid, labeled with a symbol.
A negative symbol marks the end of a `mirror()`'s source (see `GrammarMirrorState`),
  as everything pushed after it gets expanded before it's popped.
"
struct GrammarRegion{N}
    symbol::Int
    min::CellIdx{N}
    max::CellIdx{N}
    depth::Int
end
grammar_region_size(r::GrammarRegion) = r.max - r.min + one(Int32)

"A region which finished expanding, and the rule which turned it into a leaf"
struct GrammarLeafRegion{N}
    symbol::Int
    rule::Int # `GRAMMAR_ERROR_RULE` if the region went too deep, `GRAMMAR_MIRROR_RULE` for a pixel mirror
    min::CellIdx{N}
    max::CellIdx{N}
    # A bit per axis the leaf was mirrored along, by layout mirrors.
    # For a pixel mirror, it's instead the axes whose source is the high half.
    flips::UInt32
    mirror_axes::UInt32 # For a pixel mirror, a bit per axis it mirrors along
end
GrammarLeafRegion{N}(symbol, rule, min, max) where {N} = GrammarLeafRegion{N}(symbol, rule, min, max,
                                                                            zero(UInt32), zero(UInt32))
"A leaf's rule index meaning it copies one part of the region onto the rest, for a `mirror()`"
const GRAMMAR_MIRROR_RULE = -1

"A `mirror()` whose source part is being expanded"
struct GrammarMirrorState{N}
    symbol::Int
    min::CellIdx{N}
    max::CellIdx{N}
    axes::UInt32 # A bit per axis
    source_high::UInt32 # A bit per axis whose source is the high half
    layout::Bool
    first_leaf::Int # The index of the source's first leaf
end

"Temporary buffers for expanding a grammar"
struct GrammarScratch{N}
    pending::Vector{GrammarRegion{N}} # A stack, so the tree gets expanded depth-first
    leaves::Vector{GrammarLeafRegion{N}}
    weights::Vector{Float64}
    sizes::Vector{Float64}
    lengths::Vector{Int}
    symbols::Vector{Int} # The symbol of each child in `sizes`, for `repeat()`
    face_symbols::Vector{Int} # The symbol of each face (low then high side, for each axis)
    mirrors::Vector{GrammarMirrorState{N}}
end

"Everything needed to expand one region"
struct GrammarContext{N}
    op::MarkovOpGrammar
    rng::PRNG
    region::GrammarRegion{N}
    grid_size::CellIdx{N}
    scratch::GrammarScratch{N}
    calc_macros::OrderedDict{Symbol, CalcMacro}
end

"Evaluates a calculation for the current region; errors are logged and replaced with `fallback`"
eval_grammar_calc(p::CalcProgram, ctx::GrammarContext, fallback::CalcValue; len::Int = 0) = calc_evaluate(
    p, ctx.rng,
    e -> begin
        grammar_log_error(e.msg, "; using ", fallback)
        fallback
    end,
    ctx.calc_macros
    ;
    size = grammar_region_size(ctx.region).data,
    grid_size = ctx.grid_size.data,
    ndims = length(ctx.grid_size),
    depth = ctx.region.depth,
    len = len
)
"Evaluates a size or weight, treating negative or non-finite values as 0"
function eval_grammar_amount(p::CalcProgram, ctx::GrammarContext; len::Int = 0)::Float64
    v = convert(Float64, eval_grammar_calc(p, ctx, 0; len=len))
    return isfinite(v) ? max(0.0, v) : 0.0
end
"Evaluates an axis, returning it and whether it's reversed, or `nothing` if it's invalid"
function eval_grammar_axis(p::CalcProgram, ctx::GrammarContext{N})::Optional{Tuple{Int, Bool}} where {N}
    v = eval_grammar_calc(p, ctx, 0)
    if !isa(v, Int) || !in(abs(v), 1:N)
        grammar_log_error("Axis `", p, "` should be an Int from 1 to ", N, " (or negated); got `", v,
                            "`. Leaving the `", ctx.op.symbols[ctx.region.symbol], "` region untouched")
        return nothing
    end
    return (abs(v), v < 0)
end

function markov_algo_run(op::MarkovOpGrammar,
                         algo::MarkovAlgorithm, algo_state::AlgoState,
                         inherited_biases::Tuple{Vararg{AbstractMarkovBias}},
                         inherited_bias_states::Tuple,
                         grid::CellGrid{N} = algo_state.grid
                        )::OpChangedGrid where {N}
    alloc = algo_state.allocator
    scratch = GrammarScratch{N}(
        markov_allocator_acquire_array(alloc, (), GrammarRegion{N}),
        markov_allocator_acquire_array(alloc, (), GrammarLeafRegion{N}),
        markov_allocator_acquire_array(alloc, (), Float64),
        markov_allocator_acquire_array(alloc, (), Float64),
        markov_allocator_acquire_array(alloc, (), Int),
        markov_allocator_acquire_array(alloc, (), Int),
        markov_allocator_acquire_array(alloc, (), Int),
        markov_allocator_acquire_array(alloc, (), GrammarMirrorState{N})
    )
    try
        # Work out the whole layout first, then paint it, then run the Op leaves.
        expand_grammar!(op, grid, algo_state.rng, scratch, algo.calc_macros)
        @logic_logln("Grammar expanded into ", length(scratch.leaves), " leaves")

        result::OpChangedGrid = paint_grammar_leaves!(op, algo, algo_state, grid, scratch.leaves,
                                                      inherited_biases, inherited_bias_states)
        if result
            markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY)
        end

        # Pixel mirrors are in here too, in tree order, so they copy their source after its Ops ran.
        for leaf in scratch.leaves
            (leaf.rule == GRAMMAR_ERROR_RULE) && continue
            if leaf.rule == GRAMMAR_MIRROR_RULE
                if copy_grammar_mirror!(algo, algo_state, grid, leaf, inherited_biases, inherited_bias_states)
                    result = combine_op_results(result, true)
                    markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY)
                end
                continue
            end
            production = op.rules[leaf.symbol][leaf.rule].production
            if production isa GrammarOps
                leaf_result = run_grammar_ops_leaf(op, production, algo, algo_state, grid, leaf,
                                                   inherited_biases, inherited_bias_states)
                result = combine_op_results(result, leaf_result)
                markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY)
            end
        end

        markov_opalgo_tick(algo_state, STANDARD_END_OF_OP_TICK_PRIORITY + 1)
        return result
    finally
        markov_allocator_release_array(alloc, scratch.pending)
        markov_allocator_release_array(alloc, scratch.leaves)
        markov_allocator_release_array(alloc, scratch.weights)
        markov_allocator_release_array(alloc, scratch.sizes)
        markov_allocator_release_array(alloc, scratch.lengths)
        markov_allocator_release_array(alloc, scratch.symbols)
        markov_allocator_release_array(alloc, scratch.face_symbols)
        markov_allocator_release_array(alloc, scratch.mirrors)
    end
end

"Expands the whole grid into leaf regions, in depth-first order"
function expand_grammar!(op::MarkovOpGrammar, grid::CellGrid{N}, rng::PRNG,
                         scratch::GrammarScratch{N},
                         calc_macros::OrderedDict{Symbol, CalcMacro}) where {N}
    grid_size = convert(CellIdx{N}, vsize(grid))
    push!(scratch.pending, GrammarRegion{N}(1, one(CellIdx{N}), grid_size, 0))
    n_too_deep = 0
    too_deep_symbol = 0
    while !isempty(scratch.pending)
        region = pop!(scratch.pending)
        if region.symbol < 0
            finish_grammar_mirror!(scratch, scratch.mirrors[-region.symbol])
        elseif region.depth > GRAMMAR_MAX_DEPTH
            push!(scratch.leaves, GrammarLeafRegion{N}(region.symbol, GRAMMAR_ERROR_RULE, region.min, region.max))
            n_too_deep += 1
            too_deep_symbol = region.symbol
        else
            expand_grammar_region!(GrammarContext{N}(op, rng, region, grid_size, scratch, calc_macros))
        end
    end

    # Log once per run, as a runaway rule usually produces many regions.
    if n_too_deep > 0
        grammar_log_error(n_too_deep, " region(s) went past the maximum depth of ", GRAMMAR_MAX_DEPTH,
                          " (e.g. a `", op.symbols[too_deep_symbol], "`); do some rules recurse forever? ",
                          "They're painted '", CELL_TYPES[GRAMMAR_ERROR_COLOR + 1].char, "'")
    end
    return nothing
end

"
Picks one of the region's alternatives and applies it.
Alternatives whose guard fails are skipped, and so are alternatives which don't fit, without backtracking.
If none are left, the region is untouched.
"
function expand_grammar_region!(ctx::GrammarContext)
    rules = ctx.op.rules[ctx.region.symbol]
    weights = ctx.scratch.weights
    resize!(weights, length(rules))
    for (i, rule) in enumerate(rules)
        weights[i] = if exists(rule.guard) && !check_grammar_guard(rule.guard, ctx)
            0.0
        else
            eval_grammar_amount(rule.weight, ctx)
        end
    end

    while true
        total_weight = foldl(+, weights; init=0.0) # Unlike `sum`, always adds in the same order
        if !(total_weight > 0)
            @logic_logln("No alternative fits a `", ctx.op.symbols[ctx.region.symbol], "` region of size ",
                           grammar_region_size(ctx.region), "; leaving it untouched")
            return nothing
        end

        # Pick an alternative.
        t = rand(ctx.rng) * total_weight
        rule_i = findlast(>(0), weights)
        for (i, w) in enumerate(weights)
            t -= w
            if (w > 0) && (t < 0)
                rule_i = i
                break
            end
        end

        if try_expand_grammar!(rules[rule_i].production, rule_i, ctx)
            return nothing
        else
            weights[rule_i] = 0.0
        end
    end
end

function check_grammar_guard(guard::CalcProgram, ctx::GrammarContext)::Bool
    v = eval_grammar_calc(guard, ctx, false)
    if !isa(v, Bool)
        grammar_log_error("Guard `", guard, "` should give true or false; got `", v, "`. Treating it as false")
        return false
    end
    return v
end

"
Applies a production to the context's region, by pushing new regions and/or leaves.
Returns whether it fit the region.
"
try_expand_grammar!(p::GrammarProduction, rule_i::Int, ctx::GrammarContext)::Bool = error("Unhandled: ", typeof(p))

function try_expand_grammar!(p::Union{GrammarColor, GrammarBlock, GrammarOps}, rule_i::Int,
                             ctx::GrammarContext{N})::Bool where {N}
    r = ctx.region
    push!(ctx.scratch.leaves, GrammarLeafRegion{N}(r.symbol, rule_i, r.min, r.max))
    return true
end
try_expand_grammar!(::GrammarKeep, rule_i::Int, ctx::GrammarContext)::Bool = true

function try_expand_grammar!(p::GrammarRename, rule_i::Int, ctx::GrammarContext{N})::Bool where {N}
    r = ctx.region
    push!(ctx.scratch.pending, GrammarRegion{N}(p.symbol, r.min, r.max, r.depth + 1))
    return true
end

function try_expand_grammar!(p::GrammarSplit, rule_i::Int, ctx::GrammarContext)::Bool
    axis_data = eval_grammar_axis(p.axis, ctx)
    isnothing(axis_data) && return false
    (axis, reversed) = axis_data
    len = convert(Int, grammar_region_size(ctx.region)[axis])

    sizes = ctx.scratch.sizes
    resize!(sizes, length(p.children))
    for (i, child) in enumerate(p.children)
        sizes[i] = eval_grammar_amount(child.size, ctx; len=len)
    end
    if !layout_grammar_split!(sizes, p.children, len)
        return false
    end

    round_grammar_sizes!(ctx.scratch.lengths, sizes, len, ctx.rng)
    push_grammar_children!(ctx, axis, reversed, i -> p.children[i].symbol)
    return true
end

function try_expand_grammar!(p::GrammarRepeat, rule_i::Int, ctx::GrammarContext)::Bool
    axis_data = eval_grammar_axis(p.axis, ctx)
    isnothing(axis_data) && return false
    (axis, reversed) = axis_data
    len = convert(Int, grammar_region_size(ctx.region)[axis])

    # Sizes which aren't evaluated per-copy are evaluated once, up front.
    eval_size(child::GrammarChild, min_size::Float64) = begin
        v = eval_grammar_amount(child.size, ctx; len=len)
        if v < min_size
            grammar_log_error("`repeat()` size `", child.size, "` should be ",
                                iszero(min_size) ? "at least 0" : "positive", "; got ", v,
                                ". Leaving the `", ctx.op.symbols[ctx.region.symbol], "` region untouched")
            return nothing
        end
        return v
    end
    fixed_child_size = p.child_per_copy ? 0.0 : eval_size(p.child, eps(Float64))
    isnothing(fixed_child_size) && return false
    fixed_separator_size = (isnothing(p.separator) || p.separator_per_copy) ? 0.0 : eval_size(p.separator, 0.0)
    isnothing(fixed_separator_size) && return false

    # Add copies (and separators between them) until the next one doesn't fit.
    # Every copy should get at least one pixel, so there are never more copies than that.
    sizes = empty!(ctx.scratch.sizes)
    symbols = empty!(ctx.scratch.symbols)
    used = 0.0
    n_copies = 0
    while n_copies < len
        separator_size = 0.0
        if (n_copies > 0) && exists(p.separator)
            separator_size = p.separator_per_copy ? eval_size(p.separator, 0.0) : fixed_separator_size
            isnothing(separator_size) && return false
        end
        child_size = p.child_per_copy ? eval_size(p.child, eps(Float64)) : fixed_child_size
        isnothing(child_size) && return false

        # Allow for a little floating-point error, e.g. a size of `len/3`.
        if (used + separator_size + child_size) > (len + 1e-6)
            break
        end
        if (n_copies > 0) && exists(p.separator)
            push!(sizes, separator_size)
            push!(symbols, p.separator.symbol)
        end
        push!(sizes, child_size)
        push!(symbols, p.child.symbol)
        used += separator_size + child_size
        n_copies += 1
    end
    if n_copies < 1
        return false
    end

    # Hand out the leftover space.
    leftover = max(0.0, len - used)
    mode = p.leftover.mode
    if (mode == :spread) && (n_copies < 2)
        mode = :stretch # Nothing to spread
    end
    if mode in (:stretch, :spread)
        # Copies are at the odd indices, separators at the even ones.
        indices = (mode == :stretch) ? (1:2:length(sizes)) : (2:2:length(sizes))
        if isnothing(p.separator)
            indices = 1:length(sizes)
        end
        total = foldl((t, i) -> t + sizes[i], indices; init=0.0) # Not `sum()`, which may use SIMD
        for i in indices
            sizes[i] += leftover * ((total > 0) ? (sizes[i] / total) : (1 / length(indices)))
        end
    elseif mode == :before
        pushfirst!(sizes, leftover)
        pushfirst!(symbols, p.leftover.symbol)
    elseif mode == :after
        push!(sizes, leftover)
        push!(symbols, p.leftover.symbol)
    elseif mode == :around
        pushfirst!(sizes, leftover / 2)
        pushfirst!(symbols, p.leftover.symbol)
        push!(sizes, leftover / 2)
        push!(symbols, p.leftover.symbol)
    else
        error("Unhandled: ", mode)
    end

    round_grammar_sizes!(ctx.scratch.lengths, sizes, len, ctx.rng)
    push_grammar_children!(ctx, axis, reversed, i -> symbols[i])
    return true
end

function try_expand_grammar!(p::GrammarInset, rule_i::Int, ctx::GrammarContext{N})::Bool where {N}
    r = ctx.region
    region_size = grammar_region_size(r)

    # Get the inset along each axis, on each side.
    eval_amount(c::CalcProgram) = convert(Int32, round(Int, eval_grammar_amount(c, ctx)))
    single_amount = p.per_axis ? zero(Int32) : eval_amount(p.amounts[1].min)
    min_amounts = CellIdx{N}() do axis
        if !p.per_axis
            single_amount
        elseif axis <= length(p.amounts)
            eval_amount(p.amounts[axis].min)
        else
            zero(Int32)
        end
    end
    max_amounts = CellIdx{N}() do axis
        if !p.per_axis
            single_amount
        elseif axis <= length(p.amounts)
            # A single amount applies to both sides; it's only evaluated once, in case it's random.
            a = p.amounts[axis]
            isnothing(a.max) ? min_amounts[axis] : eval_amount(a.max)
        else
            zero(Int32)
        end
    end
    # The core can't be empty.
    if any(a -> (min_amounts[a] + max_amounts[a]) >= region_size[a], 1:N)
        return false
    end

    push_grammar_shell!(ctx, r.min + min_amounts, r.max - max_amounts, p.shell, p.faces, p.core)
    return true
end

function try_expand_grammar!(p::GrammarPlace, rule_i::Int, ctx::GrammarContext{N})::Bool where {N}
    box = eval_grammar_box(p.box, ctx)
    isnothing(box) && return false
    push_grammar_shell!(ctx, box[1], box[2], p.outer, p.faces, p.inner)
    return true
end

"
Finds the pixels a box covers within the current region (see `eval_calc_box()`),
  or `nothing` if it's empty
"
function eval_grammar_box(b::CalcBox, ctx::GrammarContext{N})::Optional{Tuple{CellIdx{N}, CellIdx{N}}} where {N}
    r = ctx.region
    region_size = grammar_region_size(r)
    local_box = eval_calc_box(grammar_log_error, b, region_size, ctx.calc_macros) do calc, axis_inputs
        calc_evaluate(
            calc, ctx.rng,
            e -> begin
                grammar_log_error(e.msg, "; using 0")
                0
            end,
            ctx.calc_macros
            ;
            size = region_size.data,
            grid_size = ctx.grid_size.data,
            ndims = N,
            depth = r.depth,
            axis_inputs...
        )
    end
    isnothing(local_box) && return nothing
    return (r.min + min_inclusive(local_box) - one(Int32), r.min + max_inclusive(local_box) - one(Int32))
end

grammar_has_axis(mask::UInt32, axis::Integer) = !iszero(mask & (one(UInt32) << (axis - 1)))

"
Gets the source part of a mirrored region: along each mirrored axis, the low half
  (or the high half if that axis is in `source_high`), rounded up so it includes any middle slice.
"
function grammar_mirror_source(min::CellIdx{N}, max::CellIdx{N}, axes::UInt32, source_high::UInt32
                              )::Tuple{CellIdx{N}, CellIdx{N}} where {N}
    half(a) = convert(Int32, cld(max[a] - min[a] + 1, 2))
    return (
        CellIdx{N}(a -> (grammar_has_axis(axes, a) && grammar_has_axis(source_high, a)) ?
                          (max[a] - half(a) + one(Int32)) :
                          min[a]),
        CellIdx{N}(a -> (grammar_has_axis(axes, a) && !grammar_has_axis(source_high, a)) ?
                          (min[a] + half(a) - one(Int32)) :
                          max[a])
    )
end

function try_expand_grammar!(p::GrammarMirror, rule_i::Int, ctx::GrammarContext{N})::Bool where {N}
    r = ctx.region
    axes = zero(UInt32)
    source_high = zero(UInt32)
    for axis_calc in p.axes
        axis_data = eval_grammar_axis(axis_calc, ctx)
        isnothing(axis_data) && return false
        (axis, reversed) = axis_data
        bit = one(UInt32) << (axis - 1)
        axes |= bit
        source_high = reversed ? (source_high | bit) : (source_high & ~bit)
    end

    # Expand the source, followed by a marker which finishes the mirror once the source is fully expanded.
    scratch = ctx.scratch
    push!(scratch.mirrors, GrammarMirrorState{N}(r.symbol, r.min, r.max, axes, source_high, p.layout,
                                                 length(scratch.leaves) + 1))
    push!(scratch.pending, GrammarRegion{N}(-length(scratch.mirrors), r.min, r.max, r.depth))
    if p.symbol != 0
        (source_min, source_max) = grammar_mirror_source(r.min, r.max, axes, source_high)
        push!(scratch.pending, GrammarRegion{N}(p.symbol, source_min, source_max, r.depth + 1))
    end
    return true
end

"Called once a mirror's source is fully expanded, to mirror it onto the rest of the region"
function finish_grammar_mirror!(scratch::GrammarScratch{N}, m::GrammarMirrorState{N}) where {N}
    # Pixel mirrors copy the source once its Ops have run, so they're a leaf after the source's leaves.
    if !m.layout
        push!(scratch.leaves, GrammarLeafRegion{N}(m.symbol, GRAMMAR_MIRROR_RULE, m.min, m.max,
                                                   m.source_high, m.axes))
        return nothing
    end

    # Layout mirrors copy the source's leaves onto each other part of the region,
    #   where each part is the reflection along a different subset of the axes.
    (source_min, source_max) = grammar_mirror_source(m.min, m.max, m.axes, m.source_high)
    last_leaf = length(scratch.leaves)
    reflected_axes = m.axes
    while !iszero(reflected_axes)
        for leaf_i in m.first_leaf:last_leaf
            leaf = scratch.leaves[leaf_i]
            # Reflect the leaf, then cut off any part that overlaps the source (e.g. a middle slice).
            reflect(a, p) = m.min[a] + m.max[a] - p
            new_min = CellIdx{N}() do a
                if !grammar_has_axis(reflected_axes, a)
                    leaf.min[a]
                elseif grammar_has_axis(m.source_high, a)
                    reflect(a, leaf.max[a])
                else
                    max(reflect(a, leaf.max[a]), source_max[a] + one(Int32))
                end
            end
            new_max = CellIdx{N}() do a
                if !grammar_has_axis(reflected_axes, a)
                    leaf.max[a]
                elseif grammar_has_axis(m.source_high, a)
                    min(reflect(a, leaf.min[a]), source_min[a] - one(Int32))
                else
                    reflect(a, leaf.min[a])
                end
            end
            if all(a -> new_min[a] <= new_max[a], 1:N)
                push!(scratch.leaves, GrammarLeafRegion{N}(leaf.symbol, leaf.rule, new_min, new_max,
                                                           leaf.flips ⊻ reflected_axes, leaf.mirror_axes))
            end
        end
        # Move on to the next subset of the axes.
        reflected_axes = (reflected_axes - one(UInt32)) & m.axes
    end
    return nothing
end

"Copies a pixel mirror's source onto the rest of its region, returning whether any pixels changed"
function copy_grammar_mirror!(algo::MarkovAlgorithm, algo_state::AlgoState,
                              grid::CellGrid{N}, leaf::GrammarLeafRegion{N},
                              biases::Tuple, bias_states::Tuple)::Bool where {N}
    (source_min, source_max) = grammar_mirror_source(leaf.min, leaf.max, leaf.mirror_axes, leaf.flips)
    box_ranges = ntuple(a -> leaf.min[a]:leaf.max[a], Val(N))
    previous_pixels = if isempty(biases)
        nothing
    else
        previous::Array{UInt8, N} = markov_allocator_acquire_array(algo_state.allocator, length.(box_ranges), UInt8)
        previous .= @view grid[box_ranges...]
        previous
    end
    try
        changed = false
        for pos::CellIdx{N} in leaf.min:leaf.max
            # Pixels outside the source get reflected into it.
            src = CellIdx{N}() do a
                if (pos[a] < source_min[a]) || (pos[a] > source_max[a])
                    leaf.min[a] + leaf.max[a] - pos[a]
                else
                    pos[a]
                end
            end
            if grid[pos] != grid[src]
                grid[pos] = grid[src]
                changed = true
            end
        end
        if changed && exists(previous_pixels)
            markov_opbias_update_all!(biases, bias_states, algo, algo_state,
                                      BoxI{N}(min=leaf.min, max=leaf.max), previous_pixels)
        end
        return changed
    finally
        exists(previous_pixels) && markov_allocator_release_array(algo_state.allocator, previous_pixels)
    end
end

"
Pushes a core box, and a shell around it split into one box per face.
The edges and corners of the shell go to the faces along the lowest axis.
Faces which are empty get skipped.
"
function push_grammar_shell!(ctx::GrammarContext{N}, core_min::CellIdx{N}, core_max::CellIdx{N},
                             shell::Int, faces::Vector{GrammarFace}, core::Int) where {N}
    r = ctx.region

    # Get each face's symbol: for axis `a`, the low face is at `2a - 1` and the high face at `2a`.
    face_symbols = ctx.scratch.face_symbols
    resize!(face_symbols, 2 * N)
    fill!(face_symbols, shell)
    for face in faces
        axis_data = eval_grammar_axis(face.face, ctx)
        isnothing(axis_data) && continue
        (axis, is_min_side) = axis_data
        face_symbols[(2 * axis) - (is_min_side ? 1 : 0)] = face.symbol
    end

    # Push in reverse order, so that the faces get expanded before the core.
    if core != 0
        push!(ctx.scratch.pending, GrammarRegion{N}(core, core_min, core_max, r.depth + 1))
    end
    face_range(axis, face_axis) = if axis < face_axis
        (core_min[axis], core_max[axis])
    else
        (r.min[axis], r.max[axis])
    end
    for face_axis in N:-1:1
        for is_max_side in (true, false)
            symbol = face_symbols[(2 * face_axis) - (is_max_side ? 0 : 1)]
            (symbol == 0) && continue
            face_min = CellIdx{N}(a -> if a == face_axis
                is_max_side ? (core_max[a] + one(Int32)) : r.min[a]
            else
                face_range(a, face_axis)[1]
            end)
            face_max = CellIdx{N}(a -> if a == face_axis
                is_max_side ? r.max[a] : (core_min[a] - one(Int32))
            else
                face_range(a, face_axis)[2]
            end)
            if face_min[face_axis] <= face_max[face_axis]
                push!(ctx.scratch.pending, GrammarRegion{N}(symbol, face_min, face_max, r.depth + 1))
            end
        end
    end
    return nothing
end

"
Turns the evaluated sizes of a split's children (pixels or weights)
  into pixel sizes which add up to `total`, in place.
Weighted children share the space left over by the pixel-sized ones;
  if there are no weighted children, the pixel sizes are stretched to fill the space.
Returns false if the pixel sizes are too big to fit.
"
function layout_grammar_split!(sizes::Vector{Float64}, children::Vector{GrammarChild}, total::Int)::Bool
    total_pixels = 0.0
    total_weight = 0.0
    n_weighted = 0
    for (size, child) in zip(sizes, children)
        if child.is_weight
            total_weight += size
            n_weighted += 1
        else
            total_pixels += size
        end
    end

    leftover = total - total_pixels
    if leftover < -1e-6
        return false
    end
    leftover = max(0.0, leftover)

    if n_weighted > 0
        for (i, child) in enumerate(children)
            if child.is_weight
                sizes[i] = leftover * ((total_weight > 0) ? (sizes[i] / total_weight) : (1 / n_weighted))
            end
        end
    elseif total_pixels > 0
        sizes .*= total / total_pixels
    else
        return false
    end
    return true
end

"
Rounds the sizes into integers which add up to `total`,
  using the largest-remainder method with random tie-breaks.
The sizes are modified in the process.
"
function round_grammar_sizes!(out::Vector{Int}, sizes::Vector{Float64}, total::Int, rng::PRNG)
    resize!(out, length(sizes))
    for i in 1:length(sizes)
        out[i] = floor(Int, sizes[i])
        sizes[i] -= out[i]
    end

    for _ in 1:(total - sum(out))
        best_i = 0
        best_remainder = -Inf
        n_tied = 0
        for (i, remainder) in enumerate(sizes)
            if remainder > best_remainder
                (best_i, best_remainder, n_tied) = (i, remainder, 1)
            elseif remainder == best_remainder
                # Reservoir sampling, so every tied size is equally likely.
                n_tied += 1
                if rand(rng, 1:n_tied) == 1
                    best_i = i
                end
            end
        end
        out[best_i] += 1
        sizes[best_i] = -Inf
    end
    return nothing
end

"
Pushes children laid out along an axis, using the lengths in `ctx.scratch.lengths`.
They're pushed in reverse, so they get expanded in order.
Empty children, and children which are `_`, are skipped.
"
function push_grammar_children!(ctx::GrammarContext{N}, axis::Int, reversed::Bool, get_symbol) where {N}
    r = ctx.region
    lengths = ctx.scratch.lengths
    # Walk backwards through the children, starting from the end of the axis where the last child goes.
    pos::Int32 = reversed ? r.min[axis] : r.max[axis]
    for i in length(lengths):-1:1
        len = convert(Int32, lengths[i])
        (child_start, child_end) = if reversed
            (pos, pos + len - one(Int32))
        else
            (pos - len + one(Int32), pos)
        end
        pos = reversed ? (child_end + one(Int32)) : (child_start - one(Int32))

        symbol = get_symbol(i)
        if (len > 0) && (symbol != 0)
            child_min = CellIdx{N}(a -> (a == axis) ? child_start : r.min[a])
            child_max = CellIdx{N}(a -> (a == axis) ? child_end : r.max[a])
            push!(ctx.scratch.pending, GrammarRegion{N}(symbol, child_min, child_max, r.depth + 1))
        end
    end
    return nothing
end

"Paints every color and block leaf, then updates the biases. Returns whether any pixels changed."
function paint_grammar_leaves!(op::MarkovOpGrammar, algo::MarkovAlgorithm, algo_state::AlgoState,
                               grid::CellGrid{N}, leaves::Vector{GrammarLeafRegion{N}},
                               biases::Tuple, bias_states::Tuple)::Bool where {N}
    previous_pixels = if isempty(biases)
        nothing
    else
        a::Array{UInt8, N} = markov_allocator_acquire_array(algo_state.allocator, size(grid), UInt8)
        copyto!(a, grid)
        a
    end
    try
        changed = false
        for leaf in leaves
            if leaf.rule == GRAMMAR_ERROR_RULE
                for pos::CellIdx{N} in leaf.min:leaf.max
                    if grid[pos] != GRAMMAR_ERROR_COLOR
                        grid[pos] = GRAMMAR_ERROR_COLOR
                        changed = true
                    end
                end
                continue
            end
            (leaf.rule == GRAMMAR_MIRROR_RULE) && continue # Done after the Op leaves run
            production = op.rules[leaf.symbol][leaf.rule].production
            if production isa GrammarColor
                for pos::CellIdx{N} in leaf.min:leaf.max
                    if grid[pos] != production.color
                        grid[pos] = production.color
                        changed = true
                    end
                end
            elseif production isa GrammarBlock
                changed |= paint_grammar_block!(grid, leaf.min, leaf.max, leaf.flips, production.block,
                                                algo_state.rng)
            end
        end

        if changed && exists(previous_pixels)
            markov_opbias_update_all!(biases, bias_states, algo, algo_state,
                                      BoxI{N}(min=one(CellIdx{N}), max=convert(CellIdx{N}, vsize(grid))),
                                      previous_pixels)
        end
        return changed
    finally
        exists(previous_pixels) && markov_allocator_release_array(algo_state.allocator, previous_pixels)
    end
end

"
Stretches the block over the given area of the grid, returning whether any pixels changed.
`flips` has a bit for each axis the block should be reflected along.
"
function paint_grammar_block!(grid::CellGrid{N}, area_min::CellIdx{N}, area_max::CellIdx{N}, flips::UInt32,
                              block::Array{UpscaleCellDest, NBlock}, rng::PRNG
                             )::Bool where {N, NBlock}
    area_size = area_max - area_min + one(Int32)
    changed = false
    for pos::CellIdx{N} in area_min:area_max
        local_pos = CellIdx{N}() do a
            grammar_has_axis(flips, a) ? (area_max[a] - pos[a]) : (pos[a] - area_min[a])
        end
        # Block axes past the grid's dimensions have size 1, which is enforced while parsing.
        block_idx = Vec{NBlock, Int32}() do a
            if a > N
                one(Int32)
            else
                convert(Int32, 1 + div(convert(Int, local_pos[a]) * size(block, a), convert(Int, area_size[a])))
            end
        end
        old_value = grid[pos]
        new_value = pick_upscale_value(block[block_idx], old_value, rng)
        if new_value != old_value
            grid[pos] = new_value
            changed = true
        end
    end
    return changed
end

"Runs an Op leaf within its region"
function run_grammar_ops_leaf(op::MarkovOpGrammar, leaf_ops::GrammarOps,
                              algo::MarkovAlgorithm, algo_state::AlgoState,
                              grid::CellGrid{N}, leaf::GrammarLeafRegion{N},
                              biases::Tuple, bias_states::Tuple)::OpChangedGrid where {N}
    @logic_logln("Running ", length(leaf_ops.ops), " ops for `", op.symbols[leaf.symbol], "` within ",
                   leaf.min, " to ", leaf.max)
    sub_grid = view(grid, ntuple(i -> leaf.min[i]:leaf.max[i], Val(N))...)
    return markov_opalgo_on_subgrid(algo, algo_state,
                                    sub_grid, BoxI{N}(min=leaf.min, max=leaf.max),
                                    op.inherit_biases, biases, bias_states
                                   ) do inner_biases, inner_bias_states
        # The grid can't be reallocated in here, so the bias states never get replaced.
        inner_result::OpChangedGrid = nothing
        for inner_op in leaf_ops.ops
            inner_result = combine_op_results(inner_result, markov_algo_run(
                inner_op, algo, algo_state,
                inner_biases, inner_bias_states
            ))
        end
        return inner_result
    end
end

# The inner Ops are rejected at parse-time if they reallocate.
markov_op_reallocates_grid(op::MarkovOpGrammar) = any(op.rules) do symbol_rules
    any(symbol_rules) do rule
        (rule.production isa GrammarOps) && any(markov_op_reallocates_grid, rule.production.ops)
    end
end


##################
#  DSL

dsl_format_grammar_symbol(op::MarkovOpGrammar, i::Int) = iszero(i) ? "_" : string(op.symbols[i])
dsl_format_grammar_child(op::MarkovOpGrammar, c::GrammarChild; per_copy::Bool = false) = string(
    dsl_format_grammar_symbol(op, c.symbol),
    if c.is_weight && !per_copy && is_default_grammar_weight(c.size)
        "" # A bare child has a weight of 1
    else
        "($(c.is_weight ? "w=" : "")$(per_copy ? "each" : "")(" * string(c.size) * "))"
    end
)
dsl_format_grammar_faces(op::MarkovOpGrammar, faces::Vector{GrammarFace}) = string(
    (", ($(f.face)) => $(dsl_format_grammar_symbol(op, f.symbol))" for f in faces)...
)
"Formats one value per axis (as a list), or one value for all of them"
dsl_format_grammar_per_axis(values::Vector, is_list::Bool) = if is_list
    string("[", iter_join(Iterators.map(dsl_format_grammar_value, values), ", ")..., "]")
else
    dsl_format_grammar_value(values[1])
end
dsl_format_grammar_value(c::CalcProgram) = "($c)"
dsl_format_grammar_value(a::GrammarInsetAmount) = if isnothing(a.max)
    "($(a.min))"
else
    "[($(a.min)), ($(a.max))]"
end
dsl_format_grammar(op::MarkovOpGrammar, p::GrammarSplit) = string(
    "split((", p.axis, "), ",
    iter_join(Iterators.map(c -> dsl_format_grammar_child(op, c), p.children), ", ")...,
    ")"
)
dsl_format_grammar(op::MarkovOpGrammar, p::GrammarRepeat) = string(
    "repeat((", p.axis, "), ", dsl_format_grammar_child(op, p.child; per_copy=p.child_per_copy),
    if isnothing(p.separator)
        ""
    else
        string(", ", dsl_format_grammar_child(op, p.separator; per_copy=p.separator_per_copy))
    end,
    if p.leftover.mode == :stretch
        ""
    elseif p.leftover.mode in GRAMMAR_LEFTOVER_CHILD_MODES
        "; leftover=$(p.leftover.mode)($(dsl_format_grammar_symbol(op, p.leftover.symbol)))"
    else
        "; leftover=$(p.leftover.mode)"
    end,
    ")"
)
dsl_format_grammar(op::MarkovOpGrammar, p::GrammarInset) = string(
    "inset(",
    dsl_format_grammar_per_axis(p.amounts, p.per_axis),
    ", ", dsl_format_grammar_symbol(op, p.shell),
    ", ", dsl_format_grammar_symbol(op, p.core),
    dsl_format_grammar_faces(op, p.faces),
    ")"
)
dsl_format_grammar(op::MarkovOpGrammar, p::GrammarPlace) = string(
    "place(",
    dsl_format(p.box),
    ", ", dsl_format_grammar_symbol(op, p.inner),
    ", ", dsl_format_grammar_symbol(op, p.outer),
    dsl_format_grammar_faces(op, p.faces),
    ")"
)
dsl_format_grammar(op::MarkovOpGrammar, p::GrammarMirror) = string(
    "mirror(",
    dsl_format_grammar_per_axis(p.axes, p.axes_is_list),
    ", ", dsl_format_grammar_symbol(op, p.symbol),
    p.layout ? ", layout" : "",
    ")"
)
dsl_format_grammar(op::MarkovOpGrammar, p::GrammarRename) = dsl_format_grammar_symbol(op, p.symbol)
dsl_format_grammar(op::MarkovOpGrammar, p::GrammarColor) = "'$(dsl_format(p.color))'"
dsl_format_grammar(op::MarkovOpGrammar, p::GrammarKeep) = "_"
dsl_format_grammar(op::MarkovOpGrammar, p::GrammarBlock) =
    dsl_format_rewrite_md_array(dsl_format_rewrite_dest, p.block)
dsl_format_grammar(op::MarkovOpGrammar, p::GrammarOps) = string(
    "begin\n        ",
    iter_join(Iterators.map(dsl_format, p.ops), "\n        ")...,
    "\n    end"
)

is_default_grammar_weight(w::CalcProgram) = calc_is_constant(w) && (calc_constant_value(w) == 1)

dsl_format_grammar(op::MarkovOpGrammar, symbol_i::Int, rule::GrammarRule) = string(
    op.symbols[symbol_i],
    isnothing(rule.guard) ? "" : "(($(rule.guard)))",
    " => ",
    dsl_format_grammar(op, rule.production),
    is_default_grammar_weight(rule.weight) ? "" : " *($(rule.weight))"
)

dsl_format(op::MarkovOpGrammar) = string(
    "@grammar ",
    op.inherit_biases ? "inherit_biases " : "",
    "begin\n    ",
    iter_join(
        (dsl_format_grammar(op, symbol_i, rule)
           for (symbol_i, symbol_rules) in enumerate(op.rules)
             for rule in symbol_rules),
        "\n    "
    )...,
    "\nend"
)


"
Separates a rule's weight (a trailing `*N` or `/N`) from its production.
Productions are never arithmetic, so any top-level multiplication or division is a weight.
"
function peel_grammar_weight(rhs)::Tuple{Any, Any}
    if Base.isexpr(rhs, :call) && (rhs.args[1] == :*) && (length(rhs.args) >= 3)
        factors = rhs.args[3:end]
        return (rhs.args[2], (length(factors) == 1) ? factors[1] : Expr(:call, :*, factors...))
    elseif Base.isexpr(rhs, :call) && (rhs.args[1] == :/) && (length(rhs.args) == 3)
        return (rhs.args[2], Expr(:call, :/, 1, rhs.args[3]))
    else
        return (rhs, 1)
    end
end

"Finds the index of a symbol, or 0 for `_`"
function parse_markovjunior_grammar_symbol(inputs::MacroParserInputs, loc, expr,
                                           symbol_lookup::Dict{Symbol, Int})::Int
    if expr == :_
        return 0
    elseif !isa(expr, Symbol)
        raise_parse_error(loc, inputs, "Expected a symbol (or `_`); got `", expr, "`")
    elseif !haskey(symbol_lookup, expr)
        raise_parse_error(loc, inputs, "Unknown symbol `", expr, "`; every symbol needs at least one rule")
    else
        return symbol_lookup[expr]
    end
end

"Parses a child like `Room` or `Room()` (a weight of 1), `Room(3)`, or `Room(w=2)`"
function parse_markovjunior_grammar_child(inputs::MacroParserInputs, loc, expr,
                                          symbol_lookup::Dict{Symbol, Int}
                                          ;
                                          allow_weight::Bool)::GrammarChild
    (name, size_expr) = if expr isa Symbol
        (expr, nothing)
    elseif Base.isexpr(expr, :call) && (length(expr.args) == 1)
        (expr.args[1], nothing)
    elseif Base.isexpr(expr, :call) && (length(expr.args) == 2)
        (expr.args[1], expr.args[2])
    else
        raise_parse_error(loc, inputs, "Expected a child like `Room`, `Room(3)`, or `Room(w=2)`; got `", expr, "`")
    end
    symbol = parse_markovjunior_grammar_symbol(inputs, loc, name, symbol_lookup)
    if Base.isexpr(size_expr, :call) && (size_expr.args[1] == :each)
        raise_parse_error(loc, inputs, "Only the children of `repeat()` can use `each(...)`; got `", expr, "`")
    end

    if isnothing(size_expr)
        if !allow_weight
            raise_parse_error(loc, inputs, "Child `", expr, "` needs a size in pixels, like `", name, "(3)`")
        end
        return GrammarChild(symbol, calc_compile(1, GRAMMAR_CALC_SIZE_ENV), true)
    elseif Base.isexpr(size_expr, :kw)
        if size_expr.args[1] != :w
            raise_parse_error(loc, inputs, "A child's weight is written `w=...`; got `", expr, "`")
        elseif !allow_weight
            raise_parse_error(loc, inputs, "Child `", expr, "` needs a size in pixels, not a weight")
        end
        return GrammarChild(symbol,
                            parse_markovjunior_calc(inputs, loc, size_expr.args[2], GRAMMAR_CALC_SIZE_ENV,
                                                    "Weight"),
                            true)
    else
        return GrammarChild(symbol,
                            parse_markovjunior_calc(inputs, loc, size_expr, GRAMMAR_CALC_SIZE_ENV, "Size"),
                            false)
    end
end

"
Parses a `repeat()` child like `Room(3)` or `Room(each(3:5))`,
  returning whether its size should be evaluated once per copy
"
function parse_markovjunior_grammar_repeat_child(inputs::MacroParserInputs, loc, expr,
                                                 symbol_lookup::Dict{Symbol, Int}
                                                )::Tuple{GrammarChild, Bool}
    per_copy = Base.isexpr(expr, :call) && (length(expr.args) == 2) &&
               Base.isexpr(expr.args[2], :call) && (expr.args[2].args[1] == :each)
    if per_copy
        if length(expr.args[2].args) != 2
            raise_parse_error(loc, inputs, "Expected `each(size)`; got `", expr.args[2], "`")
        end
        expr = Expr(:call, expr.args[1], expr.args[2].args[2])
    end
    child = parse_markovjunior_grammar_child(inputs, loc, expr, symbol_lookup; allow_weight=false)
    return (child, per_copy)
end

"
Separates a production's arguments into positional and keyword ones,
  checking that the keywords are allowed and not repeated
"
function parse_markovjunior_grammar_call_args(inputs::MacroParserInputs, loc, expr::Expr,
                                              allowed_keywords::Tuple{Vararg{Symbol}}
                                             )::Tuple{Vector{Any}, Dict{Symbol, Any}}
    positional = Any[ ]
    keywords = Dict{Symbol, Any}()
    add_keyword(kw) = begin
        if !Base.isexpr(kw, :kw) || !isa(kw.args[1], Symbol)
            raise_parse_error(loc, inputs, "Expected a keyword argument like `name=value`; got `", kw, "`")
        elseif !in(kw.args[1], allowed_keywords)
            raise_parse_error(loc, inputs, "`", expr.args[1], "()` has no `", kw.args[1], "` parameter",
                                isempty(allowed_keywords) ? "" : string(
                                    "; it takes ", iter_join(Iterators.map(k -> "`$k`", allowed_keywords), ", ")...
                                ))
        elseif haskey(keywords, kw.args[1])
            raise_parse_error(loc, inputs, "`", kw.args[1], "` is given more than once")
        end
        keywords[kw.args[1]] = kw.args[2]
    end
    for arg in expr.args[2:end]
        if Base.isexpr(arg, :parameters)
            foreach(add_keyword, arg.args)
        elseif Base.isexpr(arg, :kw)
            add_keyword(arg)
        else
            push!(positional, arg)
        end
    end
    return (positional, keywords)
end

"Parses face overrides like `-y => Floor`"
function parse_markovjunior_grammar_faces(inputs::MacroParserInputs, loc, exprs,
                                          symbol_lookup::Dict{Symbol, Int})::Vector{GrammarFace}
    parse_face(expr) = begin
        if !Base.isexpr(expr, :call) || (expr.args[1] != :(=>)) || (length(expr.args) != 3)
            raise_parse_error(loc, inputs, "Expected a face like `-y => Symbol`; got `", expr, "`")
        end
        GrammarFace(
            parse_markovjunior_grammar_axis(inputs, loc, expr.args[2]),
            parse_markovjunior_grammar_symbol(inputs, loc, expr.args[3], symbol_lookup)
        )
    end
    return GrammarFace[ parse_face(e) for e in exprs ]
end

"Whether the expression is a list inside a `@grammar`, written like `[a, b]` or `(a, b)`"
is_markovjunior_grammar_list(expr) = Base.isexpr(expr, :vect) || Base.isexpr(expr, :tuple)

"
Parses either one value for every axis, or a list with one value per axis,
  returning the values and whether it was a list
"
function parse_markovjunior_grammar_per_axis(parse_value, ::Type{T}, inputs::MacroParserInputs, loc, expr
                                            )::Tuple{Vector{T}, Bool} where {T}
    if is_markovjunior_grammar_list(expr)
        return (T[ parse_value(e) for e in expr.args ], true)
    else
        return (T[ parse_value(expr) ], false)
    end
end

"Parses an axis, which is a calculation; constant axes also put a floor on the grid's dimensionality"
function parse_markovjunior_grammar_axis(inputs::MacroParserInputs, loc, expr)::CalcProgram
    axis = parse_markovjunior_calc(inputs, loc, expr, GRAMMAR_CALC_ENV, "Axis")
    if calc_is_constant(axis)
        v = calc_constant_value(axis)
        if !isa(v, Int) || iszero(v)
            raise_parse_error(loc, inputs,
                              "Axis should be a non-zero Int, like `x`, `-y`, or `3`; got `", expr, "` = ", v)
        end
        inputs.min_dims = max(inputs.min_dims, abs(v))
    end
    return axis
end

function parse_markovjunior_grammar_production(inputs::MacroParserInputs, loc, expr,
                                               symbol_lookup::Dict{Symbol, Int}
                                              )::GrammarProduction
    if expr isa Char
        if !haskey(CELL_CODE_BY_CHAR, expr)
            raise_parse_error(loc, inputs, "Unsupported color value '", expr, "'")
        end
        return GrammarColor(CELL_CODE_BY_CHAR[expr])
    elseif expr == :_
        return GrammarKeep()
    elseif expr isa Symbol
        return GrammarRename(parse_markovjunior_grammar_symbol(inputs, loc, expr, symbol_lookup))
    elseif Base.isexpr(expr, :call) && (expr.args[1] == :split)
        if length(expr.args) < 3
            raise_parse_error(loc, inputs, "Expected `split(axis, children...)` with at least one child; got `",
                                expr, "`")
        end
        return GrammarSplit(
            parse_markovjunior_grammar_axis(inputs, loc, expr.args[2]),
            GrammarChild[
                parse_markovjunior_grammar_child(inputs, loc, c, symbol_lookup; allow_weight=true)
                  for c in expr.args[3:end]
            ]
        )
    elseif Base.isexpr(expr, :call) && (expr.args[1] == :repeat)
        (args, keywords) = parse_markovjunior_grammar_call_args(inputs, loc, expr, (:leftover, ))
        if !in(length(args), 2:3)
            raise_parse_error(loc, inputs, "Expected `repeat(axis, Child(size)[, Separator(size)][; leftover=...])`",
                                "; got `", expr, "`")
        end
        (child, child_per_copy) = parse_markovjunior_grammar_repeat_child(inputs, loc, args[2], symbol_lookup)
        (separator, separator_per_copy) = if length(args) >= 3
            parse_markovjunior_grammar_repeat_child(inputs, loc, args[3], symbol_lookup)
        else
            (nothing, false)
        end

        leftover_expr = get(keywords, :leftover, :stretch)
        leftover = if (leftover_expr isa Symbol) && (leftover_expr in GRAMMAR_LEFTOVER_MODES) &&
                      !in(leftover_expr, GRAMMAR_LEFTOVER_CHILD_MODES)
            GrammarLeftover(leftover_expr, 0)
        elseif Base.isexpr(leftover_expr, :call) && (length(leftover_expr.args) == 2) &&
               (leftover_expr.args[1] in GRAMMAR_LEFTOVER_CHILD_MODES)
            GrammarLeftover(leftover_expr.args[1],
                            parse_markovjunior_grammar_symbol(inputs, loc, leftover_expr.args[2], symbol_lookup))
        else
            raise_parse_error(loc, inputs, "Expected `leftover=` to be `stretch`, `spread`, `before(Symbol)`, ",
                                "`after(Symbol)`, or `around(Symbol)`; got `", leftover_expr, "`")
        end
        if (leftover.mode == :spread) && isnothing(separator)
            raise_parse_error(loc, inputs, "`leftover=spread` grows the separators, so it needs a separator")
        end

        return GrammarRepeat(
            parse_markovjunior_grammar_axis(inputs, loc, args[1]),
            child, child_per_copy,
            separator, separator_per_copy,
            leftover
        )
    elseif Base.isexpr(expr, :call) && (expr.args[1] == :inset)
        (args, _) = parse_markovjunior_grammar_call_args(inputs, loc, expr, ())
        if length(args) < 3
            raise_parse_error(loc, inputs, "Expected `inset(amount, Shell, Core, faces...)`; got `", expr, "`")
        end
        parse_amount(a) = parse_markovjunior_calc(inputs, loc, a, GRAMMAR_CALC_ENV, "Inset amount")
        (amounts, per_axis) = parse_markovjunior_grammar_per_axis(GrammarInsetAmount, inputs, loc, args[1]) do a
            if is_markovjunior_grammar_list(a) && (length(a.args) == 2)
                GrammarInsetAmount(parse_amount(a.args[1]), parse_amount(a.args[2]))
            elseif is_markovjunior_grammar_list(a)
                raise_parse_error(loc, inputs, "Expected an inset amount, or a `[low, high]` pair of them; got `",
                                    a, "`")
            else
                GrammarInsetAmount(parse_amount(a), nothing)
            end
        end
        if isempty(amounts)
            raise_parse_error(loc, inputs, "`inset()` needs an amount")
        end
        return GrammarInset(
            amounts, per_axis,
            parse_markovjunior_grammar_symbol(inputs, loc, args[2], symbol_lookup),
            parse_markovjunior_grammar_symbol(inputs, loc, args[3], symbol_lookup),
            parse_markovjunior_grammar_faces(inputs, loc, args[4:end], symbol_lookup)
        )
    elseif Base.isexpr(expr, :call) && (expr.args[1] == :place)
        (args, _) = parse_markovjunior_grammar_call_args(inputs, loc, expr, ())
        if length(args) < 3
            raise_parse_error(loc, inputs, "Expected `place(box, Inner, Outer, faces...)`; got `", expr, "`")
        end
        return GrammarPlace(
            parse_markovjunior_calc_box(inputs, loc, args[1], GRAMMAR_CALC_BOX_ENV),
            parse_markovjunior_grammar_symbol(inputs, loc, args[2], symbol_lookup),
            parse_markovjunior_grammar_symbol(inputs, loc, args[3], symbol_lookup),
            parse_markovjunior_grammar_faces(inputs, loc, args[4:end], symbol_lookup)
        )
    elseif Base.isexpr(expr, :call) && (expr.args[1] == :mirror)
        (args, _) = parse_markovjunior_grammar_call_args(inputs, loc, expr, ())
        if !in(length(args), 2:3) || ((length(args) == 3) && (args[3] != :layout))
            raise_parse_error(loc, inputs, "Expected `mirror(axes, Source[, layout])`; got `", expr, "`")
        end
        (axes, axes_is_list) = parse_markovjunior_grammar_per_axis(CalcProgram, inputs, loc, args[1]) do a
            parse_markovjunior_grammar_axis(inputs, loc, a)
        end
        if isempty(axes)
            raise_parse_error(loc, inputs, "`mirror()` needs at least one axis")
        end
        return GrammarMirror(
            axes, axes_is_list,
            parse_markovjunior_grammar_symbol(inputs, loc, args[2], symbol_lookup),
            length(args) == 3
        )
    elseif any(h -> Base.isexpr(expr, h), (:hcat, :vcat, :ncat, :vect))
        block = parse_markovjunior_rewrite_block(inputs, loc, expr,
                                                 allow_wildcard=true, allow_set=true,
                                                 allow_list=false, allow_lookup=false)
        n_block_dims = scale_ops_block_dims!(inputs, (block, ))
        return GrammarBlock(scale_ops_normalize_block(UpscaleCellDest, block, n_block_dims))
    elseif Base.isexpr(expr, :macrocall) || Base.isexpr(expr, :block)
        ops = parse_markovjunior_sequence_main_block(expr, loc, inputs)
        check_markovjunior_scoped_ops(inputs, loc, "@grammar", ops)
        return GrammarOps(ops)
    else
        raise_parse_error(loc, inputs,
                          "Expected `split()`, `repeat()`, `inset()`, `place()`, `mirror()`, another symbol, a 'C' color, `_`, ",
                            "a `[ ... ]` block, or Ops; got `", expr, "`")
    end
end

function parse_markovjunior_op(::Val{Symbol("@grammar")}, inputs::MacroParserInputs,
                               loc, expr_args, expr_line)
    args = collect(Any, expr_args)
    inherit_biases = !isempty(args) && (args[1] == :inherit_biases)
    inherit_biases && popfirst!(args)
    if (length(args) != 1) || !Base.isexpr(args[1], :block)
        raise_parse_error(loc, inputs, "Expected `@grammar [inherit_biases] begin ... end`")
    end

    # Find every rule's symbol first, so rules can refer to symbols that are defined later.
    # The root symbol always comes first, wherever its rules are.
    symbols = Symbol[ GRAMMAR_ROOT_SYMBOL ]
    symbol_lookup = Dict{Symbol, Int}(GRAMMAR_ROOT_SYMBOL => 1)
    rule_sources = Tuple{Any, Int, Any, Any}[ ] # Location, symbol, guard, and right-hand side
    parse_markovjunior_block(args[1].args) do rule_loc, line
        rule_loc = something(rule_loc, loc)
        try_parse_markovjunior_calc!(inputs, rule_loc, line) && return nothing
        with_parser_stacktrace(inputs, "Rule $(length(rule_sources) + 1)") do
            if !@capture(line, lhs_ => rhs_)
                raise_parse_error(rule_loc, inputs, "Expected a rule like `Name => production`; got `", line, "`")
            end
            (name, guard) = if lhs isa Symbol
                (lhs, nothing)
            elseif Base.isexpr(lhs, :call) && (length(lhs.args) == 2) && (lhs.args[1] isa Symbol)
                (lhs.args[1], lhs.args[2])
            else
                raise_parse_error(rule_loc, inputs,
                                  "Expected a symbol, optionally with a guard like `Room(size_min() > 3)`; got `",
                                    lhs, "`")
            end
            if name == :_
                raise_parse_error(rule_loc, inputs, "`_` means 'leave untouched', and can't have rules")
            end

            symbol_i = get!(symbol_lookup, name) do
                push!(symbols, name)
                length(symbols)
            end
            push!(rule_sources, (rule_loc, symbol_i, guard, rhs))
        end
    end
    if !any(((_, symbol_i, _, _), ) -> (symbol_i == 1), rule_sources)
        raise_parse_error(loc, inputs, "`@grammar` needs a `", GRAMMAR_ROOT_SYMBOL, "` rule, for the whole grid")
    end

    # Parse each rule.
    # Without inherited biases, Op leaves shouldn't be validated against the outer biases either.
    rules = [ GrammarRule[ ] for _ in symbols ]
    outer_bias_stack = inputs.bias_stack
    inherit_biases || (inputs.bias_stack = Stack{Vector{AbstractMarkovBias}}(16))
    try
        for (rule_i, (rule_loc, symbol_i, guard_expr, rhs)) in enumerate(rule_sources)
            with_parser_stacktrace(inputs, "Rule $rule_i `$(symbols[symbol_i])`") do
                guard = if isnothing(guard_expr)
                    nothing
                else
                    parse_markovjunior_calc(inputs, rule_loc, guard_expr, GRAMMAR_CALC_ENV, "Guard")
                end
                (production_expr, weight_expr) = peel_grammar_weight(rhs)
                weight = parse_markovjunior_calc(inputs, rule_loc, weight_expr, GRAMMAR_CALC_ENV, "Weight")
                production = parse_markovjunior_grammar_production(inputs, rule_loc, production_expr,
                                                                   symbol_lookup)
                push!(rules[symbol_i], GrammarRule(guard, weight, production))
            end
        end
    finally
        inputs.bias_stack = outer_bias_stack
    end

    return MarkovOpGrammar(inherit_biases, symbols, rules)
end
