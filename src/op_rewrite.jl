####################################
#  Rewrite rule: individual Cells

const RewriteRuleCell_Set = CellTypeSet
const RewriteRuleCell_List = UpTo{N_CELL_TYPES, UInt8}
struct RewriteRuleCell_Wildcard end
struct RewriteRuleLookup
    source_idx::Int
end

"The source entries for a rewrite rule can be a single value, ordered set, or a wildcard"
const RewriteRuleCellSource = Union{UInt8, RewriteRuleCell_Set, RewriteRuleCell_Wildcard}
"The destination entries for a rewrite rule can be a single value, unordered set, list (matching source set size), or a wildcard"
const RewriteRuleCellDest = Union{UInt8, RewriteRuleCell_Set,
                                  RewriteRuleCell_List, RewriteRuleCell_Wildcard,
                                  RewriteRuleLookup}

const RewriteCell = Tuple{RewriteRuleCellSource, RewriteRuleCellDest}

match_rewrite_source(rule::UInt8, value::UInt8) = (rule == value)
match_rewrite_source(rule::RewriteRuleCell_Set, value::UInt8) = (value in rule)
match_rewrite_source(::RewriteRuleCell_Wildcard, ::UInt8) = true

pick_rewrite_value(dest::UInt8,                    src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = dest
pick_rewrite_value(dest::RewriteRuleCell_Set,      src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = rand(rng, dest)
pick_rewrite_value(dest::RewriteRuleCell_List,     src::RewriteRuleCell_Set,   src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = dest[cell_set_index_of(src, src_values[self_idx])]
pick_rewrite_value(dest::RewriteRuleLookup,        src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = src_values[dest.source_idx]
pick_rewrite_value(dest::RewriteRuleCell_Wildcard, src::RewriteRuleCellSource, src_values::Tuple{Vararg{UInt8}}, self_idx::Int, rng::PRNG) = src_values[self_idx]


###########################
#  Rewrite rule: 1D strip

struct RewriteRule_Strip{NCells, TCells<:NTuple{NCells, RewriteCell}}
    cells::TCells
    mask::Union{Nothing, Float32, NTuple{2, Float32}} # Disabled, constant, or random range
    weight::Float32
    explicit_symmetries::Vector{GridDir}
    unlimited_symmetries_after_axis::Optional{Int} # If <1, then all symmetries are allowed and 'explicit_symmetries' is empty
end

const MaskGrid{N} = Array{Float32, N}

pick_mask(r::RewriteRule_Strip, rng::PRNG)::Float32 = pick_mask_impl(r.mask, rng)
pick_mask_impl(::Nothing, ::PRNG) = 1.0f0
pick_mask_impl(f::Float32, ::PRNG) = f
pick_mask_impl((a, b)::NTuple{2, Float32}, rng::PRNG) = lerp(a, b, rand(rng, Float32))

"
Checks for a rule matching against the given grid when applied to the given cell.
Assumes for performance that the rule does fit into the grid.
"
function rule_matches(r::RewriteRule_Strip{NCells, TCells}, grid::CellGrid{NDims},
                      c_start::CellIdx{NDims}, dir::GridDir
                     )::Bool where {NCells, TCells, NDims}
    @markovjunior_assert begin
        c_end = grid_dir_pos_along(dir, c_start, NCells - 1)
        grid_range = one(VecI{NDims}):vsize(grid)
        (c_start in grid_range) && (c_end in grid_range)
    end "Rule is outside grid bounds! $c_start along $dir by $(NCells-1) vs $(vsize(grid))"

    return all(ntuple(Val(NCells)) do i
        c_next = grid_dir_pos_along(dir, c_start, i-1)
        return match_rewrite_source(r.cells[i][1], grid[c_next])
    end)
end

"
Finds every match for a given rule that touches a given subset of the grid,
  excluding rule placements that start at cells which fail the given mask.

For each potential match your lambda receives `(start_cell, grid_dir, match_succeeded)`.
Note this does not include cells which are outside the mask.

If your lambda returns anything other than `nothing`,
  the function immediately stops checking for matches and returns that value;
  otherwise we return `nothing` at the end of the search.
"
function visit_rule_match_data(process_candidate::TLambda,
                               r::RewriteRule_Strip{NCells, TCells}, grid::CellGrid{NDims},
                               grid_mask::MaskGrid{NDims}, rule_chosen_mask::Float32,
                               inner_subset::Bplus.Math.BoxI{NDims}
                              )::Optional where {NCells, TCells, NDims, TLambda}
    function process_dir(dir::GridDir)
        clamp_range = (
            if dir.sign > 0
                1
            else
                NCells
            end,
            size(grid, dir.axis) + if dir.sign > 0
                -(NCells - 1)
            else
                0
            end
        )

        first_pos = min_inclusive(inner_subset)
        if dir.sign > 0
            first_pos = grid_dir_pos_along(dir, first_pos, -(NCells - 1))
        end
        @set! first_pos[dir.axis] = clamp(first_pos[dir.axis], clamp_range...)

        last_pos = max_inclusive(inner_subset)
        if dir.sign < 0
            last_pos = grid_dir_pos_along(dir, last_pos, -(NCells - 1))
        end
        @set! last_pos[dir.axis] = clamp(last_pos[dir.axis], clamp_range...)

        @markovjunior_assert(
            (first_pos, last_pos) == minmax(first_pos, last_pos),
            "Mismatch! $((first_pos, last_pos))"
        )

        for rule_start_pos::CellIdx{NDims} in first_pos:last_pos
            if grid_mask[rule_start_pos] <= rule_chosen_mask
                user_out = process_candidate(rule_start_pos, dir,
                                             rule_matches(r, grid, rule_start_pos, dir))
                exists(user_out) && return user_out
            end
        end
    end
    for dir in r.explicit_symmetries
        process_dir(dir)
    end
    if exists(r.unlimited_symmetries_after_axis)
        for axis in (r.unlimited_symmetries_after_axis+1):NDims
            for sign in (-1, 1)
                process_dir(GridDir(axis, sign))
            end
        end
    end
    return nothing
end


##################################
#  Rewrite rule: cache of matches

"Efficiently tracks all possible applications of rewrite rules to a grid"
struct RewriteCache{NDims, NRules, TGrid<:CellGrid{NDims}, TRules<:NTuple{NRules, RewriteRule_Strip}}
    grid::TGrid
    rules::TRules

    mask_grid::MaskGrid{NDims}
    rule_masks::NTuple{NRules, Float32}

    applications::Vector{OrderedSet{Tuple{CellIdx{NDims}, GridDir}}}
end

function RewriteCache(grid::CellGrid{NDims}, mask_grid::MaskGrid{NDims},
                      rules::NTuple{NRules, RewriteRule_Strip},
                      mask_for_each_rule::NTuple{NRules, Float32},
                      context::MarkovOpContext
                     ) where {NDims, NRules}
    whole_grid_range = Box(min=one(CellIdx{NDims}), max=vsize(grid))

    logic_logln("Caching valid rule applications...")
    logic_tab_in()
    applications_tuple = map(rules, mask_for_each_rule) do rule, mask
        set = markov_allocator_acquire_ordered_set(context.allocator, Tuple{CellIdx{NDims}, GridDir})
        empty!(set)

        logic_logln("Masked ", mask, ", rule ", rules)
        logic_tab_in()
        visit_rule_match_data(rule, grid, mask_grid, mask, whole_grid_range) do cell, dir, is_matching
            if is_matching
                logic_log("   ", cell, " along ", dir.axis, "|", (dir.sign > 0 ? "+" : ""), dir.sign)
                push!(set, (cell, dir))
            end
            return nothing
        end
        logic_logln()
        logic_tab_out()
        return set
    end
    logic_logln("Total candidates: ", sum(map(length, applications_tuple), init=0))
    logic_tab_out()

    applications = markov_allocator_acquire_array(
        context.allocator,
        tuple(length(applications_tuple)),
        OrderedSet{Tuple{CellIdx{NDims}, GridDir}}
    )
    empty!(applications)
    append!(applications, applications_tuple)

    return RewriteCache{NDims, NRules, typeof(grid), typeof(rules)}(
        grid, rules,
        mask_grid, mask_for_each_rule,
        applications
    )
end

"Updates a rewrite cache, given an area of the grid that potentially changed"
function update_rewrite_cache!(cache::RewriteCache{NDims, NRules, TGrid, TRules},
                               range_to_invalidate::BoxI{NDims}
                              )::Nothing where {NDims, NRules, TGrid, TRules}
    foreach(cache.rules, cache.rule_masks, cache.applications) do rule, mask, set
        visit_rule_match_data(rule, cache.grid, cache.mask_grid, mask, range_to_invalidate
                             ) do cell, dir, is_matching
            used_to_match = ((cell, dir) in set)
            if is_matching && !used_to_match
                push!(set, (cell, dir))
            elseif !is_matching && used_to_match
                delete!(set, (cell, dir))
            end
            return nothing
        end
    end
end

function close_rewrite_cache(cache::RewriteCache, context::MarkovOpContext)
    foreach(app -> markov_allocator_release_ordered_set(context.allocator, app), cache.applications)
    markov_allocator_release_array(context.allocator, cache.applications)
    return nothing
end


#########################
#  Rewrite Op: 1D strip

"A simple MarkovJunior rewrite op, affecting a 1D strip of pixels"
struct MarkovOpRewrite1D{TRules <: Tuple{Vararg{RewriteRule_Strip}},
                         TBias <: Tuple{Vararg{AbstractMarkovBias}}
                        } <: AbstractMarkovOp
    rules::TRules
    threshold::Optional{Threshold}
    biases::TBias
end

mutable struct MarkovOpRewrite1D_State{NDims, NRules, TGrid, TRules<:NTuple{NRules, RewriteRule_Strip},
                                       NBiases, TFullBias<:NTuple{NBiases, AbstractMarkovBias},
                                                TBiasStates<:NTuple{NBiases, Any}}
    applications_left::Optional{Int} # Infinite if not given
    rewrite_cache::RewriteCache{NDims, NRules, TGrid, TRules}
    biases::TFullBias
    bias_states::TBiasStates
    weighted_options_buffer::Vector{Tuple{Int, Float32, CellIdx{NDims}, GridDir}}
    weight_data_buffer::Base.RefValue{NTuple{3, Float32}} # Min, Max, Sum
end

function markov_op_initialize(r::MarkovOpRewrite1D{<:NTuple{NRules, RewriteRule_Strip}, TBias},
                              grid::CellGrid{NDims}, rng::PRNG,
                              context::MarkovOpContext
                             ) where {NDims, NRules, TBias}
    mask_grid = markov_allocator_acquire_array(context.allocator, size(grid), Float32)
    rand!(rng, mask_grid)

    cache = RewriteCache(grid, mask_grid, r.rules,
                         map(ru -> pick_mask(ru, rng), r.rules),
                         context)
    biases = (
        context.inherited_biases_type_stable...,
        r.biases...
    )
    TBiasStates = Tuple{map(markov_bias_state_type, typeof.(biases))...}
    bias_states::TBiasStates = map(b -> markov_bias_initialize(b, grid, rng, context.bias_context), biases)

    threshold = isnothing(r.threshold) ? nothing : get_threshold(r.threshold, grid, rng)

    out_state = MarkovOpRewrite1D_State{NDims, NRules, typeof(grid), typeof(r.rules),
                                        length(biases), typeof(biases), TBiasStates}(
        threshold, cache,
        biases, bias_states,
        markov_allocator_acquire_array(context.allocator, tuple(128),
                                       Tuple{Int, Float32, CellIdx{NDims}, GridDir}),
        Ref{NTuple{3, Float32}}()
    )
    if all(isempty, cache.applications)
        logic_logln("MarkovOpRewrite1D has no options at the start; canceling...")
        markov_op_cancel(r, out_state, context)
        return nothing
    else
        logic_logln("MarkovOpRewrite1D has been initialized; there are ",
                      sum(map(length, cache.applications), init=0), " initial options")
        logic_logln("The actual threshold is ",
                      typeof(out_state.applications_left), "(", out_state.applications_left, ")")
        return out_state
    end
end
function markov_op_iterate(r::MarkovOpRewrite1D{TRules, TSelfBiases},
                           state::MarkovOpRewrite1D_State{NDims, NRules, TGrid, TRules, NBiases, TFullBias, TBiasStates},
                           grid::CellGrid{NDims}, rng::PRNG,
                           context::MarkovOpContext,
                           ticks_left::Ref{Optional{Int}}
                          ) where {NDims, NRules, TGrid, TRules,
                                   NBiases, TFullBias, TBiasStates, TSelfBiases}
    @markovjunior_assert(get_something(state.applications_left, 1) > 0)
    logic_logln("Remaining threshold for this rewrite op: ",
                get_something(state.applications_left, "infinity"))

    while isnothing(ticks_left[]) || (ticks_left[] > 0)
        # Get all the options and their weights.
        empty!(state.weighted_options_buffer)
        state.weight_data_buffer[] = (Inf32, -Inf32, 0.0f0)
        foreach(r.rules, 1:NRules) do rule, rule_i
            for (start_cell, dir) in state.rewrite_cache.applications[rule_i]
                cell_line = CellLine(start_cell, dir, convert(Int32, length(rule.cells)))
                biases::NTuple{NBiases, Optional{Float32}} = markov_bias_calculate.(
                    state.biases, state.bias_states,
                    Ref(grid), Ref(cell_line), Ref(rng)
                )
                if all(exists, biases)
                    bias = sum(biases, init=zero(Float32))
                    push!(state.weighted_options_buffer, (
                        rule_i, bias,
                        start_cell, dir
                    ))

                    state.weight_data_buffer[] = let d = state.weight_data_buffer[]
                        (
                            min(d[1], bias),
                            max(d[2], bias),
                            d[3] + bias
                        )
                    end
                end
            end
        end

        logic_logln("There are ", length(state.weighted_options_buffer),
                      " options with biases ranging from ",
                      state.weight_data_buffer[][1], " to ", state.weight_data_buffer[][2])

        # If no options are left, we're done.
        if isempty(state.weighted_options_buffer)
            markov_op_cancel(r, state, context)
            return nothing
        end

        # Pick an option -- if all weights are equal, then use trivial uniform-random selection.
        (pick_rule_i, pick_start_cell, pick_dir) =
          if state.weight_data_buffer[][1] == state.weight_data_buffer[][2]
            logic_logln("Luckily we can use uniform-random selection, which is much faster")
            pick_option_i = rand(rng, 1:length(state.weighted_options_buffer))
            choice = state.weighted_options_buffer[pick_option_i]
            (choice[1], choice[3], choice[4])
        else
            # Normalize the bias weights by subtracting the min then dividing by their sum.
            # Adding the min to N weights means increasing the sum by min*N.
            weight_norm_offset = -state.weight_data_buffer[][1]
            weight_norm_scale = 1.0f0 / (state.weight_data_buffer[][3] *
                                        length(state.weighted_options_buffer))
            logic_logln("Bias weights transform: (+", weight_norm_offset, ")*", weight_norm_scale)

            # Pick an option using weighted randomness.
            pick_t = rand(rng, Float32)
            (() -> begin
                picked_f = 0.0f0
                for (rule_i, raw_bias, start_cell, dir) in state.weighted_options_buffer
                    bias = (raw_bias + weight_norm_offset) * weight_norm_scale
                    if picked_f <= pick_t + bias
                        return (rule_i, start_cell, dir)
                    else
                        picked_f += bias
                    end
                end
                # If we get here, then probably we picked a value near 1.0
                #    and floating-point error kept us from choosing the last element.
                return let answer = state.weighted_options_buffer[end]
                    (answer[1], answer[3], answer[4])
                end
            end)()
        end
        logic_logln("Picked rule ", pick_rule_i, " at ", pick_start_cell, " along ", pick_dir)

        # Apply the rule.
        # Because each rule is a different type but known at compile-time,
        #    we should add a layer of dispatch when executing it.
        #TODO: Include grid dir in the dispatch data; profile
        rule_len::Int32 = ((rule::RewriteRule_Strip) -> begin
            source_values = Tuple(
                grid[grid_dir_pos_along(pick_dir, pick_start_cell, i-1)]
                  for i in 1:length(rule.cells)
            )
            # Each rule's rewrite cell is also a different type known at compile-time.
            foreach(rule.cells, 1:length(rule.cells)) do (rewrite_source, rewrite_dest), cell_i
                cell_pos = grid_dir_pos_along(pick_dir, pick_start_cell, cell_i-1)
                grid[cell_pos] = pick_rewrite_value(rewrite_dest, rewrite_source,
                                                    source_values, cell_i,
                                                    rng)
            end

            return convert(Int32, length(rule.cells))
        end)(r.rules[pick_rule_i])

        # Update counters.
        pick_end_cell = grid_dir_pos_along(pick_dir, pick_start_cell, rule_len-1)
        (pick_start_cell, pick_end_cell) = minmax(pick_start_cell, pick_end_cell)
        update_rewrite_cache!(state.rewrite_cache, Box(pick_start_cell:pick_end_cell))
        if exists(ticks_left[])
            ticks_left[] -= 1
        end
        if exists(state.applications_left)
            state.applications_left -= 1
            if state.applications_left < 1
                markov_op_cancel(r, state, context)
                return nothing
            end
        end
    end

    return state
end

function markov_op_cancel(op::MarkovOpRewrite1D, s::MarkovOpRewrite1D_State,
                          context::MarkovOpContext)
    markov_allocator_release_array(context.allocator, s.rewrite_cache.mask_grid)
    markov_allocator_release_array(context.allocator, s.weighted_options_buffer)
    foreach((b, s) -> markov_bias_cleanup(b, s), s.biases, s.bias_states)
    close_rewrite_cache(s.rewrite_cache, context)
end


#####################
#  DSL integration

dsl_string_rewrite_source(rule::UInt8) = dsl_string(rule)
dsl_string_rewrite_source(rule::RewriteRuleCell_Set) = "[$(dsl_string(rule))]"
dsl_string_rewrite_source(rule::RewriteRuleCell_Wildcard) = "_"

dsl_string_rewrite_dest(rule::UInt8) = dsl_string(rule)
dsl_string_rewrite_dest(rule::RewriteRuleCell_Set) = "{$(dsl_string(rule))}"
dsl_string_rewrite_dest(rule::RewriteRuleCell_List) = "[$(string(dsl_string.(rule)...))]"
dsl_string_rewrite_dest(rule::RewriteRuleCell_Wildcard) = "_"
dsl_string_rewrite_dest(rule::RewriteRuleLookup) = "[$(rule.source_idx)]"

dsl_string(strip::RewriteRule_Strip) = string(
    dsl_string_rewrite_source.(t[1] for t in strip.cells)...,
    " => ",
    dsl_string_rewrite_dest.(t[2] for t in strip.cells)...,
    if isnothing(strip.mask)
        ()
    elseif strip.mask isa Float32
        (" %", strip.mask)
    elseif strip.mask isa NTuple{2, Float32}
        (" %(", strip.mask[1], ":", strip.mask[2], ")")
    else
        error("Unhandled: ", typeof(strip.mask))
    end...,
    (isone(strip.weight) ? () : (" *", strip.weight))...,
    #TODO: write symmetries, update DSL for '...' at end of symmetry list (multidimensional?)
)

dsl_string(r::MarkovOpRewrite1D)