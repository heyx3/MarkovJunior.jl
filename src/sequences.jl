"
Executes some rules/inner sequences on a grid, according to some logic.

Works much like an iterator:
 * `start_sequence()` returns some state.
 * `execute_sequence()` takes that state and returns the next one.
"
abstract type AbstractSequence end
"Returns the state used for the first iteration of the sequence."
start_sequence(  seq::AbstractSequence, grid::CellGrid, outer_inference::AllInference, prng       )::Any           = error("Unimplemented: ", typeof(seq))
"Executes the given iteration of the sequence and returns the next iteration (or `nothing` if finished)"
execute_sequence(seq::AbstractSequence, grid::CellGrid, prng, state)::Optional{Any} = error("Unimplemented: ", typeof(seq))
"
Tries to broadcast a sequence to an n-dimensional version of itself.
Mainly used for dimension-specific sequences to allow 1D versions to apply to any number of dimensions.

Default behavior: returns itself.
"
broadcast_sequence(seq::AbstractSequence, new_dims::Int) = seq
"Asks an executing sequence for the current state of Inference logic"
inference_of_sequence(seq::AbstractSequence, current_state)::Optional{AllInference_State} = error("Unimplemented: ", typeof.((seq, current_state)))


"Executes a set of rules exactly N times"
struct Sequence_DoN <: AbstractSequence
    rules::Vector{CellRule}
    count::Int # Secret internal behavior: if count is -666, then it is infinite.
    sequential::Bool # If true, then earlier rules always have priority
    inference::AllInference

    Sequence_DoN(rules, count, sequential::Bool, inference=AllInference()) = new(
        if rules isa Vector{CellRule}
            rules
        else
            collect(CellRule, rules)
        end,
        if (count < 0) && (count != -666)
            error("Invalid count: ", count)
        else
            convert(Int, count)
        end,
        sequential,
        inference
    )
end
start_sequence(d::Sequence_DoN, grid::CellGrid{N}, outer_inference::AllInference, rng::PRNG) where {N} = (
    1,
    RuleCache(grid, d.rules),
    AllInference_State(AllInference(outer_inference, d.inference), grid),
    Vector{Tuple{RuleApplication{N}, Float32}}()
)
function execute_sequence(d::Sequence_DoN, grid::CellGrid{N}, rng::PRNG,
                          (next_i, cache, inference, applications_buffer
                          )::Tuple{Int, RuleCache{N},
                                   AllInference_State{N},
                                   Vector{Tuple{RuleApplication{N}, Float32}}}
                         ) where {N}
    use_inference::Bool = inference_exists(inference.source)
    if SKIP_CACHE
        if (d.count >= 0) && (next_i > d.count)
            return nothing
        end

        # Get all possible applications.
        empty!(applications_buffer)
        min_inference_weight = Ref(Inf32)
        first_matching_rule_idx = Ref{Int32}(0)
        find_rule_matches(grid, d.rules) do i, c
            inference_weight = use_inference ? infer_weight(inference, d.rules[i], c, rng) : 0.0f0

            if first_matching_rule_idx[] <= 0
                first_matching_rule_idx[] = i
                min_inference_weight[] = Inf32
            end

            if exists(inference_weight)
                min_inference_weight[] = min(min_inference_weight[], inference_weight)
                #TODO: Exit early (requires rewriting the helper functions) once we hit the second rule idx

                push!(applications_buffer, (RuleApplication{N}(i, c), inference_weight))
            end
        end
        if isempty(applications_buffer)
            return nothing
        end

        # Pick an applicable option.
        filter!(t -> (!d.sequential || t[1].rule_idx == first_matching_rule_idx[]) &&
                     (t[2] == min_inference_weight[]),
                applications_buffer)
        option = rand(rng, applications_buffer)
        chosen_rule_idx = option[1].rule_idx
        chosen_cell_line = option[1].at
        chosen_inference_weight = option[2]

        rule_execute!(grid, d.rules[chosen_rule_idx], chosen_cell_line)
        recalculate!(inference)
        return (next_i + 1, cache, inference, applications_buffer)
    end

    # Are we out of options/time? Then give up.
    n_options = n_cached_rule_applications(cache)
    if (n_options < 1) || ((d.count >= 0) && (next_i > d.count))
        return nothing
    end

    # Pick the next rule application to execute.
    application::RuleApplication{N} = if d.sequential
        # If running sequentially, we need to work with only the first rule that has applications.
        try_rule_set = get_first_nonempty_cached_rule_application_set(cache)
        @markovjunior_assert exists(try_rule_set)
        (rule_idx, applications_set) = try_rule_set
        rule = d.rules[rule_idx]

        empty!(applications_buffer)
        min_inference_weight = Ref(Inf32)
        for (at, dir) in applications_set
            line = CellLine(at, dir, rule.length)
            app = RuleApplication(rule_idx, line)

            inference_weight = use_inference ? infer_weight(inference, rule, line, rng) : 0.0f0
            if exists(inference_weight)
                if inference_weight < min_inference_weight[]
                    empty!(applications_buffer)
                    min_inference_weight[] = inference_weight
                end
                if inference_weight == min_inference_weight[]
                    push!(applications_buffer, (app, inference_weight))
                end
            end
        end

        @markovjunior_assert !isempty(applications_buffer)
        rand(rng, applications_buffer)[1]
    else
        # Otherwise we must check *all* rules.
        if use_inference
            # Randomly select from the options with the lowest inference weight.
            min_inference_weight = Ref(Inf32)
            empty!(applications_buffer)
            for_each_cached_rule_application(cache) do app::RuleApplication{N}
                next_inference_weight = infer_weight(inference, d.rules[app.rule_idx], app.at, rng)
                if exists(next_inference_weight)
                    # We only care about applications with the largest inference weighting.
                    if next_inference_weight < min_inference_weight[]
                        min_inference_weight[] = next_inference_weight
                        empty!(applications_buffer)
                    end
                    if next_inference_weight == min_inference_weight[]
                        push!(applications_buffer, (app, next_inference_weight))
                    end
                end
            end

            @markovjunior_assert !isempty(applications_buffer)
            rand(rng, applications_buffer)[1]
        else
            # Choose an option with uniform-random odds.
            get_cached_rule_application(
                cache,
                rand(rng, 1:n_options)
            )
        end
    end

    # Apply the rule and update secondary data accordingly.
    rule = d.rules[application.rule_idx]
    @markovjunior_assert(
        rule_applies(grid, rule, application.at),
        "Rule $(application.rule_idx) ($rule) at $(application.at), starting with $(grid[application.at.start_cell])"
    )
    rule_execute!(grid, rule, application.at)
    update_cache!(cache, grid, application.rule_idx, application.at)
    recalculate!(inference)

    return (next_i + 1, cache, inference, applications_buffer)
end
dsl_string(dn::Sequence_DoN) = string("@do_n ", dn.count, " begin
    ",
    dn.sequential                  ?                 "@sequential\n    " : "",
    inference_exists(dn.inference) ? "$(dsl_string(dn.inference))\n    " : "",
    (iter_join(dsl_string.(dn.rules), "\n    "))..., "
end")
inference_of_sequence(::Sequence_DoN, (next_i, cache, inference, applications_buffer
                                      )::Tuple{Int, RuleCache{N},
                                               AllInference_State{N},
                                               Vector{Tuple{RuleApplication{N}, Float32}}}
                     ) where {N} = inference


"Executes a set of rules `round(N*C)` times, where C is the number of cells in the grid"
struct Sequence_DoNRelative <: AbstractSequence
    rules::Vector{CellRule}
    count_per_cell::Float64
    count_min::Int
    sequential::Bool # If true, then earlier rules always have priority
    inference::AllInference

    Sequence_DoNRelative(rules, count_per_cell, count_min,
                         sequential::Bool,
                         inference=AllInference()) = new(
        if rules isa Vector{CellRule}
            rules
        else
            collect(CellRule, rules)
        end,
        if count_per_cell < 0
            error("Invalid count_per_cell: ", count_per_cell)
        else
            convert(Int, count_per_cell)
        end,
        if count_min < 0
            error("Invalid count_min: ", count_min)
        else
            convert(Int, count_min)
        end,
        sequential,
        inference
    )
end
function start_sequence(d::Sequence_DoNRelative, grid::CellGrid{N}, outer_inference::AllInference, rng::PRNG) where {N}
    n_cells = prod(size(state.grid))
    count = max(d.count_min, Int(round(n_cells * state.count_per_cell)))

    inner_sequence = Sequence_DoN(d.rules, count, d.sequential, d.inference)
    inner_state = start_sequence(inner_sequence, grid, outer_inference, rng)

    return (inner_sequence, inner_state)
end
function execute_sequence(d::Sequence_DoNRelative, grid::CellGrid{N}, rng::PRNG,
                          (inner_sequence, inner_state)::Tuple{Sequence_DoN, Any}
                         ) where {N}
    new_inner_state = execute_sequence(inner_sequence, grid, rng, inner_state)
    if isnothing(new_inner_state)
        return nothing
    else
        return (inner_sequence, new_inner_state)
    end
end
inference_of_sequence(::Sequence_DoNRelative, (inner_sequence, inner_state)::Tuple{Sequence_DoN, Any}) =
    inference_of_sequence(inner_sequence, inner_state)

"Executes a set of rules until there are no more matches"
struct Sequence_DoAll <: AbstractSequence
    rules::Vector{CellRule}
    sequential::Bool # If true, then earlier rules always have priority
    inference::AllInference

    Sequence_DoAll(rules, sequential::Bool, inference=AllInference()) = new(
        if rules isa Vector{CellRule}
            rules
        else
            collect(CellRule, rules)
        end,
        sequential,
        inference
    )
end
function start_sequence(d::Sequence_DoAll, grid::CellGrid{N}, outer_inference::AllInference, rng::PRNG) where {N}
    inner_sequence = Sequence_DoN(d.rules, -666, d.sequential, d.inference)
    inner_state = start_sequence(inner_sequence, grid, outer_inference, rng)
    return (inner_sequence, inner_state)
end
function execute_sequence(d::Sequence_DoAll, grid::CellGrid{N}, rng::PRNG,
                          (inner_sequence, inner_state)::Tuple{Sequence_DoN, Any}) where {N}
    new_inner_state = execute_sequence(inner_sequence, grid, rng, inner_state)
    if isnothing(new_inner_state)
        return nothing
    else
        return (inner_sequence, new_inner_state)
    end
end
dsl_string(da::Sequence_DoAll) = string("@do_all begin
    ",
    da.sequential                  ?                 "@sequential\n    " : "",
    inference_exists(da.inference) ? "$(dsl_string(da.inference))\n    " : "",
    (iter_join(dsl_string.(da.rules), "\n    "))..., "
end")
inference_of_sequence(::Sequence_DoAll, (inner_sequence, inner_state)::Tuple{Sequence_DoN, Any}) =
    inference_of_sequence(inner_sequence, inner_state)

"
Executes a list of inner sequences, in order.

If put in repeating mode, it will keep executing the list
  for as long as the first element of it can find a match.
Note that a repeating Ordered sequence could exit without even running once!
"
struct Sequence_Ordered <: AbstractSequence
    list::Vector{AbstractSequence}
    repeats::Bool
    inference::AllInference

    Sequence_Ordered(sequences, repeats=false, inference=AllInference()) = new(
        if sequences isa Vector{AbstractSequence}
            sequences
        else
            collect(AbstractSequence, sequences)
        end,
        convert(Bool, repeats),
        inference
    )
end
struct Sequence_Ordered_State
    current_i::Int
    current_state::Any
    current_any_matches_yet::Bool
    inference::AllInference
end
function start_sequence(s::Sequence_Ordered, grid::CellGrid{N}, outer_inference::AllInference, rng::PRNG) where {N}
    if isempty(s.list)
        return nothing
    else
        inf = AllInference(outer_inference, s.inference)
        return Sequence_Ordered_State(
            1, start_sequence(s.list[1], grid, inf, rng),
            false, inf
        )
    end
end
execute_sequence(s::Sequence_Ordered, grid::CellGrid, rng::PRNG, ::Nothing) = nothing
function execute_sequence(s::Sequence_Ordered, grid::CellGrid{N}, rng::PRNG,
                          iter::Sequence_Ordered_State) where {N}
    next_from_element = execute_sequence(s.list[iter.current_i], grid, rng, iter.current_state)
    if exists(next_from_element)
        return Sequence_Ordered_State(iter.current_i, next_from_element, true, iter.inference)
    elseif s.repeats && (iter.current_i == 1) && !iter.current_any_matches_yet
        # We failed to restart, so give up.
        return nothing
    elseif iter.current_i < length(s.list)
        # Recursively move on to the next iteration.
        return execute_sequence(s, grid, rng, Sequence_Ordered_State(
            iter.current_i + 1,
            start_sequence(s.list[iter.current_i + 1], grid, iter.inference, rng),
            false, iter.inference
        ))
    elseif s.repeats
        # Start anew.
        return execute_sequence(s, grid, rng, Sequence_Ordered_State(
            1, start_sequence(s.list[1], grid, iter.inference, rng),
            false, iter.inference
        ))
    else
        # We ran out of things to do.
        return nothing
    end
end
dsl_string(o::Sequence_Ordered) = string(
    "@block ",
    o.repeats ? "repeat " : "",
    """begin
    """,
    inference_exists(o.inference) ? string("    ", dsl_string(o.inference), "\n") : "",
    (string.(Ref("    "), dsl_string.(o.list), Ref("\n")))...,
    "    end"
)
inference_of_sequence(o::Sequence_Ordered, ::Nothing) = nothing
inference_of_sequence(o::Sequence_Ordered, iter::Sequence_Ordered_State) = inference_of_sequence(
    o.list[iter.current_i],
    iter.current_state
)

"Draws a rectangle of pixels, specified in UV space"
struct Sequence_DrawBox{N} <: AbstractSequence
    area::Bplus.BoxF{N}
    output_type::UInt8

    Sequence_DrawBox(area::Bplus.Box{N, <:Real}, output) where {N} = new{N}(
        convert(BoxF{N}, area),
        if output isa UInt8
            output
        elseif output isa Char
            CELL_CODE_BY_CHAR[output]
        elseif output isa Integer
            CELL_TYPES[output].code
        else
            error("Unsupported cell type: ", typeof(output))
        end
    )
end
function start_sequence(d::Sequence_DrawBox{N}, grid::CellGrid{M},
                        outer_inference::AllInference, rng::PRNG) where {N, M}
    if N != M
        error("Passed an ", M, "-dimensional grid into a Sequence_DrawBox{", N, "}")
    elseif N < 1
        error("N is ", N, "; must be at least 1-dimensional")
    end

    resolution = vsize(grid)
    resolutionF = convert(VecF{N}, resolution)
    b_min = convert(v2i, floor(min_inclusive(d.area) * resolutionF))
    b_max = convert(v2i, floor(max_exclusive(d.area) * resolutionF))

    (b_min, b_max) = minmax(b_min, b_max)
    b_min = clamp(b_min, one(v2i), resolution)
    b_max = clamp(b_max, one(v2i), resolution)
    return (
        b_min, b_max,
        b_min
    )
end
function execute_sequence(d::Sequence_DrawBox{N}, grid::CellGrid{N}, rng::PRNG,
                          (b_min, b_max, b_current)::NTuple{3, VecI{N}}) where {N}
    if b_current.x > b_max.x
        return nothing
    end

    grid[b_current] = d.output_type

    b_prev = b_current
    incremented_axis::Int = 1
    while incremented_axis <= N
        @set! b_current[incremented_axis] += 1
        if b_current[incremented_axis] > b_max[incremented_axis]
            @set! b_current[incremented_axis] = b_min[incremented_axis]
            incremented_axis += 1
        else
            break
        end
    end
    if incremented_axis > N
        return nothing
    end

    return (b_min, b_max, b_current)
end
function broadcast_sequence(seq::Sequence_DrawBox{N}, new_dims::Int) where {N}
    if N == new_dims
        return seq
    elseif N == 1
        return Sequence_DrawBox(
            BoxF{new_dims}(
                min=VecF{new_dims}(i->min_inclusive(seq.area).x),
                size=VecF{new_dims}(i->size(seq.area).x)
            ),
            seq.output_type
        )
    else
        error("Can't broadcast a ", N, "D box to ", new_dims, "D")
    end
end
dsl_string(b::Sequence_DrawBox{N}) where {N} = string(
    "@draw_box('", dsl_string(b.output_type), "', ",
    "min=", min_inclusive(b.area).data, ", ",
    "size=", size(b.area).data,
    ")"
)
inference_of_sequence(::Sequence_DrawBox, ::NTuple{3, VecI{N}}) where {N} = nothing

#TODO: Slice a grid to individual M-dimensional areas, and execute a sequence on each such area