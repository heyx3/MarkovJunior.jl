"
A no-op that emits a tagged event, pausing the algorithm so users can inspect the grid.
If mutable (tag name has a leading '!'), users may also modify the grid before resuming.
"
struct MarkovOpEvent <: AbstractMarkovOp
    tag::Symbol
end

function markov_algo_run(op::MarkovOpEvent,
                         algo::MarkovAlgorithm, algo_state::AlgoState,
                         inherited_bias_tuple::Tuple{Vararg{AbstractMarkovBias}},
                         inherited_bias_state_tuple::Tuple,
                         grid::CellGrid{NGrid} = algo_state.grid
                        )::Tuple{OpChangedGrid, typeof(inherited_bias_state_tuple)} where {NGrid}
    # Read-only events are cosmetic.
    if !is_mutable_event_tag(op.tag)
        markov_algo_tick(algo_state, op.tag)
        return (nothing, inherited_bias_state_tuple)
    # Skipped mutable events count as failing to change anything.
    elseif !is_tagged_event_emitted(algo_state.ticking, op.tag)
        markov_algo_tick(algo_state, op.tag)
        return (false, inherited_bias_state_tuple)
    end

    # Snapshot the grid so we can see what the user changed.
    return markov_allocator_with_array(algo_state.allocator, size(grid), UInt8) do _previous_grid
        previous_grid::Array{UInt8, NGrid} = _previous_grid
        copyto!(previous_grid, grid)
        markov_algo_tick(algo_state, op.tag)

        # Find the bounds of the user's changes.
        grid_size = convert(Vec{NGrid, Int32}, vsize(grid))
        changed_min = grid_size + one(Int32)
        changed_max = zero(Vec{NGrid, Int32})
        for pixel in one(Vec{NGrid, Int32}):grid_size
            if grid[pixel] != previous_grid[pixel]
                changed_min = min(changed_min, pixel)
                changed_max = max(changed_max, pixel)
            end
        end
        changed_area = Box(min=changed_min, max=changed_max)
        if is_empty(changed_area)
            return (false, inherited_bias_state_tuple)
        end

        changed_slice = ntuple(i -> changed_min[i]:changed_max[i], Val(NGrid))
        return (true, markov_bias_update.(
            inherited_bias_tuple, inherited_bias_state_tuple,
            Ref(algo), Ref(algo_state),
            Ref(changed_area), Ref(@view previous_grid[changed_slice...])
        ))
    end
end


dsl_format(e::MarkovOpEvent) = "@event $(e.tag)"

function parse_markovjunior_op(::Val{Symbol("@event")}, inputs::MacroParserInputs,
                               loc, expr_args, expr_line)
    if length(expr_args) != 1
        raise_parse_error(loc, inputs, "Expected exactly one argument (the event's name); got ", length(expr_args))
    end

    (name, is_mutable) = if expr_args[1] isa Symbol
        (expr_args[1], false)
    elseif @capture(expr_args[1], !n_Symbol)
        (n, true)
    else
        raise_parse_error(loc, inputs, "Expected an event name like `MyEvent` or `!MyEvent`; got `", expr_args[1], "`")
    end
    if name in BUILTIN_TAGS
        raise_parse_error(loc, inputs, "`", name, "` is a built-in tag and can't be used by @event")
    end

    return MarkovOpEvent(is_mutable ? Symbol(MUTABLE_EVENT_TAG_PREFIX, name) : name)
end
