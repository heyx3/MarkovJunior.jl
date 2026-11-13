#######################
#  Logging

"Enable detailed algorithm logging by redefining this to return `true`"
log_logic() = false

LOG_LOGIC_TAB::String = ""
macro logic_tab_in()
    MJ = @__MODULE__
    return quote
        if $log_logic()
            $MJ.LOG_LOGIC_TAB::String *= "    "
        end
        nothing
    end
end
macro logic_tab_out()
    MJ = @__MODULE__
    return quote
        if $log_logic()
            $MJ.LOG_LOGIC_TAB = $MJ.LOG_LOGIC_TAB[5:end]
        end
        nothing
    end
end

macro logic_log(args...)
    MJ = @__MODULE__
    return quote
        if $log_logic()
            foreach(tuple( $(esc.(args)...)) ) do arg
                if arg isa String
                    arg = replace(arg, "\n"=>"\n$($MJ.LOG_LOGIC_TAB)")
                elseif arg == '\n'
                    arg = "\n$($MJ.LOG_LOGIC_TAB)"
                else
                    arg = replace(string(arg), "\n"=>"\n$($MJ.LOG_LOGIC_TAB)")
                end
                print(stderr, arg)
            end
        end
        nothing
    end
end
macro logic_logln(args...)
    MJ = @__MODULE__
    return :(
        $MJ.@logic_log(
            $(esc.(args)...),
            '\n'
        )
    )
end


##################
#  Types

abstract type AbstractMarkovOp end
abstract type AbstractMarkovBias end

# Equality between ops and between biases is important for testing.
# We set a standard implementation of checking their precise types
#   and comparing their individual fields with == (by default Julia would use === on the fields).
Base.:(==)(a::AbstractMarkovOp, b::AbstractMarkovOp) =
    (typeof(a) == typeof(b)) && all(
        f -> getfield(a, f) == getfield(b, f),
        fieldnames(typeof(a))
    )
#


struct MarkovAlgorithm
    initial_fill::UInt8

    fixed_dimension::Optional{Int}
    min_dimension::Int # Will be 1 if there aren't any restrictions;
                       #   will be fixed_dimension if it exists

    sequence::Vector{AbstractMarkovOp}

    # Pragmas should be added to the top of the algorithm block,
    #   and are stored here both chronologically and in a Dict lookup.
    pragmas_chronological::Vector{Pair{Symbol, Int}} # Each key to its index -- `args = pragmas_map[key][value]`
    pragmas_map::Dict{Symbol, Vector{Vector{Any}}} # Each key maps to a list of its arguments each time it was invoked, in chronological order

    # User functions for calculations, from `@calc` statements anywhere in the algorithm.
    # Calls look these up by name whenever they run, so users may add or replace them between runs
    #    (e.g. `push!(algo.calc_macros, calc_macro_parse("room_width() = 8"))`) to parameterize the algorithm.
    calc_macros::OrderedDict{Symbol, CalcInterpreter.CalcMacro}

    # "Add-ons" can be added by the user after parsing is finished,
    #    to inject per-instance parameters.
    add_ons::Dict{Symbol, Any}
end
Base.:(==)(a::MarkovAlgorithm, b::MarkovAlgorithm) = (
    a.initial_fill == b.initial_fill &&
    a.fixed_dimension == b.fixed_dimension &&
    a.min_dimension == b.min_dimension &&
    a.sequence == b.sequence &&
    a.pragmas_chronological == b.pragmas_chronological &&
    a.pragmas_map == b.pragmas_map &&
    a.calc_macros == b.calc_macros &&
    a.add_ons == b.add_ons
)


"The tag that is emitted right after the algorithm starts, allowing you to decide when to actually let it run"
const TAG_ALGO_STARTING = :ALGORITHM_STARTING
"The tag that's used whenever an Op replaces the grid (including before the first tick); this tag is always followed by a reference to the new grid"
const TAG_NEW_GRID = :ALGORITHM_REPLACED_GRID
"The tag for a completed algorithm run"
const TAG_ALGO_COMPLETED = :ALGORITHM_COMPLETED
"The tag for a canceled algorithm run"
const TAG_ALGO_CANCELED = :ALGORITHM_CANCELED

const BUILTIN_TAGS = Set{Symbol}([
    TAG_ALGO_STARTING, TAG_NEW_GRID,
    TAG_ALGO_COMPLETED, TAG_ALGO_CANCELED
])
const IMPORTANT_BUILTIN_TAGS = Set{Symbol}([
    TAG_NEW_GRID,
    TAG_ALGO_COMPLETED, TAG_ALGO_CANCELED
])

"The leading character of a tag which marks it as a mutable event (users may modify the grid while it's paused)"
const MUTABLE_EVENT_TAG_PREFIX = '!'
"Checks whether the given tag represents a mutable event, without allocating a String"
@inline is_mutable_event_tag(tag::Symbol)::Bool =
    unsafe_load(Base.unsafe_convert(Ptr{UInt8}, tag)) == UInt8(MUTABLE_EVENT_TAG_PREFIX)

"
Ticks level 1 and 2 are very minor events
  that we usually don't want to waste any CPU cycles on checking.
"
const STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY = 3
"The tick priority for the end of a normal Op, like @rewrite"
const STANDARD_END_OF_OP_TICK_PRIORITY = 4

"
Internal error type which signals that the algorithm has been canceled.

This is caught at the top of the algorithm loop; all Ops/Biases/etc.
  must not prevent it from propagating upward!
"
struct ErrorCancelAlgo <: Exception end
Base.showerror(io::IO, ::ErrorCancelAlgo) = error("ErrorCancelAlgo should never go uncaught!")

"""
Which tagged events get skipped instead of being sent to users.
# The important built-in ones (new grid, and canceled/completed) are never skipped.
#   * `none_skippable`: all tagged events happen.
#   * `immutable_skippable`: only mutable events (see `is_mutable_event_tag()`) happen.
#   * `all_skippable`: no other tagged events happen.
"""
@bp_enum(SkipTaggedEvents,
    none_skippable,
    immutable_skippable,
    all_skippable
)

"
Controls ticking, tagged events, and cancellation.
You may edit any of these settings mid-run.
"
mutable struct MarkovTickSettings{MinCompileTimePriority}
    # The above type parameter removes low-priority ticks at compile time;
    #   this removes further ones with runtime checks.
    min_runtime_priority::Int

    skip_tagged_events::E_SkipTaggedEvents
    cancel_algo::Bool

    # A hint to algorithm logic that the visual output matters;
    #   without this there are some clever optimizations that re-order events.
    #
    # Implemented as a stack so that specific sequences can push and pop their own hint.
    # An empty stack represents 'false'.
    animated::Vector{Bool}

    MarkovTickSettings(min_tick_priority::Integer = STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY
                       ;
                       skip_tagged_events::E_SkipTaggedEvents = SkipTaggedEvents.none_skippable
                      ) = new{min(convert(Int, min_tick_priority), STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY)}(
        convert(Int, min_tick_priority),
        skip_tagged_events,
        false, preallocated_vector(Bool, 2)
    )
end
is_animated(t::MarkovTickSettings)::Bool = !isempty(t.animated) && t.animated[end]

"Checks whether the given tagged event would be sent to users, rather than skipped"
function is_tagged_event_emitted(t::MarkovTickSettings, tag::Symbol)::Bool
    if (t.skip_tagged_events == SkipTaggedEvents.none_skippable) || in(tag, IMPORTANT_BUILTIN_TAGS)
        return true
    elseif t.skip_tagged_events == SkipTaggedEvents.immutable_skippable
        return is_mutable_event_tag(tag)
    elseif t.skip_tagged_events == SkipTaggedEvents.all_skippable
        return false
    else
        error("Unhandled: ", t.skip_tagged_events)
    end
end

const AlgoCommsChannel = Channel{Union{Int, Symbol, CellGrid}}

mutable struct AlgoState{MinCompileTimeTickPriority}
    grid::CellGrid
    rng::PRNG
    allocator::AbstractMarkovAllocator

    # Custom data for the whole run (see docs/add_ons.md).
    # It's never reset, even when the grid is reallocated or a scoped Op swaps in a view;
    #   grid-dependent data should check `grid_generation` and `grid_view_id`.
    data_store::Dict{Symbol, Any}

    channel::AlgoCommsChannel
    # (the inclusion of Symbol means that the Ints get boxed, however
    #    Julia keeps an array of pre-boxed small Ints so we shouldn't be affected --
    #    see https://github.com/JuliaLang/julia/issues/62598)

    ticking::MarkovTickSettings{MinCompileTimeTickPriority}

    # The state of every active bias, in the same order as the inherited-bias tuples.
    bias_state_stack::Vector{Any}

    # Incremented every time the grid is reallocated (which may also replace bias states).
    grid_generation::Int
    # How many scoped Ops (like `@within_box`) deep we currently are.
    # As long as it's >0, the grid is a slice of a larger array (which also means it can't be resized).
    scope_depth::Int
    # Unique identifier for the grid view currently being used within a larger grid.
    # A value of 0 means the whole original grid is being modified.
    # You can uniquely identify an algorithm scope with (grid_generation, grid_view_id).
    grid_view_id::Int

    next_grid_view_id::Int

    # How many cells of the current grid are actually being worked on.
    # Usually it's the whole grid, but some ops (`@within_regions`) work on a subset of that.
    grid_area::Int
end


#################################
##  Helper functions for Ops

"
Implements our internal coroutine logic,
  potentially pausing the calling thread until the user resumes us.

Also does other bookkeeping like checking for algorithm cancellation.
"
@inline markov_opalgo_tick(s::AlgoState{MinP}, p::Int) where {MinP} = if p < MinP
    #NOTE: it's force-inlined so that this first condition disappears at compile-time
    nothing
elseif s.ticking.cancel_algo
    throw(ErrorCancelAlgo())
elseif p < s.ticking.min_runtime_priority
    nothing
else
    put!(s.channel, p)
    @bp_check(take!(s.channel) === 0, "You should write 0 in response to every tick")
    nothing
end

"
Implements our internal coroutine logic for tagged events,
  usually pausing the calling thread until the user resumes us.

Also does some other bookkeeping like checking for algorithm cancellation.
"
function markov_opalgo_tick(s::AlgoState, tag::Symbol)
    if is_tagged_event_emitted(s.ticking, tag)
        put!(s.channel, tag)
        if tag == TAG_NEW_GRID
            put!(s.channel, s.grid)
        end
        @bp_check(take!(s.channel) === 0, "You should write 0 in response to every tick")
    end
    if s.ticking.cancel_algo
        throw(ErrorCancelAlgo())
    end
    return nothing
end


"
Reallocates the algorithm grid, lets the current biases react, then issues a tagged event for it.
Returns the new grid, a.k.a. `state.grid`,
  as well as the new bias states (which may have different types!).

The `bias_states` input is only passed for type-stability
  and must match the bias states stored in `AlgoState`.
"
function markov_opalgo_new_grid(new_state_setup_fn,
                                algo::MarkovAlgorithm, state::AlgoState,
                                new_size::Union{NTuple{NGrid, Integer}, Vec{NGrid, <:Integer}},
                                biases::NTuple{NBiases, AbstractMarkovBias},
                                bias_states::NTuple{NBiases, Any}
                               )::Tuple{CellGrid{NGrid}, NTuple{NBiases, Any}} where {NGrid, NBiases}
    if state.scope_depth > 0
        error("The grid can't be reallocated inside a scoped Op (`@within_box`, `@within_region`, etc)! ",
                "This is normally caught at parse-time using `markov_op_reallocates_grid()`, ",
                "so there may be a bug with the current Op.")
    end
    @markovjunior_assert(length(state.bias_state_stack) == NBiases &&
                           all(i -> state.bias_state_stack[i] === bias_states[i], 1:NBiases),
                         "You must pass all inherited biases to markov_opalgo_new_grid()")

    # Dispatch on the grid and allocator type:
    return ((old_grid, alloc) -> begin
        # On failure, the old grid should remain the current one.
        new_grid::CellGrid{NGrid} = markov_allocator_acquire_array(alloc, Tuple(new_size), UInt8)
        try
            new_state_setup_fn(old_grid, new_grid)
        catch
            markov_allocator_release_array(alloc, new_grid)
            rethrow()
        end

        state.grid = new_grid
        state.grid_area = length(new_grid)
        state.grid_generation += 1
        try
            # Replace each bias's state as soon as it's ready, so cleanup always sees the live one.
            foreach(1:NBiases) do i
                state.bias_state_stack[i] = markov_bias_new_grid(biases[i], bias_states[i],
                                                                 algo, state, old_grid)
            end

            markov_opalgo_tick(state, TAG_NEW_GRID)
            return (new_grid, markov_opbias_current_states(state))
        finally
            markov_allocator_release_array(alloc, old_grid)
        end
    end)(state.grid, state.allocator)
end


"
Sets up the algorithm state to run on a view within the original grid,
  executes your lambda, then restores the full grid.
Manages the biases properly, and ends by returning your lambda's return value.
Your lambda signature must be
  `(inner_biases::NTuple{N, AbstractMarkovBias}, inner_bias_states::NTuple{N, Any}) -> OpChangedGrid`.

You set up the view, tell us its `area` within the larger one,
  and tell us about biases (e.g. whether the inner logic inherits them).

The outer bias states are left alone while your lambda runs,
  then updated all at once with the area's changes (if it reports any).
If the inner logic inherits the biases, they get separate states which only see the view,
  and which are cleaned up afterwards.
Either way, nothing replaces bias states while inside the scope, so `grid_generation` stays the same.

You may also provide callbacks to prepare and restore data within the grid slice,
  in a way that the outer biases can't see, e.g. to mask out certain pixels from the scoped run.
Inherited biases are initialized after `prepare_view()` and cleaned up before `restore_view()`.

Within your lambda, Ops are forbidden from resizing the grid.
This is enforced already, but you should implement `markov_op_reallocates_grid()`
  to catch the problem at parse-time.

Note that `restore_view` is called even if the inner Ops throw,
  though the exception is rethrown afterward.

`n_cells` is how many of the view's cells are actually being worked on (see `AlgoState.grid_area`),
  if it's not all of them (e.g. a `@within_regions` region, whose bounding box also holds dummy cells).
"
function markov_opalgo_on_subgrid(to_do,
                                  algo::MarkovAlgorithm, algo_state::AlgoState,
                                  sub_grid::CellGrid{N}, area::BoxI{N},
                                  inherit_biases::Bool,
                                  outer_biases::Tuple{Vararg{AbstractMarkovBias}},
                                  outer_bias_states::Tuple,
                                  prepare_view = () -> nothing,
                                  restore_view = () -> nothing
                                  ;
                                  n_cells::Int = length(sub_grid)
                                 ) where {N}
    @markovjunior_assert(vsize(sub_grid) == size(area),
                         "sub_grid size is ", vsize(sub_grid), " but reported box is ", area)
    @markovjunior_assert(length(algo_state.bias_state_stack) == length(outer_biases),
                         "Bias stack has ", length(algo_state.bias_state_stack),
                           " states, but there are ", length(outer_biases), " outer biases")
    outer_grid = algo_state.grid
    alloc = algo_state.allocator

    # Hide the outer bias states, so the inner Ops get their own stack.
    # The outer states are never replaced by the scope, only updated afterwards;
    #    so `grid_generation` doesn't change and nested scopes are invisible to the Ops running them.
    hidden_states::Vector{Any} = markov_allocator_acquire_array(alloc, (), Any)
    append!(hidden_states, algo_state.bias_state_stack)
    empty!(algo_state.bias_state_stack)

    # Remember the scope's original pixels, so the outer biases can be updated afterwards.
    previous_pixels = if isempty(outer_biases)
        nothing
    else
        a::Array{UInt8, N} = markov_allocator_acquire_array(alloc, size(sub_grid), UInt8)
        copyto!(a, sub_grid)
        a
    end

    outer_view_id = algo_state.grid_view_id
    outer_grid_area = algo_state.grid_area
    try
        algo_state.grid = sub_grid
        algo_state.scope_depth += 1
        algo_state.next_grid_view_id += 1
        algo_state.grid_view_id = algo_state.next_grid_view_id
        algo_state.grid_area = n_cells

        result::OpChangedGrid = try
            prepare_view()
            if inherit_biases
                # The outer biases get separate states which only see the view.
                inner_bias_states = map(b -> markov_bias_initialize(b, algo, algo_state), outer_biases)
                append!(algo_state.bias_state_stack, inner_bias_states)
                to_do(outer_biases, inner_bias_states)
            else
                to_do((), ())
            end
        finally
            # Clean up any inner states while they can still see the view.
            # Inner Ops' own biases have already been popped off the stack by now.
            if inherit_biases
                for (bias, inner_state) in zip(outer_biases, algo_state.bias_state_stack)
                    markov_bias_cleanup(bias, inner_state, algo, algo_state)
                end
            end
            restore_view()

            algo_state.scope_depth -= 1
            algo_state.grid_view_id = outer_view_id
            algo_state.grid_area = outer_grid_area
            algo_state.grid = outer_grid

            empty!(algo_state.bias_state_stack)
            append!(algo_state.bias_state_stack, hidden_states)
        end

        if exists(previous_pixels) && (result === true)
            markov_opbias_update_all!(outer_biases, outer_bias_states,
                                      algo, algo_state,
                                      area, previous_pixels)
        end
        return result
    finally
        markov_allocator_release_array(alloc, hidden_states)
        exists(previous_pixels) && markov_allocator_release_array(alloc, previous_pixels)
    end
end

######################


#######################
##  Public interface

"
Runs the algorithm as a coroutine on a separate task,
   returning the Channel you use to communicate with it.
`seeds` should be a `Real` number or small enumeration of `Real` numbers.

You can use high-level tick helpers with `markov_algo_next()` and `markov_algo_complete()`,
  or follow the below low-level protocol:

1. `put!(channel, zero(Int))` to advance the algorithm
2. `take!(channel)` to get the next event
3. Event is either an integer for the Priority of the tick that just happened, or a `Symbol` tagged event.
    * If the event is `TAG_NEW_GRID`, read the new grid with `take!(channel)::CellGrid`.
    * If the event is `TAG_ALGO_COMPLETED` or `TAG_ALGO_CANCELED`,
save the most recent grid data somewhere and then write `zero(Int)` to clean up the algorithm.
4. Go back to 1 to dispatch the next 'tick'.
5. Cancel the algorithm anytime with `close(channel)`, or by setting `ticking.cancel_algo = true`.
"
function markov_algo_run(algo::MarkovAlgorithm,
                         initial::Union{CellGrid, Tuple{Vararg{Integer}}, VecT{<:Integer}},
                         ticking::MarkovTickSettings = MarkovTickSettings()
                         ;
                         seeds = rand(UInt32),
                         allocator::AbstractMarkovAllocator = MarkovAllocatorHeapReused()
                        )::AlgoCommsChannel
    initial_size::Tuple{Vararg{Int}} = if initial isa CellGrid
        size(initial)
    elseif initial isa Vec
        convert.(Ref(Int), initial.data)
    elseif initial isa Tuple
        convert.(Ref(Int), initial)
    else
        error("Unhandled: ", typeof(initial))
    end

    if exists(algo.fixed_dimension) && (length(initial_size) != algo.fixed_dimension)
        error("Can't start a ", length(initial_size), "D MarkovJunior run with a ",
              algo.fixed_dimension, "D algorithm")
    end
    if length(initial_size) < algo.min_dimension
        error("Can't start a ", length(initial_size), "D MarkovJunior run with a ",
              algo.min_dimension, "+D algorithm")
    end

    rng::PRNG = (seeds isa Real) ? PRNG(seeds) : PRNG(seeds...)
    data_store = Dict{Symbol, Any}()
    return AlgoCommsChannel(0) do ch
        initial_grid = markov_allocator_acquire_array(allocator, initial_size, UInt8)
        if initial isa CellGrid
            copyto!(initial_grid, initial)
        else
            fill!(initial_grid, algo.initial_fill)
        end
        algo_state = AlgoState(initial_grid, rng, allocator, data_store, ch, ticking, Any[ ], 0, 0, 0, 0,
                               length(initial_grid))

        try
            @logic_logln "Starting algorithm of size " size(algo_state.grid)
            @bp_check(take!(ch) === 0)
            markov_opalgo_tick(algo_state, TAG_NEW_GRID)
            markov_opalgo_tick(algo_state, TAG_ALGO_STARTING)

            @logic_tab_in()
            for op in algo.sequence
                @logic_logln "Starting a top-level op..."
                @logic_tab_in()
                markov_algo_run(op, algo, algo_state, (), ())
                @logic_tab_out()
            end
            @logic_logln "Finishing algorithm of size " size(algo_state.grid)
            @logic_tab_out()
            markov_opalgo_tick(algo_state, TAG_ALGO_COMPLETED)
        catch e
            if isa(e, ErrorCancelAlgo) || isa(e, InvalidStateException)
                @logic_logln "User canceled algorithm!"
                # If canceled cleanly (without closing the channel),
                #   issue a 'Canceled' tagged event through the channel.
                (e isa ErrorCancelAlgo) && isopen(ch) && markov_opalgo_tick(algo_state, TAG_ALGO_CANCELED)
            else
                rethrow()
            end
        finally
            markov_allocator_release_array(algo_state.allocator, algo_state.grid)
        end
    end
end

"
Runs one 'step' of the algorithm, and returns what that step was:

* The priority of the tick that just completed
* The tagged event that just happened
* The newly-allocated grid, on a `TAG_NEW_GRID` event

If the output is `TAG_ALGO_COMPLETED` or `TAG_ALGO_CANCELED`,
  you should call `markov_algo_cleanup()` to close out the algorithm run.
"
function markov_algo_next(channel::AlgoCommsChannel)::Union{Int, Symbol, CellGrid}
    put!(channel, zero(Int))
    msg = take!(channel)
    if msg isa Int
        return msg
    elseif msg == TAG_NEW_GRID
        return take!(channel)::CellGrid
    elseif msg isa Symbol
        return msg
    else
        error("Unhandled: ", typeof(msg), "(", msg, ")")
    end
end
function markov_algo_cleanup(channel::AlgoCommsChannel)
    if isopen(channel)
        put!(channel, 0)
        close(channel)
    end
    return nothing
end

"
Executes the rest of the given algorithm (using its comms channel, returned from `markov_algo_run()`).
Runs your lambda on the final `CellGrid` representing the state after the algorithm is done,
  keeping in mind that this grid is deallocated as soon as your lambda exits.
Returns the output of your lambda.

Optionally uses another lambda for processing each tick and tagged event,
   of the form `(tick::Union{Int, Symbol}, current_grid::CellGrid) -> nothing`.

If you've already partially run this algorithm, you should supply the most recent allocated grid,
   or else the grid passed to your lambdas will be `Nothing` until the next grid reallocation.

If you for some reason cancel the algorithm yourself while this function is executing,
   it will invoke your lambda with `nothing` for the grid instance.
"
function markov_algo_complete(process_final_state,
                              ch::AlgoCommsChannel,
                              current_grid::Optional{CellGrid} = nothing,
                              process_ticks = ((::Union{Int, Symbol}, ::CellGrid) -> nothing)
                             )
    try
        while true
            put!(ch, zero(Int))
            next_msg = take!(ch)

            if next_msg isa Int
                process_ticks(next_msg, current_grid)
            elseif next_msg isa Symbol
                if next_msg == TAG_NEW_GRID
                    current_grid = take!(ch)
                elseif next_msg in (TAG_ALGO_COMPLETED, TAG_ALGO_CANCELED)
                    result = process_final_state(current_grid)
                    put!(ch, zero(Int))
                    return result
                else
                    process_ticks(next_msg, current_grid)
                end
            else
                error("Unhandled: ", typeof(next_msg))
            end
        end
    catch e
        # If the user canceled the algo, catch and process that.
        # The grid is already potentially deallocated, so pass `nothing` in place of it.
        if isa(e, InvalidStateException)
            return process_final_state(nothing)
        else
            rethrow()
        end
    end

    return nothing
end

markov_algo_to_string(a::MarkovAlgorithm) = dsl_format(a)
