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
Base.:(==)(a::T, b::T) where {T<:Union{AbstractMarkovOp, AbstractMarkovBias}} =
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

"
Ticks level 1 and 2 are very minor events
  that we usually don't want to waste any CPU cycles on checking.
"
const STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY = 3

struct ErrorCancelAlgo <: Exception end
Base.showerror(io::IO, ::ErrorCancelAlgo) = error("ErrorCancelAlgo should never go uncaught!")

"
Controls how an algorithm a) time-slices itself and b) keeps the user updated on tagged events.
You may edit these settings mid-run.
"
Base.@kwdef mutable struct TickSettings{MinCompileTimePriority}
    # The above type parameter removes low-priority ticks at compile time;
    #   this removes further ones with runtime checks.
    min_runtime_tick_priority::Int = STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY + 1

    # Only "new grid" and "canceled/completed" will happen; the rest will not.
    skip_most_tagged_events::Bool = false
    cancel_algo::Bool = false

    # A hint to algorithm logic that the visual output matters;
    #   without this there are some clever optimizations that re-order events.
    #
    # Implemented as a stack so that specific sequences can push and pop their own hint.
    # An empty stack represents 'false'.
    animated::Vector{Bool} = preallocated_vector(Bool, 2)
end
TickSettings(; kw...) = TickSettings{STANDARD_MIN_COMPILE_TIME_TICK_PRIORITY}(; kw...)
is_animated(t::TickSettings)::Bool = !isempty(t.animated) && t.animated[end]

const AlgoCommsChannel = Channel{Union{Int, Symbol, CellGrid{NGrid}}}

mutable struct AlgoState{NGrid, TGrid<:CellGrid{NGrid},
                         MinCompileTimeTickPriority}
    grid::TGrid
    rng::PRNG
    allocator::AbstractMarkovAllocator

    data_store::Dict{Symbol, Any}()

    channel::AlgoCommsChannel
    # (the inclusion of Symbol means that the Ints get boxed, however
    #    Julia keeps an array of pre-boxed small Ints so we shouldn't be affected --
    #    see https://github.com/JuliaLang/julia/issues/62598)

    ticking::TickSettings{MinCompileTimeTickPriority}
end

@inline markov_algo_tick(s::AlgoState{NG, TG, MinP}, p::Int) where {NG, TG, MinP} = if p < MinP
    nothing
elseif s.ticking.cancel_algo
    throw(ErrorCancelAlgo())
elseif (p <= s.min_runtime_priority) || s.insta_finish
    nothing
else
    put!(channel, M)
    @bp_check(take!(channel) == 0, "You should write 0 in response to every tick")
    nothing
end
function markov_algo_tick(s::AlgoState, tag::Symbol)
    if !s.ticking.skip_most_tagged_events || in(tag, IMPORTANT_BUILTIN_TAGS)
        put!(s.channel, tag)
        if tag == TAG_NEW_GRID
            put!(s.channel, s.grid)
        end
        @bp_check(take!(channel) == 0, "You should write 0 in response to every tick")
    end
    if s.ticking.cancel_algo
        throw(ErrorCancelAlgo())
    end
    return nothing
end

"Reallocates the algorithm grid and issues a tagged event for it"
function markov_algo_new_grid(new_state_setup_fn,
                              state::AlgoState{NGrid, TGrid},
                              new_size::NTuple{NGrid, Int}
                             )::CellGrid{NGrid}
    new_grid = markov_allocator_acquire_array(state.allocator, new_size, UInt8)
    new_state_setup_fn(state.grid, new_grid)

    markov_allocator_release_array(state.allocator, state.grid)
    state.grid = new_grid

    markov_algo_tick(state, TAG_NEW_GRID)
    return new_grid
end

"
Runs the algorithm as a coroutine on a separate task,
   writing ticks/events through a Channel and waiting for you to write `zero(Int)`
   in response to each one.
Returns the Channel you must use to communicate with it.

`seeds` should be a `Real` number or small enumeration of `Real` numbers.
"
function markov_algo_run(algo::MarkovAlgorithm,
                         initial::Union{CellGrid, Tuple{Vararg{Integer}}, VecT{<:Integer}},
                         ticking::TickSettings = TickSettings()
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

    rng::PRNG = (seeds isa Real) ? PRNG(seeds) : PRNG(seeds...)
    data_store = Dict{Symbol, Any}()
    return AlgoCommsChannel(0) do ch
        initial_grid = markov_allocator_acquire_array(allocator, initial_size, UInt8)
        if initial isa CellGrid
            copyto!(initial_grid, initial)
        else
            fill!(initial_grid, algo.initial_fill)
        end
        algo_state = AlgoState(initial_grid, rng, allocator, data_store, ch, ticking)

        try
            @logic_logln "Starting algorithm of size " size(algo_state.grid)
            markov_algo_tick(algo_state, TAG_NEW_GRID)
            markov_algo_tick(algo_state, TAG_ALGO_STARTING)

            @logic_tab_in()
            for op in algo.sequence
                @logic_logln "Starting a top-level op..."
                @logic_tab_in()
                markov_algo_run(op, algo, algo_state, (), ())
                @logic_tab_out()
            end
            @logic_logln "Finishing algorithm of size " size(algo_state.grid)
            @logic_tab_out()
            markov_algo_tick(algo_state, TAG_ALGO_COMPLETED)
        catch e
            if e isa ErrorCancelAlgo
                @logic_logln "User canceled algorithm!"
                markov_algo_tick(algo_state, TAG_ALGO_CANCELED)
            else
                rethrow()
            end
        finally
            markov_allocator_release_array(algo_state.allocator, algo_state.grid)
        end
    end
end

"
Runs through the given algorithm (as its channel returned from `markov_algo_run()`).
On each tagged event, sends it along with the current grid into your lambda.
  You can also provide a lambda for normal ticks.
"
function markov_algo_complete(process_tagged_event,
                              ch::AlgoCommsChannel,
                              current_grid::Optional{CellGrid{NGrid}} = nothing,
                              process_tick = ((priority::Int, grid::CellGrid{NGrid}) -> nothing)
                             )::Nothing where {NGrid}
    while true
        next_msg = take!(ch)
        if next_msg isa Int
            process_tick(next_msg, current_grid)
        elseif next_msg isa Symbol
            if next_msg == TAG_NEW_GRID
                current_grid = take!(ch)
            end
            process_tagged_event(next_msg, current_grid)
            if next_msg == TAG_ALGO_COMPLETED
                break
            end
        else
            error("Unhandled: ", typeof(next_msg))
        end
    end

    return nothing
end

markov_algo_to_string(a::MarkovAlgorithm) = dsl_string(a)