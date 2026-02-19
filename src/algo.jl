#######################
#  Logging

"Redefine this to return `true` to enable detailed logging"
log_logic() = false

LOG_LOGIC_TAB::String = ""
@inline logic_tab_in() = log_logic() && (global LOG_LOGIC_TAB *= "  ")
@inline logic_tab_out() = log_logic() && (global LOG_LOGIC_TAB = LOG_LOGIC_TAB[2:end])

@inline logic_log(args...) = if log_logic()
    foreach(args) do arg
        if arg isa String
            arg = replace(arg, "\n"=>"\n$LOG_LOGIC_TAB")
        elseif arg == '\n'
            arg = "\n$LOG_LOGIC_TAB"
        else
            arg = replace(string(arg), "\n"=>"\n$LOG_LOGIC_TAB")
        end
        print(stderr, arg)
    end
else
    nothing
end
@inline logic_logln(args...) = logic_log(args..., '\n')


##################
#  Allocators

"
Some kind of memory allocator, mainly used for the grids.
Default behavior of its interface is to allocate things from the heap
   and let the GC handle them on release.
"
abstract type AbstractMarkovAllocator end
markov_allocator_acquire_array(allocator::AbstractMarkovAllocator,
                               size::Tuple{Vararg{Integer}},
                               T::Type
                              )::AbstractArray{T, length(size)} where {T} = Array{T, length(size)}(undef, size)
markov_allocator_acquire_ordered_set(allocator::AbstractMarkovAllocator,
                                     T::Type
                                    )::OrderedSet{T} = OrderedSet{T}()
markov_allocator_release_array(allocator::AbstractMarkovAllocator, data::AbstractArray) = nothing
markov_allocator_release_ordered_set(allocator::AbstractMarkovAllocator, s::OrderedSet) = nothing
Base.close(::AbstractMarkovAllocator) = nothing

"An allocator for MarkovJunior that just goes through the heap"
struct MarkovAllocatorHeap <: AbstractMarkovAllocator end


##################
#  Types

abstract type AbstractMarkovOp end
abstract type AbstractMarkovBias end

struct MarkovAlgorithm
    initial_fill::UInt8

    fixed_dimension::Optional{Int}
    min_dimension::Int # Will be 1 if there aren't any restrictions

    pragmas::Vector{Pair{Symbol, Vector{Any}}} # Nonstandard commands meant to be parsed by the user of this package
    sequence::Vector{AbstractMarkovOp}
end


"Information about the algorithm state relevant to a Bias struct"
struct MarkovBiasContext
    allocator::AbstractMarkovAllocator
end

"
Information about the algorithm state relevant to an Op.
Inherits properties from `MarkovBiasContext`.
"
struct MarkovOpContext
    inherited_biases_type_stable::Tuple{Vararg{<:AbstractMarkovBias}}
    bias_context::MarkovBiasContext
end
MarkovOpContext(algo_state) = MarkovOpContext(  #NOTE: algo state can't be typed because it's defined below
    markov_algo_flatten_biases(algo_state),
    algo_state.bias_context
)
Base.propertynames(::MarkovOpContext) = propertynames(MarkovOpContext)
Base.propertynames(::Type{MarkovOpContext}) = (
    propertynames(MarkovBiasContext)...,
    fieldnames(MarkovOpContext)...
)
@inline Base.getproperty(c::MarkovOpContext, n::Symbol) = if n in fieldnames(MarkovOpContext)
    getfield(c, n)
else
    getproperty(getfield(c, :bias_context), n)
end

"Inherits properties from `MarkovOpContext`, though they aren't directly settable"
mutable struct MarkovAlgoState{N}
    grid::Base.RefValue{CellGridConcrete{N}} # Wrapped in a Ref for implementation reasons
    n_iterations::Int

    op_idx::Int
    op_state::Any
    op_context::MarkovOpContext # Inner sequences can modify this field

    # Sequences can push new biases onto this stack, and they all go into any ops inside that stack.
    bias_chain::Stack{Vector{AbstractMarkovBias}}

    rng::PRNG

    buffer_iter_count::Base.RefValue{Optional{Int}}
end
Base.propertynames(c::MarkovAlgoState) = propertynames(typeof(c))
Base.propertynames(T::Type{<:MarkovAlgoState}) = (
    propertynames(MarkovOpContext)...,
    fieldnames(T)...
)
@inline Base.getproperty(s::MarkovAlgoState, n::Symbol) = if n in fieldnames(typeof(s))
    getfield(s, n)
else
    getproperty(getfield(s, :op_context), n)
end


##################
#  Interface

function markov_algo_start(algo::MarkovAlgorithm,
                           initial_size::Union{Tuple{Vararg{Integer}}, VecT{<:Integer}},
                           seeds::Union{Real, Tuple{Vararg{Real}}, Vec} = rand(UInt32)
                           ;
                           allocator::AbstractMarkovAllocator = MarkovAllocatorHeap()
                          )::MarkovAlgoState
    if exists(algo.fixed_dimension) && (length(initial_size) != algo.fixed_dimension)
        error("Can't start a ", length(initial_size), "D MarkovJunior run with a ",
              algo.fixed_dimension, "D algorithm")
    end

    grid = markov_allocator_acquire_array(allocator, initial_size, UInt8)
    fill!(grid, algo.initial_fill)

    rng = (seeds isa Real) ? PRNG(seeds) : PRNG(seeds...)
    bias_chain = Stack{Vector{AbstractMarkovBias}}()

    state = MarkovAlgoState(
        Ref(grid), 0,
        0, nothing,
        MarkovOpContext((), MarkovBiasContext(allocator)),
        bias_chain, rng,
        Ref{Optional{Int}}()
    )
    logic_logln("Stated an algorithm run with seeds ", seeds,
                "\n Running first step of the algorithm...")
    logic_tab_in()
    markov_algo_step(algo, state)
    logic_tab_out()
    return state
end

function markov_algo_step(algo::MarkovAlgorithm, state::MarkovAlgoState, n_iterations::Integer = 1)
    # Internal behavior: if n_iterations is negative, run all the way to the end.
    if n_iterations < 0
        state.buffer_iter_count[] = nothing
    else
        state.buffer_iter_count[] = convert(Int, n_iterations)
    end

    logic_logln("Next tick of the algorithm will run ",
                  exists(state.buffer_iter_count[]) ?
                      state.buffer_iter_count[] :
                      "infinite",
                  " iterations")

    logic_tab_in()
    while !markov_algo_is_finished(algo, state) && (isnothing(state.buffer_iter_count[]) || state.buffer_iter_count[] > 0)
        # Note that when starting the algorithm, op_idx is 0 and op_state is initialized to nothing.
        if exists(state.op_state)
            n_iters_before::Int = get_something(state.buffer_iter_count[], 0)
            state.op_state = markov_op_iterate(algo.sequence[state.op_idx], state.op_state, state.grid,
                                               state.rng, state.op_context, state.buffer_iter_count)
            n_iters_after::Int = get_something(state.buffer_iter_count[], 0)

            n_op_ticks::Int = n_iters_before - n_iters_after
            @markovjunior_assert(n_op_ticks >= 0,
                                 "Your op ", typeof(algo.sequence[state.op_idx]),
                                   "*added* ", -n_op_ticks, " to n_ticks_left??")
            state.n_iterations += n_op_ticks

            logic_logln("Completed ", n_op_ticks, " iterations, totaling ", state.n_iterations)
        end

        if isnothing(state.op_state)
            logic_logln("Moving on to the next top-level op!")
            state.op_idx += 1
            logic_tab_in()
            if state.op_idx <= length(algo.sequence)
                state.op_context = MarkovOpContext(state)
                state.op_state = markov_op_initialize(algo.sequence[state.op_idx], state.grid,
                                                      state.rng, state.op_context)
                logic_logln(typeof(algo.sequence[state.op_idx]))
                logic_logln("Initial state:\n  ", state.op_state)
            else
                logic_logln("No more ops left. This algo run is about to end")
            end
            logic_tab_out()
        end
    end
    logic_tab_out()

    return nothing
end
markov_algo_finish(algo::MarkovAlgorithm, state::MarkovAlgoState) = markov_algo_step(algo, state, -1)

"
Releases all allocations of the MarkovAlgoState back to its allocator;
  remember to call this when you're done with the algorithm state!
"
function Base.close(s::MarkovAlgoState, owning_algo::MarkovAlgorithm)
    logic_logln("Closing algo for ", vsize(s.grid[]))
    logic_tab_in()
    if !markov_algo_is_finished(owning_algo, s)
        markov_op_cancel(owning_algo.sequence[s.op_idx], s.op_state, s.op_context)
    end
    markov_allocator_release_array(s.allocator, s.grid[])
    logic_tab_out()
    return nothing
end


markov_algo_is_finished(algo::MarkovAlgorithm, state::MarkovAlgoState) = (state.op_idx > length(algo.sequence))

"Operations can require a minimum number of dimensions for the grid"
markov_op_min_dimension(::AbstractMarkovOp)::Int = 1

markov_algo_grid(s::MarkovAlgoState) = s.grid[]
markov_algo_n_iterations(s::MarkovAlgoState) = s.n_iterations


###############
#  Implementation

"Collects the current `bias_chain` into a type-stable tuple for ops to process"
markov_algo_flatten_biases(mas::MarkovAlgoState) = Tuple(Iterators.flatten(mas.bias_chain))


#################
#  Thresholds

struct ThresholdByArea
    scale::Float32
end
struct ThresholdByLength
    scale::Float32
end
const ThresholdScalar = Union{Int, ThresholdByArea, ThresholdByLength}

"A randomly chosen threshold between two (unordered) values"
struct ThresholdRange
    a::ThresholdScalar
    b::ThresholdScalar
end

const Threshold = Union{ThresholdScalar, ThresholdRange}


struct ThresholdInputs
    area::Float32
    length::Float32
    rng::PRNG
end

get_threshold(th::Threshold, in::ThresholdInputs)::Int = max(1, round(Int, get_raw_threshold(th, in)))
get_threshold(th::Threshold, grid::CellGrid, rng::PRNG)::Int = get_threshold(th, ThresholdInputs(
    convert(Float32, prod(vsize(grid))),
    convert(Float32, sum(vsize(grid)) / ndims(grid)),
    rng
))

get_raw_threshold(th::Int, in::ThresholdInputs) = convert(Float32, th)
get_raw_threshold(th::ThresholdByArea, in::ThresholdInputs) = convert(Float32, in.area) * th.scale
get_raw_threshold(th::ThresholdByLength, in::ThresholdInputs) = convert(Float32, in.length) * th.scale
function get_raw_threshold(th::ThresholdRange, in::ThresholdInputs)
    a::Float32 = get_raw_threshold(th.a, in)
    b::Float32 = get_raw_threshold(th.b, in)
    return lerp(a, b, rand(rng, Float32))
end
