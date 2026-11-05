##   Ops  ##

"Whether an Op changed the grid, using `nothing` to mean 'cosmetic event which could never change the grid'"
const OpChangedGrid = Union{Bool, Nothing}

"Combines the results of two Ops, where `nothing` doesn't count either way"
combine_op_results(a::OpChangedGrid, b::OpChangedGrid)::OpChangedGrid = if isnothing(a)
    b
elseif isnothing(b)
    a
else
    a | b
end


"
Executes the given op within an algorithm, returning whether it changed the grid (see `OpChangedGrid`).

Inherited bias states are mutated in place as the grid changes.
If this Op (or a nested one) reallocates the grid, the states may be replaced with new objects.
So Ops that run other Ops should watch `algo_state.grid_generation`,
  and call `markov_bias_current_states()` when it changes.
"
function markov_algo_run(op::AbstractMarkovOp, algo, algo_state,
                         inherited_bias_tuple, inherited_bias_state_tuple
                        )::OpChangedGrid
    throw(MethodError(markov_algo_run, (
        op, algo, algo_state, inherited_bias_tuple, inherited_bias_state_tuple
    )))
end


##  Biases  ##

"
Initialize the bias and return some state object for it (by default, `nothing`).
If the state needs to change after initialization, it must be a mutable object
  which `markov_bias_update()` modifies in place.
The state's type must stay the same until the grid is reallocated (see `markov_bias_new_grid()`).
"
markov_bias_initialize(bias::AbstractMarkovBias, algo, algo_state) = nothing

"
Updates this bias's state in place, after the given area of the grid has changed.
Returns nothing.

Note that some Ops can't be laser-precise, and will include pixels that didn't actually change.
In some circumstances the box will cover the entire grid!

By default, does nothing.
"
markov_bias_update(bias::AbstractMarkovBias, bias_state,
                   algo, algo_state,
                   grid_subset::BoxI{N}, old_subset_values::CellGrid{N}
                  ) where {N} = nothing

"
Updates this bias after an Op reallocated the grid (e.g. resizing it),
  returning its new state -- which may be a different type than the old one.
`algo_state.grid` is already the new grid; `old_grid` is the previous which will be released right after this.

By default, cleans up the old state and then initializes a new one.
"
function markov_bias_new_grid(bias::AbstractMarkovBias, bias_state,
                              algo, algo_state, old_grid::CellGrid)
    markov_bias_cleanup(bias, bias_state, algo, algo_state)
    return markov_bias_initialize(bias, algo, algo_state)
end

"Cleans up this bias after it's done being used (by default: does nothing)"
markov_bias_cleanup(bias::AbstractMarkovBias, bias_state, algo, algo_state) = nothing

"
Calculates the desirability of changing the given pixels, at the current moment.
Returns `nothing` if the action *must* not happen.
"
function markov_bias_calculate(bias::AbstractMarkovBias, bias_state,
                               algo, algo_state,
                               at::Union{CellLine{N}, CellRegion{N}}
                              )::Optional{Float32} where {N}
    throw(MethodError(markov_bias_calculate, (
        bias, bias_state, algo, algo_state, at
    )))
end


##  Bias helpers for Ops  ##

"
A scoped helper that initializes some new biases on top of the existing ones,
  runs your code, then removes those biases.

Your lambda should have the signature `(full_biases::NTuple{NB}, full_bias_states::NTuple{NB}) -> Any`,
  taking type-stable tuples of the biases and their state objects.

Ensures bias cleanup happens correctly even in weird circumstances like a reallocated grid.
"
function markov_bias_with_states(to_do, new_biases::Tuple{Vararg{AbstractMarkovBias}},
                                 algo, algo_state,
                                 inherited_biases::Tuple{Vararg{AbstractMarkovBias}},
                                 inherited_states::Tuple)
    n_inherited = length(inherited_states)
    @markovjunior_assert(length(algo_state.bias_state_stack) == n_inherited,
                         "Bias stack has ", length(algo_state.bias_state_stack),
                           " states, but there are ", n_inherited, " inherited biases")

    new_states = map(b -> markov_bias_initialize(b, algo, algo_state), new_biases)
    append!(algo_state.bias_state_stack, new_states)
    try
        return to_do(tuple(inherited_biases..., new_biases...),
                     tuple(inherited_states..., new_states...))
    finally
        for (i, new_bias) in enumerate(new_biases)
            markov_bias_cleanup(new_bias, algo_state.bias_state_stack[n_inherited + i],
                                algo, algo_state)
        end
        resize!(algo_state.bias_state_stack, n_inherited)
    end
end

"Updates every given bias after the given area of the grid changed (see `markov_bias_update()`)"
@inline function markov_bias_update_all!(biases::Tuple, bias_states::Tuple, algo, algo_state,
                                         grid_subset, old_subset_values)
    # 'map' unrolls over tuples, unlike 'foreach'.
    map(biases, bias_states) do bias, state
        markov_bias_update(bias, state, algo, algo_state, grid_subset, old_subset_values)
    end
    return nothing
end

"Extracts a type-stable tuple of the current state of every active bias."
markov_bias_current_states(algo_state)::Tuple = Tuple(algo_state.bias_state_stack)
