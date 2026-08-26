##   Ops  ##

"
Executes the given op within an algorithm.
Returns two things:

1. Whether it was able to make any changes to the grid
2. The new states of the inherited biases
"
function markov_algo_run(op::AbstractMarkovOp, algo, algo_state,
                         inherited_bias_tuple, inherited_bias_state_tuple
                        )::Tuple{Bool, typeof(inherited_bias_state_tuple)}
    throw(MethodError(markov_algo_run, (
        op, algo, algo_state, inherited_bias_tuple, inherited_bias_state_tuple
    )))
end


##  Biases  ##

"
Initialize the bias and return some state object for it (by default, `nothing`).
The state object needs to remain type-stable through this particular algorithm run.
"
markov_bias_initialize(bias::AbstractMarkovBias, algo, algo_state) = nothing

"
Updates this bias after the given area of the grid has changed.
Returns the new state object for this bias, exactly the same type as the input state.

Note that some Ops can't be laser-precise, and will include pixels that didn't actually change.
In some circumstances the box will cover the entire grid!

By default, does nothing and returns the same state.
"
markov_bias_update(bias::AbstractMarkovBias, bias_state,
                   algo, algo_state,
                   grid_subset::BoxI{N}, old_subset_values::CellGrid{N}
                  ) where {N} = bias_state

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
