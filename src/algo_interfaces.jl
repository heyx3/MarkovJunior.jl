##   Ops  ##

"Executes the given op within an algorithm and returns whether it was able to make any changes to the grid"
function markov_algo_run(op::AbstractMarkovOp, algo, algo_state,
                         inherited_bias_tuple, inherited_bias_state_tuple
                        )::Bool
    throw(MethodError(markov_algo_run, typeof.((
        op, algo, algo_state, inherited_bias_tuple, inherited_bias_state_tuple, flag_had_no_matches
    ))))
end


##  Biases  ##

"
Initialize the bias and return some state object for it (by default, `nothing`).
The state object needs to remain type-stable through this particular algorithm run.
"
markov_bias_initialize(bias::AbstractMarkovBias, algo, algo_state) = nothing

"
Updates this bias after the given area of the grid has changed.
Returns the new state object for this bias.

By default, does nothing and returns the same state.
"
markov_bias_update(bias::AbstractMarkovBias, bias_state,
                   algo, algo_state,
                   grid_subset::BoxI{N}, old_subset_values::CellGrid{N}
                  ) where {N} = bias_state

"Cleans up this bias after it's done being used (by default: does nothing)"
markov_bias_cleanup(bias::AbstractMarkovBias, bias_state, algo, algo_state) = nothing

"
Calculates the desirability of a rewrite action, at this moment, in this grid.
Returns `nothing` if the action must not happen.
"
function markov_bias_calculate(bias::AbstractMarkovBias, bias_state,
                               algo, algo_state,
                               grid::CellGrid{N}, at::Union{CellLine{N}, CellRegion{N}},
                               rng::PRNG
                              )::Optional{Float32} where {N}
    throw(MethodError(markov_bias_calculate, typeof.((
        bias, bias_state, algo, algo_state, grid, at, rng
    ))))
end
