###############################
#  AbstractMarkovOp interface

"
Initializes the given operation, and either returns some 'state' object
  or `nothing` if the op can't do anything.

The grid is wrapped in a `Ref` so that it can be reallocated as desired,
  but by default this will call through to a non-`Ref` version.
"
markov_op_initialize(op::AbstractMarkovOp, grid::Ref{<:CellGrid}, rng::PRNG, context::MarkovOpContext)::Optional = markov_op_initialize(op, grid[], rng, context)

"
Iterate on an operation N times and return either the new 'state' object, or `nothing` if the operation finished.

The `n_ticks_left` will be modified to remove however many ticks were performed;
  unless its value is `nothing` in which case the op will execute everything ASAP.

The grid is wrapped in a `Ref` so that it can be reallocated as desired.

By default this will call through to a version with an unwrapped grid,
  then to a version which perfoms exactly one tick at a time.
"
markov_op_iterate(op::AbstractMarkovOp, state,
                  grid::Ref{<:CellGrid}, rng::PRNG,
                  context::MarkovOpContext,
                  n_ticks_left::Ref{Optional{Int}}
                 )::Optional = markov_op_iterate(op, state, grid[], rng, context, n_ticks_left)

function markov_op_iterate(op::AbstractMarkovOp, state,
                           grid::CellGrid, rng::PRNG,
                           context::MarkovOpContext, n_ticks_left::Ref{Optional{Int}}
                          )::Optional
    logic_logln("Running ", exists(n_ticks_left[]) ? n_ticks_left[] : "infinite",
                "individual ticks on ", typeof(op), " this iteration...")
    logic_tab_in()
    if isnothing(n_ticks_left[])
        while exists(state)
            state = markov_op_iterate(op, state, grid, rng, context)
        end
        logic_tab_out()
        return nothing
    else
        while n_ticks_left[] > 0 && exists(state)
            state = markov_op_iterate(op, state, grid, rng, context)
            n_ticks_left[] -= 1
        end
        logic_tab_out()
        return state
    end
end

"Called if the algorithm stops running early. This gives you a chance to clean up allocations"
markov_op_cancel(op::AbstractMarkovOp, state, context::MarkovOpContext) = begin
    logic_logln("Default-Canceling ", typeof(op))
    nothing
end


###############################
#  AbstractMarkovBias interface

"Initialize the bias and return some state object for it (by default, `nothing`)"
markov_bias_initialize(bias::AbstractMarkovBias,
                       grid::CellGrid, rng::PRNG,
                       context::MarkovBiasContext) = nothing
"
Updates this bias after the given area of the grid has changed.
Returns the new state object for this bias.

Default: does nothing, returns the same state.
"
markov_bias_update(bias::AbstractMarkovBias, bias_state,
                   grid::CellGrid{N}, subset::BoxI{N},
                   rng::PRNG
                  ) where {N} = bias_state
"
Cleans up this bias after it's done being used.
Default: does nothing.
"
markov_bias_cleanup(bias::AbstractMarkovBias, bias_state) = nothing

"Adds type-stability to the bias's state object (by default, returns `Nothing`)"
markov_bias_state_type(::Type{<:AbstractMarkovBias})::Type = Nothing

"
Calculates the desirability of this action, at this moment, in this grid.
Returns `nothing` if the action must not happen.
"
function markov_bias_calculate(bias::AbstractMarkovBias, bias_state,
                               grid::CellGrid{N}, at::CellLine{N},
                               rng::PRNG
                              )::Optional{Float32} where {N}
    error("Unhandled: ", typeof(bias))
end
