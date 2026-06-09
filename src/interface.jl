###############################
#  AbstractMarkovOp interface

"
Gets the kind of object representing this op's state during execution;
  it must be stable within one full run of the op but can otherwise vary based on the state of the algorithm.

You can also treat this as a 'pre-initialize' step, but be aware
  it may be called some time before the Op really starts (e.g. at the beginning of a `@sequence`).
"
markov_op_state_type(op::AbstractMarkovOp, ::Type{<:CellGrid{NDims}}, ::PRNG, ::MarkovOpContext) where {NDims} = markov_op_state_type(typeof(op), Val(NDims))
markov_op_state_type(T::DataType, NDims::Val) = error("Not implemented for ", T, " in ", val_type(NDims), "D")

"
Initializes the given operation, and either returns some 'state' object
  or `nothing` if the op can't do anything.
The state object must match the output of `markov_op_state_type()`,
  which you can take as the second parameter.

The grid is wrapped in a `Ref` so that it can be reallocated as desired,
  but by default this will call through to a non-`Ref` version.
"
markov_op_initialize(op::AbstractMarkovOp, expected_state_type::DataType,
                     grid::Ref{<:CellGrid}, rng::PRNG, context::MarkovOpContext
                    )::Optional{expected_state_type} = markov_op_initialize(op, expected_state_type, grid[], rng, context)
markov_op_initialize(op::AbstractMarkovOp, expected_state_type::DataType,
                     grid::CellGrid, rng::PRNG, context::MarkovOpContext
                    )::Optional{expected_state_type} = markov_op_initialize(op, grid, rng, context)

"
Iterate on an operation N times and return either the new 'state' object, or `nothing` if the operation finished.
The state object must be type-stable (the same type that was passed in).

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
                 )::Optional{typeof(state)} = markov_op_iterate(op, state, grid[], rng, context, n_ticks_left)

function markov_op_iterate(op::AbstractMarkovOp, state,
                           grid::CellGrid, rng::PRNG,
                           context::MarkovOpContext, n_ticks_left::Ref{Optional{Int}}
                          )::Optional
    @logic_logln("Running ", exists(n_ticks_left[]) ? n_ticks_left[] : "infinite",
                "individual ticks on ", typeof(op), " this iteration...")
    @logic_tab_in()
    if isnothing(n_ticks_left[])
        while exists(state)
            state = markov_op_iterate(op, state, grid, rng, context)
        end
        @logic_tab_out()
        return nothing
    else
        while n_ticks_left[] > 0 && exists(state)
            state = markov_op_iterate(op, state, grid, rng, context)
            n_ticks_left[] -= 1
        end
        @logic_tab_out()
        return state
    end
end

"Called if the algorithm stops running early. This gives you a chance to clean up allocations"
markov_op_cancel(op::AbstractMarkovOp, state, context::MarkovOpContext) = begin
    @logic_logln("Default-Canceling ", typeof(op))
    nothing
end


###############################
#  AbstractMarkovBias interface

"Required for type-stability; may be different per-instance; defaults to `Nothing`"
function markov_bias_state_type(a::AbstractMarkovBias, ::Type{<:CellGrid{NDims}},
                                ::PRNG, ::MarkovBiasContext
                                )::Type where {NDims}
    return markov_bias_state_type(typeof(a), Val(NDims))
end
markov_bias_state_type(T::Type{<:AbstractMarkovBias}, ::Val)::Type = markov_bias_state_type(T)
markov_bias_state_type(::Type{<:AbstractMarkovBias}) = Nothing

"
Initialize the bias and return some state object for it (by default, `nothing`).
You can optionally take the output of `markov_bias_state_type()` as the second parameter.
"
markov_bias_initialize(bias::AbstractMarkovBias, expected_state_type,
                       grid::CellGrid, rng::PRNG,
                       context::MarkovBiasContext
                      )::expected_state_type = markov_bias_initialize(bias, grid, rng, context)
markov_bias_initialize(bias::AbstractMarkovBias,
                       grid::CellGrid, rng::PRNG,
                       context::MarkovBiasContext
                      ) = nothing

"
Updates this bias after the given area of the grid has changed.
Returns the new state object for this bias.

Default: does nothing, returns the same state.
"
markov_bias_update(bias::AbstractMarkovBias, bias_state,
                   grid::CellGrid{N}, subset::BoxI{N}, old_subset_values::CellGrid{N},
                   rng::PRNG
                  ) where {N} = bias_state
"
Cleans up this bias after it's done being used.
Default: does nothing.
"
markov_bias_cleanup(bias::AbstractMarkovBias, bias_state) = nothing

"
Calculates the desirability of a rewrite action, at this moment, in this grid.
Returns `nothing` if the action must not happen.
"
function markov_bias_calculate(bias::AbstractMarkovBias, bias_state,
                               grid::CellGrid{N}, at::Union{CellLine{N}, CellRegion{N}},
                               rng::PRNG
                              )::Optional{Float32} where {N}
    error("Unhandled: ", typeof(bias))
end
