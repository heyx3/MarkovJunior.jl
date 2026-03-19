"
A reimplementation of [https://github.com/mxgmn/MarkovJunior](mxgmn/MarkovJunior).

To run the GUI editor/player, simply call `markovjunior_run_tool()`.


# Algorithm Definition

Use `@markovjunior` to generate an algorithm instance.
The following functions allow you to work with an instance:

* `markov_algo_start` to start running it
* `markov_algo_step` to update it
* `markov_algo_is_finished` to check if it's done
* `markov_algo_finish` to run all the way to the end (no timeout for infinite loops)
* `markov_algo_grid` to read the current state of the grid

## Custom Ops

To define a new `AbstractMarkovOp` (like `@rewrite`), implement the following
  un-exported interface (see their doc-strings for more info):

* `markov_op_initialize`
* `markov_op_iterate`
* `markov_op_cancel` if you have any resources/allocations to release;
you may want to call this when your op finishes as well
* `markov_op_min_dimension` if your op implies the grid must have a certain number of dimensions
* `dsl_string`; note that ops are represented with macro call syntax
* `parse_markovjunior_op` is the inverse of `dsl_string`

## Custom Biases

To define a new `AbstractMarkovBias`, implement the following
  un-exported interface (see their doc-strings for more info):

* `markov_bias_initialize`
* `markov_bias_update`
* `markov_bias_cleanup`
* `markov_bias_calculate`
* `markov_bias_state_type` is important to provide type-stability
* `dsl_string`; note that biases are represented with function call syntax
* `parse_markovjunior_bias` is the inverse of `dsl_string`
* `check_markovjunior_biases` if your new bias has rules
  that need to be validated after parsing (e.g. only one instance allowed)

# GUI runner Tool

Normally an OpenGL GUI tool is also available to test/debug the algorithm,
  which you can run with `markovjunior_run_tool()`.

To eliminate this tool and its assets and the use of folders at `pwd()`,
  set the environment variable `JMJ_BUILD_WITH_TOOL` to 0.
"
module MarkovJunior

using Random, Setfield, Profile, Printf
const System = Base.Sys

using MacroTools, NamedTupleTools
using OrderedCollections, StaticArrays, DataStructures, Combinatorics
using StructTypes, JSON3
using GLFW, CImGui, CSyntax

using Bplus; @using_bplus

include("compatibility.jl")


@make_toggleable_asserts markovjunior_
@decentralized_module_init
const BUILT_WITH_TOOL = convert(Bool, get(ENV, "JMJ_BUILD_WITH_TOOL", 1))


# Core library:
include("utils.jl")
include("cells.jl")
include("algo.jl")
include("dsl.jl")
include("interface.jl")
include("op_rewrite.jl")
include("op_draw_box.jl")
include("op_sequence.jl")
include("bias_temperatue.jl")
export AbstractMarkovAllocator, AbstractMarkovBias, AbstractMarkovOp,
       MarkovAlgorithm, MarkovAlgoState, MarkovOpContext,
       markov_op_min_dimension, markov_algo_grid, markov_algo_n_iterations,
       markov_algo_start, markov_algo_step, markov_algo_is_started, markov_algo_is_finished
#

# Tooling:
if BUILT_WITH_TOOL
    include("tool_files.jl")
    include("tool_render_3D.jl")
    include("tool_gui.jl")
    include("tool_main.jl")
    export markovjunior_run_tool
end

end # module