"
A reimplementation of [https://github.com/mxgmn/MarkovJunior](mxgmn/MarkovJunior).

To run the GUI editor/player, simply call `markovjunior_run_gui()`.

# Usage within Julia

An algorithm instance can be defined with the macro `@markovjunior [options] begin ... end`.
The specific syntax is already [documented elsewhere](docs/dsl.md).
To parse the algorithm from a string (or existing AST), call `markov_algo_parse(source)`.

The algorithm runs as a coroutine within its own Task.
Start the task with `commsChannel::AlgoCommsChannel = markov_algo_run(algo, initial_state_or_resolution, [tick_settings]; ...)`
The `MarkovTickSettings` struct allows you to control how the algorithm runs, and can even be changed during the run.
For example you can signal the algorithm to cancel itself with `ticking.cancel_algo = true`,
  or fast-forward to the end by setting `ticking.skip_most_tagged_events = true` and
  `ticking.min_runtime_tick_priority = 99999999`.

The returned `Channel` is used to communicate with the coroutine;
  it's recommended to use `markov_algo_complete()` or `markov_algo_next()` to do that for you,
  but you can also control it manually through the protocol outlined in `markov_algo_run()`.
You can cancel the algorithm by simply closing the channel,
  or for a faster/cleaner finish (avoiding `catch` blocks) update its `MarkovTickSettings` as mentioned above.

In between some ticks are special tagged events.
There are a few important built-in events, and you can dispatch custom ones at certain points in your algorithm.
Built-in tags are stored in the following global variables:

* `TAG_ALGO_STARTING` for when the algorithm is about to start running.
* `TAG_NEW_GRID` for when an Op has reallocated (resized) the output grid.
This is also raised for the initial grid, immediately before `TAG_ALGO_STARTING`.
* `TAG_ALGO_COMPLETED` for when the algorithm completes successfully and is about to deallocate the grid.
* `TAG_ALGO_CANCELED` for when you canceled the algorithm and it's about to deallocate the grid.

# GUI runner Tool

Normally an OpenGL GUI tool is also available to test/debug the algorithm,
  which you can run with `markovjunior_run_gui()`.

To eliminate this tool and its assets and the use of folders at `pwd()`,
  set the environment variable `JMJ_BUILD_WITH_TOOL` to 0.

# IPC runner Tool

To let game engines and other tools use this library, run this package as its own process
  executing `markovjunior_run_ipc(...)`.
Then connect to its *named pipe* and implement the communication protocol outlined in the readme.
"
module MarkovJunior

using Random, Setfield, Profile, Printf, Sockets
const System = Base.Sys

using MacroTools, HybridStructs, NamedTupleTools, ConcurrentUtilities
using OrderedCollections, StaticArrays, DataStructures, Statistics, Combinatorics
using StructTypes, JSON3, ArgParse
using GLFW, CImGui, CSyntax
const LibCImGui = CImGui.LibCImGui

using Bplus; @using_bplus

include("compatibility.jl")


@make_toggleable_asserts markovjunior_
@decentralized_module_init
const BUILT_WITH_TOOL = convert(Bool, get(ENV, "JMJ_BUILD_WITH_TOOL", 1))


# Core library:
include("utils.jl")
include("cells.jl")
include("allocator.jl")
include("algo.jl")
include("dsl.jl")
include("thresholds.jl")
include("algo_interfaces.jl")
include("op_rewrite.jl")
include("op_draw_box.jl")
include("op_sequence.jl")
include("bias_temperatue.jl")
include("bias_field.jl")
export @markovjunior, markov_algo_parse,
       MarkovAlgorithm, MarkovTickSettings,
       markov_algo_run, markov_algo_next, markov_algo_complete,
       markov_algo_to_string

# Tooling:
if BUILT_WITH_TOOL
    include("tool_files.jl")
    include("tool_render_3D.jl")
    include("tool_gui.jl")
    include("tool_main.jl")
    export markovjunior_run_gui
end

# Library:
# include("lib_interface.jl")
include("ipc_interface.jl")
export markovjunior_run_ipc

end # module