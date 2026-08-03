"
A reimplementation of [https://github.com/mxgmn/MarkovJunior](mxgmn/MarkovJunior).

To run the GUI editor/player, simply call `markovjunior_run_gui()`.

# Usage within Julia

An algorithm instance can be defined with the macro `@markovjunior [options] begin ... end`.
The specific syntax is already [documented elsewhere](docs/dsl.md).
To parse the algorithm from a string (or existing AST), call `markov_algo_parse(source)`.

The algorithm runs as a coroutine within its own Task (i.e. worker thread).
Start the task with `commsChannel::AlgoCommsChannel = markov_algo_run(algo, initial_state_or_resolution, [tick_settings]; ...)`
The `TickSettings` struct allows you to control how the algorithm runs, and can even be changed during the run.
For example you can signal the algorithm to cancel itself with `ticking.cancel_algo = true`,
  or fast-forward to the end by setting `ticking.skip_most_tagged_events = true` and
  `ticking.min_runtime_tick_priority = 99999999`.

To control the algorithm, it's recommended to call `markov_algo_complete()` with a lambda for tagged events,
  and optionally a lambda for ticks as well.
Your lambdas are *blocking calls*, so you could do time-slicing by calling `sleep(x)` within them.
There are special built-in tagged events, as well as custom ones you can set up within the algorithm.
Built-in tags are stored in the following global variables:

* `TAG_ALGO_STARTING` for when the algorithm is about to start running.
* `TAG_NEW_GRID` for when an Op has reallocated (resized) the output grid.
This is also raised for the initial grid, immediately before `TAG_ALGO_STARTING`.
* `TAG_ALGO_COMPLETED` for when the algorithm completes successfully and is about to deallocate the grid.
* `TAG_ALGO_CANCELED` for when you canceled the algorithm and it's about to deallocate the grid.

If you want manual control, take the channel returned by `markov_algo_run()` (of type `AlgoCommsChannel`)
  and follow this protocol loop:

1. Read one element from the channel.
   1. If it's an `Int`, then this is a normal algorithm tick where the int is the tick's 'priority'.
   2. If it's a `Symbol`, then it's a tagged event.
      1. If the tag is `MarkovJunior.TAG_NEW_GRID`, then read one more channel element
which is a reference to the newly-allocated grid, of type `CellGrid{N}`.
      2. If the tag is `MarkovJunior.TAG_ALGO_COMPLETED` or `MarkovJunior.TAG_ALGO_CANCELED`,
then no more messages will come in and the current grid will be deallocated as soon as you let the algorithm resume.
      3. Otherwise it is an internal bookkeeping event (like `MarkovJunior.TAG_ALGO_STARTING`)
or a custom event you set up within the algorithm source. Do as you please.
2. Let the algorithm resume by writing `zero(Int)` to the channel.

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
error("#TODO: New export statements")

# Tooling:
if BUILT_WITH_TOOL
    include("tool_files.jl")
    include("tool_render_3D.jl")
    include("tool_gui.jl")
    include("tool_main.jl")
    export markovjunior_run_gui
end

# Library:
include("lib_interface.jl")
include("ipc_interface.jl")
export markovjunior_run_ipc

end # module