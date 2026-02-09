"Based on [https://github.com/mxgmn/MarkovJunior](mxgmn/MarkovJunior)"
module JMarkovJunior

using Random, Setfield, Profile, Printf
const System = Base.Sys

using OrderedCollections, GLFW, CImGui
using CSyntax # Simplifies CImGui calls
using Bplus; @using_bplus

"
Retrieves the key from an ordered dictionary (or element from an ordered set),
    by its index.

The OrderedCollections package doesnt have an official way to do this yet;
    see https://github.com/JuliaCollections/OrderedCollections.jl/issues/93.
"
function ordered_collection_get_idx(oc::Union{OrderedDict, OrderedSet}, idx::Integer)
    # Make sure cached ops have all been applied.
    iterate(oc)

    inner = if oc isa OrderedDict
        oc
    elseif oc isa OrderedSet
        oc.dict
    else
        error("Unhandled: ", typeof(oc))
    end

    return inner.keys[idx]
end


Bplus.@make_toggleable_asserts markovjunior_
#DEBUG: Enable debug mode.
# markovjunior_asserts_enabled() = true

"
A edit-time flag that disables accelerated lookups,
  falling back to simpler behavior that's much less likely to have bugs.
"
const SKIP_CACHE = false
"
An edit-time flag that enables verbose logging of Path Inference potential fields computation.
"
const LOG_INFERENCE_PATH = false


"Converts various pieces of MarkovJunior algorithm definitions into their DSL representation"
function dsl_string end

# Basics:
include("cells.jl")
include("rules.jl")
include("inference.jl")
include("sequences.jl")
include("dsl.jl")

# High-level:
scenes_path() = joinpath(@__DIR__, "..", "scenes")
assets_path() = joinpath(@__DIR__, "..", "assets")
include("renderer.jl")
include("gui.jl")


const ASSET_BYTES_EDITOR_FONT::Vector{UInt8} = read(joinpath(
    assets_path(), "FiraCode-VariableFont_wght.ttf"
))


const DEFAULT_SEQUENCE_STR = read(joinpath(scenes_path(), "MazeRandomWalk.jl"), String)

function main(; sequence_str::String = DEFAULT_SEQUENCE_STR,
                seed = @markovjunior_debug(0x1a2a3b4b5c6c7d8d, rand(UInt32))
             )::Int
    @game_loop begin
        INIT(
            v2i(800, 800), "Markov Junior Playground",
            vsync=VsyncModes.on
        )

        SETUP = begin
            LOOP.max_fps = nothing

            gui_editor_font = gui_add_font_from_memory_ttf(
                ASSET_BYTES_EDITOR_FONT, [ 19 ]
            )[1]
            gui = GuiRunner(sequence_str, gui_editor_font, seed)
        end
        LOOP = begin
            if GLFW.WindowShouldClose(LOOP.context.window)
                break
            end

            gui_main(gui, LOOP.delta_seconds)

            clear_screen(v4f(0.2, 0.2, 0.2, 0))
            clear_screen(1.0f0)
        end
        TEARDOWN = begin
            close(gui)
        end
    end

    return 0
end

function run_game(; sequence::ParsedMarkovAlgorithm = eval(Meta.parse(DEFAULT_SEQUENCE_STR)),
                    resolution::NTuple{2, Int} = get_something(
                        markov_fixed_resolution(sequence),
                        (50, 50)
                    ))
    @bp_check markov_fixed_dimension(sequence) in (nothing, 2)

    grid = fill(sequence.initial_fill, resolution)
    rng = PRNG(0x23423411)

    seq_state = start_sequence(sequence.main_sequence, grid, AllInference(), rng)
    while exists(seq_state)
        seq_state = execute_sequence(sequence.main_sequence, grid, rng, seq_state)
    end

    return grid
end

# Precompile code using a basic sequence.
# The sequence logic is super type-unstable, so this is very important to performance.
@info "Running once for precompilation..."
run_game(resolution = (25, 25))


end # module