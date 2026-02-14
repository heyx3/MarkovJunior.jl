"A reimplementation of [https://github.com/mxgmn/MarkovJunior](mxgmn/MarkovJunior)"
module MarkovJunior

using Random, Setfield, Profile, Printf, DataStructures, MacroTools
const System = Base.Sys

using OrderedCollections, GLFW, CImGui
using CSyntax # Simplifies CImGui calls
using StructTypes, JSON3
using Bplus; @using_bplus

# B+ and DataStructures both define `update! :(
const update! = Bplus.update!

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

"A scene file that we can assume is there and compiles"
const FALLBACK_SCENE = "MazeRandomWalk.jl"


"Converts various pieces of MarkovJunior algorithm definitions into their DSL representation"
function dsl_string end

# Basics:
include("cells.jl")

include("new_runner.jl")

if false

include("rules.jl")
include("inference.jl")
include("sequences.jl")
include("dsl.jl")

"Scenes are julia files containing a single `@markovjunior` statement"
scenes_path() = joinpath(@__DIR__, "..", "scenes")
scenes_path(name) = joinpath(scenes_path(), name)
"Assets are various internal files"
assets_path() = joinpath(@__DIR__, "..", "assets")
assets_path(name) = joinpath(assets_path(), name)
"Locals are temp/user files, excluded from version control"
locals_path() = joinpath(@__DIR__, "..", "locals")
locals_path(name) = joinpath(locals_path(), name)

# High-level:
include("renderer.jl")
include("gui.jl")


const ASSET_BYTES_EDITOR_FONT::Vector{UInt8} = read(assets_path("FiraCode-VariableFont_wght.ttf"))
include_dependency(assets_path("FiraCode-VariableFont_wght.ttf"))

const MEMORY_FILE_NAME = "UserSession.json"


function main()
    @game_loop begin
        INIT(
            v2i(1210, 800), "Markov Junior Playground",
            vsync=VsyncModes.on
        )

        SETUP = begin
            LOOP.max_fps = nothing

            mkpath(locals_path())

            # Put the Dear ImGUI config inside the locals directory.
            # We must store the path string in a variable so the underlying char array can't be GC-ed.
            imgui_config_path = locals_path("imgui.ini")
            imgui_config_path_c = Base.unsafe_convert(Ptr{Int8}, imgui_config_path)
            unsafe_store!(CImGui.GetIO().IniFilename, imgui_config_path_c)

            # Load user session memory from file, if available.
            memory_path = locals_path(MEMORY_FILE_NAME)
            memory = if isfile(memory_path)
                try
                    open(io -> JSON3.read(io, GuiMemory), memory_path, "r")
                catch e
                    @error "User session file is corrupt, at $memory_path! It will be wiped clean"
                    m =  GuiMemory()
                    open(io -> JSON3.write(io, m), memory_path, "w")
                    m
                end
            else
                m =  GuiMemory()
                open(io -> JSON3.write(io, m), memory_path, "w")
                m
            end
            last_written_memory = copy(memory)
            function update_memory()
                if memory != last_written_memory
                    open(io -> JSON3.write(io, memory), memory_path, "w")
                    last_written_memory = copy(memory)
                end
            end

            gui_editor_font = gui_add_font_from_memory_ttf(
                ASSET_BYTES_EDITOR_FONT, [ 19 ]
            )[1]
            gui = GuiRunner(memory, gui_editor_font)
        end
        LOOP = begin
            if GLFW.WindowShouldClose(LOOP.context.window)
                break
            end

            gui_main(gui, LOOP.delta_seconds)
            if (LOOP.frame_idx % 60) == 0
                update_memory()
            end

            clear_screen(v4f(0.2, 0.2, 0.2, 0))
            clear_screen(1.0f0)
        end
        TEARDOWN = begin
            close(gui)
        end
    end
end
function julia_main()::Cint
    try
        main()
        return 0
    catch e
        println(stderr, "ERROR running tool: ")
        showerror(stderr, e)
        println(stderr, "\n")
        return 1
    end
end

"Instantly runs a full 'game' of the given algorithm"
function run_game(; sequence::ParsedMarkovAlgorithm = parse_markovjunior(
                        read(scenes_path(FALLBACK_SCENE), String)
                    ),
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

# @info "Running once for precompilation..."
#TODO: Pick a comprehensive scene
run_game(resolution = (25, 25))

end # if false, commenting out the old code


end # module