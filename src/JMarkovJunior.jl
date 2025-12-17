"Based on [https://github.com/mxgmn/MarkovJunior](mxgmn/MarkovJunior)"
module JMarkovJunior

using Random, Setfield, Profile
using OrderedCollections, GLFW
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

"
A global debug flag that disables accelerated lookups,
falling back to simpler behavior that's much less likely to have bugs
"
const SKIP_CACHE = false # Useful for debugging the cache


include("cells.jl")
include("rules.jl")
include("inference.jl")
include("sequences.jl")
include("renderer.jl")


const DEFAULT_SEQUENCE = Sequence_Ordered(
    if false # Short-paths maze generator:
        [
            Sequence_DoN([
                CellRule("b", "w")
            ], 1, false),
            Sequence_DoAll([
                CellRule("bbw", "wEw")
            ], false),
            Sequence_DoAll([
                CellRule("E", "w")
            ], false)
        ]
    elseif true # Random walk maze generator
        [
            Sequence_DoN([
                CellRule("b", "R")
            ], 1, false),
            Sequence_DoAll([
                CellRule("Rbb", "GGR"),
                CellRule("GGR", "Rww"),
                CellRule("R", "w")
            ], true)
        ]
    elseif true # Custom
        [
            Sequence_DrawBox(
                Box2Df(
                    min=v2f(0.5, 0.5),
                    max=v2f(0.5, 0.5)
                ),
                'w'
            ),
            Sequence_DrawBox(
                Box2Df(
                    min=v2f(0, 1),
                    max=v2f(1, 1)
                ),
                'B'
            ),
            Sequence_DrawBox(
                Box2Df(
                    min=v2f(0, 0),
                    max=v2f(1, 0)
                ),
                'N'
            ),
            Sequence_DoAll(
                [
                    CellRule("wbb", "wGw")
                ],
                true,
                AllInference(
                    # [
                    #     InferPath('w', 'N', 'b', recompute_each_time=false)
                    # ],
                    # 3.0f0
                )
            ),
            Sequence_DoAll(
                [
                    CellRule("G", "w")
                ],
                false
            )
        ]
    else
        [ ]
    end
)

function main(; resolution::v2i = v2i(100, 100),
                sequence::AbstractSequence = DEFAULT_SEQUENCE,
                seed = @markovjunior_debug(0x1a2a3b4b5c6c7d8d, rand(UInt32)),
                profile::Bool = false
             )::Int
    println("Seed: ", seed)
    rng = PRNG(seed)

    grid = fill(CELL_CODE_BY_CHAR['b'], resolution.data)

    @game_loop begin
        # Bplus.GL.Context constructor parameters:
        INIT(v2i(800, 800), "Markov Junior Playground",
             vsync=VsyncModes.off)

        # Initialization:
        SETUP = begin
            profile && Profile.start_timer()

            grid_tex_ref = Ref{Texture}()
            grid_pixels_buffer = Ref{Matrix{v3f}}()

            sequence_state = execute_sequence(
                sequence, grid, rng,
                start_sequence(sequence, grid, AllInference(), rng)
            )
            LOOP.max_fps = nothing
        end

        # Core loop:
        LOOP = begin
            if GLFW.WindowShouldClose(LOOP.context.window)
                break # Cleanly ends the game loop
            end

            # Tick the algorithm.
            if exists(sequence_state)
                sequence_state = execute_sequence(sequence, grid, rng, sequence_state)
            end
            render_markov_2d(grid, v3f(0.4, 0, 0.4), grid_pixels_buffer, grid_tex_ref)

            # Render the algorithm state.
            clear_screen(vRGBAf(1, 0, 1, 1))
            clear_screen(@f32(1)) # Depth
            simple_blit(grid_tex_ref[], output_curve=@f32(2.2))
        end

        TEARDOWN = begin
            profile && Profile.stop_timer()
        end
    end

    return 0
end

# Precompile code using a basic sequence.
# The sequence logic is super type-unstable, so this is very important to performance.
function run_game(; sequence::AbstractSequence = DEFAULT_SEQUENCE,
                    resolution::Int = 50)
    grid = fill(CELL_CODE_BY_CHAR['b'], (resolution, resolution))
    rng = PRNG(0x23423411)

    seq_state = start_sequence(sequence, grid, AllInference(), rng)
    while exists(seq_state)
        seq_state = execute_sequence(sequence, grid, rng, seq_state)
    end

    return grid
end
@info "Running once for precompilation..."
run_game()

end