# This is a complement to the determinism unit-test,
#   which takes longer but catches other kinds of bugs.

# Checks that scenes generate bit-identical grids under different CPU targets,
#   by running them in separate Julia processes compiled for each target (`-C`).
# A target without FMA (like `sandybridge`) is the important comparison against `native`.
#
#   julia scripts/check_determinism.jl                       # Default targets and scenes
#   julia scripts/check_determinism.jl native sandybridge -- Fields Terrain
#
# The first run for a new target recompiles every dependency for it, which takes a while;
#   these caches are kept alongside the normal ones.
# Results only need to match within one Julia version (see the README).

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

# Only the worker processes need the engine; it's loaded at top level so its methods are callable afterwards.
const IS_WORKER = !isempty(ARGS) && (ARGS[1] == "--worker")
IS_WORKER && @eval using MarkovJunior

const DEFAULT_TARGETS = [ "native", "sandybridge" ]
const DEFAULT_SCENES = [
    # Randomized field biases (`random`, `randomness=`), ranges, `curve`-free grammar sizes...
    ("Fields", (48, 48)),
    ("Terrain", (48, 48)),
    ("MinimapGenerator", (48, 48)),
    ("Floorplan", (64, 36)),
    ("Tower2D", (60, 40)),
    ("Tower3D", (24, 24, 24))
]
const SEEDS = (1, 2, 3)

"Worker mode: runs the scenes and prints `scene seed hash` lines"
function run_worker(scenes)
    MJ = Main.MarkovJunior
    println("have_fma=", Core.Intrinsics.have_fma(Float64))
    for (name, size) in scenes
        algo = MJ.markov_algo_parse(read(joinpath(@__DIR__, "..", "scenes", "$name.jl"), String))
        for seed in SEEDS
            channel = MJ.markov_algo_run(algo, size,
                                      MJ.MarkovTickSettings(skip_tagged_events = MJ.SkipTaggedEvents.all_skippable);
                                      seeds = seed)
            grid_hash = MJ.markov_algo_complete(g -> hash(copy(g)), channel)
            println(name, " ", seed, " ", repr(grid_hash))
        end
    end
end

function main(args)
    if IS_WORKER
        scene_names = args[2:end]
        run_worker(filter(s -> s[1] in scene_names, DEFAULT_SCENES))
        return 0
    end

    split_i = findfirst(==("--"), args)
    targets = isnothing(split_i) ? (isempty(args) ? DEFAULT_TARGETS : args) : args[1:(split_i - 1)]
    scene_names = isnothing(split_i) ? first.(DEFAULT_SCENES) : args[(split_i + 1):end]
    isempty(targets) && (targets = DEFAULT_TARGETS)

    outputs = Dict{String, Vector{String}}()
    for target in targets
        println(stderr, "Running ", length(scene_names), " scenes with -C ", target, "...")
        cmd = `$(Base.julia_cmd()) --startup-file=no -C $target $(@__FILE__) --worker $scene_names`
        lines = readlines(pipeline(cmd; stdin=devnull))
        println(stderr, "  ", lines[1])
        outputs[target] = lines[2:end]
    end

    reference = outputs[targets[1]]
    all_match = true
    for target in targets[2:end]
        for (a, b) in zip(reference, outputs[target])
            if a != b
                all_match = false
                println("MISMATCH between ", targets[1], " and ", target, ": `", a, "` vs `", b, "`")
            end
        end
    end
    println(all_match ? "All $(length(reference)) runs match across $(join(targets, ", "))" : "Results differ!")
    return all_match ? 0 : 1
end

exit(main(ARGS))
