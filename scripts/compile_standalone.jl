#=
    Arguments:
        -exe to build as executable
        -dll to build as library
=#

println()

using Pkg, PackageCompiler

# Figure out what kind of thing we're building.
if ("-exe" in ARGS) + ("-dll" in ARGS) != 1
    error("Must give exactly one of '-exe' and '-dll'")
end
const MODE = if "-exe" in ARGS
    :exe
elseif "-dll" in ARGS
    :dll
else
    error("Unhandled: ", ARGS)
end
per_mode(if_exe, if_dll) = if MODE == :exe
    if_exe
elseif MODE == :dll
    if_dll
else
    error("Unhandled: ", MODE)
end

const PROJECT_DIR = abspath(joinpath(@__DIR__, ".."))

# Figure out where the product is going.
const OUTPUT_FOLDER_NAME = per_mode("build-exe", "build-dll")
const OUTPUT_DIR = abspath(joinpath(@__DIR__, "..", "locals", OUTPUT_FOLDER_NAME))

const COPY_TO_BUILD_PROJECT_RELATIVE = String[
    "assets",
    "docs",
    "scenes",
    "README.md",
    "LICENSE"
]

println("Compiling project at ", PROJECT_DIR,
        "\n    into ", OUTPUT_DIR,
        "\n    ...")
rm(OUTPUT_DIR, force=true, recursive=true)
if MODE == :exe
    create_app(PROJECT_DIR, OUTPUT_DIR,
        force=true,
        include_transitive_dependencies=false,
        executables=[ "JMarkovJunior" => "julia_main" ]
    )
elseif MODE == :dll
    create_library(PROJECT_DIR, OUTPUT_DIR,
        force=true,
        include_transitive_dependencies=false,
        lib_name="JMarkovJunior"
    )
end
println()
println()

# Copy our files to the build.
println("Copying project files to buld...")
const COPY_DEST = per_mode(joinpath(OUTPUT_DIR, "bin"), joinpath(OUTPUT_DIR, "gui_tool_files"))
for c_path in COPY_TO_BUILD_PROJECT_RELATIVE
    print("\t", c_path, " .")
    full_src = joinpath(PROJECT_DIR, c_path)
    full_dest = joinpath(COPY_DEST, c_path)
    cp(full_src, full_dest, force=true, follow_symlinks=true)
    println("..Done")
end
println()

# The preferred way to tell graphics drivers to prefer discrete GPU's
#    is to export certain symbols in the executable.
# Unfortunately PackageCompiler.jl doesn't have a way to do that yet,
#    so we need to use an external tool.
if MODE == :exe
    println("Using `nvpatch` to hint to graphics drivers to use discrete cards.")
    if isnothing(Sys.which("nvpatch"))
        @error "Unable to use `nvpatch` as it's not installed! Users will have to manually force discrete GPU in NVidia Control Panel (and AMD equivalent)"
    else
        println("If the file can't be patched, it's really bad luck and you need to jiggle the code until it compiles differently.")
        for f_name in [ "JMarkovJunior.exe", "Julia.exe" ]
            f_path = joinpath(OUTPUT_DIR, "bin", f_name)
            cmd = `nvpatch --enable "$f_path"`
            print("\t", cmd, " .")
            run(cmd)
            println("..Done")
        end
    end
    println()
end

dir_size_bytes::UInt = 0
for (path, dirs, files) in walkdir(OUTPUT_DIR)
    for file in files
        global dir_size_bytes += filesize(joinpath(path, file))
    end
end
println("Done! Total size is ", Base.format_bytes(dir_size_bytes))