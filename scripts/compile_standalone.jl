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
        lib_name="JMarkovJunior",
        header_files = [
            "src/lib_interface.h"
        ]
    )
end
println()
println()

# Copy our files to the build.
println("Copying project files to buld...")
const COPY_DEST = per_mode(joinpath(OUTPUT_DIR, "bin"), joinpath(OUTPUT_DIR, "gui_tool_files"))
mkpath(COPY_DEST)
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

# PackageCompiler.jl doesn't create a lib file for us to link into the dll,
#    so we either need to dynamically load all the functions or generate the .lib ourselves.
# We chose the latter.
if MODE == :dll
    const DLL_PATH = joinpath(OUTPUT_DIR, "bin", "JMarkovJunior.dll")
    const LIB_PATH = joinpath(OUTPUT_DIR, "bin", "JMarkovJunior.lib")
    const LIB_FINAL_PATH = joinpath(OUTPUT_DIR, "lib", "JMarkovJunior.lib")
    @assert(isfile(DLL_PATH), "Dll missing! At $DLL_PATH")
    @assert(!isfile(LIB_PATH), "We didn't create the lib yet but it exists!? At $LIB_PATH")
    @assert(!isfile(LIB_FINAL_PATH), "We didn't create the lib yet but it exists!? At $LIB_FINAL_PATH")

    "
    `vswhere` is a program provided by Microsoft which always exist at a particular path,
    and tells us where to find the rest of the VS toolchain.
    "
    const VSWHERE_PATH = "C:\\Program Files (x86)\\Microsoft Visual Studio\\Installer\\vswhere.exe"
    const GEN_LIB_PATH = joinpath(@__DIR__, "GenLibFromDll.exe")
    @assert(isfile(GEN_LIB_PATH))
    const MSG_CANT_GENERATE_LIB = "Without a generated .lib file, users can't easily link the compiled dll with the provided header."
    const MSG_NO_VS_TOOLS = "Install Visual Studio Build Tools from https://visualstudio.microsoft.com/visual-cpp-build-tools/."
    if isfile(VSWHERE_PATH)
        vs_tools_path = readchomp(`$VSWHERE_PATH -latest -property installationPath`)
        if isempty(vs_tools_path) || !isdir(vs_tools_path)
            @error "`vswhere` failed to report any VS tools on this computer! $MSG_NO_VS_TOOLS $MSG_CANT_GENERATE_LIB"
        else
            vs_cmd_line_path = joinpath(vs_tools_path, "VC", "Auxiliary", "Build", "vcvarsall.bat")
            if !isfile(vs_cmd_line_path)
                @error "VS toolchain at $vs_tools_path doesn't have a vcvarsall.bat (expected at $vs_cmd_line_path)! You could try reinstalling? $MSG_NO_VS_TOOLS $MSG_CANT_GENERATE_LIB"
            else
                println("Using VS build tools to generate a .lib file...")

                run_cmd = `cmd /C "$vs_cmd_line_path" x64 "&&" "$GEN_LIB_PATH" "$DLL_PATH"`
                run_success = try
                    run(run_cmd)
                    @assert(isfile(LIB_PATH), "Command succeeded but the lib file is missing, at $LIB_PATH")
                    true
                catch e
                    @error "Command failed: $run_cmd\n\t$(sprint(showerror, e))\n\t$MSG_CANT_GENERATE_LIB"
                    false
                end

                if run_success
                    mkpath(LIB_FINAL_PATH)
                    mv(LIB_PATH, LIB_FINAL_PATH, force=true)
                    println("\n\n\t.lib file created! At ", LIB_FINAL_PATH)

                    # Remove intermediate files created by .lib generation.
                    rm(joinpath(OUTPUT_DIR, "bin", "JMarkovJunior.def"))
                    rm(joinpath(OUTPUT_DIR, "bin", "JMarkovJunior.exp"))

                    # PackageCompiler.jl adds an extra header for Julia's startup/shutdown functions,
                    #    however we have our own version in our header to handle .lib linkage.
                    rm(joinpath(OUTPUT_DIR, "include", "julia_init.h"))
                end
            end
        end
    elseif Sys.iswindows()
        @error "`vswhere` not found! Should have been in $VSWHERE_PATH. $MSG_CANT_GENERATE_LIB"
    else
        @error "We don't yet have a solution for generating .lib files outside of Windows; you won't be able to easily link the compiled dll to the provided header."
    end
end

dir_size_bytes::UInt = 0
for (path, dirs, files) in walkdir(OUTPUT_DIR)
    for file in files
        global dir_size_bytes += filesize(joinpath(path, file))
    end
end
println("Done! Total size is ", Base.format_bytes(dir_size_bytes))