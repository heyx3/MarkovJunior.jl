#=
    Arguments:
        -exe to build as executable
        -dll to build as library
=#

println()

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
insert!(LOAD_PATH, 1, @__DIR__)

using Printf, PackageCompiler
using MarkovJunior; const MJ = MarkovJunior

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
        executables=[
            "JMarkovJunior_GuiTool" => "markovjunior_run_gui_main",
            "JMarkovJunior_IPC" => "markovjunior_run_ipc_main"
        ]
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
println("Copying project files to build...")
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

# Generate headers containing project constants.
println("Generating general headers...")
const INCLUDE_DIR = joinpath(OUTPUT_DIR, "include")
const CONSTS_FILE_NAME = "jmj_consts"
const HEADER_GUARD_NAME_1 = "JMARKOVJUNIOR_GENERATED_CONSTS_H"
const BASIC_ALGO_ESC = escape_string("""@markovjunior begin
    @rewrite 1 b=>w
    @rewrite wbb=>wgw
    @rewrite g=>w
end""")
mkpath(INCLUDE_DIR)
open(joinpath(INCLUDE_DIR, "$CONSTS_FILE_NAME.h"), "w") do file
    print(file, """
    #ifndef $HEADER_GUARD_NAME_1
    #define $HEADER_GUARD_NAME_1

    #define JMJ_N_GRID_VALUES $(MJ.N_CELL_TYPES)

    //The color for each possible grid value, as RGB triplets.
    static const float JMJ_GRID_COLORS[][3] = {
    $(map(MJ.CELL_TYPES) do c
        separator = (c.code == MJ.N_CELL_TYPES - 1) ? "" : ","
        return "\t{ $(@sprintf("%.2f", c.color.x))f, $(@sprintf("%.2f", c.color.y))f, $(@sprintf("%.2f", c.color.z))f }$separator\n"
    end...)
    };

    //The char representing each possible grid value.
    #define JMJ_GRID_CHARS "$(map(c -> c.char, MJ.CELL_TYPES)...)"

    //The full name of each possible grid value.
    static const char* const JMJ_GRID_NAMES[] = {
    $(map(MJ.CELL_TYPES) do c
        separator = (c.code == MJ.N_CELL_TYPES - 1) ? "" : ","
        return "\t\"$(escape_string(c.name))\"$separator\n"
    end...)
    };

    //A simple maze-generator algorithm that works in any number of dimensions, useful to verify your IPC code.
    #define JMJ_BASIC_MAZE "$BASIC_ALGO_ESC"

    #endif
    """)
end
open(joinpath(INCLUDE_DIR, "$CONSTS_FILE_NAME.hpp"), "w") do file
    print(file, """
    #ifndef $HEADER_GUARD_NAME_1
    #define $HEADER_GUARD_NAME_1

    #include <cstdint>
    #include <array>
    #include <string_view>

    namespace jmj {
        constexpr uint8_t NGridValues = $(MJ.N_CELL_TYPES);

        inline constexpr std::array<std::array<float, 3>, NGridValues> GridColors = {{
        $(map(MJ.CELL_TYPES) do c
            separator = (c.code == MJ.N_CELL_TYPES - 1) ? "" : ","
            return "\t{ $(@sprintf("%.2f", c.color.x))f, $(@sprintf("%.2f", c.color.y))f, $(@sprintf("%.2f", c.color.z))f }$separator\n\t"
        end...)
        }};

        inline constexpr std::string_view GridChars = "$(map(c->c.char, MJ.CELL_TYPES)...)";

        inline constexpr std::array<std::string_view, NGridValues> GridNames = {
        $(map(MJ.CELL_TYPES) do c
            separator = (c.code == MJ.N_CELL_TYPES - 1) ? "" : ","
            return "\t\"$(escape_string(c.name))\"$separator\n"
        end...)
        };

        //A simple maze-generator algorithm that works in any number of dimensions, useful to verify your IPC code.
        inline constexpr std::string_view BasicMaze = "$BASIC_ALGO_ESC";
    }

    #endif
    """)
end
open(joinpath(INCLUDE_DIR, "$CONSTS_FILE_NAME.json"), "w") do file
    print(file, """{
        "n_grid_values": $(MJ.N_CELL_TYPES),
        "grid_colors": [
        $(map(MJ.CELL_TYPES) do c
            separator = (c.code == MJ.N_CELL_TYPES - 1) ? "" : ","
            return "\t[ $(c.color.x), $(c.color.y), $(c.color.z) ]$separator\n"
        end...)
        ],
        "grid_chars": "$(map(c -> c.char, MJ.CELL_TYPES))",
        "grid_names": [
        $(map(MJ.CELL_TYPES) do c
            separator = (c.code == MJ.N_CELL_TYPES - 1) ? "" : ","
            return "\t\"$(escape_string(c.name))\"$separator\n"
        end...)
        ],
        "basic_maze": "$BASIC_ALGO_ESC"
    }""")
end

if MODE == :exe

    # The preferred way to tell graphics drivers to prefer discrete GPU's
    #    is to export certain symbols in the executable.
    # Unfortunately PackageCompiler.jl doesn't have a way to do that yet,
    #    so we need to use an external tool.
    println("Using `nvpatch` to hint to graphics drivers to use discrete cards.")
    if isnothing(Sys.which("nvpatch"))
        @error "Unable to use `nvpatch` as it's not installed! Users will have to manually force discrete GPU in NVidia Control Panel (and AMD equivalent)"
    else
        println("If the file can't be patched, it's really bad luck and you need to jiggle the code until it compiles differently.")
        for f_name in [ "JMarkovJunior_GuiTool.exe", "Julia.exe" ]
            f_path = joinpath(OUTPUT_DIR, "bin", f_name)
            cmd = `nvpatch --enable "$f_path"`
            print("\t", cmd, " .")
            run(cmd)
            println("..Done")
        end
    end
    println()

    # Generate headers containing important constants for the IPC server.
    println("Generating IPC headers...")
    const HEADER_GUARD_NAME_2 = "JMARKOVJUNIOR_IPC_GENERATED_H"
    const NAMED_PIPE_ESCAPED = escape_string(MJ.IPC_PIPE_PATH)
    open(joinpath(INCLUDE_DIR, "jmj_ipc.h"), "w") do file
        print(file, """
        #ifndef $HEADER_GUARD_NAME_2
        #define $HEADER_GUARD_NAME_2

        #define JMJ_IPC_NAMED_PIPE "$NAMED_PIPE_ESCAPED"
        #define JMJ_IPC_DEFAULT_MAX_GRID_BYTES $(MJ.IPC_DEFAULT_MAX_GRID_BYTE_SIZE)
        #define JMJ_IPC_DEFAULT_MAX_CLIENT_NAME $(MJ.IPC_DEFAULT_MAX_CLIENT_NAME_BYTES)
        #define JMJ_IPC_STDOUT_START_CODE $(MJ.IPC_MAIN_START_CODE)
        #define JMJ_IPC_STDOUT_STOP_CODE $(MJ.IPC_MAIN_STOP_CODE)

        #endif
        """)
    end
    open(joinpath(INCLUDE_DIR, "jmj_ipc.hpp"), "w") do file
        print(file, """
        #ifndef $HEADER_GUARD_NAME_2
        #define $HEADER_GUARD_NAME_2

        #include <string_view>
        #include <array>

        namespace jmj { namespace ipc {
            //Note: the null-terminator *is* there, but the string_view does not include it
            inline constexpr std::string_view NamedPipe = "$NAMED_PIPE_ESCAPED";

            constexpr size_t DefaultMaxGridBytes = $(MJ.IPC_DEFAULT_MAX_GRID_BYTE_SIZE);
            constexpr size_t DefaultMaxClientName = $(MJ.IPC_DEFAULT_MAX_CLIENT_NAME_BYTES);

            constexpr uint32_t StdoutStartCode = $(MJ.IPC_MAIN_START_CODE);
            constexpr uint32_t StdoutStopCode = $(MJ.IPC_MAIN_STOP_CODE);
        } }

        #endif
        """)
    end
    open(joinpath(INCLUDE_DIR, "jmj_ipc.json"), "w") do file
        print(file, """{
            "named_pipe": "$NAMED_PIPE_ESCAPED",
            "default_max_grid_bytes": $(MJ.IPC_DEFAULT_MAX_GRID_BYTE_SIZE),
            "default_max_client_name": $(MJ.IPC_DEFAULT_MAX_CLIENT_NAME_BYTES),
            "stdout_start_code": $(MJ.IPC_MAIN_START_CODE),
            "stdout_stop_code": $(MJ.IPC_MAIN_STOP_CODE)
        }""")
    end
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