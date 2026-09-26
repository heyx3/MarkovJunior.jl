"A scene file that we can assume is there and well-formed"
const FALLBACK_SCENE_NAME = "MazeRandomWalk.jl"

# Disk paths are set on initialization.
"Scenes are julia files containing a single `@markovjunior` statement"
SCENES_PATH::String = ""
"Assets are various internal files"
ASSETS_PATH::String = ""
"Locals are temp/user files, excluded from version control"
LOCALS_PATH::String = ""
push!(RUN_ON_INIT, () -> begin
    root_dir = pwd()
    global SCENES_PATH = joinpath(root_dir, "scenes")
    global ASSETS_PATH = joinpath(root_dir, "assets")
    global LOCALS_PATH = joinpath(root_dir, "locals")
    mkpath(LOCALS_PATH)
end)

path_scene(name) = joinpath(SCENES_PATH, name)
path_asset(name) = joinpath(ASSETS_PATH, name)
path_local(name) = joinpath(LOCALS_PATH, name)

const FONT_FILE_NAME = "FiraCode-VariableFont_wght.ttf"
ASSET_BYTES_EDITOR_FONT_BUFFER::Vector{UInt8} = preallocated_vector(UInt8, 300*1024)
function get_editor_font_bytes()
    if isempty(ASSET_BYTES_EDITOR_FONT_BUFFER)
        file_path = path_asset(FONT_FILE_NAME)

        resize!(ASSET_BYTES_EDITOR_FONT_BUFFER, filesize(file_path))
        open(io -> read!(io, ASSET_BYTES_EDITOR_FONT_BUFFER),
             file_path, "r")
    end
    return ASSET_BYTES_EDITOR_FONT_BUFFER
end

"Contains the state of the UI, serialized to disk"
const MEMORY_FILE_NAME = "UserSession.json"


function show_file_to_user(path::AbstractString)
    path = abspath(path)

    cmd = if Sys.iswindows()
        `cmd /C start "" $path`
    elseif Sys.isapple()
        `open $path`
    else
        `xdg-open $path`
    end

    try
        run(cmd; wait=false)
    catch e
        @warn "Couldn't open $path in a default application" exception=e
    end

    return nothing
end