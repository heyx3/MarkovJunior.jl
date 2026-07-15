const GUI_LEGEND_DATA = map(enumerate(CELL_TYPES)) do (i,t)
    return (
        t.color,
        Float32(i)/Float32(length(CELL_TYPES)),
        " - $(t.char)  $(i-1)"
    )
end

const gVec2 = Bplus.GUI.gVec2
const gVec4 = Bplus.GUI.gVec4

"A color that indicates a GUI widget is debug-only"
const GUI_DEBUG_COLOR = v3f(0.2, 0.2, 0.2)

# Set up a Dear ImGUI dropdown for materials by name.
const MATERIAL_NAMES = [ string(m) for m in Render3D.UboMaterialMode.instances() ]
const MATERIAL_NAMES_C_ARRAY_RUNTIME = [ Ptr{UInt8}(C_NULL) for _ in MATERIAL_NAMES ]
push!(RUN_ON_INIT, () -> begin
    for (i, name) in enumerate(MATERIAL_NAMES)
        MATERIAL_NAMES_C_ARRAY_RUNTIME[i] = pointer(name)
    end
end)


const CELL_MAT_PRAGMA = :GuiMaterial
"""
Interpret a pragma statement about cell materials, updating the given array of them.
If something goes wrong, an error/warning is logged to stderr.
"""
function interpret_pragma_gui_material!(materials::AbstractVector{<:Render3D.UboData_Material}, p_args)
    # Grab the first argument, the color being configured.
    if isempty(p_args)
        @error "No arguments given for @pragma $CELL_MAT_PRAGMA"
        return nothing
    end
    cell_code = first(p_args)
    # Parse the color.
    if !isa(cell_code, Symbol)
        @error "First argument to @pragma $CELL_MAT_PRAGMA must be a cell color (e.g. R, w, L)! Got $cell_code"
        return nothing
    elseif !haskey(CELL_CODE_BY_SYMBOL, cell_code)
        @error "Unexpected cell color in `@pragma $CELL_MAT_PRAGMA $cell_code ...`! Expected one of $(Tuple(t.name for t in CELL_TYPES))"
        return nothing
    end
    cell_i = CELL_CODE_BY_SYMBOL[cell_code] + 1
    cell_mat = materials[cell_i]

    # Grab the second argument, the material type.
    MaterialEnum = Render3D.UboMaterialMode
    if length(p_args) < 2
        @error "No material type given for @pragma $CELL_MAT_PRAGMA $cell_code! Expected one of $(Tuple(string(m) for m in MaterialEnum.instances()))"
        return nothing
    end
    material_mode_symbol = first(Iterators.drop(p_args, 1))
    material_mode = MaterialEnum.from(material_mode_symbol)
    if isnothing(material_mode)
        @error "Invalid material type given for @pragma $CELL_MAT_PRAGMA $cell_code $material_mode_symbol! Expected one of $(Tuple(string(m) for m in MaterialEnum.instances()))"
        return nothing
    end
    cell_mat.mode = material_mode

    # Set up parsing of subsequent arguments.
    msg_header() = "@pragma $CELL_MAT_PRAGMA $cell_code $material_mode_symbol"
    mode_args = Iterators.drop(p_args, 2)
    n_mode_args = count(_->true, mode_args)
    function try_parse_scalar(idx::Int, name, default)::Float32
        @bp_check(default isa Real)

        raw = if n_mode_args < idx
            @warn "$(msg_header()) is missing the '$name' parameter, which will default to $default"
            default
        else
            first(Iterators.drop(mode_args, idx-1))
        end

        return if !isa(raw, Real)
            @warn "$(msg_header()) parameter '$name' is not a number literal! It will default to $default"
            convert(Float32, default)
        else
            convert(Float32, raw)
        end
    end
    function try_parse_vector(idx::Int, name, default::Vec)::typeof(default)
        raw = if n_mode_args < idx
            @warn "$(msg_header()) is missing the '$name' parameter, which will default to $default"
            default
        else
            first(Iterators.drop(mode_args, idx-1))
        end

        return if raw isa Real
            # Broadcast literal values across all components.
            typeof(default)(i->raw)
        elseif Base.isexpr(raw, :tuple)
            if any(x -> !isa(x, Real), raw.args)
                @warn "$(msg_header()) has a garbled '$name' parameter, which will now default to $(default.data)"
                default
            else
                typeof(default)(raw.args...)
            end
        else
            @warn "$(msg_header()) has a garbled '$name' parameter, which will now default to $(default.data)"
            default
        end
    end
    if material_mode == MaterialEnum.empty
        if n_mode_args > 0
            @warn "@pragma $CELL_MAT_PRAGMA $cell_code $material_mode_symbol has unnecessary arguments which will be ignored"
        end
    elseif material_mode == MaterialEnum.glass
        cell_mat.opacity = try_parse_scalar(1, "opacity", 0.25)
        cell_mat.roughness = try_parse_scalar(2, "roughness", 0.25)
        cell_mat.albedo = try_parse_vector(3, "albedo", v3f(1, 0, 1))
    elseif material_mode == MaterialEnum.dielectric
        cell_mat.roughness = try_parse_scalar(1, "roughness", 0.75)
        cell_mat.albedo = try_parse_vector(2, "albedo", v3f(1, 0, 1))
    elseif material_mode == MaterialEnum.metal
        cell_mat.roughness = try_parse_scalar(1, "roughness", 0.1)
        cell_mat.albedo = try_parse_vector(2, "albedo", v3f(1, 0, 1))
    elseif material_mode == MaterialEnum.light_source
        cell_mat.albedo = try_parse_vector(1, "albedo", v3f(4, 1, 4))
    else
        error("Unhandled case: ", material_mode)
    end
    
    return nothing
end


"User state that should persist between program runs; serializable to/from JSON"
mutable struct GuiMemory
    #TODO: Window size/fullscreen

    next_dimensionality::Int
    next_resolution::Vector{Int32}
    legend_edits_materials::Bool

    current_scene_file_name::String
    current_scene_src::String
    current_seed_src::String
    new_scene_file_name::String

    is_playing::Bool
    prefers_running_to_end::Bool
    ticks_per_second::Float32
    ticks_per_jump::Int32
    ticks_for_profile::Int32
    max_seconds_for_run_to_end::Float32

    rendering_dim::Int

    render3D_sun_color::v3f
    render3D_sun_dir::v3f
    render3D_ambient_light::v3f
end
GuiMemory() = GuiMemory(
    2, [ 64, 64 ], false,
    FALLBACK_SCENE_NAME, read(path_scene(FALLBACK_SCENE_NAME), String),
    "0x1234567890abcdef",
    "MyScene",
    false, false, 150, 10, 1000, 10,
    2,
    v3f(1, 1, 1), vnorm(v3f(1, -1, -1)), v3f(0.05, 0.05, 0.05)
)
StructTypes.StructType(::Type{GuiMemory}) = StructTypes.Mutable()

# Equality and copying is important, to see if the memory needs updating on-disk.
Base.:(==)(a::GuiMemory, b::GuiMemory) = all(
    (getfield(a, f) == getfield(b, f))
      for f in fieldnames(GuiMemory)
)
Base.copy(m::GuiMemory) = GuiMemory((
    Iterators.map(getfield(m, f) for f in fieldnames(GuiMemory)) do f
        if f isa Vector
            copy(f)
        else
            f
        end
    end
)...)


"The state of the GUI for our MarkovJunior tool"
mutable struct GuiRunner
    memory::GuiMemory
    elapsed_seconds::Float32

    algorithm::MarkovAlgorithm
    algorithm_state::Optional{MarkovAlgoState}

    current_seed::UInt64
    current_seed_display::String
    next_seed::GuiText

    next_algorithm::GuiText
    next_algorithm_font::Ptr{CImGui.LibCImGui.ImFont}
    algorithm_error_msg::String

    available_scenes::Vector{String}
    next_scene_name::GuiText
    current_scene_idx::Int
    current_scene_has_changes::Bool

    render_3D_assets::Render3D.App
    rendering::Union{
        Tuple{Val{2}, Array{v3f, 2}, Texture}, # 2D rendering
        Tuple{Val{3}, Render3D.Scene, Render3D.FullViewport} # 3D rendering
    }
    render_settings_are_open::Bool
    reset_render_settings_window::Bool
    is_dragging_camera::Bool

    visualize_shadowmap_instead::Bool
    visualized_shadowmap::Texture
    visualized_shadowmap_target::Target

    time_till_next_tick::Float32
    textures_to_destroy::Vector{Texture} # Textures are deleted the frame *after* they're no longer used,
                                         #  in case Dear ImGUI still has instructions to render them
end

function GuiRunner(memory::GuiMemory,
                   algorithm_editor_font::Ptr{CImGui.LibCImGui.ImFont})
    initial_error_msg::String = ""

    available_scenes = Vector{String}()
    query_available_scenes!(available_scenes)
    scene_idx = get_something(findfirst(s -> s == memory.current_scene_file_name, available_scenes), 0)

    parsed_algo::MarkovAlgorithm = try
        markov_algo_parse(memory.current_scene_src)
    catch e
        initial_error_msg = string(
            "Unable to compile initial scene (fallback used instead):\n",
            sprint(showerror, e)
        )
        markov_algo_parse(read(path_scene(FALLBACK_SCENE_NAME), String))
    end
    
    app = Render3D.App()
    runner = GuiRunner(
        memory,
        0.0f0,

        parsed_algo, nothing,

        # Current seed will be parsed+stringified in a moment; we just need the textbox to be set up.
        zero(UInt64), "[UNINITIALIZED]",
        GuiText(string(memory.current_seed_src)),

        GuiText(memory.current_scene_src,
            is_multiline=true,
            imgui_flags=CImGui.LibCImGui.ImGuiInputTextFlags_AllowTabInput
        ),
        algorithm_editor_font,

        initial_error_msg,

        available_scenes,
        GuiText(memory.new_scene_file_name),
        scene_idx,
        # 'current_scene_has_changes' will be computed next.
        false,

        app,
        if memory.rendering_dim == 2
            resolution_2D = ntuple(2) do i
                if i > length(memory.next_resolution)
                    1
                else
                    memory.next_resolution[i]
                end
            end
            (
                Val(2),
                fill(zero(v3f), resolution_2D...),
                # Make a 1D texture so the lazy-init logic catches and creates it for real.
                Texture(SpecialFormats.rgb10_a2, 1)
            )
        elseif memory.rendering_dim == 3
            (
                Val(3),
                Render3D.Scene(
                    sun_dir=memory.render3D_sun_dir,
                    sun_color=memory.render3D_sun_color
                ),
                Render3D.FullViewport(v2i(1200, 1200), v3f(3, 3, 3), memory.render3D_ambient_light)
            )
        else
            error("Unhandled: ", memory.rendering_dim)
        end,
        false, false, false,

        false,
        begin
            t = Texture(
                SpecialFormats.rgb10_a2,
                v2u(1, 1),
                sampler = TexSampler{2}(
                    pixel_filter = PixelFilters.smooth
                ),
                n_mips = 1
            )
            (t, Target(TargetOutput(tex=t), nothing))
        end...,

        -1.0f0,
        Vector{Texture}()
    )

    reset_gui_runner_algo(runner, true, false, true, true)
    if runner.rendering[1] isa Val{2}
        update_gui_runner_texture_2D(runner)
    elseif runner.rendering[1] isa Val{3}
        update_gui_runner_render_3D(runner, true)
    else
        error("Unhandled: ", typeof(runner.rendering))
    end
    interpret_all_pragma_gui_materials(runner, parsed_algo)
    update_gui_runner_scenes!(runner)

    return runner
end
function Base.close(runner::GuiRunner)
    exists(runner.algorithm_state) && close(runner.algorithm_state, runner.algorithm)
    foreach(close, runner.textures_to_destroy)
    close(runner.render_3D_assets)
    close(runner.visualized_shadowmap_target)
    close(runner.visualized_shadowmap)

    if runner.rendering[1] isa Val{2}
        close(runner.rendering[3])
    elseif runner.rendering[1] isa Val{3}
        close(runner.rendering[2])
        close(runner.rendering[3])
    else
        error("Unhandled: ", typeof(runner.rendering))
    end
end

function query_available_scenes!(output::Vector{String}, empty_output_first::Bool=true)
    empty_output_first && empty!(output)
    append!(output, (
        name for name in readdir(SCENES_PATH)
          if endswith(name, ".jl")
    ))
    return nothing
end
function update_gui_runner_scenes!(runner::GuiRunner)
    query_available_scenes!(runner.available_scenes)
    next_scene_idx = findfirst(s -> s == runner.memory.current_scene_file_name, runner.available_scenes)

    if isnothing(next_scene_idx)
        push!(runner.available_scenes, runner.memory.current_scene_file_name)
        sort!(runner.available_scenes)
        next_scene_idx = findfirst(s -> s == runner.memory.current_scene_file_name, runner.available_scenes)

        open(io -> write(io, string(runner.next_algorithm)),
             path_scene(runner.memory.current_scene_file_name),
             "w")
        runner.current_scene_has_changes = false

        runner.algorithm_error_msg = string(
            "Scene '", runner.memory.current_scene_file_name, "' wasn't found, ",
              "so it was recreated on-disk from this editor"
        )
    else
        runner.current_scene_has_changes = (runner.memory.current_scene_src !=
            read(path_scene(runner.memory.current_scene_file_name), String)
        )
    end

    runner.current_scene_idx = next_scene_idx
    return nothing
end

function update_gui_runner_texture_2D(runner::GuiRunner)
    (_, array::Matrix{v3f}, tex::Texture) = runner.rendering

    # Generate the pixel buffer for upload to the GPU.
    #  * Figure out the final 2D slice resolution:
    resolution_2D = if ndims(runner.algorithm_state.grid[]) == 1
        Vec(length(runner.algorithm_state.grid[]), 1)
    else
        Vec(size(runner.algorithm_state.grid[])[1:2]...)
    end
    #  * Reallocate if necessary:
    if vsize(array) != resolution_2D
        array = fill(zero(v3f), resolution_2D...)
    end
    #  * Update the contents:
    convert_pixel(u::UInt8) = if u == CELL_CODE_INVALID
        v3f(1, 0, 1)
    else
        CELL_TYPES[u+1].color
    end
    N = ndims(runner.algorithm_state.grid[])
    if N < 3
        array .= convert_pixel.(runner.algorithm_state.grid[])
    else
        array .= convert_pixel.(runner.algorithm_state.grid[])[
            :, :,
            (1 for i in 3:N)...
        ]
    end
 
    # Update the GPU texture.
    if (tex.type != TexTypes.twoD) || (tex.size.xy != vsize(runner.algorithm_state.grid[]).xy)
        push!(runner.textures_to_destroy, tex)
        tex = Texture(
            SimpleFormat(
                FormatTypes.normalized_uint,
                SimpleFormatComponents.RGB,
                SimpleFormatBitDepths.B8
            ),
            array,
            sampler = TexSampler{2}(
                pixel_filter = PixelFilters.rough
            ),
            n_mips = 1
        )
    else
        set_tex_pixels(tex, array)
    end

    runner.rendering = (Val(2), array, tex)
    return nothing
end
function update_gui_runner_render_3D(runner::GuiRunner, rerender_view::Bool)
    (_, scene::Render3D.Scene, viewport::Render3D.FullViewport) = runner.rendering

    # Get a 3D view of the grid.
    grid_slice = runner.algorithm_state.grid[]
    if ndims(grid_slice) < 3
        grid_slice = reshape(
            grid_slice,
            size(grid_slice)...,
            (1 for i in (ndims(grid_slice)+1):3)...
        )
    elseif ndims(grid_slice) > 3
        grid_slice = @view grid_slice[
            :, :, :,
            (1 for i in 4:ndims(grid_slice))...
        ]
    end

    grid_is_new = (size(grid_slice) != tuple(scene.grid_tex_3D.size.xyz...))
    Render3D.update_scene_grid!(scene, grid_slice)
    grid_is_new && Render3D.on_new_grid!(viewport, Vec(size(grid_slice)...))
    if rerender_view
        Render3D.render(runner.render_3D_assets, scene, viewport)
    end

    runner.memory.rendering_dim = 3
    runner.rendering = (Val(3), scene, viewport)
end

function interpret_all_pragma_gui_materials(runner::GuiRunner, algo::MarkovAlgorithm)
    if runner.rendering[1] isa Val{3}
        materials = runner.rendering[2].cell_materials
        if haskey(runner.algorithm.pragmas_map, CELL_MAT_PRAGMA)
            for statement_args in runner.algorithm.pragmas_map[CELL_MAT_PRAGMA]
                interpret_pragma_gui_material!(materials, statement_args)
            end
        end
    end
    return nothing
end
function reset_gui_runner_algo(runner::GuiRunner,
                               parse_new_seed::Bool, parse_new_algorithm::Bool, update_resolution::Bool,
                               allow_run_to_end::Bool)
    old_grid_size::Tuple = if exists(runner.algorithm_state)
        size(runner.algorithm_state.grid[])
    else
        (1, 1)
    end
    new_resolution::Vector{Int32} = [ Int32.(old_grid_size)... ]
    new_dims = length(old_grid_size)

    # Re-parse the algorithm if requested.
    if parse_new_algorithm
        runner.algorithm_error_msg = ""

        @markovjunior_assert(runner.memory.current_scene_src == string(runner.next_algorithm))
        runner.algorithm = try
            a = markov_algo_parse(runner.memory.current_scene_src)
            interpret_all_pragma_gui_materials(runner, a) # We only want to do this if parsing succeeds
            a
        catch e
            runner.algorithm_error_msg = string(
                "Failed to parse: ", sprint(showerror, e),
                "\n\nFalling back to previous successful algorithm"
            )
            runner.algorithm
        end
    end
    # Update resolution from the GUI editor, if requested.
    if update_resolution
        if !exists(runner.algorithm.fixed_dimension)
            new_dims = runner.memory.next_dimensionality
        end
        new_dims = max(new_dims, runner.algorithm.min_dimension)

        new_resolution = copy(runner.memory.next_resolution)
    end

    # Finalize the dimensionality and resolution.
    if exists(runner.algorithm.fixed_dimension)
        new_dims = runner.algorithm.fixed_dimension
    end
    new_dims = max(new_dims, runner.algorithm.min_dimension)
    new_resolution = map(1:new_dims) do i
        if i > length(new_resolution)
            one(Int32)
        else
            new_resolution[i]
        end
    end
    runner.memory.next_dimensionality = new_dims
    runner.memory.next_resolution = new_resolution

    # Update the RNG seed if requested.
    if parse_new_seed
        @markovjunior_assert(runner.memory.current_seed_src == string(runner.next_seed))

        as_int = tryparse(UInt64, runner.memory.current_seed_src)
        runner.current_seed = if exists(as_int)
            as_int
        else
            hash(string(runner.next_seed))
        end

        runner.current_seed_display = "Seed: 0x$(string(runner.current_seed, base=16))"
    end

    # Finally, we can start the algorithm!
    runner.algorithm_state = markov_algo_start(runner.algorithm, Tuple(new_resolution), runner.current_seed)

    if allow_run_to_end && runner.memory.prefers_running_to_end
        finish_gui_runner_algo(runner, Ref(false))
    end

    # Configure rendering.
    if runner.rendering[1] isa Val{2}
        (_, array, tex) = runner.rendering
        if size(runner.algorithm_state.grid[]) != size(array)
            array = let N = ndims(runner.algorithm_state.grid[])
                if N == 1
                    fill(zero(v3f), length(runner.algorithm_state.grid[]), 1)
                elseif N == 2
                    fill(zero(v3f), size(runner.algorithm_state.grid[]))
                else
                    fill(zero(v3f), size(runner.algorithm_state.grid[])[1:2])
                end
            end
        end
        runner.rendering = (Val(2), array, tex)
        runner.memory.rendering_dim = 2
    elseif runner.rendering[1] isa Val{3}
        update_gui_runner_render_3D(runner, true)
    else
        error("Unhandled: ", typeof(runner.rendering))
    end

    return nothing
end
function step_gui_runner_algo(runner::GuiRunner, n_iterations::Int)
    if exists(runner.algorithm_state) && !markov_algo_is_finished(runner.algorithm, runner.algorithm_state)
        markov_algo_step(runner.algorithm, runner.algorithm_state, n_iterations)
    end

    return nothing
end
function finish_gui_runner_algo(runner::GuiRunner, grid_has_changed::Ref{Bool})
    start_t = time()
    while !gui_runner_is_finished(runner)
        step_gui_runner_algo(runner, 50) #TODO: gradually increase tick count and time it
        grid_has_changed[] = true

        if (time() - start_t) > runner.memory.max_seconds_for_run_to_end
            runner.algorithm_error_msg = string(
                "ENDLESS RUN DETECTED: took longer than ",
                 runner.memory.max_seconds_for_run_to_end,
                 " seconds to complete the grid!",
                "\n  If this was the first run, try again now that the JIT is warmed up.",
                "\n  Otherwise increase the timeout."
            )
            break
        end
    end

    return nothing
end

gui_runner_is_finished(runner::GuiRunner)::Bool = isnothing(runner.algorithm_state) ||
                                                  markov_algo_is_finished(runner.algorithm, runner.algorithm_state)

function gui_main(runner::GuiRunner, delta_seconds::Float32, cam_input::Cam3D_Input{Float32})
    print_wnd_sizes::Bool = false && @markovjunior_debug(rand(Float32) < 0.01, false)
    pane_flags = |(
        CImGui.LibCImGui.ImGuiWindowFlags_HorizontalScrollbar,
        CImGui.LibCImGui.ImGuiWindowFlags_NoDecoration & (~CImGui.LibCImGui.ImGuiWindowFlags_NoScrollbar),
        CImGui.LibCImGui.ImGuiWindowFlags_NoBringToFrontOnFocus
    )

    foreach(close, runner.textures_to_destroy)
    empty!(runner.textures_to_destroy)

    runner.elapsed_seconds += delta_seconds

    # Update the 3D renderer if applicable.
    if runner.rendering[1] isa Val{3}
        (_, scene_3D, viewport_3D) = runner.rendering

        if runner.is_dragging_camera
            (viewport_3D.cam, viewport_3D.cam_settings) = cam_update(
                viewport_3D.cam, viewport_3D.cam_settings,
                cam_input, delta_seconds
            )
        end
        Render3D.tick_scene!(scene_3D, delta_seconds, runner.render_3D_assets)
        Render3D.render(runner.render_3D_assets, scene_3D, viewport_3D)
    end

    function render_settings_window_fn()
        if runner.rendering[1] isa Val{2}
            CImGui.Text("No 2D settings")
        elseif runner.rendering[1] isa Val{3}
            Render3D.scene_settings_gui!(runner.rendering[2])
            Render3D.viewport_settings_gui!(runner.rendering[3], runner.rendering[2])

            # Remember the current settings.
            runner.memory.render3D_sun_dir = runner.rendering[2].sun_dir
            runner.memory.render3D_sun_color = runner.rendering[2].sun_color_hdr
            runner.memory.render3D_ambient_light = runner.rendering[3].ambient_light
        else
            error("Unhandled: ", typeof(runner.rendering))
        end
    end
    if runner.render_settings_are_open
        if runner.reset_render_settings_window
            runner.reset_render_settings_window = false
            CImGui.SetNextWindowPos((30, 30))
            CImGui.SetNextWindowFocus()
        end
        @c gui_window(render_settings_window_fn, "Render Settings", &runner.render_settings_are_open, |(
            CImGui.LibCImGui.ImGuiWindowFlags_AlwaysAutoResize,
            CImGui.LibCImGui.ImGuiWindowFlags_MenuBar,
            CImGui.LibCImGui.ImGuiWindowFlags_NoCollapse
        ))
    end

    gui_next_window_space(
        Box2Df(
            min=v2f(0, 0),
            max=v2f(0.45, 1)
        ),
        min_pixel_size = v2i(323, -1)
    )
    gui_window("Runner", C_NULL, pane_flags) do
        content_size = convert(v2f, CImGui.GetContentRegionAvail())
        print_wnd_sizes && println("Rnn: ", CImGui.GetWindowSize())

        grid_has_changed = Ref(false)

        # Render settings:
        if runner.rendering[1] isa Val{2}
            if CImGui.Button("Switch to 3D", v2f(100, 27))
                # Close the 2D renderer.
                (_, pixel_buf, tex_2D) = runner.rendering
                push!(runner.textures_to_destroy, tex_2D)

                # Allocate a 3D renderer.
                runner.rendering = (
                    Val(3),
                    Render3D.Scene(
                        sun_dir = runner.memory.render3D_sun_dir,
                        sun_color = runner.memory.render3D_sun_color
                    ),
                    Render3D.FullViewport(v2i(1200, 1200), v3f(3, 3, 3), runner.memory.render3D_ambient_light)
                )
                (_, scene::Render3D.Scene, viewport::Render3D.FullViewport) = runner.rendering

                update_gui_runner_render_3D(runner, true)
            end
        elseif runner.rendering[1] isa Val{3}
            @markovjunior_debug(gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, GUI_DEBUG_COLOR) do
                @c CImGui.Checkbox("Display sun shadowmap", &runner.visualize_shadowmap_instead)
            end)
            if CImGui.Button("Switch to 2D", v2f(100, 27))
                # Close the 3D renderer.
                (_, scene, viewport) = runner.rendering
                close(viewport)
                close(scene)

                # Allocate a 2D renderer.
                # Use a BS format for the texture so that the update function notices and lazy-initializes it.
                runner.rendering = (
                    Val(2),
                    fill(zero(v3f), 1, 1),
                    Texture(SpecialFormats.rgb10_a2, v3u(1, 1, 1))
                )

                # Initialize the 2D renderer.
                update_gui_runner_texture_2D(runner)
                runner.memory.rendering_dim = 2
            end

            # Offer a popup for render settings.
            CImGui.SameLine(0, 20)
            if CImGui.Button("Render Settings")
                if runner.render_settings_are_open
                    runner.reset_render_settings_window = true
                else
                    runner.render_settings_are_open = true
                end
            end
        else
            error("Unhandled: ", typeof(runner.rendering))
        end

        # Current state:
        (img_size, display_tex, flip_uv_y) =
            if runner.rendering[1] isa Val{2}
                (_, array, tex) = runner.rendering
                (convert(v2f, tex.size.xy), tex, false)
            elseif runner.rendering[1] isa Val{3}
                (_, scene_3D, viewport_3D) = runner.rendering
                if runner.visualize_shadowmap_instead
                    # Update the visualization texture if needed.
                    if runner.visualized_shadowmap.size.xy != scene_3D.sun_pov.view_depth.size.xy
                        close(runner.visualized_shadowmap_target)
                        close(runner.visualized_shadowmap)
                        runner.visualized_shadowmap = Texture(
                            runner.visualized_shadowmap.format,
                            scene_3D.sun_pov.view_depth.size.xy,
                            sampler = convert(TexSampler{2}, runner.visualized_shadowmap.sampler),
                            n_mips = 1
                        )
                        runner.visualized_shadowmap_target = Target(
                            TargetOutput(tex=runner.visualized_shadowmap),
                            nothing
                        )
                    end
                    # Regenerate the visualization.
                    target_activate(runner.visualized_shadowmap_target)
                    target_clear(runner.visualized_shadowmap_target, v4f(0, 0, 0, 1))
                    Render3D.visualize_depth_tex(runner.render_3D_assets,
                                                 scene_3D.sun_pov.view_depth,
                                                 Box2Df(min=zero(v2f), size=one(v2f)),
                                                 runner.elapsed_seconds)
                    target_activate(nothing)

                    dp_out = runner.visualized_shadowmap
                    (convert(v2f, dp_out.size.xy), dp_out, true)
                else
                    (convert(v2f, viewport_3D.view_color.size.xy), viewport_3D.view_color, true)
                end
            else
                error("Unhandled: ", typeof(runner.rendering))
            end
        img_size = if runner.rendering[1] isa Val{2}
            # Scale up but never down; don't lose pixels
            min_img_size::Float32 = content_size.x - 10
            img_size * max(1.0f0, (min_img_size / img_size)...)
        elseif runner.rendering[1] isa Val{3}
            # Scale up/down to fit the content window.
            sz = v2f(content_size.x - 10, content_size.y - 300)
            v2f(i -> max(256.0f0, min(sz...)))
            #TODO: Recreate the 3D viewport whenever the display size changes
        else
            error("Unhandled: ", typeof(runner.rendering))
        end
        img_cursor_pos = CImGui.GetCursorScreenPos();
        CImGui.Image(gui_tex_handle(display_tex),
                     convert(gVec2, img_size),
                     gVec2(0, flip_uv_y ? 1 : 0), gVec2(1, flip_uv_y ? 0 : 1),
                     gVec4(1, 1, 1, 1), gVec4(0, 0, 0, 0))
        # Draw an invisible button over the image, and detect when the user click-drags it.
        CImGui.SetCursorScreenPos(img_cursor_pos)
        CImGui.InvisibleButton("RendererMouseCapture", convert(gVec2, img_size))
        if CImGui.IsItemActive()
            if !runner.is_dragging_camera
                GLFW.SetInputMode(
                    get_context().window, GLFW.CURSOR,
                    GLFW.CURSOR_DISABLED
                )
                runner.is_dragging_camera = true
            end
        else
            if runner.is_dragging_camera
                GLFW.SetInputMode(
                    get_context().window, GLFW.CURSOR,
                    GLFW.CURSOR_NORMAL
                )
                runner.is_dragging_camera = false
            end
        end

        CImGui.Separator()

        # Execute Play logic.
        if runner.memory.is_playing && !gui_runner_is_finished(runner)
            runner.time_till_next_tick -= delta_seconds
            n_ticks = 0
            while runner.time_till_next_tick <= 0
                n_ticks += 1
                grid_has_changed[] = true
                runner.time_till_next_tick += 1.0f0 / runner.memory.ticks_per_second
            end
            step_gui_runner_algo(runner, n_ticks)
        end

        # Runner control panel, diagram below:
        #=   C1    | C2|   C3  |     C4       |  C5|   C6
            [Step]                [Reset]
            [Jump]   N  Ticks     [Profile]     500  ticks
            [Pause]  N  Ticks     [Run to End]  45   max seconds
        =#
        LCIG = CImGui.LibCImGui
        CImGui.BeginTable("##RunButtons", 6, |(
            LCIG.ImGuiTableFlags_NoSavedSettings,
            LCIG.ImGuiTableFlags_RowBg,
            LCIG.ImGuiTableFlags_BordersInnerH,
            LCIG.ImGuiTableFlags_SizingFixedFit
        ))
        BUTTON_SIZE_RUN_PLAIN = v2f(80, 30)
        BUTTON_SIZE_RUN_SPECIAL = v2f(94, 45)
        BUTTON_VPAD_RUN_PLAIN = max(0.0, (BUTTON_SIZE_RUN_SPECIAL.y - BUTTON_SIZE_RUN_PLAIN.y) / 2)
        UNITS_VPAD = 12
        BUTTON_COLOR_RUN_SPECIAL = v3f(0.3, 0.1, 0.0)
        BUTTON_COLUMN_FLAGS = |(
            LCIG.ImGuiTableColumnFlags_NoResize,
            LCIG.ImGuiTableColumnFlags_NoReorder,
            LCIG.ImGuiTableColumnFlags_NoHide,
            LCIG.ImGuiTableColumnFlags_NoSort
        )
        CImGui.TableSetupColumn("##PlainButtons", BUTTON_COLUMN_FLAGS, 84)
        CImGui.TableSetupColumn("##PlainButtonParams", BUTTON_COLUMN_FLAGS, 50)
        CImGui.TableSetupColumn("##PlainButtonUnits", BUTTON_COLUMN_FLAGS, 85)
        CImGui.TableSetupColumn("##SpecialButtons", BUTTON_COLUMN_FLAGS, 98)
        CImGui.TableSetupColumn("##SpecialButtonParams", BUTTON_COLUMN_FLAGS, 50)
        CImGui.TableSetupColumn("##SpecialButtonUnits", BUTTON_COLUMN_FLAGS, 60)
        # First row
        CImGui.TableNextColumn()
            CImGui.Dummy(0, BUTTON_VPAD_RUN_PLAIN)
            CImGui.SetNextItemWidth(-1)
            if CImGui.Button("Step", BUTTON_SIZE_RUN_PLAIN)
                step_gui_runner_algo(runner, 1)
                grid_has_changed[] = true
            end
        CImGui.TableNextColumn()
        CImGui.TableNextColumn()
        CImGui.TableNextColumn()
            CImGui.SetNextItemWidth(-1)
            gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, BUTTON_COLOR_RUN_SPECIAL) do
                if CImGui.Button("Reset", BUTTON_SIZE_RUN_SPECIAL)
                    reset_gui_runner_algo(runner, false, false, false, true)
                    grid_has_changed[] = true
                end
            end
        CImGui.TableNextColumn()
        CImGui.TableNextColumn()
        # Second row:
        CImGui.TableNextColumn()
            CImGui.Dummy(0, BUTTON_VPAD_RUN_PLAIN)
            CImGui.SetNextItemWidth(-1)
            if CImGui.Button("Jump", BUTTON_SIZE_RUN_PLAIN)
                step_gui_runner_algo(runner, convert(Int, runner.memory.ticks_per_jump))
                grid_has_changed[] = true
            end
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragInt("##TicksPerJump", &runner.memory.ticks_per_jump, 1.0, 1, 999999999, "%d")
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("Ticks")
        CImGui.TableNextColumn()
            CImGui.SetNextItemWidth(-1)
            if gui_with_style(() -> CImGui.Button("Profile", BUTTON_SIZE_RUN_SPECIAL),
                              CImGui.LibCImGui.ImGuiCol_Button, BUTTON_COLOR_RUN_SPECIAL)
            #begin
                Profile.start_timer()
                step_gui_runner_algo(runner, runner.memory.ticks_for_profile)
                Profile.stop_timer()

                prof_text_path = path_local("ProfileResult.txt")
                open(prof_text_path, "w") do io::IO
                    ctx = IOContext(io, :displaysize=>(5000, 999999))

                    println(ctx, "========================================\n==    Flat")
                    Profile.print(ctx, format=:flat, mincount=100, sortedby=:count)

                    println(ctx, "\n\n\n")

                    println(ctx, "========================================\n==    Tree")
                    Profile.print(ctx, format=:tree, noisefloor=2.0)
                end

                #TODO: modal display of profiler results?
                if Sys.iswindows()
                    run(`cmd /C $(abspath(prof_text_path))`)
                else
                    run(`$(abspath(prof_text_path))`)
                end
            end
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragInt("##TicksForProfile", &runner.memory.ticks_for_profile, 1.0, 1, 999999999, "%d")
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("ticks")
        # Third row:
        CImGui.TableNextColumn()
            CImGui.Dummy(0, BUTTON_VPAD_RUN_PLAIN)
            CImGui.SetNextItemWidth(-1)
            new_is_playing = Ref(runner.memory.is_playing)
            gui_with_style(LCIG.ImGuiCol_Button, new_is_playing[] ? v3f(0.3, 0.8, 0.5) : v3f(0.8, 0.2, 0.05)) do
                if CImGui.Button(new_is_playing[] ? "Pause" : "Play", BUTTON_SIZE_RUN_PLAIN)
                    new_is_playing[] = !new_is_playing[]
                end
            end
            if new_is_playing[] && !runner.memory.is_playing
                runner.time_till_next_tick = 1.0f0 / runner.memory.ticks_per_second
            end
            runner.memory.is_playing = new_is_playing[]
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragFloat("##TicksPerSecond", &runner.memory.ticks_per_second,
                                0.1, 0.00001, 999999999.0, "%.0f", 1.0)
            runner.time_till_next_tick = min(runner.time_till_next_tick,
                                             1.0f0 / runner.memory.ticks_per_second)
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("ticks/s")
        CImGui.TableNextColumn()
            CImGui.SetNextItemWidth(-1)
            new_runs_to_end = Ref(runner.memory.prefers_running_to_end)
            if gui_with_style(() -> CImGui.Button("Run to End", BUTTON_SIZE_RUN_SPECIAL),
                              LCIG.ImGuiCol_Button,
                              new_runs_to_end[] ? v3f(0.3, 0.8, 0.5) : BUTTON_COLOR_RUN_SPECIAL)
            #begin
                new_runs_to_end[] = !new_runs_to_end[]
                if new_runs_to_end[]
                    finish_gui_runner_algo(runner, grid_has_changed)
                end
            end
            runner.memory.prefers_running_to_end = new_runs_to_end[]
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragFloat(
                "##MaxSecondsRunningToEnd",
                &runner.memory.max_seconds_for_run_to_end,
                0.1,
                0.1, 999999.9,
                "%.1f", 1.0
            )
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("s timeout")
        # All done!
        CImGui.EndTable()
        #TODO: B+ helpers for tables

        CImGui.Separator()

        # Seed data:
        new_seed = tryparse(UInt64, string(runner.next_seed))
        new_seed_is_int = exists(new_seed)
        if !new_seed_is_int
            new_seed = hash(string(runner.next_seed))
        end
        gui_within_group() do
            CImGui.Dummy(1, 10)
            CImGui.Text(runner.current_seed_display)
            CImGui.Dummy(10, 0); CImGui.SameLine()
            gui_with_item_width(150) do
                gui_text!(runner.next_seed)
                runner.memory.current_seed_src = string(runner.next_seed)
            end
            CImGui.SameLine(0, 5)
            CImGui.Text(new_seed_is_int ? "(is number)" : "(is string)")
        end
        seed_has_changed::Bool = (new_seed != runner.current_seed)
        CImGui.SameLine(0, 40)
        gui_within_group() do
            if gui_with_style(() -> CImGui.Button("Reset with new seed", v2f(150, 32)) && seed_has_changed,
                            LCIG.ImGuiCol_Button, v3f(0.2, 0.1, 0.1),
                            unchanged = seed_has_changed)
            #begin
                reset_gui_runner_algo(runner, true, false, false, true)
                grid_has_changed[] = true
            end
            if CImGui.Button("Reset with rnd seed", v2f(150, 32))
                runner.memory.current_seed_src = string(rand(UInt64))
                update!(runner.next_seed, runner.memory.current_seed_src)

                reset_gui_runner_algo(runner, true, false, false, true)
                grid_has_changed[] = true
            end
        end

        CImGui.Separator()

        # Resolution data:
        gui_with_item_width(60) do
            if exists(runner.algorithm.fixed_dimension)
                CImGui.LabelText("Dimension", string(runner.memory.next_dimensionality))
            else
                d = convert(Int32, runner.memory.next_dimensionality)
                @c CImGui.InputInt("Dimension", &d)
                runner.memory.next_dimensionality = clamp(
                    d,
                    max(2, runner.algorithm.min_dimension),
                    max(3, runner.algorithm.min_dimension)
                )
                d = runner.memory.next_dimensionality
            end
            CImGui.SameLine(0, 30)
            # Update the resolution to match the dimensions.
            while length(runner.memory.next_resolution) < runner.memory.next_dimensionality
                push!(runner.memory.next_resolution, one(Int32))
            end
            while length(runner.memory.next_resolution) > runner.memory.next_dimensionality
                deleteat!(runner.memory.next_resolution, length(runner.memory.next_resolution))
            end
            if runner.memory.next_dimensionality == 2
                CImGui.InputInt2("Resolution", Ref(runner.memory.next_resolution, 1))
            elseif runner.memory.next_dimensionality == 3
                CImGui.InputInt3("Resolution", Ref(runner.memory.next_resolution, 1))
            else
                CImGui.Text("UNHANDLED dim $(runner.memory.next_dimensionality)")
            end
        end
        resolution_is_different::Bool =
            (runner.memory.next_dimensionality != ndims(runner.algorithm_state.grid[])) ||
            any(t -> t[1]!=t[2],
                zip(runner.memory.next_resolution, size(runner.algorithm_state.grid[])))
        CImGui.SameLine(0, 10)
        gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, v3f(0.2, 0.1, 0.1), unchanged=resolution_is_different) do
         gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonHovered, v3f(0.2, 0.1, 0.1), unchanged=resolution_is_different) do
          gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonActive, v3f(0.2, 0.1, 0.1), unchanged=resolution_is_different) do
            if CImGui.Button("Reset with new size", v2f(150, 25)) && resolution_is_different
                reset_gui_runner_algo(runner, false, false, true, true)
                grid_has_changed[] = true
            end
        end end end

        # Debug widgets:
        @markovjunior_debug begin
            CImGui.Separator()
            CImGui.Dummy(1, 50)
            CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("DEBUG")

            gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, GUI_DEBUG_COLOR) do
                render_tex = if runner.rendering[1] isa Val{2}
                    runner.rendering[3]
                elseif runner.rendering[1] isa Val{3}
                    runner.rendering[3].view_color
                else
                    error("Unhandled: ", typeof(runner.rendering))
                end

                if CImGui.Button("Recompile Shaders", (200, 30))
                    try
                        Render3D.recompile_shaders!(runner.render_3D_assets)
                    catch e
                        runner.algorithm_error_msg = sprint(showerror, e)
                        nothing
                    end
                end
                if CImGui.GetContentRegionAvailWidth() > 460
                    CImGui.SameLine(0, 40)
                end
                if CImGui.Button("Log GUI draw calls", (200, 30))
                    println(stderr, "LOGGING WITH RENDER TEX ID ",
                            gui_tex_handle(render_tex))
                    service_GUI().debug_log_render_commands = true
                end

                if CImGui.Button("Log center-pixel of render", (230, 30))
                    values = fill(zero(v4f), 1, 1)
                    pixel_pos = round(v2u, render_tex.size.xy / 2)
                    get_tex_color(render_tex, values, TexSubset(
                        Box(center=pixel_pos, size=one(v2u))
                    ))
                    println(stderr, "Render at pixel ", pixel_pos, ": ", values[1, 1])
                end
                if CImGui.GetContentRegionAvailWidth() > 500
                    CImGui.SameLine(0, 20)
                end
                if (runner.rendering[1] isa Val{3}) && CImGui.Button("Log some sun-shadowmap pixels", (230, 30))
                    println(stderr, "Depth at shadowmap pixels:")
                    values = fill(zero(Float32), 1, 1)
                    for p in (render_tex.size.xy * 0.1, render_tex.size.xy * 0.5, render_tex.size.xy * 0.9)
                        pixel_pos = round(v2u, p)
                        get_tex_depth(
                            runner.rendering[2].sun_pov.view_depth,
                            values,
                            TexSubset(Box(center=pixel_pos, size=one(v2u)))
                        )
                        println(stderr, "\t", pixel_pos, " = ", values[1, 1])
                    end
                end
            end

            if runner.rendering[1] isa Val{3}
                CImGui.Text("Sun Shadowmap matrix: World=>UV")
                gui_with_indentation() do
                    m = runner.rendering[2].sun_shadowmap_world_to_texel
                    CImGui.InputFloat4("##SunShadowmapRow1", Ref(v4f(m[1, :]...)))
                    CImGui.InputFloat4("##SunShadowmapRow2", Ref(v4f(m[2, :]...)))
                    CImGui.InputFloat4("##SunShadowmapRow3", Ref(v4f(m[3, :]...)))
                    CImGui.InputFloat4("##SunShadowmapRow4", Ref(v4f(m[4, :]...)))
                end
                CImGui.Text("Sun Shadowmap: Sample coordinates")
                gui_with_indentation() do
                    m = runner.rendering[2].sun_shadowmap_world_to_texel
                    for (j, p_world) in enumerate((v3f(0.9, 0.9, 1), v3f(0.95, 0.95, 1), v3f(1.05, 1.05, 1)))
                        CImGui.Text("$j: $p_world")
                        gui_with_nested_id(j) do
                            p_sun_texel = m * vappend(p_world, 1.0f0)
                            p_biased_sun_texel = m * vappend(p_world + runner.rendering[3].shadowmap_world_bias * runner.rendering[2].sun_dir, 1.0f0)
                            CImGui.InputFloat4("Sun Texel", Ref(p_sun_texel))
                            CImGui.InputFloat3("Sun Texel (dehom)", Ref(p_sun_texel.xyz / p_sun_texel.w))
                            CImGui.InputFloat4("Sun Texel with Bias", Ref(p_biased_sun_texel))
                            CImGui.InputFloat3("Sun Texel with Bias (dehom)", Ref(p_biased_sun_texel.xyz / p_biased_sun_texel.w))
                        end
                    end
                end
            end
        end

        # Update the state texture, if any above code changed the state.
        if grid_has_changed[]
            if runner.rendering[1] isa Val{2}
                update_gui_runner_texture_2D(runner)
            elseif runner.rendering[1] isa Val{3}
                update_gui_runner_render_3D(runner, true)
            else
                error("Unhandled: ", typeof(runner.rendering))
            end
        end
    end

    gui_next_window_space(
        Box2Df(
            min=v2f(0.45, 0),
            max=v2f(0.65, 0.5)
        ),
        max_pixel_size = v2i(-1, 316)
    )
    gui_window("Legend", C_NULL, pane_flags) do
        print_wnd_sizes && println("Legend wnd: ", CImGui.GetWindowSize())
        CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Legend")
        if runner.rendering[1] isa Val{3}
            CImGui.SameLine(85)
            @c CImGui.Checkbox("Edit Materials", &runner.memory.legend_edits_materials)
        end
        gui_within_group() do
            # If rendering 3D, this is also where we edit each cell type's Material.
            scene_materials = if (runner.rendering[1] isa Val{3}) && runner.memory.legend_edits_materials
                (runner.rendering[2]::Render3D.Scene).cell_materials
            else
                Iterators.repeated(nothing)
            end

            for (i, (color, greyscale, text), scene_material) in zip(1:N_CELL_TYPES, GUI_LEGEND_DATA, scene_materials)
                gui_with_nested_id(i) do
                    gui_draw_rect(
                        GuiDrawCursorRelative(Box2Df(
                            min=v2f(0, 0),
                            size=v2f(15, 15)
                        ), true),
                        GuiDrawColorFilled(color)
                    )
                    CImGui.SameLine()
                    CImGui.Text(text)

                    if exists(scene_material)
                        CImGui.Dummy(20, 1)
                        CImGui.SameLine()
                        gui_with_item_width(100) do
                            @c CImGui.Combo("##Type", &scene_material.mode,
                                            MATERIAL_NAMES_C_ARRAY_RUNTIME,
                                            length(MATERIAL_NAMES_C_ARRAY_RUNTIME),
                                            length(MATERIAL_NAMES_C_ARRAY_RUNTIME))
                        end
                        gui_with_item_width(50) do
                            if scene_material.mode in tuple(Render3D.UboMaterialMode.dielectric, Render3D.UboMaterialMode.metal, Render3D.UboMaterialMode.glass)
                                CImGui.SameLine(0, 10)
                                @c CImGui.SliderFloat("Rough", &scene_material.roughness, 0, 1,
                                                      "%.2f", LibCImGui.ImGuiSliderFlags_AlwaysClamp)
                            end
                        end
                        gui_with_item_width(50) do
                            if scene_material.mode in tuple(Render3D.UboMaterialMode.glass)
                                CImGui.SameLine(0, 10)
                                @c CImGui.SliderFloat("Alpha", &scene_material.opacity, 0, 1,
                                                      "%.2f", LibCImGui.ImGuiSliderFlags_AlwaysClamp)
                            end
                        end
                        gui_with_item_width(200) do
                            if scene_material.mode in tuple(Render3D.UboMaterialMode.dielectric,
                                                            Render3D.UboMaterialMode.metal,
                                                            Render3D.UboMaterialMode.glass,
                                                            Render3D.UboMaterialMode.light_source)
                            #begin
                                CImGui.Dummy(20, 1)
                                CImGui.SameLine(0, 10)
                                @c CImGui.ColorEdit3("##Albedo", &scene_material.albedo,
                                                     LibCImGui.ImGuiColorEditFlags_HDR |
                                                       #LibCImGui.ImGuiColorEditFlags_NoSidePreview |
                                                       LibCImGui.ImGuiColorEditFlags_Float |
                                                       LibCImGui.ImGuiColorEditFlags_InputRGB)
                            end
                        end
                        CImGui.Dummy(1, 5)
                    end
                end
            end
        end

        #TODO: Also display a rules legend
    end

    gui_next_window_space(
        Box2Df(
            min=v2f(0.45, 0.5),
            max=v2f(0.65, 1.0)
        )
    )
    gui_window("Scenes", C_NULL, pane_flags) do
        CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Scenes")
        CImGui.Dummy(40, 1);

        BUTTON_SIZE = v2f(100, 35)
        if CImGui.Button("Refresh", BUTTON_SIZE)
            update_gui_runner_scenes!(runner)
        end
        gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
         gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonHovered, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
          gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonActive, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
            if CImGui.Button("Reset changes", BUTTON_SIZE) && runner.current_scene_has_changes
                update!(runner.next_algorithm,
                        read(path_scene(runner.memory.current_scene_file_name), String))
                runner.memory.current_scene_src = string(runner.next_algorithm)
                runner.current_scene_has_changes = false
            end
            CImGui.SameLine(0, 20)
            if CImGui.Button("Save changes", BUTTON_SIZE) && runner.current_scene_has_changes
                open(io -> print(io, string(runner.next_algorithm)),
                     path_scene(runner.memory.current_scene_file_name),
                     "w")
                runner.current_scene_has_changes = false
            end
        end end end

        next_scene_idx_c = convert(Int32, runner.current_scene_idx - 1)
        n_available_lines = max(2, round(Int,
            (CImGui.GetContentRegionAvail().y - 100) / CImGui.GetTextLineHeightWithSpacing()
        ))
        # If there are any unsaved changes, draw the scene list disabled.
        gui_with_style(CImGui.LibCImGui.ImGuiCol_FrameBg, v3f(0.2, 0.1, 0.1), unchanged=!runner.current_scene_has_changes) do
          gui_with_style(CImGui.LibCImGui.ImGuiCol_FrameBgHovered, v3f(0.2, 0.1, 0.1), unchanged=!runner.current_scene_has_changes) do
            @c CImGui.ListBox("##Scenes", &next_scene_idx_c,
                              runner.available_scenes, length(runner.available_scenes),
                              n_available_lines)
        end end
        if (next_scene_idx_c+1 != runner.current_scene_idx) && !runner.current_scene_has_changes
            runner.current_scene_idx = next_scene_idx_c+1
            runner.memory.current_scene_file_name = runner.available_scenes[runner.current_scene_idx]
            update!(runner.next_algorithm, read(path_scene(runner.memory.current_scene_file_name), String))
        end

        was_new_scene_name_edited::Bool = gui_text!(runner.next_scene_name)
        CImGui.SameLine(0, 20)
        runner.memory.new_scene_file_name = string(runner.next_scene_name)
        new_scene_file_and_extension = runner.memory.new_scene_file_name
        if !endswith(new_scene_file_and_extension, ".jl")
            new_scene_file_and_extension *= ".jl"
        end
        can_use_new_scene_name = none(c -> c in runner.memory.new_scene_file_name,
                                      ('/', '\\', '*'))
        if !can_use_new_scene_name && isempty(runner.algorithm_error_msg)
            runner.algorithm_error_msg = "Invalid path chars in new scene name!"
        end
        if was_new_scene_name_edited && can_use_new_scene_name
            can_use_new_scene_name = !isfile(path_scene(new_scene_file_and_extension))
            if !can_use_new_scene_name && isempty(runner.algorithm_error_msg)
                runner.algorithm_error_msg = "File with that name already exists!"
            end
        end
        if can_use_new_scene_name && !runner.current_scene_has_changes && CImGui.Button("Create scene")
            # Check if the file already exists.
            # This is a bit redundant with the above, but it's the only easy way to avoid
            #    constantly checking for the file every X seconds.
            if isfile(path_scene(new_scene_file_and_extension))
                can_use_new_scene_name = false
                runner.algorithm_error_msg = "File with that name already exists!"
            else try
                new_scene_str = "@markovjunior begin\n\t@rewrite 1 b=>w\n\t@rewrite wb=>ww\nend"
                open(io -> write(io, new_scene_str),
                     path_scene(new_scene_file_and_extension),
                     "w")

                # Switch to editing the new scene.
                runner.memory.current_scene_file_name = new_scene_file_and_extension
                runner.memory.current_scene_src = new_scene_str
                update!(runner.next_algorithm, new_scene_str)

                update_gui_runner_scenes!(runner)
            catch e
                runner.algorithm_error_msg = "Unable to create file: $(sprint(showerror, e))"
            end end
        end
    end

    gui_next_window_space(Box2Df(
        min=v2f(0.65, 0),
        max=v2f(1.0, isempty(runner.algorithm_error_msg) ? 1.0 : 0.8)
    ))
    gui_window("Editor", C_NULL, pane_flags) do
        content_size = convert(v2f, CImGui.GetContentRegionAvail())
        print_wnd_sizes && println("Editor wnd: ", CImGui.GetWindowSize())

        runner.next_algorithm.multiline_requested_size = round.(Ref(Int),
            (content_size - v2f(20, 50)).data
        )
        gui_with_font(runner.next_algorithm_font) do
            runner.current_scene_has_changes |= gui_text!(runner.next_algorithm)
            runner.memory.current_scene_src = string(runner.next_algorithm)
        end

        if CImGui.Button("Restart##WithNewAlgorithm")
            reset_gui_runner_algo(runner, false, true, true, true)
            if runner.rendering[1] isa Val{2}
                update_gui_runner_texture_2D(runner)
            elseif runner.rendering[1] isa Val{3}
                update_gui_runner_render_3D(runner, true)
            else
                error("Unhandled: ", typeof(runner.rendering))
            end
        end
    end

    if !isempty(runner.algorithm_error_msg)
        gui_next_window_space(Box2Df(
            min=v2f(0.65, 0.8),
            max=v2f(1.0, 1.0)
        ))
        gui_with_style(CImGui.LibCImGui.ImGuiCol_WindowBg, v3f(0.6, 0.2, 0.26)) do
            gui_window("ErrorMsg", C_NULL, CImGui.LibCImGui.ImGuiWindowFlags_NoDecoration) do
                print_wnd_sizes && println("Error Msg wnd:", CImGui.GetWindowSize())
                gui_with_text_wrap(Inf) do
                    CImGui.Text(runner.algorithm_error_msg)
                end
            end
        end
    end

    print_wnd_sizes && println()
end