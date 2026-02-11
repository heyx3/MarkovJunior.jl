const GUI_LEGEND_DATA = map(enumerate(CELL_TYPES)) do (i,t)
    return (
        t.color,
        Float32(i)/Float32(length(CELL_TYPES)),
        " - $(t.char)"
    )
end

const gVec2 = Bplus.GUI.gVec2
const gVec4 = Bplus.GUI.gVec4

const RENDER_2D_MODE_NAMES = collect(string.(Render2DMode.instances()))


"The state of the GUI for our MarkovJunior tool"
mutable struct GuiRunner
    draw_monochrome::Bool

    next_dimensionality::Int
    next_resolution::Vector{Int32}

    state_grid::CellGrid
    state_grid_tex2D_slice::CellGrid{2}
    state_texture::Texture

    algorithm::ParsedMarkovAlgorithm
    algorithm_state::Optional
    algorithm_rng::PRNG

    render2D::Render2DData

    current_seed::UInt64
    current_seed_display::String
    next_seed::GuiText

    next_algorithm::GuiText
    next_algorithm_font::Ptr{CImGui.LibCImGui.ImFont}
    algorithm_error_msg::String

    available_scenes::Vector{String}
    current_scene::String
    current_scene_idx::Int
    current_scene_has_changes::Bool

    is_playing::Bool
    ticks_per_second::Float32
    time_till_next_tick::Float32

    ticks_per_jump::Int32
    ticks_for_profile::Int32
    max_seconds_for_run_to_end::Float32
end

function GuiRunner(initial_algorithm_str::String,
                   algorithm_editor_font::Ptr{CImGui.LibCImGui.ImFont},
                   seed="0x1234567890abcdef")
    # Create a half-baked initial instance, then "restart" the algorithm.
    fake_size = 64
    runner = GuiRunner(
        false,

        2, [ fake_size, fake_size ],
        fill(zero(UInt8), fake_size, fake_size),
        fill(zero(UInt8), fake_size, fake_size),
        Texture(SpecialFormats.rgb10_a2, v2u(1, 2)), # Size must not match fake_size
                                                     #   or else it won't be properly reallocated

        @markovjunior(begin end), nothing, PRNG(1),
        Render2DData(),

        zero(UInt64), "[NULL]",
        GuiText(string(seed)),

        GuiText(initial_algorithm_str,
            is_multiline=true,
            imgui_flags=CImGui.LibCImGui.ImGuiInputTextFlags_AllowTabInput
        ),
        algorithm_editor_font,
        "Uh-OH: UNINITIALIZED!!!!",

        String[ ], "", 1, false,

        false, 150.0f0, -1.0f0,

        10, 1000, 10.0f0
    )

    reset_gui_runner_algo(runner, true, true, true)
    update_gui_runner_texture_2D(runner)
    update_gui_runner_scenes(runner)

    return runner
end
function Base.close(runner::GuiRunner)
    close(runner.state_texture)
    isassigned(runner.render2D.output) && close(runner.render2D.output[])
end

function update_gui_runner_scenes(runner::GuiRunner)
    empty!(runner.available_scenes)
    append!(runner.available_scenes,
        (name for name in readdir(scenes_path())
          if endswith(name, ".jl"))
    )
    next_scene_idx = findfirst(s -> s == runner.current_scene, runner.available_scenes)

    if isnothing(next_scene_idx)
        # On initial construction, the "current scene" is left empty.
        if isempty(runner.current_scene)
            next_scene_idx = 1
            runner.current_scene = runner.available_scenes[1]
            runner.current_scene_has_changes = false
        else
            push!(runner.available_scenes, runner.current_scene)
            sort!(runner.available_scenes)
            next_scene_idx = findfirst(s -> s == runner.current_scene, runner.available_scenes)

            open(io -> write(io, string(runner.next_algorithm)),
                joinpath(scenes_path(), runner.current_scene),
                "w")
            runner.current_scene_has_changes = false
        end
    end

    runner.current_scene_idx = next_scene_idx
    return nothing
end

function update_gui_runner_texture_2D(runner::GuiRunner)
    if (runner.state_texture.type != TexTypes.twoD) || (runner.state_texture.size.xy != vsize(runner.state_grid).xy)
        runner.state_texture = Texture(
            SimpleFormat(
                FormatTypes.uint,
                SimpleFormatComponents.R,
                SimpleFormatBitDepths.B8 # 4-bit isn't valid because there's only one component
            ),
            runner.state_grid,
            sampler = TexSampler{ndims(runner.state_grid)}(
                pixel_filter = PixelFilters.rough
            ),
            n_mips = 1
        )
        runner.state_grid_tex2D_slice = fill(zero(UInt8), runner.state_texture.size.xy...)
    else
        runner.state_grid_tex2D_slice .= @view runner.state_grid[
            ntuple(i -> (i>2) ? 1 : Colon(),
                   ndims(runner.state_grid))...
        ]
        set_tex_pixels(runner.state_texture, runner.state_grid_tex2D_slice)
    end

    render_markov_2d(runner.state_grid_tex2D_slice, v3f(1, 0, 1),
                     runner.render2D,
                     runner.algorithm.main_sequence, runner.algorithm_state)

    return nothing
end

function reset_gui_runner_algo(runner::GuiRunner,
                               parse_new_seed::Bool, parse_new_algorithm::Bool, update_resolution::Bool)
    # Re-parse the algorithm textbox, if requested.
    if parse_new_algorithm
        runner.algorithm_error_msg = ""
        try
            algorithm_ast = Meta.parse(string(runner.next_algorithm), filename=runner.current_scene[1:end-3])
            if !Base.isexpr(algorithm_ast, :macrocall) || algorithm_ast.args[1] != Symbol("@markovjunior")
                runner.algorithm_error_msg = string(
                    "Invalid header: Expected `@markovjunior ... begin ... end`",
                    "\n\nFalling back to previous successfully-parsed algorithm"
                )
            else
                runner.algorithm = eval(algorithm_ast)
            end
        catch e
            runner.algorithm_error_msg = string(
                "Failed to parse: ", sprint(io -> showerror(io, e)),
                "\n\nFalling back to previous successfully-parsed algorithm"
            )
        end

        # If the algorithm has a fixed dimensionality, trim 'next_resolution' to fit.
        n_dims = markov_fixed_dimension(runner.algorithm)
        if exists(n_dims)
            while length(runner.next_resolution) < n_dims
                push!(runner.next_resolution, one(Int32))
            end
            while length(runner.next_resolution) > n_dims
                deleteat!(runner.next_resolution, convert(Int32, length(runner.next_resolution)))
            end
        end
    end

    # Initialize the grid.
    if update_resolution
        dimensions = let d = markov_fixed_dimension(runner.algorithm)
            if exists(d)
                d
            else
                runner.next_dimensionality
            end
        end
        resolution = let r = markov_fixed_resolution(runner.algorithm)
            if exists(r)
                r
            else
                Tuple(runner.next_resolution)
            end
        end
        runner.state_grid = fill(runner.algorithm.initial_fill, resolution)
    else
        fill!(runner.state_grid, runner.algorithm.initial_fill)
    end

    # Initialize the RNG.
    if parse_new_seed
        as_int = tryparse(UInt64, string(runner.next_seed))
        runner.current_seed = if exists(as_int)
            as_int
        else
            hash(string(runner.next_seed))
        end

        runner.current_seed_display = "Seed: 0x$(string(runner.current_seed, base=16))"
    end
    runner.algorithm_rng = PRNG(PrngStrength.strong, runner.current_seed)

    # Initialize the algorithm sequence.
    runner.algorithm_state = start_sequence(
        runner.algorithm.main_sequence,
        runner.state_grid, AllInference(),
        runner.algorithm_rng
    )

    return nothing
end
function step_gui_runner_algo(runner::GuiRunner)
    if exists(runner.algorithm_state)
        runner.algorithm_state = execute_sequence(
            runner.algorithm.main_sequence,
            runner.state_grid, runner.algorithm_rng,
            runner.algorithm_state
        )
    end

    return nothing
end

gui_runner_is_finished(runner::GuiRunner)::Bool = isnothing(runner.algorithm_state)

function gui_main(runner::GuiRunner, delta_seconds::Float32)
    print_wnd_sizes::Bool = false && @markovjunior_debug(rand(Float32) < 0.01, false)
    pane_flags = CImGui.LibCImGui.ImGuiWindowFlags_HorizontalScrollbar |
                 (CImGui.LibCImGui.ImGuiWindowFlags_NoDecoration & (~CImGui.LibCImGui.ImGuiWindowFlags_NoScrollbar))

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

        should_update_texture = Ref(false)

        # Display a 'render mode' selection box.
        current_render2D_mode_idx = convert(Int32, Render2DMode.to_index(runner.render2D.mode)-1)
        next_render2D_mode_idx = current_render2D_mode_idx
        @c CImGui.ListBox("Mode", &current_render2D_mode_idx, RENDER_2D_MODE_NAMES, length(RENDER_2D_MODE_NAMES))
        if next_render2D_mode_idx != current_render2D_mode_idx
            should_update_texture[] = true
        end
        runner.render2D.mode = Render2DMode.from_index(convert(Int, current_render2D_mode_idx+1))
        # Display info about the current render mode.
        CImGui.SameLine(0, 15)
        gui_within_group() do
            if runner.render2D.mode == Render2DMode.normal
                CImGui.Dummy(1, 1)
            elseif runner.render2D.mode == Render2DMode.potentials
                CImGui.Text("Potential range: ")
                CImGui.Text(@sprintf("%f to %f",
                    min_inclusive(runner.render2D.inference_potentials_range),
                    max_inclusive(runner.render2D.inference_potentials_range)
                ))
            else
                error("Unhandled: ", runner.render2D.mode)
            end
        end

        CImGui.Separator()

        # Current state:
        img_size = convert(v2f, runner.render2D.output[].size.xy)
        min_img_size::Float32 = content_size.x - 10
        scale::Float32 = max(1.0f0, (min_img_size / img_size)...)
        @set! img_size *= scale
        CImGui.Image(gui_tex_handle(runner.render2D.output[]),
                     convert(gVec2, img_size),
                     gVec2(0, 0), gVec2(1, 1),
                     gVec4(1, 1, 1, 1), gVec4(0, 0, 0, 0))

        CImGui.Separator()

        # Execute Play logic.
        if runner.is_playing && !gui_runner_is_finished(runner)
            runner.time_till_next_tick -= delta_seconds
            while runner.time_till_next_tick <= 0
                step_gui_runner_algo(runner)
                should_update_texture[] = true
                runner.time_till_next_tick += 1.0f0 / runner.ticks_per_second
            end
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
                step_gui_runner_algo(runner)
                should_update_texture[] = true
            end
        CImGui.TableNextColumn()
        CImGui.TableNextColumn()
        CImGui.TableNextColumn()
            CImGui.SetNextItemWidth(-1)
            gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, BUTTON_COLOR_RUN_SPECIAL) do
                if CImGui.Button("Reset", BUTTON_SIZE_RUN_SPECIAL)
                    reset_gui_runner_algo(runner, false, false, false)
                    should_update_texture[] = true
                end
            end
        CImGui.TableNextColumn()
        CImGui.TableNextColumn()
        # Second row:
        CImGui.TableNextColumn()
            CImGui.Dummy(0, BUTTON_VPAD_RUN_PLAIN)
            CImGui.SetNextItemWidth(-1)
            if CImGui.Button("Jump", BUTTON_SIZE_RUN_PLAIN)
                for i in 1:runner.ticks_per_jump
                    step_gui_runner_algo(runner)
                    should_update_texture[] = true
                end
            end
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragInt("##TicksPerJump", &runner.ticks_per_jump, 1.0, 1, 0, "%d")
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("Ticks")
        CImGui.TableNextColumn()
            CImGui.SetNextItemWidth(-1)
            if gui_with_style(() -> CImGui.Button("Profile", BUTTON_SIZE_RUN_SPECIAL),
                              CImGui.LibCImGui.ImGuiCol_Button, BUTTON_COLOR_RUN_SPECIAL)
            #begin
                Profile.start_timer()
                for i in 1:runner.ticks_for_profile
                    step_gui_runner_algo(runner)
                    if gui_runner_is_finished(runner)
                        break
                    end
                end
                Profile.stop_timer()

                prof_text_path = joinpath(@__DIR__, "..", "ProfileResult.txt")
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
            @c CImGui.DragInt("##TicksForProfile", &runner.ticks_for_profile, 1.0, 1, 0, "%d")
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("ticks")
        # Third row:
        CImGui.TableNextColumn()
            CImGui.Dummy(0, BUTTON_VPAD_RUN_PLAIN)
            CImGui.SetNextItemWidth(-1)
            new_is_playing = Ref(runner.is_playing)
            gui_with_style(LCIG.ImGuiCol_Button, new_is_playing[] ? v3f(0.3, 0.8, 0.5) : v3f(0.8, 0.2, 0.05)) do
                if CImGui.Button(new_is_playing[] ? "Pause" : "Play", BUTTON_SIZE_RUN_PLAIN)
                    new_is_playing[] = !new_is_playing[]
                end
            end
            if new_is_playing[] && !runner.is_playing
                runner.time_till_next_tick = 1.0f0 / runner.ticks_per_second
            end
            runner.is_playing = new_is_playing[]
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragFloat("##TicksPerSecond", &runner.ticks_per_second,
                                0.1, 0.00001, 0, "%.0f", 1.0)
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.Text("ticks/s")
        CImGui.TableNextColumn()
            CImGui.SetNextItemWidth(-1)
            if gui_with_style(() -> CImGui.Button("Run to End", BUTTON_SIZE_RUN_SPECIAL),
                              LCIG.ImGuiCol_Button, BUTTON_COLOR_RUN_SPECIAL)
            #begin
                start_t = time()
                while !gui_runner_is_finished(runner)
                    step_gui_runner_algo(runner)
                    should_update_texture[] = true
                    if (time() - start_t) > runner.max_seconds_for_run_to_end
                        runner.algorithm_error_msg = string(
                            "ENDLESS RUN DETECTED: took longer than ",
                            runner.max_seconds_for_run_to_end,
                            " seconds to complete the grid!"
                        )
                        break
                    end
                end
            end
        CImGui.TableNextColumn()
            CImGui.Dummy(0, UNITS_VPAD)
            CImGui.SetNextItemWidth(-1)
            @c CImGui.DragFloat(
                "##MaxSecondsRunningToEnd",
                &runner.max_seconds_for_run_to_end,
                0.1,
                0.0, 0.0,
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
            CImGui.Text(runner.current_seed_display)
            CImGui.Dummy(10, 0); CImGui.SameLine()
            gui_with_item_width(150) do
                gui_text!(runner.next_seed)
            end
            CImGui.SameLine(0, 5)
            CImGui.Text(new_seed_is_int ? "(is number)" : "(is string)")
        end
        seed_has_changed::Bool = (new_seed != runner.current_seed)
        CImGui.SameLine(0, 40)
        if gui_with_style(() -> CImGui.Button("Reset with new seed", v2f(150, 32)) && seed_has_changed,
                          LCIG.ImGuiCol_Button, v3f(0.2, 0.1, 0.1),
                          unchanged = seed_has_changed)
        #begin
            reset_gui_runner_algo(runner, true, false, false)
            should_update_texture[] = true
        end

        CImGui.Separator()

        # Resolution data:
        fixed_dims = markov_fixed_dimension(runner.algorithm)
        fixed_resolution::Optional{Tuple} = markov_fixed_resolution(runner.algorithm)
        dims = convert(Int32, get_something(fixed_dims, runner.next_dimensionality))
        resolution::Vector = if exists(fixed_resolution)
            collect(fixed_resolution)
        else
            runner.next_resolution
        end
        gui_with_item_width(60) do
            if exists(fixed_dims)
                CImGui.LabelText("Dimension", string(dims))
            else
                @c CImGui.InputInt("Dimension", &dims)
                runner.next_dimensionality = clamp(dims, 2, 3)
                dims = runner.next_dimensionality
            end
            CImGui.SameLine(0, 30)
            # Update the resolution to match the dimensions.
            while length(runner.next_resolution) < dims
                push!(runner.next_resolution, 1)
            end
            while length(runner.next_resolution) > dims
                deleteat!(runner.next_resolution, length(runner.next_resolution))
            end
            if exists(fixed_resolution)
                CImGui.LabelText("Resolution", string(fixed_resolution))
            elseif dims == 2
                CImGui.InputInt2("Resolution", Ref(resolution, 1))
            elseif dims == 3
                CImGui.InputInt3("Resolution", Ref(resolution, 1))
            else
                CImGui.Text("ERROR: Unhandled: $dims/$resoution")
            end
        end
        resolution_is_different::Bool = (dims != ndims(runner.state_grid)) ||
                                        any(t->t[1]!=t[2], zip(resolution, size(runner.state_grid)))
        CImGui.SameLine(0, 10)
        gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, v3f(0.2, 0.1, 0.1), unchanged=resolution_is_different) do
         gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonHovered, v3f(0.2, 0.1, 0.1), unchanged=resolution_is_different) do
          gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonActive, v3f(0.2, 0.1, 0.1), unchanged=resolution_is_different) do
            if CImGui.Button("Reset with new size", v2f(150, 25)) && resolution_is_different
                reset_gui_runner_algo(runner, false, false, true)
                should_update_texture[] = true
            end
        end end end

        # Update the state texture, if any above code changed the state.
        if should_update_texture[]
            update_gui_runner_texture_2D(runner)
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
        gui_within_group() do
            for (color, greyscale, text) in GUI_LEGEND_DATA
                gui_draw_rect(
                    GuiDrawCursorRelative(Box2Df(
                        min=v2f(0, 0),
                        size=v2f(15, 15)
                    ), true),
                    GuiDrawFilled(color)
                )
                CImGui.SameLine()
                CImGui.Text(text)
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
    gui_window("Files", C_NULL, pane_flags) do
        BUTTON_SIZE = v2f(100, 35)
        CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Files")
        if CImGui.Button("Refresh", BUTTON_SIZE)
            update_gui_runner_scenes(runner)
        end
        CImGui.SameLine(0, 20)
        gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
         gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonHovered, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
          gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonActive, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
            if CImGui.Button("Reset changes", BUTTON_SIZE) && runner.current_scene_has_changes
                update!(runner.next_algorithm,
                        read(joinpath(scenes_path(), runner.current_scene), String))
                runner.current_scene_has_changes = false
            end
            CImGui.SameLine(0, 20)
            if CImGui.Button("Save changes", BUTTON_SIZE) && runner.current_scene_has_changes
                open(io -> print(io, string(runner.next_algorithm)),
                      joinpath(scenes_path(), runner.current_scene),
                      "w")
                runner.current_scene_has_changes = false
            end
        end end end

        next_scene_idx_c = convert(Int32, runner.current_scene_idx - 1)
        @c CImGui.ListBox("Scenes", &next_scene_idx_c,
                          runner.available_scenes, length(runner.available_scenes))
        if (next_scene_idx_c+1 != runner.current_scene_idx) && !runner.current_scene_has_changes
            runner.current_scene_idx = next_scene_idx_c+1
            runner.current_scene = runner.available_scenes[runner.current_scene_idx]
            update!(runner.next_algorithm, read(joinpath(scenes_path(), runner.current_scene), String))
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
        end

        if CImGui.Button("Restart##WithNewAlgorithm")
            reset_gui_runner_algo(runner, false, true, true)
            update_gui_runner_texture_2D(runner)
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