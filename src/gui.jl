@bp_bitflag(MarkovJuniorGuiRenderFlags,
    NONE,
    potentials, rules
)
const GUI_LEGEND_DATA = map(enumerate(CELL_TYPES)) do (i,t)
    return (
        t.color,
        Float32(i)/Float32(length(CELL_TYPES)),
        " - $(t.char)"
    )
end

const gVec2 = Bplus.GUI.gVec2
const gVec4 = Bplus.GUI.gVec4


"The state of the GUI for our MarkovJunior tool"
mutable struct GuiRunner
    draw_monochrome::Bool
    draw_flags::E_MarkovJuniorGuiRenderFlags

    next_dimensionality::Int
    next_resolution::Vector{Int}

    state_grid::CellGrid
    state_grid_tex2D_slice::CellGrid{2}
    state_texture::Texture

    algorithm::ParsedMarkovAlgorithm
    algorithm_state::Optional
    algorithm_rng::PRNG

    algorithm_render2D::Ref{Texture}
    algorithm_render2D_buffer::Ref{Matrix{v3f}}

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
        false, MarkovJuniorGuiRenderFlags.NONE,

        2, [ fake_size, fake_size ],
        fill(zero(UInt8), fake_size, fake_size),
        fill(zero(UInt8), fake_size, fake_size),
        Texture(SpecialFormats.rgb10_a2, v2u(1, 2)), # Size must not match fake_size
                                                     #   or else it won't be properly reallocated

        @markovjunior(begin end), nothing, PRNG(1),
        Ref{Texture}(), Ref{Matrix{v3f}}(),
        zero(UInt64), "[NULL]",
        GuiText(string(seed)),

        GuiText(initial_algorithm_str,
            is_multiline=true,
            imgui_flags=CImGui.LibCImGui.ImGuiInputTextFlags_AllowTabInput
        ),
        algorithm_editor_font,
        "Uh-OH: UNINITIALIZED!!!!",

        String[ ], "", 1, false,

        false, 5.0f0, -1.0f0,

        10, 1000, 10.0f0
    )

    reset_gui_runner_algo(runner, true, true)
    update_gui_runner_texture_2D(runner)
    update_gui_runner_scenes(runner)

    return runner
end
function Base.close(runner::GuiRunner)
    close(runner.state_texture)
    close(runner.algorithm_render2D[])
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

    render_markov_2d(runner.state_grid_tex2D_slice, v3f(0.2, 0.2, 0.2),
                     runner.algorithm_render2D_buffer,
                     runner.algorithm_render2D)

    return nothing
end

function reset_gui_runner_algo(runner::GuiRunner,
                               parse_new_seed::Bool, parse_new_algorithm::Bool)
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
                push!(runner.next_resolution, 1)
            end
            while length(runner.next_resolution) > n_dims
                deleteat!(runner.next_resolution, length(runner.next_resolution))
            end
        end
    end

    # Initialize the grid.
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

    # Always stop Playing if the algorithm is finished.
    if isnothing(runner.algorithm_state)
        runner.is_playing = false
    end

    return nothing
end

gui_runner_is_finished(runner::GuiRunner)::Bool = isnothing(runner.algorithm_state)

function gui_main(runner::GuiRunner, delta_seconds::Float32)
    print_wnd_sizes::Bool = @markovjunior_debug(rand(Float32) < 0.01, false)

    gui_next_window_space(
        Box2Df(
            min=v2f(0, 0),
            max=v2f(0.45, 1)
        ),
        min_pixel_size = v2i(323, -1)
    )
    gui_window("Runner", C_NULL, CImGui.LibCImGui.ImGuiWindowFlags_NoDecoration) do
        content_size = convert(v2f, CImGui.GetContentRegionAvail())
        print_wnd_sizes && println("Rnn: ", CImGui.GetWindowSize())

        # Render settings:
        @c CImGui.Selectable("Monochrome", &runner.draw_monochrome)
        CImGui.SameLine(0, 20)
        for (name, flag) in [ ("Potentials", MarkovJuniorGuiRenderFlags.potentials),
                                ("Rules", MarkovJuniorGuiRenderFlags.rules) ]
        #begin
            if runner.draw_monochrome
                if CImGui.Selectable(name, contains(flag, runner.draw_flags))
                    runner.draw_flags |= flag
                else
                    runner.draw_flags -= flag
                end
            else
                # Draw a disabled version of the widget.
                CImGui.Text(name)
            end
            CImGui.SameLine(0, 5)
        end
        CImGui.Dummy(0, 0) # To cancel the last SameLine() call

        # Current state:
        CImGui.BeginChild(CImGui.GetID("StateDisplayArea"),
                          gVec2(content_size.x - 20,
                                content_size.y - 200))
            img_size = convert(v2f, runner.algorithm_render2D[].size.xy)
            min_img_size::Float32 = content_size.x - 4
            scale::Float32 = max(1.0f0, (min_img_size / img_size)...)
            @set! img_size *= scale

            CImGui.Image(gui_tex_handle(runner.algorithm_render2D[]),
                         convert(gVec2, img_size),
                         gVec2(0, 0), gVec2(1, 1),
                         gVec4(1, 1, 1, 1), gVec4(0, 0, 0, 0))
        CImGui.EndChild()
        #TODO: Add B+ helper for scroll regions once this is verified working

        # Below actions may invalidate the algorithm state.
        should_update_texture = Ref(false)

        # Run buttons:
        #   * Step
        if CImGui.Button("Step")
            step_gui_runner_algo(runner)
            should_update_texture[] = true
        end
        CImGui.SameLine(0, 40)
        #   * Jump
        gui_with_item_width(20) do
            @c CImGui.DragInt("##TicksPerJump", &runner.ticks_per_jump, 1.0, 1, 0, "%d")
        end
        CImGui.SameLine()
        if CImGui.Button("Jump")
            for i in 1:runner.ticks_per_jump
                step_gui_runner_algo(runner)
                should_update_texture[] = true
            end
        end
        CImGui.SameLine(0, 40)
        #    * Play/Pause
        if runner.is_playing
            if gui_runner_is_finished(runner)
                runner.is_playing = false
            else
                runner.time_till_next_tick -= delta_seconds
                while runner.time_till_next_tick <= 0
                    step_gui_runner_algo(runner)
                    should_update_texture[] = true
                    runner.time_till_next_tick += 1.0f0 / runner.ticks_per_second
                end
            end
        end
        gui_with_item_width(30) do
            @c CImGui.DragFloat("##TicksPerSecond", &runner.ticks_per_second, 0.1, 0.00001, 0, "%.1f", 1.0)
        end
        CImGui.SameLine()
        new_is_playing = Ref(runner.is_playing)
        gui_with_item_width(60) do
            CImGui.Selectable(new_is_playing[] ? "Pause" : "Play",
                              new_is_playing)
        end
        if new_is_playing[] && !runner.is_playing
            runner.time_till_next_tick = 1.0f0 / runner.ticks_per_second
        end
        runner.is_playing = new_is_playing[]

        CImGui.Dummy(0, 20)

        # Special control buttons:
        gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, v3f(0.5, 0.3, 0.35)) do
            gui_with_item_width(30) do
                @c CImGui.DragFloat(
                    "##MaxSecondsRunningToEnd",
                    &runner.max_seconds_for_run_to_end,
                    0.1,
                    0.0, 0.0,
                    "%.1f", 1.0
                )
            end
            CImGui.SameLine()
            gui_tooltip("Max seconds, before canceling the run-to-end")
            if CImGui.Button("Run to End")
                start_t = time()
                while !gui_runner_is_finished(runner)
                    step_gui_runner_algo(runner)
                    should_update_texture[] = true

                    if (time() - start_t) > runner.max_seconds_for_run_to_end
                        runner.algorithm_error_msg = string(
                            "ENDLESS RUN DETECTED: took longer than ",
                            runner.max_seconds_for_run_to_end, " seconds to run to the end!",
                            "\n\nYou may increase this cutoff time if you like."
                        )
                        break
                    end
                end
            end
            CImGui.SameLine(0, 40)
            if CImGui.Button("Reset")
                reset_gui_runner_algo(runner, false, false)
                should_update_texture[] = true
            end
            CImGui.SameLine(0, 20)
            gui_with_item_width(20) do
                @c CImGui.DragInt("##TicksForProfile", &runner.ticks_for_profile, 1.0, 1, 0, "%d")
            end
            CImGui.SameLine()
            if CImGui.Button("Profile")
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
                    Profile.print(ctx, format=:flat, mincount=100)

                    println(ctx, "\n\n\n")

                    println(ctx, "========================================\n==    Tree")
                    Profile.print(ctx, format=:tree, noisefloor=2.0)
                end
                if Sys.iswindows()
                    run(`cmd /C $(abspath(prof_text_path))`)
                else
                    run(`$(abspath(prof_text_path))`)
                end
            end
            CImGui.SameLine(0, 10)
            CImGui.Text("#TODO: profiled modal view")
        end

        CImGui.Dummy(0, 20)

        # Seed data:
        CImGui.Text(runner.current_seed_display)
        CImGui.Dummy(20, 0); CImGui.SameLine()
        gui_with_item_width(100) do
            gui_text!(runner.next_seed)
        end
        CImGui.SameLine(0, 10)
        if CImGui.Button("Restart##WithNewSeed")
            reset_gui_runner_algo(runner, true, false)
            should_update_texture[] = true
        end
        CImGui.SameLine(0, 20)
        CImGui.Text(isnothing(tryparse(UInt64, string(runner.next_seed))) ?
                      "as String" :
                      "as number")

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
    gui_window("Legend", C_NULL, CImGui.LibCImGui.ImGuiWindowFlags_NoDecoration) do
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
    gui_window("Files", C_NULL, CImGui.LibCImGui.ImGuiWindowFlags_NoDecoration) do
        CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Files")
        if CImGui.Button("Refresh")
            update_gui_runner_scenes(runner)
        end
        CImGui.SameLine(0, 20)
        gui_with_style(CImGui.LibCImGui.ImGuiCol_Button, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
         gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonHovered, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
          gui_with_style(CImGui.LibCImGui.ImGuiCol_ButtonActive, v3f(0.2, 0.1, 0.1), unchanged=runner.current_scene_has_changes) do
            if CImGui.Button("Reset changes") && runner.current_scene_has_changes
                update!(runner.next_algorithm,
                        read(joinpath(scenes_path(), runner.current_scene), String))
                runner.current_scene_has_changes = false
            end
            CImGui.SameLine(0, 20)
            if CImGui.Button("Save changes") && runner.current_scene_has_changes
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
    gui_window("Editor", C_NULL, CImGui.LibCImGui.ImGuiWindowFlags_NoDecoration) do
        content_size = convert(v2f, CImGui.GetContentRegionAvail())
        print_wnd_sizes && println("Editor wnd: ", CImGui.GetWindowSize())

        runner.next_algorithm.multiline_requested_size = round.(Ref(Int),
            (content_size - v2f(20, 50)).data
        )
        gui_with_font(runner.next_algorithm_font) do
            runner.current_scene_has_changes |= gui_text!(runner.next_algorithm)
        end

        if CImGui.Button("Restart##WithNewAlgorithm")
            reset_gui_runner_algo(runner, false, true)
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

    #TODO: File management window
end