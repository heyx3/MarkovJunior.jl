function setup_gui_tool_inputs()
    create_axis("CAM: forward", AxisInput([
        ButtonAsAxis(GLFW.KEY_W),
        ButtonAsAxis_Negative(GLFW.KEY_S)
    ]))
    create_axis("CAM: upward", AxisInput([
        ButtonAsAxis(GLFW.KEY_E),
        ButtonAsAxis_Negative(GLFW.KEY_Q)
    ]))
    create_axis("CAM: rightward", AxisInput([
        ButtonAsAxis(GLFW.KEY_D),
        ButtonAsAxis_Negative(GLFW.KEY_A)
    ]))
    create_axis("CAM: pitch", AxisInput(
        MouseAxes.y, AxisModes.delta;
        value_scale=-0.05
    ))
    create_axis("CAM: yaw", AxisInput(
        MouseAxes.x, AxisModes.delta;
        value_scale=0.05
    ))
end
get_gui_tool_cam_input() = Cam3D_Input{Float32}(
    controlling_rotation = gui_tool_imgui_denying_mouse_control(),
    forward = get_axis("CAM: forward"),
    right = get_axis("CAM: rightward"),
    up = get_axis("CAM: upward"),
    pitch = get_axis("CAM: pitch"),
    yaw = get_axis("CAM: yaw")
)

gui_tool_imgui_denying_mouse_control()::Bool = let io = CImGui.GetIO()
    unsafe_load(io.WantCaptureKeyboard) || unsafe_load(io.WantCaptureMouse)
end


function markovjunior_run_gui()
    @game_loop begin
        INIT(
            v2i(1210, 800), "Markov Junior Playground",
            vsync=VsyncModes.on,
            glfw_hints = Dict(
                GLFW.FOCUSED => 1
            )
        )

        SETUP = begin
            LOOP.max_fps = 165
            LOOP.max_frame_duration = 0.1 # Prevent large camera movements during hitches

            # Put the Dear ImGUI config inside the locals directory.
            # We must store the path string in a variable so the underlying char array can't be GC-ed.
            imgui_config_path = path_local("imgui.ini")
            imgui_config_path_c = Base.unsafe_convert(Ptr{Int8}, imgui_config_path)
            unsafe_store!(CImGui.GetIO().IniFilename, imgui_config_path_c)

            # Load user session memory from file, if available.
            memory_path = path_local(MEMORY_FILE_NAME)
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
                get_editor_font_bytes(), [ 19 ]
            )[1]
            gui = GuiRunner(memory, gui_editor_font)

            setup_gui_tool_inputs()
        end

        LOOP = begin
            if GLFW.WindowShouldClose(LOOP.context.window)
                break
            end

            gui_main(gui, LOOP.delta_seconds, get_gui_tool_cam_input())
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
    return nothing
end
"The entry-point for standalone executables"
function markovjunior_run_gui_main()::Cint
    try
        markovjunior_run_gui()
        return 0
    catch e
        showerror(stderr, e, catch_backtrace())
        return 1
    end
end