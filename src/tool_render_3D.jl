module Render3D

using Setfield, CSyntax
using CImGui
using Bplus
@using_bplus

using ..MarkovJunior: path_asset,
                      CellTypeSet, N_CELL_TYPES,
                      CellGrid, CellGridConcrete


"All the data available to shaders, in one UBO"
@std140 struct UboData
    grid_3D::UInt64

    matrix_view_proj::fmat4x4
    cam_pos::v3f

    sun_dir::v3f
    sun_color::v3f
    sun_shadowmap::UInt64
    sun_shadowmap_mat_viewproj::fmat4x4
    sun_shadowmap_mat_world_to_texel::fmat4x4

    shadowmap_world_bias::Float32

    cell_air_lookup::StaticBlockArray{N_CELL_TYPES, Bool}
end
const UBO_DATA_SHADER_SRC = """
    layout(std140, binding=0) uniform UniformBlock {
        $(glsl_decl(UboData))
    } u_data;
"""


###########
#   App

"Singleton rendering resources, across the entire app"
mutable struct App
    shader_src_utils::String
    shader_src_lighting::String

    render_cubes_depth::Program
    render_cubes_forward::Program
    visualize_depth_tex::Program

    dummy_depth_tex::Texture
end
function App()
    dummy_depth = Texture(DepthStencilFormats.depth_16u, fill(1.0f0, 1, 1),
                          sampler=TexSampler{2}(
                            depth_comparison_mode = ValueTests.greater_than_or_equal
                          ))

    make_dummy_shader() = GL.bp_glsl_str("""
        #START_VERTEX
            void main() { gl_Position = vec4(0); }
        #START_FRAGMENT
            void main() { }
    """)

    app = App(
        # Put in dummy shader data, then immediately call recompile_shaders!() below.
        "", "",
        make_dummy_shader(), make_dummy_shader(),
        make_dummy_shader(),
        dummy_depth
    )

    recompile_shaders!(app)
    return app
end
Base.close(a::App) = close.((
    a.render_cubes_depth,
    a.render_cubes_forward,
    a.visualize_depth_tex,
    a.dummy_depth_tex
))

function recompile_shaders!(a::App)
    # Only close the old ones if the new ones compile successfully.
    function replace_shader(field::Symbol, new_src::String)
        new_prog = try
            GL.bp_glsl_str(new_src)
        catch e
            error("Failed to recompile app.", field, ": ", sprint(showerror, e))
        end
        close(getproperty(a, field))
        setproperty!(a, field, new_prog)
    end

    a.shader_src_utils = read(path_asset("render3D/utils.glsl"), String)
    a.shader_src_lighting = read(path_asset("render3D/lighting.glsl"), String)

    shader_src_render_cubes = """
        #line 10000
        $(a.shader_src_utils)
        #line 20000
        $(a.shader_src_lighting)

        #line 30000
        $UBO_DATA_SHADER_SRC

        #line 1
        $(read(path_asset("render3D/render_cubes.glsl"), String))
    """
    replace_shader(:render_cubes_depth, """
        #define DEPTH_ONY 1
        $shader_src_render_cubes
    """)
    replace_shader(:render_cubes_forward, """
        #define DEPTH_ONY 0
        $shader_src_render_cubes
    """)

    shader_src_viz_depth_tex = read(path_asset("render3D/render_depth_map.glsl"), String)
    shader_src_viz_depth_tex = replace(shader_src_viz_depth_tex,
        "@blit_vs" => """
            uniform vec2 u_ndcPos;
            uniform vec2 u_ndcSize;
            $(make_vertex_shader_blit(
                "pos4.xy = (pos4.xy * u_ndcSize) + u_ndcPos;"
            ))
        """
    )
    replace_shader(:visualize_depth_tex, """
        #line 10000
        $(a.shader_src_utils)

        #line 1
        $shader_src_viz_depth_tex
    """)

    return nothing
end

"
Draws the given depth-texture using the depth-visualization shader,
  taking a time parameter for the shader.
"
function visualize_depth_tex(app::App, tex::Texture, out_target_uv::Box2Df, total_seconds::Float32)
    set_uniform(app.visualize_depth_tex, "u_time", total_seconds)
    set_uniform(app.visualize_depth_tex, "u_tex", get_view(tex, TexSampler{2}(
        pixel_filter = PixelFilters.rough,
        wrapping = WrapModes.clamp
    )))
    set_uniform(app.visualize_depth_tex, "u_ndcPos", convert(v2f, (center(out_target_uv) - 0.5) * 2))
    set_uniform(app.visualize_depth_tex, "u_ndcSize", convert(v2f, size(out_target_uv)  * 2))
    custom_blit(tex, app.visualize_depth_tex)
end


#############
#   Scene

"
Bare minimum viewport into the world (e.g. for a depth-only view).
Used directly to represent shadow-casters, or pseudo-inherited from for real renderable viewports.
"
mutable struct BasicViewport
    cam::Cam3D{Float32}
    view_depth::Texture
    view_target::Target

    # If true, backfaces are drawn instead of front-faces.
    flip_face_culling::Bool

    render_buffer_cpu::UboData
    render_buffer_gpu::Buffer
end
function BasicViewport(resolution::Vec2{<:Integer}, start_pos::Vec3
                       ;
                       color_outputs::Vector{TargetOutput} = TargetOutput[ ],
                       projection::Union{PerspectiveProjection{Float32}, OrthographicProjection{Float32}} =
                         PerspectiveProjection{Float32}(),
                       flip_face_culling::Bool = false)
    return BasicViewport(
        Cam3D{Float32}(
            pos = convert(v3f, start_pos),
            forward = -vnorm(convert(v3f, start_pos)),
            projection = projection
        ),

        begin
            view_depth = Texture(
                DepthStencilFormats.depth_24u, resolution,
                n_mips=1
            )
            view_target = Target(color_outputs, TargetOutput(tex=view_depth))
            (view_depth, view_target)
        end...,

        flip_face_culling,

        UboData(), Buffer(true, UboData)
    )
end
Base.close(v::BasicViewport) = close.((
    v.view_target, v.view_depth,
    v.render_buffer_gpu
))


"Settings and resources for a single MarkovJunior scene"
mutable struct Scene
    grid_tex_3D::Texture
    grid_buffer_3D::CellGridConcrete{3}

    air_cells::CellTypeSet
    sun_color_hdr::v3f
    sun_dir::v3f
    sun_dir_gui_buffer::Ref{Float32}
    #TODO: Floor, Sky, Fog, Per-pixel Materials

    sun_pov::BasicViewport
    sun_shadowmap_viewproj::fmat4x4
    sun_shadowmap_world_to_texel::fmat4x4

    old_textures::Vector{Pair{Texture, Int}} # Elements are cleaned up after X frames,
                                             #   to prevent crashes from Dear ImGUI trying to draw with it.
end
function Scene(; sun_dir::v3f = norm(v3f(1, -1, -1)),
                 sun_color::v3f = v3f(1, 1, 1),
                 sun_shadow_resolution::Int = 2048)
    return Scene(
        Texture(
            SimpleFormat(FormatTypes.uint, SimpleFormatComponents.R, SimpleFormatBitDepths.B8),
            v3u(1, 1, 1),
            sampler=TexSampler{3}(
                pixel_filter=PixelFilters.rough
            ),
            n_mips=1
        ),
        fill(zero(UInt8), 1, 1, 1),

        CellTypeSet('b'),
        sun_color, sun_dir,
        Ref{Float32}(0),

        # Note that none of the arguments to the sun's viewport are important here
        #    except for its resolution.
        BasicViewport(v2i(i->sun_shadow_resolution), zero(v3f),
                      flip_face_culling=true),
        m_identityf(4, 4), m_identityf(4, 4),

        Vector{Pair{Texture, Int}}()
    )
end
Base.close(s::Scene) = close.((
    s.grid_tex_3D, s.sun_pov,
    (p[1] for p in s.old_textures)...
))

function update_scene_grid!(scene::Scene, new_grid_view::CellGrid{3})
    if scene.grid_tex_3D.size.xyz == vsize(new_grid_view)
        scene.grid_buffer_3D .= new_grid_view
        set_tex_color(scene.grid_tex_3D, scene.grid_buffer_3D)
    else
        scene.grid_buffer_3D = fill(zero(UInt8), size(new_grid_view)...)
        scene.grid_buffer_3D .= new_grid_view

        push!(scene.old_textures, scene.grid_tex_3D => 2)
        scene.grid_tex_3D = Texture(
            scene.grid_tex_3D.format,
            scene.grid_buffer_3D,
            sampler=scene.grid_tex_3D.sampler,
            n_mips=1
        )
    end
    return nothing
end

function tick_scene!(scene::Scene, delta_seconds::Float32, app::App)
    # Clean up old resources.
    i = length(scene.old_textures)
    while i > 0
        if (scene.old_textures[i][2] < 1)
            close(scene.old_textures[i][1])
            deleteat!(scene.old_textures, i)
        else
            scene.old_textures[i] = scene.old_textures[i][1] => (scene.old_textures[i][2] - 1)
        end
        i -= 1
    end

    # Update the sun shadow-map.
    #  1. Get the view matrix for the sun's POV.
    scene.sun_pov.cam = Cam3D{Float32}(
        pos = zero(v3f),
        forward = scene.sun_dir,
        projection = Box3Df(center=zero(v3f), size=zero(v3f)) # Filled in later
    )
    m_sun_view = cam_view_mat(scene.sun_pov.cam )
    #  2. Transform the level bounds to the sun's view-space.
    level_corners_world = map(v -> v * scene.grid_tex_3D.size.xyz, (
        v3f(0, 0, 0),
        v3f(1, 0, 0),
        v3f(0, 1, 0),
        v3f(0, 0, 1),
        v3f(1, 1, 0),
        v3f(0, 1, 1),
        v3f(1, 0, 1),
        v3f(1, 1, 1)
    ))
    level_corners_view = map(v -> m_apply_point_affine(m_sun_view, v), level_corners_world)
    #  3. Get an ortho projection for the sun that fits the level bounds inside its view.
    scene.sun_pov.cam = let c = scene.sun_pov.cam
        @set! c.projection = Box3Df(
            min=min(level_corners_view...),
            max=max(level_corners_view...)
        )
        c
    end
    m_sun_projection = cam_projection_mat(scene.sun_pov.cam)
    #  4. Get the final matrices for sun shadowmapping.
    scene.sun_shadowmap_viewproj = m_combine(m_sun_view, m_sun_projection)
    scene.sun_shadowmap_world_to_texel = m_combine(
        scene.sun_shadowmap_viewproj,
        # Map NDC space (-1, +1) to texel space (0, 1).
        m_scale(v4f(0.5, 0.5, 1, 1)),
        m4_translate(v3f(0.5, 0.5, 0))
    )
    #  5. Render into the shadowmap depth buffer.
    render(app, scene, scene.sun_pov)
end

"Runs some Dear ImGUI widgets to edit this scene's settings"
function scene_settings_gui!(scene::Scene)
    CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Sun")
    @c CImGui.ColorEdit3("Color##Sun", &scene.sun_color_hdr,
        CImGui.LibCImGui.ImGuiColorEditFlags_HDR |
          CImGui.LibCImGui.ImGuiColorEditFlags_Float
    )
    gui_with_item_width(100) do
        scene.sun_dir = gui_spherical_vector("Direction##Sun", scene.sun_dir,
                                             stays_normalized=true,
                                             fallback_yaw=scene.sun_dir_gui_buffer)
    end
end


##############
#   Viewport

"
A fully-featured viewport into the world, pseudo-inheriting from `BasicViewport`
  (takes its fields and is automatically convertible to it).

Explicitly convert it to a `BasicViewport` with `v.base`.
"
mutable struct FullViewport
    base::BasicViewport
    cam_settings::Cam3D_Settings{Float32}

    shadowmap_world_bias::Float32

    view_color::Texture
    view_target_outputs_depth_only::Vector{Optional{Int}}
    view_target_outputs_forward::Vector{Optional{Int}}
end
function FullViewport(resolution::Vec2{<:Integer}, start_pos::Vec3)
    view_color = Texture(
        SpecialFormats.rgb10_a2, resolution,
        sampler = TexSampler{2}(
            pixel_filter=PixelFilters.smooth,
            wrapping=WrapModes.clamp
        ),
        n_mips=1
    )
    return FullViewport(
        BasicViewport(
            resolution, start_pos,
            color_outputs = [ TargetOutput(tex=view_color) ]
        ),
        Cam3D_Settings{Float32}(),

        -0.1f0,

        view_color,
        Optional{Int}[ ],
        Optional{Int}[ 1 ]
    )
end
Base.close(v::FullViewport) = close.((
    v.base, v.view_color
))
Base.convert(::Type{BasicViewport}, v::FullViewport) = v.base
Base.propertynames(v::FullViewport) = (
    propertynames(v.base)...,
    fieldnames(FullViewport)...
)
@inline Base.getproperty(v::FullViewport, n::Symbol) = if n in fieldnames(BasicViewport)
    getproperty(getfield(v, :base), n)
else
    getfield(v, n)
end
@inline Base.setproperty!(fv::FullViewport, n::Symbol, v) = if n in fieldnames(BasicViewport)
    setproperty!(getfield(fv, :base), n, v)
else
    setfield!(fv, n, v)
end



"Resets this viewport's position to oversee a new grid of the given resolution"
function on_new_grid!(viewport::FullViewport, new_grid_size::Vec3{<:Integer})
    viewport.cam = let c = viewport.cam
        @set! c.pos = vappend(maximum(new_grid_size.xy) * v2f(i->1.25f0),
                              convert(Float32, (sum(new_grid_size.data)/3) * 1.5f0))
        @set! c.forward = -vnorm(c.pos)
        c
    end
    return nothing
end

"Runs some Dear ImGUI widgets to edit this viewport's settings"
function viewport_settings_gui!(viewport::FullViewport, scene::Scene)
    CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Camera")
    viewport.cam = let c = viewport.cam
        p = c.pos
        @c CImGui.DragFloat3("Pos##Camera", &p, 1.0)
        @set! c.pos = p

        f = c.forward
        f = gui_with_item_width(100) do
            gui_spherical_vector("Direction##Camera", f,
                                 stays_normalized=true)
        end
        @set! c.forward = f

        c
    end

    CImGui.Separator(); CImGui.SameLine(30); CImGui.Text("Viewport")
    @c CImGui.DragFloat("Shadow World Bias", &viewport.shadowmap_world_bias,
                        0.1, -10, 10, "%.2f")
end

function render(app::App, scene::Scene, view::Union{BasicViewport, FullViewport})
    # Generate the UBO data.
    view.render_buffer_cpu.cam_pos = view.cam.pos
    view.render_buffer_cpu.grid_3D = GL.gl_type(get_ogl_handle(get_view(scene.grid_tex_3D)))
    view.render_buffer_cpu.matrix_view_proj = m_combine(
        cam_view_mat(view.cam),
        cam_projection_mat(view.cam)
    )
    view.render_buffer_cpu.sun_dir = scene.sun_dir
    view.render_buffer_cpu.sun_color = scene.sun_color_hdr
    # If this viewport *is* the sun's shadowmap, use a dummy shadowmap here.
    sun_shadowmap_tex = if view == scene.sun_pov
        app.dummy_depth_tex
    else
        scene.sun_pov.view_depth
    end
    sun_shadowmap_view = get_view(sun_shadowmap_tex, TexSampler{2}(
        pixel_filter = PixelFilters.smooth,
        wrapping = WrapModes.clamp,
        depth_comparison_mode = ValueTests.less_than_or_equal
    ))
    view.render_buffer_cpu.sun_shadowmap = GL.gl_type(get_ogl_handle(sun_shadowmap_view))
    view.render_buffer_cpu.sun_shadowmap_mat_viewproj = scene.sun_shadowmap_viewproj
    view.render_buffer_cpu.sun_shadowmap_mat_world_to_texel = scene.sun_shadowmap_world_to_texel
    view.render_buffer_cpu.shadowmap_world_bias = if view isa FullViewport
        view.shadowmap_world_bias
    else
        0.0f0
    end
    for i in 1:N_CELL_TYPES
        view.render_buffer_cpu.cell_air_lookup[i] = (i-1) in scene.air_cells
    end

    # Upload the UBO data.
    set_buffer_data(view.render_buffer_gpu, view.render_buffer_cpu)
    set_uniform_block(view.render_buffer_gpu, 1) #NOTE: GLSL index is 0-based, while this is 1-based!

    # Set up rendering state.
    set_render_state(RenderState(
        depth_test = ValueTests.less_than,
        cull_mode = view.flip_face_culling ? FaceCullModes.backwards : FaceCullModes.on
    ))
    view_activate(scene.grid_tex_3D)
    view_activate(sun_shadowmap_view)
    target_activate(view.view_target)
    n_cubes = prod(scene.grid_tex_3D.size.xyz)
    n_faces = n_cubes * 6
    n_tris = n_faces * 2
    n_verts = n_tris * 3
    empty_mesh = service_BasicGraphics().empty_mesh

    # Do the depth pre-pass.
    if view isa FullViewport
        target_configure_fragment_outputs(view.view_target, view.view_target_outputs_depth_only)
    end
    target_clear(view.view_target, @f32(1))
    render_mesh(
        empty_mesh, app.render_cubes_depth,
        shape=PrimitiveTypes.triangle,
        elements = IntervalU(
            min=1,
            size=n_verts
        )
    )

    # Do the forward pass.
    if view isa FullViewport
        target_configure_fragment_outputs(view.view_target, view.view_target_outputs_forward)
        set_depth_test(ValueTests.less_than_or_equal)
        target_clear(view.view_target, v4f(0.7, 0.7, 1, 1))
        render_mesh(
            empty_mesh, app.render_cubes_forward,
            shape=PrimitiveTypes.triangle,
            elements = IntervalU(
                min=1,
                size=n_verts
            )
        )
    end

    # Clean up.
    view_deactivate(sun_shadowmap_view)
    view_deactivate(scene.grid_tex_3D)
    target_activate(nothing)
    return nothing
end


end # module