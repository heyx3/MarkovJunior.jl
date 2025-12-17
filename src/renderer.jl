"
Rendering a 2D markov-junior scene is easy: copy the cell colors into an array.
You may pass a null ref the first time you call this.
This function also takes a reusable buffer to eliminate heap allocations.
"
function render_markov_2d(grid::CellGrid{2}, null_color::v3f,
                          buffer::Ref{Matrix{v3f}}, output::Ref{Texture})
    # Resize/initialize the output texture if necessary.
    if !isassigned(output) || (output[].size.xy != vsize(grid))
        output[] = Texture(
            SimpleFormat(FormatTypes.normalized_uint,
                         SimpleFormatComponents.RGB,
                         SimpleFormatBitDepths.B8),
            vsize(grid)::Vec2,
            sampler = Bplus.GL.TexSampler{2}(
                pixel_filter = PixelFilters.rough
            )
        )
    end

    # Build the texture's array.
    if !isassigned(buffer) || (vsize(buffer[]) != vsize(grid))
        buffer[] = fill(zero(v3f), size(grid))
    end
    buffer_memory::Matrix{v3f} = buffer[]
    for pixel::v2i in one(v2i):vsize(buffer_memory)
        cell::UInt8 = grid[pixel]
        color::v3f = if cell == CELL_CODE_INVALID
            null_color
        else
            CELL_TYPES[cell + 1].color
        end
        buffer_memory[pixel] = color
    end

    # Upload.
    Bplus.GL.set_tex_color(output[], buffer_memory)
end


Bplus.GL.@std140 struct CellColorLookupUBO
    colors::StaticBlockArray{length(CELL_TYPES), v3f}
end
const UBO_IDX_CELL_COLOR_LOOKUP = 1

"Static process-wide resources for rendering markov grids"
struct MarkovRendererResources
    cell_color_lookup::Bplus.GL.Buffer
    program_render_cube::Bplus.GL.Program
end

# The shader source references some constants so it must be evaluated within this module.
const SHADER_SOURCE_RENDER_CUBES = eval(string(
    "\"\"\"",
    open(io -> read(io, String),
         joinpath(@__FILE__, "..", "CubeRender.glsl")),
    "\"\"\""
))

function MarkovRendererResources()
    r = MarkovRendererResources(
        CellColorLookupUBO(map(c -> c.color, CELL_TYPES)),
        Bplus.GL.bp_glsl_str(SHADER_SOURCE_RENDER_CUBES)
    )
    Bplus.GL.set_uniform_block(r.cell_color_lookup, UBO_IDX_CELL_COLOR_LOOKUP)
    return r
end


struct MarkovRenderer3D{TGrid <: CellGrid{3}}
    grid::TGrid
    state_tex::Bplus.GL.Texture

    cam::Cam3D{Float32}

    output_color::Bplus.GL.Texture
    output_target::Bplus.GL.Target
end

function MarkovRenderer3D(grid::CellGrid{3}, resolution::v2u)
    grid_size = vsize(grid)
    output_color = Texture(SpecialFormats.rgb_tiny_ufloats, resolution)
    cam_back_up_scale = @f32(2)
    return MarkovRenderer3D{typeof(grid)}(
        grid,
        Texture(
            SimpleFormat(
                FormatTypes.uint,
                SimpleFormatComponents.R,
                SimpleFormatBitDepths.B8
            ),
            convert(v3u, grid_size),
            sampler=TexSampler{2}(wrapping = WrapModes.clamp)
        ),

        Cam3D{Float32}(
            pos = convert(v3f, grid_size) * cam_back_up_scale,
            forward = vnorm(v3f(-1, -1, -1)),
            up = v3f(0, 0, 1),
            projection = OrthographicProjection{Float32}(
                center=v3f(
                    0, 0,
                    convert(v3f, grid_size) * cam_back_up_scale / 2
                ),
                size=v3f(
                    (one(v2f) * @f32(max(grid_size) * 1.4))...,
                    vlength(grid_size * cam_back_up_scale)
                )
            )
        ),

        output_color,
        Target(TargetOutput(tex=output_color),
               DepthStencilFormats.depth_32u)
    )
end

Base.close(mr::MarkovRenderer3D) = close.((
    mr.grid_state,
    mr.output_color,
    mr.output_target
))

function render_markov_3d(grid::CellGrid{3}, renderer::MarkovRenderer3D, assets::MarkovRendererResources)
    set_tex_color(renderer.state_tex, grid)

    target_activate(renderer.output_target)
        clear_screen(vRGBAf(0.6, 0.7, 1.0, 1.0))
        clear_screen(@f32(1))

        set_uniform(assets.program_render_cube, "u_transformWorld", m_identityf(4, 4))
        set_uniform(assets.program_render_cube, "u_transformWorldNormals", m_identityf(3, 3))
        set_uniform(assets.program_render_cube, "u_transformViewProj", Bplus.Math.m_combine(
            Bplus.BplusTools.cam_view_mat(renderer.cam),
            Bplus.BplusTools.cam_projection_mat(renderer.cam)
        ))
        set_uniform(assets.program_render_cube, "u_gridTex", renderer.state_tex)
        set_uniform_block(assets.program_render_cube, "CellLookupUBO", UBO_IDX_CELL_COLOR_LOOKUP)

        view_activate(renderer.state_tex)
            render_mesh(
                service_BasicGraphics().empty_mesh,
                assets.program_render_cube,
                shape = PrimitiveTypes.point,
                elements = IntervalU(min=1, size=prod(size(grid)))
            )
        view_deactivate(renderer.state_tex)

    target_activate(nothing)
end