//Cubes will be rendered by dispatching 36 vertices per grid cell
//  and discarding ones for the cells which represent air.

#if (RENDER_PASS == RENDER_PASS_FORWARD) || (RENDER_PASS == RENDER_PASS_TRANSPARENT)
    #define DEPTH_ONLY 0
#else
    #define DEPTH_ONLY 1
#endif


#START_VERTEX

#if !DEPTH_ONLY
out ivec3 o_gridCell;
out int o_cellValue;
out vec3 o_gridPosF;
out vec3 o_normal;
out vec2 o_faceUV;
#endif

//You can discard primitives by positioning their vertices at NaN.
#define MAKE_NAN32 (intBitsToFloat(int(0xFFC00000u)))

void main() {
    //Figure out what vertex this is.
    int cubeIdx = (gl_VertexID / 36),    //Which cube?
        faceIdx = (gl_VertexID / 6) % 6, //Which face of cube?
        triIdx  = (gl_VertexID / 3) % 2, //Which tri of face?
        vertIdx = (gl_VertexID % 3),     //Which vert of tri?
        faceVertIdx = (gl_VertexID % 6); //Which vert of face?
    int faceAxis = (faceIdx / 2),
        faceDirMask = (faceIdx % 2),
        faceDir = -1 + (2 * faceDirMask);

    //Read this grid cell.
    usampler3D gridTex = usampler3D(u_data.grid_3D);
    ivec3 gridResolution = textureSize(gridTex, 0);
    ivec3 gridCell = ivec3(
        cubeIdx % gridResolution.x,
        (cubeIdx / gridResolution.x) % gridResolution.y,
        (cubeIdx / (gridResolution.x * gridResolution.y)) % gridResolution.z
    );
    int cellValue = int(texelFetch(gridTex, gridCell, 0).r);
    UboData_Material cellMat = u_data.cell_materials[cellValue];

    //If this cell is meant to be empty, discard the primitive.
    if (cellMat.mode == CELL_MATERIAL_EMPTY ||
        ((RENDER_PASS == RENDER_PASS_TRANSPARENT) != (cellMat.mode == CELL_MATERIAL_GLASS)))
    {
        #if !DEPTH_ONLY
            o_cellValue = cellValue;
            o_gridCell = gridCell;
            o_faceUV = vec2(0, 0);
            o_gridPosF = vec3(0, 0, 0);
            o_normal = vec3(0, 0, 0);
        #endif

        gl_Position = vec4(MAKE_NAN32, MAKE_NAN32, MAKE_NAN32, -1);
        return;
    }

    //Calculate vertex data.
    bool flipWinding[6] = {
        false, true,
        true, false,
        false, true
    };
    int faceIdxFlipLookup[6] = {
        0, 2, 1,
        3, 5, 4
    };
    int faceIdxIdx = flipWinding[faceIdx] ? faceIdxFlipLookup[faceVertIdx] : faceVertIdx;
    int faceIdcLookup[6] = {
        0, 1, 2,
        0, 2, 3
    };
    vec2 faceUVOptions[4] = {
        vec2(0, 0),
        vec2(0, 1),
        vec2(1, 1),
        vec2(1, 0)
    };
    ivec3 faceTangentAxesOptions[3] = {
        ivec3(1, 2, 0),
        ivec3(0, 2, 1),
        ivec3(0, 1, 2)
    };
    int elementIdx = faceIdcLookup[faceIdxIdx];
    ivec3 faceTangentAxes = faceTangentAxesOptions[faceAxis];

    //Generate the fragment inputs.
    vec2 faceUV = faceUVOptions[elementIdx];
    vec3 gridPosF = vec3(gridCell);
    gridPosF[faceTangentAxes.x] += faceUV.x;
    gridPosF[faceTangentAxes.y] += faceUV.y;
    gridPosF[faceTangentAxes.z] += faceDirMask;
    #if !DEPTH_ONLY
        o_normal = vec3(0, 0, 0);
        o_normal[faceAxis] = float(faceDir);

        o_cellValue = cellValue;
        o_gridCell = gridCell;
        o_faceUV = faceUV;
        o_gridPosF = gridPosF;
    #endif

    gl_Position = u_data.matrix_view_proj * vec4(gridPosF, 1);
}


#START_FRAGMENT

in flat ivec3 o_gridCell;
in flat int o_cellValue;
in vec3 o_gridPosF;
in vec3 o_normal;
in vec2 o_faceUV;

#if !DEPTH_ONLY
    out vec4 outColor;
#endif

void main() {
#if !DEPTH_ONLY

    //Compute coordinate stuff.
    vec3 fragToCamera = u_data.cam_pos - o_gridPosF;
    float distToCamera = length(fragToCamera);
    vec3 dirTowardsCamera = fragToCamera / max(0.000001, distToCamera);

    //Compute lighting.
    UboData_Material mat = u_data.cell_materials[o_cellValue];
    if (!u_data.lighting_enabled || mat.mode == CELL_MATERIAL_LIGHT_SOURCE)
    {
        outColor = vec4(mat.albedo, mat.opacity);
    }
    else
    {
        float shadowMask = computeShadows(
            o_gridPosF,
            u_data.sun_dir,
            sampler2DShadow(u_data.sun_shadowmap),
            u_data.sun_shadowmap_mat_world_to_texel,
            u_data.shadowmap_world_bias
        );
        vec3 litColor = microfacetLighting(
            o_normal, dirTowardsCamera,
            -u_data.sun_dir, u_data.sun_color * shadowMask, u_data.ambient_light,
            mat.albedo, (mat.mode == CELL_MATERIAL_DIELECTRIC) ? 0.0 : 1.0, mat.roughness
        );
        outColor = vec4(litColor, mat.opacity);
    }


#endif
}