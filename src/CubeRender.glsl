//NOTE: this file is evaluated as the inside of a Julia string,
//    so we can reference constants in the codebase using Julia string interpolation.

#define MISSING_CELL $CELL_CODE_INVALID
#define N_CELL_TYPES $(length(CELL_TYPES))

uniform mat4 u_transformWorld, u_transformViewProj;
uniform mat3 u_transformWorldNormals;
uniform usampler3D u_gridTex;
uniform CellLookupUBO {
    vec3 u_colors[N_CELL_TYPES];
};


#START_VERTEX

out ivec3 gIn_voxelIdx;

void main() {
    //Convert from primitive index to voxel grid cell.
    ivec3 nVoxels = textureSize(u_gridTex, 0);
    gIn_voxelIdx = ivec3(
        gl_VertexID % nVoxels.x,
        (gl_VertexID / u_nVoxels.x) % nVoxels.y,
        gl_VertexID / (nVoxels.x * nVoxels.y)
    );
}


#START_GEOMETRY

layout(points) in;
in ivec3 gIn_voxelIdx[];

layout(triangle_strip, max_vertices=24) out;
out vec2 fIn_uv;
out vec3 fIn_localPos;
out vec3 fIn_worldPos;
flat out vec3 fIn_localNormal;
flat out vec3 fIn_worldNormal;
flat out vec3 fIn_color;

void main() {
    ivec3 voxelIdx = gIn_voxelIdx[0];

    uint colorCode = texelFetch(u_gridTex, voxelIdx, 0).r;
    if (colorCode == MISSING_CELL)
        return;

    vec3 color = u_colors[colorCode];
    ivec3 nVoxels = textureSize(u_gridTex, 0);

    //Generate each face.
    for (int axis = 0; axis < 3; ++axis)
    {
        ivec2 otherAxesChoices[3] = {
            ivec2(1, 2),
            ivec2(0, 2),
            ivec2(0, 1)
        };
        ivec2 otherAxes = otherAxesChoices[axis];

        for (int bDir = 0; bDir < 2; ++bDir)
        {
            int dir = (bDir * 2) - 1;

            //Get the neighbor on this face.
            ivec3 neighborPos = voxelIdx;
            neighborPos[axis] += dir;
            //If the neighbor is past the edge of the voxel grid, assume it's empty space.
            uint neighborVoxel = MISSING_CELL;
            if (neighborPos[axis] >= 0 && neighborPos[axis] < nVoxels[axis])
                neighborVoxel = texelFetch(u_gridTex, neighborPos, 0).r;

            //If the neighbor is empty, then this face of our voxel should be rendered.
            if (neighborVoxel == MISSING_CELL)
            {
                //Compute the 4 corners of this face,
                //    and emit a triangle strip for them.
                vec3 minLocalCorner = vec3(voxelIdx);
                const vec2 cornerFaceOffsets[4] = {
                    vec2(0, 0),
                    vec2(1, 0),
                    vec2(0, 1),
                    vec2(1, 1)
                };
                for (int cornerI = 0; cornerI < 4; ++cornerI)
                {
                    vec2 uv = cornerFaceOffsets[cornerI];

                    vec3 localCorner = minLocalCorner;
                    localCorner[axis] += bDir;
                    localCorner[otherAxes.x] += uv.x;
                    localCorner[otherAxes.y] += uv.y;

                    vec3 localNormal = vec3(0, 0, 0);
                    localNormal[axis] = dir;

                    vec4 worldCorner4 = u_transformWorld * vec4(localCorner, 1);
                    vec3 worldNormal = normalize(u_transformWorldNormals * localNormal);

                    fIn_uv = uv;
                    fIn_localPos = localCorner;
                    fIn_worldPos = worldCorner4.xyz / worldCorner4.w;
                    fIn_localNormal = localNormal;
                    fIn_worldNormal = worldNormal;
                    fIn_color = color;
                    gl_Position = u_transformViewProj * fIn_worldPos;
                    EmitVertex();
                }
                EndPrimitive();
            }
        }
    }
}


#START_FRAGMENT

in vec2 fIn_uv;
in vec3 fIn_localPos;
in vec3 fIn_worldPos;
flat in vec3 fIn_localNormal;
flat in vec3 fIn_worldNormal;
flat in vec3 fIn_color;

out vec4 fOut_color;

void main() {
    //TODO: Add more interesting lighting.
    float diffuse = max(0.2, dot(fIn_worldNormal, normalize(vec3(1, 1, -1))));
    fIn_color = fOut_color * diffuse;
}