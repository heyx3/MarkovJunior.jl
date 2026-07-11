#START_VERTEX

//Code will paste in a basic vertex-shader here.
@blit_vs
#line 6

#START_FRAGMENT

in vec2 vOut_uv;
out vec4 fOut_color;

uniform float u_time;
uniform sampler2D u_tex;

void main() {
    float rawDepth = textureLod(u_tex, vOut_uv, 0).r;
    fOut_color = vec4(vec3(
        pow(rawDepth, 1.0)
    ), 1);
}