#define PI (3.1415926535897932384626433832795)
#define PI2 (PI * 2.0)

#define OSCILLATE(a, b, inp) (mix((a), (b), 0.5 + (0.5 * sin(PI2 * (inp)))))

#define INV_LERP(a, b, x) (((x)-(a)) / ((b)-(a)))
#define SATURATE(x) clamp(x, 0.0, 1.0)
#define SHARPEN(t) smoothstep(0.0, 1.0, (t))
#define SHARPENER(t) SMOOTHERSTEP(t)

//Converts a 0-1 value to the (-1, +1) range.
#define NORMALIZED_TO_SIGNED(uv) (-1.0 + (2.0 * (uv)))
//Converts a (-1, +1) range value to the 0-1 range.
#define SIGNED_TO_NORMALIZED(normal) (0.5 + (0.5 * (normal)))

//To prevent a comma from being noticed by a macro invocation.
#define COMMA ,

#define RAND_IN_ARRAY(array, t) (array[int(mix(0.0, float(array.length()) - 0.00001, (t)))])

//A higher-quality smoothstep(), with a zero second-derivative at the edges.
#define SMOOTHERSTEP(t) clamp(t * t * t * (t * (t*6.0 - 15.0) + 10.0), \
                            0.0, 1.0)