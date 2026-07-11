//BRDF-related equations, using the "micro-facet" model.
//Reference: https://learnopengl.com/PBR/Lighting

//Approximates the light reflected from a surface, given its glancing angle.
vec3 fresnelSchlick(float diffuseStrength, vec3 F0) {
    return F0 + ((1.0 - F0) * pow(1.0 - diffuseStrength, 5.0));
}
//Approximates the proportion of micro-facets
//    which are facing the right way to reflect light into the camera.
float distributionGGX(float specularStrength, float roughness) {
    float a   = roughness*roughness,
          a2  = a*a;

    float num   = a2;
    float denom = (specularStrength * specularStrength * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;

    return num / denom;
}
//Approximates the proportion of micro-facets which are visible
//    to both the light and the camera.
float geometrySchlickGGX(float diffuseStrength, float roughness) {
    float num   = diffuseStrength;
    float denom = diffuseStrength * (1.0 - roughness) + roughness;

    return diffuseStrength / ((diffuseStrength * (1.0 - roughness)) + roughness);
}
float geometrySmith(float diffuseNormalAndCamera,
                    float diffuseNormalAndLight,
                    float roughness) {
    float r = (roughness + 1.0);
    float k = (r*r) / 8.0;
    return geometrySchlickGGX(diffuseNormalAndCamera, k) *
            geometrySchlickGGX(diffuseNormalAndLight, k);
}

//Implements a microfacet lighting model, using approximations for various factors
//    (see the functions above).
vec3 microfacetLighting(vec3 normal, vec3 towardsCameraN, vec3 towardsLightN,
                        vec3 lightIrradiance,
                        vec3 albedo, float metallic, float roughness) {
    vec3 idealNormal = normalize(towardsLightN + towardsCameraN);

    float diffuseStrength = SATURATE(dot(normal, towardsLightN)),
          specularStrength = SATURATE(dot(idealNormal, normal)),
          normalClosenessToCamera = SATURATE(dot(normal, towardsCameraN));

    vec3 F0 = mix(vec3(0.04), albedo, metallic),
         F = fresnelSchlick(SATURATE(dot(idealNormal, towardsCameraN)), F0);

    vec3 energyOfReflection = F,
         energyOfDiffuse = (1.0 - energyOfReflection) * (1.0 - metallic);

    float NDF = distributionGGX(specularStrength, roughness),
          G = geometrySmith(normalClosenessToCamera, diffuseStrength, roughness);

    vec3 specular = F * (NDF * G / max(0.0001, 4.0 * normalClosenessToCamera * diffuseStrength));
    vec3 totalLight = (((energyOfDiffuse / PI) * albedo) + specular) *
                      lightIrradiance * diffuseStrength;

    return totalLight;
}

//Computes the amount of global height-fog between the camera and the fragment.
vec3 computeFoggedColor(float camHeight,
                        float fragWorldHeight,
                        float fragDist3D, float fragDistVertical,
                        vec3 surfaceColor,
                        float fogHeightOffset, float fogHeightScale,
                        float fogDensity, float fogDropoff,
                        vec3 fogColor) {
    //Height-fog density is only a function of vertical position.
    //As long as that function can be integrated analytically,
    //    then total fog density itself can be integrated analytically.

    //Using a function of 'f(z) = exp(z)', the integral is 'if(z) = exp(z) + C'.
    //The C cancels out in the definite integral, which is "if(z1) - if(z2)".

    float fogStartHeight = fogHeightScale * (camHeight - fogHeightOffset),
          fogEndHeight = fogHeightScale * (fragWorldHeight - fogHeightOffset),
          fogIntegralScale = (fragDist3D / max(0.00001, fragDistVertical)),
          fogDensityIntegral = abs(fogIntegralScale * (exp(-fogEndHeight) - exp(-fogStartHeight)));

    float fogThickness = SATURATE(fogDensity * fogDensityIntegral);
    fogThickness = pow(fogThickness, fogDropoff);

    return mix(surfaceColor, fogColor, fogThickness);
}

//Computes ambient lighting.
vec3 computeAmbient(vec3 surfacePos, vec3 normal, vec3 albedo)
{
    return mix(vec3(0.03), albedo, 0.05);
}

//Computes shadow-maps.
//Returns a 0-1 mask (0 is total shadow, 1 is fully-lit).
float computeShadows(vec3 worldPos, vec3 worldLightDir,
                     sampler2DShadow shadowmap, mat4 worldToShadowTexel,
                     float shadowWorldBias) {
    worldPos += (worldLightDir * shadowWorldBias);

    vec4 texel4 = worldToShadowTexel * vec4(worldPos, 1);
    vec3 texel = texel4.xyz;
    texel /= texel4.w;
    //texel.y = 1.0 - texel.y;
    texel.z = 0.5 + (0.5 * texel.z);

    //The shadowmap bounds exactly cover the level bounds.
    //Therefore, if the position is outside the view of the shadow-map then it is not in shadow.
    if (any(lessThan(texel.xy, vec2(0))) || any(greaterThan(texel.xy, vec2(1))))
        return 1.0;
    float shadowMask = textureLod(shadowmap, texel, 0);
    return shadowMask;
}

//Computes sky color.
vec3 computeSkyColor(vec3 dirToSky, vec3 dirToSun)
{
    vec3 atmosphereColor = vec3(0.8, 0.825, 1.0);
    float sunSharpness = 256.0;

    float sunCloseness = max(0.0, dot(dirToSky, dirToSun));
    return atmosphereColor + pow(sunCloseness, sunSharpness);
}