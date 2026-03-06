#version 330

#pragma optionNV(unroll all)

#ifndef AO_PERSPECTIVE
#define AO_PERSPECTIVE 1
#endif

#define PI 3.14159265359

// UBO uniform
struct GTAOData {
    vec4  projInfo;
    vec2  invRes;
    float projScl;
    float zMax;
};

uniform sampler2D uDepth;
uniform sampler2D uNormal;
uniform sampler2D uRandom;

layout(location = 0) out vec4 outColor;

layout(std140) uniform uUbo {
    GTAOData ubo;
};

vec2 rand2(vec2 n) {
    return fract(sin(dot(n.xy, vec2(12.9898, 78.233)))* vec2(43758.5453, 28001.8384));
}

vec3 UvToView(vec2 uv, float eye_z) {
#if AO_PERSPECTIVE
    return vec3((uv * ubo.projInfo.xy + ubo.projInfo.zw) * eye_z, eye_z);
#else
    return vec3((uv * ubo.projInfo.xy + ubo.projInfo.zw), eye_z);
#endif
}

vec3 reconstructViewPos(vec2 uv, float depth) {
    return UvToView(uv, depth);
}

vec3 decodeViewNormal(vec2 enc) {
    vec2 fenc = enc*4-2;
    float f = dot(fenc,fenc);
    float g = sqrt(1-f/4.0);
    vec3  n = vec3(fenc*g, 1-f/2.0);
    return n * vec3(1,1,-1);
}

float pixelFootprint(float viewZ) {
    return viewZ * ubo.projInfo.x * ubo.invRes.y;
}

float logStepDistance(int i, int steps, float Rmin, float Rmax) {
    float t = float(i) / float(steps - 1); return Rmin * pow(Rmax / Rmin, t);
}

float integrateSliceGTAO_Log( vec2 uv, vec3 viewPos, vec3 viewNormal, vec2 dir, float RminPixels, float RmaxPixels, int steps) {
    float pixelSize = pixelFootprint(viewPos.z);
    float RminWorld = RminPixels * pixelSize;
    float RmaxWorld = RmaxPixels * pixelSize;

    vec2 jitter = rand2(gl_FragCoord.xy);

    // View direction (point -> camera)
    vec3 V = normalize(-viewPos);

    // Build a stable slice basis:
    // Convert screen-space direction 'dir' into a view-space direction using derivatives.
    vec3 dPdx = dFdx(viewPos);
    vec3 dPdy = dFdy(viewPos);

    vec3 T = dir.x * dPdx + dir.y * dPdy;
    float t2 = dot(T, T);
    if (t2 < 1e-10) {
        // Fallback (should be rare)
        T = vec3(dir, 0.0);
    } else {
        T *= inversesqrt(t2);
    }

    // Orthonormalize T against V (important for stability)
    T = T - V * dot(T, V);
    float t2o = dot(T, T);
    if (t2o < 1e-10) {
        // Another fallback: pick any vector perpendicular to V
        vec3 a = (abs(V.z) < 0.999) ? vec3(0,0,1) : vec3(0,1,0);
        T = normalize(cross(a, V));
    } else {
        T *= inversesqrt(t2o);
    }

    vec3 B = cross(T, V);
    float b2 = dot(B, B);
    B *= inversesqrt(max(b2, 1e-10));

    // Project normal into slice plane
    vec3 Nproj = viewNormal - B * dot(viewNormal, B);
    float n2 = dot(Nproj, Nproj);
    Nproj *= inversesqrt(max(n2, 1e-10));

    // Signed elevation of normal (same angle space as horizons)
    float N_up   = dot(Nproj, V);                 // [-1, +1]
    float N_side = length(cross(Nproj, V));       // [0,  1]
    N_side = max(N_side, 1e-6);                   // avoid atan(?, 0) sensitivity
    float thetaN = atan(N_up, N_side);            // [-pi/2, +pi/2]

    float maxHorizonPos = -PI * 0.5;
    float maxHorizonNeg = -PI * 0.5;

    const float lodBias = -4.7;

    for (int i = 0; i < steps; ++i) {
        float distWorld = logStepDistance(i, steps, RminWorld, RmaxWorld);
        float distPixels = distWorld / pixelSize;
        distPixels *= jitter.y * 0.5 + 0.5;

        // Compute LOD from projected pixel radius
        float lod = clamp(log2(distPixels) + lodBias, 0.0, 0.0);
        vec2 uvOffset = dir * (distPixels * ubo.invRes);

        // + direction
        vec2 uvPos = uv + uvOffset;
        float depthPos = textureLod(uDepth, uvPos, lod).r;
        vec3 posPos = reconstructViewPos(uvPos, depthPos);
        vec3 deltaPos = posPos - viewPos;

        float dvPos = dot(deltaPos, V);
        if (dvPos > 0.0) {
            float perpPos = length(deltaPos - V * dvPos);
            float elevPos = atan(dvPos, perpPos);
            maxHorizonPos = max(maxHorizonPos, elevPos);
        }

        // - direction
        vec2 uvNeg = uv - uvOffset;
        float depthNeg = textureLod(uDepth, uvNeg, lod).r;
        vec3 posNeg = reconstructViewPos(uvNeg, depthNeg);
        vec3 deltaNeg = posNeg - viewPos;

        float dvNeg = dot(deltaNeg, V);
        if (dvNeg > 0.0) {
            float perpNeg = length(deltaNeg - V * dvNeg);
            float elevNeg = atan(dvNeg, perpNeg);
            maxHorizonNeg = max(maxHorizonNeg, elevNeg);
        }

        // Early termination
        if (maxHorizonPos >= thetaN && maxHorizonNeg >= thetaN) break;
    }

    // Clamp
    maxHorizonPos = clamp(maxHorizonPos, -PI*0.5, PI*0.5);
    maxHorizonNeg = clamp(maxHorizonNeg, -PI*0.5, PI*0.5);

    // Cosine-weighted analytic visibility
    float visibilityPos = sin(thetaN) - sin(maxHorizonPos);
    float visibilityNeg = sin(thetaN) - sin(maxHorizonNeg);
    float sliceAO = 0.5 * (visibilityPos + visibilityNeg);
    return clamp(sliceAO, 0.0, 1.0);
}

void main() {
    const float RminPixels = 3.0;
    const float RmaxPixels = 128.0;
    const int   num_directions = 4;
    const int   num_steps = 4;

    vec2 uv = gl_FragCoord.xy * ubo.invRes;
    ivec2 randomUv = ivec2(gl_FragCoord.xy) & 3;

    float depth  = textureLod(uDepth,  uv, 0).r;
    vec2  normal = textureLod(uNormal, uv, 0).rg;
    float random = texelFetch(uRandom, randomUv, 0).r;

    if (depth >= ubo.zMax) {
        outColor = vec4(1.0); // No AO for far plane
        return;
    }

    vec3 viewPos = reconstructViewPos(uv, depth);
    vec3 viewNormal = decodeViewNormal(normal);

    float angleOffset = random * 2.0 * PI;

    float ao = 0.0;
    for (int d = 0; d < num_directions; ++d) {
        float angle = angleOffset + (2.0 * PI * float(d) / float(num_directions));
        vec2 dir = vec2(cos(angle), sin(angle));
        ao += integrateSliceGTAO_Log( uv, viewPos, viewNormal, dir, RminPixels, RmaxPixels, num_steps);
    }
    
    ao /= float(num_directions);

    // Output AO as grayscale
    outColor = vec4(vec3(ao), 1.0);
}