// Shared-ray overlays: extinction is combined in each slab, surface events are sorted
// by ray distance, not by the user's layer list. Single volumes retain their exact solver.
static constexpr str_t overlay_fragment_source = STR_LIT(R"(
#version 410 core
#define LAYERS 6
#define EVENTS 48
uniform int u_count;
uniform sampler3D u_volume[LAYERS];
uniform sampler2D u_tf[LAYERS];
uniform sampler2D u_depth;
uniform mat4 u_clip_to_model[LAYERS];
uniform mat3 u_normal_to_view[LAYERS];
uniform mat4 u_clip_to_view;
uniform vec3 u_min[LAYERS], u_max[LAYERS];
uniform vec2 u_range[LAYERS];
uniform int u_dvr[LAYERS], u_iso_count[LAYERS];
uniform float u_values[EVENTS];
uniform vec4 u_colors[EVENTS];
uniform vec2 u_inv_res;
uniform vec3 u_env;
layout(location=0) out vec4 out_color;

vec3 unproject(mat4 m, vec2 xy, float z) {
    vec4 p = m * vec4(xy, z * 2.0 - 1.0, 1.0);
    return p.xyz / p.w;
}
vec3 linear_rgb(vec3 c) {
    return mix(c / 12.92, pow((c + 0.055) / 1.055, vec3(2.4)), greaterThan(c, vec3(0.04045)));
}
// Same plastic Cook-Torrance surface model as the single-volume and deferred passes.
vec4 shade(vec4 color, vec3 N, vec3 V) {
    float nv = max(dot(N,V), 0.0);
    vec3 light = normalize(vec3(1.0));
    float nl = max(dot(N,light), 0.0);
    vec3 H = normalize(light + V);
    float nh = max(dot(N,H), 0.0), hv = max(dot(H,V), 0.0);
    float a2 = pow(0.3,4.0);
    float d = nh*nh*(a2-1.0)+1.0;
    float D = a2/(3.14159265*d*d);
    float k = 1.3*1.3/8.0;
    float G = nv/(nv*(1.0-k)+k) * nl/(nl*(1.0-k)+k);
    float F = 0.04 + 0.96*pow(1.0-hv,5.0);
    float Fe = 0.04 + 0.66*pow(1.0-nv,5.0);
    float Fv = 0.04 + 0.96*pow(1.0-nv,5.0);
    vec3 diffuse = color.rgb / 3.14159265 * ((1.0-F)*10.0*nl + (1.0-Fe)*u_env);
    vec3 specular = vec3(D*G*F/(4.0*nv*nl+0.0001)*10.0*nl) + Fe*u_env;
    return vec4(color.a*diffuse + specular, 1.0-(1.0-color.a)*(1.0-Fv));
}
void slab(inout vec3 L, inout float T, vec3 emission, float extinction, float length) {
    if (extinction <= 0.0 || length <= 0.0) return;
    float a = 1.0-exp(-extinction*length);
    L += T*a*emission/extinction;
    T *= 1.0-a;
}
void main() {
    vec2 xy = gl_FragCoord.xy*u_inv_res*2.0-1.0;
    float depth = texelFetch(u_depth, ivec2(gl_FragCoord.xy), 0).r;
    vec3 pn[LAYERS], seg[LAYERS];
    vec2 spans[LAYERS];
    float first = 1.0, last = 0.0, dt = 1.0;
    for (int j=0;j<u_count;++j) {
        pn[j] = unproject(u_clip_to_model[j],xy,0.0);
        seg[j] = unproject(u_clip_to_model[j],xy,depth)-pn[j];
        vec3 inv = 1.0/mix(seg[j],vec3(1e-20),equal(seg[j],vec3(0.0)));
        vec3 a = (u_min[j]-pn[j])*inv, b = (u_max[j]-pn[j])*inv;
        vec3 lo=min(a,b), hi=max(a,b);
        float beg=max(max(lo.x,lo.y),max(lo.z,0.0));
        float end=min(min(hi.x,hi.y),min(hi.z,1.0));
        for (int axis=0;axis<3;++axis)
            if (abs(seg[j][axis])<1e-20 && (pn[j][axis]<u_min[j][axis] || pn[j][axis]>u_max[j][axis])) end=beg;
        spans[j]=vec2(beg,end);
        if (beg<end) {
            first=min(first,beg); last=max(last,end);
            dt=min(dt,0.5/max(length(seg[j]*vec3(textureSize(u_volume[j],0))),1e-6));
        }
    }
    if (first>=last) discard;
    vec3 V = normalize(unproject(u_clip_to_view,xy,0.0)-unproject(u_clip_to_view,xy,depth));
    vec3 L=vec3(0.0); float T=1.0;
    float t=first;
    while (t<last && T>0.005) {
        float end=min(t+dt,last);
        float event_t[EVENTS]; vec4 event_c[EVENTS]; int count=0;
        vec3 emission=vec3(0.0); float extinction=0.0;
        bool occupied=false; float next_entry=last;
        for (int j=0;j<u_count;++j) {
            if (spans[j].x>t) next_entry=min(next_entry,spans[j].x);
            float a=max(t,spans[j].x), b=min(end,spans[j].y);
            if (a>=b) continue;
            occupied=true;
            if (u_dvr[j]!=0) {
                vec3 p=pn[j]+seg[j]*(a+b)*0.5;
                float density=texture(u_volume[j],p).r;
                float mapped=clamp((density-u_range[j].x)*u_range[j].y,0.0,1.0);
                vec4 tf=texture(u_tf[j],vec2(mapped,0.5));
                float tau=-log(max(1.0-clamp(tf.a,0.0,1.0),1e-6))*150.0*length(seg[j])*(b-a)/(end-t);
                extinction+=tau; emission+=tau*linear_rgb(tf.rgb);
            } else {
                float va=texture(u_volume[j],pn[j]+seg[j]*a).r;
                float vb=texture(u_volume[j],pn[j]+seg[j]*b).r;
                for (int k=0;k<u_iso_count[j];++k) {
                    int idx=j*8+k;
                    float value=u_values[idx];
                    if ((va<value)==(vb<value)) continue;
                    float lo=a,hi=b;
                    for (int refine=0;refine<8;++refine) {
                        float mid=(lo+hi)*0.5;
                        float vm=texture(u_volume[j],pn[j]+seg[j]*mid).r;
                        if ((vm<value)==(va<value)) lo=mid; else hi=mid;
                    }
                    float crossing=(lo+hi)*0.5;
                    vec3 p=pn[j]+seg[j]*crossing;
                    vec3 eps=1.0/vec3(textureSize(u_volume[j],0));
                    vec3 grad=vec3(
                        texture(u_volume[j],p+vec3(eps.x,0,0)).r-texture(u_volume[j],p-vec3(eps.x,0,0)).r,
                        texture(u_volume[j],p+vec3(0,eps.y,0)).r-texture(u_volume[j],p-vec3(0,eps.y,0)).r,
                        texture(u_volume[j],p+vec3(0,0,eps.z)).r-texture(u_volume[j],p-vec3(0,0,eps.z)).r)/eps;
                    vec3 N=u_normal_to_view[j]*grad;
                    N=length(N)>1e-10 ? normalize(N) : V;
                    if (dot(N,V)<0.0) N=-N;
                    vec4 c=shade(u_colors[idx],N,V);
                    int pos=count++;
                    while (pos>0 && event_t[pos-1]>crossing) {
                        event_t[pos]=event_t[pos-1]; event_c[pos]=event_c[pos-1]; --pos;
                    }
                    event_t[pos]=crossing; event_c[pos]=c;
                }
            }
        }
        float prev=t;
        for (int e=0;e<count;++e) {
            slab(L,T,emission,extinction,event_t[e]-prev);
            L+=T*event_c[e].rgb; T*=1.0-event_c[e].a;
            prev=event_t[e];
        }
        slab(L,T,emission,extinction,end-prev);
        t=occupied ? end : max(end,next_entry);
    }
    out_color=vec4(L,1.0-T);
}
)");
