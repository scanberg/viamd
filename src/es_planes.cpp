#include <es_planes.h>

#include <viamd.h>
#include <representation.h>
#include <gfx/gl_utils.h>

#include <imgui.h>
#include <imgui_internal.h>
#include <implot.h>
#include <ImGuizmo.h>

#include <math.h>
#include <float.h>

void plane_frame(vec3_t* u, vec3_t* v, const Plane& p, vec3_t right) {
    const vec3_t n = p.normal;
    // 'right' without its part along n; the next of x, y when it is (nearly) along n: of x and y one is always at least
    // 45 degrees off n (n.x^2 + n.y^2 <= 1), so the last is never too short
    const vec3_t candidates[3] = {right, vec3_set(1, 0, 0), vec3_set(0, 1, 0)};
    vec3_t r = vec3_set(1, 0, 0);
    for (int i = 0; i < 3; ++i) {
        const float len = vec3_length(candidates[i]);
        if (!(len > 0.0f)) continue;
        const vec3_t c = vec3_mul1(candidates[i], 1.0f / len);
        r = vec3_sub(c, vec3_mul1(n, vec3_dot(c, n)));
        if (vec3_length(r) > 0.1f) break;
    }
    *u = vec3_normalize(r);
    *v = vec3_cross(n, *u);
}

static inline vec3_t transform_point(const mat4_t& m, vec3_t p) {
    const vec4_t r = mat4_mul_vec4(m, vec4_set(p.x, p.y, p.z, 1.0f));
    return vec3_set(r.x, r.y, r.z);
}

int plane_box_polygon(vec3_t out[6], const Plane& plane, const mat4_t& box_to_world) {
    vec3_t c[8];
    float dist[8];
    for (int i = 0; i < 8; ++i) {
        c[i] = transform_point(box_to_world, vec3_set((float)(i & 1), (float)((i >> 1) & 1), (float)((i >> 2) & 1)));
        dist[i] = plane_distance(plane, c[i]);
    }
    // the 12 edges as pairs of corners that differ in one bit
    static const int edge[12][2] = {{0,1},{2,3},{4,5},{6,7},{0,2},{1,3},{4,6},{5,7},{0,4},{1,5},{2,6},{3,7}};
    vec3_t pts[12];
    int n = 0;
    for (int e = 0; e < 12; ++e) {
        const float da = dist[edge[e][0]], db = dist[edge[e][1]];
        if ((da > 0.0f) == (db > 0.0f)) continue;
        const float t = da / (da - db);
        pts[n++] = vec3_lerp(c[edge[e][0]], c[edge[e][1]], t);
    }
    if (n < 3) return 0;
    // in order around their centre, by angle in the plane
    vec3_t centre = vec3_zero();
    for (int i = 0; i < n; ++i) centre = vec3_add(centre, pts[i]);
    centre = vec3_mul1(centre, 1.0f / (float)n);
    vec3_t u, v;
    plane_frame(&u, &v, plane, vec3_set(1, 0, 0));
    float ang[12];
    for (int i = 0; i < n; ++i) {
        const vec3_t d = vec3_sub(pts[i], centre);
        ang[i] = atan2f(vec3_dot(d, v), vec3_dot(d, u));
    }
    for (int i = 1; i < n; ++i) {
        for (int j = i; j > 0 && ang[j] < ang[j - 1]; --j) {
            const float ta = ang[j]; ang[j] = ang[j - 1]; ang[j - 1] = ta;
            const vec3_t tp = pts[j]; pts[j] = pts[j - 1]; pts[j - 1] = tp;
        }
    }
    // an edge through a corner gives the corner twice
    int m = 0;
    for (int i = 0; i < n && m < 6; ++i) {
        if (m > 0 && vec3_distance(pts[i], out[m - 1]) < 1.0e-5f) continue;
        out[m++] = pts[i];
    }
    if (m > 1 && vec3_distance(out[0], out[m - 1]) < 1.0e-5f) --m;
    return m >= 3 ? m : 0;
}

bool es_volume_box(mat4_t* box_to_world, const ElectronicStructureRepresentation& es) {
    if (!es.density_vol.tex_id) return false;
    *box_to_world = es.density_vol.texture_to_world;
    return true;
}

int es_planes_add(ApplicationState* state, ElectronicStructureRepresentation* es, uint32_t rep_id) {
    mat4_t box;
    if (es->num_planes >= ES_MAX_PLANES || !es_volume_box(&box, *es)) return -1;
    const vec3_t centre = transform_point(box, vec3_set(0.5f, 0.5f, 0.5f));
    const vec3_t forward = quat_mul_vec3(state->view.camera.orientation, vec3_set(0, 0, -1));
    const int i = es->num_planes++;
    es->planes[i].plane = plane_from_point_normal(centre, forward);
    es->planes[i].enabled = true;
    state->representation.gizmo.rep_id = rep_id;
    state->representation.gizmo.plane = i;
    return i;
}

vec4_t es_volume_plane(const ElectronicStructureRepresentation& es, int i) {
    // n . (A p + t) + d = (A^T n) . p + (n . t + d), with p in the texture's unit cube
    const mat4_t& m = es.density_vol.texture_to_world;
    const vec3_t n = es.planes[i].plane.normal;
    const vec3_t c0 = vec3_set(m.elem[0][0], m.elem[0][1], m.elem[0][2]);
    const vec3_t c1 = vec3_set(m.elem[1][0], m.elem[1][1], m.elem[1][2]);
    const vec3_t c2 = vec3_set(m.elem[2][0], m.elem[2][1], m.elem[2][2]);
    const vec3_t t  = vec3_set(m.elem[3][0], m.elem[3][1], m.elem[3][2]);
    return vec4_set(vec3_dot(c0, n), vec3_dot(c1, n), vec3_dot(c2, n), vec3_dot(n, t) + es.planes[i].plane.d);
}

int es_volume_clip_planes(vec4_t* out, const ElectronicStructureRepresentation& es) {
    int k = 0;
    for (int i = 0; i < es.num_planes && i < ES_MAX_PLANES; ++i) {
        if (es.planes[i].enabled) out[k++] = es_volume_plane(es, i);
    }
    return k;
}

uint32_t es_volume_transfer_function(ElectronicStructureRepresentation* es) {
    const auto& vol = es->volume;
    uint64_t h = md_hash64(&vol.map.colormap, sizeof(vol.map.colormap), 0);
    h = md_hash64(&vol.map.symmetric, sizeof(vol.map.symmetric), h);
    const float ramp[4] = {vol.opacity, vol.ramp_beg, vol.ramp_end, vol.ramp_exponent};
    h = md_hash64(ramp, sizeof(ramp), h);
    if (vol.tf_tex && h == vol.tf_hash) return vol.tf_tex;

    enum { RES = 256 };
    float px[RES * 4];
    const float beg = CLAMP(vol.ramp_beg, 0.0f, 1.0f);
    const float end = MAX(CLAMP(vol.ramp_end, 0.0f, 1.0f), beg + 1.0e-4f);
    for (int i = 0; i < RES; ++i) {
        const float t = (float)i / (float)(RES - 1);
        const ImVec4 col = ImPlot::SampleColormap(t, vol.map.colormap);
        // the ramp over the value, or over its magnitude on a symmetric range
        const float x = vol.map.symmetric ? fabsf(2.0f * t - 1.0f) : t;
        const float r = CLAMP((x - beg) / (end - beg), 0.0f, 1.0f);
        const float a = vol.opacity * powf(r, MAX(vol.ramp_exponent, 0.01f));
        px[4 * i + 0] = col.x;
        px[4 * i + 1] = col.y;
        px[4 * i + 2] = col.z;
        px[4 * i + 3] = CLAMP(a, 0.0f, 1.0f);
    }
    if (!es->volume.tf_tex) gl::init_texture_2D(&es->volume.tf_tex, RES, 1, GL_RGBA32F);
    gl::set_texture_2D_data(es->volume.tf_tex, 0, px, GL_RGBA32F);
    es->volume.tf_hash = h;
    return es->volume.tf_tex;
}

// ---------------------------------------------------------------------------------------------------------------
// The gizmo
// ---------------------------------------------------------------------------------------------------------------

// The frame the gizmo was given when a drag began, and what it made of it since: passed back to it while the drag
// lasts, as it turns what it is given (a frame made again from the plane each frame could turn about the normal
// under it)
static struct {
    bool     active = false;
    uint32_t rep_id = 0;
    int      plane = -1;
    mat4_t   frame = {};
} gizmo_drag;

static inline bool project(ImVec2* out, const mat4_t& world_to_clip, vec3_t p, vec2_t canvas_min, vec2_t canvas_size) {
    const vec4_t c = mat4_mul_vec4(world_to_clip, vec4_set(p.x, p.y, p.z, 1.0f));
    if (c.w <= 1.0e-6f) return false;
    const float x = c.x / c.w, y = c.y / c.w;
    out->x = canvas_min.x + (x * 0.5f + 0.5f) * canvas_size.x;
    out->y = canvas_min.y + (1.0f - (y * 0.5f + 0.5f)) * canvas_size.y;
    return true;
}

static inline float segment_distance(ImVec2 p, ImVec2 a, ImVec2 b) {
    const ImVec2 ab = b - a, ap = p - a;
    const float l2 = ab.x * ab.x + ab.y * ab.y;
    const float t = l2 > 0.0f ? CLAMP((ap.x * ab.x + ap.y * ab.y) / l2, 0.0f, 1.0f) : 0.0f;
    const ImVec2 d = ap - ab * t;
    return sqrtf(d.x * d.x + d.y * d.y);
}

static Representation* gizmo_target(ApplicationState* state, uint32_t rep_id, int plane) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        Representation& rep = state->representation.reps[i];
        if (rep.node.id != rep_id) continue;
        const ElectronicStructureRepresentation& es = rep.electronic_structure;
        if (rep.type != RepresentationType::ElectronicStructure || !representation_shown(*state, rep) ||
            es.display == ElectronicStructureDisplay::Isosurface || plane < 0 || plane >= es.num_planes) {
            return nullptr;
        }
        return &rep;
    }
    return nullptr;
}

bool es_planes_gizmo(ApplicationState* state, vec2_t canvas_min, vec2_t canvas_size) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const mat4_t& view = state->view.param.matrix.curr.view;
    const mat4_t& proj = state->view.param.matrix.curr.proj_no_jitter;
    const mat4_t world_to_clip = proj * view;
    const ImVec2 mouse = ImGui::GetMousePos();
    const bool canvas_hovered = ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);

    // Something else has the mouse (the view being turned, a slider being dragged): what started first keeps it, so the
    // planes neither take it nor pick, though the mouse passes over them. ActiveId is that of the frame before here,
    // as the view's surface is submitted after this.
    const bool mouse_taken = ImGui::GetActiveID() != 0 && !gizmo_drag.active;

    auto& sel = state->representation.gizmo;
    Representation* target = gizmo_target(state, sel.rep_id, sel.plane);
    if (!target) {
        sel.plane = -1;
        gizmo_drag.active = false;
    }

    // The outlines, and the one the mouse is on
    uint32_t hover_rep = 0;
    int hover_plane = -1;
    float hover_dist = 6.0f;
    vec3_t gizmo_centre = vec3_zero();
    bool gizmo_centre_valid = false;
    for (size_t r = 0; r < md_array_size(state->representation.reps); ++r) {
        const Representation& rep = state->representation.reps[r];
        if (rep.type != RepresentationType::ElectronicStructure || !representation_shown(*state, rep) || !rep.type_is_valid) continue;
        const ElectronicStructureRepresentation& es = rep.electronic_structure;
        if (es.display == ElectronicStructureDisplay::Isosurface) continue;
        mat4_t box;
        if (!es_volume_box(&box, es)) continue;
        const int count = es.display == ElectronicStructureDisplay::Plane ? MIN(es.num_planes, 1) : es.num_planes;
        for (int i = 0; i < count; ++i) {
            vec3_t poly[6];
            const int n = plane_box_polygon(poly, es.planes[i].plane, box);
            if (n == 0) continue;
            const bool selected = target == &rep && sel.plane == i;
            if (selected) {
                for (int k = 0; k < n; ++k) gizmo_centre = vec3_add(gizmo_centre, poly[k]);
                gizmo_centre = vec3_mul1(gizmo_centre, 1.0f / (float)n);
                gizmo_centre_valid = true;
            }
            ImVec2 sp[6];
            bool all = true;
            for (int k = 0; k < n; ++k) all &= project(&sp[k], world_to_clip, poly[k], canvas_min, canvas_size);
            if (!all) continue;
            const bool enabled = es.planes[i].enabled;
            // a light line on a dark halo: seen over a light background and a dark one alike
            const ImU32 col = selected ? IM_COL32(255, 200, 40, 255) : enabled ? IM_COL32(255, 255, 255, 220) : IM_COL32(200, 200, 200, 130);
            const ImU32 halo = enabled || selected ? IM_COL32(0, 0, 0, 110) : IM_COL32(0, 0, 0, 60);
            const float width = selected ? 2.0f : 1.5f;
            dl->AddPolyline(sp, n, halo, ImDrawFlags_Closed, width + 2.0f);
            dl->AddPolyline(sp, n, col, ImDrawFlags_Closed, width);
            for (int k = 0; k < n; ++k) {
                const float d = segment_distance(mouse, sp[k], sp[(k + 1) % n]);
                if (d < hover_dist) {
                    hover_dist = d;
                    hover_rep = rep.node.id;
                    hover_plane = i;
                }
            }
        }
    }

    // The gizmo on the picked plane, at the middle of where it cuts the box: it slides along the normal and turns about
    // the two axes in the plane
    bool capture = false;
    static bool gizmo_over_last = false;
    if (!(target && gizmo_centre_valid)) gizmo_over_last = false;
    if (target && gizmo_centre_valid) {
        ElectronicStructureRepresentation& es = target->electronic_structure;
        Plane& plane = es.planes[sel.plane].plane;
        mat4_t frame;
        if (gizmo_drag.active && gizmo_drag.rep_id == sel.rep_id && gizmo_drag.plane == sel.plane) {
            frame = gizmo_drag.frame;
        } else {
            // its rings along the screen: the first axis in the plane is the camera's right as near as the plane allows
            vec3_t u, v;
            plane_frame(&u, &v, plane, quat_mul_vec3(state->view.camera.orientation, vec3_set(1, 0, 0)));
            const vec3_t n = plane.normal;
            const vec3_t o = plane_project(plane, gizmo_centre);
            frame = mat4_ident();
            frame.elem[0][0] = u.x; frame.elem[0][1] = u.y; frame.elem[0][2] = u.z;
            frame.elem[1][0] = v.x; frame.elem[1][1] = v.y; frame.elem[1][2] = v.z;
            frame.elem[2][0] = n.x; frame.elem[2][1] = n.y; frame.elem[2][2] = n.z;
            frame.elem[3][0] = o.x; frame.elem[3][1] = o.y; frame.elem[3][2] = o.z;
        }
        ImGuizmo::SetOrthographic(false);
        ImGuizmo::SetDrawlist(dl);
        ImGuizmo::SetRect(canvas_min.x, canvas_min.y, canvas_size.x, canvas_size.y);
        ImGuizmo::PushID((int)(sel.rep_id * ES_MAX_PLANES + (uint32_t)sel.plane));
        const ImGuizmo::OPERATION op = (ImGuizmo::OPERATION)(ImGuizmo::TRANSLATE_Z | ImGuizmo::ROTATE_X | ImGuizmo::ROTATE_Y);
        // The gizmo takes a click only when no item was hovered the frame before (ImGuizmo's CanActivate), and the view's
        // interaction surface, disabled while the gizmo has the mouse, is still hovered (ImGui's disabled items keep
        // HoveredId). So once the gizmo had the mouse last frame, and nothing is active, that hover is the surface's
        // and is let go of.
        if (gizmo_over_last && canvas_hovered && ImGui::GetActiveID() == 0) {
            GImGui->HoveredIdPreviousFrame = 0;
        }
        const bool changed = ImGuizmo::Manipulate(&view.elem[0][0], &proj.elem[0][0], op, ImGuizmo::LOCAL, &frame.elem[0][0]);
        const bool using_it = ImGuizmo::IsUsing();
        capture = using_it || (ImGuizmo::IsOver() && !mouse_taken);
        gizmo_over_last = capture;
        ImGuizmo::PopID();
        if (changed) {
            const vec3_t n = vec3_set(frame.elem[2][0], frame.elem[2][1], frame.elem[2][2]);
            const vec3_t o = vec3_set(frame.elem[3][0], frame.elem[3][1], frame.elem[3][2]);
            if (vec3_length(n) > 1.0e-6f) plane = plane_from_point_normal(o, n);
        }
        gizmo_drag.active = using_it;
        gizmo_drag.rep_id = sel.rep_id;
        gizmo_drag.plane = sel.plane;
        gizmo_drag.frame = frame;
    }

    // A click on an outline picks its plane
    if (!capture && !mouse_taken && canvas_hovered && hover_plane >= 0) {
        ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            sel.rep_id = hover_rep;
            sel.plane = hover_plane;
            gizmo_drag.active = false;
            capture = true;
        }
    }
    return capture;
}
