#include <es_planes.h>

#include <viamd.h>
#include <representation.h>
#include <gfx/gl_utils.h>

#include <imgui.h>
#include <imgui_internal.h>
#include <implot.h>

#include <math.h>
#include <float.h>

static inline vec3_t transform_point(const mat4_t& m, vec3_t p) {
    const vec4_t r = mat4_mul_vec4(m, vec4_set(p.x, p.y, p.z, 1.0f));
    return vec3_set(r.x, r.y, r.z);
}


bool es_volume_box(mat4_t* box_to_world, const ElectronicStructureRepresentation& es) {
    if (!es.has_box) return false;
    *box_to_world = es.box_to_world;
    return true;
}

bool es_plane_image_grid(md_grid_t* grid, const Plane& plane_in, const mat4_t& box, float spacing, int max_dim) {
    // Which way the normal points is of no matter to the image (a slice's is turned to the camera as it is picked):
    // one of the two, so that turning it does not make the image again
    Plane plane = plane_in;
    if (vec3_dot(plane.normal, vec3_set(0.5773f, 0.6172f, 0.5345f)) < 0.0f) plane_flip(&plane);
    vec3_t poly[6];
    const int n = plane_box_polygon(poly, plane, box);
    if (n < 3 || !(spacing > 0.0f) || max_dim < 1) return false;
    // A frame of the plane by the box's first axis, not the view: the image stays put as the view turns
    const vec3_t ax = vec3_set(box.elem[0][0], box.elem[0][1], box.elem[0][2]);
    vec3_t u, v;
    plane_frame(&u, &v, plane, ax);
    float s0 = FLT_MAX, s1 = -FLT_MAX, t0 = FLT_MAX, t1 = -FLT_MAX;
    for (int k = 0; k < n; ++k) {
        const float s = vec3_dot(poly[k], u), t = vec3_dot(poly[k], v);
        s0 = MIN(s0, s); s1 = MAX(s1, s);
        t0 = MIN(t0, t); t1 = MAX(t1, t);
    }
    const float ext = MAX(s1 - s0, t1 - t0);
    if (ext > spacing * (float)max_dim) spacing = ext / (float)max_dim;
    const int w = CLAMP((int)ceilf((s1 - s0) / spacing), 1, max_dim);
    const int h = CLAMP((int)ceilf((t1 - t0) / spacing), 1, max_dim);
    // Pixel centres on the plane: the one layer's centre (the kernels sample at voxel centres) a half spacing along n
    // from the corner
    const vec3_t on = vec3_mul1(plane.normal, -plane.d);   // the plane's point nearest the origin
    const vec3_t corner = vec3_add(vec3_add(vec3_add(on, vec3_mul1(u, s0)), vec3_mul1(v, t0)), vec3_mul1(plane.normal, -0.5f * spacing));
    const float to_bohr = (float)(1.0 / BOHR_TO_ANGSTROM);
    *grid = md_grid_t{};
    grid->orientation = mat3_t{};
    for (int k = 0; k < 3; ++k) {
        grid->orientation.elem[0][k] = u.elem[k];
        grid->orientation.elem[1][k] = v.elem[k];
        grid->orientation.elem[2][k] = plane.normal.elem[k];
    }
    grid->origin  = vec3_mul1(corner, to_bohr);
    grid->spacing = vec3_set1(spacing * to_bohr);
    grid->dim[0] = w;
    grid->dim[1] = h;
    grid->dim[2] = 1;
    return true;
}

static bool faces_camera(const ApplicationState* state, const Plane& p);

int es_planes_add(ApplicationState* state, ElectronicStructureRepresentation* es, uint32_t rep_id) {
    mat4_t box;
    ElectronicStructurePlaneSet& set = electronic_structure_planes(*es);
    if (set.count >= ES_MAX_PLANES || !es_volume_box(&box, *es)) return -1;
    const vec3_t centre = transform_point(box, vec3_set(0.5f, 0.5f, 0.5f));
    const vec3_t forward = quat_mul_vec3(state->view.camera.orientation, vec3_set(0, 0, -1));
    const int i = set.count++;
    // a cut facing away from the camera keeps the half behind it; a slice faces the camera
    set.items[i] = ElectronicStructurePlane{};
    set.items[i].plane = plane_from_point_normal(centre, es->display == ElectronicStructureDisplay::Planes ? vec3_mul1(forward, -1.0f) : forward);
    if (es->display == ElectronicStructureDisplay::Planes && !faces_camera(state, set.items[i].plane)) plane_flip(&set.items[i].plane);
    state->representation.gizmo.rep_id = rep_id;
    state->representation.gizmo.plane = i;
    return i;
}

Plane es_plane_fit_selection(const ApplicationState* state, const Plane& plane, PlaneFitInfo* info) {
    const md_system_state_t& st = state->mold.state;
    const size_t count = md_bitmask_count(&state->selection.selection_mask);
    if (count == 0 || !st.xyz) {
        if (info) *info = PlaneFitInfo{};
        return plane;
    }
    md_temp_scope_t temp = md_temp_begin();
    vec3_t* pts = (vec3_t*)md_temp_alloc(temp, sizeof(vec3_t) * count);
    size_t n = 0;
    md_bitset_iter_t it = md_bitmask_iter(&state->selection.selection_mask);
    while (md_bitset_iter_next(&it) && n < count) {
        if (it.idx < st.num_atoms) pts[n++] = st.xyz[it.idx];
    }
    if (n > 1 && md_unitcell_flags(&st.unitcell) != 0) {
        // across the boundary: each the image nearest to the first
        for (size_t i = 1; i < n; ++i) pts[i] = vec3_sub(pts[i], pts[0]);
        md_util_min_image_vec3(pts + 1, n - 1, &st.unitcell);
        for (size_t i = 1; i < n; ++i) pts[i] = vec3_add(pts[i], pts[0]);
    }
    const Plane fitted = plane_fit(plane, pts, n, info);
    md_temp_end(temp);
    return fitted;
}

vec4_t es_volume_plane(const ElectronicStructureRepresentation& es, const Plane& plane) {
    // n . (A p + t) + d = (A^T n) . p + (n . t + d), with p in the texture's unit cube
    const mat4_t& m = es.box_to_world;
    const vec3_t n = plane.normal;
    const vec3_t c0 = vec3_set(m.elem[0][0], m.elem[0][1], m.elem[0][2]);
    const vec3_t c1 = vec3_set(m.elem[1][0], m.elem[1][1], m.elem[1][2]);
    const vec3_t c2 = vec3_set(m.elem[2][0], m.elem[2][1], m.elem[2][2]);
    const vec3_t t  = vec3_set(m.elem[3][0], m.elem[3][1], m.elem[3][2]);
    return vec4_set(vec3_dot(c0, n), vec3_dot(c1, n), vec3_dot(c2, n), vec3_dot(n, t) + plane.d);
}

int es_volume_clip_planes(vec4_t* out, const ElectronicStructureRepresentation& es) {
    int k = 0;
    for (int i = 0; i < es.cuts.count && i < ES_MAX_PLANES; ++i) {
        if (es.cuts.items[i].enabled) out[k++] = es_volume_plane(es, es.cuts.items[i].plane);
    }
    return k;
}

int es_volume_slices(vec4_t* out, float* out_opacity, const ElectronicStructureRepresentation& es, int* out_index) {
    int k = 0;
    for (int i = 0; i < es.slices.count && i < ES_MAX_PLANES; ++i) {
        const ElectronicStructurePlane& p = es.slices.items[i];
        if (!p.enabled || !(p.opacity > 0.0f)) continue;
        out_opacity[k] = p.opacity > 1.0f ? 1.0f : p.opacity;
        if (out_index) out_index[k] = i;
        out[k++] = es_volume_plane(es, p.plane);
    }
    return k;
}

uint32_t es_volume_transfer_function(ElectronicStructureRepresentation* es) {
    const uint64_t h = transfer_function_hash(es->volume.tf);
    if (es->volume.tf_tex && h == es->volume.tf_hash) return es->volume.tf_tex;
    enum { RES = 256 };
    float px[RES * 4];
    transfer_function_bake_rgba(es->volume.tf, px, RES);
    if (!es->volume.tf_tex) gl::init_texture_2D(&es->volume.tf_tex, RES, 1, GL_RGBA32F);
    gl::set_texture_2D_data(es->volume.tf_tex, 0, px, GL_RGBA32F);
    es->volume.tf_hash = h;
    return es->volume.tf_tex;
}

// ---------------------------------------------------------------------------------------------------------------
// The gizmo
// ---------------------------------------------------------------------------------------------------------------

// The gizmo of the picked plane, in screen space at the middle of where the plane cuts its box (its anchor):
//
//   - a ring: dragging within it turns the normal as a trackball turns, about the axis in the view across the drag,
//     by the drag's length over the ring's radius; from any view, the plane seen face on or edge on alike. The
//     normal has two degrees of freedom and this turns only those: there is no turn about the normal to keep.
//     With shift, the normal snaps to the nearest axis within 10 degrees.
//   - the normal, an arrow from the anchor, and a handle at the anchor: dragging either slides the plane along the
//     normal, by where the mouse's ray passes closest to the normal's line; when the normal points (nearly) at the
//     viewer, by the mouse's vertical movement instead. Drawn faint when it points away from the viewer.
//
// The anchor is held while a drag lasts: the plane turns about it, and slides along its normal through it. After,
// it stays where the drag left it (as long as it is in the box, projected onto the plane as it moves), rather than
// jump to the middle of the plane's new cut.

enum GizmoPart { GizmoPart_None, GizmoPart_Rotate, GizmoPart_Slide };

static struct {
    GizmoPart part = GizmoPart_None;    // of the drag under way, None when there is none
    uint32_t  rep_id = 0;
    int       plane = -1;
    vec3_t    anchor = {};              // on the plane, held while the drag lasts
    ImVec2    last_mouse = {};
    // a slide: the anchor, the normal and the plane's d when it began; where along the normal from that anchor the
    // mouse took it (NaN: by the mouse's vertical movement, the normal pointing at the viewer), and that movement
    vec3_t    slide_a0 = {};
    vec3_t    slide_n = {};
    float     slide_d0 = 0.0f;
    float     slide_t = 0.0f;
    float     slide_off = 0.0f;
} gizmo_drag;

bool es_planes_dragging(uint32_t rep_id, int plane) {
    return gizmo_drag.part != GizmoPart_None && gizmo_drag.rep_id == rep_id && gizmo_drag.plane == plane;
}

// Where the last drag left the anchor, of this plane
static struct {
    uint32_t rep_id = 0;
    int      plane = -1;
    vec3_t   point = {};
} gizmo_anchor;

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

// Which way p faces the camera: true when the camera is on the side its normal points to. In perspective, by where the
// camera is (the side of the plane it is on), so that it turns exactly as the plane is seen edge on; in orthographic,
// by the direction it looks along, as every point sees the plane alike. In the space the planes are in (the view
// matrix's world).
static bool faces_camera(const ApplicationState* state, const Plane& p) {
    const mat4_t inv = mat4_inverse(state->view.param.matrix.curr.view);
    if (state->view.mode == CameraMode::Orthographic) {
        const vec4_t back = mat4_mul_vec4(inv, vec4_set(0, 0, 1, 0));
        return vec3_dot(p.normal, vec3_set(back.x, back.y, back.z)) >= 0.0f;
    }
    const vec4_t eye = mat4_mul_vec4(inv, vec4_set(0, 0, 0, 1));
    return plane_distance(p, vec3_set(eye.x, eye.y, eye.z)) >= 0.0f;
}

// The outlines es_planes_gizmo found this frame, as segments in normalized device coordinates (z for the depth test),
// for es_planes_draw_outlines; the gizmo's own parts are drawn by ImGui, over everything
struct OutlineSegment {
    vec3_t a, b;
    ImU32  col, halo;
    float  width;      // points
};
static OutlineSegment outline_segs[256];
static int            outline_count = 0;

static inline bool project_ndc(vec3_t* out, const mat4_t& world_to_clip, vec3_t p) {
    const vec4_t c = mat4_mul_vec4(world_to_clip, vec4_set(p.x, p.y, p.z, 1.0f));
    if (c.w <= 1.0e-6f) return false;
    *out = vec3_set(c.x / c.w, c.y / c.w, c.z / c.w);
    return true;
}

static Representation* gizmo_target(ApplicationState* state, uint32_t rep_id, int plane) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        Representation& rep = state->representation.reps[i];
        if (rep.node.id != rep_id) continue;
        const ElectronicStructureRepresentation& es = rep.electronic_structure;
        if (rep.type != RepresentationType::ElectronicStructure || !representation_shown(*state, rep) ||
            es.display == ElectronicStructureDisplay::Isosurface || plane < 0 || plane >= electronic_structure_planes(es).count) {
            return nullptr;
        }
        return &rep;
    }
    return nullptr;
}

// The mouse's ray in world space
static void mouse_ray(vec3_t* origin, vec3_t* dir, const mat4_t& world_to_clip, ImVec2 mouse, vec2_t canvas_min, vec2_t canvas_size) {
    const mat4_t inv = mat4_inverse(world_to_clip);
    const float x = (mouse.x - canvas_min.x) / canvas_size.x * 2.0f - 1.0f;
    const float y = 1.0f - (mouse.y - canvas_min.y) / canvas_size.y * 2.0f;
    vec4_t a = mat4_mul_vec4(inv, vec4_set(x, y, -1.0f, 1.0f));
    vec4_t b = mat4_mul_vec4(inv, vec4_set(x, y,  1.0f, 1.0f));
    const vec3_t pa = vec3_set(a.x / a.w, a.y / a.w, a.z / a.w);
    const vec3_t pb = vec3_set(b.x / b.w, b.y / b.w, b.z / b.w);
    *origin = pa;
    *dir = vec3_normalize(vec3_sub(pb, pa));
}

// Along the line c + t n (n unit), where the ray o + s r passes closest: false when the two are (nearly) parallel
static bool line_ray_closest(float* t, vec3_t c, vec3_t n, vec3_t o, vec3_t r) {
    const vec3_t w = vec3_sub(c, o);
    const float b = vec3_dot(n, r);
    const float den = 1.0f - b * b;
    if (den < 1.0e-3f) return false;
    *t = (b * vec3_dot(r, w) - vec3_dot(n, w)) / den;
    return true;
}

// v turned about the unit axis by angle (Rodrigues)
static vec3_t rotate_about(vec3_t v, vec3_t axis, float angle) {
    const float c = cosf(angle), s = sinf(angle);
    return vec3_add(vec3_add(vec3_mul1(v, c), vec3_mul1(vec3_cross(axis, v), s)), vec3_mul1(axis, vec3_dot(axis, v) * (1.0f - c)));
}

// The nearest axis to n within 10 degrees, or n
static vec3_t snap_to_axis(vec3_t n) {
    int k = 0;
    for (int i = 1; i < 3; ++i) if (fabsf(n.elem[i]) > fabsf(n.elem[k])) k = i;
    if (fabsf(n.elem[k]) < cosf(10.0f * 3.14159265f / 180.0f)) return n;
    vec3_t a = vec3_zero();
    a.elem[k] = n.elem[k] < 0.0f ? -1.0f : 1.0f;
    return a;
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
    const bool dragging = gizmo_drag.part != GizmoPart_None;
    const bool mouse_taken = ImGui::GetActiveID() != 0 && !dragging;

    auto& sel = state->representation.gizmo;
    Representation* target = gizmo_target(state, sel.rep_id, sel.plane);
    if (!target) {
        sel.plane = -1;
        gizmo_drag.part = GizmoPart_None;
    }

    // The outlines, and the one the mouse is on
    outline_count = 0;
    uint32_t hover_rep = 0;
    int hover_plane = -1;
    float hover_dist = 6.0f;
    vec3_t centre = vec3_zero();
    bool centre_valid = false;
    for (size_t r = 0; r < md_array_size(state->representation.reps); ++r) {
        const Representation& rep = state->representation.reps[r];
        if (rep.type != RepresentationType::ElectronicStructure || !representation_shown(*state, rep) || !rep.type_is_valid) continue;
        const ElectronicStructureRepresentation& es = rep.electronic_structure;
        if (es.display == ElectronicStructureDisplay::Isosurface) continue;
        mat4_t box;
        if (!es_volume_box(&box, es)) continue;
        const ElectronicStructurePlaneSet& set = electronic_structure_planes(es);
        for (int i = 0; i < set.count; ++i) {
            vec3_t poly[6];
            const int n = plane_box_polygon(poly, set.items[i].plane, box);
            if (n == 0) continue;
            const bool selected = target == &rep && sel.plane == i;
            if (selected) {
                for (int k = 0; k < n; ++k) centre = vec3_add(centre, poly[k]);
                centre = vec3_mul1(centre, 1.0f / (float)n);
                centre_valid = true;
            }
            ImVec2 sp[6];
            vec3_t ndc[6];
            bool all = true;
            for (int k = 0; k < n; ++k) all &= project(&sp[k], world_to_clip, poly[k], canvas_min, canvas_size) && project_ndc(&ndc[k], world_to_clip, poly[k]);
            if (!all) continue;
            const bool enabled = set.items[i].enabled;
            // a light line on a dark halo: seen over a light background and a dark one alike; hidden behind what is in
            // front of it (es_planes_draw_outlines)
            const ImU32 col = selected ? IM_COL32(255, 200, 40, 255) : enabled ? IM_COL32(255, 255, 255, 220) : IM_COL32(200, 200, 200, 130);
            const ImU32 halo = enabled || selected ? IM_COL32(0, 0, 0, 110) : IM_COL32(0, 0, 0, 60);
            const float width = selected ? 2.0f : 1.5f;
            for (int k = 0; k < n && outline_count < (int)ARRAY_SIZE(outline_segs); ++k) {
                outline_segs[outline_count++] = OutlineSegment{ndc[k], ndc[(k + 1) % n], col, halo, width};
            }
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

    bool capture = false;
    GizmoPart hover_part = GizmoPart_None;
    if (target && (centre_valid || dragging)) {
        Plane& plane = electronic_structure_planes(target->electronic_structure).items[sel.plane].plane;
        const bool ours = dragging && gizmo_drag.rep_id == sel.rep_id && gizmo_drag.plane == sel.plane;
        if (dragging && !ours) gizmo_drag.part = GizmoPart_None;
        // A slice's normal is of no matter which way it points: it is turned to the camera, so that the arrow is
        // seen and taken hold of (not while it is dragged, which reads it)
        if (target->electronic_structure.display == ElectronicStructureDisplay::Planes && !ours && !faces_camera(state, plane)) {
            plane_flip(&plane);
        }
        vec3_t anchor = ours ? gizmo_drag.anchor : plane_project(plane, centre);
        if (!ours && gizmo_anchor.rep_id == sel.rep_id && gizmo_anchor.plane == sel.plane) {
            mat4_t box;
            if (es_volume_box(&box, target->electronic_structure)) {
                const vec3_t p = plane_project(plane, gizmo_anchor.point);
                const vec4_t u = mat4_mul_vec4(mat4_inverse(box), vec4_set(p.x, p.y, p.z, 1.0f));
                if (u.x >= 0.0f && u.x <= 1.0f && u.y >= 0.0f && u.y <= 1.0f && u.z >= 0.0f && u.z <= 1.0f) anchor = p;
            }
        }
        const vec3_t n = plane.normal;

        // its size on screen, and the normal's length in the world that makes the arrow as long as the ring's radius
        const float ring_radius = 64.0f;
        const vec3_t right = quat_mul_vec3(state->view.camera.orientation, vec3_set(1, 0, 0));
        ImVec2 sc, sr;
        const bool visible = project(&sc, world_to_clip, anchor, canvas_min, canvas_size) &&
                             project(&sr, world_to_clip, vec3_add(anchor, right), canvas_min, canvas_size);
        const float px_per_unit = visible ? sqrtf((sr.x - sc.x) * (sr.x - sc.x) + (sr.y - sc.y) * (sr.y - sc.y)) : 0.0f;
        if (visible && px_per_unit > 0.0f) {
            const float arrow_len = ring_radius / px_per_unit;
            ImVec2 tip;
            project(&tip, world_to_clip, vec3_add(anchor, vec3_mul1(n, arrow_len)), canvas_min, canvas_size);
            const float arrow_px = sqrtf((tip.x - sc.x) * (tip.x - sc.x) + (tip.y - sc.y) * (tip.y - sc.y));
            const bool facing = faces_camera(state, plane);

            // what the mouse is over: the handle, the arrow, then the ring's disc
            const float md = sqrtf((mouse.x - sc.x) * (mouse.x - sc.x) + (mouse.y - sc.y) * (mouse.y - sc.y));
            if (!dragging && canvas_hovered && !mouse_taken) {
                if (md <= 8.0f || (arrow_px > 12.0f && segment_distance(mouse, sc, tip) <= 6.0f)) hover_part = GizmoPart_Slide;
                else if (md <= ring_radius + 4.0f) hover_part = GizmoPart_Rotate;
            }
            const GizmoPart lit = ours ? gizmo_drag.part : hover_part;

            // the ring
            const ImU32 ring_col = lit == GizmoPart_Rotate ? IM_COL32(255, 220, 90, 255) : IM_COL32(255, 200, 40, 170);
            if (lit == GizmoPart_Rotate) dl->AddCircleFilled(sc, ring_radius, IM_COL32(255, 200, 40, 28), 48);
            dl->AddCircle(sc, ring_radius, IM_COL32(0, 0, 0, 90), 48, 4.0f);
            dl->AddCircle(sc, ring_radius, ring_col, 48, 2.0f);
            // the normal: an arrow, faint behind the plane (pointing away)
            const ImU32 arrow_col = lit == GizmoPart_Slide ? IM_COL32(120, 200, 255, 255) : facing ? IM_COL32(60, 150, 255, 235) : IM_COL32(60, 150, 255, 110);
            if (arrow_px > 2.0f) {
                const ImVec2 d = ImVec2((tip.x - sc.x) / arrow_px, (tip.y - sc.y) / arrow_px);
                const ImVec2 o = ImVec2(-d.y, d.x);
                const float head = MIN(12.0f, 0.4f * arrow_px);
                const ImVec2 base = ImVec2(tip.x - d.x * head, tip.y - d.y * head);
                dl->AddLine(sc, base, IM_COL32(0, 0, 0, 110), 5.0f);
                dl->AddLine(sc, base, arrow_col, 3.0f);
                dl->AddTriangleFilled(tip, ImVec2(base.x + o.x * head * 0.5f, base.y + o.y * head * 0.5f),
                                      ImVec2(base.x - o.x * head * 0.5f, base.y - o.y * head * 0.5f), arrow_col);
            }
            dl->AddCircleFilled(sc, 6.0f, IM_COL32(0, 0, 0, 120), 16);
            dl->AddCircleFilled(sc, 4.5f, arrow_col, 16);

            // a drag begun
            if (!dragging && hover_part != GizmoPart_None && ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
                gizmo_drag.part = hover_part;
                gizmo_drag.rep_id = sel.rep_id;
                gizmo_drag.plane = sel.plane;
                gizmo_drag.anchor = anchor;
                gizmo_drag.last_mouse = mouse;
                gizmo_drag.slide_a0 = anchor;
                gizmo_drag.slide_n = n;
                gizmo_drag.slide_d0 = plane.d;
                gizmo_drag.slide_off = 0.0f;
                vec3_t o, r;
                mouse_ray(&o, &r, world_to_clip, mouse, canvas_min, canvas_size);
                float t = 0.0f;
                gizmo_drag.slide_t = line_ray_closest(&t, anchor, n, o, r) && arrow_px > 12.0f ? t : NAN;
            }
            // a drag going on
            if (gizmo_drag.part != GizmoPart_None && gizmo_drag.rep_id == sel.rep_id && gizmo_drag.plane == sel.plane) {
                if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                    gizmo_drag.part = GizmoPart_None;
                    gizmo_anchor.rep_id = sel.rep_id;
                    gizmo_anchor.plane = sel.plane;
                    gizmo_anchor.point = gizmo_drag.anchor;
                } else if (gizmo_drag.part == GizmoPart_Rotate) {
                    const ImVec2 dm = ImVec2(mouse.x - gizmo_drag.last_mouse.x, mouse.y - gizmo_drag.last_mouse.y);
                    const float len = sqrtf(dm.x * dm.x + dm.y * dm.y);
                    if (len > 0.0f) {
                        // the front of a ball under the mouse moves with it: about (dy, dx, 0) in the view, by len / radius
                        const vec3_t axis = vec3_normalize(quat_mul_vec3(state->view.camera.orientation, vec3_set(dm.y, dm.x, 0.0f)));
                        vec3_t nn = vec3_normalize(rotate_about(plane.normal, axis, len / ring_radius));
                        plane = plane_from_point_normal(gizmo_drag.anchor, nn);
                    }
                    if (ImGui::GetIO().KeyShift) plane = plane_from_point_normal(gizmo_drag.anchor, snap_to_axis(plane.normal));
                } else {
                    // along the normal it was taken with, from where it was taken
                    const vec3_t sn = gizmo_drag.slide_n;
                    float t = 0.0f;
                    vec3_t o, r;
                    mouse_ray(&o, &r, world_to_clip, mouse, canvas_min, canvas_size);
                    if (gizmo_drag.slide_t == gizmo_drag.slide_t && line_ray_closest(&t, gizmo_drag.slide_a0, sn, o, r)) {
                        t -= gizmo_drag.slide_t;
                    } else {
                        // seen along the normal: up the screen is along it
                        gizmo_drag.slide_off += (gizmo_drag.last_mouse.y - mouse.y) / px_per_unit;
                        t = gizmo_drag.slide_off;
                    }
                    gizmo_drag.anchor = vec3_add(gizmo_drag.slide_a0, vec3_mul1(sn, t));
                    plane.normal = sn;
                    plane.d = gizmo_drag.slide_d0 - t;
                }
                gizmo_drag.last_mouse = mouse;
                capture = true;
            }
            if (hover_part != GizmoPart_None) {
                capture = true;
                ImGui::SetMouseCursor(hover_part == GizmoPart_Slide ? ImGuiMouseCursor_ResizeNS : ImGuiMouseCursor_Hand);
            }
        } else if (ours) {
            // off screen while dragged: the drag ends with the button
            if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) gizmo_drag.part = GizmoPart_None;
            capture = true;
        }
    }

    // A click on an outline picks its plane
    if (!capture && !mouse_taken && canvas_hovered && hover_plane >= 0) {
        ImGui::SetMouseCursor(ImGuiMouseCursor_Hand);
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) {
            sel.rep_id = hover_rep;
            sel.plane = hover_plane;
            gizmo_drag.part = GizmoPart_None;
            capture = true;
        }
    }
    return capture;
}

// ---------------------------------------------------------------------------------------------------------------------
// The outlines drawn: each segment a quad of its width in pixels, the halos under the lines, over the image presented,
// each fragment kept where nothing of the scene is in front of it (the G-buffer's depth). They write no depth and are
// not part of the scene: they change nothing in it, as the gizmo does not.

static struct {
    GLuint program = 0, vao = 0, vbo = 0;
    GLint  loc_depth = -1, loc_scale = -1;
    bool   failed = false;
} outline_gl;

struct OutlineVertex {
    float    x, y, z;
    uint32_t col;
    float    across;    // pixels from the line's middle, signed
    float    half;      // half its width, pixels
};

static const char* outline_vs = R"(
#version 330 core
layout(location = 0) in vec3 a_pos;
layout(location = 1) in vec4 a_col;
layout(location = 2) in vec2 a_edge;    // across, half
out vec4 v_col;
out vec2 v_edge;
void main() {
    v_col = a_col;
    v_edge = a_edge;
    gl_Position = vec4(a_pos, 1.0);
}
)";

static const char* outline_fs = R"(
#version 330 core
uniform sampler2D u_depth;      // the scene's, window depth
uniform vec2 u_scale;           // its texels per framebuffer pixel
in vec4 v_col;
in vec2 v_edge;
out vec4 out_frag;
void main() {
    float scene = texelFetch(u_depth, ivec2(gl_FragCoord.xy * u_scale), 0).r;
    // behind the scene (a margin for the scene's jitter and the depth's rounding)
    if (gl_FragCoord.z > scene + 2.0e-6) discard;
    // its edges feathered over a pixel
    float cover = clamp(v_edge.y + 0.5 - abs(v_edge.x), 0.0, 1.0);
    out_frag = vec4(v_col.rgb, v_col.a * cover);
}
)";

static bool outline_gl_init() {
    if (outline_gl.program || outline_gl.failed) return outline_gl.program != 0;
    GLuint sh[2] = {
        gl::compile_shader_from_source(str_from_cstr(outline_vs), GL_VERTEX_SHADER),
        gl::compile_shader_from_source(str_from_cstr(outline_fs), GL_FRAGMENT_SHADER),
    };
    if (!sh[0] || !sh[1]) {
        outline_gl.failed = true;
    } else {
        outline_gl.program = glCreateProgram();
        if (!gl::attach_link_detach(outline_gl.program, sh, 2)) {
            glDeleteProgram(outline_gl.program);
            outline_gl.program = 0;
            outline_gl.failed = true;
        }
    }
    if (sh[0]) glDeleteShader(sh[0]);
    if (sh[1]) glDeleteShader(sh[1]);
    if (outline_gl.failed) return false;
    outline_gl.loc_depth = glGetUniformLocation(outline_gl.program, "u_depth");
    outline_gl.loc_scale = glGetUniformLocation(outline_gl.program, "u_scale");
    glGenVertexArrays(1, &outline_gl.vao);
    glGenBuffers(1, &outline_gl.vbo);
    glBindVertexArray(outline_gl.vao);
    glBindBuffer(GL_ARRAY_BUFFER, outline_gl.vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(OutlineVertex), (const GLvoid*)offsetof(OutlineVertex, x));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(OutlineVertex), (const GLvoid*)offsetof(OutlineVertex, col));
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(OutlineVertex), (const GLvoid*)offsetof(OutlineVertex, across));
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    return true;
}

void es_planes_draw_outlines(uint32_t depth_tex, int depth_width, int depth_height, int fb_width, int fb_height, float px_per_point) {
    if (outline_count == 0 || fb_width <= 0 || fb_height <= 0 || !depth_tex || !outline_gl_init()) return;
    // two quads (halo, line) of 6 vertices per segment; the halos first
    static OutlineVertex verts[ARRAY_SIZE(outline_segs) * 12];
    int nv = 0;
    const vec2_t half = vec2_set(0.5f * (float)fb_width, 0.5f * (float)fb_height);
    for (int pass = 0; pass < 2; ++pass) {
        for (int i = 0; i < outline_count; ++i) {
            const OutlineSegment& s = outline_segs[i];
            // in pixels: along it and across it, half the width out each way, and as far past each end (square caps
            // that close the corners)
            const vec2_t a = vec2_set(s.a.x * half.x, s.a.y * half.y), b = vec2_set(s.b.x * half.x, s.b.y * half.y);
            vec2_t d = vec2_sub(b, a);
            const float len = vec2_length(d);
            if (!(len > 1.0e-3f)) continue;
            d = vec2_mul1(d, 1.0f / len);
            const float w = 0.5f * (pass == 0 ? s.width + 2.0f : s.width) * px_per_point;
            const float we = w + 1.0f;     // a pixel more, feathered (outline_fs)
            const vec2_t o = vec2_set(-d.y * we, d.x * we), e = vec2_mul1(d, w);
            const vec2_t pa = vec2_sub(a, e), pb = vec2_add(b, e);
            const vec2_t q[4] = { vec2_add(pa, o), vec2_sub(pa, o), vec2_sub(pb, o), vec2_add(pb, o) };
            const float z[4] = { s.a.z, s.a.z, s.b.z, s.b.z };
            const float across[4] = { we, -we, -we, we };
            const uint32_t col = pass == 0 ? s.halo : s.col;
            static const int idx[6] = {0, 1, 2, 0, 2, 3};
            for (int k = 0; k < 6; ++k) {
                const int j = idx[k];
                verts[nv++] = OutlineVertex{q[j].x / half.x, q[j].y / half.y, z[j], col, across[j], w};
            }
        }
    }
    if (nv == 0) return;

    glUseProgram(outline_gl.program);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, depth_tex);
    glUniform1i(outline_gl.loc_depth, 0);
    glUniform2f(outline_gl.loc_scale, (float)depth_width / (float)fb_width, (float)depth_height / (float)fb_height);
    glBindVertexArray(outline_gl.vao);
    glBindBuffer(GL_ARRAY_BUFFER, outline_gl.vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(OutlineVertex) * nv, verts, GL_STREAM_DRAW);
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glDisable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDrawArrays(GL_TRIANGLES, 0, nv);
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glUseProgram(0);
}

void es_planes_shutdown() {
    if (outline_gl.program) glDeleteProgram(outline_gl.program);
    if (outline_gl.vbo) glDeleteBuffers(1, &outline_gl.vbo);
    if (outline_gl.vao) glDeleteVertexArrays(1, &outline_gl.vao);
    outline_gl = {};
}
