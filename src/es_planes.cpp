#include <es_planes.h>

#include <viamd.h>
#include <representation.h>
#include <gfx/gl_utils.h>

#include <imgui.h>
#include <imgui_internal.h>
#include <implot.h>

#include <math.h>
#include <float.h>
#include <core/md_vec_math.h>
#include <md_util_geom.h>

static inline md_float3_t transform_point(const md_float4x4_t& m, md_float3_t p) {
    const md_float4_t r = md_float4x4_mul_float4(m, md_float4_set(p.x, p.y, p.z, 1.0f));
    return md_float3_set(r.x, r.y, r.z);
}


bool es_volume_box(md_float4x4_t* box_to_world, const ElectronicStructureRepresentation& es) {
    if (!es.has_box) return false;
    *box_to_world = es.box_to_world;
    return true;
}

// The rectangle around where the plane cuts the box, with its first side along 'right' (as near as the plane allows);
// where the plane misses the box, a square about the box's middle, as wide as the box's widest side
static void slice_fit_box(ElectronicStructurePlane* p, const md_float4x4_t& box, md_float3_t right) {
    md_float3_t u, v;
    plane_frame(&u, &v, p->plane, right);
    md_float3_t poly[6];
    const int n = plane_box_polygon(poly, p->plane, box);
    float r[4] = {FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX};
    if (n >= 3) {
        for (int k = 0; k < n; ++k) {
            const float s = md_float3_dot(poly[k], u), t = md_float3_dot(poly[k], v);
            r[0] = MD_MIN(r[0], s); r[1] = MD_MAX(r[1], s);
            r[2] = MD_MIN(r[2], t); r[3] = MD_MAX(r[3], t);
        }
    } else {
        const md_float3_t c = plane_project(p->plane, transform_point(box, md_float3_set(0.5f, 0.5f, 0.5f)));
        float side = 0.0f;
        for (int k = 0; k < 3; ++k) side = MD_MAX(side, md_float3_length(md_float3_set(box.elem[k][0], box.elem[k][1], box.elem[k][2])));
        r[0] = md_float3_dot(c, u) - 0.5f * side; r[1] = r[0] + side;
        r[2] = md_float3_dot(c, v) - 0.5f * side; r[3] = r[2] + side;
    }
    // the plane's point nearest the origin, and from it along u and v to the rectangle's middle
    const md_float3_t on = md_float3_mul1(p->plane.normal, -p->plane.d);
    p->center = md_float3_add(on, md_float3_add(md_float3_mul1(u, 0.5f * (r[0] + r[1])), md_float3_mul1(v, 0.5f * (r[2] + r[3]))));
    p->axis = u;
    p->extent = md_float2_set(MD_MAX(0.5f * (r[1] - r[0]), 0.1f), MD_MAX(0.5f * (r[3] - r[2]), 0.1f));
}

void es_slices_conform(ElectronicStructureRepresentation* es) {
    md_float4x4_t box;
    const bool have_box = es_volume_box(&box, *es);
    for (int i = 0; i < es->slices.count && i < ES_MAX_PLANES; ++i) {
        ElectronicStructurePlane& p = es->slices.items[i];
        if (!(p.extent.x > 0.0f && p.extent.y > 0.0f)) {
            // by the box's first axis, not the view: as it was before slices had rectangles of their own
            if (have_box) slice_fit_box(&p, box, md_float3_set(box.elem[0][0], box.elem[0][1], box.elem[0][2]));
            continue;
        }
        // only what is off by more than rounding: a slice already on its plane is left exactly as it is (what the
        // scene's hash and its image's key are of)
        const float dist = plane_distance(p.plane, p.center);
        if (fabsf(dist) > 1.0e-5f) p.center = md_float3_sub(p.center, md_float3_mul1(p.plane.normal, dist));
        const float along = md_float3_dot(p.axis, p.plane.normal);
        const float len = md_float3_length(p.axis);
        if (fabsf(along) > 1.0e-5f || fabsf(len - 1.0f) > 1.0e-5f) {
            md_float3_t u, v;
            plane_frame(&u, &v, p.plane, p.axis);
            p.axis = u;
        }
    }
}

bool es_slice_frame(EsSliceFrame* f, const ElectronicStructurePlane& p) {
    if (!(p.extent.x > 0.0f && p.extent.y > 0.0f)) return false;
    // Which way the normal points is of no matter to the image (a slice's is turned to the camera as it is seen): one
    // of the two, so that turning it does not make the image again
    md_float3_t n = p.plane.normal;
    if (md_float3_dot(n, md_float3_set(0.5773f, 0.6172f, 0.5345f)) < 0.0f) n = md_float3_mul1(n, -1.0f);
    f->n = n;
    f->u = p.axis;
    f->v = md_float3_cross(n, p.axis);
    f->c = p.center;
    f->hu = p.extent.x;
    f->hv = p.extent.y;
    return true;
}

md_float4x4_t es_slice_rect_from_world(const EsSliceFrame& f) {
    // md_float4x4_t's elem[column][row]
    md_float4x4_t m = {};
    const md_float3_t ax[3] = {md_float3_mul1(f.u, 0.5f / f.hu), md_float3_mul1(f.v, 0.5f / f.hv), f.n};
    for (int r = 0; r < 3; ++r) {
        for (int k = 0; k < 3; ++k) m.elem[k][r] = ax[r].elem[k];
        m.elem[3][r] = -md_float3_dot(ax[r], f.c) + (r < 2 ? 0.5f : 0.0f);
    }
    m.elem[3][3] = 1.0f;
    return m;
}

void es_slice_image_grid(md_grid_t* grid, const EsSliceFrame& f, float spacing) {
    const int w = MD_CLAMP((int)ceilf(2.0f * f.hu / spacing), 1, GPU_PLANE_DIM);
    const int h = MD_CLAMP((int)ceilf(2.0f * f.hv / spacing), 1, GPU_PLANE_DIM);
    // the pixels' edges on the rectangle's: each axis its own spacing, near the one asked
    const float su = 2.0f * f.hu / (float)w, sv = 2.0f * f.hv / (float)h, sn = MD_MIN(su, sv);
    // Pixel centres on the plane: the one layer's centre (the kernels sample at voxel centres) a half spacing along n
    // from the corner
    const md_float3_t corner = md_float3_sub(md_float3_sub(md_float3_sub(f.c, md_float3_mul1(f.u, f.hu)), md_float3_mul1(f.v, f.hv)), md_float3_mul1(f.n, 0.5f * sn));
    const float to_bohr = (float)(1.0 / BOHR_TO_ANGSTROM);
    *grid = md_grid_t{};
    grid->orientation = md_float3x3_t{};
    for (int k = 0; k < 3; ++k) {
        grid->orientation.elem[0][k] = f.u.elem[k];
        grid->orientation.elem[1][k] = f.v.elem[k];
        grid->orientation.elem[2][k] = f.n.elem[k];
    }
    grid->origin  = md_float3_mul1(corner, to_bohr);
    grid->spacing = md_float3_set(su * to_bohr, sv * to_bohr, sn * to_bohr);
    grid->dim[0] = w;
    grid->dim[1] = h;
    grid->dim[2] = 1;
}

static bool faces_camera(const ApplicationState* state, const Plane& p);

int es_planes_add(ApplicationState* state, ElectronicStructureRepresentation* es, uint32_t rep_id) {
    md_float4x4_t box;
    ElectronicStructurePlaneSet& set = electronic_structure_planes(*es);
    if (set.count >= ES_MAX_PLANES || !es_volume_box(&box, *es)) return -1;
    const md_float3_t centre = transform_point(box, md_float3_set(0.5f, 0.5f, 0.5f));
    const md_float3_t forward = md_quat_mul_float3(state->view.camera.orientation, md_float3_set(0, 0, -1));
    const md_float3_t right = md_quat_mul_float3(state->view.camera.orientation, md_float3_set(1, 0, 0));
    const bool slice = es->display == ElectronicStructureDisplay::Planes;
    const int i = set.count++;
    // a cut facing away from the camera keeps the half behind it; a slice faces the camera, its rectangle around where
    // it cuts the box, square to the view
    ElectronicStructurePlane& p = set.items[i];
    p = ElectronicStructurePlane{};
    p.plane = plane_from_point_normal(centre, slice ? md_float3_mul1(forward, -1.0f) : forward);
    if (slice) {
        if (!faces_camera(state, p.plane)) plane_flip(&p.plane);
        slice_fit_box(&p, box, right);
    }
    state->representation.gizmo.rep_id = rep_id;
    state->representation.gizmo.plane = i;
    return i;
}

Plane es_plane_fit_selection(const ApplicationState* state, const Plane& plane, PlaneFitInfo* info, md_float3_t* out_centre) {
    const md_system_state_t& st = state->mold.state;
    const size_t count = md_bitmask_count(&state->selection.selection_mask);
    if (count == 0 || !st.xyz) {
        if (info) *info = PlaneFitInfo{};
        return plane;
    }
    md_temp_scope_t temp = md_temp_begin();
    md_float3_t* pts = (md_float3_t*)md_temp_alloc(temp, sizeof(md_float3_t) * count);
    size_t n = 0;
    md_bitset_iter_t it = md_bitmask_iter(&state->selection.selection_mask);
    while (md_bitset_iter_next(&it) && n < count) {
        if (it.idx < st.num_atoms) pts[n++] = st.xyz[it.idx];
    }
    if (n > 1 && md_unitcell_flags(&st.unitcell) != 0) {
        // across the boundary: each the image nearest to the first
        for (size_t i = 1; i < n; ++i) pts[i] = md_float3_sub(pts[i], pts[0]);
        md_util_min_image((float*)(pts + 1), 0, n - 1, &st.unitcell);
        for (size_t i = 1; i < n; ++i) pts[i] = md_float3_add(pts[i], pts[0]);
    }
    const Plane fitted = plane_fit(plane, pts, n, info);
    if (out_centre && n > 0) {
        md_float3_t c = md_float3_zero();
        for (size_t i = 0; i < n; ++i) c = md_float3_add(c, pts[i]);
        *out_centre = md_float3_mul1(c, 1.0f / (float)n);
    }
    md_temp_end(temp);
    return fitted;
}

md_float4_t es_volume_plane(const ElectronicStructureRepresentation& es, const Plane& plane) {
    // n . (A p + t) + d = (A^T n) . p + (n . t + d), with p in the texture's unit cube
    const md_float4x4_t& m = es.box_to_world;
    const md_float3_t n = plane.normal;
    const md_float3_t c0 = md_float3_set(m.elem[0][0], m.elem[0][1], m.elem[0][2]);
    const md_float3_t c1 = md_float3_set(m.elem[1][0], m.elem[1][1], m.elem[1][2]);
    const md_float3_t c2 = md_float3_set(m.elem[2][0], m.elem[2][1], m.elem[2][2]);
    const md_float3_t t  = md_float3_set(m.elem[3][0], m.elem[3][1], m.elem[3][2]);
    return md_float4_set(md_float3_dot(c0, n), md_float3_dot(c1, n), md_float3_dot(c2, n), md_float3_dot(n, t) + plane.d);
}

int es_volume_clip_planes(md_float4_t* out, const ElectronicStructureRepresentation& es) {
    int k = 0;
    for (int i = 0; i < es.cuts.count && i < ES_MAX_PLANES; ++i) {
        if (es.cuts.items[i].enabled) out[k++] = es_volume_plane(es, es.cuts.items[i].plane);
    }
    return k;
}

int es_volume_slices(md_float4_t* out_plane, float* out_opacity, md_float4x4_t* out_rect_from_model, const ElectronicStructureRepresentation& es) {
    const int count = MD_MIN(es.slices.count, ES_MAX_PLANES);
    for (int i = 0; i < count; ++i) {
        const ElectronicStructurePlane& p = es.slices.items[i];
        out_plane[i] = es_volume_plane(es, p.plane);
        EsSliceFrame f;
        if (es_slice_frame(&f, p)) {
            out_opacity[i] = p.enabled ? MD_CLAMP(p.opacity, 0.0f, 1.0f) : 0.0f;
            out_rect_from_model[i] = es_slice_rect_from_world(f) * es.box_to_world;
        } else {
            // no rectangle yet: nowhere
            out_opacity[i] = 0.0f;
            out_rect_from_model[i] = md_float4x4_t{};
            out_rect_from_model[i].elem[3][0] = -1.0f;
            out_rect_from_model[i].elem[3][1] = -1.0f;
            out_rect_from_model[i].elem[3][3] = 1.0f;
        }
    }
    return count;
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

// The gizmo of the picked plane, in screen space at its anchor: a cut's, the middle of where it cuts its box; a slice's,
// the middle of its rectangle.
//
//   - a ring: dragging within it turns the normal as a trackball turns, about the axis in the view across the drag,
//     by the drag's length over the ring's radius (a slice's rectangle turns with it, about its middle). With shift,
//     the normal snaps to the nearest axis within 10 degrees.
//   - the normal, an arrow from the anchor: dragging it slides the plane along the normal, by where the mouse's ray
//     passes closest to the normal's line; when the normal points (nearly) at the viewer, by the mouse's vertical
//     movement instead. Drawn faint when it points away from the viewer.
//   - a cut's handle at the anchor slides it as the arrow does. A slice's moves its rectangle along the plane (with
//     ctrl: slides it as the arrow does, which is how a slice seen face on, its arrow too short to take, is slid).
//   - a slice's rectangle: its sides and corners resize it, the opposite side staying (with shift: about its middle);
//     the knob past its first side turns it about the normal (with shift: by 15 degrees).
//
// What moves along the plane is not offered where the plane is seen (nearly) edge on: the mouse's ray meets it nowhere
// near where it is pointed. A cut's anchor is held while a drag lasts: the plane turns about it, and slides along its
// normal through it. After, it stays where the drag left it (as long as it is in the box, projected onto the plane as
// it moves), rather than jump to the middle of the plane's new cut.

enum GizmoPart { GizmoPart_None, GizmoPart_Rotate, GizmoPart_Slide, GizmoPart_Move, GizmoPart_Resize, GizmoPart_Twist };

static struct {
    GizmoPart part = GizmoPart_None;    // of the drag under way, None when there is none
    uint32_t  rep_id = 0;
    int       plane = -1;
    md_float3_t    anchor = {};              // on the plane, held while the drag lasts
    ImVec2    last_mouse = {};
    // a slide: the anchor, the normal and the plane's d when it began; where along the normal from that anchor the
    // mouse took it (NaN: by the mouse's vertical movement, the normal pointing at the viewer), and that movement
    md_float3_t    slide_a0 = {};
    md_float3_t    slide_n = {};
    float     slide_d0 = 0.0f;
    float     slide_t = 0.0f;
    float     slide_off = 0.0f;
    // a slice's plane and rectangle when it began (u0 and v0 as its handles were: v0 = n x u0 of the normal then), and
    // where on that plane the mouse took it: what a move, a resize and a turn about the normal go from
    Plane     plane0 = {};
    md_float3_t    c0 = {}, u0 = {}, v0 = {};
    float     hu0 = 0.0f, hv0 = 0.0f;
    md_float3_t    grab = {};
    int       su = 0, sv = 0;           // the sides resized: -1, +1 the one along -u, +u (-v, +v); 0 neither
    float     twist0 = 0.0f;            // the angle about the centre, from u0 to v0, the mouse took it at
} gizmo_drag;

bool es_planes_dragging(uint32_t rep_id, int plane) {
    return gizmo_drag.part != GizmoPart_None && gizmo_drag.rep_id == rep_id && gizmo_drag.plane == plane;
}

// Where the last drag left a cut's anchor, of this plane
static struct {
    uint32_t rep_id = 0;
    int      plane = -1;
    md_float3_t   point = {};
} gizmo_anchor;

void es_planes_remove(ApplicationState* state, Representation* rep, int i) {
    ElectronicStructureRepresentation& es = rep->electronic_structure;
    ElectronicStructurePlaneSet& set = electronic_structure_planes(es);
    if (i < 0 || i >= set.count) return;
    for (int k = i; k + 1 < set.count; ++k) set.items[k] = set.items[k + 1];
    set.items[set.count - 1] = ElectronicStructurePlane{};
    set.count -= 1;
    if (&set == &es.slices) {
        // its image with it: what is being evaluated into it first let finish (it would land in a texture freed)
        ElectronicStructurePlaneImage& im = es.plane_images[i];
        bool any = false;
        for (Volume& v : im.img) {
            if (v.tex_id) {
                gpu_volume_job_cancel(state, v.tex_id);
                any = true;
            }
        }
        if (any) gpu_volume_jobs_drain(state);
        for (Volume& v : im.img) if (v.tex_id) gl::free_texture(&v.tex_id);
        for (int k = i; k + 1 < ES_MAX_PLANES; ++k) es.plane_images[k] = es.plane_images[k + 1];
        es.plane_images[ES_MAX_PLANES - 1] = ElectronicStructurePlaneImage{};
    }
    const uint32_t id = rep->node.id;
    auto& sel = state->representation.gizmo;
    if (sel.rep_id == id) {
        if (sel.plane == i) sel.plane = -1;
        else if (sel.plane > i) sel.plane -= 1;
    }
    if (gizmo_drag.rep_id == id) gizmo_drag.part = GizmoPart_None;
    if (gizmo_anchor.rep_id == id) {
        if (gizmo_anchor.plane == i) gizmo_anchor.plane = -1;
        else if (gizmo_anchor.plane > i) gizmo_anchor.plane -= 1;
    }
}

static inline bool project(ImVec2* out, const md_float4x4_t& world_to_clip, md_float3_t p, md_float2_t canvas_min, md_float2_t canvas_size) {
    const md_float4_t c = md_float4x4_mul_float4(world_to_clip, md_float4_set(p.x, p.y, p.z, 1.0f));
    if (c.w <= 1.0e-6f) return false;
    const float x = c.x / c.w, y = c.y / c.w;
    out->x = canvas_min.x + (x * 0.5f + 0.5f) * canvas_size.x;
    out->y = canvas_min.y + (1.0f - (y * 0.5f + 0.5f)) * canvas_size.y;
    return true;
}

static inline float segment_distance(ImVec2 p, ImVec2 a, ImVec2 b) {
    const ImVec2 ab = b - a, ap = p - a;
    const float l2 = ab.x * ab.x + ab.y * ab.y;
    const float t = l2 > 0.0f ? MD_CLAMP((ap.x * ab.x + ap.y * ab.y) / l2, 0.0f, 1.0f) : 0.0f;
    const ImVec2 d = ap - ab * t;
    return sqrtf(d.x * d.x + d.y * d.y);
}

static inline float point_distance(ImVec2 a, ImVec2 b) {
    return sqrtf((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
}

// Which way p faces the camera: true when the camera is on the side its normal points to. In perspective, by where the
// camera is (the side of the plane it is on), so that it turns exactly as the plane is seen edge on; in orthographic,
// by the direction it looks along, as every point sees the plane alike. In the space the planes are in (the view
// matrix's world).
static bool faces_camera(const ApplicationState* state, const Plane& p) {
    const md_float4x4_t inv = md_float4x4_inverse(state->view.param.matrix.curr.view);
    if (state->view.mode == CameraMode::Orthographic) {
        const md_float4_t back = md_float4x4_mul_float4(inv, md_float4_set(0, 0, 1, 0));
        return md_float3_dot(p.normal, md_float3_set(back.x, back.y, back.z)) >= 0.0f;
    }
    const md_float4_t eye = md_float4x4_mul_float4(inv, md_float4_set(0, 0, 0, 1));
    return plane_distance(p, md_float3_set(eye.x, eye.y, eye.z)) >= 0.0f;
}

// The outlines es_planes_gizmo found this frame, as segments in normalized device coordinates (z for the depth test),
// for es_planes_draw_outlines; the gizmo's own parts are drawn by ImGui, over everything. Every cut's (they are not
// drawn otherwise), and the picked slice's: the others' edges are drawn with the slices (the volume rendering), behind
// what of them is in front.
struct OutlineSegment {
    md_float3_t a, b;
    ImU32  col, halo;
    float  width;      // points
};
static OutlineSegment outline_segs[256];
static int            outline_count = 0;

static inline bool project_ndc(md_float3_t* out, const md_float4x4_t& world_to_clip, md_float3_t p) {
    const md_float4_t c = md_float4x4_mul_float4(world_to_clip, md_float4_set(p.x, p.y, p.z, 1.0f));
    if (c.w <= 1.0e-6f) return false;
    *out = md_float3_set(c.x / c.w, c.y / c.w, c.z / c.w);
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
static void mouse_ray(md_float3_t* origin, md_float3_t* dir, const md_float4x4_t& world_to_clip, ImVec2 mouse, md_float2_t canvas_min, md_float2_t canvas_size) {
    const md_float4x4_t inv = md_float4x4_inverse(world_to_clip);
    const float x = (mouse.x - canvas_min.x) / canvas_size.x * 2.0f - 1.0f;
    const float y = 1.0f - (mouse.y - canvas_min.y) / canvas_size.y * 2.0f;
    md_float4_t a = md_float4x4_mul_float4(inv, md_float4_set(x, y, -1.0f, 1.0f));
    md_float4_t b = md_float4x4_mul_float4(inv, md_float4_set(x, y,  1.0f, 1.0f));
    const md_float3_t pa = md_float3_set(a.x / a.w, a.y / a.w, a.z / a.w);
    const md_float3_t pb = md_float3_set(b.x / b.w, b.y / b.w, b.z / b.w);
    *origin = pa;
    *dir = md_float3_normalize(md_float3_sub(pb, pa));
}

// Along the line c + t n (n unit), where the ray o + s r passes closest: false when the two are (nearly) parallel
static bool line_ray_closest(float* t, md_float3_t c, md_float3_t n, md_float3_t o, md_float3_t r) {
    const md_float3_t w = md_float3_sub(c, o);
    const float b = md_float3_dot(n, r);
    const float den = 1.0f - b * b;
    if (den < 1.0e-3f) return false;
    *t = (b * md_float3_dot(r, w) - md_float3_dot(n, w)) / den;
    return true;
}

// Where the ray o + s r (r unit) meets the plane ahead of o: false where it meets it at a grazing angle (within about
// 3 degrees: where it does, the point runs away as the mouse moves a pixel) or behind
static bool ray_plane(md_float3_t* out, const Plane& p, md_float3_t o, md_float3_t r) {
    const float b = md_float3_dot(p.normal, r);
    if (fabsf(b) < 0.05f) return false;
    const float s = -plane_distance(p, o) / b;
    if (s < 0.0f) return false;
    *out = md_float3_add(o, md_float3_mul1(r, s));
    return true;
}

// v turned about the unit axis by angle (Rodrigues)
static md_float3_t rotate_about(md_float3_t v, md_float3_t axis, float angle) {
    const float c = cosf(angle), s = sinf(angle);
    return md_float3_add(md_float3_add(md_float3_mul1(v, c), md_float3_mul1(md_float3_cross(axis, v), s)), md_float3_mul1(axis, md_float3_dot(axis, v) * (1.0f - c)));
}

// The nearest axis to n within 10 degrees, or n
static md_float3_t snap_to_axis(md_float3_t n) {
    int k = 0;
    for (int i = 1; i < 3; ++i) if (fabsf(n.elem[i]) > fabsf(n.elem[k])) k = i;
    if (fabsf(n.elem[k]) < cosf(10.0f * 3.14159265f / 180.0f)) return n;
    md_float3_t a = md_float3_zero();
    a.elem[k] = n.elem[k] < 0.0f ? -1.0f : 1.0f;
    return a;
}

// A slice's rectangle's corners, in order around it: (-u -v), (+u -v), (+u +v), (-u +v)
static void slice_corners(md_float3_t out[4], const ElectronicStructurePlane& p) {
    const md_float3_t u = md_float3_mul1(p.axis, p.extent.x);
    const md_float3_t v = md_float3_mul1(md_float3_cross(p.plane.normal, p.axis), p.extent.y);
    out[0] = md_float3_sub(md_float3_sub(p.center, u), v);
    out[1] = md_float3_sub(md_float3_add(p.center, u), v);
    out[2] = md_float3_add(md_float3_add(p.center, u), v);
    out[3] = md_float3_add(md_float3_sub(p.center, u), v);
}

static inline bool slice_has_rect(const ElectronicStructurePlane& p) { return p.extent.x > 0.0f && p.extent.y > 0.0f; }

// The handles of a slice's rectangle: its 4 corners, the middles of its 4 sides, and the knob that turns it
struct SliceHandle {
    md_float3_t world;
    int    su, sv;
};

bool es_planes_gizmo(ApplicationState* state, md_float2_t canvas_min, md_float2_t canvas_size) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImGuiIO& io = ImGui::GetIO();
    const md_float4x4_t& view = state->view.param.matrix.curr.view;
    const md_float4x4_t& proj = state->view.param.matrix.curr.proj_no_jitter;
    const md_float4x4_t world_to_clip = proj * view;
    const ImVec2 mouse = ImGui::GetMousePos();
    const bool canvas_hovered = ImGui::IsWindowHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
    // Something else has the mouse (the view being turned, a slider being dragged): what started first keeps it, so the
    // planes neither take it nor pick, though the mouse passes over them. ActiveId is that of the frame before here,
    // as the view's surface is submitted after this.
    const bool dragging = gizmo_drag.part != GizmoPart_None;
    const bool mouse_taken = ImGui::GetActiveID() != 0 && !dragging;

    // every slice with its rectangle, on its plane
    for (size_t r = 0; r < md_array_size(state->representation.reps); ++r) {
        Representation& rep = state->representation.reps[r];
        if (rep.type == RepresentationType::ElectronicStructure) es_slices_conform(&rep.electronic_structure);
    }

    auto& sel = state->representation.gizmo;
    Representation* target = gizmo_target(state, sel.rep_id, sel.plane);
    if (!target) {
        sel.plane = -1;
        gizmo_drag.part = GizmoPart_None;
    }

    // Delete (backspace, as macOS calls it delete) removes the picked plane, where no text is being typed
    if (target && !dragging && !io.WantTextInput && ImGui::GetActiveID() == 0 &&
        (ImGui::IsKeyPressed(ImGuiKey_Delete, false) || ImGui::IsKeyPressed(ImGuiKey_Backspace, false))) {
        es_planes_remove(state, target, sel.plane);
        target = nullptr;
    }

    // The outlines, and the one the mouse is on
    outline_count = 0;
    uint32_t hover_rep = 0;
    int hover_plane = -1;
    float hover_dist = 6.0f;
    md_float3_t centre = md_float3_zero();
    bool centre_valid = false;
    for (size_t r = 0; r < md_array_size(state->representation.reps); ++r) {
        const Representation& rep = state->representation.reps[r];
        if (rep.type != RepresentationType::ElectronicStructure || !representation_shown(*state, rep) || !rep.type_is_valid) continue;
        const ElectronicStructureRepresentation& es = rep.electronic_structure;
        if (es.display == ElectronicStructureDisplay::Isosurface) continue;
        md_float4x4_t box;
        if (!es_volume_box(&box, es)) continue;
        const bool slices = es.display == ElectronicStructureDisplay::Planes;
        const ElectronicStructurePlaneSet& set = electronic_structure_planes(es);
        for (int i = 0; i < set.count; ++i) {
            const bool selected = target == &rep && sel.plane == i;
            // a slice's rectangle, a cut's polygon where it cuts the box
            md_float3_t poly[6];
            int n = 0;
            if (slices) {
                if (!slice_has_rect(set.items[i])) continue;
                slice_corners(poly, set.items[i]);
                n = 4;
            } else {
                n = plane_box_polygon(poly, set.items[i].plane, box);
            }
            if (n == 0) continue;
            if (selected) {
                if (slices) {
                    centre = set.items[i].center;
                } else {
                    for (int k = 0; k < n; ++k) centre = md_float3_add(centre, poly[k]);
                    centre = md_float3_mul1(centre, 1.0f / (float)n);
                }
                centre_valid = true;
            }
            ImVec2 sp[6];
            md_float3_t ndc[6];
            bool all = true;
            for (int k = 0; k < n; ++k) all &= project(&sp[k], world_to_clip, poly[k], canvas_min, canvas_size) && project_ndc(&ndc[k], world_to_clip, poly[k]);
            if (!all) continue;
            if (!slices || selected) {
                const bool enabled = set.items[i].enabled;
                // a light line on a dark halo: seen over a light background and a dark one alike; hidden behind what
                // is in front of it (es_planes_draw_outlines)
                const ImU32 col = selected ? IM_COL32(255, 200, 40, 255) : enabled ? IM_COL32(255, 255, 255, 220) : IM_COL32(200, 200, 200, 130);
                const ImU32 halo = enabled || selected ? IM_COL32(0, 0, 0, 110) : IM_COL32(0, 0, 0, 60);
                const float width = selected ? 2.0f : 1.5f;
                for (int k = 0; k < n && outline_count < (int)MD_COUNTOF(outline_segs); ++k) {
                    outline_segs[outline_count++] = OutlineSegment{ndc[k], ndc[(k + 1) % n], col, halo, width};
                }
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
    int hover_su = 0, hover_sv = 0;
    if (target && (centre_valid || dragging)) {
        ElectronicStructurePlane& item = electronic_structure_planes(target->electronic_structure).items[sel.plane];
        Plane& plane = item.plane;
        const bool slice = target->electronic_structure.display == ElectronicStructureDisplay::Planes;
        const bool ours = dragging && gizmo_drag.rep_id == sel.rep_id && gizmo_drag.plane == sel.plane;
        if (dragging && !ours) gizmo_drag.part = GizmoPart_None;
        // A slice's normal is of no matter which way it points: it is turned to the camera, so that the arrow is
        // seen and taken hold of (not while it is dragged, which reads it)
        if (slice && !ours && !faces_camera(state, plane)) plane_flip(&plane);
        md_float3_t anchor = slice ? item.center : ours ? gizmo_drag.anchor : plane_project(plane, centre);
        if (!slice && !ours && gizmo_anchor.rep_id == sel.rep_id && gizmo_anchor.plane == sel.plane) {
            md_float4x4_t box;
            if (es_volume_box(&box, target->electronic_structure)) {
                const md_float3_t p = plane_project(plane, gizmo_anchor.point);
                const md_float4_t u = md_float4x4_mul_float4(md_float4x4_inverse(box), md_float4_set(p.x, p.y, p.z, 1.0f));
                if (u.x >= 0.0f && u.x <= 1.0f && u.y >= 0.0f && u.y <= 1.0f && u.z >= 0.0f && u.z <= 1.0f) anchor = p;
            }
        }
        const md_float3_t n = plane.normal;

        // its size on screen, and the normal's length in the world that makes the arrow as long as the ring's radius
        const float ring_radius = 64.0f;
        const md_float3_t right = md_quat_mul_float3(state->view.camera.orientation, md_float3_set(1, 0, 0));
        ImVec2 sc, sr;
        const bool visible = project(&sc, world_to_clip, anchor, canvas_min, canvas_size) &&
                             project(&sr, world_to_clip, md_float3_add(anchor, right), canvas_min, canvas_size);
        const float px_per_unit = visible ? point_distance(sr, sc) : 0.0f;
        if (visible && px_per_unit > 0.0f) {
            const float arrow_len = ring_radius / px_per_unit;
            ImVec2 tip;
            project(&tip, world_to_clip, md_float3_add(anchor, md_float3_mul1(n, arrow_len)), canvas_min, canvas_size);
            const float arrow_px = point_distance(tip, sc);
            const bool facing = faces_camera(state, plane);
            md_float3_t mo, mr;
            mouse_ray(&mo, &mr, world_to_clip, mouse, canvas_min, canvas_size);

            // A slice's handles, where its plane is not seen edge on (along the ray to its middle)
            SliceHandle handles[9];
            ImVec2 hs[9];
            int num_handles = 0;
            bool in_plane = false;
            if (slice && slice_has_rect(item)) {
                md_float3_t co, cr;
                mouse_ray(&co, &cr, world_to_clip, sc, canvas_min, canvas_size);
                in_plane = fabsf(md_float3_dot(n, cr)) > 0.15f;
                const md_float3_t u = md_float3_mul1(item.axis, item.extent.x);
                const md_float3_t v = md_float3_mul1(md_float3_cross(n, item.axis), item.extent.y);
                static const int sides[8][2] = {{-1,-1},{1,-1},{1,1},{-1,1},{1,0},{-1,0},{0,1},{0,-1}};
                bool all = in_plane;
                for (int k = 0; k < 8 && all; ++k) {
                    handles[k].world = md_float3_add(item.center, md_float3_add(md_float3_mul1(u, (float)sides[k][0]), md_float3_mul1(v, (float)sides[k][1])));
                    handles[k].su = sides[k][0];
                    handles[k].sv = sides[k][1];
                    all &= project(&hs[k], world_to_clip, handles[k].world, canvas_min, canvas_size);
                }
                if (all) {
                    num_handles = 8;
                    // the knob past the +u side's middle, out from the middle on screen
                    const ImVec2 d = hs[4] - sc;
                    const float len = point_distance(hs[4], sc);
                    if (len > 1.0f) {
                        hs[8] = hs[4] + d * (22.0f / len);
                        handles[8] = SliceHandle{handles[4].world, 0, 0};
                        num_handles = 9;
                    }
                }
            }

            // what the mouse is over: a slice's middle, its knob, its corners, the arrow's head, its sides; the arrow,
            // then the ring's disc
            const float md = point_distance(mouse, sc);
            const bool slide_by_ctrl = io.KeyCtrl;
            if (!dragging && canvas_hovered && !mouse_taken) {
                if (md <= 8.0f) {
                    hover_part = num_handles > 0 && !slide_by_ctrl ? GizmoPart_Move : GizmoPart_Slide;
                } else if (num_handles == 9 && point_distance(mouse, hs[8]) <= 7.0f) {
                    hover_part = GizmoPart_Twist;
                }
                for (int k = 0; k < 4 && num_handles > 0 && hover_part == GizmoPart_None; ++k) {
                    if (point_distance(mouse, hs[k]) <= 7.0f) {
                        hover_part = GizmoPart_Resize;
                        hover_su = handles[k].su;
                        hover_sv = handles[k].sv;
                    }
                }
                if (hover_part == GizmoPart_None && arrow_px > 12.0f && point_distance(mouse, tip) <= 9.0f) hover_part = GizmoPart_Slide;
                if (hover_part == GizmoPart_None && num_handles > 0) {
                    // the sides: corners 0-1 the -v side, 1-2 +u, 2-3 +v, 3-0 -u
                    static const int side_of[4][2] = {{0,-1},{1,0},{0,1},{-1,0}};
                    float best = 5.0f;
                    for (int k = 0; k < 4; ++k) {
                        const float d = segment_distance(mouse, hs[k], hs[(k + 1) % 4]);
                        if (d <= best) {
                            best = d;
                            hover_part = GizmoPart_Resize;
                            hover_su = side_of[k][0];
                            hover_sv = side_of[k][1];
                        }
                    }
                }
                if (hover_part == GizmoPart_None && arrow_px > 12.0f && segment_distance(mouse, sc, tip) <= 6.0f) hover_part = GizmoPart_Slide;
                if (hover_part == GizmoPart_None && md <= ring_radius + 4.0f) hover_part = GizmoPart_Rotate;
            }
            const GizmoPart lit = ours ? gizmo_drag.part : hover_part;
            const int lit_su = ours ? gizmo_drag.su : hover_su, lit_sv = ours ? gizmo_drag.sv : hover_sv;

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
                const float head = MD_MIN(12.0f, 0.4f * arrow_px);
                const ImVec2 base = ImVec2(tip.x - d.x * head, tip.y - d.y * head);
                dl->AddLine(sc, base, IM_COL32(0, 0, 0, 110), 5.0f);
                dl->AddLine(sc, base, arrow_col, 3.0f);
                dl->AddTriangleFilled(tip, ImVec2(base.x + o.x * head * 0.5f, base.y + o.y * head * 0.5f),
                                      ImVec2(base.x - o.x * head * 0.5f, base.y - o.y * head * 0.5f), arrow_col);
            }
            const ImU32 handle_col = IM_COL32(255, 200, 40, 255), handle_lit = IM_COL32(255, 240, 160, 255), handle_edge = IM_COL32(0, 0, 0, 150);
            auto square = [&](ImVec2 p, float h, ImU32 col) {
                dl->AddRectFilled(ImVec2(p.x - h - 1.0f, p.y - h - 1.0f), ImVec2(p.x + h + 1.0f, p.y + h + 1.0f), handle_edge);
                dl->AddRectFilled(ImVec2(p.x - h, p.y - h), ImVec2(p.x + h, p.y + h), col);
            };
            if (num_handles > 0) {
                for (int k = 0; k < 8; ++k) {
                    const bool on = lit == GizmoPart_Resize && handles[k].su == lit_su && handles[k].sv == lit_sv;
                    square(hs[k], on ? 5.0f : 3.5f, on ? handle_lit : handle_col);
                }
                if (num_handles == 9) {
                    const bool on = lit == GizmoPart_Twist;
                    dl->AddLine(hs[4], hs[8], handle_edge, 3.0f);
                    dl->AddLine(hs[4], hs[8], handle_col, 1.5f);
                    dl->AddCircleFilled(hs[8], on ? 7.0f : 6.0f, handle_edge, 16);
                    dl->AddCircleFilled(hs[8], on ? 6.0f : 5.0f, on ? handle_lit : handle_col, 16);
                }
                // the middle: moves it (with ctrl, slides it as the arrow does)
                const bool on = lit == GizmoPart_Move || (lit == GizmoPart_Slide && slide_by_ctrl);
                square(sc, on ? 6.0f : 5.0f, lit == GizmoPart_Slide ? arrow_col : on ? handle_lit : handle_col);
            } else {
                dl->AddCircleFilled(sc, 6.0f, IM_COL32(0, 0, 0, 120), 16);
                dl->AddCircleFilled(sc, 4.5f, arrow_col, 16);
            }

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
                float t = 0.0f;
                gizmo_drag.slide_t = line_ray_closest(&t, anchor, n, mo, mr) && arrow_px > 12.0f ? t : NAN;
                gizmo_drag.plane0 = plane;
                gizmo_drag.c0 = item.center;
                gizmo_drag.u0 = item.axis;
                gizmo_drag.v0 = md_float3_cross(n, item.axis);
                gizmo_drag.hu0 = item.extent.x;
                gizmo_drag.hv0 = item.extent.y;
                gizmo_drag.su = hover_su;
                gizmo_drag.sv = hover_sv;
                // where on the plane it was taken: the handle itself where the mouse's ray misses the plane (it does
                // not: the handles are not there then)
                if (!ray_plane(&gizmo_drag.grab, plane, mo, mr)) gizmo_drag.grab = anchor;
                const md_float3_t g = md_float3_sub(gizmo_drag.grab, item.center);
                gizmo_drag.twist0 = atan2f(md_float3_dot(g, gizmo_drag.v0), md_float3_dot(g, gizmo_drag.u0));
            }
            // a drag going on
            bool moved = false;
            if (gizmo_drag.part != GizmoPart_None && gizmo_drag.rep_id == sel.rep_id && gizmo_drag.plane == sel.plane) {
                if (!ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                    gizmo_drag.part = GizmoPart_None;
                    if (!slice) {
                        gizmo_anchor.rep_id = sel.rep_id;
                        gizmo_anchor.plane = sel.plane;
                        gizmo_anchor.point = gizmo_drag.anchor;
                    }
                } else if (gizmo_drag.part == GizmoPart_Rotate) {
                    const ImVec2 dm = ImVec2(mouse.x - gizmo_drag.last_mouse.x, mouse.y - gizmo_drag.last_mouse.y);
                    const float len = sqrtf(dm.x * dm.x + dm.y * dm.y);
                    const md_float3_t pivot = slice ? item.center : gizmo_drag.anchor;
                    if (len > 0.0f) {
                        // the front of a ball under the mouse moves with it: about (dy, dx, 0) in the view, by len / radius
                        const md_float3_t axis = md_float3_normalize(md_quat_mul_float3(state->view.camera.orientation, md_float3_set(dm.y, dm.x, 0.0f)));
                        const float angle = len / ring_radius;
                        plane = plane_from_point_normal(pivot, md_float3_normalize(rotate_about(plane.normal, axis, angle)));
                        // the rectangle turned with it
                        if (slice) item.axis = rotate_about(item.axis, axis, angle);
                    }
                    if (io.KeyShift) plane = plane_from_point_normal(pivot, snap_to_axis(plane.normal));
                    moved = true;
                } else if (gizmo_drag.part == GizmoPart_Slide) {
                    // along the normal it was taken with, from where it was taken
                    const md_float3_t sn = gizmo_drag.slide_n;
                    float t = 0.0f;
                    if (gizmo_drag.slide_t == gizmo_drag.slide_t && line_ray_closest(&t, gizmo_drag.slide_a0, sn, mo, mr)) {
                        t -= gizmo_drag.slide_t;
                    } else {
                        // seen along the normal: up the screen is along it
                        gizmo_drag.slide_off += (gizmo_drag.last_mouse.y - mouse.y) / px_per_unit;
                        t = gizmo_drag.slide_off;
                    }
                    gizmo_drag.anchor = md_float3_add(gizmo_drag.slide_a0, md_float3_mul1(sn, t));
                    plane.normal = sn;
                    plane.d = gizmo_drag.slide_d0 - t;
                    if (slice) item.center = md_float3_add(gizmo_drag.c0, md_float3_mul1(sn, t));
                    moved = true;
                } else {
                    // along the plane as it was taken: where the mouse's ray meets it, from where it was taken
                    md_float3_t p;
                    if (ray_plane(&p, gizmo_drag.plane0, mo, mr)) {
                        const md_float3_t d = md_float3_sub(p, gizmo_drag.grab);
                        const md_float3_t u0 = gizmo_drag.u0, v0 = gizmo_drag.v0;
                        if (gizmo_drag.part == GizmoPart_Move) {
                            item.center = md_float3_add(gizmo_drag.c0, d);
                        } else if (gizmo_drag.part == GizmoPart_Resize) {
                            // each side taken follows the mouse, the opposite one staying (with shift, both moving)
                            const float min_half = 0.05f;
                            const float ds = md_float3_dot(d, u0), dt = md_float3_dot(d, v0);
                            float hu = gizmo_drag.hu0, hv = gizmo_drag.hv0, ou = 0.0f, ov = 0.0f;
                            if (gizmo_drag.su != 0) {
                                const float grow = (float)gizmo_drag.su * ds;
                                hu = MD_MAX(min_half, io.KeyShift ? gizmo_drag.hu0 + grow : gizmo_drag.hu0 + 0.5f * grow);
                                ou = io.KeyShift ? 0.0f : (float)gizmo_drag.su * (hu - gizmo_drag.hu0);
                            }
                            if (gizmo_drag.sv != 0) {
                                const float grow = (float)gizmo_drag.sv * dt;
                                hv = MD_MAX(min_half, io.KeyShift ? gizmo_drag.hv0 + grow : gizmo_drag.hv0 + 0.5f * grow);
                                ov = io.KeyShift ? 0.0f : (float)gizmo_drag.sv * (hv - gizmo_drag.hv0);
                            }
                            item.extent = md_float2_set(hu, hv);
                            item.center = md_float3_add(gizmo_drag.c0, md_float3_add(md_float3_mul1(u0, ou), md_float3_mul1(v0, ov)));
                        } else if (gizmo_drag.part == GizmoPart_Twist) {
                            const md_float3_t g = md_float3_sub(p, gizmo_drag.c0);
                            float a = atan2f(md_float3_dot(g, v0), md_float3_dot(g, u0)) - gizmo_drag.twist0;
                            if (io.KeyShift) {
                                const float step = 15.0f * 3.14159265f / 180.0f;
                                a = roundf(a / step) * step;
                            }
                            item.axis = md_float3_normalize(rotate_about(u0, gizmo_drag.plane0.normal, a));
                        }
                        moved = true;
                    }
                }
                gizmo_drag.last_mouse = mouse;
                capture = true;
            }
            if (moved && slice) es_slices_conform(&target->electronic_structure);
            if (hover_part != GizmoPart_None) {
                capture = true;
                ImGuiMouseCursor cursor = ImGuiMouseCursor_Hand;
                if (hover_part == GizmoPart_Slide) cursor = ImGuiMouseCursor_ResizeNS;
                if (hover_part == GizmoPart_Move) cursor = ImGuiMouseCursor_ResizeAll;
                if (hover_part == GizmoPart_Resize) {
                    // across the side (or along the corner's diagonal) on screen
                    int k = 0;
                    while (k < 8 && !(handles[k].su == hover_su && handles[k].sv == hover_sv)) ++k;
                    const ImVec2 d = (k < 8 ? hs[k] : mouse) - sc;
                    if (hover_su != 0 && hover_sv != 0) cursor = d.x * d.y > 0.0f ? ImGuiMouseCursor_ResizeNWSE : ImGuiMouseCursor_ResizeNESW;
                    else cursor = fabsf(d.x) > fabsf(d.y) ? ImGuiMouseCursor_ResizeEW : ImGuiMouseCursor_ResizeNS;
                }
                ImGui::SetMouseCursor(cursor);
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
        gl::compile_shader_from_source(md_str_from_cstr(outline_vs), GL_VERTEX_SHADER),
        gl::compile_shader_from_source(md_str_from_cstr(outline_fs), GL_FRAGMENT_SHADER),
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
    static OutlineVertex verts[MD_COUNTOF(outline_segs) * 12];
    int nv = 0;
    const md_float2_t half = md_float2_set(0.5f * (float)fb_width, 0.5f * (float)fb_height);
    for (int pass = 0; pass < 2; ++pass) {
        for (int i = 0; i < outline_count; ++i) {
            const OutlineSegment& s = outline_segs[i];
            // in pixels: along it and across it, half the width out each way, and as far past each end (square caps
            // that close the corners)
            const md_float2_t a = md_float2_set(s.a.x * half.x, s.a.y * half.y), b = md_float2_set(s.b.x * half.x, s.b.y * half.y);
            md_float2_t d = md_float2_sub(b, a);
            const float len = md_float2_length(d);
            if (!(len > 1.0e-3f)) continue;
            d = md_float2_mul1(d, 1.0f / len);
            const float w = 0.5f * (pass == 0 ? s.width + 2.0f : s.width) * px_per_point;
            const float we = w + 1.0f;     // a pixel more, feathered (outline_fs)
            const md_float2_t o = md_float2_set(-d.y * we, d.x * we), e = md_float2_mul1(d, w);
            const md_float2_t pa = md_float2_sub(a, e), pb = md_float2_add(b, e);
            const md_float2_t q[4] = { md_float2_add(pa, o), md_float2_sub(pa, o), md_float2_sub(pb, o), md_float2_add(pb, o) };
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
