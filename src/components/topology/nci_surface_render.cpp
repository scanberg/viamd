#include "nci_surface_render.h"

#include <color_scale.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <gfx/gl_utils.h>
#include <gfx/isomesh_utils.h>
#include <implot.h>
#include <imgui.h>
#include <iso_meshes.h>
#include <viamd.h>

#include <array>
#include <cstdio>

namespace topology {
namespace {

uint32_t surface_colormap(NciState& state, NciSurface& n) {
    const ColorScale scale = nci_effective_scale(state, n);
    if (n.colormap && n.uploaded_colormap == scale.colormap) return n.colormap;
    uint32_t pixels[256];
    for (int i = 0; i < 256; ++i) {
        ImVec4 color = ImPlot::SampleColormap((float)i / 255.0f, scale.colormap);
        color.w = 1.0f;
        pixels[i] = ImGui::ColorConvertFloat4ToU32(color);
    }
    if ((!n.colormap && !gl::init_texture_2D(&n.colormap, 256, 1, GL_RGBA8)) ||
        !gl::set_texture_2D_data(n.colormap, 0, pixels, GL_RGBA8)) {
        MD_LOG_ERROR("NCI surface: could not upload the colour map");
        return 0;
    }
    n.uploaded_colormap = scale.colormap;
    return n.colormap;
}

void append_surface(NciState& state, MeshSurfaceRequest& request, NciSurface& n, bool& overflow) {
    if (!nci_surface_visible(state, n)) return;
    if (request.count >= request.cap) {
        overflow = true;
        return;
    }
    const uint32_t colormap = surface_colormap(state, n);
    if (!colormap) return;
    isomesh::Surface& surface = request.surfaces[request.count++];
    surface = {};
    surface.mesh = &n.mesh;
    surface.index_count = n.mesh.index_count;
    surface.model = iso_meshes_model_matrix();
    surface.color = vec4_set(1.0f, 1.0f, 1.0f, n.opacity);
    surface.coloring = isomesh::Coloring_VertexField;
    surface.field_colormap = colormap;
    const ColorScale scale = nci_effective_scale(state, n);
    surface.field_range_beg = scale.range_beg;
    surface.field_range_end = scale.range_end;
    request.hash = md_hash64(&scale.colormap, sizeof(scale.colormap), request.hash);
}

}  // namespace

void draw_nci_legends(const NciState& state, int first_slot) {
    if (!state.scale.show_legend) return;
    bool any = nci_surface_visible(state, state.primary);
    for (const auto& n : state.full_system) any |= nci_surface_visible(state, *n);
    for (const auto& n : state.previews) any |= nci_surface_visible(state, *n);
    for (const auto& n : state.objects) any |= nci_surface_visible(state, *n);
    if (!any) return;
    color_scale_draw_legend(state.scale, md_unit_pow(md_unit_bohr_radius(), -3), "sign(\xCE\xBB\xE2\x82\x82)\xCF\x81",
                            "NCI", -0x4E4349, first_slot, nullptr, true);
}

void nci_append_mesh_surfaces(NciState& state, MeshSurfaceRequest& request) {
    bool overflow = false;
    append_surface(state, request, state.primary, overflow);
    for (auto& n : state.previews) append_surface(state, request, *n, overflow);
    for (auto& n : state.full_system) append_surface(state, request, *n, overflow);
    for (auto& n : state.objects) append_surface(state, request, *n, overflow);
    if (overflow && !state.capacity_reported)
        MD_LOG_ERROR("NCI surfaces: renderer capacity exceeded; hide some surfaces to display the remaining objects");
    state.capacity_reported = overflow;
}

void nci_hash_render_state(const NciState& state, RenderKey& render_key) {
    const auto hash_surface = [&](const NciSurface& n) {
        const bool visible = nci_surface_visible(state, n);
        render_key.key = md_hash64(&visible, sizeof(visible), render_key.key);
        if (!visible) return;
        const ColorScale scale = nci_effective_scale(state, n);
        const float look[3] = {scale.range_beg, scale.range_end, n.opacity};
        render_key.key = md_hash64(look, sizeof(look), render_key.key);
        render_key.key = md_hash64(&n.id, sizeof(n.id), render_key.key);
        render_key.key = md_hash64(&scale.colormap, sizeof(scale.colormap), render_key.key);
        render_key.key = md_hash64(&n.mesh.version, sizeof(n.mesh.version), render_key.key);
        render_key.key = md_hash64(&n.want_key, sizeof(n.want_key), render_key.key);
        render_key.key = md_hash64(&n.have_key, sizeof(n.have_key), render_key.key);
    };
    hash_surface(state.primary);
    for (const auto& n : state.previews) hash_surface(*n);
    for (const auto& n : state.full_system) hash_surface(*n);
    for (const auto& n : state.objects) hash_surface(*n);
    render_key.key = md_hash64(&state.enabled, sizeof(state.enabled), render_key.key);
    render_key.key = md_hash64(&state.view, sizeof(state.view), render_key.key);
    render_key.key = md_hash64(&state.classified_overview, sizeof(state.classified_overview), render_key.key);
    render_key.key = md_hash64(&state.active_id, sizeof(state.active_id), render_key.key);
}

}  // namespace topology
