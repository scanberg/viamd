#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>

#include <viamd_event.h>
#include <viamd.h>
#include <serialization_utils.h>
#include <representation.h>
#include <property_grid.h>
#include <filter_input.h>

#include <gfx/gl_utils.h>
#include <gfx/volumerender_utils.h>
#include <gfx/immediate_draw_utils.h>
#include <color_utils.h>

#include <imgui_internal.h>
#include <imgui_widgets.h>

#include <implot_internal.h>
#include <implot_widgets.h>

#include <string.h>

constexpr uint64_t interaction_surface_density_vol = HASH_STR_LIT64("interaction surface density volume");

enum LegendColorMapMode_ {
    LegendColorMapMode_Opaque,
    LegendColorMapMode_Transparent,
    LegendColorMapMode_Split,
};

// Of an ensemble of reference structures, at most this many are drawn: each is a subset of the molecule per
// representation (md_gl_subset_t), drawn under its own transform
constexpr size_t DENSITY_VOLUME_MAX_SHOWN_STRUCTURES = 256;

// A representation of the reference structures. Its settings, name, filter and colours are those of a representation of
// the Representations window: the same Representation, edited and drawn by the same code (representation.h). Its
// filter is evaluated within each reference structure, which is its context as a group's filter is the context of what
// is in the group: 'all' is the structure, residue(1) its first residue, and 'out' reaches outside of it. What it picks
// within each structure shown is a subset of the molecule, drawn under the transform of that structure. A Cartoon draws
// the cartoon of the lipids and nucleic acids within it too: the glyphs whose molecule the subset holds the reference
// atom of (md_gl_subset_set_derived).
struct DensityVolumeRep {
    Representation rep = {};
    md_array(md_gl_subset_t) subsets = nullptr;       // [structures shown] What the filter picks within each
    md_array(md_gl_subset_t) glyph_subsets = nullptr; // [structures shown] Of a Cartoon: the glyphs of what each subset holds
    md_gl_rep_t members = {};                         // The colours of the structures after the first, faded (ensemble_fade)
    md_gl_rep_t members_glyphs = {};                  // And of their glyphs
    bool dirty = true;                            // What its filter picks within the structures to be found again, then its colours
    bool colors_dirty = true;                     // Its colours to be made again

    // What the subsets and colours were made for
    uint64_t built_version = 0;                   // DensityVolume::structures_version
    size_t   built_shown = 0;
    float    built_fade = 0.0f;
    bool     built_cartoon = false;               // With the glyphs of a Cartoon
};

// An empty filter, or 'all', picks all of a structure: nothing to evaluate
static bool filter_is_all(const char* text) {
    str_t s = str_trim(str_from_cstr(text));
    return str_empty(s) || str_eq(s, STR_LIT("all"));
}

struct DensityVolume : viamd::EventHandler {
    bool show_window = false;
    bool enabled = false;

    // The volume property shown, by path; empty for none. One at a time. Of an array of volumes, one of them: v[member + 1].
    SeriesKey volume_key = {};
    uint32_t  volume_member = 0;
    ApplicationState* app = nullptr;

    md_allocator_i* arena = nullptr;

    struct {
        bool enabled = true;
        struct {
            uint32_t id = 0;
            float alpha_scale = 1.f;
            int colormap = DEFAULT_COLORMAP;
            bool dirty = true;
            float min_val = 0.0f;
            float max_val = 1.0f;
        } tf;
    } dvr;

    struct {
        bool enabled = false;
        float values[8] = {};
        vec4_t colors[8] = {};
        size_t count = 0;
    } iso;

    struct {
        uint32_t id = 0;
        bool dirty = false;
        int  dim[3] = {0};
    } volume_texture;


    struct {
        vec3_t min = {0, 0, 0};
        vec3_t max = {1, 1, 1};
    } clip_volume;

    struct {
        bool enabled = true;
        bool checkerboard = true;
        int  colormap_mode = LegendColorMapMode_Split;
    } legend;

    vec3_t voxel_spacing = {1.0f, 1.0f, 1.0f};
    float resolution_scale = 2.0f;

    vec4_t clip_volume_color = {1,0,0,1};
    vec4_t bounding_box_color = {0,0,0,1};

    bool show_bounding_box = true;
    bool show_reference_structures = true;
    bool show_reference_ensemble = false;
    bool show_density_volume = false;
    bool show_coordinate_system_widget = true;

    // The settings beside the view, and how wide they are
    bool  show_panel = true;
    float panel_width = 320.0f;

    // The structures after the first of an ensemble, faded towards the background by this much: the first stands out
    float ensemble_fade = 0.5f;

    bool dirty_rep = false;     // The reference structures are to be evaluated again
    bool dirty_vol = false;
    bool view_pending = true;   // Another volume is shown: the view goes to it once its extent is known

    // How the volume shown stands against the script (series_state), what it is called ("v", "v[2]"), and how many
    // volumes its property has, as of the last update
    SeriesState volume_state = {};
    char volume_label[64] = "";
    uint32_t volume_count = 0;

    // The representations of the reference structures, in the order they are listed and drawn
    md_array(DensityVolumeRep) reps = nullptr;
    uint32_t next_rep_id = 1;

    // The reference structures as last evaluated: the atoms of each (kept, up to DENSITY_VOLUME_MAX_SHOWN_STRUCTURES),
    // where each is moved to align with the reference, and how many the script gives in all
    md_array(md_bitset_t) structures = nullptr;
    md_array(mat4_t) structure_mats = nullptr;
    size_t   structure_count = 0;
    uint64_t structures_hash = 0;
    uint64_t structures_version = 1;   // Changes with the atoms of the structures

    // Reading a workspace: how many representations it says there are (-1 for a workspace from before it said)
    int workspace_rep_count = -1;

    // Base model matrix to position the volume correctly in the world.
    mat4_t model_mat = mat4_ident();

    GBuffer gbuf = {};
    PickingSurface picking_surface = {};
    Camera camera = {};
	ViewTransform target = {};
    ViewTransform default_view = {};

    DensityVolume() {
        viamd::event_system_register_handler(*this);
    }

    // ## Representations of the reference structures

    DensityVolumeRep* add_rep(RepresentationType type = RepresentationType::BallAndStick, ColorMapping mapping = ColorMapping::Type) {
        DensityVolumeRep r = {};
        r.rep.type = type;
        r.rep.color_mapping = mapping;
        r.rep.node.id = next_rep_id++;
        snprintf(r.rep.name, sizeof(r.rep.name), "%s", "structure");
        md_array_push(reps, r, md_get_heap_allocator());
        return md_array_last(reps);
    }

    // What a representation holds on the GPU and of its evaluation, which belong to the molecule of the system: it
    // is made again from its settings when it is drawn next
    static void release_rep_resources(DensityVolumeRep& r) {
        for (size_t i = 0; i < md_array_size(r.subsets); ++i) {
            md_gl_subset_destroy(r.subsets[i]);
        }
        md_array_free(r.subsets, md_get_heap_allocator());
        r.subsets = nullptr;
        release_glyph_subsets(r);
        md_gl_rep_destroy(r.rep.md_rep);
        md_gl_rep_destroy(r.members);
        md_gl_rep_destroy(r.rep.glyph_rep);
        md_gl_rep_destroy(r.members_glyphs);
        r.rep.md_rep = {};
        r.members = {};
        r.rep.glyph_rep = {};
        r.members_glyphs = {};
        r.rep.glyphs_shown = false;
        md_bitset_free(&r.rep.atom_mask, md_get_heap_allocator());
        r.dirty = true;
        r.colors_dirty = true;
    }

    static void release_glyph_subsets(DensityVolumeRep& r) {
        for (size_t i = 0; i < md_array_size(r.glyph_subsets); ++i) {
            md_gl_subset_destroy(r.glyph_subsets[i]);
        }
        md_array_free(r.glyph_subsets, md_get_heap_allocator());
        r.glyph_subsets = nullptr;
    }

    void remove_rep(size_t idx) {
        ASSERT(idx < md_array_size(reps));
        release_rep_resources(reps[idx]);
        for (size_t i = idx; i + 1 < md_array_size(reps); ++i) {
            reps[i] = reps[i + 1];
        }
        md_array_pop(reps);
    }

    void clone_rep(size_t idx) {
        ASSERT(idx < md_array_size(reps));
        DensityVolumeRep copy = {};
        copy.rep = reps[idx].rep;
        // Its own GPU side and atoms, made when it is drawn
        copy.rep.md_rep = {};
        copy.rep.glyph_rep = {};
        copy.rep.atom_mask = {};
        copy.rep.node.id = next_rep_id++;
        md_array_push(reps, copy, md_get_heap_allocator());
        // Just after the original
        for (size_t i = md_array_size(reps) - 1; i > idx + 1; --i) {
            const DensityVolumeRep tmp = reps[i];
            reps[i] = reps[i - 1];
            reps[i - 1] = tmp;
        }
    }

    void remove_all_reps() {
        while (md_array_size(reps) > 0) {
            remove_rep(md_array_size(reps) - 1);
        }
    }

    void release_all_rep_resources() {
        for (size_t i = 0; i < md_array_size(reps); ++i) {
            release_rep_resources(reps[i]);
        }
    }

    // The atoms of the structures and where each goes. Those of the atoms are kept only when they changed: the
    // transforms change with every frame, the atoms seldom.
    void set_structures(const md_bitset_t* structs, const mat4_t* mats, size_t count) {
        md_allocator_i* heap = md_get_heap_allocator();
        structure_count = count;
        const size_t kept = MIN(count, DENSITY_VOLUME_MAX_SHOWN_STRUCTURES);

        md_array_resize(structure_mats, kept, heap);
        if (kept) MEMCPY(structure_mats, mats, kept * sizeof(mat4_t));

        uint64_t hash = md_hash64(&kept, sizeof(kept), 0);
        for (size_t i = 0; i < kept; ++i) {
            hash = md_bitset_hash64(structs[i], hash);
        }
        if (hash != structures_hash || md_array_size(structures) != kept) {
            for (size_t i = 0; i < md_array_size(structures); ++i) {
                md_bitset_free(&structures[i], heap);
            }
            md_array_resize(structures, kept, heap);
            for (size_t i = 0; i < kept; ++i) {
                structures[i] = md_bitset_copy(structs[i], heap);
            }
            structures_hash = hash;
            structures_version += 1;
        }
    }

    void clear_structures() {
        md_allocator_i* heap = md_get_heap_allocator();
        for (size_t i = 0; i < md_array_size(structures); ++i) {
            md_bitset_free(&structures[i], heap);
        }
        md_array_shrink(structures, 0);
        md_array_shrink(structure_mats, 0);
        structure_count = 0;
        structures_hash = 0;
        structures_version += 1;
    }

    // The structures drawn: the first, or all that are kept of an ensemble
    size_t shown_structures() const {
        const size_t n = md_array_size(structures);
        return show_reference_ensemble ? n : MIN(n, (size_t)1);
    }

    // The subsets of a representation: what its filter picks within each structure shown
    void build_subsets(ApplicationState* state, DensityVolumeRep& r, size_t shown) {
        md_allocator_i* heap = md_get_heap_allocator();
        md_allocator_i* frame = state->allocator.frame;
        md_temp_scope_t temp = md_temp_begin_in(frame);
        defer { md_temp_end(temp); };

        Representation& rep = r.rep;
        rep.type_is_valid = representation_atom_type_valid(*state, rep.type);

        // A subset per structure shown. Should md_gl run out of them, those that are left are not drawn: a subset of
        // id 0 would draw all of the molecule.
        while (md_array_size(r.subsets) < shown) {
            const md_gl_subset_t subset = md_gl_subset_create(state->mold.gl_mol);
            if (!subset.id) break;
            md_array_push(r.subsets, subset, heap);
        }
        while (md_array_size(r.subsets) > shown) {
            md_gl_subset_destroy(*md_array_last(r.subsets));
            md_array_pop(r.subsets);
        }

        // Within each structure what the filter picks, the structure its context. The first that fails says why.
        const bool all = filter_is_all(rep.filter.text);
        md_filter_result_t failure = {};
        bool failed = false;
        bool dynamic = false;
        md_bitset_t picked_all = md_bitset_empty_set();
        const size_t made = MIN(shown, md_array_size(r.subsets));
        for (size_t i = 0; i < made; ++i) {
            md_bitset_t picked = structures[i];
            if (!all) {
                md_filter_context_t ctx = filter_context(state);
                ctx.atoms = &structures[i];
                md_bitset_t out = md_bitset_empty_set();
                filter_evaluate(&rep.filter, &out, ctx, frame);
                dynamic |= rep.filter.result.is_dynamic;
                if (rep.filter.result.success) {
                    picked = out;
                } else {
                    picked = md_bitset_empty_set();
                    if (!failed) {
                        failure = rep.filter.result;
                        failed = true;
                    }
                }
            }
            const size_t count = (size_t)md_bitset_count(picked);
            uint32_t* indices = (uint32_t*)md_temp_alloc(temp, MAX(count, (size_t)1) * sizeof(uint32_t));
            md_bitset_extract_indices(indices, count, picked);
            md_gl_subset_set_atoms(r.subsets[i], indices, (uint32_t)count);
            picked_all = md_bitset_or(picked_all, picked, frame);
        }

        // Of a Cartoon, the glyphs of the lipids and nucleic acids of each subset
        const CartoonGlyphs& glyphs = state->mold.cartoon_glyphs;
        if (rep.type == RepresentationType::Cartoon && cartoon_glyphs_valid(glyphs)) {
            while (md_array_size(r.glyph_subsets) < made) {
                const md_gl_subset_t subset = md_gl_subset_create(glyphs.gl_mol);
                if (!subset.id) break;
                md_array_push(r.glyph_subsets, subset, heap);
            }
            while (md_array_size(r.glyph_subsets) > made) {
                md_gl_subset_destroy(*md_array_last(r.glyph_subsets));
                md_array_pop(r.glyph_subsets);
            }
            for (size_t i = 0; i < md_array_size(r.glyph_subsets); ++i) {
                md_gl_subset_set_derived(r.glyph_subsets[i], r.subsets[i]);
            }
        } else {
            release_glyph_subsets(r);
        }
        if (all) {
            rep.filter.result = {};
            rep.filter.result.success = true;
        } else if (made == 0) {
            // No structure to evaluate it within: whether it is valid at all, within the whole system
            md_bitset_t out = md_bitset_empty_set();
            filter_evaluate(&rep.filter, &out, filter_context(state), frame);
        } else if (failed) {
            rep.filter.result = failure;
            rep.filter.result.is_dynamic = dynamic;
        } else {
            rep.filter.result.is_dynamic = dynamic;
        }
        // What it shows in all: the span of the colour scale of an attribute
        md_bitset_assign(&rep.atom_mask, picked_all, heap);

        r.built_version = structures_version;
        r.built_shown = shown;
        r.built_cartoon = rep.type == RepresentationType::Cartoon;
        r.dirty = false;
        r.colors_dirty = true;
    }

    // The colours of a representation, and those of the structures after the first of an ensemble, faded
    void build_colors(ApplicationState* state, DensityVolumeRep& r, size_t shown, float fade) {
        md_allocator_i* frame = state->allocator.frame;
        md_temp_scope_t temp = md_temp_begin_in(frame);
        defer { md_temp_end(temp); };

        const size_t num_atoms = md_system_atom_count(&state->mold.sys);
        Representation& rep = r.rep;
        if (!rep.md_rep.id) rep.md_rep = md_gl_rep_create(state->mold.gl_mol);

        // Opaque: what is drawn is what the subsets hold
        uint32_t* colors = (uint32_t*)md_temp_alloc(temp, MAX(num_atoms, (size_t)1) * sizeof(uint32_t));
        representation_atom_colors(colors, num_atoms, state, &rep, &rep.atom_mask);
        for (size_t i = 0; i < num_atoms; ++i) {
            colors[i] |= 0xFF000000u;
        }
        md_gl_rep_set_atom_colors(rep.md_rep, 0, (uint32_t)num_atoms, colors, 0);

        // Of a Cartoon, its glyphs in the colours of their atoms
        const CartoonGlyphs& glyphs = state->mold.cartoon_glyphs;
        const bool with_glyphs = md_array_size(r.glyph_subsets) > 0;
        rep.glyphs_shown = false;
        if (with_glyphs) {
            if (!rep.glyph_rep.id) rep.glyph_rep = cartoon_glyphs_rep_create(glyphs);
            rep.glyphs_shown = cartoon_glyphs_rep_set_colors(glyphs, rep.glyph_rep, colors, num_atoms, rep.two_head_spheres);
        }

        // The other structures of an ensemble, towards the background
        if (shown > 1 && fade > 0.0f) {
            if (!r.members.id) r.members = md_gl_rep_create(state->mold.gl_mol);
            const vec3_t bg = state->visuals.background.color;
            const float t = CLAMP(fade, 0.0f, 1.0f);
            for (size_t i = 0; i < num_atoms; ++i) {
                vec4_t c = convert_color(colors[i]);
                c.x = c.x + (bg.x - c.x) * t;
                c.y = c.y + (bg.y - c.y) * t;
                c.z = c.z + (bg.z - c.z) * t;
                colors[i] = convert_color(c);
            }
            md_gl_rep_set_atom_colors(r.members, 0, (uint32_t)num_atoms, colors, 0);
            if (rep.glyphs_shown) {
                if (!r.members_glyphs.id) r.members_glyphs = cartoon_glyphs_rep_create(glyphs);
                cartoon_glyphs_rep_set_colors(glyphs, r.members_glyphs, colors, num_atoms, rep.two_head_spheres);
            }
        }

        r.built_fade = fade;
        r.colors_dirty = false;
    }

    // Those that need it, of those shown
    void update_reps(ApplicationState* state, bool frame_changed) {
        if (!state->mold.gl_mol.id || !show_reference_structures) return;
        const size_t shown = shown_structures();
        const float fade = shown > 1 ? ensemble_fade : 0.0f;
        for (size_t i = 0; i < md_array_size(reps); ++i) {
            DensityVolumeRep& r = reps[i];
            if (r.built_version != structures_version || r.built_shown != shown) r.dirty = true;
            if (r.built_cartoon != (r.rep.type == RepresentationType::Cartoon)) r.dirty = true;
            if (r.built_fade != fade) r.colors_dirty = true;
            if (frame_changed) {
                if (r.rep.filter.result.is_dynamic) r.dirty = true;
                if (r.rep.color_mapping == ColorMapping::SecondaryStructure || r.rep.color_mapping == ColorMapping::Attribute) r.colors_dirty = true;
            }
            if (!r.rep.enabled) continue;
            if (r.dirty) build_subsets(state, r, shown);
            if (r.colors_dirty) build_colors(state, r, shown, fade);
        }
    }

    // ## Workspace: which volume is shown, and how

    void reset_workspace_settings() {
        volume_key = {};
        volume_member = 0;
        dvr.enabled = true;
        dvr.tf.alpha_scale = 1.0f;
        dvr.tf.colormap = DEFAULT_COLORMAP;
        dvr.tf.min_val = 0.0f;
        dvr.tf.max_val = 1.0f;
        dvr.tf.dirty = true;
        iso.enabled = false;
        iso.count = 0;
        MEMSET(iso.values, 0, sizeof(iso.values));
        MEMSET(iso.colors, 0, sizeof(iso.colors));
        clip_volume.min = {0, 0, 0};
        clip_volume.max = {1, 1, 1};
        legend.enabled = true;
        legend.checkerboard = true;
        legend.colormap_mode = LegendColorMapMode_Split;
        resolution_scale = 2.0f;
        clip_volume_color = {1,0,0,1};
        bounding_box_color = {0,0,0,1};
        show_bounding_box = true;
        show_reference_structures = true;
        show_reference_ensemble = false;
        show_coordinate_system_widget = true;
        show_panel = true;
        panel_width = 320.0f;
        ensemble_fade = 0.5f;
        // The workspace says which, or one as there always was (deserialize_end)
        remove_all_reps();
        workspace_rep_count = -1;
        dirty_rep = true;
        dirty_vol = true;
    }

    void serialize(viamd::serialization_state_t& state) {
        viamd::write_section_header(state, STR_LIT("DensityVolume"));
        if (volume_key.path[0] != '\0') {
            viamd::write_str(state, STR_LIT("PropertySource"), str_from_cstr(series_source_name(volume_key.source)));
            viamd::write_str(state, STR_LIT("PropertyPath"), str_from_cstr(volume_key.path));
            if (volume_key.period != 0) {
                viamd::write_int(state, STR_LIT("PropertyPeriod"), (int)volume_key.period);
            }
            if (volume_member != 0) {
                viamd::write_int(state, STR_LIT("PropertyMember"), (int)volume_member);
            }
        }
        viamd::write_bool(state, STR_LIT("DvrEnabled"), dvr.enabled);
        viamd::write_int (state, STR_LIT("DvrColormap"), dvr.tf.colormap);
        viamd::write_flt (state, STR_LIT("DvrAlphaScale"), dvr.tf.alpha_scale);
        const float tf_range[2] = { dvr.tf.min_val, dvr.tf.max_val };
        viamd::write_flt_vec(state, STR_LIT("DvrRange"), tf_range, 2);
        viamd::write_bool(state, STR_LIT("IsoEnabled"), iso.enabled);
        viamd::write_int (state, STR_LIT("IsoCount"), (int64_t)iso.count);
        for (size_t i = 0; i < iso.count && i < ARRAY_SIZE(iso.values); ++i) {
            char key[32];
            snprintf(key, sizeof(key), "IsoValue%zu", i);
            viamd::write_flt(state, str_from_cstr(key), iso.values[i]);
            snprintf(key, sizeof(key), "IsoColor%zu", i);
            viamd::write_vec4(state, str_from_cstr(key), iso.colors[i]);
        }
        viamd::write_vec3(state, STR_LIT("ClipMin"), clip_volume.min);
        viamd::write_vec3(state, STR_LIT("ClipMax"), clip_volume.max);
        viamd::write_bool(state, STR_LIT("LegendEnabled"), legend.enabled);
        viamd::write_bool(state, STR_LIT("LegendCheckerboard"), legend.checkerboard);
        viamd::write_int (state, STR_LIT("LegendColormapMode"), legend.colormap_mode);
        viamd::write_flt (state, STR_LIT("ResolutionScale"), resolution_scale);
        viamd::write_vec4(state, STR_LIT("ClipVolumeColor"), clip_volume_color);
        viamd::write_vec4(state, STR_LIT("BoundingBoxColor"), bounding_box_color);
        viamd::write_bool(state, STR_LIT("ShowBoundingBox"), show_bounding_box);
        viamd::write_bool(state, STR_LIT("ShowReferenceStructures"), show_reference_structures);
        viamd::write_bool(state, STR_LIT("ShowReferenceEnsemble"), show_reference_ensemble);
        viamd::write_flt (state, STR_LIT("EnsembleFade"), ensemble_fade);
        viamd::write_bool(state, STR_LIT("ShowCoordinateSystem"), show_coordinate_system_widget);
        viamd::write_bool(state, STR_LIT("ShowPanel"), show_panel);
        viamd::write_flt (state, STR_LIT("PanelWidth"), panel_width);
        viamd::write_int (state, STR_LIT("RepresentationCount"), (int64_t)md_array_size(reps));

        // Each representation of the reference structures in a section of its own after this one, with the entries
        // of a representation of the Representations window
        for (size_t i = 0; i < md_array_size(reps); ++i) {
            viamd::write_section_header(state, STR_LIT("DensityVolumeRepresentation"));
            representation_write_entries(state, app, reps[i].rep);
        }
    }

    // The representation the entries of a workspace from before there were several are of
    DensityVolumeRep* legacy_rep() {
        if (md_array_size(reps) == 0) add_rep();
        return reps;
    }

    void deserialize(viamd::deserialization_state_t& state) {
        SeriesSource source = SeriesSource_Script;
        uint32_t period = 0;
        uint32_t member = 0;
        str_t path = {};
        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            if (str_eq_cstr(ident, "PropertySource")) {
                str_t name;
                viamd::extract_str(name, arg);
                series_source_from_name(&source, name, &period);
            } else if (str_eq_cstr(ident, "PropertyPeriod")) {
                int i = 0;
                if (viamd::extract_int(i, arg)) period = (uint32_t)MAX(i, 0);
            } else if (str_eq_cstr(ident, "PropertyMember")) {
                int i = 0;
                if (viamd::extract_int(i, arg)) member = (uint32_t)MAX(i, 0);
            } else if (str_eq_cstr(ident, "PropertyPath")) {
                viamd::extract_str(path, arg);
            } else if (str_eq_cstr(ident, "DvrEnabled")) {
                viamd::extract_bool(dvr.enabled, arg);
            } else if (str_eq_cstr(ident, "DvrColormap")) {
                viamd::extract_int(dvr.tf.colormap, arg);
            } else if (str_eq_cstr(ident, "DvrAlphaScale")) {
                viamd::extract_flt(dvr.tf.alpha_scale, arg);
            } else if (str_eq_cstr(ident, "DvrRange")) {
                float r[2];
                if (viamd::extract_flt_vec(r, 2, arg)) { dvr.tf.min_val = r[0]; dvr.tf.max_val = r[1]; }
            } else if (str_eq_cstr(ident, "IsoEnabled")) {
                viamd::extract_bool(iso.enabled, arg);
            } else if (str_eq_cstr(ident, "IsoCount")) {
                int n;
                if (viamd::extract_int(n, arg)) iso.count = (size_t)CLAMP(n, 0, (int)ARRAY_SIZE(iso.values));
            } else if (str_eq_cstr(ident, "ClipMin")) {
                viamd::extract_vec3(clip_volume.min, arg);
            } else if (str_eq_cstr(ident, "ClipMax")) {
                viamd::extract_vec3(clip_volume.max, arg);
            } else if (str_eq_cstr(ident, "LegendEnabled")) {
                viamd::extract_bool(legend.enabled, arg);
            } else if (str_eq_cstr(ident, "LegendCheckerboard")) {
                viamd::extract_bool(legend.checkerboard, arg);
            } else if (str_eq_cstr(ident, "LegendColormapMode")) {
                viamd::extract_int(legend.colormap_mode, arg);
            } else if (str_eq_cstr(ident, "ResolutionScale")) {
                viamd::extract_flt(resolution_scale, arg);
            } else if (str_eq_cstr(ident, "ClipVolumeColor")) {
                viamd::extract_vec4(clip_volume_color, arg);
            } else if (str_eq_cstr(ident, "BoundingBoxColor")) {
                viamd::extract_vec4(bounding_box_color, arg);
            } else if (str_eq_cstr(ident, "ShowBoundingBox")) {
                viamd::extract_bool(show_bounding_box, arg);
            } else if (str_eq_cstr(ident, "ShowReferenceStructures")) {
                viamd::extract_bool(show_reference_structures, arg);
            } else if (str_eq_cstr(ident, "ShowReferenceEnsemble")) {
                viamd::extract_bool(show_reference_ensemble, arg);
            } else if (str_eq_cstr(ident, "EnsembleFade")) {
                viamd::extract_flt(ensemble_fade, arg);
            } else if (str_eq_cstr(ident, "ShowCoordinateSystem")) {
                viamd::extract_bool(show_coordinate_system_widget, arg);
            } else if (str_eq_cstr(ident, "ShowPanel")) {
                viamd::extract_bool(show_panel, arg);
            } else if (str_eq_cstr(ident, "PanelWidth")) {
                viamd::extract_flt(panel_width, arg);
            } else if (str_eq_cstr(ident, "RepresentationCount")) {
                int n = 0;
                if (viamd::extract_int(n, arg)) workspace_rep_count = MAX(n, 0);
            } else if (str_eq_cstr(ident, "RepType")) {
                // From before there were several representations: the one there was
                viamd::extract_enum(legacy_rep()->rep.type, arg, (int)RepresentationType::Cartoon + 1);
            } else if (str_eq_cstr(ident, "RepColorMapping")) {
                viamd::extract_enum(legacy_rep()->rep.color_mapping, arg, (int)ColorMapping::Count);
            } else if (str_eq_cstr(ident, "RepParam")) {
                // Its parameters were the scales, in the order of a representation's
                viamd::extract_flt_vec(legacy_rep()->rep.scale.elem, 4, arg);
            } else if (str_eq_cstr(ident, "RepColor")) {
                viamd::extract_vec4(legacy_rep()->rep.base_color, arg);
            } else {
                for (size_t i = 0; i < ARRAY_SIZE(iso.values); ++i) {
                    char key[32];
                    snprintf(key, sizeof(key), "IsoValue%zu", i);
                    if (str_eq_cstr(ident, key)) { viamd::extract_flt(iso.values[i], arg); break; }
                    snprintf(key, sizeof(key), "IsoColor%zu", i);
                    if (str_eq_cstr(ident, key)) { viamd::extract_vec4(iso.colors[i], arg); break; }
                }
            }
        }
        if (!str_empty(path)) {
            volume_key = series_key(source, path, SeriesVariant_Values, period);
            volume_member = member;
        }
        // Workspaces from before the two were exclusive can have both on: DVR wins, as it did by default
        if (dvr.enabled && iso.enabled) {
            iso.enabled = false;
        }
        ensemble_fade = CLAMP(ensemble_fade, 0.0f, 1.0f);
        panel_width = MAX(panel_width, 100.0f);
        dvr.tf.dirty = true;
        dirty_rep = true;
        dirty_vol = true;
    }

    void deserialize_rep(viamd::deserialization_state_t& state) {
        DensityVolumeRep* r = add_rep();
        str_t ident, arg;
        while (viamd::next_entry(ident, arg, state)) {
            representation_read_entry(&r->rep, ident, arg);
        }
        // Only representations of atoms are drawn here
        if (r->rep.type > RepresentationType::Cartoon) {
            r->rep.type = RepresentationType::BallAndStick;
        }
    }

    void deserialize_end() {
        // A workspace from before there were several, or from before there was a choice: one as there always was
        if (workspace_rep_count < 0 && md_array_size(reps) == 0) {
            add_rep();
        }
    }

    void update(ApplicationState* state) {
        md_allocator_i* temp_arena = state->allocator.frame;
        md_temp_scope_t temp_scope = md_temp_begin_in(temp_arena);
        defer { md_temp_end(temp_scope); };

        if (dvr.tf.dirty) {
            dvr.tf.dirty = false;
            // Update colormap texture
            volume::compute_transfer_function_texture_simple(&dvr.tf.id, dvr.tf.colormap, dvr.tf.alpha_scale);
        }

        // Resolved by path every frame: the evaluation's table is rebuilt when the script recompiles
        const bool selected = volume_key.path[0] != '\0';
        SeriesVolumeView vol_view = {};
        const bool resolved = selected && series_resolve_volume(&vol_view, state, volume_key, volume_member);
        if (resolved) {
            // Of an array which became shorter, the last of it
            volume_member = vol_view.member;
        }

        const float* vol_data = resolved ? vol_view.data : NULL;
        const uint32_t* vol_dim = vol_view.dim;
        volume_state = resolved ? vol_view.state : SeriesState{};
        volume_count = resolved ? vol_view.count : 0;
        str_copy_to_char_buf(volume_label, sizeof(volume_label), str_from_cstr(vol_view.label));
        const md_script_vis_ref_t vis_ref = vol_view.vis_ref;
        const uint64_t data_version = vol_view.version;

        static SeriesKey s_volume_key = {};
        static bool s_first = true;
        if (s_first || !series_key_equal(s_volume_key, volume_key)) {
            // Another volume, or the first: looked at whole when its extent is known, wherever the view was
            if (selected) {
                view_pending = true;
            }
            s_first = false;
            s_volume_key = volume_key;
            dirty_vol = true;
            dirty_rep = true;
        }
        show_density_volume = selected;

        static uint64_t s_script_fingerprint = 0;
        if (s_script_fingerprint != md_script_ir_fingerprint(state->script.eval_ir)) {
            s_script_fingerprint = md_script_ir_fingerprint(state->script.eval_ir);
            dirty_vol = true;
            dirty_rep = true;
        }

        static uint64_t s_data_version = 0;
        if (s_data_version != data_version) {
            s_data_version = data_version;
            dirty_vol = true;
        }

        // Another volume of the array: other data, and other reference structures
        static uint32_t s_member = 0;
        if (s_member != volume_member) {
            s_member = volume_member;
            dirty_vol = true;
            dirty_rep = true;
        }

        static double s_frame = 0;
        bool frame_changed = false;
        if (s_frame != state->animation.frame) {
            s_frame = state->animation.frame;
            dirty_rep = true;
            frame_changed = true;
        }

        if (dirty_rep) {
            if (vol_data && !md_script_vis_ref_empty(vis_ref)) {
                dirty_rep = false;
                bool result = false;
                md_script_vis_t vis = {};

                if (md_script_ir_valid(state->script.eval_ir)) {
                    md_script_vis_init(&vis, state->allocator.frame);
                    md_script_vis_ctx_t ctx = {
                        .ir = state->script.eval_ir,
                        .sys = &state->mold.sys,
                        .state = &state->mold.state
                    };
                    // The structures of the member shown: subidx i is the i-th volume of an array, 0 the one volume
                    result = md_script_vis_eval_ref(&vis, vis_ref, (int)vol_view.member, &ctx, MD_SCRIPT_VISUALIZE_SDF);
                }

                if (result) {
                    if (vis.sdf.extent) {
                        const float s = vis.sdf.extent;
                        vec3_t min_aabb = vec3_set1(-s);
                        vec3_t max_aabb = vec3_set1( s);
                        model_mat = volume::compute_model_to_world_matrix(min_aabb, max_aabb);
                        voxel_spacing = vec3_t{2*s / vol_dim[0], 2*s / vol_dim[1], 2*s / vol_dim[2]};
                        default_view = compute_optimal_view((min_aabb + max_aabb) * 0.5f, (max_aabb - min_aabb) * 0.5f);
                        if (view_pending) {
                            view_pending = false;
                            target = default_view;
                            camera = default_view;
                        }
                    }
                    const size_t count = MIN(md_array_size(vis.sdf.structures), md_array_size(vis.sdf.matrices));
                    set_structures(vis.sdf.structures, vis.sdf.matrices, count);
                } else {
                    clear_structures();
                }
            }
        }

        update_reps(state, frame_changed);

        if (dirty_vol) {
            if (vol_data) {
                dirty_vol = false;
                if (!volume_texture.id) {
                    int dim[3] = { (int)vol_dim[0], (int)vol_dim[1], (int)vol_dim[2] };
                    gl::init_texture_3D(&volume_texture.id, dim[0], dim[1], dim[2], GL_R16F);
                    MEMCPY(volume_texture.dim, dim, sizeof(dim));
                }
                gl::set_texture_3D_data(volume_texture.id, 0, vol_data, GL_R32F);
                volume::notify_data_changed(volume_texture.id);
            }
        }
    }

    // ## Settings

    // The volume properties of the script, one to pick: picking one replaces the other, picking it again hides it
    void draw_volume_choices(ApplicationState* state) {
        int candidate_count = 0;
        auto list = [&](uint32_t period) {
            series_for_each_script_property(state, SeriesSource_Script, MD_SCRIPT_PROPERTY_FLAG_VOLUME, [&](const SeriesKey& over_all) {
                SeriesKey key = over_all;
                key.period = period;
                char label[96];
                series_label(label, sizeof(label), state, key);
                const uint32_t count = series_volume_count(state, key);
                if (count > 1) {
                    const size_t len = strnlen(label, sizeof(label));
                    snprintf(label + len, sizeof(label) - len, " \xC2\xB7 %u volumes", count);
                }
                const bool is_selected = series_key_equal(key, volume_key);
                ImGui::PushID((int)period);
                ImPlot::ItemIcon(series_default_color(state, key)); ImGui::SameLine();
                if (ImGui::Selectable(label, is_selected)) {
                    volume_key = is_selected ? SeriesKey{} : key;
                    volume_member = 0;
                }
                if (ImGui::IsItemHovered()) {
                    const str_t ident = series_script_ident(key);
                    script_visualize_ref(state, md_script_ir_property_vis_ref(state->script.eval_ir, ident), -1, MD_SCRIPT_VISUALIZE_DEFAULT);
                    script_set_hovered_property(state, ident);
                }
                ImGui::PopID();
                candidate_count += 1;
            });
        };
        // Over every frame, and over each timeline period: "v (Production)"
        list(0);
        for (int i = 0; i < state->timeline.periods.count && candidate_count > 0; ++i) {
            if (i == 0) ImGui::SeparatorText("Over a period");
            list(state->timeline.periods.items[i].id);
        }
        if (candidate_count == 0) {
            ImGui::TextDisabled("No volume properties in the script.");
        }
    }

    // Of an array of volumes ('sdf(...) in residue(1:n)'), the one shown, as the script names them
    void draw_member_slider(ApplicationState* state, float width) {
        const str_t ident = series_script_ident(volume_key);
        int m = (int)volume_member + 1;
        char fmt[96];
        snprintf(fmt, sizeof(fmt), STR_FMT "[%%d] of %u", STR_ARG(ident), volume_count);
        ImGui::SetNextItemWidth(width);
        if (ImGui::SliderInt("##member", &m, 1, (int)volume_count, fmt, ImGuiSliderFlags_AlwaysClamp)) {
            volume_member = (uint32_t)(m - 1);
        }
        if (ImGui::IsItemHovered() || ImGui::IsItemActive()) {
            // Whose volume it is: the reference structures of that one
            script_visualize_ref(state, md_script_ir_property_vis_ref(state->script.eval_ir, ident), (int)volume_member, MD_SCRIPT_VISUALIZE_DEFAULT);
            script_set_hovered_property(state, ident);
        }
    }

    void draw_section_volume(ApplicationState* state) {
        if (!prop::section("Volume")) return;
        char preview[96] = "none";
        if (volume_key.path[0] != '\0') {
            series_label(preview, sizeof(preview), state, volume_key);
        }
        prop::row("property", "A volume of the script (sdf(...) and the like), which can also be dragged here from its mark in the script editor");
        if (ImGui::BeginCombo("##property", preview, ImGuiComboFlags_HeightLarge)) {
            draw_volume_choices(state);
            ImGui::EndCombo();
        }
        if (volume_count > 1) {
            prop::row("volume", "Of the array of volumes, the one shown");
            draw_member_slider(state, -FLT_MIN);
        }
        if (volume_key.path[0] != '\0' && volume_state.status != SeriesStatus_Current) {
            char tip[256];
            series_state_describe(tip, sizeof(tip), volume_state);
            if (tip[0]) prop::hint("%s", tip);
        }
        prop::section_end();
    }

    void draw_section_rendering() {
        if (!prop::section("Rendering")) return;
        // One or the other (or neither): the two are separate renderers and never mixed
        const char* modes[] = {"Off", "Direct volume rendering", "Isosurfaces"};
        int mode = dvr.enabled ? 1 : iso.enabled ? 2 : 0;
        prop::row("mode");
        if (ImGui::Combo("##mode", &mode, modes, IM_ARRAYSIZE(modes))) {
            dvr.enabled = mode == 1;
            iso.enabled = mode == 2;
        }
        if (dvr.enabled) {
            prop::row("colormap");
            const ImVec2 size = {prop::control_width(), ImGui::GetFrameHeight()};
            if (ImPlot::ColormapButton(ImPlot::GetColormapName(dvr.tf.colormap), size, dvr.tf.colormap)) {
                ImGui::OpenPopup("Colormap Selector");
            }
            if (ImGui::BeginPopup("Colormap Selector")) {
                for (int map = 4; map < ImPlot::GetColormapCount(); ++map) {
                    if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), size, map)) {
                        dvr.tf.colormap = map;
                        dvr.tf.dirty = true;
                        ImGui::CloseCurrentPopup();
                    }
                }
                ImGui::EndPopup();
            }
            prop::row("opacity", "Scale of the opacity of the transfer function");
            if (ImGui::SliderFloat("##alpha_scale", &dvr.tf.alpha_scale, 0.001f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic)) {
                dvr.tf.dirty = true;
            }
            const float tf_min = 0.0f;
            const float tf_max = 1000.0f;
            prop::row("min value", "The value at the start of the colormap");
            ImGui::SliderFloat("##tf_min", &dvr.tf.min_val, tf_min, dvr.tf.max_val, "%.3f", ImGuiSliderFlags_Logarithmic);
            prop::row("max value", "The value at the end of the colormap");
            ImGui::SliderFloat("##tf_max", &dvr.tf.max_val, dvr.tf.min_val, tf_max, "%.3f", ImGuiSliderFlags_Logarithmic);
        }
        if (iso.enabled) {
            // Each isosurface: its value and colour
            const float gap = ImGui::GetStyle().ItemInnerSpacing.x;
            const float btn = ImGui::GetFrameHeight();
            int remove = -1;
            for (int i = 0; i < (int)iso.count; ++i) {
                ImGui::PushID(i);
                prop::row(i == 0 ? "isovalues" : "");
                ImGui::SetNextItemWidth(MAX(prop::control_width() - 2.0f * (btn + gap), 1.0f));
                ImGui::SliderFloat("##isovalue", &iso.values[i], 0.0f, 10.f, "%.3f", ImGuiSliderFlags_Logarithmic);
                ImGui::SameLine(0.0f, gap);
                ImGui::ColorEdit4Minimal("##color", iso.colors[i].elem);
                ImGui::SameLine(0.0f, gap);
                if (ImGui::DeleteButton(ICON_FA_XMARK, ImVec2(btn, btn))) {
                    remove = i;
                }
                ImGui::SetItemTooltip("Remove the isosurface");
                ImGui::PopID();
            }
            if (remove >= 0) {
                for (int j = remove; j < (int)iso.count - 1; ++j) {
                    iso.colors[j] = iso.colors[j+1];
                    iso.values[j] = iso.values[j+1];
                }
                iso.count -= 1;
            }
            prop::row(iso.count == 0 ? "isovalues" : "");
            ImGui::BeginDisabled(iso.count >= ARRAY_SIZE(iso.values));
            if (ImGui::Button(ICON_FA_PLUS " add")) {
                size_t idx = iso.count++;
                iso.values[idx] = 0.1f;
                iso.colors[idx] = { 0.2f, 0.1f, 0.9f, 1.0f };
            }
            ImGui::EndDisabled();
            ImGui::SameLine();
            ImGui::BeginDisabled(iso.count == 0);
            if (ImGui::Button("clear")) {
                iso.count = 0;
            }
            ImGui::EndDisabled();
        }
        prop::section_end();
    }

    void draw_section_clipping() {
        if (!prop::section("Clipping", false)) return;
        // The part of the volume drawn, as fractions of its box along each axis
        const char* axes[3] = {"x", "y", "z"};
        for (int k = 0; k < 3; ++k) {
            ImGui::PushID(k);
            prop::row(axes[k]);
            ImGui::RangeSliderFloat("##range", &clip_volume.min.elem[k], &clip_volume.max.elem[k], 0.0f, 1.0f, "(%.2f, %.2f)");
            ImGui::PopID();
        }
        prop::row("");
        const bool clipped = clip_volume.min.x > 0 || clip_volume.min.y > 0 || clip_volume.min.z > 0 ||
                             clip_volume.max.x < 1 || clip_volume.max.y < 1 || clip_volume.max.z < 1;
        ImGui::BeginDisabled(!clipped);
        if (ImGui::Button("reset")) {
            clip_volume.min = {0, 0, 0};
            clip_volume.max = {1, 1, 1};
        }
        ImGui::EndDisabled();
        prop::section_end();
    }

    void draw_section_overlays() {
        if (!prop::section("Overlays", false)) return;
        const float gap = ImGui::GetStyle().ItemInnerSpacing.x;
        prop::row("bounding box", "The box of the volume, and of the part of it drawn (Clipping)");
        ImGui::Checkbox("##bounding_box", &show_bounding_box);
        if (show_bounding_box) {
            ImGui::SameLine(0.0f, gap);
            ImGui::ColorEdit4("##bounding_box_color", bounding_box_color.elem, ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_NoLabel);
        }
        prop::row("legend", "The colormap of the direct volume rendering, in the corner of the view");
        ImGui::Checkbox("##legend", &legend.enabled);
        if (legend.enabled) {
            const char* colormap_modes[] = {"Opaque", "Transparent", "Split"};
            prop::row("legend colors", "The colormap as it is, as it is with the opacity of the transfer function, or both");
            ImGui::Combo("##legend_mode", &legend.colormap_mode, colormap_modes, IM_ARRAYSIZE(colormap_modes));
            prop::row("checkerboard", "A checkerboard behind the transparent part of the legend");
            ImGui::Checkbox("##checkerboard", &legend.checkerboard);
        }
        prop::row("axes", "The axes of the view, in the corner; a click on one looks along it");
        ImGui::Checkbox("##axes", &show_coordinate_system_widget);
        prop::section_end();
    }

    // A representation of the reference structures: its row (the arrow, its name, its buttons), and below it when it
    // is open its properties, as in the Representations window
    void draw_rep(ApplicationState* state, size_t idx, size_t* remove, size_t* clone) {
        DensityVolumeRep& r = reps[idx];
        Representation& rep = r.rep;
        ImGui::PushID((int)rep.node.id);

        const float btn = ImGui::GetFrameHeight();
        const float gap = ImGui::GetStyle().ItemInnerSpacing.x;
        const float right_x = ImGui::GetCursorPosX() + ImGui::GetContentRegionAvail().x;

        const bool open = ImGui::TreeNodeEx("##node", ImGuiTreeNodeFlags_FramePadding | ImGuiTreeNodeFlags_DefaultOpen | ImGuiTreeNodeFlags_AllowOverlap);
        ImGui::SameLine();
        const float buttons = 3.0f * btn + 3.0f * gap;
        ImGui::SetNextItemWidth(MAX(right_x - ImGui::GetCursorPosX() - buttons, 1.0f));
        const bool invalid = !rep.filter.result.success && !r.dirty && rep.filter.text[0] != '\0';
        if (invalid) ImGui::PushInvalid();
        ImGui::InputText("##name", rep.name, sizeof(rep.name));
        if (invalid) {
            ImGui::PopInvalid();
            ImGui::SetItemTooltip("%s", rep.filter.result.error);
        }
        ImGui::SameLine(right_x - buttons + gap, 0.0f);
        if (ImGui::Button(rep.enabled ? ICON_FA_EYE : ICON_FA_EYE_SLASH, ImVec2(btn, btn))) {
            rep.enabled = !rep.enabled;
        }
        ImGui::SetItemTooltip("Show/Hide");
        ImGui::SameLine(0.0f, gap);
        if (ImGui::Button(ICON_FA_COPY, ImVec2(btn, btn))) {
            *clone = idx;
        }
        ImGui::SetItemTooltip("Duplicate");
        ImGui::SameLine(0.0f, gap);
        if (ImGui::DeleteButton(ICON_FA_XMARK, ImVec2(btn, btn))) {
            *remove = idx;
        }
        ImGui::SetItemTooltip("Remove");

        if (open) {
            // Its type only says what is drawn of the subsets, its geometry how: neither is found again, but for the
            // glyphs of a Cartoon (update_reps)
            if (draw_representation_type(state, rep, true)) {
                rep.type_is_valid = representation_atom_type_valid(*state, rep.type);
            }

            // Completed within the first structure, as each is evaluated within its own
            md_filter_context_t ctx = filter_context(state);
            if (md_array_size(structures) > 0) ctx.atoms = &structures[0];
            prop::row("filter", "Evaluated within each reference structure: all is the structure, residue(1) its first residue.\n"
                                "'out' reaches outside of it, as in a script: within(5, out resname(\"LIG\")).");
            if (filter_input("##filter", &rep.filter, ctx, 0, nullptr, "all: the whole structure")) {
                r.dirty = true;
            }

            if (prop::section("Color")) {
                if (draw_representation_coloring(state, rep)) {
                    r.colors_dirty = true;
                }
                prop::section_end();
            }
            if (prop::section("Geometry", false)) {
                // Only what the args of the draws take, but whether the head of a lipid is two spheres, which its glyphs'
                // colours do
                const bool two_heads = rep.two_head_spheres;
                draw_representation_geometry(state, rep);
                if (rep.two_head_spheres != two_heads) r.colors_dirty = true;
                prop::section_end();
            }
            ImGui::Spacing();
            ImGui::TreePop();
        }
        ImGui::PopID();
    }

    void draw_section_structures(ApplicationState* state) {
        if (!prop::section_toggle("Reference structures", &show_reference_structures, nullptr, true,
                                  "The structures the volume is accumulated around, aligned with it")) {
            return;
        }
        const size_t shown = shown_structures();
        prop::row("superimposed", "Every structure the volume was accumulated over, each where it aligns with the reference, rather than the first only");
        ImGui::Checkbox("##ensemble", &show_reference_ensemble);
        if (structure_count > 0) {
            ImGui::SameLine();
            if (show_reference_ensemble && structure_count > shown) {
                ImGui::TextDisabled("%zu of %zu shown", shown, structure_count);
            } else {
                ImGui::TextDisabled("%zu structure%s", structure_count, structure_count == 1 ? "" : "s");
            }
        }
        if (show_reference_ensemble && structure_count > 1) {
            prop::row("fade others", "How far the structures after the first go towards the background, which sets the first apart");
            ImGui::SliderFloat("##fade", &ensemble_fade, 0.0f, 1.0f, "%.2f", ImGuiSliderFlags_AlwaysClamp);
        }

        size_t remove = SIZE_MAX;
        size_t clone = SIZE_MAX;
        for (size_t i = 0; i < md_array_size(reps); ++i) {
            draw_rep(state, i, &remove, &clone);
        }
        if (remove != SIZE_MAX) {
            remove_rep(remove);
        } else if (clone != SIZE_MAX) {
            clone_rep(clone);
        }
        if (ImGui::Button(ICON_FA_PLUS " representation")) {
            add_rep();
        }
        ImGui::SetItemTooltip("Another representation of the reference structures, a part of them picked by its filter");
        prop::section_end();
    }

    void draw_panel(ApplicationState* state) {
        prop::begin("density_volume");
        draw_section_volume(state);
        draw_section_rendering();
        draw_section_structures(state);
        draw_section_clipping();
        draw_section_overlays();
        prop::end();
    }

    void draw(ApplicationState* state) {
        if (!show_window) return;

        md_allocator_i* temp_arena = state->allocator.frame;
        md_temp_scope_t temp_scope = md_temp_begin_in(temp_arena);
        defer { md_temp_end(temp_scope); };

        ImGui::SetNextWindowSize(ImVec2(800, 500), ImGuiCond_FirstUseEver);
        if (ImGui::Begin("Density Volume", &show_window, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoFocusOnAppearing)) {
            if (ImGui::IsWindowFocused(ImGuiFocusedFlags_RootAndChildWindows) && ImGui::IsKeyPressed(KEY_PLAY_PAUSE, false) && !ImGui::GetIO().WantTextInput) {
                state->animation.mode = state->animation.mode == PlaybackMode::Playing ? PlaybackMode::Stopped : PlaybackMode::Playing;
            }

            if (ImGui::BeginMenuBar()) {
                if (ImGui::BeginMenu("Property")) {
                    draw_volume_choices(state);
                    if (volume_count > 1) {
                        ImGui::SeparatorText("Volume of the array");
                        draw_member_slider(state, 160.0f);
                    }
                    ImGui::EndMenu();
                }
                // The settings beside the view, at the right end of the bar
                const char* label = ICON_FA_SLIDERS;
                const float w = ImGui::CalcTextSize(label).x + 2.0f * ImGui::GetStyle().ItemSpacing.x;
                ImGui::SetCursorPosX(ImGui::GetCursorPosX() + MAX(ImGui::GetContentRegionAvail().x - w, 0.0f));
                if (ImGui::MenuItem(label, nullptr, show_panel)) {
                    show_panel = !show_panel;
                }
                ImGui::SetItemTooltip(show_panel ? "Hide the settings" : "Show the settings");
                ImGui::EndMenuBar();
            }

            update(state);

            // The view on the left, the settings on the right, a splitter between them which sets their width
            const ImVec2 origin = ImGui::GetCursorScreenPos();
            const ImVec2 avail = ImGui::GetContentRegionAvail();
            const float splitter_w = 6.0f;
            float panel_w = 0.0f;
            if (show_panel) {
                panel_w = CLAMP(panel_width, 200.0f, MAX(avail.x - 150.0f, 200.0f));
                const float canvas_w = MAX(avail.x - panel_w - splitter_w, 50.0f);

                ImGui::SetCursorScreenPos(ImVec2(origin.x + canvas_w, origin.y));
                ImGui::InvisibleButton("##splitter", ImVec2(splitter_w, MAX(avail.y, 1.0f)));
                if (ImGui::IsItemActive()) {
                    panel_width = CLAMP(panel_w - ImGui::GetIO().MouseDelta.x, 200.0f, MAX(avail.x - 150.0f, 200.0f));
                }
                if (ImGui::IsItemHovered() || ImGui::IsItemActive()) {
                    ImGui::SetMouseCursor(ImGuiMouseCursor_ResizeEW);
                }
                const float x = origin.x + canvas_w + 0.5f * splitter_w;
                ImGui::GetWindowDrawList()->AddLine(ImVec2(x, origin.y), ImVec2(x, origin.y + avail.y),
                                                    ImGui::GetColorU32(ImGui::IsItemHovered() || ImGui::IsItemActive() ? ImGuiCol_SeparatorHovered : ImGuiCol_Separator));

                ImGui::SetCursorScreenPos(ImVec2(origin.x + canvas_w + splitter_w, origin.y));
                if (ImGui::BeginChild("##panel", ImVec2(panel_w, MAX(avail.y, 1.0f)), ImGuiChildFlags_None)) {
                    draw_panel(state);
                }
                ImGui::EndChild();
                ImGui::SetCursorScreenPos(origin);
                panel_w += splitter_w;
            }

            // Animate camera towards targets
            camera_animate(&camera, target, state->app.timing.delta_s);

            const ImVec2 canvas_sz = ImMax(ImVec2(avail.x - panel_w, avail.y), ImVec2(50.0f, 50.0f));

            int width  = (int)(canvas_sz.x * ImGui::GetIO().DisplayFramebufferScale.x);
            int height = (int)(canvas_sz.y * ImGui::GetIO().DisplayFramebufferScale.y);
            if ((int)gbuf.width != width || (int)gbuf.height != height) {
                gbuffer_init(&gbuf, width, height);
            }

            const float aspect_ratio = canvas_sz.x / canvas_sz.y;
            const mat4_t world_to_view = camera_world_to_view_matrix(camera);
            const mat4_t view_to_world = camera_view_to_world_matrix(camera);

            const mat4_t view_to_clip  = camera_view_to_clip_matrix_persp(camera, aspect_ratio);
            const mat4_t clip_to_view  = camera_clip_to_view_matrix_persp(camera, aspect_ratio);

            const mat4_t world_to_clip = view_to_clip * world_to_view;
            const mat4_t clip_to_world = view_to_world * clip_to_view;

            const TrackballControllerParam trackball_param = {
                .min_distance = 1.0,
                .max_distance = 1000.0,
            };

            const InteractionSurfaceFlags surface_flags = InteractionSurfaceFlags_NoRegionSelect;
            InteractionSurfaceState surface_state = interaction_surface(interaction_surface_density_vol, vec_cast(canvas_sz), surface_flags);
            const ImVec2 canvas_p0 = ImGui::GetItemRectMin();
            const ImVec2 canvas_p1 = ImGui::GetItemRectMax();
            const ImRect canvas_rect = ImRect(canvas_p0, canvas_p1);
            PickingHit hit = {};

            // A volume dragged here (from its mark in the script editor) is the one shown
            if (ImGui::BeginDragDropTarget()) {
                if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload(VOLUME_SERIES_DND)) {
                    volume_key = ((const SeriesDragPayload*)payload->Data)->key;
                    volume_member = 0;
                }
                ImGui::EndDragDropTarget();
            }

            if (surface_state.hovered) {
                InteractionSurfaceHitArgs args = {
                    .picking_surface = &picking_surface,
                    .picking_handler = state->picking_handler,
                    .fbo = gbuf.fbo,
                    .width = gbuf.width,
                    .height = gbuf.height,
                    .clip_to_world = clip_to_world,
                };

                interaction_surface_hit_extract(&hit, surface_state, args);

                InteractionSurfaceEvent event = {};
                interaction_surface_event_extract(&event, surface_state, hit);

                event.clip_to_world = clip_to_world;
                event.world_to_clip = world_to_clip;

                if (event.kind != InteractionSurfaceEventKind::RegionSelect) {
                    viamd::event_system_broadcast_event(viamd::EventType_ViamdInteractionSurface, viamd::EventPayloadType_InteractionSurfaceEvent, &event);
                }
            }

            InteractionSurfaceViewTransformArgs view_args = {
                .camera = camera,
                .trackball_param = trackball_param,
            };

            InteractionSurfaceViewTransformResult view_result = interaction_surface_view_transform_apply(&target, surface_state, view_args);
            if (view_result.reset_requested) {
                ViewTransform reset_transform = default_view;
                if (hit.depth < 1.0f) {
                    reset_transform.distance = target.distance;
                    reset_transform.orientation = camera.orientation;
                    reset_transform.position = hit.world_pos + camera.orientation * vec3_set(0, 0, target.distance);
                }
                target = reset_transform;
            }

            // Draw border and background color
            ImDrawList* draw_list = ImGui::GetWindowDrawList();
            draw_list->AddImage((ImTextureID)(intptr_t)gbuf.tex.transparency, canvas_p0, canvas_p1, { 0,1 }, { 1,0 });
            draw_list->AddRect(canvas_p0, canvas_p1, IM_COL32(50, 50, 50, 255));

            if (volume_key.path[0] == '\0') {
                // Nothing shown: where to get something to show
                const char* text = "Pick a volume property of the script (Property),\nor drag one here from the script editor";
                const ImVec2 sz = ImGui::CalcTextSize(text);
                draw_list->AddText(canvas_rect.GetCenter() - sz * 0.5f, ImGui::GetColorU32(ImGuiCol_TextDisabled), text);
            }

            if (volume_state.status != SeriesStatus_Current || volume_count > 1) {
                // What is shown is not all of what the script says: how it stands, in the corner. Of an array of
                // volumes, which one it is ("v[2]") as well.
                char suffix[32], text[128], tip[256];
                series_state_suffix(suffix, sizeof(suffix), volume_state);
                snprintf(text, sizeof(text), "%s%s", volume_label, suffix);
                const ImVec2 pad = ImGui::GetStyle().FramePadding;
                const ImVec2 p0 = canvas_p0 + pad;
                const ImVec2 p1 = p0 + ImGui::CalcTextSize(text) + pad * 2.0f;
                const bool outdated = volume_state.status >= SeriesStatus_Outdated;
                draw_list->AddRectFilled(p0, p1, outdated ? IM_COL32(120, 80, 20, 220) : IM_COL32(40, 40, 40, 200), ImGui::GetStyle().FrameRounding);
                draw_list->AddText(p0 + pad, IM_COL32(255, 255, 255, 255), text);
                if (ImGui::IsMouseHoveringRect(p0, p1)) {
                    series_state_describe(tip, sizeof(tip), volume_state);
                    if (volume_count > 1) {
                        const size_t len = strnlen(tip, sizeof(tip));
                        snprintf(tip + len, sizeof(tip) - len, "%sVolume %u of the %u in the array: another is picked under Volume",
                                 len ? "\n" : "", volume_member + 1, volume_count);
                    }
                    ImGui::SetTooltip("%s", tip);
                }
            }

            if (dvr.enabled && legend.enabled) {
                ImVec2 canvas_ext = canvas_p1 - canvas_p0;
                ImVec2 cmap_ext = {MIN(canvas_ext.x * 0.5f, 250.0f), MIN(canvas_ext.y * 0.25f, 30.0f)};
                ImVec2 cmap_pad = {10, 10};
                ImVec2 cmap_pos = canvas_p1 - ImVec2(cmap_ext.x, cmap_ext.y) - cmap_pad;
                ImPlotColormap cmap = dvr.tf.colormap;
                ImPlotContext& gp = *ImPlot::GetCurrentContext();
                ImU32 checker_bg = IM_COL32(255, 255, 255, 255);
                ImU32 checker_fg = IM_COL32(128, 128, 128, 255);
                float checker_size = 8.0f;
                ImVec2 checker_offset = ImVec2(0,0);

                int mode = legend.colormap_mode;

                ImVec2 opaque_scl = ImVec2(1,1);
                ImVec2 transp_scl = ImVec2(0,0);

                if (mode == LegendColorMapMode_Split) {
                    opaque_scl = ImVec2(1.0f, 0.5f);
                    transp_scl = ImVec2(0.0f, 0.5f);
                }

                ImRect opaque_rect = ImRect(cmap_pos, cmap_pos + cmap_ext * opaque_scl);
                ImRect transp_rect = ImRect(cmap_pos + cmap_ext * transp_scl, cmap_pos + cmap_ext);

                // Opaque
                if (mode == LegendColorMapMode_Opaque || mode == LegendColorMapMode_Split) {
                    ImPlot::RenderColorBar(gp.ColormapData.GetKeys(cmap),gp.ColormapData.GetKeyCount(cmap),*draw_list,opaque_rect,false,false,!gp.ColormapData.IsQual(cmap));
                }

                if (mode == LegendColorMapMode_Transparent || mode == LegendColorMapMode_Split) {
                    if (legend.checkerboard) {
                        // Checkerboard
                        ImGui::DrawCheckerboard(draw_list, transp_rect.Min, transp_rect.Max, checker_bg, checker_fg, checker_size, checker_offset);
                    }
                    // Transparent
                    draw_list->AddImage((ImTextureID)(intptr_t)dvr.tf.id, transp_rect.Min, transp_rect.Max);
                }

                // Boarder
                draw_list->AddRect(cmap_pos, cmap_pos + cmap_ext, IM_COL32(0, 0, 0, 255), 0.0f, 0, 0.1f);
            }

            if (show_coordinate_system_widget) {
                float  ext = MIN(canvas_rect.GetWidth(), canvas_rect.GetHeight()) * 0.2f;
                float  pad = 0.1f * ext;
			    ImVec2 size = { ext, ext };

			    ImGui::SetCursorScreenPos(ImVec2(canvas_p0.x + pad, canvas_p1.y - ext - pad));

                quat_t out_orientation = target.orientation;
                if (ImGui::CoordinateSystemWidget(&out_orientation, camera.orientation, size)) {
                    const vec3_t look_at = camera_get_look_at(target);
                    target.orientation = quat_normalize(out_orientation);
                    target.position = camera_position_from_look_at(look_at, target.orientation, target.distance);
                }
            }

            PUSH_GPU_SECTION("RENDER DENSITY VOLUME");
            gbuffer_clear(&gbuf);

            const GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT_COLOR, GL_COLOR_ATTACHMENT_NORMAL, GL_COLOR_ATTACHMENT_VELOCITY,
                GL_COLOR_ATTACHMENT_PICKING, GL_COLOR_ATTACHMENT_TRANSPARENCY };

            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);

            glEnable(GL_DEPTH_TEST);
            glDepthMask(GL_TRUE);
            glDepthFunc(GL_LESS);

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
            glDrawBuffers((int)ARRAY_SIZE(draw_buffers), draw_buffers);
            glViewport(0, 0, gbuf.width, gbuf.height);
            glScissor(0, 0,  gbuf.width, gbuf.height);

            const bool selected_property = volume_key.path[0] != '\0';

            if (selected_property && show_reference_structures) {
                // Each representation, once per structure shown: what its filter picks within that structure, where the
                // structure aligns with the reference (a rigid transform)
                const bool has_backbone = state->mold.sys.protein_backbone.range.count > 0;
                md_array(md_gl_draw_op_t) draw_ops = nullptr;
                md_array(md_gl_draw_op_t) glyph_ops = nullptr;
                CartoonGlyphs& glyphs = state->mold.cartoon_glyphs;
                for (size_t i = 0; i < md_array_size(reps); ++i) {
                    const DensityVolumeRep& r = reps[i];
                    if (!r.rep.enabled || !r.rep.type_is_valid || r.dirty || !r.rep.md_rep.id) continue;
                    if (r.rep.glyphs_shown && md_array_size(r.glyph_subsets) > 0) {
                        // Derived from the positions of the atoms when a cartoon of the main view has not
                        if (!glyphs.positions_valid) cartoon_glyphs_update_positions(&glyphs, state->mold.state);
                        md_array(md_gl_draw_op_t) ops = nullptr;
                        representation_glyph_draw_ops(&ops, *state, r.rep, temp_arena);
                        const size_t n = MIN(md_array_size(r.glyph_subsets), md_array_size(structure_mats));
                        for (size_t k = 0; k < md_array_size(ops); ++k) {
                            for (size_t j = 0; j < n; ++j) {
                                if (!r.glyph_subsets[j].id) continue;
                                md_gl_draw_op_t op = ops[k];
                                op.rep = (j > 0 && r.members_glyphs.id) ? r.members_glyphs : r.rep.glyph_rep;
                                op.subset = r.glyph_subsets[j];
                                op.copies.transforms = &structure_mats[j].elem[0][0];
                                op.copies.count = 1;
                                md_array_push(glyph_ops, op, temp_arena);
                            }
                        }
                    }
                    if ((r.rep.type == RepresentationType::Ribbons || r.rep.type == RepresentationType::Cartoon) && !has_backbone) continue;
                    const md_gl_draw_op_t base = representation_draw_op(r.rep);
                    const size_t n = MIN(md_array_size(r.subsets), md_array_size(structure_mats));
                    for (size_t j = 0; j < n; ++j) {
                        if (!r.subsets[j].id) continue;
                        md_gl_draw_op_t op = base;
                        op.rep = (j > 0 && r.members.id) ? r.members : r.rep.md_rep;
                        op.subset = r.subsets[j];
                        op.copies.transforms = &structure_mats[j].elem[0][0];
                        op.copies.count = 1;
                        md_array_push(draw_ops, op, temp_arena);
                    }
                }

                md_gl_draw_args_t draw_args = {
                    .shaders = state->gl.shaders,
                    .draw_operations = {
                        .count = md_array_size(draw_ops),
                        .ops = draw_ops
                    },
                    .view_transform = {
                        .view_matrix = &world_to_view.elem[0][0],
                        .proj_matrix = &view_to_clip.elem[0][0],
                    },
                    .picking_offset = gl_picking_offset(*state),
                };
                if (md_array_size(draw_ops) > 0) {
                    md_gl_draw(&draw_args);
                }
                // The glyphs, picked as those of the main view are
                if (md_array_size(glyph_ops) > 0) {
                    draw_args.draw_operations = { .count = md_array_size(glyph_ops), .ops = glyph_ops };
                    draw_args.picking_offset = gl_picking_offset_glyphs(*state);
                    md_gl_draw(&draw_args);
                }

                glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
            }

            bool iso_written = false;
            if (show_density_volume) {
                if (dvr.enabled) {
                    volume::DvrRenderDesc vol_desc = {
                        .render_target = {
                            .depth  = gbuf.tex.depth,
                            .color  = gbuf.tex.transparency,
                            .width  = gbuf.width,
                            .height = gbuf.height,
                        },
                        .texture = {
                            .density_volume = volume_texture.id,
                            .transfer_function = dvr.tf.id,
                        },
                        .matrix = {
                            .model = model_mat,
                            .view = world_to_view,
                            .proj = view_to_clip,
                            .inv_proj = clip_to_view,
                        },
                        .clip_volume = {
                            .min = clip_volume.min,
                            .max = clip_volume.max,
                        },
                        .tf = {
                            .min_value = dvr.tf.min_val,
                            .max_value = dvr.tf.max_val,
                        },
                    };
                    volume::render_dvr(vol_desc);
                } else if (iso.enabled) {
                    volume::IsoRenderDesc vol_desc = {
                        .render_target = {
                            .depth  = gbuf.tex.depth,
                            .color  = gbuf.tex.transparency_hdr,
                            .width  = gbuf.width,
                            .height = gbuf.height,
                            .clear_color = true,
                        },
                        .texture = {
                            .density_volume = volume_texture.id,
                        },
                        .matrix = {
                            .model = model_mat,
                            .view = world_to_view,
                            .proj = view_to_clip,
                            .inv_proj = clip_to_view,
                        },
                        .clip_volume = {
                            .min = clip_volume.min,
                            .max = clip_volume.max,
                        },
                        .iso = {
                            .count = iso.count,
                            .values = iso.values,
                            .colors = iso.colors,
                            .exact = state->settings.exact_isosurfaces,
                        },
                        .shading = {
                            .roughness = 0.3f,
                            .ior = 1.5f,
                            .mode = state->visuals.shading.mode,
                            .ambient = state->visuals.shading.ambient,
                            .specular = state->visuals.shading.specular,
                            .colors = state->visuals.shading.colors,
                            .light = state->visuals.shading.light,
                            .environment = state->visuals.shading.environment,
                        },
                    };
                    iso_written = volume::render_isosurfaces(vol_desc);
                }
            }

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, gbuf.fbo);
            glDrawBuffer(GL_COLOR_ATTACHMENT_TRANSPARENCY);
            glViewport(0, 0, gbuf.width, gbuf.height);
            glScissor(0, 0, gbuf.width, gbuf.height);

            if (show_bounding_box && selected_property) {
                glEnable(GL_DEPTH_TEST);
                glDepthMask(GL_TRUE);
                glEnable(GL_BLEND);
                // Into the transparency buffer, which holds premultiplied colour
                glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

				immediate::Queue* queue = immediate::queue_create("Density Volume Bounds");

                immediate::set_model(queue, model_mat);

                uint32_t box_color = convert_color(bounding_box_color);
                uint32_t clip_color = convert_color(clip_volume_color);
                immediate::box_wireframe(queue, {0,0,0}, {1,1,1}, box_color);
                immediate::box_wireframe(queue, clip_volume.min, clip_volume.max, clip_color);

                immediate::RenderParams params = {};
                params.view = world_to_view;
                params.proj = view_to_clip;
                immediate::render(queue, params);

				immediate::queue_destroy(queue);

                glDisable(GL_BLEND);
            }

            PUSH_GPU_SECTION("Postprocessing")
            postprocess_pipeline::Settings postprocess_settings = {};
            postprocess_pipeline::Inputs postprocess_inputs = {};

            postprocess_settings.background_color = state->visuals.background.color;
            postprocess_settings.shading.mode = state->visuals.shading.mode;
            postprocess_settings.shading.ambient = state->visuals.shading.ambient;
            postprocess_settings.shading.specular = state->visuals.shading.specular;
            postprocess_settings.shading.colors = state->visuals.shading.colors;
            postprocess_settings.shading.light = state->visuals.shading.light;
            postprocess_settings.shading.environment = state->visuals.shading.environment;
            postprocess_settings.tonemap.enabled = state->visuals.tonemapping.enabled;
            postprocess_settings.tonemap.mode = state->visuals.tonemapping.tonemapper;
            postprocess_settings.tonemap.exposure = state->visuals.tonemapping.exposure;
            postprocess_settings.tonemap.gamma = state->visuals.tonemapping.gamma;
            postprocess_settings.shadows.enabled = false;
            postprocess_settings.ssao.enabled = false;
            postprocess_settings.dof.enabled = false;
            postprocess_settings.fxaa.enabled = true;
            postprocess_settings.taa.enabled = false;
            postprocess_settings.sharpen.enabled = false;

            postprocess_inputs.depth = gbuf.tex.depth;
            postprocess_inputs.color = gbuf.tex.color;
            postprocess_inputs.normal = gbuf.tex.normal;
            postprocess_inputs.velocity = gbuf.tex.velocity;
            postprocess_inputs.transparency = gbuf.tex.transparency;
            postprocess_inputs.transparency_hdr = iso_written ? gbuf.tex.transparency_hdr : 0;

            ViewParam view_param = {
                .matrix = {
                    .curr = {
                        .view = world_to_view,
                        .proj = view_to_clip,
                        .norm = world_to_view,
                    },
                    .inv = {
                        .proj = clip_to_view,
                    }
                },
                .clip_planes = {
                    .near = camera.near_plane,
                    .far = camera.far_plane,
                },
                .resolution = {canvas_sz.x, canvas_sz.y},
                .fov_y = camera.fov_y,
            };

            postprocess_pipeline::execute(postprocess_inputs, postprocess_settings, view_param);
            POP_GPU_SECTION()

            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
            glDrawBuffer(GL_BACK);

            POP_GPU_SECTION();
        }

        ImGui::End();
    }

    void process_events(const viamd::Event* events, size_t num_events) final {
        for (size_t i = 0; i < num_events; ++i) {
            const viamd::Event e = events[i];
            switch (e.type) {
            case viamd::EventType_ViamdInitialize: {
                // Initialize component
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState* state = (ApplicationState*)e.payload;
                app = state;
                arena = md_arena_allocator_create(state->allocator.persistent, MEGABYTES(1));
                picking_surface_init(&picking_surface, interaction_surface_density_vol);
                workspace_register_window("DensityVolume", &show_window);
                add_rep();
                break;
            }
            case viamd::EventType_ViamdSerialize:
                serialize(*(viamd::serialization_state_t*)e.payload);
                break;
            case viamd::EventType_ViamdDeserializeBegin:
                reset_workspace_settings();
                break;
            case viamd::EventType_ViamdDeserialize: {
                viamd::deserialization_state_t& state = *(viamd::deserialization_state_t*)e.payload;
                const str_t section = viamd::section_header(state);
                if (str_eq(section, STR_LIT("DensityVolume"))) {
                    deserialize(state);
                } else if (str_eq(section, STR_LIT("DensityVolumeRepresentation"))) {
                    deserialize_rep(state);
                }
                break;
            }
            case viamd::EventType_ViamdDeserializeEnd:
                // Over the timeline filter, from before there were periods: over the period it followed
                if (app) volume_key.period = series_legacy_period(app, volume_key.period);
                deserialize_end();
                break;
            case viamd::EventType_ViamdShutdown: {
                // Cleanup component. The GPU side goes with md_gl.
                clear_structures();
                md_arena_allocator_destroy(arena);
                break;
            }
            case viamd::EventType_ViamdSystemInit: {
                dirty_rep = true;
                break;
            }
            case viamd::EventType_ViamdSystemFree: {
                // What the representations hold is of the molecule that goes; their settings stay
                release_all_rep_resources();
                clear_structures();
                model_mat = {0};
                break;
            }
            case viamd::EventType_ViamdFrameTick: {
                ASSERT(e.payload_type == viamd::EventPayloadType_ApplicationState);
                ApplicationState* state = (ApplicationState*)e.payload;
                draw(state);
                break;
            }
            case viamd::EventType_ViamdWindowDrawMenu:
                ImGui::Checkbox("Density Volume", &show_window);
                break;
            default:
                break;
            }
        }
    }
};

static DensityVolume instance;
