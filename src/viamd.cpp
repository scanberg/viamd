
#include <md_util.h>
#include <md_filter.h>

#include <core/md_log.h>
#include <core/md_str_builder.h>
#include <core/md_arena_allocator.h>

#include <viamd.h>
#include <viamd_event.h>
#include <event.h>
#include <loader.h>
#include <color_utils.h>
#include <serialization_utils.h>

#include <gfx/gl_utils.h>
#include <gfx/volumerender_utils.h>

#include <imgui.h>
#include <imgui_notify.h>
#include <implot.h>
#include <implot_internal.h>

static const str_t* find_in_arr(str_t str, const str_t arr[], size_t len) {
    for (size_t i = 0; i < len; ++i) {
        if (str_eq(arr[i], str)) {
            return &arr[i];
        }
    }
    return NULL;
}

void init_volume(Volume* vol, const md_grid_t& grid, GLenum format) {
    ASSERT(vol);
    MEMCPY(vol->dim, grid.dim, sizeof(vol->dim));

    // A grid is in BOHR and a Volume's transforms are in Angstrom; this is where the two meet.
    const float scl = (float)BOHR_TO_ANGSTROM;

    vec3_t extent = md_grid_extent(&grid);
    vol->world_to_model   = compute_world_to_model_mat(grid.orientation, grid.origin * scl);
    vol->texture_to_world = compute_texture_to_world_mat(grid.orientation, grid.origin * scl, extent * scl);
    vol->voxel_size       = grid.spacing * scl;
    gl::init_texture_3D(&vol->tex_id, vol->dim[0], vol->dim[1], vol->dim[2], format);
}

static void init_all_representations(ApplicationState* state);

static void fill_picking_tooltip_text(md_strb_t* sb, const ApplicationState& state, const PickingHit& hit) {
    ASSERT(sb);
    const md_system_t& sys = state.mold.sys;
	const md_system_state_t& sys_state = state.mold.state;

    if (hit.domain == PickingDomain_Atom && hit.local_idx < sys.atom.count) {
        int atom_idx = hit.local_idx;
        int local_idx = atom_idx;
        const vec3_t pos = md_state_coord(&sys_state, atom_idx);
        str_t type = md_atom_name(&sys.atom, atom_idx);
        md_atomic_number_t z = md_atom_atomic_number(&sys.atom, atom_idx);
        str_t elem = z ? md_util_element_name(z)   : str_t{};
        str_t symb = z ? md_util_element_symbol(z) : str_t{};

        int comp_idx = md_component_find_by_atom_idx(&sys.component, atom_idx);
        str_t comp_name = {};
        int comp_seq_id = 0;
        if (comp_idx != -1) {
            comp_name   = md_component_name(&sys.component, comp_idx);
            comp_seq_id = md_component_seq_id(&sys.component, comp_idx);
            md_urange_t range = md_component_atom_range(&sys.component, comp_idx);
            local_idx = atom_idx - range.beg;
        }

        int inst_idx = md_system_instance_find_by_atom_idx(&sys, atom_idx);
        str_t inst_id = {};
		str_t auth_id = {};
        if (inst_idx != -1) {
            inst_id = md_instance_id(&sys.instance, inst_idx);
			auth_id = md_instance_auth_id(&sys.instance, inst_idx);
        }

        // @NOTE(Robin): External indices begin with 1 not 0
		if (state.selection.granularity == SelectionGranularity::Atom) {
            md_strb_fmt(sb, "atom[%i]", atom_idx + 1);
            if (comp_idx != -1) {
                md_strb_fmt(sb, "[%i]: ", local_idx + 1);
            } else {
				md_strb_push_cstr(sb, ": ");
            }
            md_strb_push_str(sb, type);
			md_strb_push_char(sb, ' ');
            if (z) {
                md_strb_fmt(sb, "%.*s %.*s ", STR_ARG(elem), STR_ARG(symb));
            }
            md_strb_fmt(sb, "(%.3f, %.3f, %.3f)\n", pos.x, pos.y, pos.z);    
        }
        
        if (comp_idx != -1 && (state.selection.granularity == SelectionGranularity::Atom || state.selection.granularity == SelectionGranularity::Component)) {
            md_strb_fmt(sb, "component[%i]", comp_idx + 1);
            if (comp_name) {
                md_strb_fmt(sb, ": " STR_FMT, STR_ARG(comp_name));
            }
			md_strb_fmt(sb, " (seq_id: %i)\n", comp_seq_id);
        }
        if (inst_idx != -1) {
            md_strb_fmt(sb, "structure[%i]", inst_idx + 1);
            if (inst_id) {
                md_strb_fmt(sb, ": " STR_FMT, STR_ARG(inst_id));
            }
			if (auth_id) {
                md_strb_fmt(sb, " (" STR_FMT ")", STR_ARG(auth_id));
            }
            md_strb_push_char(sb, '\n');
        }

        uint32_t flags = 0;

        if (state.selection.granularity == SelectionGranularity::Atom) {
			flags = md_system_atom_flags(&sys, atom_idx);
		} else if (state.selection.granularity == SelectionGranularity::Component && comp_idx != -1) {
			flags = md_system_component_flags(&sys, comp_idx);
        } else if (state.selection.granularity == SelectionGranularity::Instance && inst_idx != -1) {
			flags = md_system_instance_flags(&sys, inst_idx);
        }

		const uint32_t TERM_N = MD_FLAG_AMINO_ACID   | MD_FLAG_TERMINAL_BEG;
		const uint32_t TERM_C = MD_FLAG_AMINO_ACID   | MD_FLAG_TERMINAL_END;
		const uint32_t TERM_5 = MD_FLAG_NUCLEIC_ACID | MD_FLAG_TERMINAL_BEG;
		const uint32_t TERM_3 = MD_FLAG_NUCLEIC_ACID | MD_FLAG_TERMINAL_END;

        if (flags) {
            *sb += "flags: ";
            if (flags & MD_FLAG_HETERO)         { *sb += "HETERO "; }
			if (flags & MD_FLAG_POLYPEPTIDE)    { *sb += "POLYPEPTIDE "; }
            if (flags & MD_FLAG_AMINO_ACID)     { *sb += "AMINO-ACID "; }
            if (flags & MD_FLAG_SIDE_CHAIN)     { *sb += "SIDE-CHAIN "; }
			if (flags & MD_FLAG_NUCLEIC_ACID)   { *sb += "NUCLEIC-ACID "; }
            if (flags & MD_FLAG_NUCLEOTIDE)     { *sb += "NUCLEOTIDE "; }
            if (flags & MD_FLAG_NUCLEOSIDE)     { *sb += "NUCLEOSIDE "; }
            if (flags & MD_FLAG_NUCLEOBASE)     { *sb += "NUCLEOBASE "; }
            if (flags & MD_FLAG_WATER)          { *sb += "WATER "; }
            if (flags & MD_FLAG_ION)            { *sb += "ION "; }
            if (flags & MD_FLAG_BACKBONE)       { *sb += "BACKBONE "; }
            if ((flags & TERM_N) == TERM_N)     { *sb += "N-TERMINUS "; }
            if ((flags & TERM_C) == TERM_C)     { *sb += "C-TERMINUS "; }
			if ((flags & TERM_5) == TERM_5)     { *sb += "5'-TERMINUS "; }
			if ((flags & TERM_3) == TERM_3)     { *sb += "3'-TERMINUS "; }
            if (flags & MD_FLAG_SP)             { *sb += "SP "; }
            if (flags & MD_FLAG_SP2)            { *sb += "SP2 "; }
            if (flags & MD_FLAG_SP3)            { *sb += "SP3 "; }
            if (flags & MD_FLAG_AROMATIC)       { *sb += "AROMATIC "; }
            *sb += "\n";
        }
        /*
        // @TODO: REIMPLEMENT THIS
        if (res_idx < sys.backbone.segment.angleangles.size() && res_idx < sys.backbone.segments.size() && valid_backbone_atoms(sys.backbone.segments[res_idx])) {
        const auto angles = RAD_TO_DEG((vec2)sys.backbone.angles[res_idx]);
        len += snprintf(buff + len, 256 - len, u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0\n", angles.x, angles.y);
        }
        */
    }
    else if (hit.domain == PickingDomain_Bond) {
        int bond_idx = hit.local_idx;
        if (0 <= bond_idx && bond_idx < (int)sys.bond.count) {
            md_atom_pair_t   pair = sys.bond.pairs[bond_idx];
            md_bond_flags_t flags = sys.bond.flags[bond_idx];
            char bond_flags_buf[256] = {};
            int  len = 0;

            typedef struct {
                md_bond_flags_t flag;
                const char* label;
            } bond_flag_label_t;

            bond_flag_label_t bond_flag_map[] = {
                {MD_BOND_FLAG_COVALENT,     "COVALENT"},
                {MD_BOND_FLAG_DOUBLE,       "DOUBLE"},
                {MD_BOND_FLAG_TRIPLE,       "TRIPLE"},
                {MD_BOND_FLAG_QUADRUPLE,    "QUADRUPLE"},
                {MD_BOND_FLAG_AROMATIC,     "AROMATIC"},
                {MD_BOND_FLAG_COORDINATE,   "COORD"},
				{MD_BOND_FLAG_METAL,        "METAL"},
				{MD_BOND_FLAG_USER_DEFINED, "USER"},
            };

            for (size_t i = 0; i < ARRAY_SIZE(bond_flag_map); ++i) {
                if (flags & bond_flag_map[i].flag) {
                    len += snprintf(bond_flags_buf + len, sizeof(bond_flags_buf) - len, "%s ", bond_flag_map[i].label);
                }
            }
            
            char bond_type = '-';
            if (flags & MD_BOND_FLAG_DOUBLE) bond_type = '=';
            if (flags & MD_BOND_FLAG_TRIPLE) bond_type = '#';
            if (flags & MD_BOND_FLAG_QUADRUPLE) bond_type = '$';
            if (flags & MD_BOND_FLAG_AROMATIC) bond_type = ':';

            vec3_t p0 = md_state_coord(&sys_state, pair.idx[0]);
            vec3_t p1 = md_state_coord(&sys_state, pair.idx[1]);
            float d = vec3_distance(p0, p1);

            str_t type0 = md_atom_name(&sys.atom, pair.idx[0]);
            str_t type1 = md_atom_name(&sys.atom, pair.idx[1]);

            md_strb_fmt(sb, "bond: " STR_FMT "%c" STR_FMT "\n", STR_ARG(type0), bond_type, STR_ARG(type1));
            md_strb_fmt(sb, "flags: %.*s\n", len, bond_flags_buf);
            md_strb_fmt(sb, "length: %.3f\n", d);
        }
    } else if (hit.domain == PickingDomain_Dipole) {
        DipoleMoment dipoles[64];
        size_t num_dipoles = MIN(dipole_moments_gather(dipoles, ARRAY_SIZE(dipoles), state.mold.sys), ARRAY_SIZE(dipoles));
        int dipole_idx = hit.local_idx;
        if (0 <= dipole_idx && dipole_idx < (int)num_dipoles) {
            const DipoleMoment& d = dipoles[dipole_idx];

            char label[64];
            int label_len = dipole_entry_label(label, sizeof(label), d);
            md_strb_fmt(sb, "%.*s\n", label_len, label);

            // Debye is the readable unit for a dipole. The conversion is the attribute's own
            // business now, and it refuses rather than rescaling if the producer published
            // something which is not a dipole moment at all. One element of the group, because a
            // transition dipole attribute holds one per excited state.
            const md_attribute_t* attr = md_attributes_get(&state.mold.sys.attributes, d.key);
            float debye[3];
            const md_attribute_slice_t slice = (attr && attr->format.rank > 0) ? md_attribute_slice_1(d.index) : md_attribute_slice_all();
            if (attr && md_attribute_extract_slice_f32(debye, ARRAY_SIZE(debye), attr, &slice, md_unit_debye()) == 3) {
                md_strb_fmt(sb, "(%.3f %.3f %.3f) Debye\n", debye[0], debye[1], debye[2]);
            } else {
                char unit_buf[32];
                size_t unit_len = md_unit_print(unit_buf, sizeof(unit_buf), d.unit);
                md_strb_fmt(sb, "(%.3f %.3f %.3f) %.*s\n", d.vec.x, d.vec.y, d.vec.z, (int)unit_len, unit_buf);
            }
        }
    }
}

void draw_picking_tooltip_window(const PickingHit& hit, const ApplicationState& state) {
    if (hit.raw_idx == INVALID_PICKING_IDX) return;
    
	md_temp_scope_t temp = md_temp_begin_in(state.allocator.frame);
    defer { md_temp_end(temp); };

    PickingTooltipTextRequest tooltip_request = {
        .app = state,
        .hit = hit,
        .sb = md_strb_create(state.allocator.frame),
    };

    viamd::event_system_broadcast_event(viamd::EventType_ViamdPickingTooltipTextRequest, viamd::EventPayloadType_PickingTooltipTextRequest, &tooltip_request);

    if (!md_strb_empty(tooltip_request.sb)) {
        const ImVec2 offset = { 10.f, 18.f };
        const ImVec2 new_pos = {ImGui::GetMousePos().x + offset.x, ImGui::GetMousePos().y + offset.y};
        ImGui::SetNextWindowPos(new_pos);
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
        ImGui::Begin("##Picking Tooltip Window", 0,
            ImGuiWindowFlags_Tooltip | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoDocking);
        ImGui::Text("%s", md_strb_to_cstr(tooltip_request.sb));
        ImGui::End();
        ImGui::PopStyleColor();
    }
}

void interrupt_async_tasks(ApplicationState* state) {
    task_system::pool_interrupt_running_tasks();

    if (state->script.full_eval) md_script_eval_interrupt(state->script.full_eval);
    if (state->script.filt_eval) md_script_eval_interrupt(state->script.filt_eval);

    task_system::pool_wait_for_completion();
}

static inline void clear_frame_cache(FrameCache* cache) {
#if FRAME_CACHE_SIZE == 4
    md_mm_storeu_epi32(cache->frame_idx, md_mm_set1_epi32(-1));
    md_lru_cache4_init(&cache->lru);
#elif FRAME_CACHE_SIZE == 8
    md_mm256_storeu_epi32(cache->frame_idx, md_mm256_set1_epi32(-1));
    md_lru_cache8_init(&cache->lru);
#endif
}

static inline void init_frame_cache(FrameCache* cache, size_t num_atoms, md_allocator_i* alloc) {
    clear_frame_cache(cache);
    size_t capacity = ALIGN_TO(num_atoms, 16);
    for (size_t i = 0; i < FRAME_CACHE_SIZE; ++i) {
        md_array_resize(cache->states[i].x, capacity, alloc);
        md_array_resize(cache->states[i].y, capacity, alloc);
        md_array_resize(cache->states[i].z, capacity, alloc);
    }
}

static inline void free_frame_cache(FrameCache* cache, md_allocator_i* alloc) {
    for (size_t i = 0; i < FRAME_CACHE_SIZE; ++i) {
        md_array_free(cache->states[i].x, alloc);
        md_array_free(cache->states[i].y, alloc);
        md_array_free(cache->states[i].z, alloc);
    }
    clear_frame_cache(cache);
}

static inline bool find_frame_in_cache(int* out_slot_idx, int64_t frame_idx, const FrameCache* cache) {
    ASSERT(out_slot_idx);
#if FRAME_CACHE_SIZE == 4
    md_128i frame_indices = md_mm_loadu_epi32(cache->frame_idx);
    md_128i cmp_mask = md_mm_cmpeq_epi32(frame_indices, md_mm_set1_epi32((int32_t)frame_idx));
    int mask = md_mm_movemask_epi8(cmp_mask);
    *out_slot_idx = ctz32(mask) >> 2; // Each int32 comparison results in 4 bytes in the mask
    return mask != 0;
#elif FRAME_CACHE_SIZE == 8
    md_256i frame_indices = md_mm256_loadu_epi32(cache->frame_idx);
    md_256i cmp_mask = md_mm256_cmpeq_epi32(frame_indices, md_mm256_set1_epi32((int32_t)frame_idx));
    int mask = md_mm256_movemask_epi8(cmp_mask);
    *out_slot_idx = ctz32(mask) >> 2; // Each int32 comparison results in 4 bytes in the mask
    return mask != 0;
#endif
}

static inline int find_lru_cache_slot(const FrameCache* cache) {
#if FRAME_CACHE_SIZE == 4
    return md_lru_cache4_get_lru(cache->lru);
#elif FRAME_CACHE_SIZE == 8
    return md_lru_cache8_get_lru(cache->lru);
#endif
}

static inline void set_mru_cache_slot(FrameCache* cache, int slot_idx) {
#if FRAME_CACHE_SIZE == 4
    md_lru_cache4_set_mru(&cache->lru, slot_idx);
#elif FRAME_CACHE_SIZE == 8
    md_lru_cache8_set_mru(&cache->lru, slot_idx);
#endif
}

void clear_system_frame_cache(ApplicationState* state) {
    ASSERT(state);
    clear_frame_cache(&state->mold.frame_cache);
}

enum class SecondaryStructureRenderClass : uint8_t {
    Coil,
    Helix,
    Sheet,
};

static SecondaryStructureRenderClass secondary_structure_render_class(md_secondary_structure_t ss) {
    switch (ss) {
    case MD_SECONDARY_STRUCTURE_HELIX_310:
    case MD_SECONDARY_STRUCTURE_HELIX_ALPHA:
    case MD_SECONDARY_STRUCTURE_HELIX_PI:
        return SecondaryStructureRenderClass::Helix;
    case MD_SECONDARY_STRUCTURE_BETA_SHEET:
    case MD_SECONDARY_STRUCTURE_BETA_BRIDGE:
        return SecondaryStructureRenderClass::Sheet;
    case MD_SECONDARY_STRUCTURE_COIL:
    case MD_SECONDARY_STRUCTURE_TURN:
    case MD_SECONDARY_STRUCTURE_BEND:
    case MD_SECONDARY_STRUCTURE_UNKNOWN:
    default:
        return SecondaryStructureRenderClass::Coil;
    }
}

static md_secondary_structure_t secondary_structure_from_render_class(SecondaryStructureRenderClass cls) {
    switch (cls) {
    case SecondaryStructureRenderClass::Helix:
        return MD_SECONDARY_STRUCTURE_HELIX_ALPHA;
    case SecondaryStructureRenderClass::Sheet:
        return MD_SECONDARY_STRUCTURE_BETA_SHEET;
    case SecondaryStructureRenderClass::Coil:
    default:
        return MD_SECONDARY_STRUCTURE_COIL;
    }
}

static void secondary_structure_render_denoise(md_secondary_structure_t* dst, const md_secondary_structure_t* src, size_t num_frames, size_t stride) {
    ASSERT(dst);
    ASSERT(src);

    if (num_frames == 0 || stride == 0) {
        return;
    }

    md_temp_scope_t temp_scope = md_temp_begin();
    defer { md_temp_end(temp_scope); };

    SecondaryStructureRenderClass* classes  = md_temp_alloc_array(temp_scope, SecondaryStructureRenderClass, num_frames);
    SecondaryStructureRenderClass* filtered = md_temp_alloc_array(temp_scope, SecondaryStructureRenderClass, num_frames);

    constexpr int window_radius = 2;

    for (size_t seg_idx = 0; seg_idx < stride; ++seg_idx) {
        for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
            classes[frame_idx] = secondary_structure_render_class(src[frame_idx * stride + seg_idx]);
        }

        for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
            int counts[3] = {0, 0, 0};
            size_t beg = frame_idx > window_radius ? frame_idx - window_radius : 0;
            size_t end = MIN(frame_idx + window_radius + 1, num_frames);
            for (size_t j = beg; j < end; ++j) {
                counts[(int)classes[j]] += 1;
            }

            SecondaryStructureRenderClass majority = classes[frame_idx];
            int majority_count = 0;
            for (int j = 0; j < 3; ++j) {
                if (counts[j] > majority_count) {
                    majority = (SecondaryStructureRenderClass)j;
                    majority_count = counts[j];
                }
            }

            filtered[frame_idx] = majority_count > (int)((end - beg) / 2) ? majority : classes[frame_idx];
        }

        size_t run_beg = 0;
        while (run_beg < num_frames) {
            size_t run_end = run_beg + 1;
            while (run_end < num_frames && filtered[run_end] == filtered[run_beg]) {
                ++run_end;
            }

            size_t run_len = run_end - run_beg;
            if (run_beg > 0 && run_end < num_frames && filtered[run_beg - 1] == filtered[run_end]) {
                SecondaryStructureRenderClass replacement = filtered[run_beg - 1];
                bool replace_single = run_len <= 1;
                bool replace_short_structured = run_len <= 2 && replacement != SecondaryStructureRenderClass::Coil;
                if (replace_single || replace_short_structured) {
                    for (size_t j = run_beg; j < run_end; ++j) {
                        filtered[j] = replacement;
                    }
                }
            }

            run_beg = run_end;
        }

        for (size_t frame_idx = 0; frame_idx < num_frames; ++frame_idx) {
            dst[frame_idx * stride + seg_idx] = secondary_structure_from_render_class(filtered[frame_idx]);
        }
    }
}

// #trajectorydata
void free_trajectory_data(ApplicationState* state) {
    ASSERT(state);

    if (state->mold.sys.trajectory) {
        md_trajectory_free(state->mold.sys.trajectory);
        state->mold.sys.trajectory = nullptr;
    }
    state->files.trajectory[0] = '\0';

    md_array_free(state->timeline.x_values,  state->allocator.persistent);
    for (size_t i = 0; i < md_array_size(state->display_properties); ++i) {
        md_array_free(state->display_properties[i].hist.bins, state->display_properties[i].hist.alloc);
        state->display_properties[i].hist.bins = nullptr;
    }
    md_array_free(state->display_properties, state->allocator.persistent);

    // backbone_angles.data and secondary_structure.data are VIEWS into sys.attributes and
    // are not md_arrays - freeing them here would hand the wrong pointer to the array allocator. The
    // table owns them and releases them with the system; dropping the views is all that is owed.
    state->trajectory_data.backbone_angles.data       = nullptr;
    state->trajectory_data.backbone_angles.stride     = 0;
    state->trajectory_data.backbone_angles.count      = 0;
    state->trajectory_data.secondary_structure.data   = nullptr;
    state->trajectory_data.secondary_structure.stride = 0;
    state->trajectory_data.secondary_structure.count  = 0;

    // The denoised render copy is still an ordinary array owned here.
    md_array_free(state->trajectory_data.secondary_structure_render.data,    state->allocator.persistent);

    free_frame_cache(&state->mold.frame_cache, state->allocator.persistent);
}

void init_trajectory_data(ApplicationState* data) {
    size_t num_frames = md_trajectory_num_frames(data->mold.sys.trajectory);
    if (num_frames > 0) {
        size_t min_frame = 0;
        size_t max_frame = num_frames - 1;
        md_trajectory_header_t header = {};
        md_trajectory_get_header(data->mold.sys.trajectory, &header);

        init_frame_cache(&data->mold.frame_cache, data->mold.sys.atom.count, data->allocator.persistent);

        ASSERT(header.frame_times);
        double min_time = header.frame_times[0];
        double max_time = header.frame_times[num_frames - 1];

        data->timeline.view_range = {min_time, max_time};
        data->timeline.filter.beg_frame = (double)min_frame;
        data->timeline.filter.end_frame = (double)max_frame;

        md_array_resize(data->timeline.x_values, num_frames, data->allocator.persistent);
        for (size_t i = 0; i < num_frames; ++i) {
            data->timeline.x_values[i] = (float)header.frame_times[i];
        }

        // The COORDINATE for the frame axis, published as an ordinary attribute rather than as
        // axis metadata on every temporal quantity. This is the ANCHORING idea applied to an index
        // axis: a dipole publishes its origin as a sibling, and a quantity indexed by frame gets
        // its frame times the same way - one attribute serving every temporal attribute in the
        // table, since they all share that axis.
        //
        // It also retires a convention. md_trajectory_header_t documents that an EMPTY time_unit
        // means the format could not tell us real time and the values are ordinals standing in for
        // it. Here that is structural: the unit is on the attribute, and a reader that cannot
        // express time simply publishes none.
        {
            data->mold.sys.attributes.num_frames = (uint32_t)num_frames;

            const md_attribute_desc_t time_desc = {
                .path   = STR_LIT("frame/time"),
                .format = {
                    .type = MD_ATTRIBUTE_TYPE_F64, .components = 1,
                    .rank = 1, .shape = { (uint32_t)num_frames },
                },
                // Temporal AND resident: the frame axis is its only axis, and a whole extract is a
                // copy of a few hundred kilobytes, which is exactly what a plot of the time axis
                // wants. That is why the whole-extract rule is about cost rather than the axis.
                .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
                .unit   = header.time_unit,
                .label  = STR_LIT("Time"),
                .data   = header.frame_times,
                .byte_size = num_frames * sizeof(double),
            };
            md_attributes_replace(&data->mold.sys.attributes, &time_desc);

            if (header.frame_steps) {
                const md_attribute_desc_t step_desc = {
                    .path   = STR_LIT("frame/step"),
                    .format = {
                        .type = MD_ATTRIBUTE_TYPE_I64, .components = 1,
                        .rank = 1, .shape = { (uint32_t)num_frames },
                    },
                    .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
                    .unit   = md_unit_none(),
                    .label  = STR_LIT("Step"),
                    .description = STR_LIT("The file's own notion of where a frame sits in the run, not the frame ordinal"),
                    .data   = header.frame_steps,
                    .byte_size = num_frames * sizeof(int64_t),
                };
                md_attributes_replace(&data->mold.sys.attributes, &step_desc);
            }
        }

        data->animation.frame = CLAMP(data->animation.frame, (double)min_frame, (double)max_frame);
        int64_t frame_idx = CLAMP((int64_t)(data->animation.frame + 0.5), 0, (int64_t)max_frame);

        md_trajectory_load_frame(data->mold.sys.trajectory, frame_idx, &data->mold.state);

        if (data->mold.sys.protein_backbone.segment.count > 0) {
            // The angles and the per frame secondary structure are TEMPORAL attributes: the table
            // owns the storage and the trajectory_data fields below are views onto it. That is the
            // point of the move - 'stride' and 'count' were a hand rolled shape {F,S}, and the one
            // place they could disagree with the buffer was here, in three lines repeated per
            // quantity. Now the shape IS the declaration, and md_attributes_query answers "what
            // varies over time in this dataset" without anybody maintaining a second list.
            //
            // Both are declared together so the frame axis and the segment axis are stated once for
            // the pair. Resident rather than computed on demand, because the ramachandran density
            // reads EVERY frame at once - a per plane provider would recompute the whole trajectory
            // each time the window redraws. Reserved-and-filled-lazily is the interesting third
            // option, and it only starts paying at trajectory lengths where this allocation hurts.
            const size_t num_segments = data->mold.sys.protein_backbone.segment.count;
            md_attributes_t* attributes = &data->mold.sys.attributes;

            // What the TEMPORAL tag is checked against. Set before anything temporal is published,
            // so a declaration whose outermost extent disagrees is refused at the point the mistake
            // is made rather than surviving until something reads past the end of it.
            attributes->num_frames = (uint32_t)num_frames;

            // md_secondary_structure_t is a 4 byte enum, so I32 is the storage and the PATH is what
            // tells a consumer these are labels rather than numbers to average - which is the rule
            // md_system.h already states for integral attributes.
            const md_attribute_desc_t ss_desc = {
                .path   = STR_LIT("backbone/secondary_structure"),
                .format = {
                    .type = MD_ATTRIBUTE_TYPE_I32, .components = 1,
                    .rank = 2, .shape = { (uint32_t)num_frames, (uint32_t)num_segments },
                },
                .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
                .unit   = md_unit_none(),
                .label  = STR_LIT("Secondary Structure"),
            };
            md_attribute_id_t ss_id = md_attributes_replace(attributes, &ss_desc);

            data->trajectory_data.secondary_structure.data   = (md_secondary_structure_t*)md_attributes_data(attributes, ss_id, MD_ATTRIBUTE_TYPE_I32);
            data->trajectory_data.secondary_structure.stride = num_segments;
            data->trajectory_data.secondary_structure.count  = num_segments * num_frames;

            const md_attribute_format_t angle_format = {
                .type = MD_ATTRIBUTE_TYPE_F32, .components = 2,   // phi, psi - one value, two parts
                .rank = 2, .shape = { (uint32_t)num_frames, (uint32_t)num_segments },
            };
            const md_attribute_desc_t angle_desc = {
                .path   = STR_LIT("backbone/angle"),
                .format = angle_format,
                .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
                .unit   = md_unit_radian(),
                .label  = STR_LIT("Backbone Angles"),
            };
            md_attribute_id_t angle_id = md_attributes_replace(attributes, &angle_desc);

            // md_backbone_angles_t is exactly {float phi; float psi;}, so the table's storage IS a
            // md_backbone_angles_t array - no repacking, and the existing consumers keep their type.
            data->trajectory_data.backbone_angles.data   = (md_backbone_angles_t*)md_attributes_data(attributes, angle_id, MD_ATTRIBUTE_TYPE_F32);
            data->trajectory_data.backbone_angles.stride = num_segments;
            data->trajectory_data.backbone_angles.count  = num_segments * num_frames;

            // The denoised copy stays an ordinary array: it is a PRESENTATION smoothing of the one
            // above, read only by the renderer, and publishing it would put two answers to "what is
            // the secondary structure at frame f" in the same table.
            data->trajectory_data.secondary_structure_render.stride = data->mold.sys.protein_backbone.segment.count;
            data->trajectory_data.secondary_structure_render.count = data->mold.sys.protein_backbone.segment.count * num_frames;
            md_array_resize(data->trajectory_data.secondary_structure_render.data, data->mold.sys.protein_backbone.segment.count * num_frames, data->allocator.persistent);
            MEMSET(data->trajectory_data.secondary_structure_render.data, 0, md_array_bytes(data->trajectory_data.secondary_structure_render.data));

            // Launch work to compute the values
            task_system::task_interrupt_and_wait_for(data->tasks.backbone_computations);

            data->tasks.backbone_computations = task_system::create_pool_task(STR_LIT("Backbone Operations"), (uint32_t)num_frames, [data](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                // Create copy here of molecule since we use the full structure as input
                md_system_t sys = data->mold.sys;
                md_trajectory_i* traj = data->mold.sys.trajectory;

                md_temp_scope_t temp = md_temp_begin();
                defer { md_temp_end(temp); };

                md_system_state_t frame_state = { .alloc = temp.arena };
				md_system_state_init(&frame_state, sys.atom.count);

                md_trajectory_reader_i reader;
                if (md_trajectory_reader_init(&reader, traj)) {
                    for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                        md_backbone_angles_t* bb_dst = data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx;
                        md_secondary_structure_t* ss_dst = data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx;
                        
                        md_trajectory_reader_load_frame(reader, frame_idx, &frame_state);
                        md_util_backbone_angles_compute(bb_dst, data->trajectory_data.backbone_angles.stride, frame_state.x, frame_state.y, frame_state.z, &frame_state.unitcell, &sys.protein_backbone);
                        md_util_backbone_secondary_structure_infer(ss_dst, data->trajectory_data.secondary_structure.stride, frame_state.x, frame_state.y, frame_state.z, &frame_state.unitcell, &sys.protein_backbone);
                    }
					md_trajectory_reader_free(&reader);
                } else {
                    for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                        md_backbone_angles_t* bb_dst = data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx;
                        md_secondary_structure_t* ss_dst = data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx;
                        
                        md_trajectory_load_frame(traj, frame_idx, &frame_state);
                        md_util_backbone_angles_compute(bb_dst, data->trajectory_data.backbone_angles.stride, frame_state.x, frame_state.y, frame_state.z, &frame_state.unitcell, &sys.protein_backbone);
                        md_util_backbone_secondary_structure_infer(ss_dst, data->trajectory_data.secondary_structure.stride, frame_state.x, frame_state.y, frame_state.z, &frame_state.unitcell, &sys.protein_backbone);
                    }
                }
            });

            uint64_t time = (uint64_t)md_time_now();
            task_system::ID main_task = task_system::create_main_task(STR_LIT("Update Trajectory Data"), [data, t0 = time, num_frames]() {
                secondary_structure_render_denoise(
                    data->trajectory_data.secondary_structure_render.data,
                    data->trajectory_data.secondary_structure.data,
                    num_frames,
                    data->trajectory_data.secondary_structure.stride);

                uint64_t t1 = (uint64_t)md_time_now();
                double elapsed = md_time_as_seconds(t1 - t0);
                MD_LOG_INFO("Finished computing trajectory data (%.2fs)", elapsed);
                // The fill is finished, so the two temporal attributes have changed - which is the
                // one moment a consumer caching something derived from them needs to hear about.
                // The table carries that now; the hand rolled fingerprints beside the data are gone,
                // because a stamp that lives next to storage the table owns is a second answer
                // waiting to disagree with it.
                md_attributes_t* attributes = &data->mold.sys.attributes;
                md_attributes_touch(attributes, md_attributes_id_from_path(STR_LIT("backbone/angle")));
                md_attributes_touch(attributes, md_attributes_id_from_path(STR_LIT("backbone/secondary_structure")));
                data->trajectory_data.secondary_structure_render.fingerprint = generate_fingerprint();

				data->mold.interpolate_system_state = true;
                data->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
                flag_all_representations_as_dirty(data);
            });

            task_system::set_task_dependency(main_task, data->tasks.backbone_computations);
            task_system::enqueue_task(data->tasks.backbone_computations);
        }

        data->mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
        data->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;

        // Prefetch frames
        //launch_prefetch_job(data);
    }
}

void init_system_data(ApplicationState* data) {
    if (data->mold.sys.atom.count) {
        md_bitfield_clear(&data->operations.recenter_query.mask);
        md_bitfield_clear(&data->operations.selection_mask);
        recenter_mark_selection_dirty(data);
        data->operations.recenter_query.valid = false;
        data->operations.recenter_query.dynamic = false;
        data->operations.recenter_query.evaluated_version = 0;
        data->operations.recenter_query.ir_fingerprint = 0;
        data->operations.initial_frame.target_version = 0;
        recenter_mark_query_dirty(data);

        data->mold.gl_mol = md_gl_mol_create(&data->mold.sys);
        if (data->mold.sys.protein_backbone.segment.count > 0) {
            data->mold.interpolated_properties.secondary_structure = md_array_create(md_gl_secondary_structure_t, data->mold.sys.protein_backbone.segment.count, data->mold.sys.alloc);
        }

        mat3_t A;
		md_unitcell_A_extract_float(A.elem, &data->mold.state.unitcell);
		vec3_t c = mat3_mul_vec3(A, vec3_set(0.5f, 0.5f, 0.5f));
        data->mold.unitcell_transform = mat4_translate(-c.x, -c.y, -c.z);

        vec3_t aabb_min = {};
        vec3_t aabb_max = {};
        md_util_aabb_compute(aabb_min.elem, aabb_max.elem, data->mold.state.x, data->mold.state.y, data->mold.state.z, nullptr, nullptr, data->mold.state.num_atoms);

        const vec3_t cell_ext = mat3_mul_vec3(A, vec3_set1(1.0f));
        const float max_cell_ext = vec3_reduce_max(cell_ext);
        const float max_aabb_ext = vec3_reduce_max(vec3_sub(aabb_max, aabb_min));

        // Calculate a default view transform to use later as a reset target
        ViewTransform default_view = {};
        reset_view(&default_view, data->mold.state);

        data->view.camera.near_plane = 1.0f;
        data->view.camera.far_plane = 100000.0f;
        data->view.trackball_param.max_distance = MAX(max_cell_ext, max_aabb_ext) * 10.0f;

        data->view.target = default_view;
        data->view.camera = default_view;
        

#if EXPERIMENTAL_GFX_API
        const md_system_t& mol = data->mold.sys;
        vec3_t& aabb_min = data->mold.sys_aabb_min;
        vec3_t& aabb_max = data->mold.sys_aabb_max;
        md_util_compute_aabb_soa(&aabb_min, &aabb_max, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.radius, mol.atom.count);

        data->mold.gfx_structure = md_gfx_structure_create(mol.atom.count, mol.covalent.count, mol.backbone.count, mol.backbone.range_count, mol.residue.count, mol.instance.count);
        md_gfx_structure_set_atom_position(data->mold.gfx_structure, 0, mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
        md_gfx_structure_set_atom_radius(data->mold.gfx_structure, 0, mol.atom.count, mol.atom.radius, 0);
        md_gfx_structure_set_aabb(data->mold.gfx_structure, &data->mold.sys_aabb_min, &data->mold.sys_aabb_max);
        if (mol.instance.count > 0) {
            md_gfx_structure_set_instance_atom_ranges(data->mold.gfx_structure, 0, mol.instance.count, (md_gfx_range_t*)mol.instance.atom_range, 0);
            md_gfx_structure_set_instance_transforms(data->mold.gfx_structure, 0, mol.instance.count, mol.instance.transform, 0);
        }
#endif
    }
    viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemInit, viamd::EventPayloadType_ApplicationState, data);

    update_representation_info(data);
    init_all_representations(data);
    data->script.compile_ir = true;
}

void free_system_data(ApplicationState* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    // The arena OUTLIVES the system it backs, and the zeroing below is why there used to be a
    // second handle to it: MEMSET over md_system_t clears sys.alloc along with everything else, so
    // the only pointer to the arena has to be carried across by hand. Rewound, not destroyed - the
    // next load reuses it.
    md_allocator_i* sys_alloc = data->mold.sys.alloc;
    md_arena_allocator_reset(sys_alloc);
    MEMSET(&data->mold.sys, 0, sizeof(data->mold.sys));
    MEMSET(&data->mold.state, 0, sizeof(data->mold.state));
    data->mold.sys.alloc   = sys_alloc;
    data->mold.state.alloc = sys_alloc;

    md_array_free(data->operations.initial_frame.rel_xyzw, data->allocator.persistent);
    data->operations.initial_frame.rel_xyzw = nullptr;

    md_gl_mol_destroy(data->mold.gl_mol);

    // The dataset's GPU data goes with the dataset, for the same reason gl_mol does.
    system_gpu_data_free(data);

    MEMSET(data->files.molecule, 0, sizeof(data->files.molecule));

    data->mold.interpolated_properties.secondary_structure = nullptr;

    MEMSET(data->mold.frame_cache.states, 0, sizeof(data->mold.frame_cache.states));
    clear_frame_cache(&data->mold.frame_cache);

    md_bitfield_clear(&data->selection.selection_mask);
    md_bitfield_clear(&data->selection.highlight_mask);
    md_script_ir_t* ir = data->script.ir;
    md_script_ir_t* eval_ir = data->script.eval_ir;
    data->script.ir = nullptr;
    data->script.eval_ir = nullptr;
    if (ir) {
        md_script_ir_free(ir);
    }
    if (eval_ir && eval_ir != ir) {
        md_script_ir_free(eval_ir);
    }
    if (data->script.full_eval) {
        md_script_eval_free(data->script.full_eval);
        data->script.full_eval = nullptr;
    }
    if (data->script.filt_eval) {
        md_script_eval_free(data->script.filt_eval);
        data->script.filt_eval = nullptr;
    }

    viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemFree, viamd::EventPayloadType_ApplicationState, data);
}

bool load_data_from_file(ApplicationState* state, str_t filepath, const loader::LoaderState& load_state) {
    ASSERT(state);

    bool success = false;
    str_t path_to_file = md_path_make_canonical(filepath, state->allocator.frame);
    if (path_to_file) {
        if (load_state.flags & LoaderFlag_Supplemental) {
            // Do nothing here to hinder system or trajectory paths
        } else if (load_state.flags & LoaderFlag_System) {
            interrupt_async_tasks(state);
            free_trajectory_data(state);
            free_system_data(state);

            // sys.alloc and state.alloc are already the dataset arena - set at startup and put
            // back by free_system_data above. The loaded coordinates are the initial current state;
            // inference copies them into sys.reference, so the two stay independent.
            ASSERT(state->mold.sys.alloc && state->mold.state.alloc == state->mold.sys.alloc);
            if (!loader::load(&state->mold.sys, &state->mold.state, path_to_file, load_state)) {
                VIAMD_LOG_ERROR("Failed to load molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));
                return false;
            }
            success = true;
            VIAMD_LOG_SUCCESS("Successfully loaded molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));

            str_copy_to_char_buf(state->files.molecule, sizeof(state->files.molecule), path_to_file);
            state->files.coarse_grained = load_state.flags & LoaderFlag_CoarseGrained;
            // @NOTE: If the dataset is coarse-grained, then postprocessing must be aware
            md_infer_flags_t flags = state->files.coarse_grained ? MD_UTIL_INFER_NONE : MD_UTIL_INFER_ALL;
            md_util_system_infer(&state->mold.sys, &state->mold.state, flags);
            init_system_data(state);

            init_trajectory_data(state);
        } else if (load_state.flags & LoaderFlag_Trajectory) {
            if (!state->mold.sys.atom.count) {
                VIAMD_LOG_ERROR("Before loading a trajectory, molecular data needs to be present");
                return false;
            }
            interrupt_async_tasks(state);
            free_trajectory_data(state);
            state->animation.frame = 0;

            success = loader::load(&state->mold.sys, &state->mold.state, path_to_file, load_state);
            if (success) {
                init_trajectory_data(state);
                str_copy_to_char_buf(state->files.trajectory, sizeof(state->files.trajectory), path_to_file);
                VIAMD_LOG_SUCCESS("Successfully opened trajectory from file '" STR_FMT "'", STR_ARG(path_to_file));
            } else {
                VIAMD_LOG_ERROR("Failed to open trajectory from file '" STR_FMT "'", STR_ARG(path_to_file));
            }
        }
        LoadDataPayload data = {
            .app_state = state,
            .loader_state = load_state,
            .path_to_file = path_to_file,
        };
        viamd::event_system_broadcast_event(viamd::EventType_ViamdLoadData, viamd::EventPayloadType_LoadData, &data);
        update_representation_info(state);
    }

    return success;
}

void load_workspace(ApplicationState* data, str_t filename) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    str_t txt = load_textfile(filename, temp_alloc);

    if (str_empty(txt)) {
        VIAMD_LOG_ERROR("Could not open workspace file: '" STR_FMT "'", STR_ARG(filename));
        return;
    }

    // Reset and clear things
    remove_all_selections(data);
    remove_all_representations(data);
    data->editor.SetText("");
    data->files.workspace[0]  = '\0';

    data->animation = {};

    md_array(md_atom_pair_t) user_bonds = 0;    

    str_t new_molecule_file   = {};
    str_t new_trajectory_file = {};
    bool  new_coarse_grained  = false;
    double new_frame = 0;

    str_t folder = {};
    extract_folder_path(&folder, filename);

    viamd::deserialization_state_t state = {
        .filename = filename,
        .text = txt,
    };
    str_t section;
    while (viamd::next_section_header(section, state)) {
        if (str_eq(section, STR_LIT("Files")) || str_eq(section, STR_LIT("File"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("MoleculeFile"))) {
                    str_t file;
                    viamd::extract_str(file, arg);
                    if (!str_empty(file)) {
                        md_strb_t path = md_strb_create(temp_alloc);
                        // Paths are stored relative to the workspace, except when no relative
                        // path exists between the two (a separate volume on windows)
                        if (!md_path_is_absolute(file)) {
                            path += folder;
                        }
                        path += file;
                        new_molecule_file = md_path_make_canonical(path, temp_alloc);
                    }
                } else if (str_eq(ident, STR_LIT("TrajectoryFile"))) {
                    str_t file;
                    viamd::extract_str(file, arg);
                    if (!str_empty(file)) {
                        md_strb_t path = md_strb_create(temp_alloc);
                        if (!md_path_is_absolute(file)) {
                            path += folder;
                        }
                        path += file;
                        new_trajectory_file = md_path_make_canonical(path, temp_alloc);
                    }
                } else if (str_eq(ident, STR_LIT("CoarseGrained"))) {
                    viamd::extract_bool(new_coarse_grained, arg);
                }
            }
        } else if (str_eq(section, STR_LIT("Animation"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Frame"))) {
                    viamd::extract_dbl(new_frame, arg);
                } else if (str_eq(ident, STR_LIT("Fps"))) {
                    viamd::extract_flt(data->animation.fps, arg);
                } else if (str_eq(ident, STR_LIT("Interpolation"))) {
                    int mode;
                    viamd::extract_int(mode, arg);
                    data->animation.interpolation = (InterpolationMode)mode;
                }
            }
        } else if (str_eq(section, STR_LIT("RenderSettings"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("SsaoEnabled"))) {
                    viamd::extract_bool(data->visuals.ssao.enabled, arg);
                } else if (str_eq(ident, STR_LIT("SsaoIntensity"))) {
                    viamd::extract_flt(data->visuals.ssao.intensity, arg);
                } else if (str_eq(ident, STR_LIT("SsaoRadius"))) {
                    viamd::extract_flt(data->visuals.ssao.radius, arg);
                } else if (str_eq(ident, STR_LIT("SsaoBias"))) {
                    viamd::extract_flt(data->visuals.ssao.bias, arg);
                } else if (str_eq(ident, STR_LIT("DofEnabled"))) {
                    viamd::extract_bool(data->visuals.dof.enabled, arg);
                } else if (str_eq(ident, STR_LIT("DofFocusScale"))) {
                    viamd::extract_flt(data->visuals.dof.focus_scale, arg);
                }
            }
        } else if (str_eq(section, STR_LIT("Camera"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Position"))) {
                    viamd::extract_vec3(data->view.camera.position, arg);
                } else if (str_eq(ident, STR_LIT("Orientation"))) {
                    viamd::extract_quat(data->view.camera.orientation, arg);
                } else if (str_eq(ident, STR_LIT("Distance"))) {
                    viamd::extract_flt(data->view.camera.distance, arg);
                } else if (str_eq(ident, STR_LIT("Rotation"))) {
                    // DEPRECATED
                    viamd::extract_quat(data->view.camera.orientation, arg);
                }
            }
        } else if (str_eq(section, STR_LIT("Representation"))) {
            Representation* rep = create_representation(data);
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Name"))) {
                    viamd::extract_to_char_buf(rep->name, sizeof(rep->name), arg);
                } else if (str_eq(ident, STR_LIT("Filter"))) {
                    viamd::extract_to_char_buf(rep->filt, sizeof(rep->filt), arg);
                } else if (str_eq(ident, STR_LIT("Enabled"))) {
                    viamd::extract_bool(rep->enabled, arg);
                } else if (str_eq(ident, STR_LIT("Type"))) {
                    int type;
                    viamd::extract_int(type, arg);
                    type = CLAMP(type, 0, (int)RepresentationType::Count);
                    rep->type = (RepresentationType)type;
                } else if (str_eq(ident, STR_LIT("ColorMapping"))) {
                    int mapping;
                    viamd::extract_int(mapping, arg);
                    mapping = CLAMP(mapping, 0, (int)ColorMapping::Count);
                    rep->color_mapping = (ColorMapping)mapping;
                } else if (str_eq(ident, STR_LIT("StaticColor")) || str_eq(ident, STR_LIT("BaseColor"))) {
                    viamd::extract_flt_vec(rep->base_color.elem, 4, arg);
                } else if (str_eq(ident, STR_LIT("Saturation"))) {
                    viamd::extract_flt(rep->saturation, arg);
                } else if (str_eq(ident, STR_LIT("Radius"))) {
                    // DEPRECATED
                    viamd::extract_flt(rep->scale.x, arg);
                } else if (str_eq(ident, STR_LIT("Tension"))) {
                    // DEPRECATED
                } else if (str_eq(ident, STR_LIT("Width"))) {
                    viamd::extract_flt(rep->scale.x, arg);
                } else if (str_eq(ident, STR_LIT("Thickness"))) {
                    viamd::extract_flt(rep->scale.y, arg);
                } else if (str_eq(ident, STR_LIT("Param"))) {
                    viamd::extract_flt_vec(rep->scale.elem, 4, arg);
                } else if (str_eq(ident, STR_LIT("DynamicEval"))) {
                    viamd::extract_bool(rep->dynamic_evaluation, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureMoIdx")) || str_eq(ident, STR_LIT("ElectronicStructureOrbitalIdx"))) {
                    viamd::extract_int(rep->electronic_structure.orbital_idx, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureNtoIdx")) || str_eq(ident, STR_LIT("ElectronicStructureExcitedStateIdx"))) {
                    viamd::extract_int(rep->electronic_structure.excited_state_idx, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureNtoLambdaIdx"))) {
                    viamd::extract_int(rep->electronic_structure.nto_lambda_idx, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureDensityPropertyPath"))) {
                    // The id is a function of the path, so this resolves without the dataset being
                    // loaded yet - and names the same property after a reload, which the index it
                    // replaces did not.
                    rep->electronic_structure.density_property_key = md_attributes_id_from_path(arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureDensityPropertyIdx"))) {
                    // DEPRECATED. A position in a list that only existed while a particular file was
                    // open; there is nothing to resolve it against. Workspaces written before the
                    // path was stored fall back to the first property available.
                } else if (str_eq(ident, STR_LIT("ElectronicStructureRes"))) {
                    int res;
                    viamd::extract_int(res, arg);
                    rep->electronic_structure.resolution = (VolumeResolution)res;
                } else if (str_eq(ident, STR_LIT("ElectronicStructureType"))) {
                    int type;
                    viamd::extract_int(type, arg);
                    electronic_structure_set_legacy_type(&rep->electronic_structure, type);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureSource"))) {
                    int source;
                    viamd::extract_int(source, arg);
                    rep->electronic_structure.source = (ElectronicStructureSource)source;
                    electronic_structure_set_source_defaults(&rep->electronic_structure);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureField"))) {
                    int field;
                    viamd::extract_int(field, arg);
                    electronic_structure_set_legacy_field(&rep->electronic_structure, field);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureUseMagnitude"))) {
                    viamd::extract_bool(rep->electronic_structure.use_magnitude, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureSpin"))) {
                    int spin;
                    viamd::extract_int(spin, arg);
                    rep->electronic_structure.spin = (ElectronicStructureSpin)spin;
                } else if (str_eq(ident, STR_LIT("ElectronicStructureNtoComponent"))) {
                    int component;
                    viamd::extract_int(component, arg);
                    rep->electronic_structure.nto_component = (ElectronicStructureNtoComponent)component;
                } else if (str_eq(ident, STR_LIT("ElectronicStructureTransitionDensityComponent"))) {
                    int component;
                    viamd::extract_int(component, arg);
                    rep->electronic_structure.transition_density_component = (ElectronicStructureTransitionDensityComponent)component;
                } else if (str_eq(ident, STR_LIT("ElectronicStructureIso"))) {
                    viamd::extract_dbl(rep->electronic_structure.iso_value, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureColPos"))) {
                    viamd::extract_vec4(rep->electronic_structure.col_psi_pos, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureColNeg"))) {
                    viamd::extract_vec4(rep->electronic_structure.col_psi_neg, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureColDen"))) {
                    viamd::extract_vec4(rep->electronic_structure.col_den, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureColAtt"))) {
                    viamd::extract_vec4(rep->electronic_structure.col_att, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureColDet"))) {
                    viamd::extract_vec4(rep->electronic_structure.col_det, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureDensityPropertyIsoCount"))) {
                    viamd::extract_int(rep->electronic_structure.density_property.num_isos, arg);
                } else {
                    for (int i = 0; i < (int)ARRAY_SIZE(rep->electronic_structure.density_property.values); ++i) {
                        char key[64];
                        snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoValue%d", i);
                        if (str_eq(ident, str_from_cstr(key))) {
                            viamd::extract_dbl(rep->electronic_structure.density_property.values[i], arg);
                            break;
                        }
                        snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoColor%d", i);
                        if (str_eq(ident, str_from_cstr(key))) {
                            viamd::extract_vec4(rep->electronic_structure.density_property.colors[i], arg);
                            break;
                        }
                    }
                }
            }
        }
        else if (str_eq(section, STR_LIT("UserBonds"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("atoms"))) {
                    int atom_indices[2] = {-1, -1};
                    viamd::extract_int_vec(atom_indices, 2, arg);
                    if (atom_indices[0] >= 0 && atom_indices[1] >= 0) {
                        md_atom_pair_t pair = { atom_indices[0], atom_indices[1] };
                        md_array_push(user_bonds, pair, temp_alloc);
                    }
                }
            }
        }
        else if (str_eq(section, STR_LIT("Script"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Text"))) {
                    str_t str;
                    viamd::extract_str(str, arg);
                    data->editor.SetText(std::string(str.ptr, str.len));
                }
            }
        }
        else if (str_eq(section, STR_LIT("Selection"))) {
            str_t ident, arg;
            str_t label = {};
            str_t mask_base64 = {};
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Label"))) {
                    viamd::extract_str(label, arg);
                } else if (str_eq(ident, STR_LIT("Mask"))) {
                    viamd::extract_str(mask_base64, arg);
                }
            }
            if (!str_empty(label) && !str_empty(mask_base64)) {
                Selection* sel = create_selection(data, label);
                md_bitfield_deserialize(&sel->atom_mask, mask_base64.ptr, mask_base64.len);
            }
        } else {
            viamd::event_system_broadcast_event(viamd::EventType_ViamdDeserialize, viamd::EventPayloadType_DeserializationState, &state);
        }
    }

    data->view.target.position    = data->view.camera.position;
    data->view.target.orientation = data->view.camera.orientation;
    data->view.target.distance    = data->view.camera.distance;
    
    str_copy_to_char_buf(data->files.workspace, sizeof(data->files.workspace), filename);
    
    data->files.coarse_grained  = new_coarse_grained;

    loader::LoaderState loader_state = {};
    loader::init(&loader_state, new_molecule_file);

    if (new_coarse_grained) {
        loader_state.flags |= LoaderFlag_CoarseGrained;
    }

    if (new_molecule_file && load_data_from_file(data, new_molecule_file, loader_state)) {
        str_copy_to_char_buf(data->files.molecule, sizeof(data->files.molecule), new_molecule_file);
    } else {
        data->files.molecule[0]   = '\0';
    }

    if (new_trajectory_file) {
        loader::init(&loader_state, new_trajectory_file);
        if (load_data_from_file(data, new_trajectory_file, loader_state)) {
            str_copy_to_char_buf(data->files.trajectory, sizeof(data->files.trajectory), new_trajectory_file);
        }
        data->animation.frame = new_frame;
    } else {
        data->files.trajectory[0] = '\0';
    }

    // Add user bonds
    size_t num_user_bonds = md_array_size(user_bonds);
    if (num_user_bonds > 0) {
        for (size_t i = 0; i < num_user_bonds; ++i) {
            const md_atom_pair_t& pair = user_bonds[i];
            if (pair.idx[0] >= 0 && pair.idx[1] >= 0 && (size_t)pair.idx[0] < data->mold.sys.atom.count && (size_t)pair.idx[1] < data->mold.sys.atom.count) {
                md_system_bond_insert(&data->mold.sys, pair.idx[0], pair.idx[1], MD_BOND_FLAG_USER_DEFINED);
            }
        }
        data->mold.dirty_gpu_buffers |= MolBit_DirtyBonds;
    }

    //apply_atom_elem_mappings(data);
}

void save_workspace(ApplicationState* app_state, str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        VIAMD_LOG_ERROR("Could not open workspace file for writing: '%.*s", (int)filename.len, filename.ptr);
        return;
    }
    defer { md_file_close(&file); };

    md_allocator_i* temp_alloc = app_state->allocator.frame;

    viamd::serialization_state_t state {
        .filename = filename,
        .sb = md_strb_create(temp_alloc),
    };

    constexpr str_t header_snippet = STR_LIT(
        R"(
        #01010110#01001001#01000001#01001101#01000100#01001101#01000001#01001001#01010110#
        #                                                                                #
        #            VIAMD — Visual Interactive Analysis of Molecular Dynamics           #
        #                                                                                #
        #                    github: https://github.com/scanberg/viamd                   #
        #                 manual: https://github.com/scanberg/viamd/wiki                 #
        #                    youtube playlist: https://bit.ly/4aRsPrh                    #
        #                                twitter: @VIAMD_                                #
        #                                                                                #
        #                If you use VIAMD in your research, please cite:                 #
        #   "VIAMD: a Software for Visual Interactive Analysis of Molecular Dynamics"    #
        #       Robin Skånberg, Ingrid Hotz, Anders Ynnerman, and Mathieu Linares        #
        #                 J. Chem. Inf. Model. 2023, 63, 23, 7382–7391                   #
        #                   https://doi.org/10.1021/acs.jcim.3c01033                     #
        #                                                                                #
        #01010110#01001001#01000001#01001101#01000100#01001101#01000001#01001001#01010110#
        )");

    // Write big ass header
    state.sb += header_snippet;
    state.sb += '\n';

    // Files are stored relative to the workspace so that a dataset can be moved as a whole.
    // If no relative path exists (a separate volume on windows) we fall back to the absolute
    // path, and a slot which holds no file is written as an empty string.
    auto workspace_relative_path = [filename, temp_alloc](const char* file) -> str_t {
        if (!file || file[0] == '\0') {
            return {};
        }
        str_t path = str_from_cstr(file);
        str_t rel  = md_path_make_relative(filename, path, temp_alloc);
        return str_empty(rel) ? md_path_make_canonical(path, temp_alloc) : rel;
    };

    str_t mol_file  = workspace_relative_path(app_state->files.molecule);
    str_t traj_file = workspace_relative_path(app_state->files.trajectory);

    viamd::write_section_header(state, STR_LIT("Files"));
    viamd::write_str(state, STR_LIT("MoleculeFile"), mol_file);
    viamd::write_str(state, STR_LIT("TrajectoryFile"), traj_file);
    viamd::write_int(state, STR_LIT("CoarseGrained"), app_state->files.coarse_grained);

    viamd::write_section_header(state, STR_LIT("Animation"));
    viamd::write_dbl(state, STR_LIT("Frame"), app_state->animation.frame);
    viamd::write_flt(state, STR_LIT("Fps"), app_state->animation.fps);
    viamd::write_int(state, STR_LIT("Interpolation"), (int)app_state->animation.interpolation);

    viamd::write_section_header(state, STR_LIT("RenderSettings"));
    viamd::write_bool(state, STR_LIT("SsaoEnabled"), app_state->visuals.ssao.enabled);
    viamd::write_flt(state, STR_LIT("SsaoIntensity"), app_state->visuals.ssao.intensity);
    viamd::write_flt(state, STR_LIT("SsaoRadius"), app_state->visuals.ssao.radius);
    viamd::write_bool(state, STR_LIT("DofEnabled"), app_state->visuals.dof.enabled);
    viamd::write_flt(state, STR_LIT("DofFocusScale"), app_state->visuals.dof.focus_scale);

    viamd::write_section_header(state, STR_LIT("Camera"));
    viamd::write_vec3(state, STR_LIT("Position"), app_state->view.camera.position);
    viamd::write_quat(state, STR_LIT("Orientation"), app_state->view.camera.orientation);
    viamd::write_flt(state,  STR_LIT("Distance"), app_state->view.camera.distance);
    viamd::write_int(state,  STR_LIT("Mode"), (int)app_state->view.mode);


    for (size_t i = 0; i < md_array_size(app_state->representation.reps); ++i) {
        const Representation& rep = app_state->representation.reps[i];
        viamd::write_section_header(state, STR_LIT("Representation"));
        viamd::write_str(state,  STR_LIT("Name"), str_from_cstr(rep.name));
        viamd::write_str(state,  STR_LIT("Filter"), str_from_cstr(rep.filt));
        viamd::write_bool(state, STR_LIT("Enabled"), rep.enabled);
        viamd::write_int(state,  STR_LIT("Type"), (int)rep.type);
        viamd::write_int(state,  STR_LIT("ColorMapping"), (int)rep.color_mapping);
        viamd::write_vec4(state, STR_LIT("BaseColor"), rep.base_color);
        viamd::write_flt(state,  STR_LIT("Saturation"), rep.saturation);

        viamd::write_vec4(state, STR_LIT("Param"), rep.scale);
        viamd::write_bool(state, STR_LIT("DynamicEval"), rep.dynamic_evaluation);

        if (rep.type == RepresentationType::ElectronicStructure) {
            viamd::write_int(state,  STR_LIT("ElectronicStructureMoIdx"),    rep.electronic_structure.orbital_idx);
            viamd::write_int(state,  STR_LIT("ElectronicStructureNtoIdx"),   rep.electronic_structure.excited_state_idx);
            viamd::write_int(state,  STR_LIT("ElectronicStructureNtoLambdaIdx"), rep.electronic_structure.nto_lambda_idx);
            if (const md_attribute_t* prop = md_attributes_get(&app_state->mold.sys.attributes, rep.electronic_structure.density_property_key)) {
                viamd::write_str(state, STR_LIT("ElectronicStructureDensityPropertyPath"), prop->path);
            }
            viamd::write_int(state,  STR_LIT("ElectronicStructureType"),     electronic_structure_legacy_type(rep.electronic_structure));
            viamd::write_int(state,  STR_LIT("ElectronicStructureSource"),   (int)rep.electronic_structure.source);
            viamd::write_int(state,  STR_LIT("ElectronicStructureField"),    (int)electronic_structure_legacy_field(rep.electronic_structure));
            viamd::write_bool(state, STR_LIT("ElectronicStructureUseMagnitude"), rep.electronic_structure.use_magnitude);
            viamd::write_int(state,  STR_LIT("ElectronicStructureSpin"),     (int)rep.electronic_structure.spin);
            viamd::write_int(state,  STR_LIT("ElectronicStructureNtoComponent"), (int)rep.electronic_structure.nto_component);
            viamd::write_int(state,  STR_LIT("ElectronicStructureTransitionDensityComponent"), (int)rep.electronic_structure.transition_density_component);
            viamd::write_int(state,  STR_LIT("ElectronicStructureRes"), (int)rep.electronic_structure.resolution);
            viamd::write_dbl(state,  STR_LIT("ElectronicStructureIso"),      rep.electronic_structure.iso_value);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColPos"),   rep.electronic_structure.col_psi_pos);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColNeg"),   rep.electronic_structure.col_psi_neg);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColDen"),   rep.electronic_structure.col_den);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColAtt"),   rep.electronic_structure.col_att);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColDet"),   rep.electronic_structure.col_det);
            viamd::write_int(state,  STR_LIT("ElectronicStructureDensityPropertyIsoCount"), rep.electronic_structure.density_property.num_isos);
            for (int j = 0; j < rep.electronic_structure.density_property.num_isos; ++j) {
                char key[64];
                snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoValue%d", j);
                viamd::write_dbl(state, str_from_cstr(key), rep.electronic_structure.density_property.values[j]);
                snprintf(key, sizeof(key), "ElectronicStructureDensityPropertyIsoColor%d", j);
                viamd::write_vec4(state, str_from_cstr(key), rep.electronic_structure.density_property.colors[j]);
            }
        }
    }

    // TODO: Move atom element remappings to dataset component
    /*
    for (size_t i = 0; i < md_array_size(app_state->dataset.atom_element_remappings); ++i) {
        const AtomElementMapping& mapping = app_state->dataset.atom_element_remappings[i];
        viamd::write_section_header(state, STR_LIT("AtomElementMapping"));
        viamd::write_str(state, STR_LIT("Label"), str_from_cstr(mapping.lbl));
        viamd::write_int(state, STR_LIT("Element"), mapping.elem);
    }
    */

    {
        std::string text = app_state->editor.GetText();
        viamd::write_section_header(state, STR_LIT("Script"));
        viamd::write_str(state, STR_LIT("Text"), str_t{text.c_str(), text.size()});
    }

    for (size_t i = 0; i < md_array_size(app_state->selection.stored_selections); ++i) {
        const Selection& sel = app_state->selection.stored_selections[i];
        size_t cap = md_bitfield_serialize_size_in_bytes(&sel.atom_mask);
        char* buf = (char*)md_vm_arena_push(temp_alloc, cap);
        size_t len = md_bitfield_serialize(buf, &sel.atom_mask);
        str_t encoded_mask = {buf, len};
        viamd::write_section_header(state, STR_LIT("Selection"));
        viamd::write_str(state, STR_LIT("Label"), str_from_cstr(sel.name));
        viamd::write_str(state, STR_LIT("Mask"),  encoded_mask);
    }

    // Save user defined bonds
    md_array(md_atom_pair_t) user_bonds = 0;
    for (size_t i = 0; i < app_state->mold.sys.bond.count; ++i) {
        md_bond_flags_t flags = app_state->mold.sys.bond.flags[i];
        if (flags & MD_BOND_FLAG_USER_DEFINED) {
            md_array_push(user_bonds, app_state->mold.sys.bond.pairs[i], temp_alloc);
        }
    }
    if (md_array_size(user_bonds) > 0) {
        viamd::write_section_header(state, STR_LIT("UserBonds"));
        for (size_t i = 0; i < md_array_size(user_bonds); ++i) {
            const md_atom_pair_t& pair = user_bonds[i];
            viamd::write_int_vec(state, STR_LIT("atoms"), pair.idx, 2);
        }
    }

    viamd::event_system_broadcast_event(viamd::EventType_ViamdSerialize, viamd::EventPayloadType_SerializationState, &state);

    str_t text = md_strb_to_str(state.sb);
    md_file_write(file, str_ptr(text), str_len(text));
}

// --- SELECTION ---
Selection* create_selection(ApplicationState* state, str_t name, md_bitfield_t* atom_mask) {
    ASSERT(state);
    Selection sel;
    str_copy_to_char_buf(sel.name, sizeof(sel.name), name);
    md_bitfield_init(&sel.atom_mask, state->allocator.persistent);
    if (atom_mask) {
        md_bitfield_copy(&sel.atom_mask, atom_mask);
    }
    md_array_push(state->selection.stored_selections, sel, state->allocator.persistent);
    return md_array_last(state->selection.stored_selections);
}

void remove_all_selections(ApplicationState* state) {
    md_array_shrink(state->selection.stored_selections, 0);
}

void remove_selection(ApplicationState* state, size_t idx) {
    ASSERT(state);
    if (md_array_size(state->selection.stored_selections) <= idx) {
        VIAMD_LOG_ERROR("Index [%zu] out of range when trying to remove selection", idx);
    }
    auto item = &state->selection.stored_selections[idx];
    md_bitfield_free(&item->atom_mask);

    state->selection.stored_selections[idx] = *md_array_last(state->selection.stored_selections);
    md_array_pop(state->selection.stored_selections);
}

// --- Representation ---

static void init_representation(ApplicationState* state, Representation* rep) {
#if EXPERIMENTAL_GFX_API
    rep->gfx_rep = md_gfx_rep_create(state->mold.sys.atom.count);
#endif
    rep->md_rep = md_gl_rep_create(state->mold.gl_mol);
    md_bitfield_init(&rep->atom_mask, state->allocator.persistent);

    // Default to the first per atom field the system offers, if it offers any. There is no list to
    // consult: the attribute table is the list.
    md_attribute_id_t first_property = MD_ATTRIBUTE_INVALID;
    if (atom_property_query(&first_property, 1, state->mold.sys) > 0) {
        atom_property_select(&rep->atomic_property, first_property, state->mold.sys);
    }

    flag_representation_as_dirty(rep);
}

Representation* create_representation(ApplicationState* state, RepresentationType type, ColorMapping color_mapping, str_t filter) {
    ASSERT(state);
    md_array_push(state->representation.reps, Representation(), state->allocator.persistent);
    Representation* rep = md_array_last(state->representation.reps);
    rep->type = type;
    rep->color_mapping = color_mapping;
    if (!str_empty(filter)) {
        str_copy_to_char_buf(rep->filt, sizeof(rep->filt), filter);
    }
    rep->electronic_structure.orbital_idx = (int)state->representation.info.alpha.homo_idx;
    init_representation(state, rep);
    return rep;
}

Representation* clone_representation(ApplicationState* state, const Representation& rep) {
    ASSERT(state);
    md_array_push(state->representation.reps, rep, state->allocator.persistent);
    Representation* clone = md_array_last(state->representation.reps);
    clone->md_rep = {0};
    clone->atom_mask = {0};
    init_representation(state, clone);
    return clone;
}

void remove_representation(ApplicationState* state, size_t idx) {
    ASSERT(state);
    ASSERT(idx < md_array_size(state->representation.reps));
    auto& rep = state->representation.reps[idx];
    md_bitfield_free(&rep.atom_mask);
    md_gl_rep_destroy(rep.md_rep);
    // A readback queued for this representation's volume would land in a texture that no longer
    // exists, so let the queue run out first.
    gpu_volume_jobs_drain(state);
    if (rep.electronic_structure.density_vol.tex_id) gl::free_texture(&rep.electronic_structure.density_vol.tex_id);
    if (rep.electronic_structure.color_vol.tex_id)   gl::free_texture(&rep.electronic_structure.color_vol.tex_id);
    if (rep.electronic_structure.dvr.tf_tex)         gl::free_texture(&rep.electronic_structure.dvr.tf_tex);
    md_array_swap_back_and_pop(state->representation.reps, idx);
    recompute_atom_visibility_mask(state);
}

void recompute_atom_visibility_mask(ApplicationState* state) {
    ASSERT(state);
    auto& mask = state->representation.visibility_mask;

    md_bitfield_clear(&mask);
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        auto& rep = state->representation.reps[i];
        if (!rep.enabled) continue;
        md_bitfield_or_inplace(&mask, &rep.atom_mask);
    }
    state->representation.visibility_mask_hash = md_bitfield_hash64(&mask, 0);
}

// "ground_state" reads as "Ground State" in a menu. The path segment is the identity; this is
// presentation only, which is why it lives here and not in mdlib. Writes into a caller buffer so
// gathering stays allocation free.
int dipole_label_pretty(char* buf, size_t cap, str_t group) {
    if (!buf || cap == 0) return 0;
    size_t n = MIN(group.len, cap - 1);
    bool boundary = true;
    for (size_t i = 0; i < n; ++i) {
        char c = group.ptr[i];
        if (c == '_' || c == '-') {
            buf[i] = ' ';
            boundary = true;
            continue;
        }
        buf[i] = (boundary && c >= 'a' && c <= 'z') ? (char)(c - 'a' + 'A') : c;
        boundary = false;
    }
    buf[n] = '\0';
    return (int)n;
}

int dipole_entry_label(char* buf, size_t cap, const DipoleMoment& dipole) {
    if (!buf || cap == 0) return 0;

    int len = dipole_label_pretty(buf, cap, dipole.label);
    if (dipole.count > 1 && (size_t)len + 1 < cap) {
        len += snprintf(buf + len, cap - (size_t)len, " %u", dipole.index + 1);
    }
    return len;
}

size_t dipole_moments_gather(DipoleMoment out[], size_t cap, const md_system_t& sys) {
    str_t groups[64];
    size_t num_groups = md_attributes_query_children(groups, ARRAY_SIZE(groups), &sys.attributes, STR_LIT("dipole"));
    num_groups = MIN(num_groups, ARRAY_SIZE(groups));

    size_t count = 0;
    for (size_t i = 0; i < num_groups; ++i) {
        char path[256];
        int len = snprintf(path, sizeof(path), "dipole/" STR_FMT "/vector", STR_ARG(groups[i]));
        const md_attribute_t* vec = md_attributes_find(&sys.attributes, str_from_cstrn(path, len));
        len = snprintf(path, sizeof(path), "dipole/" STR_FMT "/origin", STR_ARG(groups[i]));
        const md_attribute_t* org = md_attributes_find(&sys.attributes, str_from_cstrn(path, len));

        // A group carrying only one half is not a dipole anyone can draw.
        if (!vec || !org) continue;
        if (md_attribute_components(&vec->format) != 3 || md_attribute_components(&org->format) != 3) continue;

        // The vector decides how many dipoles the group holds. The origin is addressed by the same
        // index space but need not have the same shape: an origin with fewer index axes is constant
        // over the ones it lacks, which is how one shared anchor serves every excited state.
        const size_t num_elem = md_attribute_value_count(&vec->format);
        const size_t num_org  = md_attribute_value_count(&org->format);
        if (num_elem == 0 || (num_org != num_elem && num_org != 1)) continue;

        for (size_t e = 0; e < num_elem; ++e) {
            if (out && count < cap) {
                // One slice, clamped to each half's own rank. A shared origin is rank 0 and takes
                // no index; a per state one is rank 1 and takes the same index as the vector. That
                // clamp is the whole cost of not requiring the two to have equal shapes.
                const md_attribute_slice_t vec_slice = vec->format.rank > 0 ? md_attribute_slice_1((uint32_t)e) : md_attribute_slice_all();
                const md_attribute_slice_t org_slice = org->format.rank > 0 ? vec_slice : md_attribute_slice_all();

                // The vector comes back as stored, whatever the producer chose; the origin must be a
                // length in system space, and extraction refuses if it is not. A refusal leaves the
                // entry at zero rather than dropping it, so that counting and writing agree on how
                // many there are - a picking index into this array depends on that.
                float v[3] = {0, 0, 0};
                float o[3] = {0, 0, 0};
                md_attribute_extract_slice_f32(v, ARRAY_SIZE(v), vec, &vec_slice, md_unit_none());
                md_attribute_extract_slice_f32(o, ARRAY_SIZE(o), org, &org_slice, md_unit_angstrom());

                out[count] = {
                    .key    = vec->id,
                    .index  = (uint32_t)e,
                    .count  = (uint32_t)num_elem,
                    .label  = groups[i],
                    .vec    = vec3_set(v[0], v[1], v[2]),
                    .origin = vec3_set(o[0], o[1], o[2]),
                    .unit   = vec->unit,
                };
            }
            count += 1;
        }
    }

    return count;
}


// A per atom scalar field: values one component wide, the atom axis LAST, and at most one axis of
// variants ahead of it. mdlib deliberately does not know that an "atom/..." path is over this
// system's atoms - categories are not predeclared - so this is where that convention is checked.
static bool atom_property_qualifies(const md_attribute_t* attr, const md_system_t& sys) {
    const md_attribute_format_t& fmt = attr->format;
    if (fmt.rank < 1 || fmt.rank > 2) return false;
    if (md_attribute_components(&fmt) != 1) return false;
    return fmt.shape[fmt.rank - 1] == (uint32_t)sys.atom.count;
}

size_t atom_property_query(md_attribute_id_t out_ids[], size_t cap, const md_system_t& sys) {
    md_attribute_id_t ids[128];
    size_t num_ids = md_attributes_query(ids, ARRAY_SIZE(ids), &sys.attributes, STR_LIT("atom"));
    num_ids = MIN(num_ids, ARRAY_SIZE(ids));

    size_t count = 0;
    for (size_t i = 0; i < num_ids; ++i) {
        const md_attribute_t* attr = md_attributes_get(&sys.attributes, ids[i]);
        if (!attr || !atom_property_qualifies(attr, sys)) continue;

        if (out_ids && count < cap) {
            out_ids[count] = ids[i];
        }
        count += 1;
    }

    return count;
}

str_t atom_property_label(const md_attribute_t* attr) {
    if (!attr) return str_t{};
    // An empty label is a valid state, and the leaf is what the path spells for itself.
    return str_empty(attr->label) ? md_attribute_leaf(attr) : attr->label;
}

int atom_property_variant_count(const md_attribute_t* attr) {
    if (!attr) return 0;
    return attr->format.rank > 1 ? (int)attr->format.shape[0] : 1;
}

bool atom_property_value_range(float* out_min, float* out_max, const md_attribute_t* attr) {
    if (!attr) return false;

    const size_t num_values = md_attribute_element_count(&attr->format);
    if (num_values == 0) return false;

    md_temp_scope_t temp = md_temp_begin();
    float* values = (float*)md_temp_alloc(temp, sizeof(float) * num_values);

    // Deliberately the whole attribute and not one variant: a span recomputed per variant would
    // make the colours shift as the index slider moves, which reads as the data changing.
    bool result = values && md_attribute_extract_f32(values, num_values, attr, md_unit_none()) == num_values;
    if (result) {
        float value_min =  FLT_MAX;
        float value_max = -FLT_MAX;
        for (size_t i = 0; i < num_values; ++i) {
            value_min = MIN(value_min, values[i]);
            value_max = MAX(value_max, values[i]);
        }
        if (out_min) *out_min = value_min;
        if (out_max) *out_max = value_max;
    }
    md_temp_end(temp);

    return result;
}

void atom_property_select(AtomicPropertyRepresentation* prop, md_attribute_id_t key, const md_system_t& sys) {
    ASSERT(prop);

    prop->key = key;
    prop->variant_idx = 0;

    float value_min = 0.0f;
    float value_max = 1.0f;
    atom_property_value_range(&value_min, &value_max, md_attributes_get(&sys.attributes, key));

    prop->value_min = value_min;
    prop->value_max = value_max;
    prop->range_beg = value_min;
    prop->range_end = value_max;
}

// ---------------------------------------------------------------------------
// Per dataset GPU data
// ---------------------------------------------------------------------------
// Built from the system's own basis/ attributes, so it needs no loader and no component - which is
// the point of publishing the basis in the first place. Lives beside gl_mol because it is the same
// kind of thing: derived FROM the system so that something can draw or evaluate it, and never read
// back by mdlib.
//
// The device scratch is grown here rather than allocated per dataset. A wider basis loaded later
// grows it; a narrower one leaves it alone, because it is scratch and only the maximum matters.

void system_gpu_data_free(ApplicationState* state) {
    ASSERT(state);
#if MD_ENABLE_GPU
    if (state->mold.gpu_basis) {
        md_gto_gpu_basis_destroy(state->mold.gpu_basis);
        state->mold.gpu_basis = nullptr;
    }
    if (state->mold.gpu_atoms) {
        md_gpu_free(state->mold.gpu_atoms, state->gpu_stream);
        state->mold.gpu_atoms = nullptr;
    }
    state->mold.gpu_atoms_hash = 0;
#else
    (void)state;
#endif
}

bool system_gpu_data_update(ApplicationState* state, double cutoff) {
    ASSERT(state);
#if MD_ENABLE_GPU
    system_gpu_data_free(state);

    if (!state->gpu_device || !state->gpu_pool) {
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &state->mold.sys.attributes, md_temp_allocator(temp))) {
        // A system with no basis published is the normal case, not a failure.
        return false;
    }

    md_gto_gpu_basis_desc_t desc = { .basis = &basis, .cutoff = cutoff };
    state->mold.gpu_basis = md_gto_gpu_basis_create(state->gpu_pool, state->gpu_stream, &desc);
    if (!state->mold.gpu_basis) {
        MD_LOG_ERROR("Failed to upload the GTO basis to the device");
        return false;
    }

    const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
    const size_t num_atoms = md_gto_gpu_basis_num_atoms(state->mold.gpu_basis);

    state->mold.gpu_atoms = md_gpu_malloc(state->gpu_pool, md_gto_gpu_atom_buffer_size(num_atoms), state->gpu_stream);
    state->mold.gpu_atoms_hash = 0;

    // Density coefficients are the larger of the two packings, so one size covers both the density
    // and the MO evaluation paths.
    const size_t coeff_size = md_gto_gpu_coeff_size_density(num_cgtos);
    if (coeff_size > state->gpu_coeff_capacity) {
        if (state->gpu_coeff) {
            md_gpu_free(state->gpu_coeff, state->gpu_stream);
        }
        state->gpu_coeff = md_gpu_malloc(state->gpu_pool, coeff_size, state->gpu_stream);
        state->gpu_coeff_capacity = state->gpu_coeff ? coeff_size : 0;
    }

    return state->mold.gpu_atoms != nullptr && state->gpu_coeff != nullptr;
#else
    (void)state; (void)cutoff;
    return false;
#endif
}

// ---------------------------------------------------------------------------
// Orbital evaluation
// ---------------------------------------------------------------------------
// Nothing here holds a reader. The basis and the AO coefficients are attributes on the system, the
// atom positions are the system's own state, and the grid, the destination texture and the
// evaluation parameters belong to the caller. There is no vlx pointer in scope, no event, and no
// component to ask - which is the whole point of publishing the coefficients.
//
// The basis is rebuilt per call for now. That is one interleave over a few hundred shells, but it
// is the obvious thing for a representation to cache, keyed on the ids of the basis/ attributes it
// was built from.
double* orbital_coefficients_extract(size_t* out_num_ao, md_temp_scope_t temp, const md_system_t& sys, str_t coefficient_path, const md_attribute_slice_t* slice) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, coefficient_path);
    if (!attr) {
        MD_LOG_DEBUG("No orbital coefficients published at '" STR_FMT "'", STR_ARG(coefficient_path));
        return nullptr;
    }

    // Whatever the attribute's own rank, the slice has to leave exactly ONE row of AO coefficients:
    // {M,A} sliced by the orbital, {S,L,A} sliced by state and lambda. Asking the slice for its
    // format rather than the attribute is what makes those the same call.
    md_attribute_format_t format = {};
    if (!md_attribute_slice_format(&format, attr, slice)) {
        MD_LOG_ERROR("The slice does not address '" STR_FMT "'", STR_ARG(coefficient_path));
        return nullptr;
    }
    if (format.rank != 1 || md_attribute_components(&format) != 1) {
        MD_LOG_ERROR("'" STR_FMT "' does not slice down to one row of orbital coefficients", STR_ARG(coefficient_path));
        return nullptr;
    }

    // Ask the slice how big it is, then allocate for exactly that. The size comes from the format
    // alone, so this same shape works whether the coefficients are stored or worked out on demand.
    const size_t num_ao = md_attribute_slice_count(attr, slice);
    if (num_ao == 0) {
        return nullptr;
    }

    double* dst = (double*)md_temp_alloc(temp, sizeof(double) * num_ao);
    if (!dst) {
        return nullptr;
    }

    // f64, because that is what md_gto takes: the coefficients are double at this boundary to keep
    // the QM code's precision, and extracting them through floats would spend it here.
    if (md_attribute_extract_slice_f64(dst, num_ao, attr, slice, md_unit_none()) != num_ao) {
        return nullptr;
    }

    if (out_num_ao) *out_num_ao = num_ao;
    return dst;
}

// ---------------------------------------------------------------------------
// Which system atom each basis atom is
// ---------------------------------------------------------------------------
// A QM calculation may cover only PART of a loaded system - a chromophore inside a protein - and
// then the basis' atom indices are its own and not the system's. That map was the last thing an
// evaluation still needed a loader for, so it is published beside the basis and read from there.
//
// Absent means the identity, which is what a standalone load is, so nothing publishes it in the
// common case and every consumer needs the same one line to handle both.
//
// It is VIAMD that publishes it and not the reader, because the file alone cannot decide: the same
// h5 carries a local-to-global map whether it is opened on its own - where the map must NOT be
// applied, since the system IS the QM atoms - or against a larger system, where it must. Only the
// side holding both knows which, and that is here.
static const str_t BASIS_ATOM_MAP_PATH = STR_LIT("basis/atom/system_index");

void basis_atom_map_publish(md_system_t* sys, const int32_t* system_index, size_t count) {
    ASSERT(sys);

    // Clearing is as important as setting: a stale map from a previous load would silently send
    // every evaluation to the wrong atoms.
    if (const md_attribute_t* existing = md_attributes_find(&sys->attributes, BASIS_ATOM_MAP_PATH)) {
        md_attributes_remove(&sys->attributes, existing->id);
    }
    if (!system_index || count == 0) {
        return;
    }

    const md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_U32, .components = 1, .rank = 1, .shape = { (uint32_t)count },
    };
    const md_attribute_desc_t desc = {
        .path   = BASIS_ATOM_MAP_PATH,
        .format = format,
        .unit   = md_unit_none(),
        .label  = STR_LIT("System Atom Index"),
    };
    md_attribute_id_t id = md_attributes_create(&sys->attributes, &desc);

    uint32_t* dst = (uint32_t*)md_attributes_data(&sys->attributes, id, MD_ATTRIBUTE_TYPE_U32);
    if (!dst) {
        if (id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&sys->attributes, id);
        return;
    }
    for (size_t i = 0; i < count; ++i) {
        dst[i] = (uint32_t)system_index[i];
    }
}

// NULL when there is none, which the callers read as the identity. Same explicitness as
// gto_attr_column in md_gto.c and for the same reason: the table is open, so a path with the wrong
// shape is a mistake to refuse rather than to reinterpret.
static const uint32_t* basis_atom_map_find(size_t* out_count, const md_system_t& sys) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, BASIS_ATOM_MAP_PATH);
    if (!attr) {
        return nullptr;
    }
    if (attr->format.type != MD_ATTRIBUTE_TYPE_U32 || attr->format.rank != 1 || md_attribute_components(&attr->format) != 1 || !attr->data) {
        MD_LOG_ERROR("'" STR_FMT "' is not a plain column of atom indices", STR_ARG(BASIS_ATOM_MAP_PATH));
        return nullptr;
    }
    if (out_count) *out_count = attr->format.shape[0];
    return (const uint32_t*)attr->data;
}

// Gathers the BASIS atoms' positions out of the state, in the Bohr and the interleaved layout
// md_gto works in, in basis atom order. Returns the number written, 0 on failure.
//
// This conversion is the one input to an evaluation which is neither an attribute nor a caller
// parameter, and that is the right shape: the basis deliberately stores no coordinates, so that it
// survives a geometry change and the positions come from wherever the current ones live.
static size_t basis_atom_positions_gather(vec3_t* dst, size_t cap, const md_system_t& sys, const md_system_state_t& state, size_t num_basis_atoms) {
    if (!dst || num_basis_atoms == 0 || num_basis_atoms > cap) {
        return 0;
    }
    // Checked here rather than at each call site: this is the only place the state is read, so this
    // is the only place that can be wrong about it.
    if (state.num_atoms == 0 || !state.x || !state.y || !state.z) {
        return 0;
    }

    size_t map_count = 0;
    const uint32_t* map = basis_atom_map_find(&map_count, sys);
    if (map && map_count < num_basis_atoms) {
        MD_LOG_ERROR("The basis spans %zu atoms and '" STR_FMT "' maps %zu", num_basis_atoms, STR_ARG(BASIS_ATOM_MAP_PATH), map_count);
        return 0;
    }

    for (size_t i = 0; i < num_basis_atoms; ++i) {
        const size_t idx = map ? (size_t)map[i] : i;
        if (idx >= state.num_atoms) {
            MD_LOG_ERROR("Basis atom %zu is system atom %zu and the state holds %zu", i, idx, state.num_atoms);
            return 0;
        }
        dst[i] = vec3_set(state.x[idx], state.y[idx], state.z[idx]) * (float)ANGSTROM_TO_BOHR;
    }
    return num_basis_atoms;
}

// ---------------------------------------------------------------------------
// Volume readback
// ---------------------------------------------------------------------------
// Getting an evaluated volume out of the device scratch texture and into the GL texture a
// representation draws. Nothing about it is specific to what was evaluated or to who asked, and
// every handle it touches is on the ApplicationState above - which is why it lives here and not in
// whichever component happened to need it first.

#if MD_ENABLE_GPU

// Prefers a pixel unpack buffer: one write into memory the GPU already sees, and the transfer
// overlaps instead of blocking. Falls back to the plain client pointer upload when no buffer is
// available.
static void gpu_volume_upload_to_gl(uint32_t vol_tex, const void* src, size_t size) {
    if (!src) return;
    if (void* dst = gl::pbo_upload_begin(size)) {
        MEMCPY(dst, src, size);
        if (gl::pbo_upload_end_texture_3D(vol_tex, 0, GL_R32F)) return;
    }
    gl::set_texture_3D_data(vol_tex, 0, src, GL_R32F);
}

static ApplicationState::GpuVolumeJob* gpu_volume_job_acquire(ApplicationState* state) {
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        if (!state->gpu_volume_jobs[i].in_flight) {
            state->gpu_volume_jobs[i] = ApplicationState::GpuVolumeJob{};
            state->gpu_volume_jobs[i].in_flight = true;
            state->gpu_volume_jobs[i].owner = state;
            return &state->gpu_volume_jobs[i];
        }
    }
    return nullptr;
}

static bool gpu_volume_job_any_in_flight(const ApplicationState* state) {
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        if (state->gpu_volume_jobs[i].in_flight) return true;
    }
    return false;
}

static bool gpu_volume_job_in_flight_for(const ApplicationState* state, uint32_t tex_id) {
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        if (state->gpu_volume_jobs[i].in_flight && state->gpu_volume_jobs[i].tex_id == tex_id) return true;
    }
    return false;
}

// Runs on the GL thread, from md_gpu_device_poll() in the frame loop.
static void gpu_volume_job_complete(void* user) {
    ApplicationState::GpuVolumeJob* job = (ApplicationState::GpuVolumeJob*)user;
    ApplicationState* self = job->owner;
    if (self && job->tex_id) {
        gpu_volume_upload_to_gl(job->tex_id, md_gpu_host_ptr(job->rb), job->size);
    }
    if (self && job->rb) md_gpu_free(job->rb, self->gpu_stream);
    job->rb        = nullptr;
    job->in_flight = false;
}

// Used when no job slot is free, and when a caller genuinely needs the data before it returns.
static bool gpu_volume_readback_blocking(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, size_t size) {
    md_gpu_ptr_t rb = md_gpu_malloc(state->gpu_rb_pool, size, state->gpu_stream);
    if (!rb) return false;
    const md_gpu_tex_region_t region = {
        .offset = {0, 0, 0},
        .extent = { (uint32_t)grid.dim[0], (uint32_t)grid.dim[1], (uint32_t)grid.dim[2] },
    };
    bool ok = md_gpu_memcpy_from_tex_async(rb, state->gpu_volume, &region, size, state->gpu_stream);
    md_gpu_stream_sync(state->gpu_stream);
    if (ok) gpu_volume_upload_to_gl(vol_tex, md_gpu_host_ptr(rb), size);
    md_gpu_free(rb, state->gpu_stream);
    return ok;
}

// Queues a readback of the evaluated region of gpu_volume. Returns once the copy is recorded; the
// GL texture is filled later, from gpu_volume_job_complete(), which md_gpu_device_poll() calls on
// the GL thread.
//
// Returns true when the work was QUEUED -- not that the texture holds data.
static bool gpu_volume_readback(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid) {
    if (!state->gpu_stream || !state->gpu_rb_pool || !state->gpu_volume) {
        return false;
    }
    const size_t size = sizeof(float) * (size_t)grid.dim[0] * (size_t)grid.dim[1] * (size_t)grid.dim[2];

    // Never allow two outstanding readbacks for the same texture. Two reasons, both of which
    // produce a wrong image rather than a slow one:
    //
    //  - a staging block can be most of a gigabyte, so N in flight is N times that;
    //  - md_gpu_stream_sync does not fire user callbacks, only md_gpu_device_poll does. So a
    //    blocking fallback taken while an older job is still pending would upload new data now and
    //    let the older callback overwrite it with stale data next frame.
    //
    // Draining costs the stall we are trying to avoid, but only when the user outruns the GPU, and
    // it keeps uploads strictly ordered.
    if (gpu_volume_job_in_flight_for(state, vol_tex)) {
        gpu_volume_jobs_drain(state);
    }

    ApplicationState::GpuVolumeJob* job = gpu_volume_job_acquire(state);
    if (!job) {
        gpu_volume_jobs_drain(state);
        job = gpu_volume_job_acquire(state);
    }
    if (!job) {
        return gpu_volume_readback_blocking(state, vol_tex, grid, size);
    }

    md_gpu_ptr_t rb = md_gpu_malloc(state->gpu_rb_pool, size, state->gpu_stream);
    if (!rb) {
        MD_LOG_ERROR("Failed to allocate volume readback staging (%zu bytes)", size);
        job->in_flight = false;
        return false;
    }

    const md_gpu_tex_region_t region = {
        .offset = {0, 0, 0},
        .extent = { (uint32_t)grid.dim[0], (uint32_t)grid.dim[1], (uint32_t)grid.dim[2] },
    };
    if (!md_gpu_memcpy_from_tex_async(rb, state->gpu_volume, &region, size, state->gpu_stream)) {
        MD_LOG_ERROR("Failed to record the volume readback");
        md_gpu_free(rb, state->gpu_stream);
        job->in_flight = false;
        return false;
    }

    job->tex_id = vol_tex;
    job->rb     = rb;
    job->size   = size;

    if (!md_gpu_launch_host_fn(state->gpu_stream, gpu_volume_job_complete, job)) {
        // No callback means nothing would ever free rb or fill the texture, so finish this one
        // synchronously instead of leaking it.
        MD_LOG_ERROR("Failed to queue the volume completion; falling back to a blocking readback");
        md_gpu_stream_sync(state->gpu_stream);
        gpu_volume_upload_to_gl(vol_tex, md_gpu_host_ptr(rb), size);
        md_gpu_free(rb, state->gpu_stream);
        job->in_flight = false;
        return true;
    }
    return true;
}

// Queues the atom upload on gpu_stream if it is dirty. Any evaluation launched afterwards on the
// same stream observes it. md_gpu_upload_begin writes straight into the destination when that is
// safe and stages through a transient arena otherwise, so there is one path regardless of whether
// the device is discrete.
static bool gpu_atoms_ensure_uploaded(ApplicationState* state, const md_system_t& sys, const md_system_state_t& sys_state) {
    if (!state->gpu_stream || !state->mold.gpu_basis || !state->mold.gpu_atoms) {
        return false;
    }
    const size_t num_atoms = md_gto_gpu_basis_num_atoms(state->mold.gpu_basis);
    const size_t sz = md_gto_gpu_atom_buffer_size(num_atoms);

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    // Same gather as the GL path, through the same map, so the two backends cannot disagree about
    // which atom a basis index means. Packed to vec4 because that is what the device wants.
    vec3_t* xyz = (vec3_t*)md_temp_alloc(temp, sizeof(vec3_t) * MAX(num_atoms, (size_t)1));
    if (!xyz || basis_atom_positions_gather(xyz, num_atoms, sys, sys_state, num_atoms) != num_atoms) {
        return false;
    }

    // Gather first, then decide: the gather is what the comparison is ABOUT, and it is cheap next
    // to the transfer it may save.
    const uint64_t hash = md_hash64(xyz, sizeof(vec3_t) * num_atoms, 0);
    if (hash != 0 && hash == state->mold.gpu_atoms_hash) {
        return true;
    }

    vec4_t* xyzw = (vec4_t*)md_temp_alloc(temp, sizeof(vec4_t) * MAX(num_atoms, (size_t)1));
    if (!xyzw) {
        return false;
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        xyzw[i] = vec4_set(xyz[i].x, xyz[i].y, xyz[i].z, 1.0f);
    }

    float* dst = (float*)md_gpu_upload_begin(state->gpu_stream, state->mold.gpu_atoms, sz);
    if (!dst) {
        MD_LOG_ERROR("Failed to upload the atom positions to the device");
        return false;
    }
    md_gto_gpu_atom_pack(dst, (const float*)xyzw, sizeof(vec4_t), num_atoms);
    if (!md_gpu_upload_end(state->gpu_stream)) {
        return false;
    }

    state->mold.gpu_atoms_hash = hash;
    return true;
}

#endif // MD_ENABLE_GPU

void gpu_volume_jobs_drain(ApplicationState* state) {
    ASSERT(state);
#if MD_ENABLE_GPU
    if (!state->gpu_stream || !gpu_volume_job_any_in_flight(state)) return;
    md_gpu_stream_sync(state->gpu_stream);
    md_gpu_device_poll(state->gpu_device);   // this is what actually runs the callbacks

    // Backstop: if a callback somehow did not fire, release the staging block here rather than
    // leaking it, since the pool may be about to go.
    for (int i = 0; i < ApplicationState::GPU_VOLUME_JOB_SLOTS; ++i) {
        ApplicationState::GpuVolumeJob& j = state->gpu_volume_jobs[i];
        if (j.in_flight) {
            if (j.rb) md_gpu_free(j.rb, state->gpu_stream);
            j.rb = nullptr;
            j.in_flight = false;
        }
    }
#else
    (void)state;
#endif
}

// The two things every GTO evaluation needs besides its own coefficients: the basis, rebuilt from
// the system's basis/ attributes, and the basis atoms' positions in the units md_gto works in.
// Neither depends on WHAT is being evaluated, and both come out of the system - so this is the
// whole of the context an evaluation has, and there is deliberately nothing else in it.
static bool gto_evaluation_context(md_gto_basis_t* out_basis, vec3_t** out_atom_xyz, md_temp_scope_t temp,
                                   const md_system_t& sys, const md_system_state_t& state) {
    ASSERT(out_basis);
    ASSERT(out_atom_xyz);

    MEMSET(out_basis, 0, sizeof(*out_basis));
    if (!md_gto_basis_extract_attributes(out_basis, &sys.attributes, md_temp_allocator(temp))) {
        return false;
    }

    const size_t num_basis_atoms = md_gto_basis_num_atoms(out_basis);
    vec3_t* atom_xyz = (vec3_t*)md_temp_alloc(temp, sizeof(vec3_t) * MAX(num_basis_atoms, (size_t)1));
    if (!atom_xyz || basis_atom_positions_gather(atom_xyz, num_basis_atoms, sys, state, num_basis_atoms) != num_basis_atoms) {
        return false;
    }

    *out_atom_xyz = atom_xyz;
    return true;
}

bool orbital_evaluate_gl(uint32_t vol_tex, const md_grid_t& grid, const md_system_t& sys, const md_system_state_t& state,
                         str_t coefficient_path, const md_attribute_slice_t* slice, md_gto_eval_mode_t mode, md_gto_op_t op, double cutoff) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    size_t  num_ao    = 0;
    double* ao_coeffs = orbital_coefficients_extract(&num_ao, temp, sys, coefficient_path, slice);
    if (!ao_coeffs) {
        return false;
    }

    md_gto_basis_t basis = {};
    vec3_t* atom_xyz = nullptr;
    if (!gto_evaluation_context(&basis, &atom_xyz, temp, sys, state)) {
        return false;
    }

    // The coefficients and the basis are separate attributes, so nothing guarantees they agree
    // until it is checked here.
    if (md_gto_basis_num_ao(&basis) != num_ao) {
        MD_LOG_ERROR("The basis spans %zu atomic orbitals and '" STR_FMT "' %zu", md_gto_basis_num_ao(&basis), STR_ARG(coefficient_path), num_ao);
        return false;
    }

    md_gto_grid_evaluate_mo_GL(vol_tex, &grid, &basis, (const float*)atom_xyz, sizeof(vec3_t), ao_coeffs, cutoff, mode, op);
    return true;
}

// ---------------------------------------------------------------------------
// Density evaluation
// ---------------------------------------------------------------------------
// A density in this table is an AO x AO matrix. It is either that on its own - the SCF densities,
// the density properties - or the innermost two axes of something indexed by state, as the
// transition densities are at {S,A,A}. Both are addressed the same way here: a path plus a slice
// which narrows to one matrix. The caller names what it wants and this cares neither which of the
// two shapes it came from, nor whether the values were stored or worked out on demand, because
// md_attribute_slice_format answers the first and the extract answers the second.
double* density_matrix_extract(size_t* out_dim, md_temp_scope_t temp, const md_system_t& sys, str_t density_path, const md_attribute_slice_t* slice) {
    const md_attribute_t* attr = md_attributes_find(&sys.attributes, density_path);
    if (!attr) {
        MD_LOG_DEBUG("No density published at '" STR_FMT "'", STR_ARG(density_path));
        return nullptr;
    }

    md_attribute_format_t format = {};
    if (!md_attribute_slice_format(&format, attr, slice)) {
        MD_LOG_ERROR("The slice does not address '" STR_FMT "'", STR_ARG(density_path));
        return nullptr;
    }

    // Square is not pedantry: the GL and GPU density paths both pack the upper triangle and never
    // read the lower half, so a non square matrix would be silently half consumed.
    if (format.rank != 2 || format.shape[0] != format.shape[1] || md_attribute_components(&format) != 1) {
        MD_LOG_ERROR("'" STR_FMT "' does not slice down to a square density matrix", STR_ARG(density_path));
        return nullptr;
    }

    const size_t count = md_attribute_slice_count(attr, slice);
    if (count == 0) {
        return nullptr;
    }

    double* dst = (double*)md_temp_alloc(temp, sizeof(double) * count);
    if (!dst) {
        return nullptr;
    }
    if (md_attribute_extract_slice_f64(dst, count, attr, slice, md_unit_none()) != count) {
        return nullptr;
    }

    if (out_dim) *out_dim = format.shape[0];
    return dst;
}

bool density_evaluate_gl(uint32_t vol_tex, const md_grid_t& grid, const md_system_t& sys, const md_system_state_t& state,
                         str_t density_path, const md_attribute_slice_t* slice, md_gto_op_t op) {
    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    size_t  dim = 0;
    double* density_matrix = density_matrix_extract(&dim, temp, sys, density_path, slice);
    if (!density_matrix) {
        return false;
    }

    md_gto_basis_t basis = {};
    vec3_t* atom_xyz = nullptr;
    if (!gto_evaluation_context(&basis, &atom_xyz, temp, sys, state)) {
        return false;
    }

    if (md_gto_basis_num_ao(&basis) != dim) {
        MD_LOG_ERROR("The basis spans %zu atomic orbitals and '" STR_FMT "' %zu", md_gto_basis_num_ao(&basis), STR_ARG(density_path), dim);
        return false;
    }

    md_gto_grid_evaluate_density_GL(vol_tex, &grid, &basis, (const float*)atom_xyz, sizeof(vec3_t), density_matrix, false, op);
    return true;
}

static bool es_attribute_exists(const md_system_t& sys, str_t path) {
    return md_attributes_find(&sys.attributes, path) != nullptr;
}

// Evaluates whatever the representation is currently pointed at into its own volume texture.
// Everything it needs is an attribute on the system or a field of the representation; there is no
// reader, no component and no event in the path, which is the whole point of the port.
static bool electronic_structure_evaluate(ApplicationState* state, Representation* rep) {
    const ElectronicStructureRepresentation& es = rep->electronic_structure;

    const double samples_per_angstrom = volume_resolution_samples_per_angstrom[(int)es.resolution];

    md_grid_t grid = {};
    if (!electronic_structure_grid_init(&grid, state->mold.sys, state->mold.state, samples_per_angstrom * BOHR_TO_ANGSTROM)) {
        return false;
    }
    init_volume(&rep->electronic_structure.density_vol, grid, GL_R32F);
    const uint32_t tex_id = rep->electronic_structure.density_vol.tex_id;

    switch (es.source) {
    case ElectronicStructureSource::MolecularOrbital: {
        const str_t path = (es.spin == ElectronicStructureSpin::Beta) ? es_path::beta_coefficient : es_path::alpha_coefficient;
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)es.orbital_idx);
        return orbital_evaluate(state, tex_id, grid, path, &slice, MD_GTO_EVAL_MODE_PSI, es_gto_op(es.use_magnitude), DEFAULT_GTO_CUTOFF_VALUE);
    }
    case ElectronicStructureSource::NaturalTransitionOrbital: {
        const str_t path = (es.nto_component == ElectronicStructureNtoComponent::Particle) ? es_path::nto_particle : es_path::nto_hole;
        // {S,L,A}: the excited state and then the lambda pair within it. Two indices instead of the
        // one an ordinary orbital takes, which is exactly what a slice is for.
        const md_attribute_slice_t slice = md_attribute_slice_2((uint32_t)es.excited_state_idx, (uint32_t)es.nto_lambda_idx);
        return orbital_evaluate(state, tex_id, grid, path, &slice, MD_GTO_EVAL_MODE_PSI, es_gto_op(es.use_magnitude), DEFAULT_GTO_CUTOFF_VALUE);
    }
    case ElectronicStructureSource::TransitionDensity: {
        str_t path = es_path::attachment_density;
        switch (es.transition_density_component) {
        case ElectronicStructureTransitionDensityComponent::Detachment: path = es_path::detachment_density; break;
        case ElectronicStructureTransitionDensityComponent::Difference: path = es_path::transition_diff;    break;
        case ElectronicStructureTransitionDensityComponent::Attachment:
        default: break;
        }
        // These are VIRTUAL: the slice is what tells the provider to reconstruct one state rather
        // than all of them, and it is the only reason asking for one is affordable.
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)es.excited_state_idx);
        return density_evaluate(state, tex_id, grid, path, &slice, MD_GTO_OP_SET);
    }
    case ElectronicStructureSource::ElectronDensity: {
        // A spin difference is signed, so it is the one density the magnitude toggle applies to.
        const bool magnitude = (es.spin == ElectronicStructureSpin::Difference) && es.use_magnitude;
        return density_evaluate(state, tex_id, grid, es_electron_density_path(es.spin), nullptr, es_gto_op(magnitude));
    }
    case ElectronicStructureSource::DensityProperty: {
        const md_attribute_t* attr = md_attributes_get(&state->mold.sys.attributes, es.density_property_key);
        if (!attr) {
            return false;
        }
        return density_evaluate(state, tex_id, grid, attr->path, nullptr, MD_GTO_OP_SET);
    }
    default:
        MD_LOG_ERROR("Unknown electronic structure source");
        return false;
    }
}

// The colour volume beside the density one: the atoms' own colours splatted into a downsampled 3D
// texture, so that a surface can be shaded by whatever the representation is colouring atoms with.
static void electronic_structure_color_volume_update(ApplicationState* state, Representation* rep, const uint32_t* atom_colors) {
    if (!atom_colors) {
        return;
    }

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    const md_system_t&       sys       = state->mold.sys;
    const md_system_state_t& sys_state = state->mold.state;

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, md_temp_allocator(temp))) {
        return;
    }
    const size_t num_points = md_gto_basis_num_atoms(&basis);
    if (num_points == 0) {
        return;
    }

    const int downsample_factor = 1;
    int dim[3] = {
        (int)(rep->electronic_structure.density_vol.dim[0] / downsample_factor),
        (int)(rep->electronic_structure.density_vol.dim[1] / downsample_factor),
        (int)(rep->electronic_structure.density_vol.dim[2] / downsample_factor),
    };
    MEMCPY(rep->electronic_structure.color_vol.dim, dim, sizeof(dim));
    rep->electronic_structure.color_vol.world_to_model   = rep->electronic_structure.density_vol.world_to_model;
    rep->electronic_structure.color_vol.texture_to_world = rep->electronic_structure.density_vol.texture_to_world;
    rep->electronic_structure.color_vol.voxel_size       = rep->electronic_structure.density_vol.voxel_size * (float)downsample_factor;
    gl::init_texture_3D(&rep->electronic_structure.color_vol.tex_id, dim[0], dim[1], dim[2], GL_RGBA8);

    const vec3_t& voxel_size     = rep->electronic_structure.color_vol.voxel_size;
    const mat4_t& world_to_model = rep->electronic_structure.color_vol.world_to_model;
    mat4_t index_to_world = rep->electronic_structure.color_vol.texture_to_world
        * mat4_scale(1.0f / dim[0], 1.0f / dim[1], 1.0f / dim[2])
        * mat4_translate(0.5f, 0.5f, 0.5f); // Center of the corner voxel should be at the origin

    // The colours are per SYSTEM atom and the splats are per BASIS atom, so the map is consulted
    // here too - and once more it is the same lookup in the same direction the evaluation used.
    size_t map_count = 0;
    const uint32_t* map = basis_atom_map_find(&map_count, sys);
    if (map && map_count < num_points) {
        return;
    }

    vec4_t*   point_xyzw   = (vec4_t*)md_temp_alloc(temp, sizeof(vec4_t) * num_points);
    uint32_t* point_colors = (uint32_t*)md_temp_alloc(temp, sizeof(uint32_t) * num_points);
    if (!point_xyzw || !point_colors) {
        return;
    }
    if (sys_state.num_atoms == 0 || !sys_state.x || !sys_state.y || !sys_state.z) {
        return;
    }
    for (size_t i = 0; i < num_points; ++i) {
        const size_t idx = map ? (size_t)map[i] : i;
        if (idx >= sys_state.num_atoms) {
            return;
        }
        // Angstrom and not Bohr: these are world coordinates for the splatting pass, which shares
        // the Volume's transforms, unlike the grid the density was evaluated on.
        const float radius = md_atom_radius(&sys.atom, idx);
        point_xyzw[i]   = vec4_set(sys_state.x[idx], sys_state.y[idx], sys_state.z[idx], radius);
        point_colors[i] = atom_colors[idx];
    }

    volume::compute_point_color_volume(rep->electronic_structure.color_vol.tex_id, dim, voxel_size.elem, world_to_model.elem,
                                       index_to_world.elem, point_xyzw, point_colors, num_points,
                                       rep->electronic_structure.gaussian_splatting_power);
}

// ---------------------------------------------------------------------------
// Where to evaluate
// ---------------------------------------------------------------------------
// An object aligned box around the BASIS atoms, padded, at the requested sample density. Derived
// from the system on every call rather than cached at load, for two reasons: it is a PCA over the
// tens to hundreds of atoms a QM calculation covers, which is nothing next to the evaluation it
// precedes, and computing it here is what makes the box follow the geometry instead of pinning it
// to whatever the coordinates were when the file was opened.
bool electronic_structure_grid_init(md_grid_t* grid, const md_system_t& sys, const md_system_state_t& state, double samples_per_unit_length) {
    ASSERT(grid);

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, md_temp_allocator(temp))) {
        return false;
    }

    const size_t num_atoms = md_gto_basis_num_atoms(&basis);
    if (num_atoms == 0) {
        return false;
    }

    vec3_t* xyz = (vec3_t*)md_temp_alloc(temp, sizeof(vec3_t) * num_atoms);
    if (!xyz || basis_atom_positions_gather(xyz, num_atoms, sys, state, num_atoms) != num_atoms) {
        return false;
    }

    // mat3_PCA and calculate_bounds both want a homogeneous point, and w == 1 is what makes the
    // rotation in calculate_bounds a rotation of a POINT rather than of a direction.
    vec4_t* xyzw = (vec4_t*)md_temp_alloc(temp, sizeof(vec4_t) * num_atoms);
    if (!xyzw) {
        return false;
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        xyzw[i] = vec4_set(xyz[i].x, xyz[i].y, xyz[i].z, 1.0f);
    }

    OABB oabb = {};
    oabb.orientation = mat3_PCA(xyzw, num_atoms);
    calculate_bounds(oabb.min_ext.elem, oabb.max_ext.elem, xyzw, num_atoms, oabb.orientation);

    init_grid(grid, oabb.orientation, oabb.min_ext, oabb.max_ext, samples_per_unit_length);
    return true;
}

// ---------------------------------------------------------------------------
// Evaluation: the entry points a representation calls
// ---------------------------------------------------------------------------
// One function per KIND of thing being evaluated - an orbital, a density - and the choice of
// backend inside it. A caller names the attribute it wants drawn and gets a filled texture; which
// device did the work is not its business, and pushing that choice down here is what lets the two
// backends share the basis, the atom gather and the map they all agree on.
//
// The GPU path writes the device scratch volume and queues a readback into vol_tex; the texture is
// filled from the frame loop's md_gpu_device_poll rather than before this returns. Call
// gpu_volume_jobs_drain() where the contents are needed immediately.

#if MD_ENABLE_GPU
// True when the device scratch and this dataset's uploaded basis are both present, which is what
// every GPU path below needs before it can start.
static bool gpu_evaluation_ready(const ApplicationState* state) {
    return state->gpu_stream && state->mold.gpu_basis && state->mold.gpu_atoms && state->gpu_coeff && state->gpu_volume;
}
#endif

bool orbital_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, str_t coefficient_path,
                      const md_attribute_slice_t* slice, md_gto_eval_mode_t mode, md_gto_op_t op, double cutoff) {
    ASSERT(state);

#if MD_ENABLE_GPU
    if (gpu_evaluation_ready(state)) {
        const md_system_t&       sys       = state->mold.sys;
        const md_system_state_t& sys_state = state->mold.state;

        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        size_t  num_ao    = 0;
        const double* ao_coeffs = orbital_coefficients_extract(&num_ao, temp, sys, coefficient_path, slice);
        if (!ao_coeffs) {
            return false;
        }

        const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
        if (num_cgtos != num_ao) {
            MD_LOG_ERROR("The uploaded basis spans %zu atomic orbitals and '" STR_FMT "' %zu", num_cgtos, STR_ARG(coefficient_path), num_ao);
            return false;
        }
        if (!gpu_atoms_ensure_uploaded(state, sys, sys_state)) {
            return false;
        }

        const double* coeff_ptrs[1] = { ao_coeffs };
        float* dst = (float*)md_gpu_upload_begin(state->gpu_stream, state->gpu_coeff, md_gto_gpu_coeff_size_mo(1, num_cgtos));
        if (!dst) {
            return false;
        }
        md_gto_gpu_coeff_pack_mo(dst, coeff_ptrs, nullptr, 1, num_cgtos);
        md_gpu_upload_end(state->gpu_stream);

        md_gto_gpu_orbital_desc_t desc = {
            .basis        = state->mold.gpu_basis,
            .atom_xyz     = state->mold.gpu_atoms,
            .coeff        = state->gpu_coeff,
            .out_tex      = state->gpu_volume,
            .grid         = &grid,
            .sample_offset = {0.5f, 0.5f, 0.5f},
            .num_orbitals = 1,
            .eval_mode    = mode,
            .op           = op,
        };
        md_gto_gpu_orbital_launch(state->gpu_stream, &desc);
        return gpu_volume_readback(state, vol_tex, grid);
    }
#endif

    return orbital_evaluate_gl(vol_tex, grid, state->mold.sys, state->mold.state, coefficient_path, slice, mode, op, cutoff);
}

#if MD_ENABLE_GPU
bool density_evaluate_to_gpu_volume(ApplicationState* state, const md_grid_t& grid, str_t density_path,
                                    const md_attribute_slice_t* slice, md_gto_op_t op) {
    ASSERT(state);
    if (!gpu_evaluation_ready(state)) {
        return false;
    }

    const md_system_t&       sys       = state->mold.sys;
    const md_system_state_t& sys_state = state->mold.state;

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    size_t dim = 0;
    const double* density_matrix = density_matrix_extract(&dim, temp, sys, density_path, slice);
    if (!density_matrix) {
        return false;
    }

    const size_t num_cgtos = md_gto_gpu_basis_num_cgtos(state->mold.gpu_basis);
    if (num_cgtos != dim) {
        MD_LOG_ERROR("The uploaded basis spans %zu atomic orbitals and '" STR_FMT "' %zu", num_cgtos, STR_ARG(density_path), dim);
        return false;
    }
    if (!gpu_atoms_ensure_uploaded(state, sys, sys_state)) {
        return false;
    }

    float* dst = (float*)md_gpu_upload_begin(state->gpu_stream, state->gpu_coeff, md_gto_gpu_coeff_size_density(num_cgtos));
    if (!dst) {
        return false;
    }
    md_gto_gpu_coeff_pack_density(dst, density_matrix, num_cgtos);
    md_gpu_upload_end(state->gpu_stream);

    md_gto_gpu_density_desc_t desc = {
        .basis         = state->mold.gpu_basis,
        .atom_xyz      = state->mold.gpu_atoms,
        .coeff         = state->gpu_coeff,
        .out_tex       = state->gpu_volume,
        .grid          = &grid,
        .sample_offset = {0.5f, 0.5f, 0.5f},
        .op            = op,
    };
    md_gto_gpu_density_launch(state->gpu_stream, &desc);
    return true;
}
#endif

bool density_evaluate(ApplicationState* state, uint32_t vol_tex, const md_grid_t& grid, str_t density_path,
                      const md_attribute_slice_t* slice, md_gto_op_t op) {
    ASSERT(state);

#if MD_ENABLE_GPU
    // Evaluate, then read back. The split exists because a consumer that runs another kernel over
    // the result - the critical point extraction does - wants the first half and not the second.
    if (density_evaluate_to_gpu_volume(state, grid, density_path, slice, op)) {
        return gpu_volume_readback(state, vol_tex, grid);
    }
    if (gpu_evaluation_ready(state)) {
        return false;   // the GPU path was available and failed; the GL one would fail the same way
    }
#endif

    return density_evaluate_gl(vol_tex, grid, state->mold.sys, state->mold.state, density_path, slice, op);
}

void update_all_representations(ApplicationState* state) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        update_representation(state, &state->representation.reps[i]);
    }
}

bool representation_uses_atom_colors(const Representation& rep) {
    switch (rep.type) {
        case RepresentationType::ElectronicStructure:
            return rep.electronic_structure.use_atom_colors;
        case RepresentationType::DipoleMoment:
            return false;
        default:
            return true;
    }
}

void update_representation(ApplicationState* state, Representation* rep) {
    ASSERT(state);
    ASSERT(rep);

    if (!rep->enabled) return;
    if (!rep->needs_update) return;

    const auto& sys = state->mold.sys;
    size_t num_atoms = md_system_atom_count(&sys);

    md_allocator_i* frame_alloc = state->allocator.frame;
    md_temp_scope_t temp = md_temp_begin_in(frame_alloc);
    defer { md_temp_end(temp); };

    const size_t bytes = num_atoms * sizeof(uint32_t);

    //md_script_property_t prop = {0};
    //if (rep->color_mapping == ColorMapping::Property) {
    //rep->prop_is_valid = md_script_compile_and_eval_property(&prop, rep->prop, &data->mold.sys, frame_allocator, &data->script.ir, rep->prop_error.beg(), rep->prop_error.capacity());
    //}

    uint32_t* colors = 0;
    if (representation_uses_atom_colors(*rep)) {
        colors = (uint32_t*)md_vm_arena_push(frame_alloc, sizeof(uint32_t) * num_atoms);

        switch (rep->color_mapping) {
        case ColorMapping::Uniform:
            color_atoms_uniform(colors, num_atoms, convert_color(rep->base_color));
            break;
        case ColorMapping::Type:
            color_atoms_type(colors, num_atoms, sys);
            break;
        case ColorMapping::Serial:
            color_atoms_idx(colors, num_atoms, sys);
            break;
        case ColorMapping::CompName:
            color_atoms_comp_name(colors, num_atoms, sys);
            break;
        case ColorMapping::CompSeqId:
            color_atoms_comp_seq_id(colors, num_atoms, sys);
            break;
        case ColorMapping::CompIndex:
            color_atoms_comp_idx(colors, num_atoms, sys);
            break;
        case ColorMapping::InstId:
            color_atoms_inst_id(colors, num_atoms, sys);
            break;
        case ColorMapping::InstIndex:
            color_atoms_inst_idx(colors, num_atoms, sys);
            break;
        case ColorMapping::SecondaryStructure: {
            SecondaryStructurePalette palette = {
                .coil  = convert_color(rep->secondary_structure.color_coil),
                .helix = convert_color(rep->secondary_structure.color_helix),
                .sheet = convert_color(rep->secondary_structure.color_sheet),
            };
            color_atoms_secondary_structure(colors, num_atoms, sys, palette);
            break;
        }
        case ColorMapping::Property:
            // @TODO: Map colors accordingly
            //color_atoms_uniform(colors, mol.atom.count, rep->uniform_color);

            {
                const md_attribute_t* attr = md_attributes_get(&sys.attributes, rep->atomic_property.key);
                size_t num_extracted = 0;
                float* values = nullptr;

                if (attr) {
                    values = (float*)md_vm_arena_push(frame_alloc, sizeof(float) * num_atoms);

                    // Fixing the variant axis hands back exactly the atom axis, so there is no
                    // offset arithmetic here to get wrong. A field with no variant axis is rank 1
                    // and takes no indices at all.
                    const uint32_t variant = (uint32_t)CLAMP(rep->atomic_property.variant_idx, 0, MAX(atom_property_variant_count(attr) - 1, 0));
                    const md_attribute_slice_t slice = attr->format.rank > 1 ? md_attribute_slice_1(variant) : md_attribute_slice_all();
                    num_extracted = md_attribute_extract_slice_f32(values, num_atoms, attr, &slice, md_unit_none());
                }

                if (num_extracted == num_atoms) {
                    float range_ext = (rep->atomic_property.range_end - rep->atomic_property.range_beg);
                    range_ext = MAX(range_ext, 0.001f);
                    for (size_t i = 0; i < num_atoms; ++i) {
                        float t = (values[i] - rep->atomic_property.range_beg) / range_ext;
                        colors[i] = ImPlot::SampleColormapU32(ImClamp(t, 0.0f, 1.0f), rep->atomic_property.colormap);
                    }
                } else {
                    if (attr) {
                        MD_LOG_DEBUG("Failed to extract values for the selected atom property");
                    }
                    MEMSET(colors, 0xFFFFFFFFu, bytes);
                }
            }
#if 0
            if (rep->prop) {
                MEMSET(colors, 0xFFFFFFFF, bytes);
                md_script_pro
                    const float* values = rep->prop->data.values;
                if (rep->prop->data.aggregate) {
                    const int dim = rep->prop->data.dim[0];
                    md_script_vis_t vis = {0};
                    bool result = false;

                    //if (md_semaphore_aquire(&data->script.ir_semaphore)) {
                    //    defer { md_semaphore_release(&data->script.ir_semaphore); };

                    if (md_script_ir_valid(state->script.eval_ir)) {
                        md_script_vis_init(&vis, frame_alloc);
                        md_script_vis_ctx_t ctx = {
                            .ir = state->script.eval_ir,
                            .mol = &state->mold.sys,
                            .traj = state->mold.sys.trajectory,
                        };
                        result = md_script_vis_eval_payload(&vis, rep->prop->vis_payload, 0, &ctx, MD_SCRIPT_VISUALIZE_ATOMS);
                    }
                    //}
                    if (result) {
                        if (dim == (int)md_array_size(vis.structure)) {
                            int i0 = CLAMP((int)state->animation.frame + 0, 0, (int)rep->prop->data.num_values / dim - 1);
                            int i1 = CLAMP((int)state->animation.frame + 1, 0, (int)rep->prop->data.num_values / dim - 1);
                            float frame_fract = fractf((float)state->animation.frame);

                            md_bitfield_t mask = {0};
                            md_bitfield_init(&mask, frame_alloc);
                            for (int i = 0; i < dim; ++i) {
                                md_bitfield_and(&mask, &rep->atom_mask, &vis.structure[i]);
                                float value = lerpf(values[i0 * dim + i], values[i1 * dim + i], frame_fract);
                                float t = CLAMP((value - rep->map_beg) / (rep->map_end - rep->map_beg), 0, 1);
                                ImVec4 color = ImPlot::SampleColormap(t, rep->color_map);
                                color_atoms_uniform(colors, mol.atom.count, vec_cast(color), &mask);
                            }
                        }
                    }
                } else {
                    int i0 = CLAMP((int)state->animation.frame + 0, 0, (int)rep->prop->data.num_values - 1);
                    int i1 = CLAMP((int)state->animation.frame + 1, 0, (int)rep->prop->data.num_values - 1);
                    float value = lerpf(values[i0], values[i1], fractf((float)state->animation.frame));
                    float t = CLAMP((value - rep->map_beg) / (rep->map_end - rep->map_beg), 0, 1);
                    ImVec4 color = ImPlot::SampleColormap(t, rep->color_map);
                    color_atoms_uniform(colors, mol.atom.count, vec_cast(color));
                }
            } else {
                color_atoms_uniform(colors, mol.atom.count, rep->uniform_color);
            }
#endif
            break;
        default:
            ASSERT(false);
            break;
        }
    }

    if (colors && (rep->tint_scale > 0.0f || rep->saturation < 1.0f)) {
        uint32_t tint_color = convert_color(rep->tint_color);
        tint_colors(colors, num_atoms, tint_color, rep->tint_scale, rep->saturation);
    }

    switch (rep->type) {
    case RepresentationType::SpaceFill:
        rep->type_is_valid = sys.atom.count > 0;
        break;
    case RepresentationType::Licorice:
        rep->type_is_valid = sys.bond.count > 0;
        break;
    case RepresentationType::BallAndStick:
        rep->type_is_valid = sys.atom.count > 0;
        break;
    case RepresentationType::Ribbons:
    case RepresentationType::Cartoon:
        rep->type_is_valid = sys.protein_backbone.range.count > 0;
        break;
    case RepresentationType::ElectronicStructure: {
        rep->type_is_valid = electronic_structure_source_supported(state->representation.info.electronic_structure_source_mask, rep->electronic_structure.source);
        if (rep->type_is_valid && rep->enabled) {
            // Re-evaluating is expensive and everything it depends on is right here, so it is
            // gated on a hash of exactly those inputs - the frame included, since the geometry
            // moves the grid.
            uint64_t vol_hash = md_hash64(&state->animation.frame, sizeof(state->animation.frame), 0);
            vol_hash = md_hash64(&rep->electronic_structure.source,                       sizeof(rep->electronic_structure.source), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.use_magnitude,                sizeof(rep->electronic_structure.use_magnitude), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.spin,                         sizeof(rep->electronic_structure.spin), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.nto_component,                sizeof(rep->electronic_structure.nto_component), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.transition_density_component, sizeof(rep->electronic_structure.transition_density_component), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.resolution,                   sizeof(rep->electronic_structure.resolution), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.orbital_idx,                  sizeof(rep->electronic_structure.orbital_idx), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.excited_state_idx,            sizeof(rep->electronic_structure.excited_state_idx), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.nto_lambda_idx,               sizeof(rep->electronic_structure.nto_lambda_idx), vol_hash);
            vol_hash = md_hash64(&rep->electronic_structure.density_property_key,         sizeof(rep->electronic_structure.density_property_key), vol_hash);

            if (vol_hash != rep->electronic_structure.vol_hash) {
                rep->electronic_structure.vol_hash = vol_hash;
                electronic_structure_evaluate(state, rep);
            }

            if (rep->electronic_structure.use_atom_colors) {
                uint64_t col_hash = md_hash64(&rep->electronic_structure.gaussian_splatting_power, sizeof(rep->electronic_structure.gaussian_splatting_power), (int)rep->color_mapping);
                col_hash = md_hash64_combine(col_hash, rep->electronic_structure.vol_hash);
                if (col_hash != rep->electronic_structure.col_hash) {
                    rep->electronic_structure.col_hash = col_hash;
                    electronic_structure_color_volume_update(state, rep, colors);
                }
            }
        }
        break;
    }
	case RepresentationType::DipoleMoment:
		rep->type_is_valid = md_attributes_get(&state->mold.sys.attributes, rep->dipole.dipole_key) != NULL;
		break;
    default:
        ASSERT(false);
        break;
    }

    if (colors) {
        if (rep->dynamic_evaluation) {
            rep->filt_is_dirty = true;
        }

        if (rep->filt_is_dirty) {
			rep->filt_is_valid = md_filter(&rep->atom_mask, str_from_cstr(rep->filt), &state->mold.sys, &state->mold.state, state->script.ir, &rep->filt_is_dynamic, rep->filt_error, sizeof(rep->filt_error));
            rep->filt_is_dirty = false;
        }

        if (rep->filt_is_valid) {
            filter_colors(colors, num_atoms, &rep->atom_mask);
            state->representation.atom_visibility_mask_dirty = true;
            md_gl_rep_set_atom_colors(rep->md_rep, 0, (uint32_t)num_atoms, colors, 0);

#if EXPERIMENTAL_GFX_API
            md_gfx_rep_attr_t attributes = {};
            attributes.spacefill.radius_scale = 1.0f;
            md_gfx_rep_set_type_and_attr(rep->gfx_rep, MD_GFX_REP_TYPE_SPACEFILL, &attributes);
            md_gfx_rep_set_color(rep->gfx_rep, 0, (uint32_t)mol.atom.count, (md_gfx_color_t*)colors, 0);
#endif
        }
    }

    rep->needs_update = false;
}

void update_representation_info(ApplicationState* state) {
    md_allocator_i* alloc = state->representation.info.alloc;
    md_arena_allocator_reset(alloc);

    // Clear info, maintain allocator
    state->representation.info = {};
    state->representation.info.alloc = alloc;

    // Broadcast event to populate info
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRepresentationInfoFill, viamd::EventPayloadType_RepresentationInfo, &state->representation.info);

    // What an electronic structure representation can SHOW is decided here and not by whoever
    // loaded the file: each source needs particular attributes, and asking the table which of them
    // exist is both the exact test and one that a second reader satisfies for free.
    const md_system_t& sys = state->mold.sys;
    ElectronicStructureSourceFlags& mask = state->representation.info.electronic_structure_source_mask;

    if (es_attribute_exists(sys, es_path::alpha_coefficient)) mask |= ElectronicStructureSourceFlag_MolecularOrbital;
    if (es_attribute_exists(sys, es_path::alpha_density))     mask |= ElectronicStructureSourceFlag_ElectronDensity;
    if (es_attribute_exists(sys, es_path::nto_particle))      mask |= ElectronicStructureSourceFlag_NaturalTransitionOrbital;
    if (es_attribute_exists(sys, es_path::attachment_density))mask |= ElectronicStructureSourceFlag_TransitionDensity;

    // Every leaf under the density property group is one property. The list is the table's, in the
    // table's own path order, and the label falls back to the last path segment when the file gave
    // the attribute none of its own.
    {
        md_temp_scope_t temp = md_temp_begin();
        defer { md_temp_end(temp); };

        const size_t num_props = md_attributes_query(nullptr, 0, &sys.attributes, es_path::density_property);
        md_attribute_id_t* ids = num_props ? (md_attribute_id_t*)md_temp_alloc(temp, sizeof(md_attribute_id_t) * num_props) : nullptr;
        if (ids) {
            md_attributes_query(ids, num_props, &sys.attributes, es_path::density_property);
            for (size_t i = 0; i < num_props; ++i) {
                const md_attribute_t* attr = md_attributes_get(&sys.attributes, ids[i]);
                if (!attr) continue;

                str_t label = attr->label;
                if (str_empty(label)) {
                    label = attr->path;
                    size_t sep = 0;
                    if (str_rfind_char(&sep, label, '/')) {
                        label = str_substr(label, sep + 1, SIZE_MAX);
                    }
                }
                DensityProperty prop = { .key = attr->id, .label = str_copy(label, alloc) };
                md_array_push(state->representation.info.density_properties, prop, alloc);
            }
            if (md_array_size(state->representation.info.density_properties) > 0) {
                mask |= ElectronicStructureSourceFlag_DensityProperty;
            }
        }
    }

    // Per atom scalar fields are deliberately NOT collected here. They live in the system's
    // attribute table under atom/, whoever loaded the data put them there, and the UI reads that
    // table directly through atom_property_query. A copy kept here could only go stale against it.
}

static void init_all_representations(ApplicationState* state) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        auto& rep = state->representation.reps[i];
        init_representation(state, &rep);
    }
}

void flag_representation_as_dirty(Representation* rep) {
    ASSERT(rep);
    rep->filt_is_dirty = true;
    rep->needs_update  = true;
    rep->electronic_structure.vol_hash = 0;
}

void flag_all_representations_as_dirty(ApplicationState* state) {
    ASSERT(state);
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        flag_representation_as_dirty(&state->representation.reps[i]);
    }
}

void remove_all_representations(ApplicationState* state) {
    while (md_array_size(state->representation.reps) > 0) {
        remove_representation(state, (int32_t)md_array_size(state->representation.reps) - 1);
    }
}


void create_default_representations(ApplicationState* state) {
    bool amino_acid_present = false;
    bool nucleic_present = false;
    bool ion_present = false;
    bool water_present = false;
    bool ligand_present = false;
    bool coarse_grained = false;
    bool orbitals_present = state->representation.info.alpha.num_orbitals > 0;

    if (state->mold.sys.atom.count > 3'000'000) {
        VIAMD_LOG_INFO("Large system detected, creating default representation for all atoms");
        Representation* rep = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    if (state->mold.sys.component.count == 0) {
        // No residues present
        Representation* rep = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Type, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    // TODO: Redo this check with entities instead of atom flags
    for (size_t i = 0; i < state->mold.sys.atom.count; ++i) {
        uint32_t flags = state->mold.sys.atom.flags[i];
        if (flags & MD_FLAG_AMINO_ACID) amino_acid_present = true;
        if (flags & MD_FLAG_NUCLEOTIDE) nucleic_present = true;
        if (flags & MD_FLAG_ION) ion_present = true;
        if (flags & MD_FLAG_WATER) water_present = true;
        if (flags & MD_FLAG_COARSE_GRAINED) coarse_grained = true;

        if (!(flags & (MD_FLAG_AMINO_ACID | MD_FLAG_NUCLEOTIDE | MD_FLAG_ION | MD_FLAG_WATER))) {
            ligand_present = true;
        }
    }

    if (coarse_grained) {
        Representation* rep = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    if (amino_acid_present) {
        RepresentationType type = RepresentationType::Cartoon;
        ColorMapping color = ColorMapping::SecondaryStructure;

        if (state->mold.sys.instance.count > 1) {
            color = ColorMapping::InstId;
        } else {
            size_t res_count = md_instance_comp_count(&state->mold.sys.instance, 0);
            if (res_count < 20) {
                type = RepresentationType::BallAndStick;
                color = ColorMapping::Type;
            }
        }

        Representation* prot = create_representation(state, type, color, STR_LIT("protein"));
        snprintf(prot->name, sizeof(prot->name), "protein");
    }
    if (nucleic_present) {
        Representation* nucl = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Type, STR_LIT("nucleic"));
        snprintf(nucl->name, sizeof(nucl->name), "nucleic");
    }
    if (ion_present) {
        Representation* ion = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("ion"));
        snprintf(ion->name, sizeof(ion->name), "ion");
    }
    if (ligand_present) {
        Representation* ligand = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Type, STR_LIT("not (protein or nucleic or water or ion)"));
        snprintf(ligand->name, sizeof(ligand->name), "ligand");
    }
    if (water_present) {
        Representation* water = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("water"));
        water->scale.x = 0.5f;
        snprintf(water->name, sizeof(water->name), "water");
        water->enabled = false;
        if (!amino_acid_present && !nucleic_present && !ligand_present) {
            water->enabled = true;
        }
    }

done:
    if (orbitals_present) {
        Representation* rep = create_representation(state, RepresentationType::ElectronicStructure);
        snprintf(rep->name, sizeof(rep->name), "electronic structure");
        rep->enabled = true;

		// ONE representation, for the ground state moment. A group each meant a VeloxChem file with
		// transition dipoles opened with a representation per group, and with the electric,
		// magnetic and velocity sets that is three nobody asked for on top of the one they wanted.
		// The ground state is what a dipole means to someone who has not said otherwise; the rest
		// are a representation away, and the index slider covers the excited states within each.
		DipoleMoment dipoles[64];
		size_t num_dipoles = MIN(dipole_moments_gather(dipoles, ARRAY_SIZE(dipoles), state->mold.sys), ARRAY_SIZE(dipoles));
        for (size_t i = 0; i < num_dipoles; ++i) {
            // The group name is the identity here, the same string the path spells.
            if (dipoles[i].index != 0 || !str_eq(dipoles[i].label, STR_LIT("ground_state"))) continue;
            if (vec3_length(dipoles[i].vec) <= 1e-3f) continue;

            Representation* dipole_rep = create_representation(state, RepresentationType::DipoleMoment);
            dipole_rep->dipole.dipole_key   = dipoles[i].key;
            dipole_rep->dipole.dipole_index = dipoles[i].index;

            // Lower case, like every other auto created representation - "protein", "water",
            // "electronic structure". dipole_label_pretty title cases for menus and tooltips,
            // which is a different job and stays as it is.
            snprintf(dipole_rep->name, sizeof(dipole_rep->name), "ground state");
            dipole_rep->enabled = true;
            break;
        }
    }

    recompute_atom_visibility_mask(state);
}

void interpolate_system_state(ApplicationState* app) {
    ASSERT(app);
	const auto& sys  = app->mold.sys;

	size_t num_atoms = app->mold.state.num_atoms;
    if (num_atoms == 0 || !md_trajectory_num_frames(sys.trajectory)) return;

    const int64_t last_frame = MAX(0LL, (int64_t)md_trajectory_num_frames(sys.trajectory) - 1);
    // This is not actually time, but the fractional frame representation
    const double time = CLAMP(app->animation.frame, 0.0, double(last_frame));

    // Scaling factor for cubic spline
    const int64_t frame = (int64_t)time;
    const int64_t nearest_frame = CLAMP((int64_t)(time + 0.5), 0LL, last_frame);

    static int64_t curr_nearest_frame = -1;
    if (app->animation.interpolation == InterpolationMode::Nearest) {
        if (curr_nearest_frame == nearest_frame) {
            return;
        }
        app->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
    }
    curr_nearest_frame = nearest_frame;

    // This represents the frames that we would like to load into memory for interpolation (worst case).
    const int64_t frames[4] = {
        MAX(0LL, frame - 1),
        MAX(0LL, frame),
        MIN(frame + 1, last_frame),
        MIN(frame + 2, last_frame),
    };

    const size_t num_threads = task_system::pool_num_threads();

    // The number of atoms to be processed per thread when divided into chunks
    const uint32_t grain_size = 1024;

    md_allocator_i* temp_arena = app->allocator.frame;
    md_temp_scope_t temp = md_temp_begin_in(temp_arena);
    defer { md_temp_end(temp); };

    struct Payload {
        ApplicationState* app;
        float s;
        float t;
        InterpolationMode mode;

        int64_t nearest_frame;
        int64_t frames[4];

        md_system_state_t* src_states[4];
		md_system_state_t* dst_state;

        vec3_t* aabb_min;
        vec3_t* aabb_max;

        mat4_t recenter_transform;
    };

    const InterpolationMode mode = (frames[1] == frames[2]) ? InterpolationMode::Nearest : app->animation.interpolation;

    Payload payload = {
        .app = app,
        .s = 1.0f - CLAMP(app->animation.tension, 0.0f, 1.0f),
        .t = (float)fract(time),
        .mode = mode,
        .nearest_frame = nearest_frame,
        .frames = { frames[0], frames[1], frames[2], frames[3]},
        .dst_state = &app->mold.state,
        .aabb_min = md_temp_alloc_array(temp, vec3_t, num_threads),
        .aabb_max = md_temp_alloc_array(temp, vec3_t, num_threads),
    };

    // Stamp the destination with the frame it is about to represent. The interpolated state is not
    // loaded through md_trajectory_reader_load_frame, so nothing else would fill this in, and a
    // stale value is worse than an absent one. Nearest snaps to a whole frame; the other modes land
    // between two, which is exactly what the fractional part is for.
    app->mold.state.frame = (mode == InterpolationMode::Nearest) ? (double)nearest_frame : time;

    int requested_frames[4] = { 0 };
    int num_requested_frames = 0;

    switch (mode) {
        case InterpolationMode::Nearest:
            requested_frames[num_requested_frames++] = (int)nearest_frame;
            break;
        case InterpolationMode::Linear:
            requested_frames[num_requested_frames++] = (int)frames[1];
            requested_frames[num_requested_frames++] = (int)frames[2];
            break;
        case InterpolationMode::CubicSpline:
            requested_frames[num_requested_frames++] = (int)frames[0];
            requested_frames[num_requested_frames++] = (int)frames[1];
            requested_frames[num_requested_frames++] = (int)frames[2];
            requested_frames[num_requested_frames++] = (int)frames[3];
            break;
        default:
            ASSERT(false);
            break;
    }

    // Represents the frame cache slot indices for the requested frames, -1 if not present in cache
    int frame_cache_slot_idx[4] = { -1, -1, -1, -1 };

    // Array of frame_cache slots which requires a load
    int frame_cache_load_slot[4] = {0};
    int num_frames_to_load = 0;

    for (int i = 0; i < num_requested_frames; ++i) {
        int slot_idx = -1;
        if (!find_frame_in_cache(&slot_idx, requested_frames[i], &app->mold.frame_cache)) {
            slot_idx = find_lru_cache_slot(&app->mold.frame_cache);
            frame_cache_load_slot[num_frames_to_load++] = slot_idx;
        }
        set_mru_cache_slot(&app->mold.frame_cache, slot_idx);
        app->mold.frame_cache.frame_idx[slot_idx] = requested_frames[i];
        frame_cache_slot_idx[i] = slot_idx;
    }

    for (int i = 0; i < num_requested_frames; ++i) {
        int slot_idx = frame_cache_slot_idx[i];
        payload.src_states[i] = &app->mold.frame_cache.states[slot_idx];
    }

    // This holds the chain of tasks we are about to submit
    task_system::ID tasks[16] = {0};
    int num_tasks = 0;

    task_system::ID load_task = task_system::create_pool_task(STR_LIT("## Load Frames"), num_frames_to_load,
        [data = &payload, frame_cache_load_slot](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            for (uint32_t i = range_beg; i < range_end; ++i) {
                int slot_idx = frame_cache_load_slot[i];
                int frame_idx = data->app->mold.frame_cache.frame_idx[slot_idx];
                md_system_state_t* state = &data->app->mold.frame_cache.states[slot_idx];
                md_trajectory_load_frame(data->app->mold.sys.trajectory, frame_idx, state);
            }
        }
    );

    tasks[num_tasks++] = load_task;

    switch (mode) {
        case InterpolationMode::Nearest: {
            task_system::ID interp_task = task_system::create_pool_task(STR_LIT("## Interpolate"), [data = &payload]() {
                data->dst_state->unitcell = data->src_states[0]->unitcell;
                MEMCPY(data->dst_state->x, data->src_states[0]->x, sizeof(float) * data->dst_state->num_atoms);
                MEMCPY(data->dst_state->y, data->src_states[0]->y, sizeof(float) * data->dst_state->num_atoms);
                MEMCPY(data->dst_state->z, data->src_states[0]->z, sizeof(float) * data->dst_state->num_atoms);
            });
            tasks[num_tasks++] = interp_task;
            break;
        }
        case InterpolationMode::Linear: {
            task_system::ID iterp_cell_task = task_system::create_pool_task(STR_LIT("## Interp Unitcell"), [data = &payload]() {
                double x  = lerp(data->src_states[0]->unitcell.x,  data->src_states[1]->unitcell.x,  data->t);
                double y  = lerp(data->src_states[0]->unitcell.y,  data->src_states[1]->unitcell.y,  data->t);
                double z  = lerp(data->src_states[0]->unitcell.z,  data->src_states[1]->unitcell.z,  data->t);
                double xy = lerp(data->src_states[0]->unitcell.xy, data->src_states[1]->unitcell.xy, data->t);
                double xz = lerp(data->src_states[0]->unitcell.xz, data->src_states[1]->unitcell.xz, data->t);
                double yz = lerp(data->src_states[0]->unitcell.yz, data->src_states[1]->unitcell.yz, data->t);
                data->dst_state->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
			});

            task_system::ID interp_coord_task = task_system::create_pool_task(STR_LIT("## Interp Coord Data"), (uint32_t)num_atoms, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                float* dst_x = data->dst_state->x + range_beg;
                float* dst_y = data->dst_state->y + range_beg;
                float* dst_z = data->dst_state->z + range_beg;
                const float* src_x[2] = { data->src_states[0]->x + range_beg, data->src_states[1]->x + range_beg};
                const float* src_y[2] = { data->src_states[0]->y + range_beg, data->src_states[1]->y + range_beg};
                const float* src_z[2] = { data->src_states[0]->z + range_beg, data->src_states[1]->z + range_beg};

                md_util_interpolate_linear(dst_x, dst_y, dst_z, src_x, src_y, src_z, count, &data->dst_state->unitcell, data->t);
            }, grain_size);

			tasks[num_tasks++] = iterp_cell_task;
            tasks[num_tasks++] = interp_coord_task;

            break;
        }
        case InterpolationMode::CubicSpline: {
            task_system::ID iterp_cell_task = task_system::create_pool_task(STR_LIT("## Interp Unitcell"), [data = &payload]() {
                double x  = cubic_spline(data->src_states[0]->unitcell.x,  data->src_states[1]->unitcell.x,  data->src_states[2]->unitcell.x,  data->src_states[3]->unitcell.x,  data->t, data->s);
                double y  = cubic_spline(data->src_states[0]->unitcell.y,  data->src_states[1]->unitcell.y,  data->src_states[2]->unitcell.y,  data->src_states[3]->unitcell.y,  data->t, data->s);
                double z  = cubic_spline(data->src_states[0]->unitcell.z,  data->src_states[1]->unitcell.z,  data->src_states[2]->unitcell.z,  data->src_states[3]->unitcell.z,  data->t, data->s);
                double xy = cubic_spline(data->src_states[0]->unitcell.xy, data->src_states[1]->unitcell.xy, data->src_states[2]->unitcell.xy, data->src_states[3]->unitcell.xy, data->t, data->s);
                double xz = cubic_spline(data->src_states[0]->unitcell.xz, data->src_states[1]->unitcell.xz, data->src_states[2]->unitcell.xz, data->src_states[3]->unitcell.xz, data->t, data->s);
                double yz = cubic_spline(data->src_states[0]->unitcell.yz, data->src_states[1]->unitcell.yz, data->src_states[2]->unitcell.yz, data->src_states[3]->unitcell.yz, data->t, data->s);
                data->dst_state->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
            });

            task_system::ID interp_coord_task = task_system::create_pool_task(STR_LIT("## Interp Coord Data"), (uint32_t)num_atoms, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                float* dst_x = data->dst_state->x + range_beg;
                float* dst_y = data->dst_state->y + range_beg;
                float* dst_z = data->dst_state->z + range_beg;
                const float* src_x[4] = { data->src_states[0]->x + range_beg, data->src_states[1]->x + range_beg, data->src_states[2]->x + range_beg, data->src_states[3]->x + range_beg};
                const float* src_y[4] = { data->src_states[0]->y + range_beg, data->src_states[1]->y + range_beg, data->src_states[2]->y + range_beg, data->src_states[3]->y + range_beg};
                const float* src_z[4] = { data->src_states[0]->z + range_beg, data->src_states[1]->z + range_beg, data->src_states[2]->z + range_beg, data->src_states[3]->z + range_beg};

                md_util_interpolate_cubic_spline(dst_x, dst_y, dst_z, src_x, src_y, src_z, count, &data->dst_state->unitcell, data->t, data->s);
            }, grain_size);

            tasks[num_tasks++] = iterp_cell_task;
            tasks[num_tasks++] = interp_coord_task;
            
            break;
        }
        default:
            ASSERT(false);
            break;
    }

    {
        // Calculate a global AABB for the molecule
        task_system::ID aabb_task = task_system::create_pool_task(STR_LIT("## Compute AABB"), (uint32_t)num_atoms, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            size_t range_len = range_end - range_beg;
            const float* x = data->dst_state->x + range_beg;
            const float* y = data->dst_state->y + range_beg;
            const float* z = data->dst_state->z + range_beg;

            md_temp_scope_t temp = md_temp_begin();
            defer { md_temp_end(temp); };
            float* r = md_temp_alloc_array(temp, float, range_len);
            md_atom_extract_radii(r, range_beg, range_len, &data->app->mold.sys.atom);

            vec3_t aabb_min = vec3_set1(FLT_MAX);
            vec3_t aabb_max = vec3_set1(-FLT_MAX);
            md_util_aabb_compute(aabb_min.elem, aabb_max.elem, x, y, z, r, 0, range_len);

            data->aabb_min[thread_num] = aabb_min;
            data->aabb_max[thread_num] = aabb_max;
        });
        tasks[num_tasks++] = aabb_task;
    }

    if (sys.protein_backbone.segment.count > 0 && sys.protein_backbone.segment.angle) {
        switch (mode) {
            case InterpolationMode::Nearest: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Compute Backbone Angles"), [data = &payload]() {
                    const md_backbone_angles_t* src_angles[2] = {
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    const md_backbone_angles_t* src_angle = data->t < 0.5f ? src_angles[0] : src_angles[1];
                    MEMCPY(data->app->mold.sys.protein_backbone.segment.angle, src_angle, data->app->mold.sys.protein_backbone.segment.count * sizeof(md_backbone_angles_t));
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            case InterpolationMode::Linear: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Compute Backbone Angles"), (uint32_t)sys.protein_backbone.segment.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_backbone_angles_t* src_angles[2] = {
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    md_system_t& sys = data->app->mold.sys;
                    for (size_t i = range_beg; i < range_end; ++i) {
                        float phi[2] = {src_angles[0][i].phi, src_angles[1][i].phi};
                        float psi[2] = {src_angles[0][i].psi, src_angles[1][i].psi};

                        phi[1] = deperiodize_orthof(phi[1], phi[0], (float)TWO_PI);
                        psi[1] = deperiodize_orthof(psi[1], psi[0], (float)TWO_PI);

                        float final_phi = lerp(phi[0], phi[1], data->t);
                        float final_psi = lerp(psi[0], psi[1], data->t);
                        sys.protein_backbone.segment.angle[i] = {deperiodize_orthof(final_phi, 0, (float)TWO_PI), deperiodize_orthof(final_psi, 0, (float)TWO_PI)};
                    }
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            case InterpolationMode::CubicSpline: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Interpolate Backbone Angles"), (uint32_t)sys.protein_backbone.segment.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_backbone_angles_t* src_angles[4] = {
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[0],
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[2],
                        data->app->trajectory_data.backbone_angles.data + data->app->trajectory_data.backbone_angles.stride * data->frames[3],
                    };
                    md_system_t& sys = data->app->mold.sys;
                    for (size_t i = range_beg; i < range_end; ++i) {
                        float phi[4] = {src_angles[0][i].phi, src_angles[1][i].phi, src_angles[2][i].phi, src_angles[3][i].phi};
                        float psi[4] = {src_angles[0][i].psi, src_angles[1][i].psi, src_angles[2][i].psi, src_angles[3][i].psi};

                        phi[0] = deperiodize_orthof(phi[0], phi[1], (float)TWO_PI);
                        phi[2] = deperiodize_orthof(phi[2], phi[1], (float)TWO_PI);
                        phi[3] = deperiodize_orthof(phi[3], phi[2], (float)TWO_PI);

                        psi[0] = deperiodize_orthof(psi[0], psi[1], (float)TWO_PI);
                        psi[2] = deperiodize_orthof(psi[2], psi[1], (float)TWO_PI);
                        psi[3] = deperiodize_orthof(psi[3], psi[2], (float)TWO_PI);

                        float final_phi = cubic_spline(phi[0], phi[1], phi[2], phi[3], data->t, data->s);
                        float final_psi = cubic_spline(psi[0], psi[1], psi[2], psi[3], data->t, data->s);
                        sys.protein_backbone.segment.angle[i] = {deperiodize_orthof(final_phi, 0, (float)TWO_PI), deperiodize_orthof(final_psi, 0, (float)TWO_PI)};
                    }
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            default:
                ASSERT(false);
                break;
        }
    }

    if (sys.protein_backbone.segment.count > 0 && sys.protein_backbone.segment.secondary_structure) {
        if (md_array_size(app->mold.interpolated_properties.secondary_structure) != sys.protein_backbone.segment.count) {
			MD_LOG_ERROR("Secondary structure array size does not match the number of segments.");
        }
        size_t num_backbone_segments = sys.protein_backbone.segment.count;
        task_system::ID ss_task = task_system::create_pool_task(STR_LIT("## Interpolate Secondary Structures"), (uint32_t)num_backbone_segments, [data = &payload, mode](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            const md_secondary_structure_t* ss_data = data->app->trajectory_data.secondary_structure_render.data ?
                data->app->trajectory_data.secondary_structure_render.data :
                data->app->trajectory_data.secondary_structure.data;
            const size_t ss_stride = data->app->trajectory_data.secondary_structure_render.data ?
                data->app->trajectory_data.secondary_structure_render.stride :
                data->app->trajectory_data.secondary_structure.stride;
            const md_secondary_structure_t* src_ss[4] = {
                ss_data + ss_stride * data->frames[0],
                ss_data + ss_stride * data->frames[1],
                ss_data + ss_stride * data->frames[2],
                ss_data + ss_stride * data->frames[3],
            };
            const md_secondary_structure_t* src_ss_nearest = data->t < 0.5f ? src_ss[1] : src_ss[2];
            auto normalize_ss = [](md_gl_secondary_structure_t ss) {
                ss.helix = CLAMP(ss.helix, 0.0f, 1.0f);
                ss.sheet = CLAMP(ss.sheet, 0.0f, 1.0f);
                float sum = ss.helix + ss.sheet;
                if (sum > 1.0f) {
                    ss.helix /= sum;
                    ss.sheet /= sum;
                }
                return ss;
            };
            auto blend_ss = [normalize_ss](md_gl_secondary_structure_t a, md_gl_secondary_structure_t b, float t) {
                md_gl_secondary_structure_t ss = {
                    .helix = lerp(a.helix, b.helix, t),
                    .sheet = lerp(a.sheet, b.sheet, t),
                };
                return normalize_ss(ss);
            };
            auto smoothstep5 = [](float t) {
                t = CLAMP(t, 0.0f, 1.0f);
                return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f);
            };
            auto is_eq = [](md_gl_secondary_structure_t a, md_gl_secondary_structure_t b) {
                return a.helix == b.helix && a.sheet == b.sheet;
            };
            switch (mode) {
            default:
                MD_LOG_DEBUG("Unsupported interpolation mode for secondary structure interpolation");
                [[fallthrough]];
            case InterpolationMode::Nearest: {
                for (size_t i = range_beg; i < range_end; ++i) {
                    md_secondary_structure_t ss = src_ss_nearest[i];
                    // Set both the analytical (nearest) and interpolated secondary structure (rendering)
                    data->app->mold.sys.protein_backbone.segment.secondary_structure[i] = ss;
                    data->app->mold.interpolated_properties.secondary_structure[i] = md_gl_secondary_structure_convert(ss);
                }
                break;
            }
            case InterpolationMode::Linear: {
                for (size_t i = range_beg; i < range_end; ++i) {
                    md_secondary_structure_t ss[2] = { src_ss[1][i], src_ss[2][i] };
                    md_gl_secondary_structure_t ss_gl[2] = { md_gl_secondary_structure_convert(ss[0]), md_gl_secondary_structure_convert(ss[1]) };
                    md_gl_secondary_structure_t ss_gl_i = blend_ss(ss_gl[0], ss_gl[1], data->t);
                    // Set both the analytical (nearest) and interpolated secondary structure (rendering)
                    data->app->mold.sys.protein_backbone.segment.secondary_structure[i] = src_ss_nearest[i];
                    data->app->mold.interpolated_properties.secondary_structure[i] = ss_gl_i;
                }
                break;
            }
            case InterpolationMode::CubicSpline: {
                
                for (size_t i = range_beg; i < range_end; ++i) {
                    md_secondary_structure_t ss[4] = { src_ss[0][i], src_ss[1][i], src_ss[2][i], src_ss[3][i] };
                    
                    md_gl_secondary_structure_t ss_gl[4] = {
                        md_gl_secondary_structure_convert(ss[0]),
                        md_gl_secondary_structure_convert(ss[1]),
                        md_gl_secondary_structure_convert(ss[2]),
                        md_gl_secondary_structure_convert(ss[3]),
                    };

                    // Cleanup isolated temporal assignments to reduce noise during transitions.
                    if (is_eq(ss_gl[0], ss_gl[2]) && !is_eq(ss_gl[1], ss_gl[0])) {
                        ss_gl[1] = ss_gl[0];
                    }
                    if (is_eq(ss_gl[1], ss_gl[3]) && !is_eq(ss_gl[2], ss_gl[1])) {
                        ss_gl[2] = ss_gl[1];
                    }

                    md_gl_secondary_structure_t ss_gl_i = blend_ss(ss_gl[1], ss_gl[2], smoothstep5(data->t));
                    // Set both the analytical (nearest) and interpolated secondary structure (rendering)
                    data->app->mold.sys.protein_backbone.segment.secondary_structure[i] = src_ss_nearest[i];
                    data->app->mold.interpolated_properties.secondary_structure[i] = ss_gl_i;
                }
                break;
            }
            }
        });
        tasks[num_tasks++] = ss_task;

#if 1
        // Task for cleaning up isolated coils to its neighbors (if the same), to reduce the noise in the secondary structure during transitions. This is a non temporal filtering step
        task_system::ID ss_cleanup_task = task_system::create_pool_task(STR_LIT("## Cleanup Secondary Structures"), [data = &payload]() {
            // Cleanup isolated coils to reduce noise during transitions
            auto is_eq = [](md_gl_secondary_structure_t a, md_gl_secondary_structure_t b) {
                return a.helix == b.helix && a.sheet == b.sheet;
            };
            const md_gl_secondary_structure_t ss_coil = { 0,0 };
            const md_gl_secondary_structure_t ss_helix = { .helix = 1.0f };
            const md_gl_secondary_structure_t ss_sheet = { .sheet = 1.0f };

            md_gl_secondary_structure_t* ss_gl = data->app->mold.interpolated_properties.secondary_structure;
            for (size_t i = 0; i < data->app->mold.sys.protein_backbone.range.count; ++i) {
                size_t range_beg = data->app->mold.sys.protein_backbone.range.offset[i];
                size_t range_end = data->app->mold.sys.protein_backbone.range.offset[i + 1];
                for (size_t j = range_beg + 1; j + 1 < range_end; ++j) {
                    // Set isolated coils between matching structured segments to reduce noise during transitions.
                    if (is_eq(ss_gl[j - 1], ss_helix) && is_eq(ss_gl[j + 1], ss_helix) && is_eq(ss_gl[j], ss_coil)) {
                        ss_gl[j] = ss_helix;
                    }
                    if (is_eq(ss_gl[j - 1], ss_sheet) && is_eq(ss_gl[j + 1], ss_sheet) && is_eq(ss_gl[j], ss_coil)) {
                        ss_gl[j] = ss_sheet;
                    }
                }
            }
        });
        tasks[num_tasks++] = ss_cleanup_task;
#endif

        app->mold.dirty_gpu_buffers |= MolBit_DirtySecondaryStructure;
    }

    if (num_tasks > 0) {
        for (int i = 1; i < num_tasks; ++i) {
            task_system::set_task_dependency(tasks[i], tasks[i-1]);
        }
        task_system::enqueue_task(tasks[0]);
        task_system::task_wait_for(tasks[num_tasks - 1]);
    }

    vec3_t aabb_min = payload.aabb_min[0];
    vec3_t aabb_max = payload.aabb_max[0];
    for (size_t i = 1; i < task_system::pool_num_threads(); ++i) {
        aabb_min = vec3_min(aabb_min, payload.aabb_min[i]);
        aabb_max = vec3_max(aabb_max, payload.aabb_max[i]);
    }
    app->mold.sys_aabb_min = aabb_min;
    app->mold.sys_aabb_max = aabb_max;

    // unitcell transform is essentially just a translation to place the center of the unitcell at the origin
    mat3_t A;
    md_unitcell_A_extract_float(A.elem, &app->mold.state.unitcell);
    vec3_t c = mat3_mul_vec3(A, vec3_set(0.5f, 0.5f, 0.5f));
    app->mold.unitcell_transform = mat4_translate(-c.x, -c.y, -c.z);

#if 0
    if (sys.unitcell.flags) {
        vec3_t c = sys.unitcell.basis * vec3_set1(0.5f);
        app->mold.model_mat = mat4_translate_vec3(-c);
    }
#endif

    app->mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
}

void recenter_mark_query_dirty(ApplicationState* state) {
    ASSERT(state);
    state->operations.recenter_query.version += 1;
    if (state->operations.recenter_query.version == 0) {
        state->operations.recenter_query.version = 1;
    }
}

void recenter_mark_selection_dirty(ApplicationState* state) {
    ASSERT(state);
    state->operations.selection_version += 1;
    if (state->operations.selection_version == 0) {
        state->operations.selection_version = 1;
    }
}

const md_bitfield_t& recenter_get_active_target_mask(const ApplicationState* state) {
    ASSERT(state);
    return state->operations.recenter_query.enabled ? state->operations.recenter_query.mask : state->operations.selection_mask;
}

uint64_t recenter_get_active_target_version(const ApplicationState* state) {
    ASSERT(state);
    const uint64_t source_version = state->operations.recenter_query.enabled ? state->operations.recenter_query.evaluated_version : state->operations.selection_version;
    const uint64_t source_idx = state->operations.recenter_query.enabled ? 1 : 0;
    return (source_version << 1) | source_idx;
}

bool recenter_update_query_mask(ApplicationState* state) {
    ASSERT(state);

    auto& query = state->operations.recenter_query;
    const uint64_t ir_fingerprint = state->script.ir ? state->script.ir_fingerprint : 0;
    if (query.ir_fingerprint != ir_fingerprint) {
        query.ir_fingerprint = ir_fingerprint;
        recenter_mark_query_dirty(state);
    }

    if (!query.enabled) {
        return false;
    }

    if (query.dynamic) {
        recenter_mark_query_dirty(state);
    }

    if (query.evaluated_version == query.version) {
        return false;
    }

    md_bitfield_clear(&query.mask);
    query.dynamic = false;
    query.valid = md_filter(&query.mask, str_from_cstr(query.query), &state->mold.sys, &state->mold.state, state->script.ir, &query.dynamic, query.error, sizeof(query.error));
    query.evaluated_version = query.version;
    return true;
}

void recenter_update(ApplicationState* state) {
    ASSERT(state);
    recenter_update_query_mask(state);
    recenter_update_target_data(state);
}

void recenter_update_target_data(ApplicationState* state) {
    if (!state->mold.sys.trajectory) return;

    const md_bitfield_t& target_mask = recenter_get_active_target_mask(state);
    const uint64_t target_version = recenter_get_active_target_version(state);
    if (target_version != state->operations.initial_frame.target_version) {
        // Need to recalculate the initial frame
        state->operations.initial_frame.target_version = target_version;
        size_t count = md_bitfield_popcount(&target_mask);

        md_array_resize(state->operations.initial_frame.rel_xyzw, count, state->allocator.persistent);

        state->operations.initial_frame.com = vec3_zero();

        if (state->operations.initial_frame.rel_xyzw && count > 0) {
            // Fetch initial frame data required for orienting the structure
            size_t num_atoms = state->mold.sys.atom.count;
            float* temp_x = (float*)md_vm_arena_push(state->allocator.frame, sizeof(float) * num_atoms);
            float* temp_y = (float*)md_vm_arena_push(state->allocator.frame, sizeof(float) * num_atoms);
            float* temp_z = (float*)md_vm_arena_push(state->allocator.frame, sizeof(float) * num_atoms);

            md_system_state_t temp_state = { 0 };
            temp_state.num_atoms = num_atoms;
            temp_state.x = temp_x;
            temp_state.y = temp_y;
            temp_state.z = temp_z;
            md_trajectory_load_frame(state->mold.sys.trajectory, 0, &temp_state);

            md_bitfield_iter_t it = md_bitfield_iter_create(&target_mask);
            int dst_idx = 0;
            while (md_bitfield_iter_next(&it)) {
                uint64_t src_idx = md_bitfield_iter_idx(&it);
                float mass = md_atom_mass(&state->mold.sys.atom, src_idx);
                state->operations.initial_frame.rel_xyzw[dst_idx++] = vec4_set(temp_x[src_idx], temp_y[src_idx], temp_z[src_idx], mass);
            }

			state->operations.initial_frame.com = md_util_com_compute_vec4(state->operations.initial_frame.rel_xyzw, NULL, count, &temp_state.unitcell);
			md_util_convert_to_relative_coordinates_vec4(state->operations.initial_frame.rel_xyzw, state->operations.initial_frame.com, count, &temp_state.unitcell);
        }
    }
}

void recenter_calculate_transform(float M[4][4], const ApplicationState* app) {
    ASSERT(M);
    ASSERT(app);

    const md_bitfield_t& target_mask = recenter_get_active_target_mask(app);
    size_t count = md_bitfield_popcount(&target_mask);
    mat4_t transform = mat4_ident();

    if (count > 0) {
        md_temp_scope_t temp = md_temp_begin_in(app->allocator.frame);
        defer { md_temp_end(temp); };

        // Extract xyzw subset of target
        vec4_t* target_xyzw = md_temp_alloc_array(temp, vec4_t, count);

		md_util_system_extract_xyzw_from_mask(target_xyzw, &target_mask, &app->mold.sys, &app->mold.state);
		vec3_t target_com = md_util_com_compute_vec4(target_xyzw, NULL, count, &app->mold.state.unitcell);

        // Calculate target
        vec3_t target = {0};
        if (md_unitcell_flags(&app->mold.state.unitcell) != 0) {
            mat3_t A = {0};
            md_unitcell_A_extract_float(A.elem, &app->mold.state.unitcell);
            target = mat3_mul_vec3(A, vec3_set1(0.5f));
        } 

        // Place the target in mutually consistent images and take its centre. When there is a
        // reference to hold the orientation against, the rotation, the centre and each point's
        // periodic image are solved for together, so the images are chosen to minimise the alignment
        // residual rather than being committed to beforehand by a criterion unrelated to it.
        //
        // Both calls report the centre in the image the STATE coordinates actually occupy, which is
        // what makes the translation below valid - it is applied to those same coordinates, untouched.
        // A centre folded into the reference cell would place a target living outside that cell one
        // lattice vector off, and with a rotation in play, off by R times a lattice vector.
        mat3_t R = mat3_ident();

        // The reference has to have been built from the SAME target that is being fitted now.
        // A size match is not sufficient: the selection can change to a different set of equal
        // size between recenter_update() and this call, which would silently pair up unrelated
        // atoms and yield a garbage rotation. The version is the identity of the target.
        const bool reference_valid =
            app->operations.initial_frame.rel_xyzw &&
            md_array_size(app->operations.initial_frame.rel_xyzw) == count &&
            app->operations.initial_frame.target_version == recenter_get_active_target_version(app);

        if (app->operations.fixate_orientation && reference_valid) {
		    md_util_convert_to_relative_coordinates_vec4(target_xyzw, target_com, count, &app->mold.state.unitcell);
            R = md_util_optimal_rotation_rel_vec4(target_xyzw, app->operations.initial_frame.rel_xyzw, count);
            R = mat3_orthonormalize(R);
        }

        const mat4_t A = app->operations.alignment_mat;
        transform = mat4_translate_vec3(target) * A * mat4_from_mat3(R) * mat4_translate_vec3(-target_com);
    }
    mat4_store((float*)M, transform);
}

bool picking_range_reserve(PickingRange* out_range, PickingSpace* space, PickingDomainID domain, size_t count) {
    ASSERT(space);

    if (count > 0 && space->num_ranges < ARRAY_SIZE(space->ranges)) {
        PickingRange* curr_range = &space->ranges[space->num_ranges++];
        PickingRange* prev_range = space->num_ranges > 1 ? &space->ranges[space->num_ranges - 2] : NULL;
        curr_range->domain = domain;
        curr_range->beg = prev_range ? prev_range->end : 0;
        curr_range->end = curr_range->beg + (uint32_t)count;
        if (out_range) {
            MEMCPY(out_range, curr_range, sizeof(PickingRange));
        }
        return true;
    }

    return false;
}

void picking_handler_new_frame(PickingHandler* handler) {
    ASSERT(handler);
    handler->frame_idx += 1;

    const uint32_t slot_idx = handler->frame_idx % ARRAY_SIZE(handler->history);
    handler->history[slot_idx].submitted_frame_idx = handler->frame_idx;
    handler->history[slot_idx].space = PickingSpace{};
}

PickingSpace* picking_handler_current_space(PickingHandler* handler) {
    ASSERT(handler);
    return &handler->history[handler->frame_idx % ARRAY_SIZE(handler->history)].space;
}

const PickingSpace* picking_handler_find_space(const PickingHandler& handler, uint32_t submitted_frame_idx) {
    for (size_t i = 0; i < ARRAY_SIZE(handler.history); ++i) {
        const auto& hist = handler.history[i];
        if (hist.submitted_frame_idx == submitted_frame_idx) {
            return &hist.space;
        }
    }

    return nullptr;
}

void picking_surface_init(PickingSurface* surface, PickingSourceID source) {
    ASSERT(surface);

    *surface = PickingSurface{};
    surface->source = source;

    for (size_t i = 0; i < ARRAY_SIZE(surface->slots); ++i) {
        auto& slot = surface->slots[i];
        glGenBuffers(1, &slot.color_pbo);
        glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.color_pbo);
        glBufferData(GL_PIXEL_PACK_BUFFER, 4, nullptr, GL_DYNAMIC_READ);

        glGenBuffers(1, &slot.depth_pbo);
        glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.depth_pbo);
        glBufferData(GL_PIXEL_PACK_BUFFER, 4, nullptr, GL_DYNAMIC_READ);
    }

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
}

void picking_surface_free(PickingSurface* surface) {
    ASSERT(surface);

    for (size_t i = 0; i < ARRAY_SIZE(surface->slots); ++i) {
        auto& slot = surface->slots[i];
        if (slot.color_pbo) glDeleteBuffers(1, &slot.color_pbo);
        if (slot.depth_pbo) glDeleteBuffers(1, &slot.depth_pbo);
    }

    *surface = PickingSurface{};
}

bool picking_surface_submit_readback(
    PickingSurface* surface,
    uint32_t fbo,
    uint32_t width,
    uint32_t height,
    uint32_t submitted_frame_idx,
    vec2_t surface_coord,
    vec2_t screen_coord,
    const mat4_t& clip_to_world
) {
    ASSERT(surface);

    if (!fbo || width == 0 || height == 0) {
        return false;
    }

    const int x = (int)surface_coord.x;
    const int y = (int)surface_coord.y;
    if (x < 0 || y < 0 || x >= (int)width || y >= (int)height) {
        return false;
    }

    const uint32_t queue_idx = surface->slot_cursor % ARRAY_SIZE(surface->slots);
    surface->slot_cursor += 1;

    auto& slot = surface->slots[queue_idx];
    if (slot.color_pbo == 0 || slot.depth_pbo == 0) {
        MD_LOG_ERROR("Invalid PBOs in picking surface slot");
        return false;
    }

    slot.submitted_frame_idx = submitted_frame_idx;
    slot.pending = true;
    slot.viewport_width = width;
    slot.viewport_height = height;
    slot.surface_coord = surface_coord;
    slot.screen_coord = screen_coord;
    slot.clip_to_world = clip_to_world;

    PUSH_GPU_SECTION("QUEUE PICKING READBACK")
    glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo);
    glReadBuffer(GL_COLOR_ATTACHMENT_PICKING);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.color_pbo);
    glReadPixels(x, y, 1, 1, GL_BGRA, GL_UNSIGNED_BYTE, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.depth_pbo);
    glReadPixels(x, y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, 0);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    POP_GPU_SECTION()

    return true;
}

static const PickingRange* find_picking_range(const PickingSpace* space, uint32_t raw_idx) {
    ASSERT(space);

    for (size_t i = 0; i < space->num_ranges; ++i) {
        const PickingRange* range = &space->ranges[i];
        if (range->beg <= raw_idx && raw_idx < range->end) {
            return range;
        }
    }

    return nullptr;
}

bool picking_surface_poll_hit(
    PickingHit* out_hit,
    PickingSurface* surface,
    const PickingHandler& handler
) {
    ASSERT(out_hit);
    ASSERT(surface);

    *out_hit = {};

    if (surface->slot_cursor < ARRAY_SIZE(surface->slots)) {
        return false;
    }

    const uint32_t read_idx = surface->slot_cursor % ARRAY_SIZE(surface->slots);
    auto& slot = surface->slots[read_idx];
    if (!slot.pending) {
        return false;
    }

    uint8_t color[4] = {};
    float depth = 1.0f;

    PUSH_GPU_SECTION("POLL PICKING READBACK")
    glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.color_pbo);
    glGetBufferSubData(GL_PIXEL_PACK_BUFFER, 0, sizeof(color), color);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, slot.depth_pbo);
    glGetBufferSubData(GL_PIXEL_PACK_BUFFER, 0, sizeof(depth), &depth);

    glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    POP_GPU_SECTION()

    slot.pending = false;

    const uint32_t raw_idx = (color[0] << 16) | (color[1] << 8) | (color[2] << 0) | (color[3] << 24);

    const PickingSpace* space = picking_handler_find_space(handler, slot.submitted_frame_idx);
    if (!space) {
        return false;
    }

    PickingRange range = {0};
    if (const PickingRange* space_range = find_picking_range(space, raw_idx)) {
        range = *space_range;
    }

    out_hit->source = surface->source;
    out_hit->domain = range.domain;
    out_hit->frame_idx = slot.submitted_frame_idx;
    out_hit->raw_idx = raw_idx;
    out_hit->local_idx = raw_idx - range.beg;
    out_hit->surface_coord = slot.surface_coord;
    out_hit->screen_coord = slot.screen_coord;
    out_hit->depth = depth;

    const vec4_t viewport = {0, 0, (float)slot.viewport_width, (float)slot.viewport_height};
    out_hit->world_pos = mat4_unproject({slot.surface_coord.x, slot.surface_coord.y, depth}, slot.clip_to_world, viewport);

    return true;
}

bool picking_surface_submit_readback_and_poll_hit(
    PickingHit* out_hit,
    PickingSurface* surface,
    const PickingHandler& handler,
    const PickingReadbackRequest& request
) {
    ASSERT(out_hit);
    ASSERT(surface);

    picking_surface_submit_readback(
        surface,
        request.fbo,
        request.width,
        request.height,
        handler.frame_idx,
        request.surface_coord,
        request.screen_coord,
        request.clip_to_world
    );

    return picking_surface_poll_hit(out_hit, surface, handler);
}

InteractionSurfaceState interaction_surface(InteractionSurfaceID id, const vec2_t& size, InteractionSurfaceFlags flags) {
    InteractionSurfaceState state = {};

    ImGuiWindow* window = ImGui::GetCurrentWindow();
    if (!window) {
        MD_LOG_ERROR("No current ImGui window for interaction surface");
        return state;
    }

    state.surface_id = id;

    static const ImGuiButtonFlags btn_flags = ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_AllowOverlap;
    ImGui::InvisibleButton("interaction surface button", vec_cast(size), btn_flags);
    state.item_id = ImGui::GetItemID();

    state.hovered = ImGui::IsItemHovered();
    state.active = ImGui::IsItemActive();
    state.activated = ImGui::IsItemActivated();
    state.deactivated = ImGui::IsItemDeactivated();

    ImDrawList* draw_list = window->DrawList;
    ASSERT(draw_list);

    const ImVec2 canvas_min = ImGui::GetItemRectMin();
    const ImVec2 canvas_max = ImGui::GetItemRectMax();
    const ImVec2 canvas_size = ImGui::GetItemRectSize();

    state.surface_size = vec_cast(canvas_size);

    const ImVec2 mouse_pos = ImGui::GetMousePos();
    const ImVec2 local_mouse = mouse_pos - canvas_min;
    state.mouse_local = vec_cast(local_mouse);

    if (state.activated || state.active || state.deactivated) {
        if (ImGui::IsKeyPressed(ImGuiMod_Shift, false)) {
            ImGui::ResetMouseDragDelta(ImGuiMouseButton_Left);
            ImGui::ResetMouseDragDelta(ImGuiMouseButton_Right);
        }

        if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
            ImGuiMouseButton button = ImGuiMouseButton_Left;
            if (ImGui::IsMouseDown(ImGuiMouseButton_Left) || ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                state.selection_mode = InteractionSelectionMode::Append;
                button = ImGuiMouseButton_Left;
            } else if (ImGui::IsMouseDown(ImGuiMouseButton_Right) || ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
                state.selection_mode = InteractionSelectionMode::Remove;
                button = ImGuiMouseButton_Right;
            }

            if (state.selection_mode != InteractionSelectionMode::None) {
                const ImVec2 ext = ImGui::GetMouseDragDelta(button);
                const ImVec2 pos = ImGui::GetMousePos() - ext;

                ImVec2 sel_min = ImClamp(ImMin(pos, pos + ext), canvas_min, canvas_max);
                ImVec2 sel_max = ImClamp(ImMax(pos, pos + ext), canvas_min, canvas_max);

                draw_list->AddRectFilled(sel_min, sel_max, 0x22222222);
                draw_list->AddRect(sel_min, sel_max, 0x88888888);

                sel_min -= canvas_min;
                sel_max -= canvas_min;

                state.region_min = vec_cast(sel_min);
                state.region_max = vec_cast(sel_max);
            }
        }
    }

    if (flags & InteractionSurfaceFlags_NoRegionSelect) {
        state.region_min = { 0, 0 };
        state.region_max = { 0, 0 };
    }

    return state;
}

bool interaction_surface_hit_extract(PickingHit* out_hit, const InteractionSurfaceState& state, const InteractionSurfaceHitArgs& args) {
    ASSERT(out_hit);
    if (state.hovered && args.fbo && args.width && args.height) {
        const ImVec2 local_coord = ImVec2(state.mouse_local.x, state.surface_size.y - state.mouse_local.y) * ImGui::GetIO().DisplayFramebufferScale;

        PickingReadbackRequest request = {
            .fbo = args.fbo,
            .width = args.width,
            .height = args.height,
            .surface_coord = {local_coord.x, local_coord.y},
            .screen_coord = {state.mouse_local.x, state.mouse_local.y},
            .clip_to_world = args.clip_to_world,
        };

        return picking_surface_submit_readback_and_poll_hit(out_hit, args.picking_surface, args.picking_handler, request);
    }
    return false;
}

InteractionSurfaceViewTransformResult interaction_surface_view_transform_apply(ViewTransform* target, const InteractionSurfaceState& state, const InteractionSurfaceViewTransformArgs& args) {
    ASSERT(target);
    InteractionSurfaceViewTransformResult result = {};
    if (state.active || state.hovered) {
        if (state.selection_mode == InteractionSelectionMode::None) {
            const vec2_t delta = vec_cast(ImGui::GetIO().MouseDelta);
            const vec2_t coord = state.mouse_local;
            const float scroll_delta = ImGui::GetIO().MouseWheel;

            TrackballControllerInput input = {};
            input.rotate_button = ImGui::IsMouseDown(ImGuiMouseButton_Left);
            input.pan_button = ImGui::IsMouseDown(ImGuiMouseButton_Right);
            input.dolly_button = ImGui::IsMouseDown(ImGuiMouseButton_Middle);
            input.mouse_coord_curr = coord;
            input.mouse_coord_prev = coord - delta;
            input.screen_size = state.surface_size;
            input.dolly_delta = scroll_delta;
            input.fov_y = args.camera.fov_y;

            TrackballFlags flags = TrackballFlags_None;
            if (state.active) {
                flags |= TrackballFlags_EnableAllInteractions;
            } else {
                flags |= TrackballFlags_DollyEnabled;
            }

            camera_controller_trackball(target, input, args.trackball_param, flags);

            if (ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left)) {
                result.reset_requested = true;
            }
        }
    }
    return result;
}

void interaction_surface_event_extract(InteractionSurfaceEvent* event, const InteractionSurfaceState& state, const PickingHit& hit) {
    ASSERT(event);
    *event = InteractionSurfaceEvent{};

    event->surface_id   = state.surface_id;
    event->item_id      = state.item_id;
    event->mouse_local  = state.mouse_local;
    event->surface_size = state.surface_size;
    event->region_min   = state.region_min;
    event->region_max   = state.region_max;
    event->hit = hit;

    if (ImGui::IsKeyDown(ImGuiMod_Shift) && state.region_max != state.region_min) {
        event->selection_mode = state.selection_mode;
        event->kind = InteractionSurfaceEventKind::RegionSelect;
        if (state.active) {
            event->region_phase = InteractionSurfaceEventPhase::Update;
        } else if (state.deactivated) {
            event->region_phase = InteractionSurfaceEventPhase::Commit;
        }
    } else if (ImGui::IsKeyDown(ImGuiMod_Shift) && (ImGui::IsMouseReleased(ImGuiMouseButton_Left) || ImGui::IsMouseReleased(ImGuiMouseButton_Right))) {
        event->selection_mode = state.selection_mode;
        event->kind = InteractionSurfaceEventKind::Click;
        if (state.deactivated) {
            event->region_phase = InteractionSurfaceEventPhase::Commit;
        }
    } else if (state.hovered) {
        if (ImGui::IsMouseReleased(ImGuiMouseButton_Right) && ImGui::GetMouseDragDelta(ImGuiMouseButton_Right) == ImVec2(0,0)) {
            event->kind = InteractionSurfaceEventKind::ContextMenu;
        } else {
            event->kind = InteractionSurfaceEventKind::Hover;
        }
    }
}

void point_set_region_mask_compute(md_bitfield_t* mask,
    const float x[],
    const float y[],
    const float z[],
    size_t count,
    const md_bitfield_t* candidate_mask,
    const mat4_t& world_to_clip,
    const vec2_t& region_min,
    const vec2_t& region_max,
    const vec2_t& surface_size)
{
    ASSERT(mask);
    ASSERT(x);
    ASSERT(y);
    ASSERT(z);

    md_bitfield_clear(mask);

    // Transform visible atoms using supplied world_to_clip and set if within region

    if (candidate_mask) {
        // Do for candidate set
        md_bitfield_iter_t it = md_bitfield_iter_create(candidate_mask);

        while (md_bitfield_iter_next(&it)) {
            size_t idx = md_bitfield_iter_idx(&it);
            vec4_t xyz1 = { x[idx], y[idx], z[idx], 1.0f };
            vec4_t coord = mat4_mul_vec4(world_to_clip, xyz1);
            vec2_t surf_coord = {( coord.x / coord.w * 0.5f + 0.5f) * surface_size.x,
                                    (-coord.y / coord.w * 0.5f + 0.5f) * surface_size.y};
            if (region_min.x <= surf_coord.x && surf_coord.x <= region_max.x &&
                region_min.y <= surf_coord.y && surf_coord.y <= region_max.y) {
                md_bitfield_set_bit(mask, idx);
            }
        }
    } else {
        // Do for full set
        for (size_t i = 0; i < count; ++i) {
            vec4_t xyz1 = { x[i], y[i], z[i], 1.0f };
            vec4_t coord = mat4_mul_vec4(world_to_clip, xyz1);
            vec2_t surf_coord = {( coord.x / coord.w * 0.5f + 0.5f) * surface_size.x,
                                    (-coord.y / coord.w * 0.5f + 0.5f) * surface_size.y};
            if (region_min.x <= surf_coord.x && surf_coord.x <= region_max.x &&
                region_min.y <= surf_coord.y && surf_coord.y <= region_max.y) {
                md_bitfield_set_bit(mask, i);
            }
        }
    }
}

bool file_queue_empty(const FileQueue* queue) {
    return queue->head == queue->tail;
}

bool file_queue_full(const FileQueue* queue) {
    return (queue->head + 1) % ARRAY_SIZE(queue->arr) == queue->tail;
}

void file_queue_push(FileQueue* queue, str_t path, FileFlags flags) {
    ASSERT(queue);
    ASSERT(!file_queue_full(queue));
    int prio = 5;


    str_t ext;
    if (extract_ext(&ext, path)) {
        LoaderType type = loader::type_from_ext(ext);
        LoaderFlags loader_flags = loader::type_flags(type);
        if (str_eq(ext, WORKSPACE_FILE_EXTENSION)) {
            prio = 1;
        } else if (loader_flags & LoaderFlag_System) {
            prio = 2;
        } else if (loader_flags & LoaderFlag_Trajectory) {
            prio = 3;
        } else if (find_in_arr(ext, SCRIPT_IMPORT_FILE_EXTENSIONS, ARRAY_SIZE(SCRIPT_IMPORT_FILE_EXTENSIONS))) {
            prio = 4;
        } else {
            flags |= FileFlags_ShowDialogue;
        }
    } else {
        // Unknown extension
        flags |= FileFlags_ShowDialogue;
    }

    uint32_t i = queue->head;
    queue->arr[queue->head] = {str_copy(path, queue->ring), flags, prio};
    queue->head = (queue->head + 1) % ARRAY_SIZE(queue->arr);

    // Sort queue based on prio
     while (i != queue->tail && queue->arr[i].prio < queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)].prio) {
        FileQueue::Entry tmp = queue->arr[i];
        queue->arr[i] = queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)];
        queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)] = tmp;
        i = (i - 1) % ARRAY_SIZE(queue->arr);
     }
}

FileQueue::Entry file_queue_front(const FileQueue* queue) {
    ASSERT(!file_queue_empty(queue));
    return queue->arr[queue->tail];
}

FileQueue::Entry file_queue_pop(FileQueue* queue) {
    ASSERT(queue);
    ASSERT(!file_queue_empty(queue));
    FileQueue::Entry front = file_queue_front(queue);
    queue->tail = (queue->tail + 1) % ARRAY_SIZE(queue->arr);
    return front;
}

void file_queue_process(ApplicationState* state) {
    ASSERT(state);
    if (!file_queue_empty(&state->file_queue) && !state->load_dataset.show_window) {
        FileQueue::Entry e = file_queue_pop(&state->file_queue);

        str_t ext;
        extract_ext(&ext, e.path);
        const str_t* res = 0;

        if (str_eq_ignore_case(ext, WORKSPACE_FILE_EXTENSION)) {
            load_workspace(state, e.path);
            reset_view(&state->view.camera, state->mold.state, &state->representation.visibility_mask);
        } else if ((res = find_in_arr(ext, SCRIPT_IMPORT_FILE_EXTENSIONS, ARRAY_SIZE(SCRIPT_IMPORT_FILE_EXTENSIONS)))) {
            char buf[1024];
            str_t base_path = {};
            if (state->files.workspace[0] != '\0') {
                base_path = str_from_cstr(state->files.workspace);
            } else if (state->files.trajectory[0] != '\0') {
                base_path = str_from_cstr(state->files.trajectory);
            } else if (state->files.molecule[0] != '\0') {
                base_path = str_from_cstr(state->files.molecule);
            } else {
                md_path_write_cwd(buf, sizeof(buf));
                base_path = str_from_cstr(buf);
            }

            str_t rel_path = md_path_make_relative(base_path, e.path, state->allocator.frame);
            MD_LOG_DEBUG("Attempting to make relative path from '" STR_FMT "' to '" STR_FMT "'", STR_ARG(base_path), STR_ARG(e.path));
            MD_LOG_DEBUG("Relative path: '" STR_FMT "'", STR_ARG(rel_path));
            if (str_empty(rel_path)) {
                // No relative path exists between the two (a separate volume on windows)
                rel_path = md_path_make_canonical(e.path, state->allocator.frame);
            }
            if (!str_empty(rel_path)) {
                snprintf(buf, sizeof(buf), "table = import(\"%.*s\");\n", STR_ARG(rel_path));
                TextEditor::Coordinates pos = state->editor.GetCursorPosition();
                pos.mLine += 1;
                state->editor.SetCursorPosition({0,0});
                state->editor.InsertText(buf);
                state->editor.SetCursorPosition(pos);
            }
        } else {
            loader::LoaderState loader_state = {};
            loader::init(&loader_state, e.path, &state->mold.sys);
                
            if ((e.flags & FileFlags_ShowDialogue) || (loader_state.flags & LoaderFlag_RequiresDialogue)) {
                state->load_dataset = LoadDatasetWindowState();
                str_copy_to_char_buf(state->load_dataset.path_buf, sizeof(state->load_dataset.path_buf), e.path);
                state->load_dataset.path_changed = true;
                state->load_dataset.show_window = true;
                state->load_dataset.coarse_grained = e.flags & FileFlags_CoarseGrained;
            } else {
                loader_state.flags |= (e.flags & FileFlags_DisableCacheWrite) ? LoaderFlag_DisableCacheWrite : 0;
                loader_state.flags |= (e.flags & FileFlags_CoarseGrained) ? LoaderFlag_CoarseGrained : 0;
                if (load_data_from_file(state, e.path, loader_state)) {
                    state->animation = {};
                    // @TODO @FIX: This is hacky, just because the loader CAN set system state does not mean it always will.
                    // This should be instead captured and performed by the Event that signals when a new system is loaded.
                    if (loader_state.flags & LoaderFlag_System) {
                        md_bitfield_reset(&state->representation.visibility_mask);

                        if (!state->settings.keep_representations) {
                            remove_all_representations(state);
                            create_default_representations(state);
                        }
                        recompute_atom_visibility_mask(state);
                        state->mold.interpolate_system_state = true;
                        state->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
                        reset_view(&state->view.camera, state->mold.state, &state->representation.visibility_mask);
                    }
                }
            }
        }
    }
}

void reset_view(ViewTransform* transform, const md_system_state_t& state, const md_bitfield_t* mask) {
    ASSERT(transform);
    if (!state.num_atoms) return;

    md_temp_scope_t temp = md_temp_begin();
    defer { md_temp_end(temp); };

    size_t popcount = 0;
    if (mask) {
        popcount = md_bitfield_popcount(mask);
    }

    int32_t* indices = nullptr;
    if (0 < popcount && popcount < state.num_atoms) {
        indices = md_temp_alloc_array(temp, int32_t, popcount);
        size_t len = md_bitfield_iter_extract_indices(indices, popcount, md_bitfield_iter_create(mask));
        if (len > popcount || len > state.num_atoms) {
            MD_LOG_DEBUG("Error: Invalid number of indices");
            len = MIN(popcount, state.num_atoms);
        }
    }

    size_t count = popcount ? popcount : state.num_atoms;
    vec3_t com = md_util_com_compute(state.x, state.y, state.z, nullptr, indices, count, &state.unitcell);

    mat3_t PCA = mat3_ident();
    if (count > 4) {
        mat3_t C = mat3_covariance_matrix(state.x, state.y, state.z, nullptr, indices, count, com);
        mat3_eigen_t eigen = mat3_eigen(C);
        PCA = mat3_orthonormalize(mat3_extract_rotation(eigen.vectors));
    }

    // Compute min and maximum extent along the PCA axes
    vec4_t min_ext = vec4_set1( FLT_MAX);
    vec4_t max_ext = vec4_set1(-FLT_MAX);
    mat4_t Ri  = mat4_from_mat3(PCA);

    // Transform the atom (x,y,z,radius) into the PCA frame to find the min and max extend within it
    for (size_t i = 0; i < count; ++i) {
        int32_t idx = indices ? indices[i] : (int32_t)i;
        vec4_t xyz1 = { state.x[idx], state.y[idx], state.z[idx], 1.0f };

        vec4_t p = mat4_mul_vec4(Ri, xyz1);
        min_ext = vec4_min(min_ext, p);
        max_ext = vec4_max(max_ext, p);
    }

    const float radius = 1.0f;
    min_ext -= vec4_set1(radius);
    max_ext += vec4_set1(radius);

    mat3_t basis = mat3_transpose(PCA);
    vec3_t half_ext = (vec3_from_vec4(max_ext) - vec3_from_vec4(min_ext)) * 0.5f;

    mat3_t A = {};
    md_unitcell_A_extract(A.elem, &state.unitcell);
    mat4_t unitcell_transform = mat4_translate_vec3(-mat3_mul_vec3(A, vec3_set1(0.5f)));
    com = mat4_mul_vec3(unitcell_transform, com, 1.0f);

    ViewTransform opt_view = compute_optimal_view(com, half_ext, basis);

    if (count <= 4) {
        // Apply same as double click camera reset, but center on COM of target, also use optimal distance but keep current orientation
        transform->distance = opt_view.distance;
        transform->position = com + transform->orientation * vec3_set(0, 0, transform->distance);
    } else {
        // Copy full optimal view
        *transform = opt_view;
    }
}

void ViamdEventHandler::process_events(const viamd::Event* events, size_t num_events) {
    for (size_t i = 0; i < num_events; ++i) {
        const viamd::Event& event = events[i];
        switch (event.type) {
        case viamd::EventType_ViamdFrameTick:
            break;
        case viamd::EventType_ViamdPickingRangeReserve: {
			ASSERT(event.payload_type == viamd::EventPayloadType_PickingSpace);
			PickingSpace* space = (PickingSpace*)event.payload;
            size_t num_atoms = state->mold.sys.atom.count;
            size_t num_bonds = state->mold.sys.bond.count;
            size_t num_dipoles = dipole_moments_gather(NULL, 0, state->mold.sys);
            picking_range_reserve(&state->picking_range_atom, space, PickingDomain_Atom, num_atoms);
            picking_range_reserve(&state->picking_range_bond, space, PickingDomain_Bond, num_bonds);
            picking_range_reserve(&state->picking_range_dipole, space, PickingDomain_Dipole, num_dipoles);
            break;
        }
        case viamd::EventType_ViamdInteractionSurface:
            if (event.payload) {
                ASSERT(event.payload_type == viamd::EventPayloadType_InteractionSurfaceEvent);
                InteractionSurfaceEvent* surf = (InteractionSurfaceEvent*)event.payload;
                switch (surf->kind) {
                case InteractionSurfaceEventKind::Hover:
                    draw_picking_tooltip_window(surf->hit, *state);
                    [[fallthrough]];
                case InteractionSurfaceEventKind::Click:
                    // Use highlight mask as intermediate mask for selection operations and to provide hover feedback
                    md_bitfield_clear(&state->selection.highlight_mask);
                    if (surf->hit.domain == PickingDomain_Atom) {
                        int32_t atom_idx = surf->hit.local_idx;
                        if (atom_idx >= 0 && (size_t)atom_idx < state->mold.sys.atom.count) {
                            md_bitfield_set_bit(&state->selection.highlight_mask, atom_idx);
                            if (surf->selection_mode == InteractionSelectionMode::Append) {
                                single_selection_sequence_push_idx(&state->selection.single_selection_sequence, atom_idx);
                            }
                            else if (surf->selection_mode == InteractionSelectionMode::Remove) {
                                single_selection_sequence_pop_idx(&state->selection.single_selection_sequence, atom_idx);
                            }
                        }
                    } else if (surf->hit.domain == PickingDomain_Bond) {
                        size_t bond_idx = surf->hit.local_idx;
                        if (bond_idx < state->mold.sys.bond.count) {
                            md_bitfield_set_bit(&state->selection.highlight_mask, state->mold.sys.bond.pairs[bond_idx].idx[0]);
                            md_bitfield_set_bit(&state->selection.highlight_mask, state->mold.sys.bond.pairs[bond_idx].idx[1]);
                        }
                    }
                    
                    // Commit to selection mask upon click release, for hover we only update the highlight mask
                    if (surf->hit.domain == PickingDomain_Atom || surf->hit.domain == PickingDomain_Bond) {
                        grow_mask_by_selection_granularity(&state->selection.highlight_mask, state->selection.granularity, state->mold.sys);
                        if (surf->selection_mode == InteractionSelectionMode::Append) {
                            md_bitfield_or_inplace(&state->selection.selection_mask, &state->selection.highlight_mask);
                        }
                        else if (surf->selection_mode == InteractionSelectionMode::Remove) {
                            md_bitfield_andnot_inplace(&state->selection.selection_mask, &state->selection.highlight_mask);
                        }   
                    } else if (surf->hit.domain == 0) {
                        if (surf->selection_mode == InteractionSelectionMode::Remove) {
                            md_bitfield_clear(&state->selection.selection_mask);
                            single_selection_sequence_clear(&state->selection.single_selection_sequence);
                        }
                    }

                    break;
                case InteractionSurfaceEventKind::RegionSelect:
                    break;
                case InteractionSurfaceEventKind::ContextMenu:
                    break;
                case InteractionSurfaceEventKind::None: [[fallthrough]];
                default:
                    break;
                }
            }
            break;
        case viamd::EventType_ViamdPickingTooltipTextRequest: {
            ASSERT(event.payload_type == viamd::EventPayloadType_PickingTooltipTextRequest);
            PickingTooltipTextRequest* req = (PickingTooltipTextRequest*)event.payload;
            if (req->hit.domain == PickingDomain_Atom || req->hit.domain == PickingDomain_Bond || req->hit.domain == PickingDomain_Dipole) {
                fill_picking_tooltip_text(&req->sb, *state, req->hit);
            }
            break;
        }
        case viamd::EventType_ViamdViewFit: {
            ASSERT(event.payload_type == viamd::EventPayloadType_ViewFitRequest);
            ViewFitRequest* req = (ViewFitRequest*)event.payload;
            if (req) {
                md_bitfield_t* bf = nullptr;
                switch (req->round) {
                case ViewFitRound_Highlight:
                    bf = &state->selection.highlight_mask; break;
                case ViewFitRound_Selection:
                    bf = &state->selection.selection_mask; break;
                case ViewFitRound_Visible:
                    bf = &state->representation.visibility_mask; break;
                default:
                    break;
                }

                size_t popcount = bf ? md_bitfield_popcount(bf) : 0;
                if (popcount > 0) {
                    vec4_t* dst_xyzw = md_array_extend(req->xyzw, popcount, req->alloc);
                    if (dst_xyzw) {
                        md_util_system_extract_xyzw_from_mask(dst_xyzw, bf, &state->mold.sys, &state->mold.state);
                    }
                }
            }
            break;
        }
        case viamd::EventType_ViamdSystemStateChanged: {
			// Apply operators if toggled on system change

            // EXTRACT STATE HERE FROM PAYLOAD
            ASSERT(event.payload_type == viamd::EventPayloadType_ApplicationState);
            ApplicationState* app = (ApplicationState*)event.payload;

            int num_tasks = 0;
            task_system::ID tasks[16];
            
            md_system_t& sys = app->mold.sys;
			md_system_state_t& sys_state = app->mold.state;
            // Identity, not the zero matrix: a transform that never gets computed must leave
            // the system where it is rather than collapse it onto the origin.
            mat4_t recenter_transform = mat4_ident();

            // Whether any operation below actually rewrote mold.state coordinates.
            bool coords_modified = false;

            if (app->operations.recalc_bonds) {
                static int64_t cur_nearest_frame = -1;

                // We cannot recalculate bonds while the full or filtered evaluation is running
                // because it would overwrite the bond data while we are reading it
                int64_t nearest_frame = (int64_t)(app->animation.frame + 0.5);
                if (!task_system::task_is_running(app->tasks.evaluate_full) && !task_system::task_is_running(app->tasks.evaluate_filt)) {
                    if (app->mold.sys.trajectory == NULL || (cur_nearest_frame != nearest_frame)) {
                        cur_nearest_frame = nearest_frame;
                        task_system::ID recalc_bond_task = task_system::create_pool_task(STR_LIT("## Recalc bond task"), [&sys, app, &nearest_frame]() {
                            md_temp_scope_t temp = md_temp_begin();
                            defer { md_temp_end(temp); };

							md_system_state_t ref_state = sys.reference;

                            if (sys.trajectory) {
                                // Use state from frame closest to the current animation time
                                md_system_state_t frame_state = { .alloc = temp.arena };
								md_system_state_init(&frame_state, sys.atom.count);
                                if (!md_trajectory_load_frame(sys.trajectory, nearest_frame, &frame_state)) {
                                    MD_LOG_ERROR("Failed to extract frame data");
                                }
                            }

                            MD_LOG_DEBUG("RECALCULATING BONDS");
                            md_util_infer_covalent_bonds(&sys.bond, &ref_state, &sys, sys.alloc);
                            md_bond_build_connectivity(&sys.bond, sys.atom.count, sys.alloc);

                            app->mold.dirty_gpu_buffers |= MolBit_DirtyBonds;
                        });
                        tasks[num_tasks++] = recalc_bond_task;
                    }
                }
            }

            if (state->operations.recenter) {
                const md_bitfield_t& target_mask = recenter_get_active_target_mask(state);
                size_t num_idx = md_bitfield_popcount(&target_mask);
                if (num_idx > 0) {
                    // Create async task to calculate transformation matrix (Its only expressed as a task to ensure that it runs after some of the previous tasks in the workflow)
                    task_system::ID calc_transform_task = task_system::create_pool_task(STR_LIT("## Calculate Recenter Transform"), [&recenter_transform, app]() {
                        recenter_calculate_transform(recenter_transform.elem, app);
                    });

                    // Batch transform all atoms
                    task_system::ID apply_transform_task = task_system::create_pool_task(STR_LIT("## Recenter"), (uint32_t)sys.atom.count, [&sys_state, &recenter_transform](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                        (void)thread_num;
                        size_t count = range_end - range_beg;
                        float* x = sys_state.x + range_beg;
                        float* y = sys_state.y + range_beg;
                        float* z = sys_state.z + range_beg;
                        mat4_batch_transform_inplace(x, y, z, 1.0f, count, recenter_transform);
                    }, 1024);

                    tasks[num_tasks++] = calc_transform_task;
                    tasks[num_tasks++] = apply_transform_task;
                    coords_modified = true;
                }
            }

            if (state->operations.apply_pbc) {
                task_system::ID pbc_task = task_system::create_pool_task(STR_LIT("## Apply PBC"), (uint32_t)sys.atom.count, [&sys_state](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    size_t count = range_end - range_beg;
                    float* x = sys_state.x + range_beg;
                    float* y = sys_state.y + range_beg;
                    float* z = sys_state.z + range_beg;
                    md_util_pbc(x, y, z, NULL, count, &sys_state.unitcell);
                });
                tasks[num_tasks++] = pbc_task;
                coords_modified = true;
            } 

            if (state->operations.unwrap_structures) {
                size_t num_structures = md_structure_count(&sys.structure);
                task_system::ID unwrap_task = task_system::create_pool_task(STR_LIT("## Unwrap Structures"), (uint32_t)num_structures, [&sys_state, &sys](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    for (uint32_t i = range_beg; i < range_end; ++i) {
                        md_structure_t structure = {};
                        md_structure_extract(&structure, &sys.structure, i);
						md_util_unwrap_structure(&sys_state, &structure);
                    }
                });
                tasks[num_tasks++] = unwrap_task;
                coords_modified = true;
            }

            if (num_tasks > 0) {
                for (int j = 1; j < num_tasks; ++j) {
                    task_system::set_task_dependency(tasks[j], tasks[j-1]);
                }
                task_system::enqueue_task(tasks[0]);
                task_system::task_wait_for(tasks[num_tasks - 1]);
            }

            // The operations above rewrote the coordinates that update_md_buffers uploads.
            // Nothing else flags them: the synchronous broadcasts of this event do not pass
            // through the interpolation step that would otherwise have set the bit.
            if (coords_modified) {
                app->mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
            }
            break;
        }
        default:
            break;
        }
    }
}

void script_visualize_payload(ApplicationState* state, const md_script_vis_payload_o* payload, int subidx, md_script_vis_flags_t flags) {
    ASSERT(state);

    md_script_vis_ctx_t ctx = {
        .ir   = state->script.eval_ir,
        .sys  = &state->mold.sys,
        .state = &state->mold.state,
    };

    if (md_script_vis_eval_payload(&state->script.vis, payload, subidx, &ctx, flags)) {
        if (!md_bitfield_empty(&state->script.vis.atom_mask)) {
            md_bitfield_copy(&state->selection.highlight_mask, &state->script.vis.atom_mask);
        }
    }
}

void script_visualize_str(ApplicationState* state, str_t str, md_script_vis_flags_t flags) {
    ASSERT(state);

    md_script_vis_ctx_t ctx = {
        .ir    = state->script.eval_ir,
        .sys   = &state->mold.sys,
        .state = &state->mold.state,
    };

    if (md_script_vis_eval_string(&state->script.vis, str, &ctx, flags)) {
        if (!md_bitfield_empty(&state->script.vis.atom_mask)) {
            md_bitfield_copy(&state->selection.highlight_mask, &state->script.vis.atom_mask);
        }
    }
}

void script_set_hovered_property(ApplicationState* state, str_t label, int population_idx) {
    state->hovered_display_property_label = label;
    state->hovered_display_property_pop_idx = population_idx;
}
