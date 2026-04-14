
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

static void init_all_representations(ApplicationState* state);

static void fill_picking_tooltip_text(md_strb_t* sb, const ApplicationState& state, const PickingHit& hit) {
    ASSERT(sb);
    const md_system_t& sys = state.mold.sys;

    if (hit.domain == PickingDomain_Atom && hit.local_idx < sys.atom.count) {
        int atom_idx = hit.local_idx;
        int local_idx = atom_idx;
        const vec3_t pos = md_atom_coord(&sys.atom, atom_idx);
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
            md_strb_fmt(sb, "comp[%i]", comp_idx + 1);
            if (comp_name) {
                md_strb_fmt(sb, ": " STR_FMT, STR_ARG(comp_name));
            }
			md_strb_fmt(sb, " (seq_id: %i)\n", comp_seq_id);
        }
        if (inst_idx != -1) {
            md_strb_fmt(sb, "inst[%i]", inst_idx + 1);
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

            vec3_t p0 = md_atom_coord(&sys.atom, pair.idx[0]);
            vec3_t p1 = md_atom_coord(&sys.atom, pair.idx[1]);
            float d = vec3_distance(p0, p1);

            str_t type0 = md_atom_name(&sys.atom, pair.idx[0]);
            str_t type1 = md_atom_name(&sys.atom, pair.idx[1]);

            md_strb_fmt(sb, "bond: " STR_FMT "%c" STR_FMT "\n", STR_ARG(type0), bond_type, STR_ARG(type1));
            md_strb_fmt(sb, "flags: %.*s\n", len, bond_flags_buf);
            md_strb_fmt(sb, "length: %.3f\n", d);
        }
    }
}

void draw_picking_tooltip_window(const PickingHit& hit, const ApplicationState& state) {
    if (hit.raw_idx == INVALID_PICKING_IDX) return;
    
	md_vm_arena_temp_t temp = md_vm_arena_temp_begin(state.allocator.frame);
    defer { md_vm_arena_temp_end(temp); };

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
        md_array_resize(cache->states[i].atom_x, capacity, alloc);
        md_array_resize(cache->states[i].atom_y, capacity, alloc);
        md_array_resize(cache->states[i].atom_z, capacity, alloc);
    }
}

static inline void free_frame_cache(FrameCache* cache, md_allocator_i* alloc) {
    for (size_t i = 0; i < FRAME_CACHE_SIZE; ++i) {
        md_array_free(cache->states[i].atom_x, alloc);
        md_array_free(cache->states[i].atom_y, alloc);
        md_array_free(cache->states[i].atom_z, alloc);
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

// #trajectorydata
void free_trajectory_data(ApplicationState* state) {
    ASSERT(state);

    if (state->mold.sys.trajectory) {
        md_trajectory_free(state->mold.sys.trajectory);
        state->mold.sys.trajectory = nullptr;
    }
    state->files.trajectory[0] = '\0';

    md_array_free(state->timeline.x_values,  state->allocator.persistent);
    md_array_free(state->display_properties, state->allocator.persistent);

    md_array_free(state->trajectory_data.backbone_angles.data,     state->allocator.persistent);
    md_array_free(state->trajectory_data.secondary_structure.data, state->allocator.persistent);

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

        data->animation.frame = CLAMP(data->animation.frame, (double)min_frame, (double)max_frame);
        int64_t frame_idx = CLAMP((int64_t)(data->animation.frame + 0.5), 0, (int64_t)max_frame);

        md_trajectory_frame_header_t frame_header = {};
        md_trajectory_load_frame(data->mold.sys.trajectory, frame_idx, &frame_header, data->mold.sys.atom.x, data->mold.sys.atom.y, data->mold.sys.atom.z);
        data->mold.sys.unitcell = frame_header.unitcell;

        if (data->mold.sys.protein_backbone.segment.count > 0) {
            data->trajectory_data.secondary_structure.stride = data->mold.sys.protein_backbone.segment.count;
            data->trajectory_data.secondary_structure.count = data->mold.sys.protein_backbone.segment.count * num_frames;
            md_array_resize(data->trajectory_data.secondary_structure.data, data->mold.sys.protein_backbone.segment.count * num_frames, data->allocator.persistent);
            MEMSET(data->trajectory_data.secondary_structure.data, 0, md_array_bytes(data->trajectory_data.secondary_structure.data));

            data->trajectory_data.backbone_angles.stride = data->mold.sys.protein_backbone.segment.count;
            data->trajectory_data.backbone_angles.count = data->mold.sys.protein_backbone.segment.count * num_frames;
            md_array_resize(data->trajectory_data.backbone_angles.data, data->mold.sys.protein_backbone.segment.count * num_frames, data->allocator.persistent);
            MEMSET(data->trajectory_data.backbone_angles.data, 0, md_array_bytes(data->trajectory_data.backbone_angles.data));

            // Launch work to compute the values
            task_system::task_interrupt_and_wait_for(data->tasks.backbone_computations);

            data->tasks.backbone_computations = task_system::create_pool_task(STR_LIT("Backbone Operations"), (uint32_t)num_frames, [data](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                // Create copy here of molecule since we use the full structure as input
                md_system_t sys = data->mold.sys;
                md_trajectory_i* traj = data->mold.sys.trajectory;

                md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(1));
                defer { md_vm_arena_destroy(temp_arena); };

                const size_t capacity = ALIGN_TO(sys.atom.count, 16);
                float* x = (float*)md_vm_arena_push(temp_arena, sizeof(float) * capacity);
                float* y = (float*)md_vm_arena_push(temp_arena, sizeof(float) * capacity);
                float* z = (float*)md_vm_arena_push(temp_arena, sizeof(float) * capacity);

                md_trajectory_reader_i reader;
                if (md_trajectory_reader_init(&reader, traj)) {
                    for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                        md_backbone_angles_t* bb_dst = data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx;
                        md_secondary_structure_t* ss_dst = data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx;
                        
                        md_trajectory_frame_header_t frame_header = {0};
                        md_trajectory_reader_load_frame(reader, frame_idx, &frame_header, x, y, z);
                        md_util_backbone_angles_compute(bb_dst, data->trajectory_data.backbone_angles.stride, x, y, z, &frame_header.unitcell, &sys.protein_backbone);
                        md_util_backbone_secondary_structure_infer(ss_dst, data->trajectory_data.secondary_structure.stride, x, y, z, &frame_header.unitcell, &sys.protein_backbone);
                    }
					md_trajectory_reader_free(&reader);
                } else {
                    for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                        md_backbone_angles_t* bb_dst = data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx;
                        md_secondary_structure_t* ss_dst = data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx;
                        
                        md_trajectory_frame_header_t frame_header = {0};
                        md_trajectory_load_frame(traj, frame_idx, &frame_header, x, y, z);
                        md_util_backbone_angles_compute(bb_dst, data->trajectory_data.backbone_angles.stride, x, y, z, &frame_header.unitcell, &sys.protein_backbone);
                        md_util_backbone_secondary_structure_infer(ss_dst, data->trajectory_data.secondary_structure.stride, x, y, z, &frame_header.unitcell, &sys.protein_backbone);
                    }
                }
            });

            uint64_t time = (uint64_t)md_time_now();
            task_system::ID main_task = task_system::create_main_task(STR_LIT("Update Trajectory Data"), [data, t0 = time]() {
                uint64_t t1 = (uint64_t)md_time_now();
                double elapsed = md_time_as_seconds(t1 - t0);
                MD_LOG_INFO("Finished computing trajectory data (%.2fs)", elapsed);
                data->trajectory_data.backbone_angles.fingerprint     = generate_fingerprint();
                data->trajectory_data.secondary_structure.fingerprint = generate_fingerprint();

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
        md_bitfield_clear(&data->operations.target_mask);

        data->mold.gl_mol = md_gl_mol_create(&data->mold.sys);
        if (data->mold.sys.protein_backbone.segment.count > 0) {
            data->interpolated_properties.secondary_structure = md_array_create(md_gl_secondary_structure_t, data->mold.sys.protein_backbone.segment.count, data->mold.sys_alloc);
        }

        mat3_t A;
		md_unitcell_A_extract_float(A.elem, &data->mold.sys.unitcell);
		vec3_t c = mat3_mul_vec3(A, vec3_set(0.5f, 0.5f, 0.5f));
        data->mold.unitcell_transform = mat4_translate(-c.x, -c.y, -c.z);

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

static inline void clear_density_volume(ApplicationState* state) {
    md_array_shrink(state->density_volume.gl_reps, 0);
    md_array_shrink(state->density_volume.rep_model_mats, 0);
    state->density_volume.model_mat = {0};
}

void free_system_data(ApplicationState* data) {
    ASSERT(data);
    interrupt_async_tasks(data);

    //md_molecule_free(&data->mold.sys, persistent_alloc);
    md_arena_allocator_reset(data->mold.sys_alloc);
    MEMSET(&data->mold.sys, 0, sizeof(data->mold.sys));

    md_array_free(data->operations.initial_frame.xyzw, data->allocator.persistent);
    data->operations.initial_frame.xyzw = nullptr;

    md_gl_mol_destroy(data->mold.gl_mol);
    MEMSET(data->files.molecule, 0, sizeof(data->files.molecule));

    data->interpolated_properties.secondary_structure = nullptr;

    MEMSET(data->mold.frame_cache.states, 0, sizeof(data->mold.frame_cache.states));
    clear_frame_cache(&data->mold.frame_cache);

    md_bitfield_clear(&data->selection.selection_mask);
    md_bitfield_clear(&data->selection.highlight_mask);
    if (data->script.ir) {
        md_script_ir_free(data->script.ir);
        data->script.ir = nullptr;
    }
    if (data->script.eval_ir && data->script.eval_ir != data->script.ir) {
        md_script_ir_free(data->script.eval_ir);
        data->script.eval_ir = nullptr;
    }
    if (data->script.full_eval) {
        md_script_eval_free(data->script.full_eval);
        data->script.full_eval = nullptr;
    }
    if (data->script.filt_eval) {
        md_script_eval_free(data->script.filt_eval);
        data->script.filt_eval = nullptr;
    }
    clear_density_volume(data);

    viamd::event_system_broadcast_event(viamd::EventType_ViamdSystemFree, viamd::EventPayloadType_ApplicationState, data);
}

bool load_data_from_file(ApplicationState* state, str_t filepath, const loader::State& load_state) {
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

            state->mold.sys.alloc = state->mold.sys_alloc;
            if (!loader::load(&state->mold.sys, path_to_file, load_state)) {
                VIAMD_LOG_ERROR("Failed to load molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));
                return false;
            }
            success = true;
            VIAMD_LOG_SUCCESS("Successfully loaded molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));

            str_copy_to_char_buf(state->files.molecule, sizeof(state->files.molecule), path_to_file);
            state->files.coarse_grained = load_state.flags & LoaderFlag_CoarseGrained;
            // @NOTE: If the dataset is coarse-grained, then postprocessing must be aware
            md_postprocess_flags_t flags = state->files.coarse_grained ? MD_UTIL_POSTPROCESS_NONE : MD_UTIL_POSTPROCESS_ALL;
            md_util_system_postprocess(&state->mold.sys, flags);
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

            success = loader::load(&state->mold.sys, path_to_file, load_state);
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
    ScopedTemp temp_reset;
    md_allocator_i* temp_alloc = md_get_temp_allocator();

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
                        path += folder;
                        path += file;
                        new_molecule_file = md_path_make_canonical(path, temp_alloc);
                    }
                } else if (str_eq(ident, STR_LIT("TrajectoryFile"))) {
                    str_t file;
                    viamd::extract_str(file, arg);
                    if (!str_empty(file)) {
                        md_strb_t path = md_strb_create(temp_alloc);
                        path += folder;
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
                } else if (str_eq(ident, STR_LIT("ElectronicStructureMoIdx"))) {
                    viamd::extract_int(rep->electronic_structure.mo_idx, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureNtoIdx"))) {
                    viamd::extract_int(rep->electronic_structure.nto_idx, arg);
                } else if (str_eq(ident, STR_LIT("ElectronicStructureRes"))) {
                    int res;
                    viamd::extract_int(res, arg);
                    rep->electronic_structure.resolution = (VolumeResolution)res;
                } else if (str_eq(ident, STR_LIT("ElectronicStructureType"))) {
                    int type;
                    viamd::extract_int(type, arg);
                    rep->electronic_structure.type = (ElectronicStructureType)type;
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
                }
            }
        }
        /*else if (str_eq(section, STR_LIT("AtomElementMapping"))) {
            str_t lbl = {};
            int elem = 0;
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("Label"))) {
                    viamd::extract_str(lbl, arg);
                } else if (str_eq(ident, STR_LIT("Element"))) {
                    viamd::extract_int(elem, arg);
                }
            }
            if (!str_empty(lbl) && elem) {
                add_atom_elem_mapping(data, lbl, (md_element_t)elem);
            }
        } */
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

    loader::State loader_state = {};
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

    str_t mol_file  = md_path_make_relative(filename, str_from_cstr(app_state->files.molecule),   temp_alloc);
    str_t traj_file = md_path_make_relative(filename, str_from_cstr(app_state->files.trajectory), temp_alloc);

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
            viamd::write_int(state,  STR_LIT("ElectronicStructureMoIdx"),    rep.electronic_structure.mo_idx);
            viamd::write_int(state,  STR_LIT("ElectronicStructureNtoIdx"),   rep.electronic_structure.nto_idx);
            viamd::write_int(state,  STR_LIT("ElectronicStructureType"),(int)rep.electronic_structure.type);
            viamd::write_int(state,  STR_LIT("ElectronicStructureRes"), (int)rep.electronic_structure.resolution);
            viamd::write_dbl(state,  STR_LIT("ElectronicStructureIso"),      rep.electronic_structure.iso_value);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColPos"),   rep.electronic_structure.col_psi_pos);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColNeg"),   rep.electronic_structure.col_psi_neg);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColDen"),   rep.electronic_structure.col_den);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColAtt"),   rep.electronic_structure.col_att);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColDet"),   rep.electronic_structure.col_det);
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

void remove_selection(ApplicationState* state, int idx) {
    ASSERT(state);
    if (idx < 0 || (int)md_array_size(state->selection.stored_selections) <= idx) {
        VIAMD_LOG_ERROR("Index [%i] out of range when trying to remove selection", idx);
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

    size_t num_props = md_array_size(state->representation.info.atom_properties);
    if (num_props > 0) {
        rep->atomic_property.idx = 0;
        rep->atomic_property.range_beg = state->representation.info.atom_properties[0].value_min;
        rep->atomic_property.range_end = state->representation.info.atom_properties[0].value_max;
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
    rep->electronic_structure.mo_idx = (int)state->representation.info.alpha.homo_idx;
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

void remove_representation(ApplicationState* state, int idx) {
    ASSERT(state);
    ASSERT(idx < md_array_size(state->representation.reps));
    auto& rep = state->representation.reps[idx];
    md_bitfield_free(&rep.atom_mask);
    md_gl_rep_destroy(rep.md_rep);
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

void update_all_representations(ApplicationState* state) {
    for (size_t i = 0; i < md_array_size(state->representation.reps); ++i) {
        update_representation(state, &state->representation.reps[i]);
    }
}

bool representation_uses_atom_colors(const Representation& rep) {
    return rep.type != RepresentationType::ElectronicStructure || rep.electronic_structure.use_atom_colors;
}

void update_representation(ApplicationState* state, Representation* rep) {
    ASSERT(state);
    ASSERT(rep);

    if (!rep->enabled) return;
    if (!rep->needs_update) return;

    const auto& sys = state->mold.sys;
    size_t num_atoms = md_system_atom_count(&sys);

    md_allocator_i* frame_alloc = state->allocator.frame;
    md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(frame_alloc);
    defer { md_vm_arena_temp_end(tmp); };

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

            if (md_array_size(state->representation.info.atom_properties) > 0) {
                float* values = (float*)md_vm_arena_push(frame_alloc, sizeof(float) * num_atoms);
                EvalAtomProperty eval = {
                    .key = state->representation.info.atom_properties[rep->atomic_property.idx].key,
                    .idx = rep->atomic_property.sub_idx,
                    .num_values = num_atoms,
                    .dst_values = values,
                    .output_written = false,
                };
                viamd::event_system_broadcast_event(viamd::EventType_ViamdRepresentationEvalAtomProperty, viamd::EventPayloadType_EvalAtomProperty, &eval);

                if (eval.output_written) {
                    float range_ext = (rep->atomic_property.range_end - rep->atomic_property.range_beg);
                    range_ext = MAX(range_ext, 0.001f);
                    for (size_t i = 0; i < num_atoms; ++i) {
                        float t = (values[i] - rep->atomic_property.range_beg) / range_ext;
                        colors[i] = ImPlot::SampleColormapU32(ImClamp(t, 0.0f, 1.0f), rep->atomic_property.colormap);
                    }
                } else {
                    MD_LOG_DEBUG("No output written for EvalAtomProperty event");
                    MEMSET(colors, 0xFFFFFFFFu, bytes);
                }
            } else {
                MEMSET(colors, 0xFFFFFFFFu, bytes);
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

    if (rep->tint_scale > 0.0f || rep->saturation < 1.0f) {
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
        size_t num_mos = state->representation.info.alpha.num_orbitals;
        rep->type_is_valid = num_mos > 0;
        if (rep->type_is_valid && rep->enabled) {
            EvalElectronicStructure data = {
                .sys = &state->mold.sys,
                .frame = state->animation.frame,
                .rep = rep,
                .atom_colors = colors,
            };
            viamd::event_system_broadcast_event(viamd::EventType_ViamdRepresentationEvalElectronicStructure, viamd::EventPayloadType_EvalElectronicStructure, &data);
        }
        break;
    }
    default:
        ASSERT(false);
        break;
    }

    if (colors) {
        if (rep->dynamic_evaluation) {
            rep->filt_is_dirty = true;
        }

        if (rep->filt_is_dirty) {
			rep->filt_is_valid = md_filter(&rep->atom_mask, str_from_cstr(rep->filt), &state->mold.sys, state->mold.sys.atom.x, state->mold.sys.atom.y, state->mold.sys.atom.z, state->script.ir, &rep->filt_is_dynamic, rep->filt_error, sizeof(rep->filt_error));
            rep->filt_is_dirty = false;
        }

        if (rep->filt_is_valid) {
            filter_colors(colors, num_atoms, &rep->atom_mask);
            state->representation.atom_visibility_mask_dirty = true;
            md_gl_rep_set_color(rep->md_rep, 0, (uint32_t)num_atoms, colors, 0);

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
    }

    recompute_atom_visibility_mask(state);
}

void interpolate_system_state(ApplicationState* state) {
    ASSERT(state);
    auto& sys = state->mold.sys;
    const auto& traj = state->mold.sys.trajectory;

    if (!sys.atom.count || !md_trajectory_num_frames(traj)) return;

    const int64_t last_frame = MAX(0LL, (int64_t)md_trajectory_num_frames(traj) - 1);
    // This is not actually time, but the fractional frame representation
    const double time = CLAMP(state->animation.frame, 0.0, double(last_frame));

    // Scaling factor for cubic spline
    const int64_t frame = (int64_t)time;
    const int64_t nearest_frame = CLAMP((int64_t)(time + 0.5), 0LL, last_frame);

    static int64_t curr_nearest_frame = -1;
    if (state->animation.interpolation == InterpolationMode::Nearest) {
        if (curr_nearest_frame == nearest_frame) {
            return;
        }
        state->mold.dirty_gpu_buffers |= MolBit_ClearVelocity;
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

    md_allocator_i* temp_arena = state->allocator.frame;
    md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(temp_arena);
    defer { md_vm_arena_temp_end(tmp); };

    struct Payload {
        ApplicationState* state;
        float s;
        float t;
        InterpolationMode mode;

        int64_t nearest_frame;
        int64_t frames[4];
        md_unitcell_t unitcell;

        size_t count;

        md_system_state_t* src_states[4];

        float* dst_x;
        float* dst_y;
        float* dst_z;

        vec3_t* aabb_min;
        vec3_t* aabb_max;

        mat4_t recenter_transform;
    };

    const InterpolationMode mode = (frames[1] == frames[2]) ? InterpolationMode::Nearest : state->animation.interpolation;

    Payload payload = {
        .state = state,
        .s = 1.0f - CLAMP(state->animation.tension, 0.0f, 1.0f),
        .t = (float)fract(time),
        .mode = mode,
        .nearest_frame = nearest_frame,
        .frames = { frames[0], frames[1], frames[2], frames[3]},
        .count = sys.atom.count,
        .dst_x = sys.atom.x,
        .dst_y = sys.atom.y,
        .dst_z = sys.atom.z,
        .aabb_min = (vec3_t*)md_vm_arena_push(temp_arena, num_threads * sizeof(vec3_t)),
        .aabb_max = (vec3_t*)md_vm_arena_push(temp_arena, num_threads * sizeof(vec3_t)),
    };

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
        if (!find_frame_in_cache(&slot_idx, requested_frames[i], &state->mold.frame_cache)) {
            slot_idx = find_lru_cache_slot(&state->mold.frame_cache);
            frame_cache_load_slot[num_frames_to_load++] = slot_idx;
        }
        set_mru_cache_slot(&state->mold.frame_cache, slot_idx);
        state->mold.frame_cache.frame_idx[slot_idx] = requested_frames[i];
        frame_cache_slot_idx[i] = slot_idx;
    }

    for (int i = 0; i < num_requested_frames; ++i) {
        int slot_idx = frame_cache_slot_idx[i];
        payload.src_states[i] = &state->mold.frame_cache.states[slot_idx];
    }

    // This holds the chain of tasks we are about to submit
    task_system::ID tasks[16] = {0};
    int num_tasks = 0;

    task_system::ID load_task = task_system::create_pool_task(STR_LIT("## Load Frames"), num_frames_to_load,
        [data = &payload, frame_cache_load_slot](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            for (uint32_t i = range_beg; i < range_end; ++i) {
                int slot_idx = frame_cache_load_slot[i];
                int frame_idx = data->state->mold.frame_cache.frame_idx[slot_idx];
                md_system_state_t* state = &data->state->mold.frame_cache.states[slot_idx];
                md_trajectory_frame_header_t header = {0};
                md_trajectory_load_frame(data->state->mold.sys.trajectory, frame_idx, &header, state->atom_x, state->atom_y, state->atom_z);
                state->unitcell = header.unitcell;
            }
        }
    );

    tasks[num_tasks++] = load_task;

    switch (mode) {
        case InterpolationMode::Nearest: {
            task_system::ID interp_task = task_system::create_pool_task(STR_LIT("## Interpolate"), [data = &payload]() {
                data->unitcell = data->src_states[0]->unitcell;
                MEMCPY(data->dst_x, data->src_states[0]->atom_x, sizeof(float) * data->count);
                MEMCPY(data->dst_y, data->src_states[0]->atom_y, sizeof(float) * data->count);
                MEMCPY(data->dst_z, data->src_states[0]->atom_z, sizeof(float) * data->count);
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
                data->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
			});

            task_system::ID interp_coord_task = task_system::create_pool_task(STR_LIT("## Interp Coord Data"), (uint32_t)sys.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                float* dst_x = data->dst_x + range_beg;
                float* dst_y = data->dst_y + range_beg;
                float* dst_z = data->dst_z + range_beg;
                const float* src_x[2] = { data->src_states[0]->atom_x + range_beg, data->src_states[1]->atom_x + range_beg};
                const float* src_y[2] = { data->src_states[0]->atom_y + range_beg, data->src_states[1]->atom_y + range_beg};
                const float* src_z[2] = { data->src_states[0]->atom_z + range_beg, data->src_states[1]->atom_z + range_beg};

                md_util_interpolate_linear(dst_x, dst_y, dst_z, src_x, src_y, src_z, count, &data->unitcell, data->t);
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
                data->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
            });

            task_system::ID interp_coord_task = task_system::create_pool_task(STR_LIT("## Interp Coord Data"), (uint32_t)sys.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                float* dst_x = data->dst_x + range_beg;
                float* dst_y = data->dst_y + range_beg;
                float* dst_z = data->dst_z + range_beg;
                const float* src_x[4] = { data->src_states[0]->atom_x + range_beg, data->src_states[1]->atom_x + range_beg, data->src_states[2]->atom_x + range_beg, data->src_states[3]->atom_x + range_beg};
                const float* src_y[4] = { data->src_states[0]->atom_y + range_beg, data->src_states[1]->atom_y + range_beg, data->src_states[2]->atom_y + range_beg, data->src_states[3]->atom_y + range_beg};
                const float* src_z[4] = { data->src_states[0]->atom_z + range_beg, data->src_states[1]->atom_z + range_beg, data->src_states[2]->atom_z + range_beg, data->src_states[3]->atom_z + range_beg};

                md_util_interpolate_cubic_spline(dst_x, dst_y, dst_z, src_x, src_y, src_z, count, &data->unitcell, data->t, data->s);
            }, grain_size);

            tasks[num_tasks++] = iterp_cell_task;
            tasks[num_tasks++] = interp_coord_task;
            
            break;
        }
        default:
            ASSERT(false);
            break;
    }

    if (state->operations.recalc_bonds) {
        // We cannot recalculate bonds while the full or filtered evaluation is running
        // because it would overwrite the bond data while we are reading it
        if (!task_system::task_is_running(state->tasks.evaluate_full) && !task_system::task_is_running(state->tasks.evaluate_filt)) {
            static int64_t cur_nearest_frame = -1;
            if (cur_nearest_frame != payload.nearest_frame) {
                cur_nearest_frame = payload.nearest_frame;
                task_system::ID recalc_bond_task = task_system::create_pool_task(STR_LIT("## Recalc bond task"), [data = &payload]() {
                    const auto& sys = data->state->mold.sys;
                    const md_system_state_t* src_state = data->src_states[0];
                    int offset = data->t < 0.5 ? 0 : 1;

                    switch (data->mode) {
                    case InterpolationMode::Nearest: break;
                    case InterpolationMode::Linear:
                        src_state = data->src_states[0 + offset];
                        break;
                    case InterpolationMode::CubicSpline:
                        src_state = data->src_states[1 + offset];
                        break;
                    default:
                        break;
                    };

                    md_util_infer_covalent_bonds(&data->state->mold.sys.bond, src_state->atom_x, src_state->atom_y, src_state->atom_z, &src_state->unitcell, &sys, sys.alloc);
                    data->state->mold.dirty_gpu_buffers |= MolBit_DirtyBonds;
                });
                tasks[num_tasks++] = recalc_bond_task;
            }
        }
    }

    if (state->operations.recenter) {
        size_t num_idx = md_bitfield_popcount(&state->operations.target_mask);
        if (num_idx > 0) {
            // Create async task to calculate transformation matrix (Its only expressed as a task to ensure that it runs after some of the previous tasks in the workflow)
            task_system::ID calc_transform_task = task_system::create_pool_task(STR_LIT("## Calculate Recenter Transform"), [data = &payload]() {
                recenter_calculate_transform(data->recenter_transform.elem, data->state);
            });
            
            // Batch transform all atoms
            task_system::ID apply_transform_task = task_system::create_pool_task(STR_LIT("## Recenter"), (uint32_t)sys.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                float* x = data->dst_x + range_beg;
                float* y = data->dst_y + range_beg;
                float* z = data->dst_z + range_beg;
                mat4_batch_transform_inplace(x, y, z, 1.0f, count, data->recenter_transform);
            }, grain_size);

            tasks[num_tasks++] = calc_transform_task;
            tasks[num_tasks++] = apply_transform_task;
        }
    }

    if (state->operations.apply_pbc) {
        task_system::ID pbc_task = task_system::create_pool_task(STR_LIT("## Apply PBC"), (uint32_t)sys.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            size_t count = range_end - range_beg;
            float* x = data->dst_x + range_beg;
            float* y = data->dst_y + range_beg;
            float* z = data->dst_z + range_beg;
            md_util_pbc(x, y, z, NULL, count, &data->unitcell);
        });
        tasks[num_tasks++] = pbc_task;
    } 

    if (state->operations.unwrap_structures) {
        size_t num_structures = md_index_data_num_ranges(&sys.structure);
        task_system::ID unwrap_task = task_system::create_pool_task(STR_LIT("## Unwrap Structures"), (uint32_t)num_structures, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            for (uint32_t i = range_beg; i < range_end; ++i) {
                int*     s_idx = md_index_range_ptr(&data->state->mold.sys.structure, i);
                size_t   s_len = md_index_range_size(&data->state->mold.sys.structure, i);
                md_util_unwrap(data->dst_x, data->dst_y, data->dst_z, s_idx, s_len, &data->state->mold.sys.bond, &data->unitcell);
            }
        });
        tasks[num_tasks++] = unwrap_task;
    }

    {
        // Calculate a global AABB for the molecule
        task_system::ID aabb_task = task_system::create_pool_task(STR_LIT("## Compute AABB"), (uint32_t)sys.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            size_t range_len = range_end - range_beg;
            const float* x = data->state->mold.sys.atom.x + range_beg;
            const float* y = data->state->mold.sys.atom.y + range_beg;
            const float* z = data->state->mold.sys.atom.z + range_beg;

            size_t temp_pos = md_temp_get_pos();
            defer { md_temp_set_pos_back(temp_pos); };
            float* r = (float*)md_temp_push(sizeof(float) * range_len);
            md_atom_extract_radii(r, range_beg, range_len, &data->state->mold.sys.atom);

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
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    const md_backbone_angles_t* src_angle = data->t < 0.5f ? src_angles[0] : src_angles[1];
                    MEMCPY(data->state->mold.sys.protein_backbone.segment.angle, src_angle, data->state->mold.sys.protein_backbone.segment.count * sizeof(md_backbone_angles_t));
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            case InterpolationMode::Linear: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Compute Backbone Angles"), (uint32_t)sys.protein_backbone.segment.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_backbone_angles_t* src_angles[2] = {
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    md_system_t& mol = data->state->mold.sys;
                    for (size_t i = range_beg; i < range_end; ++i) {
                        float phi[2] = {src_angles[0][i].phi, src_angles[1][i].phi};
                        float psi[2] = {src_angles[0][i].psi, src_angles[1][i].psi};

                        phi[1] = deperiodize_orthof(phi[1], phi[0], (float)TWO_PI);
                        psi[1] = deperiodize_orthof(psi[1], psi[0], (float)TWO_PI);

                        float final_phi = lerp(phi[0], phi[1], data->t);
                        float final_psi = lerp(psi[0], psi[1], data->t);
                        mol.protein_backbone.segment.angle[i] = {deperiodize_orthof(final_phi, 0, (float)TWO_PI), deperiodize_orthof(final_psi, 0, (float)TWO_PI)};
                    }
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            case InterpolationMode::CubicSpline: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Interpolate Backbone Angles"), (uint32_t)sys.protein_backbone.segment.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_backbone_angles_t* src_angles[4] = {
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[0],
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[1],
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[2],
                        data->state->trajectory_data.backbone_angles.data + data->state->trajectory_data.backbone_angles.stride * data->frames[3],
                    };
                    md_system_t& sys = data->state->mold.sys;
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
        if (md_array_size(state->interpolated_properties.secondary_structure) != sys.protein_backbone.segment.count) {
			MD_LOG_ERROR("Secondary structure array size does not match the number of segments.");
        }
        size_t num_backbone_segments = sys.protein_backbone.segment.count;
        task_system::ID ss_task = task_system::create_pool_task(STR_LIT("## Interpolate Secondary Structures"), (uint32_t)num_backbone_segments, [data = &payload, mode](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            const md_secondary_structure_t* src_ss[4] = {
                (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[0],
                (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[1],
                (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[2],
                (md_secondary_structure_t*)data->state->trajectory_data.secondary_structure.data + data->state->trajectory_data.secondary_structure.stride * data->frames[3],
            };
            const md_secondary_structure_t* src_ss_nearest = data->t < 0.5f ? src_ss[1] : src_ss[2];
            switch (mode) {
            default:
                MD_LOG_DEBUG("Unsupported interpolation mode for secondary structure interpolation");
                [[fallthrough]];
            case InterpolationMode::Nearest: {
                for (size_t i = range_beg; i < range_end; ++i) {
                    md_secondary_structure_t ss = src_ss_nearest[i];
                    // Set both the analytical (nearest) and interpolated secondary structure (rendering)
                    data->state->mold.sys.protein_backbone.segment.secondary_structure[i] = ss;
                    data->state->interpolated_properties.secondary_structure[i] = md_gl_secondary_structure_convert(ss);
                }
                break;
            }
            case InterpolationMode::Linear: {
                for (size_t i = range_beg; i < range_end; ++i) {
                    md_secondary_structure_t ss[2] = { src_ss[1][i], src_ss[2][i] };
                    md_gl_secondary_structure_t ss_gl[2] = { md_gl_secondary_structure_convert(ss[0]), md_gl_secondary_structure_convert(ss[1]) };
                    md_gl_secondary_structure_t ss_gl_i = {
                        .helix = lerp(ss_gl[0].helix, ss_gl[1].helix, data->t),
                        .sheet = lerp(ss_gl[0].sheet, ss_gl[1].sheet, data->t),
                    };
                    // Set both the analytical (nearest) and interpolated secondary structure (rendering)
                    data->state->mold.sys.protein_backbone.segment.secondary_structure[i] = src_ss_nearest[i];
                    data->state->interpolated_properties.secondary_structure[i] = ss_gl_i;
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
                    
                    
                    // Cleanup isolated coils temporally to reduce noise during transitions.
                    // This is a common issue with secondary structure assignment during transitions, where a segment might flip between coil and helix/sheet rapidly, creating a flickering effect.
                    // We compare indices 0, 1, 2, 3 and if idx 1 or 2 differs from the other 3 (which are the same), we set 1 to be the same as the others, effectively removing isolated coil assignments.
                    
                    #if 0
                    const md_gl_secondary_structure_t ss_coil = { 0,0 };
                    const md_gl_secondary_structure_t ss_sheet = { .sheet = 1.0f };

                    auto is_eq = [](md_gl_secondary_structure_t a, md_gl_secondary_structure_t b) {
                        return a.helix == b.helix && a.sheet == b.sheet;
                    };

                    if (is_eq(ss_gl[1], ss_coil)) {
                        if (is_eq(ss_gl[0], ss_gl[2]) && is_eq(ss_gl[0], ss_gl[3]) && !is_eq(ss_gl[1], ss_gl[0])) {
                            ss_gl[1] = ss_gl[0];
                        }
                    }
                    if (is_eq(ss_gl[2], ss_coil)) {
                        if (is_eq(ss_gl[0], ss_gl[1]) && is_eq(ss_gl[0], ss_gl[3]) && !is_eq(ss_gl[2], ss_gl[0])) {
                            ss_gl[2] = ss_gl[0];
                        }
                    }
#endif

                    md_gl_secondary_structure_t ss_gl_i = {
                        .helix = cubic_spline(ss_gl[0].helix, ss_gl[1].helix, ss_gl[2].helix, ss_gl[3].helix, data->t, data->s),
                        .sheet = cubic_spline(ss_gl[0].sheet, ss_gl[1].sheet, ss_gl[2].sheet, ss_gl[3].sheet, data->t, data->s),
                    };
                    // Set both the analytical (nearest) and interpolated secondary structure (rendering)
                    data->state->mold.sys.protein_backbone.segment.secondary_structure[i] = src_ss_nearest[i];
                    data->state->interpolated_properties.secondary_structure[i] = ss_gl_i;
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
            const md_gl_secondary_structure_t ss_sheet = { .sheet = 1.0f };

            md_gl_secondary_structure_t* ss_gl = data->state->interpolated_properties.secondary_structure;
            for (size_t i = 0; i < data->state->mold.sys.protein_backbone.range.count; ++i) {
                size_t range_beg = data->state->mold.sys.protein_backbone.range.offset[i];
                size_t range_end = data->state->mold.sys.protein_backbone.range.offset[i + 1];
                for (size_t j = range_beg + 1; j + 1 < range_end; ++j) {
                    // Set isolated coils between sheets to sheet to reduce noise during transitions, as this is likely a result of flickering during secondary structure assignment
                    if (is_eq(ss_gl[j - 1], ss_sheet) && is_eq(ss_gl[j + 1], ss_sheet) && is_eq(ss_gl[j], ss_coil)) {
                        ss_gl[j] = ss_sheet;
                    }
                }
            }
        });
        tasks[num_tasks++] = ss_cleanup_task;
#endif

        state->mold.dirty_gpu_buffers |= MolBit_DirtySecondaryStructure;
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
    state->mold.sys_aabb_min = aabb_min;
    state->mold.sys_aabb_max = aabb_max;
    sys.unitcell = payload.unitcell;

    // unitcell transform is essentially just a translation to place the center of the unitcell at the origin
    mat3_t A;
    md_unitcell_A_extract_float(A.elem, &state->mold.sys.unitcell);
    vec3_t c = mat3_mul_vec3(A, vec3_set(0.5f, 0.5f, 0.5f));
    state->mold.unitcell_transform = mat4_translate(-c.x, -c.y, -c.z);

#if 0
    if (mol.unitcell.flags) {
        vec3_t c = mol.unitcell.basis * vec3_set1(0.5f);
        state->mold.model_mat = mat4_translate_vec3(-c);
    }
#endif

    state->mold.dirty_gpu_buffers |= MolBit_DirtyPosition;
}

void recenter_update_target_data(ApplicationState* state) {
    if (!state->mold.sys.trajectory) return;

    uint64_t hash = md_bitfield_hash64(&state->operations.target_mask, 0xDEADBEEF);
    if (hash != state->operations.initial_frame.hash) {
        // Need to recalculate the initial frame
        state->operations.initial_frame.hash = hash;
        size_t count = md_bitfield_popcount(&state->operations.target_mask);

        md_array_resize(state->operations.initial_frame.xyzw, count, state->allocator.persistent);

        state->operations.initial_frame.com = vec3_zero();
        state->operations.initial_frame.alignment_mat = mat4_ident();

        if (state->operations.initial_frame.xyzw && count > 0) {
            // Fetch initial frame data required for orienting the structure
            size_t num_atoms = state->mold.sys.atom.count;
            float* temp_x = (float*)md_vm_arena_push(state->allocator.frame, sizeof(float) * num_atoms);
            float* temp_y = (float*)md_vm_arena_push(state->allocator.frame, sizeof(float) * num_atoms);
            float* temp_z = (float*)md_vm_arena_push(state->allocator.frame, sizeof(float) * num_atoms);

            md_trajectory_load_frame(state->mold.sys.trajectory, 0, NULL, temp_x, temp_y, temp_z);

            md_bitfield_iter_t it = md_bitfield_iter_create(&state->operations.target_mask);
            int dst_idx = 0;
            while (md_bitfield_iter_next(&it)) {
                uint64_t src_idx = md_bitfield_iter_idx(&it);
                float mass = md_atom_mass(&state->mold.sys.atom, src_idx);
                state->operations.initial_frame.xyzw[dst_idx++] = vec4_set(temp_x[src_idx], temp_y[src_idx], temp_z[src_idx], mass);
            }

            state->operations.initial_frame.com = md_util_com_compute_vec4(state->operations.initial_frame.xyzw, NULL, count, &state->mold.sys.unitcell);
        }
    }
}

void recenter_calculate_transform(float M[4][4], const ApplicationState* state) {
    ASSERT(M);
    ASSERT(state);

    size_t count = md_bitfield_popcount(&state->operations.target_mask);
    mat4_t transform = mat4_ident();

    if (count > 0) {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(state->allocator.frame);
        defer { md_vm_arena_temp_end(temp); };

        // Extract xyzw subset of target
        vec4_t* target_xyzw = (vec4_t*)md_vm_arena_push(state->allocator.frame, sizeof(vec4_t) * count);
        md_util_system_extract_xyzw_from_mask(target_xyzw, &state->operations.target_mask, &state->mold.sys);

        // Unwrap target structure (required for rotation)
        md_util_unwrap_vec4(target_xyzw, NULL, count, &state->mold.sys.bond, &state->mold.sys.unitcell);

        // Calculate target
        vec3_t target = {0};
        if (md_unitcell_flags(&state->mold.sys.unitcell) != 0) {
            mat3_t A = {0};
            md_unitcell_A_extract_float(A.elem, &state->mold.sys.unitcell);
            target = mat3_mul_vec3(A, vec3_set1(0.5f));
        } 

        // Calculate COM
        vec3_t target_com = md_util_com_compute_vec4(target_xyzw, NULL, count, &state->mold.sys.unitcell);

        // Calculate Rotation
        mat3_t R = mat3_ident();
        if (state->operations.rotate && state->operations.initial_frame.xyzw) {
            ASSERT(md_array_size(state->operations.initial_frame.xyzw) == count);
            const vec4_t* xyzw[2] = {
                state->operations.initial_frame.xyzw,
                target_xyzw,
            };
            const vec3_t com[2] = {
                state->operations.initial_frame.com,
                target_com,
            };
            R = mat3_optimal_rotation_vec4(xyzw, NULL, count, com);
            R = mat3_orthonormalize(R);
        }
        const mat4_t A = state->operations.initial_frame.alignment_mat;
        mat4_t T = mat4_translate_vec3(target) * A * mat4_from_mat3(R) * mat4_translate_vec3(-target_com);
        transform = T;
    }
    mat4_store((float*)M, transform);
}

bool picking_reserve_range(PickingRange* out_range, PickingSpace* space, PickingDomainID domain, size_t count) {
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
    MEMSET(&handler->history[slot_idx].space, 0, sizeof(PickingSpace));
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

    MEMSET(surface, 0, sizeof(PickingSurface));
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

    MEMSET(surface, 0, sizeof(PickingSurface));
}

bool picking_surface_submit_readback(
    PickingSurface* surface,
    uint32_t fbo,
    uint32_t width,
    uint32_t height,
    uint32_t submitted_frame_idx,
    vec2_t surface_coord,
    vec2_t screen_coord,
    const mat4_t& inv_mvp
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
    slot.inv_mvp = inv_mvp;

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
    if (raw_idx == INVALID_PICKING_IDX) {
        return false;
    }

    const PickingSpace* space = picking_handler_find_space(handler, slot.submitted_frame_idx);
    if (!space) {
        return false;
    }

    const PickingRange* range = find_picking_range(space, raw_idx);
    if (!range) {
        return false;
    }

    out_hit->source = surface->source;
    out_hit->domain = range->domain;
    out_hit->frame_idx = slot.submitted_frame_idx;
    out_hit->raw_idx = raw_idx;
    out_hit->local_idx = raw_idx - range->beg;
    out_hit->surface_coord = slot.surface_coord;
    out_hit->screen_coord = slot.screen_coord;
    out_hit->depth = depth;

    const vec4_t viewport = {0, 0, (float)slot.viewport_width, (float)slot.viewport_height};
    out_hit->world_pos = mat4_unproject({slot.surface_coord.x, slot.surface_coord.y, depth}, slot.inv_mvp, viewport);

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
        request.inv_mvp
    );

    return picking_surface_poll_hit(out_hit, surface, handler);
}

bool interaction_surface(PickingHit* out_hit, const InteractionSurfaceArgs& args) {
    ASSERT(out_hit);
    ASSERT(args.sys);
    ASSERT(args.highlight_mask);
    ASSERT(args.selection_mask);
    ASSERT(args.single_selection_sequence);
    ASSERT(args.camera);
    ASSERT(args.mvp);
    ASSERT(args.inv_mvp);
    ASSERT(args.view_target);
    ASSERT(args.trackball_param);
    ASSERT(args.picking_surface);
    ASSERT(args.picking_handler);

    *out_hit = {};

    enum class RegionMode { Append, Remove };

    const bool pressed = ImGui::InvisibleButton(
        "canvas",
        {args.size.x, args.size.y},
        ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_AllowOverlap
    );

    ImGuiWindow* window = ImGui::GetCurrentWindow();
    if (!window) {
        return pressed;
    }

    bool is_performing_region_selection = false;
    ImDrawList* draw_list = window->DrawList;
    ASSERT(draw_list);

    const ImVec2 canvas_min = ImGui::GetItemRectMin();
    const ImVec2 canvas_max = ImGui::GetItemRectMax();
    const ImVec2 canvas_size = ImGui::GetItemRectSize();

    if (ImGui::IsItemHovered() && args.fbo && args.width && args.height) {
        const ImVec2 mouse_pos = ImGui::GetMousePos();
        const ImVec2 local_mouse = mouse_pos - canvas_min;
        const ImVec2 local_coord = ImVec2(local_mouse.x, canvas_size.y - local_mouse.y) * ImGui::GetIO().DisplayFramebufferScale;

        PickingReadbackRequest request = {
            .fbo = args.fbo,
            .width = args.width,
            .height = args.height,
            .surface_coord = { local_coord.x, local_coord.y },
            .screen_coord = { local_mouse.x, local_mouse.y },
            .inv_mvp = *args.inv_mvp,
        };

        picking_surface_submit_readback_and_poll_hit(out_hit, args.picking_surface, *args.picking_handler, request);
    }

    const bool valid_hit = out_hit->raw_idx != INVALID_PICKING_IDX &&
        (args.expected_source == 0 || out_hit->source == args.expected_source);

    if (pressed || ImGui::IsItemActive() || ImGui::IsItemDeactivated()) {
        if (ImGui::IsKeyPressed(ImGuiMod_Shift, false)) {
            ImGui::ResetMouseDragDelta();
        }

        if (ImGui::IsKeyDown(ImGuiMod_Shift)) {
            RegionMode mode = RegionMode::Append;
            if (ImGui::IsMouseDown(ImGuiMouseButton_Left) || ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                mode = RegionMode::Append;
            } else if (ImGui::IsMouseDown(ImGuiMouseButton_Right) || ImGui::IsMouseReleased(ImGuiMouseButton_Right)) {
                mode = RegionMode::Remove;
            }

            const ImVec2 ext = ImGui::GetMouseDragDelta(mode == RegionMode::Append ? ImGuiMouseButton_Left : ImGuiMouseButton_Right);
            const ImVec2 pos = ImGui::GetMousePos() - ext;
            const ImU32 fill_col = 0x22222222;
            const ImU32 line_col = 0x88888888;

            ImVec2 sel_min = ImClamp(ImMin(pos, pos + ext), canvas_min, canvas_max);
            ImVec2 sel_max = ImClamp(ImMax(pos, pos + ext), canvas_min, canvas_max);

            draw_list->AddRectFilled(sel_min, sel_max, fill_col);
            draw_list->AddRect(sel_min, sel_max, line_col);

            sel_min -= canvas_min;
            sel_max -= canvas_min;

            md_allocator_i* temp_alloc = md_get_temp_allocator();
            md_bitfield_t mask = md_bitfield_create(temp_alloc);

            if (sel_min != sel_max) {
                md_bitfield_clear(args.highlight_mask);
                is_performing_region_selection = true;

                if (args.candidate_mask) {
                    md_bitfield_iter_t it = md_bitfield_iter_create(args.candidate_mask);
                    while (md_bitfield_iter_next(&it)) {
                        const uint64_t idx = md_bitfield_iter_idx(&it);
                        const vec4_t p = mat4_mul_vec4(*args.mvp, vec4_set(args.sys->atom.x[idx], args.sys->atom.y[idx], args.sys->atom.z[idx], 1.0f));
                        const vec2_t c = {
                            ( p.x / p.w * 0.5f + 0.5f) * canvas_size.x,
                            (-p.y / p.w * 0.5f + 0.5f) * canvas_size.y,
                        };

                        if (sel_min.x <= c.x && c.x <= sel_max.x && sel_min.y <= c.y && c.y <= sel_max.y) {
                            md_bitfield_set_bit(&mask, idx);
                        }
                    }
                } else {
                    for (size_t idx = 0; idx < args.sys->atom.count; ++idx) {
                        const vec4_t p = mat4_mul_vec4(*args.mvp, vec4_set(args.sys->atom.x[idx], args.sys->atom.y[idx], args.sys->atom.z[idx], 1.0f));
                        const vec2_t c = {
                            ( p.x / p.w * 0.5f + 0.5f) * canvas_size.x,
                            (-p.y / p.w * 0.5f + 0.5f) * canvas_size.y,
                        };

                        if (sel_min.x <= c.x && c.x <= sel_max.x && sel_min.y <= c.y && c.y <= sel_max.y) {
                            md_bitfield_set_bit(&mask, idx);
                        }
                    }
                }

                grow_mask_by_selection_granularity(&mask, args.selection_granularity, *args.sys);

                if (mode == RegionMode::Append) {
                    md_bitfield_or(args.highlight_mask, args.selection_mask, &mask);
                } else {
                    md_bitfield_andnot(args.highlight_mask, args.selection_mask, &mask);
                }

                if (pressed || ImGui::IsMouseReleased(ImGuiMouseButton_Left)) {
                    md_bitfield_copy(args.selection_mask, args.highlight_mask);
                }
            } else if (pressed) {
                if (valid_hit) {
                    if (mode == RegionMode::Append) {
                        if (out_hit->domain == PickingDomain_Atom) {
                            const int atom_idx = (int)out_hit->local_idx;
                            single_selection_sequence_push_idx(args.single_selection_sequence, atom_idx);
                            md_bitfield_set_bit(&mask, atom_idx);
                        } else if (out_hit->domain == PickingDomain_Bond && args.sys->bond.pairs) {
                            const md_atom_pair_t pair = args.sys->bond.pairs[out_hit->local_idx];
                            md_bitfield_set_bit(&mask, pair.idx[0]);
                            md_bitfield_set_bit(&mask, pair.idx[1]);
                        }
                        grow_mask_by_selection_granularity(&mask, args.selection_granularity, *args.sys);
                        md_bitfield_or_inplace(args.selection_mask, &mask);
                    } else {
                        if (out_hit->domain == PickingDomain_Atom) {
                            const int atom_idx = (int)out_hit->local_idx;
                            single_selection_sequence_pop_idx(args.single_selection_sequence, atom_idx);
                            md_bitfield_set_bit(&mask, atom_idx);
                        } else if (out_hit->domain == PickingDomain_Bond && args.sys->bond.pairs) {
                            const md_atom_pair_t pair = args.sys->bond.pairs[out_hit->local_idx];
                            md_bitfield_set_bit(&mask, pair.idx[0]);
                            md_bitfield_set_bit(&mask, pair.idx[1]);
                        }
                        grow_mask_by_selection_granularity(&mask, args.selection_granularity, *args.sys);
                        md_bitfield_andnot_inplace(args.selection_mask, &mask);
                    }
                } else {
                    single_selection_sequence_clear(args.single_selection_sequence);
                    md_bitfield_clear(args.selection_mask);
                    md_bitfield_clear(args.highlight_mask);
                }
            }
        }
    } else if (ImGui::IsItemHovered() && !ImGui::IsAnyItemActive()) {
        md_bitfield_clear(args.highlight_mask);
        if (valid_hit) {
            if (out_hit->domain == PickingDomain_Atom) {
                const int atom_idx = (int)out_hit->local_idx;
                if (0 <= atom_idx && atom_idx < (int)args.sys->atom.count) {
                    md_bitfield_set_bit(args.highlight_mask, atom_idx);
                    grow_mask_by_selection_granularity(args.highlight_mask, args.selection_granularity, *args.sys);
                }
            } else if (out_hit->domain == PickingDomain_Bond && args.sys->bond.pairs) {
                const int bond_idx = (int)out_hit->local_idx;
                if (0 <= bond_idx && bond_idx < (int)args.sys->bond.count) {
                    const md_atom_pair_t pair = args.sys->bond.pairs[bond_idx];
                    if (0 <= pair.idx[0] && pair.idx[0] < (int)args.sys->atom.count) {
                        md_bitfield_set_bit(args.highlight_mask, pair.idx[0]);
                    }
                    if (0 <= pair.idx[1] && pair.idx[1] < (int)args.sys->atom.count) {
                        md_bitfield_set_bit(args.highlight_mask, pair.idx[1]);
                    }
                    grow_mask_by_selection_granularity(args.highlight_mask, args.selection_granularity, *args.sys);
                }
            }
        }
    }

    if (ImGui::IsItemActive() || ImGui::IsItemHovered()) {
        if (!ImGui::IsKeyDown(ImGuiMod_Shift) && !is_performing_region_selection) {
            const ImVec2 delta = ImGui::GetIO().MouseDelta;
            const ImVec2 coord = ImGui::GetMousePos() - canvas_min;
            const float scroll_delta = ImGui::GetIO().MouseWheel;

            TrackballControllerInput input = {};
            input.rotate_button = ImGui::IsMouseDown(ImGuiMouseButton_Left);
            input.pan_button = ImGui::IsMouseDown(ImGuiMouseButton_Right);
            input.dolly_button = ImGui::IsMouseDown(ImGuiMouseButton_Middle);
            input.mouse_coord_curr = {coord.x, coord.y};
            input.mouse_coord_prev = {coord.x - delta.x, coord.y - delta.y};
            input.screen_size = {canvas_size.x, canvas_size.y};
            input.dolly_delta = scroll_delta;
            input.fov_y = args.camera->fov_y;

            TrackballFlags flags = TrackballFlags_None;
            if (ImGui::IsItemActive()) {
                flags |= TrackballFlags_EnableAllInteractions;
            } else {
                flags |= TrackballFlags_DollyEnabled;
            }

            camera_controller_trackball(args.view_target, input, *args.trackball_param, flags);
        }
    }

    return pressed;
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
            reset_view(state, &state->representation.visibility_mask, false, true);
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
            if (!str_empty(rel_path)) {
                snprintf(buf, sizeof(buf), "table = import(\"%.*s\");\n", STR_ARG(rel_path));
                TextEditor::Coordinates pos = state->editor.GetCursorPosition();
                pos.mLine += 1;
                state->editor.SetCursorPosition({0,0});
                state->editor.InsertText(buf);
                state->editor.SetCursorPosition(pos);
            }
        } else {
            loader::State loader_state = {};
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
                        reset_view(state, &state->representation.visibility_mask, true, false);
                    }
                }
            }
        }
    }
}

void reset_view(ApplicationState* state, const md_bitfield_t* target, bool move_camera, bool smooth_transition) {
    ASSERT(state);

    md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(state->allocator.frame);
    defer { md_vm_arena_temp_end(tmp); };

    if (!state->mold.sys.atom.count) return;
    const auto& mol = state->mold.sys;

    size_t popcount = 0;
    if (target) {
        popcount = md_bitfield_popcount(target);
    }

    int32_t* indices = nullptr;
    if (0 < popcount && popcount < mol.atom.count) {
        indices = (int32_t*)md_vm_arena_push_array(state->allocator.frame, int32_t, popcount);
        size_t len = md_bitfield_iter_extract_indices(indices, popcount, md_bitfield_iter_create(target));
        if (len > popcount || len > mol.atom.count) {
            MD_LOG_DEBUG("Error: Invalid number of indices");
            len = MIN(popcount, mol.atom.count);
        }
    }

    size_t count = popcount ? popcount : mol.atom.count;
    vec3_t com = md_util_com_compute(mol.atom.x, mol.atom.y, mol.atom.z, nullptr, indices, count, &mol.unitcell);
    mat3_t C = mat3_covariance_matrix(mol.atom.x, mol.atom.y, mol.atom.z, nullptr, indices, count, com);
    mat3_eigen_t eigen = mat3_eigen(C);
    mat3_t PCA = mat3_orthonormalize(mat3_extract_rotation(eigen.vectors));
    mat4_t Ri  = mat4_from_mat3(PCA);

    // Compute min and maximum extent along the PCA axes
    vec3_t min_ext = vec3_set1( FLT_MAX);
    vec3_t max_ext = vec3_set1(-FLT_MAX);
    mat3_t basis = mat3_transpose(PCA);

    const float radius = 1.0f;

    // Transform the atom (x,y,z,radius) into the PCA frame to find the min and max extend within it
    for (size_t i = 0; i < count; ++i) {
        int32_t idx = indices ? indices[i] : (int32_t)i;
        vec3_t xyz = { mol.atom.x[idx], mol.atom.y[idx], mol.atom.z[idx] };

        vec3_t p = mat4_mul_vec3(Ri, xyz, 1.0f);
        min_ext = vec3_min(min_ext, vec3_sub_f(p, radius));
        max_ext = vec3_max(max_ext, vec3_add_f(p, radius));
    }

    ViewTransform opt_view = compute_optimal_view(min_ext, max_ext, basis);
	opt_view.position = mat4_mul_vec3(state->mold.unitcell_transform, opt_view.position, 1.0f);

    if (move_camera) {
		state->view.target = opt_view;
        if (!smooth_transition) {
			state->view.camera = opt_view;
        }
    }

    mat3_t A = { 0 };
	md_unitcell_A_extract_float(A.elem, &mol.unitcell);
    const vec3_t cell_ext = mat3_mul_vec3(A, vec3_set1(1.0f));
    const float  max_cell_ext = vec3_reduce_max(cell_ext);
    const float  max_aabb_ext = vec3_reduce_max(vec3_sub(max_ext, min_ext));

    state->view.camera.near_plane = 1.0f;
    state->view.camera.far_plane = 100000.0f;
    state->view.trackball_param.max_distance = MAX(max_cell_ext, max_aabb_ext) * 10.0f;
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
            picking_reserve_range(&state->picking_range_atom, space, PickingDomain_Atom, num_atoms);
            picking_reserve_range(&state->picking_range_bond, space, PickingDomain_Bond, num_bonds);
            break;
        }
        case viamd::EventType_ViamdPickingHit: {
            break;
        }
        case viamd::EventType_ViamdPickingTooltipTextRequest: {
            ASSERT(event.payload_type == viamd::EventPayloadType_PickingTooltipTextRequest);
            PickingTooltipTextRequest* req = (PickingTooltipTextRequest*)event.payload;
            if (req->hit.domain == PickingDomain_Atom || req->hit.domain == PickingDomain_Bond) {
                fill_picking_tooltip_text(&req->sb, *state, req->hit);
            }
            break;
        }
        default:
            break;
        }
    }
}
