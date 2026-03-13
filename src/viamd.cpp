#include <viamd.h>

#include <core/md_log.h>
#include <core/md_str_builder.h>
#include <core/md_arena_allocator.h>
#include <md_util.h>
#include <md_filter.h>

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

const str_t* find_in_arr(str_t str, const str_t arr[], size_t len) {
    for (size_t i = 0; i < len; ++i) {
        if (str_eq(arr[i], str)) {
            return &arr[i];
        }
    }
    return NULL;
}

static void init_all_representations(ApplicationState* state);

Dataset* create_new_dataset(ApplicationState& state) {
    if (state.dataset.num_datasets >= ARRAY_SIZE(state.dataset.datasets)) {
        MD_LOG_ERROR("Maximum number of datasets reached");
        return nullptr;
    }
    Dataset& dataset = state.dataset.datasets[state.dataset.num_datasets++];
    snprintf(dataset.identifier, sizeof(dataset.identifier), "dataset_%i", state.dataset.num_datasets);

    dataset.alloc  = md_arena_allocator_create(state.allocator.persistent, MEGABYTES(1));
    dataset.representation.info.alloc = md_arena_allocator_create(state.allocator.persistent, MEGABYTES(1));

    md_bitfield_init(&current_dataset(state).selection.selection_mask, dataset.alloc);
    md_bitfield_init(&current_dataset(state).selection.highlight_mask, dataset.alloc);
    md_bitfield_init(&state.selection.query.mask, dataset.alloc);
    md_bitfield_init(&state.selection.grow.mask, dataset.alloc);

    md_bitfield_init(&current_dataset(state).representation.visibility_mask, dataset.alloc);

    return &dataset;
}

void file_queue_push(FileQueue* queue, str_t path, FileFlags flags) {
    ASSERT(queue);
    ASSERT(!file_queue_full(queue));
    int prio = 0;

    str_t ext;
    if (extract_ext(&ext, path) && str_eq(ext, WORKSPACE_FILE_EXTENSION)) {
        prio = 1;
    }
    else if (load::mol::loader_from_ext(ext)) {
        prio = 2;
    }
    else if (load::traj::loader_from_ext(ext)) {
        prio = 3;
    }
    else if (find_in_arr(ext, SCRIPT_IMPORT_FILE_EXTENSIONS, ARRAY_SIZE(SCRIPT_IMPORT_FILE_EXTENSIONS))) {
        prio = 4;
    }
    else {
        // Unknown extension
        prio = 5;
        flags |= FileFlags_ShowDialogue;
    }

    uint64_t i = queue->head;
    queue->arr[queue->head] = { str_copy(path, queue->ring), flags, prio };
    queue->head = (queue->head + 1) % ARRAY_SIZE(queue->arr);

    // Sort queue based on prio
    while (i != queue->tail && queue->arr[i].prio < queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)].prio) {
        FileQueue::Entry tmp = queue->arr[i];
        queue->arr[i] = queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)];
        queue->arr[(i - 1) % ARRAY_SIZE(queue->arr)] = tmp;
        i = (i - 1) % ARRAY_SIZE(queue->arr);
    }
}

void draw_info_window(const ApplicationState& state, uint32_t picking_idx) {
    const auto& sys = current_dataset(state).sys;
    if (picking_idx == INVALID_PICKING_IDX) return;
    
	md_vm_arena_temp_t temp = md_vm_arena_temp_begin(state.allocator.frame);
    defer { md_vm_arena_temp_end(temp); };

    md_strb_t sb = md_strb_create(state.allocator.frame);

    if (picking_idx < sys.atom.count) {
        int atom_idx = picking_idx;
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
		if (current_dataset(state).selection.granularity == SelectionGranularity::Atom) {
            md_strb_fmt(&sb, "atom[%i]", atom_idx + 1);
            if (comp_idx != -1) {
                md_strb_fmt(&sb, "[%i]: ", local_idx + 1);
            } else {
				md_strb_push_cstr(&sb, ": ");
            }
            md_strb_push_str(&sb, type);
			md_strb_push_char(&sb, ' ');
            if (z) {
                md_strb_fmt(&sb, "%.*s %.*s ", STR_ARG(elem), STR_ARG(symb));
            }
            md_strb_fmt(&sb, "(%.3f, %.3f, %.3f)\n", pos.x, pos.y, pos.z);    
        }
        
        if (comp_idx != -1 && (current_dataset(state).selection.granularity == SelectionGranularity::Atom || current_dataset(state).selection.granularity == SelectionGranularity::Component)) {
            md_strb_fmt(&sb, "comp[%i]", comp_idx + 1);
            if (comp_name) {
                md_strb_fmt(&sb, ": " STR_FMT, STR_ARG(comp_name));
            }
			md_strb_fmt(&sb, " (seq_id: %i)\n", comp_seq_id);
        }
        if (inst_idx != -1) {
            md_strb_fmt(&sb, "inst[%i]", inst_idx + 1);
            if (inst_id) {
                md_strb_fmt(&sb, ": " STR_FMT, STR_ARG(inst_id));
            }
			if (auth_id) {
                md_strb_fmt(&sb, " (" STR_FMT ")", STR_ARG(auth_id));
            }
            md_strb_push_char(&sb, '\n');
        }

        uint32_t flags = 0;

        if (current_dataset(state).selection.granularity == SelectionGranularity::Atom) {
			flags = md_system_atom_flags(&sys, atom_idx);
		} else if (current_dataset(state).selection.granularity == SelectionGranularity::Component && comp_idx != -1) {
			flags = md_system_component_flags(&sys, comp_idx);
        } else if (current_dataset(state).selection.granularity == SelectionGranularity::Instance && inst_idx != -1) {
			flags = md_system_instance_flags(&sys, inst_idx);
        }

		const uint32_t TERM_N = MD_FLAG_AMINO_ACID   | MD_FLAG_TERMINAL_BEG;
		const uint32_t TERM_C = MD_FLAG_AMINO_ACID   | MD_FLAG_TERMINAL_END;
		const uint32_t TERM_5 = MD_FLAG_NUCLEIC_ACID | MD_FLAG_TERMINAL_BEG;
		const uint32_t TERM_3 = MD_FLAG_NUCLEIC_ACID | MD_FLAG_TERMINAL_END;

        if (flags) {
            sb += "flags: ";
            if (flags & MD_FLAG_HETERO)         { sb += "HETERO "; }
			if (flags & MD_FLAG_POLYPEPTIDE)    { sb += "POLYPEPTIDE "; }
            if (flags & MD_FLAG_AMINO_ACID)     { sb += "AMINO-ACID "; }
            if (flags & MD_FLAG_SIDE_CHAIN)     { sb += "SIDE-CHAIN "; }
			if (flags & MD_FLAG_NUCLEIC_ACID)   { sb += "NUCLEIC-ACID "; }
            if (flags & MD_FLAG_NUCLEOTIDE)     { sb += "NUCLEOTIDE "; }
            if (flags & MD_FLAG_NUCLEOSIDE)     { sb += "NUCLEOSIDE "; }
            if (flags & MD_FLAG_NUCLEOBASE)     { sb += "NUCLEOBASE "; }
            if (flags & MD_FLAG_WATER)          { sb += "WATER "; }
            if (flags & MD_FLAG_ION)            { sb += "ION "; }
            if (flags & MD_FLAG_BACKBONE)       { sb += "BACKBONE "; }
            if ((flags & TERM_N) == TERM_N)     { sb += "N-TERMINUS "; }
            if ((flags & TERM_C) == TERM_C)     { sb += "C-TERMINUS "; }
			if ((flags & TERM_5) == TERM_5)     { sb += "5'-TERMINUS "; }
			if ((flags & TERM_3) == TERM_3)     { sb += "3'-TERMINUS "; }
            if (flags & MD_FLAG_SP)             { sb += "SP "; }
            if (flags & MD_FLAG_SP2)            { sb += "SP2 "; }
            if (flags & MD_FLAG_SP3)            { sb += "SP3 "; }
            if (flags & MD_FLAG_AROMATIC)       { sb += "AROMATIC "; }
            sb += "\n";
        }
        /*
        // @TODO: REIMPLEMENT THIS
        if (res_idx < sys.backbone.segment.angleangles.size() && res_idx < sys.backbone.segments.size() && valid_backbone_atoms(sys.backbone.segments[res_idx])) {
        const auto angles = RAD_TO_DEG((vec2)sys.backbone.angles[res_idx]);
        len += snprintf(buff + len, 256 - len, u8"\u03C6: %.1f\u00b0, \u03C8: %.1f\u00b0\n", angles.x, angles.y);
        }
        */
    }
    else if (picking_idx >= 0x80000000) {
        int bond_idx = picking_idx & 0x7FFFFFFF;
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

            md_strb_fmt(&sb, "bond: " STR_FMT "%c" STR_FMT "\n", STR_ARG(type0), bond_type, STR_ARG(type1));
            md_strb_fmt(&sb, "flags: %.*s\n", len, bond_flags_buf);
            md_strb_fmt(&sb, "length: %.3f\n", d);
        }
    }

    const ImVec2 offset = { 10.f, 18.f };
    const ImVec2 new_pos = {ImGui::GetMousePos().x + offset.x, ImGui::GetMousePos().y + offset.y};
    ImGui::SetNextWindowPos(new_pos);
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0, 0, 0, 0.5f));
    ImGui::Begin("##Atom Info", 0,
        ImGuiWindowFlags_Tooltip | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoDocking);
    ImGui::Text("%s", md_strb_to_cstr(sb));
    ImGui::End();
    ImGui::PopStyleColor();
}

void extract_picking_data(PickingData& out_picking, GBuffer& gbuffer, const vec2_t& coord, const mat4_t& inv_MVP) {
    out_picking = {};

#if MD_PLATFORM_OSX
    coord = coord * vec_cast(ImGui::GetIO().DisplayFramebufferScale);
#endif
    if (0.f < coord.x && coord.x < (float)gbuffer.width && 0.f < coord.y && coord.y < (float)gbuffer.height) {
        extract_gbuffer_picking_idx_and_depth(&out_picking.idx, &out_picking.depth, &gbuffer, (int)coord.x, (int)coord.y);
        const vec4_t viewport = {0, 0, (float)gbuffer.width, (float)gbuffer.height};
        out_picking.world_coord = mat4_unproject({coord.x, coord.y, out_picking.depth}, inv_MVP, viewport);
        out_picking.screen_coord = {coord.x, coord.y};
    }
}

md_atom_idx_t atom_idx_from_picking_idx(uint32_t picking_idx) {
    if (picking_idx < 0x80000000) {
        return (md_atom_idx_t)picking_idx;
    }
    return -1;
}

md_bond_idx_t bond_idx_from_picking_idx(uint32_t picking_idx) {
    if (picking_idx >= 0x80000000) {
        return (md_bond_idx_t)(picking_idx & 0x7FFFFFFF);
    }
    return -1;
}

void interrupt_async_tasks(ApplicationState* state) {
    task_system::pool_interrupt_running_tasks();

    if (state->script.full_eval) md_script_eval_interrupt(state->script.full_eval);
    if (state->script.filt_eval) md_script_eval_interrupt(state->script.filt_eval);

    task_system::pool_wait_for_completion();
}

// #trajectorydata
void free_trajectory_data(ApplicationState* state) {
    ASSERT(state);
    interrupt_async_tasks(state);

    Dataset& dataset = current_dataset(*state);
    if (dataset.traj_loader) {
        dataset.traj_loader->destroy(dataset.traj);
    }
    dataset.traj = nullptr;
    dataset.traj_loader = nullptr;
    dataset.traj_path = {};
    dataset_clear_recenter_target(dataset);
    dataset_clear_frame_cache(dataset);

    dataset.sys.unitcell = {};
    md_array_free(state->timeline.x_values,  state->allocator.persistent);
    md_array_free(state->display_properties, state->allocator.persistent);

    md_array_free(dataset.trajectory_data.backbone_angles.data,     state->allocator.persistent);
    md_array_free(dataset.trajectory_data.secondary_structure.data, state->allocator.persistent);
}

void init_trajectory_data(ApplicationState* data) {
    size_t num_frames = md_trajectory_num_frames(current_dataset(*data).traj);
    if (num_frames > 0) {
        size_t min_frame = 0;
        size_t max_frame = num_frames - 1;
        md_trajectory_header_t header;
        md_trajectory_get_header(current_dataset(*data).traj, &header);

        double min_time = header.frame_times[0];
        double max_time = header.frame_times[num_frames - 1];

        data->timeline.view_range = {min_time, max_time};
        data->timeline.filter.beg_frame = min_frame;
        data->timeline.filter.end_frame = max_frame;

        md_array_resize(data->timeline.x_values, num_frames, data->allocator.persistent);
        for (size_t i = 0; i < num_frames; ++i) {
            data->timeline.x_values[i] = header.frame_times[i];
        }

        current_dataset(*data).animation.frame = CLAMP(current_dataset(*data).animation.frame, (double)min_frame, (double)max_frame);
        int64_t frame_idx = CLAMP((int64_t)(current_dataset(*data).animation.frame + 0.5), 0, (int64_t)max_frame);

        md_trajectory_frame_header_t frame_header;
        md_trajectory_load_frame(current_dataset(*data).traj, frame_idx, &frame_header, current_dataset(*data).sys.atom.x, current_dataset(*data).sys.atom.y, current_dataset(*data).sys.atom.z);
        current_dataset(*data).sys.unitcell = frame_header.unitcell;

        if (current_dataset(*data).sys.protein_backbone.segment.count > 0) {
            current_dataset(*data).trajectory_data.secondary_structure.stride = current_dataset(*data).sys.protein_backbone.segment.count;
            current_dataset(*data).trajectory_data.secondary_structure.count = current_dataset(*data).sys.protein_backbone.segment.count * num_frames;
            md_array_resize(current_dataset(*data).trajectory_data.secondary_structure.data, current_dataset(*data).sys.protein_backbone.segment.count * num_frames, data->allocator.persistent);
            MEMSET(current_dataset(*data).trajectory_data.secondary_structure.data, 0, md_array_bytes(current_dataset(*data).trajectory_data.secondary_structure.data));

            current_dataset(*data).trajectory_data.backbone_angles.stride = current_dataset(*data).sys.protein_backbone.segment.count;
            current_dataset(*data).trajectory_data.backbone_angles.count = current_dataset(*data).sys.protein_backbone.segment.count * num_frames;
            md_array_resize(current_dataset(*data).trajectory_data.backbone_angles.data, current_dataset(*data).sys.protein_backbone.segment.count * num_frames, data->allocator.persistent);
            MEMSET(current_dataset(*data).trajectory_data.backbone_angles.data, 0, md_array_bytes(current_dataset(*data).trajectory_data.backbone_angles.data));

            // Launch work to compute the values
            task_system::task_interrupt_and_wait_for(data->tasks.backbone_computations);

            data->tasks.backbone_computations = task_system::create_pool_task(STR_LIT("Backbone Operations"), (uint32_t)num_frames, [data](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                // Create copy here of molecule since we use the full structure as input
                md_system_t sys = current_dataset(*data).sys;
                md_trajectory_i* traj = current_dataset(*data).traj;

                md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(1));
                defer { md_vm_arena_destroy(temp_arena); };

                const size_t capacity = ALIGN_TO(sys.atom.count, 16);
                float* x = (float*)md_vm_arena_push(temp_arena, sizeof(float) * capacity);
                float* y = (float*)md_vm_arena_push(temp_arena, sizeof(float) * capacity);
                float* z = (float*)md_vm_arena_push(temp_arena, sizeof(float) * capacity);

                for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                    md_trajectory_frame_header_t frame_header;
                    md_backbone_angles_t* bb_dst = current_dataset(*data).trajectory_data.backbone_angles.data + current_dataset(*data).trajectory_data.backbone_angles.stride * frame_idx;
                    md_secondary_structure_t* ss_dst = current_dataset(*data).trajectory_data.secondary_structure.data + current_dataset(*data).trajectory_data.secondary_structure.stride * frame_idx;

                    md_trajectory_load_frame(traj, frame_idx, &frame_header, x, y, z);
                    md_util_backbone_angles_compute(bb_dst, current_dataset(*data).trajectory_data.backbone_angles.stride, x, y, z, &frame_header.unitcell, &sys.protein_backbone);
                    md_util_backbone_secondary_structure_infer(ss_dst, current_dataset(*data).trajectory_data.secondary_structure.stride, x, y, z, &frame_header.unitcell, &sys.protein_backbone);
                }
                });

            uint64_t time = (uint64_t)md_time_current();
            task_system::ID main_task = task_system::create_main_task(STR_LIT("Update Trajectory Data"), [data, t0 = time]() {
                uint64_t t1 = (uint64_t)md_time_current();
                double elapsed = md_time_as_seconds(t1 - t0);
                MD_LOG_INFO("Finished computing trajectory data (%.2fs)", elapsed);
                current_dataset(*data).trajectory_data.backbone_angles.fingerprint     = generate_fingerprint();
                current_dataset(*data).trajectory_data.secondary_structure.fingerprint = generate_fingerprint();

				current_dataset(*data).interpolate_system_state = true;
                current_dataset(*data).dirty_gpu_buffers |= MolBit_ClearVelocity;
                flag_all_representations_as_dirty(data);
            });

            task_system::set_task_dependency(main_task, data->tasks.backbone_computations);
            task_system::enqueue_task(data->tasks.backbone_computations);
        }

        current_dataset(*data).dirty_gpu_buffers |= MolBit_DirtyPosition;
        current_dataset(*data).dirty_gpu_buffers |= MolBit_ClearVelocity;

        // Prefetch frames
        //launch_prefetch_job(data);
    }
}

bool load_trajectory_data(ApplicationState* data, str_t filename, md_trajectory_loader_i* loader, md_trajectory_flags_t flags) {
    Dataset& dataset = current_dataset(*data);
    md_trajectory_i* traj = loader->create(filename, dataset.alloc, flags);
    if (traj) {
        free_trajectory_data(data);
        dataset.traj = traj;
        dataset.traj_path = str_copy(filename, data->allocator.persistent);
        init_trajectory_data(data);
        dataset.animation.frame = 0;
        return true;
    }

    return false;
}

void init_molecule_data(ApplicationState* data) {
    if (current_dataset(*data).sys.atom.count) {

        data->picking.idx = INVALID_PICKING_IDX;
        current_dataset(*data).selection.atom_idx.hovered = -1;
        current_dataset(*data).selection.atom_idx.right_click = -1;
        current_dataset(*data).selection.bond_idx.hovered = -1;
        current_dataset(*data).selection.bond_idx.right_click = -1;

        current_dataset(*data).gl_mol = md_gl_mol_create(&current_dataset(*data).sys);
        if (current_dataset(*data).sys.protein_backbone.segment.count > 0) {
            current_dataset(*data).interpolated_properties.secondary_structure = md_array_create(md_gl_secondary_structure_t, current_dataset(*data).sys.protein_backbone.segment.count, current_dataset(*data).alloc);
        }

#if EXPERIMENTAL_GFX_API
        const md_system_t& mol = current_dataset(*data).sys;
        vec3_t& aabb_min = current_dataset(*data).sys_aabb_min;
        vec3_t& aabb_max = current_dataset(*data).sys_aabb_max;
        md_util_compute_aabb_soa(&aabb_min, &aabb_max, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.radius, mol.atom.count);

        current_dataset(*data).gfx_structure = md_gfx_structure_create(mol.atom.count, mol.covalent.count, mol.backbone.count, mol.backbone.range_count, mol.residue.count, mol.instance.count);
        md_gfx_structure_set_atom_position(current_dataset(*data).gfx_structure, 0, mol.atom.count, mol.atom.x, mol.atom.y, mol.atom.z, 0);
        md_gfx_structure_set_atom_radius(current_dataset(*data).gfx_structure, 0, mol.atom.count, mol.atom.radius, 0);
        md_gfx_structure_set_aabb(current_dataset(*data).gfx_structure, &current_dataset(*data).sys_aabb_min, &current_dataset(*data).sys_aabb_max);
        if (mol.instance.count > 0) {
            md_gfx_structure_set_instance_atom_ranges(current_dataset(*data).gfx_structure, 0, mol.instance.count, (md_gfx_range_t*)mol.instance.atom_range, 0);
            md_gfx_structure_set_instance_transforms(current_dataset(*data).gfx_structure, 0, mol.instance.count, mol.instance.transform, 0);
        }
#endif
    }
    viamd::event_system_broadcast_event(viamd::EventType_ViamdTopologyInit, viamd::EventPayloadType_ApplicationState, data);

    update_representation_info(data);
    init_all_representations(data);
    data->script.compile_ir = true;

}

static inline void clear_density_volume(ApplicationState* state) {
    md_array_shrink(state->density_volume.gl_reps, 0);
    md_array_shrink(state->density_volume.rep_model_mats, 0);
    state->density_volume.model_mat = {0};
}

void free_molecule_data(ApplicationState* data) {
    ASSERT(data);
    interrupt_async_tasks(data);
    Dataset& dataset = current_dataset(*data);

    md_arena_allocator_reset(dataset.alloc);
    MEMSET(&dataset.sys, 0, sizeof(dataset.sys));

    md_gl_mol_destroy(dataset.gl_mol);
    dataset.sys_path = {};

    dataset.interpolated_properties.secondary_structure = nullptr;

    md_bitfield_clear(&dataset.selection.selection_mask);
    md_bitfield_clear(&dataset.selection.highlight_mask);
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

    viamd::event_system_broadcast_event(viamd::EventType_ViamdTopologyFree, viamd::EventPayloadType_ApplicationState, data);
}

void process_file_queue(ApplicationState* state) {
    if (!file_queue_empty(&state->file_queue) && !state->load_dataset.show_window) {
        FileQueue::Entry e = file_queue_pop(&state->file_queue);

        str_t ext;
        extract_ext(&ext, e.path);
        const str_t* res = 0;

        if (str_eq_ignore_case(ext, WORKSPACE_FILE_EXTENSION)) {
            load_workspace(state, e.path);
            reset_view(state, &current_dataset(*state).representation.visibility_mask, false, true);
        } else if ((res = find_in_arr(ext, SCRIPT_IMPORT_FILE_EXTENSIONS, ARRAY_SIZE(SCRIPT_IMPORT_FILE_EXTENSIONS)))) {
            char buf[1024];
            str_t base_path = {};
            if (!str_empty(state->workspace_path)) {
                base_path = state->workspace_path;
            } else if (!str_empty(current_dataset(*state).traj_path)) {
                base_path = current_dataset(*state).traj_path;
            } else if (!str_empty(current_dataset(*state).sys_path)) {
                base_path = current_dataset(*state).sys_path;
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
            load::LoaderState loader = {};
            bool success = load::init_loader_state(&loader, e.path, state->allocator.frame);

            if (!success || (e.flags & FileFlags_ShowDialogue) || (loader.flags & LoaderStateFlag_RequiresDialogue)) {
                state->load_dataset = LoadDatasetWindowState();
                str_copy_to_char_buf(state->load_dataset.path_buf, sizeof(state->load_dataset.path_buf), e.path);
                state->load_dataset.path_changed = true;
                state->load_dataset.show_window = true;
                state->load_dataset.coarse_grained = e.flags & FileFlags_CoarseGrained;
            } else if (success) {
                LoadParam param = {};
                param.sys_loader  = loader.sys_loader;
                param.traj_loader = loader.traj_loader;
                param.file_path      = e.path;
                param.coarse_grained = e.flags & FileFlags_CoarseGrained;
                param.sys_loader_arg = loader.sys_loader_arg;
                param.traj_loader_flags = (e.flags & FileFlags_DisableCacheWrite) ? MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE : MD_TRAJECTORY_FLAG_NONE;
                if (load_dataset_from_file(state, param)) {
                    current_dataset(*state).animation = {};
                    if (param.sys_loader) {
                        md_bitfield_reset(&current_dataset(*state).representation.visibility_mask);

                        remove_all_representations(state);
                        create_default_representations(state);

                        recompute_atom_visibility_mask(state);
						current_dataset(*state).interpolate_system_state = true;
                        current_dataset(*state).dirty_gpu_buffers |= MolBit_ClearVelocity;
                        reset_view(state, &current_dataset(*state).representation.visibility_mask, true, false);
                    }
                }
            }
        }
    }
}

bool load_dataset_from_file(ApplicationState* data, const LoadParam& param) {
    ASSERT(data);
    Dataset& dataset = current_dataset(*data);

    str_t path_to_file = md_path_make_canonical(param.file_path, data->allocator.frame);
    if (path_to_file) {
        if (param.sys_loader) {
            interrupt_async_tasks(data);
            free_trajectory_data(data);
            free_molecule_data(data);

            if (!param.sys_loader->init_from_file(&dataset.sys, path_to_file, param.sys_loader_arg, dataset.alloc)) {
                VIAMD_LOG_ERROR("Failed to load molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));
                return false;
            }
            VIAMD_LOG_SUCCESS("Successfully loaded molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));

            dataset.sys_path = str_copy(path_to_file, data->allocator.persistent);
            // @NOTE: If the dataset is coarse-grained, then postprocessing must be aware
            md_postprocess_flags_t flags = param.coarse_grained ? MD_UTIL_POSTPROCESS_NONE : MD_UTIL_POSTPROCESS_ALL;
            md_util_system_postprocess(&dataset.sys, dataset.alloc, flags);
            init_molecule_data(data);

            // @NOTE: Some files contain both atomic coordinates and trajectory
            if (param.traj_loader) {
                VIAMD_LOG_INFO("File may also contain trajectory, attempting to load trajectory");
            } else {
                return true;
            }
        }

        if (param.traj_loader) {
            if (!dataset.sys.atom.count) {
                VIAMD_LOG_ERROR("Before loading a trajectory, molecular data needs to be present");
                return false;
            }
            interrupt_async_tasks(data);

            bool success = load_trajectory_data(data, path_to_file, param.traj_loader, param.traj_loader_flags);
            if (success) {
                VIAMD_LOG_SUCCESS("Successfully opened trajectory from file '" STR_FMT "'", STR_ARG(path_to_file));
                return true;
            } else {
                if (param.sys_loader && param.traj_loader) {
                    // Don't record this as an error, as the trajectory may be optional (In case of PDB for example)
                    return true;
                }
                VIAMD_LOG_ERROR("Failed to opened trajectory from file '" STR_FMT "'", STR_ARG(path_to_file));
            }
        }
    }

    return false;
}

void reset_state(ApplicationState* state) {
    interrupt_async_tasks(state);
    free_trajectory_data(state);
    free_molecule_data(state);
    state->workspace_path = {};

    // Reset and clear things
    remove_all_selections(state);
    remove_all_representations(state);
    state->editor.SetText("");
    state->workspace_path = {};
}

void load_workspace(ApplicationState* app_state, str_t filename) {
    ScopedTemp temp_reset;
    md_allocator_i* temp_arena = app_state->allocator.frame;

    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(temp_arena);
    defer { md_vm_arena_temp_end(temp); };

    str_t txt = load_textfile(filename, temp_arena);

    if (str_empty(txt)) {
        VIAMD_LOG_ERROR("Could not open workspace file: '" STR_FMT "'", STR_ARG(filename));
        return;
    }

    reset_state(app_state);

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
                        md_strb_t path = md_strb_create(temp_arena);
                        path += folder;
                        path += file;
                        new_molecule_file = md_path_make_canonical(path, temp_arena);
                    }
                } else if (str_eq(ident, STR_LIT("TrajectoryFile"))) {
                    str_t file;
                    viamd::extract_str(file, arg);
                    if (!str_empty(file)) {
                        md_strb_t path = md_strb_create(temp_arena);
                        path += folder;
                        path += file;
                        new_trajectory_file = md_path_make_canonical(path, temp_arena);
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
                    viamd::extract_flt(current_dataset(*data).animation.fps, arg);
                } else if (str_eq(ident, STR_LIT("Interpolation"))) {
                    int mode;
                    viamd::extract_int(mode, arg);
                    current_dataset(*data).animation.interpolation = (InterpolationMode)mode;
                }
            }
        } else if (str_eq(section, STR_LIT("RenderSettings"))) {
            str_t ident, arg;
            while (viamd::next_entry(ident, arg, state)) {
                if (str_eq(ident, STR_LIT("SsaoEnabled"))) {
                    viamd::extract_bool(current_dataset(*data).visuals.ssao.enabled, arg);
                } else if (str_eq(ident, STR_LIT("SsaoIntensity"))) {
                    viamd::extract_flt(current_dataset(*data).visuals.ssao.intensity, arg);
                } else if (str_eq(ident, STR_LIT("SsaoRadius"))) {
                    viamd::extract_flt(current_dataset(*data).visuals.ssao.radius, arg);
                } else if (str_eq(ident, STR_LIT("SsaoBias"))) {
                    viamd::extract_flt(current_dataset(*data).visuals.ssao.bias, arg);
                } else if (str_eq(ident, STR_LIT("DofEnabled"))) {
                    viamd::extract_bool(current_dataset(*data).visuals.dof.enabled, arg);
                } else if (str_eq(ident, STR_LIT("DofFocusScale"))) {
                    viamd::extract_flt(current_dataset(*data).visuals.dof.focus_scale, arg);
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
                    viamd::extract_flt(data->view.camera.focus_distance, arg);
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
                } else if (str_eq(ident, STR_LIT("StaticColor"))) {
                    viamd::extract_flt_vec(rep->uniform_color.elem, 4, arg);
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

    data->view.animation.target_position    = data->view.camera.position;
    data->view.animation.target_orientation = data->view.camera.orientation;
    data->view.animation.target_distance    = data->view.camera.focus_distance;
    
    str_free(data->workspace_path, data->allocator.persistent);
    data->workspace_path = str_copy(filename, data->allocator.persistent);

    load::LoaderState loader_state = {};
    load::init_loader_state(&loader_state, new_molecule_file, temp_arena);

    LoadParam param = {};
    param.file_path = new_molecule_file;
    param.sys_loader = loader_state.sys_loader;
    param.traj_loader = loader_state.traj_loader;
    param.coarse_grained = new_coarse_grained;
    param.sys_loader_arg = loader_state.sys_loader_arg;

    if (new_molecule_file && load_dataset_from_file(data, param)) {
        current_dataset(*data).sys_path = str_copy(new_molecule_file, data->allocator.persistent);
    } else {
        current_dataset(*data).sys_path = {};
    }

    if (new_trajectory_file) {
        load::init_loader_state(&loader_state, new_trajectory_file, data->allocator.frame);
        param.sys_loader = 0;
        param.traj_loader = loader_state.traj_loader;
        param.file_path = new_trajectory_file;
        if (load_dataset_from_file(data, param)) {
            current_dataset(*data).traj_path = str_copy(new_trajectory_file, data->allocator.persistent);
        }
        current_dataset(*data).animation.frame = new_frame;
    } else {
        current_dataset(*data).traj_path = {};
    }

    //apply_atom_elem_mappings(data);
}

void save_workspace(ApplicationState* app_state, str_t filename) {
    md_file_o* file = md_file_open(filename, MD_FILE_WRITE);
    if (!file) {
        VIAMD_LOG_ERROR("Could not open workspace file for writing: '%.*s", (int)filename.len, filename.ptr);
        return;
    }
    defer { md_file_close(file); };

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

    const Dataset& current = current_dataset(*app_state);
    str_t mol_file  = md_path_make_relative(filename, current.sys_path, temp_alloc);
    str_t traj_file = md_path_make_relative(filename, current.traj_path, temp_alloc);

    viamd::write_section_header(state, STR_LIT("Files"));
    viamd::write_str(state, STR_LIT("MoleculeFile"), mol_file);
    viamd::write_str(state, STR_LIT("TrajectoryFile"), traj_file);

    viamd::write_section_header(state, STR_LIT("Animation"));
    viamd::write_dbl(state, STR_LIT("Frame"), current.animation.frame);
    viamd::write_flt(state, STR_LIT("Fps"), current.animation.fps);
    viamd::write_int(state, STR_LIT("Interpolation"), (int)current.animation.interpolation);

    viamd::write_section_header(state, STR_LIT("RenderSettings"));
    viamd::write_bool(state, STR_LIT("SsaoEnabled"), current.visuals.ssao.enabled);
    viamd::write_flt(state, STR_LIT("SsaoIntensity"), current.visuals.ssao.intensity);
    viamd::write_flt(state, STR_LIT("SsaoRadius"), current.visuals.ssao.radius);
    viamd::write_bool(state, STR_LIT("DofEnabled"), current.visuals.dof.enabled);
    viamd::write_flt(state, STR_LIT("DofFocusScale"), current.visuals.dof.focus_scale);

    viamd::write_section_header(state, STR_LIT("Camera"));
    viamd::write_vec3(state, STR_LIT("Position"), current.view.camera.position);
    viamd::write_quat(state, STR_LIT("Orientation"), current.view.camera.orientation);
    viamd::write_flt(state,  STR_LIT("Distance"), current.view.camera.focus_distance);
    viamd::write_int(state,  STR_LIT("Mode"), (int)current.view.mode);

    for (size_t i = 0; i < md_array_size(current.representation.reps); ++i) {
        const Representation& rep = current.representation.reps[i];
        viamd::write_section_header(state, STR_LIT("Representation"));
        viamd::write_str(state,  STR_LIT("Name"), str_from_cstr(rep.name));
        viamd::write_str(state,  STR_LIT("Filter"), str_from_cstr(rep.filt));
        viamd::write_bool(state, STR_LIT("Enabled"), rep.enabled);
        viamd::write_int(state,  STR_LIT("Type"), (int)rep.type);
        viamd::write_int(state,  STR_LIT("ColorMapping"), (int)rep.color_mapping);
        viamd::write_vec4(state, STR_LIT("StaticColor"), rep.uniform_color);
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

    for (size_t i = 0; i < md_array_size(current.selection.stored_selections); ++i) {
        const Selection& sel = current.selection.stored_selections[i];
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
    md_array_push(current_dataset(*state).selection.stored_selections, sel, state->allocator.persistent);
    return md_array_last(current_dataset(*state).selection.stored_selections);
}

void remove_all_selections(ApplicationState* state) {
    md_array_shrink(current_dataset(*state).selection.stored_selections, 0);
}

void remove_selection(ApplicationState* state, int idx) {
    ASSERT(state);
    if (idx < 0 || (int)md_array_size(current_dataset(*state).selection.stored_selections) <= idx) {
        VIAMD_LOG_ERROR("Index [%i] out of range when trying to remove selection", idx);
    }
    auto item = &current_dataset(*state).selection.stored_selections[idx];
    md_bitfield_free(&item->atom_mask);

    current_dataset(*state).selection.stored_selections[idx] = *md_array_last(current_dataset(*state).selection.stored_selections);
    md_array_pop(current_dataset(*state).selection.stored_selections);
}

// --- Representation ---

static void init_representation(ApplicationState* state, Representation* rep) {
#if EXPERIMENTAL_GFX_API
    rep->gfx_rep = md_gfx_rep_create(current_dataset(*state).sys.atom.count);
#endif
    rep->md_rep = md_gl_rep_create(current_dataset(*state).gl_mol);
    md_bitfield_init(&rep->atom_mask, state->allocator.persistent);

    size_t num_props = md_array_size(current_dataset(*state).representation.info.atom_properties);
    if (num_props > 0) {
        rep->prop.idx = 0;
        rep->prop.range_beg = current_dataset(*state).representation.info.atom_properties[0].value_min;
        rep->prop.range_end = current_dataset(*state).representation.info.atom_properties[0].value_max;
    }

    flag_representation_as_dirty(rep);
}

Representation* create_representation(ApplicationState* state, RepresentationType type, ColorMapping color_mapping, str_t filter) {
    ASSERT(state);
    md_array_push(current_dataset(*state).representation.reps, Representation(), state->allocator.persistent);
    Representation* rep = md_array_last(current_dataset(*state).representation.reps);
    rep->type = type;
    rep->color_mapping = color_mapping;
    if (!str_empty(filter)) {
        str_copy_to_char_buf(rep->filt, sizeof(rep->filt), filter);
    }
    rep->electronic_structure.mo_idx = (int)current_dataset(*state).representation.info.alpha.homo_idx;
    init_representation(state, rep);
    return rep;
}

Representation* clone_representation(ApplicationState* state, const Representation& rep) {
    ASSERT(state);
    md_array_push(current_dataset(*state).representation.reps, rep, state->allocator.persistent);
    Representation* clone = md_array_last(current_dataset(*state).representation.reps);
    clone->md_rep = {0};
    clone->atom_mask = {0};
    init_representation(state, clone);
    return clone;
}

void remove_representation(ApplicationState* state, int idx) {
    ASSERT(state);
    ASSERT(idx < md_array_size(current_dataset(*state).representation.reps));
    auto& rep = current_dataset(*state).representation.reps[idx];
    md_bitfield_free(&rep.atom_mask);
    md_gl_rep_destroy(rep.md_rep);
    if (rep.electronic_structure.vol.tex_id != 0) gl::free_texture(&rep.electronic_structure.vol.tex_id);
    if (rep.electronic_structure.dvr.tf_tex != 0) gl::free_texture(&rep.electronic_structure.dvr.tf_tex);
    md_array_swap_back_and_pop(current_dataset(*state).representation.reps, idx);
    recompute_atom_visibility_mask(state);
}

void recompute_atom_visibility_mask(ApplicationState* state) {
    ASSERT(state);
    auto& mask = current_dataset(*state).representation.visibility_mask;

    md_bitfield_clear(&mask);
    for (size_t i = 0; i < md_array_size(current_dataset(*state).representation.reps); ++i) {
        auto& rep = current_dataset(*state).representation.reps[i];
        if (!rep.enabled) continue;
        md_bitfield_or_inplace(&mask, &rep.atom_mask);
    }
    current_dataset(*state).representation.visibility_mask_hash = md_bitfield_hash64(&mask, 0);
}

void update_all_representations(ApplicationState* state) {
    for (size_t i = 0; i < md_array_size(current_dataset(*state).representation.reps); ++i) {
        update_representation(state, &current_dataset(*state).representation.reps[i]);
    }
}

static inline bool rep_type_uses_atomic_colors(RepresentationType type) {
    return type != RepresentationType::ElectronicStructure;
}

void update_representation(ApplicationState* state, Representation* rep) {
    ASSERT(state);
    ASSERT(rep);

    if (!rep->enabled) return;
    if (!rep->needs_update) return;

    const auto& sys = current_dataset(*state).sys;
    size_t num_atoms = md_system_atom_count(&sys);

    md_allocator_i* frame_alloc = state->allocator.frame;
    md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(frame_alloc);
    defer { md_vm_arena_temp_end(tmp); };


    const size_t bytes = num_atoms * sizeof(uint32_t);

    //md_script_property_t prop = {0};
    //if (rep->color_mapping == ColorMapping::Property) {
    //rep->prop_is_valid = md_script_compile_and_eval_property(&prop, rep->prop, &current_dataset(*data).sys, frame_allocator, &data->script.ir, rep->prop_error.beg(), rep->prop_error.capacity());
    //}

    uint32_t* colors = 0;
    if (rep_type_uses_atomic_colors(rep->type)) {
        colors = (uint32_t*)md_vm_arena_push(frame_alloc, sizeof(uint32_t) * num_atoms);

        switch (rep->color_mapping) {
        case ColorMapping::Uniform:
            color_atoms_uniform(colors, num_atoms, rep->uniform_color);
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
        case ColorMapping::SecondaryStructure:
            color_atoms_sec_str(colors, num_atoms, sys);
            break;
        case ColorMapping::Property:
            // @TODO: Map colors accordingly
            //color_atoms_uniform(colors, mol.atom.count, rep->uniform_color);

            if (md_array_size(current_dataset(*state).representation.info.atom_properties) > 0) {
                float* values = (float*)md_vm_arena_push(frame_alloc, sizeof(float) * num_atoms);
                EvalAtomProperty eval = {
                    .property_id = current_dataset(*state).representation.info.atom_properties[rep->prop.idx].id,
                    .idx = 0,
                    .output_written = false,
                    .num_values = num_atoms,
                    .dst_values = values,
                };
                viamd::event_system_broadcast_event(viamd::EventType_ViamdRepresentationEvalAtomProperty, viamd::EventPayloadType_EvalAtomProperty, &eval);

                if (eval.output_written) {
                    float range_ext = (rep->prop.range_end - rep->prop.range_beg);
                    range_ext = MAX(range_ext, 0.001f);
                    for (size_t i = 0; i < num_atoms; ++i) {
                        float t = (values[i] - rep->prop.range_beg) / range_ext;
                        colors[i] = ImPlot::SampleColormapU32(ImClamp(t, 0.0f, 1.0f), rep->prop.colormap);
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
                            .mol = &current_dataset(*state).sys,
                            .traj = current_dataset(*state).traj,
                        };
                        result = md_script_vis_eval_payload(&vis, rep->prop->vis_payload, 0, &ctx, MD_SCRIPT_VISUALIZE_ATOMS);
                    }
                    //}
                    if (result) {
                        if (dim == (int)md_array_size(vis.structure)) {
                            int i0 = CLAMP((int)current_dataset(*state).animation.frame + 0, 0, (int)rep->prop->data.num_values / dim - 1);
                            int i1 = CLAMP((int)current_dataset(*state).animation.frame + 1, 0, (int)rep->prop->data.num_values / dim - 1);
                            float frame_fract = fractf((float)current_dataset(*state).animation.frame);

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
                    int i0 = CLAMP((int)current_dataset(*state).animation.frame + 0, 0, (int)rep->prop->data.num_values - 1);
                    int i1 = CLAMP((int)current_dataset(*state).animation.frame + 1, 0, (int)rep->prop->data.num_values - 1);
                    float value = lerpf(values[i0], values[i1], fractf((float)current_dataset(*state).animation.frame));
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

    if (rep->saturation != 1.0f) {
        scale_saturation(colors, num_atoms, rep->saturation);
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
        size_t num_mos = current_dataset(*state).representation.info.alpha.num_orbitals;
        rep->type_is_valid = num_mos > 0;
        if (num_mos > 0) {
            uint64_t orb_idx = 0;
            uint64_t sub_idx = 0;

            switch(rep->electronic_structure.type) {
            case ElectronicStructureType::MolecularOrbital:
            case ElectronicStructureType::MolecularOrbitalDensity:
                orb_idx = rep->electronic_structure.mo_idx;
                break;
            case ElectronicStructureType::NaturalTransitionOrbitalParticle:
            case ElectronicStructureType::NaturalTransitionOrbitalHole:
            case ElectronicStructureType::NaturalTransitionOrbitalDensityParticle:
            case ElectronicStructureType::NaturalTransitionOrbitalDensityHole:
                orb_idx = rep->electronic_structure.nto_idx;
                sub_idx = rep->electronic_structure.nto_lambda_idx;
                break;
            case ElectronicStructureType::AttachmentDensity:
            case ElectronicStructureType::DetachmentDensity:
                orb_idx = rep->electronic_structure.nto_idx;
                break;
            default:
                break;
            }
            uint64_t vol_hash = (uint64_t)rep->electronic_structure.type | ((uint64_t)rep->electronic_structure.resolution << 8) | (orb_idx << 24) | (sub_idx << 48);
            if (vol_hash != rep->electronic_structure.vol_hash) {
                const float samples_per_angstrom[(int)VolumeResolution::Count] = {
                    4.0f,
                    8.0f,
                    16.0f,
                };
                EvalElectronicStructure data = {
                    .system = &current_dataset(*state).sys,
                    .type = rep->electronic_structure.type,
                    .major_idx = (int)orb_idx,
                    .minor_idx = (int)sub_idx,
                    .samples_per_angstrom = samples_per_angstrom[(int)rep->electronic_structure.resolution],
                    .dst_volume = &rep->electronic_structure.vol,
                };
                viamd::event_system_broadcast_event(viamd::EventType_ViamdRepresentationEvalElectronicStructure, viamd::EventPayloadType_EvalElectronicStructure, &data);

                if (data.output_written) {
                    rep->electronic_structure.vol_hash = vol_hash;
                }
            }
        }
        uint64_t tf_hash = md_hash64(&rep->electronic_structure.dvr.colormap, sizeof(rep->electronic_structure.dvr.colormap), 0);
        if (tf_hash != rep->electronic_structure.tf_hash) {
            rep->electronic_structure.tf_hash = tf_hash;
            volume::compute_transfer_function_texture_simple(&rep->electronic_structure.dvr.tf_tex, rep->electronic_structure.dvr.colormap);
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
			rep->filt_is_valid = md_filter(&rep->atom_mask, str_from_cstr(rep->filt), &current_dataset(*state).sys, current_dataset(*state).sys.atom.x, current_dataset(*state).sys.atom.y, current_dataset(*state).sys.atom.z, state->script.ir, &rep->filt_is_dynamic, rep->filt_error, sizeof(rep->filt_error));
            rep->filt_is_dirty = false;
        }

        if (rep->filt_is_valid) {
            filter_colors(colors, num_atoms, &rep->atom_mask);
            current_dataset(*state).representation.atom_visibility_mask_dirty = true;
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
    md_allocator_i* alloc = current_dataset(*state).representation.info.alloc;
    md_arena_allocator_reset(alloc);

    // Clear info, maintain allocator
    current_dataset(*state).representation.info = {};
    current_dataset(*state).representation.info.alloc = alloc;

    // Broadcast event to populate info
    viamd::event_system_broadcast_event(viamd::EventType_ViamdRepresentationInfoFill, viamd::EventPayloadType_RepresentationInfo, &current_dataset(*state).representation.info);
}

static void init_all_representations(ApplicationState* state) {
    for (size_t i = 0; i < md_array_size(current_dataset(*state).representation.reps); ++i) {
        auto& rep = current_dataset(*state).representation.reps[i];
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
    for (size_t i = 0; i < md_array_size(current_dataset(*state).representation.reps); ++i) {
        flag_representation_as_dirty(&current_dataset(*state).representation.reps[i]);
    }
}

void remove_all_representations(ApplicationState* state) {
    while (md_array_size(current_dataset(*state).representation.reps) > 0) {
        remove_representation(state, (int32_t)md_array_size(current_dataset(*state).representation.reps) - 1);
    }
}


void create_default_representations(ApplicationState* state) {
    bool amino_acid_present = false;
    bool nucleic_present = false;
    bool ion_present = false;
    bool water_present = false;
    bool ligand_present = false;
    bool coarse_grained = false;
    bool orbitals_present = current_dataset(*state).representation.info.alpha.num_orbitals > 0;

    if (current_dataset(*state).sys.atom.count > 3'000'000) {
        VIAMD_LOG_INFO("Large system detected, creating default representation for all atoms");
        Representation* rep = create_representation(state, RepresentationType::SpaceFill, ColorMapping::Type, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    if (current_dataset(*state).sys.component.count == 0) {
        // No residues present
        Representation* rep = create_representation(state, RepresentationType::BallAndStick, ColorMapping::Type, STR_LIT("all"));
        snprintf(rep->name, sizeof(rep->name), "default");
        goto done;
    }

    // TODO: Redo this check with entities instead of atom flags
    for (size_t i = 0; i < current_dataset(*state).sys.atom.count; ++i) {
        uint32_t flags = current_dataset(*state).sys.atom.flags[i];
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

        if (current_dataset(*state).sys.instance.count > 1) {
            color = ColorMapping::InstId;
        } else {
            size_t res_count = md_instance_comp_count(&current_dataset(*state).sys.instance, 0);
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
    auto& sys = current_dataset(*state).sys;
    const auto& traj = current_dataset(*state).traj;

    if (!sys.atom.count || !md_trajectory_num_frames(traj)) return;

    const int64_t last_frame = MAX(0LL, (int64_t)md_trajectory_num_frames(traj) - 1);
    // This is not actually time, but the fractional frame representation
    const double time = CLAMP(current_dataset(*state).animation.frame, 0.0, double(last_frame));

    // Scaling factor for cubic spline
    const int64_t frame = (int64_t)time;
    const int64_t nearest_frame = CLAMP((int64_t)(time + 0.5), 0LL, last_frame);

    static int64_t curr_nearest_frame = -1;
    if (current_dataset(*state).animation.interpolation == InterpolationMode::Nearest) {
        if (curr_nearest_frame == nearest_frame) {
            return;
        }
        current_dataset(*state).dirty_gpu_buffers |= MolBit_ClearVelocity;
    }
    curr_nearest_frame = nearest_frame;

    const int64_t frames[4] = {
        MAX(0LL, frame - 1),
        MAX(0LL, frame),
        MIN(frame + 1, last_frame),
        MIN(frame + 2, last_frame),
    };

    const size_t num_threads = task_system::pool_num_threads();

    const size_t stride = ALIGN_TO(sys.atom.count, 16);    // The interploation uses SIMD vectorization without bounds, so we make sure there is no overlap between the data segments
    const size_t bytes = stride * sizeof(float) * 3 * 4;
    
    // The number of atoms to be processed per thread when divided into chunks
    const uint32_t grain_size = 1024;

    md_allocator_i* temp_arena = state->allocator.frame;
    md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(temp_arena);
    defer { md_vm_arena_temp_end(tmp); };

    void* mem = md_vm_arena_push(temp_arena, bytes);
    struct Payload {
        ApplicationState* state;
        float s;
        float t;
        InterpolationMode mode;

        int64_t nearest_frame;
        int64_t frames[4];
        md_trajectory_frame_header_t headers[4];
        md_unitcell_t unitcell;

        size_t count;

        float* src_x[4];
        float* src_y[4];
        float* src_z[4];

        float* dst_x;
        float* dst_y;
        float* dst_z;

        vec3_t* aabb_min;
        vec3_t* aabb_max;
    };

    const InterpolationMode mode = (frames[1] != frames[2]) ? current_dataset(*state).animation.interpolation : InterpolationMode::Nearest;

    Payload payload = {
        .state = state,
        .s = 1.0f - CLAMP(current_dataset(*state).animation.tension, 0.0f, 1.0f),
        .t = (float)fractf(time),
        .mode = mode,
        .nearest_frame = nearest_frame,
        .frames = { frames[0], frames[1], frames[2], frames[3]},
        .count = sys.atom.count,
        .src_x = { (float*)mem + stride * 0, (float*)mem + stride * 1, (float*)mem + stride * 2,  (float*)mem + stride * 3 },
        .src_y = { (float*)mem + stride * 4, (float*)mem + stride * 5, (float*)mem + stride * 6,  (float*)mem + stride * 7 },
        .src_z = { (float*)mem + stride * 8, (float*)mem + stride * 9, (float*)mem + stride * 10, (float*)mem + stride * 11},
        .dst_x = sys.atom.x,
        .dst_y = sys.atom.y,
        .dst_z = sys.atom.z,
        .aabb_min = (vec3_t*)md_vm_arena_push(temp_arena, num_threads * sizeof(vec3_t)),
        .aabb_max = (vec3_t*)md_vm_arena_push(temp_arena, num_threads * sizeof(vec3_t)),
    };

    // This holds the chain of tasks we are about to submit
    task_system::ID tasks[16] = {0};
    int num_tasks = 0;

    switch (mode) {
        case InterpolationMode::Nearest: {
            task_system::ID load_task = task_system::create_pool_task(STR_LIT("## Load Frame"),[data = &payload]() {
                md_trajectory_frame_header_t header;
                md_trajectory_load_frame(current_dataset(*data->state).traj, data->nearest_frame, &header, data->dst_x, data->dst_y, data->dst_z);
                MEMCPY(&data->unitcell, &header.unitcell, sizeof(md_unitcell_t));
            });

            tasks[num_tasks++] = load_task;
            break;
        }
        case InterpolationMode::Linear: {
            task_system::ID load_task = task_system::create_pool_task(STR_LIT("## Load Frame"), 2, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                for (uint32_t i = range_beg; i < range_end; ++i) {
                    md_trajectory_load_frame(current_dataset(*data->state).traj, data->frames[i+1], &data->headers[i], data->src_x[i], data->src_y[i], data->src_z[i]);
                }
            });

            task_system::ID interp_unit_cell_task = task_system::create_pool_task(STR_LIT("## Interp Unit Cell Data"), [data = &payload]() {
                if (data->headers[0].unitcell.flags) {
                    double x  = lerp(data->headers[0].unitcell.x,  data->headers[1].unitcell.x,  data->t);
                    double y  = lerp(data->headers[0].unitcell.y,  data->headers[1].unitcell.y,  data->t);
                    double z  = lerp(data->headers[0].unitcell.z,  data->headers[1].unitcell.z,  data->t);
                    double xy = lerp(data->headers[0].unitcell.xy, data->headers[1].unitcell.xy, data->t);
                    double xz = lerp(data->headers[0].unitcell.xz, data->headers[1].unitcell.xz, data->t);
                    double yz = lerp(data->headers[0].unitcell.yz, data->headers[1].unitcell.yz, data->t);
                    data->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
                }
            });

            task_system::ID interp_coord_task = task_system::create_pool_task(STR_LIT("## Interp Coord Data"), (uint32_t)sys.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                float* dst_x = data->dst_x + range_beg;
                float* dst_y = data->dst_y + range_beg;
                float* dst_z = data->dst_z + range_beg;
                const float* src_x[2] = { data->src_x[0] + range_beg, data->src_x[1] + range_beg};
                const float* src_y[2] = { data->src_y[0] + range_beg, data->src_y[1] + range_beg};
                const float* src_z[2] = { data->src_z[0] + range_beg, data->src_z[1] + range_beg};

                md_util_interpolate_linear(dst_x, dst_y, dst_z, src_x, src_y, src_z, count, &data->unitcell, data->t);
            }, grain_size);

            tasks[num_tasks++] = load_task;
            tasks[num_tasks++] = interp_unit_cell_task;
            tasks[num_tasks++] = interp_coord_task;

            break;
        }
        case InterpolationMode::CubicSpline: {
            task_system::ID load_task = task_system::create_pool_task(STR_LIT("## Load Frame"), 4, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                for (uint32_t i = range_beg; i < range_end; ++i) {
                    md_trajectory_load_frame(current_dataset(*data->state).traj, data->frames[i], &data->headers[i], data->src_x[i], data->src_y[i], data->src_z[i]);
                }
            });

            task_system::ID interp_unit_cell_task = task_system::create_pool_task(STR_LIT("## Interp Unit Cell Data"), [data = &payload]() {
                if (data->headers[0].unitcell.flags) {
                    double x  = cubic_spline(data->headers[0].unitcell.x,  data->headers[1].unitcell.x,  data->headers[2].unitcell.x,  data->headers[3].unitcell.x,  data->t, data->s);
                    double y  = cubic_spline(data->headers[0].unitcell.y,  data->headers[1].unitcell.y,  data->headers[2].unitcell.y,  data->headers[3].unitcell.y,  data->t, data->s);
                    double z  = cubic_spline(data->headers[0].unitcell.z,  data->headers[1].unitcell.z,  data->headers[2].unitcell.z,  data->headers[3].unitcell.z,  data->t, data->s);
                    double xy = cubic_spline(data->headers[0].unitcell.xy, data->headers[1].unitcell.xy, data->headers[2].unitcell.xy, data->headers[3].unitcell.xy, data->t, data->s);
                    double xz = cubic_spline(data->headers[0].unitcell.xz, data->headers[1].unitcell.xz, data->headers[2].unitcell.xz, data->headers[3].unitcell.xz, data->t, data->s);
                    double yz = cubic_spline(data->headers[0].unitcell.yz, data->headers[1].unitcell.yz, data->headers[2].unitcell.yz, data->headers[3].unitcell.yz, data->t, data->s);
                    
                    data->unitcell = md_unitcell_from_basis_parameters(x, y, z, xy, xz, yz);
                }
            });

            task_system::ID interp_coord_task = task_system::create_pool_task(STR_LIT("## Interp Coord Data"), (uint32_t)sys.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                (void)thread_num;
                size_t count = range_end - range_beg;
                float* dst_x = data->dst_x + range_beg;
                float* dst_y = data->dst_y + range_beg;
                float* dst_z = data->dst_z + range_beg;
                const float* src_x[4] = { data->src_x[0] + range_beg, data->src_x[1] + range_beg, data->src_x[2] + range_beg, data->src_x[3] + range_beg};
                const float* src_y[4] = { data->src_y[0] + range_beg, data->src_y[1] + range_beg, data->src_y[2] + range_beg, data->src_y[3] + range_beg};
                const float* src_z[4] = { data->src_z[0] + range_beg, data->src_z[1] + range_beg, data->src_z[2] + range_beg, data->src_z[3] + range_beg};

                md_util_interpolate_cubic_spline(dst_x, dst_y, dst_z, src_x, src_y, src_z, count, &data->unitcell, data->t, data->s);
            }, grain_size);

            tasks[num_tasks++] = load_task;
            tasks[num_tasks++] = interp_unit_cell_task;
            tasks[num_tasks++] = interp_coord_task;
            
            break;
        }
        default:
            ASSERT(false);
            break;
    }

    if (current_dataset(*state).operations.recalc_bonds) {
        if (!task_system::task_is_running(state->tasks.evaluate_full) && !task_system::task_is_running(state->tasks.evaluate_filt)) {
            // We cannot recalculate bonds while the full or filtered evaluation is running
            // because it would overwrite the bond data while we are reading it

            static int64_t cur_nearest_frame = -1;
            if (cur_nearest_frame != payload.nearest_frame) {
                cur_nearest_frame = payload.nearest_frame;
                task_system::ID recalc_bond_task = task_system::create_pool_task(STR_LIT("## Recalc bond task"), [data = &payload]() {
                    const auto& sys = current_dataset(*data->state).sys;
                    const float* x = sys.atom.x;
                    const float* y = sys.atom.y;
                    const float* z = sys.atom.z;
                    const md_unitcell_t* cell = &sys.unitcell;
                    int offset = data->t < 0.5 ? 0 : 1;

                    switch (data->mode) {
                    case InterpolationMode::Nearest: break;
                    case InterpolationMode::Linear:
                        x = data->src_x[0 + offset];
                        y = data->src_y[0 + offset];
                        z = data->src_z[0 + offset];
                        cell = &data->headers[0 + offset].unitcell;
                        break;
                    case InterpolationMode::CubicSpline:
                        x = data->src_x[1 + offset];
                        y = data->src_y[1 + offset];
                        z = data->src_z[1 + offset];
                        cell = &data->headers[1 + offset].unitcell;
                        break;
                    default:
                        break;
                    };

                    md_bond_data_t* bonds = &current_dataset(*data->state).sys.bond;
                    md_bond_data_clear(bonds);

                    md_util_infer_covalent_bonds(bonds, x, y, z, cell, &sys, current_dataset(*data->state).alloc);
                    current_dataset(*data->state).dirty_gpu_buffers |= MolBit_DirtyBonds;
                });
                tasks[num_tasks++] = recalc_bond_task;
            }
        }
    }

    if (current_dataset(*state).operations.apply_pbc) {
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
    if (current_dataset(*state).operations.unwrap_structures) {
        size_t num_structures = md_index_data_num_ranges(&sys.structure);
        task_system::ID unwrap_task = task_system::create_pool_task(STR_LIT("## Unwrap Structures"), (uint32_t)num_structures, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            for (uint32_t i = range_beg; i < range_end; ++i) {
                int*     s_idx = md_index_range_ptr(&current_dataset(*data->state).sys.structure, i);
                size_t   s_len = md_index_range_size(&current_dataset(*data->state).sys.structure, i);
                md_util_unwrap(data->dst_x, data->dst_y, data->dst_z, s_idx, s_len, &data->unitcell);
            }
        });
        tasks[num_tasks++] = unwrap_task;
    }

    {
        // Calculate a global AABB for the molecule
        task_system::ID aabb_task = task_system::create_pool_task(STR_LIT("## Compute AABB"), (uint32_t)sys.atom.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            size_t range_len = range_end - range_beg;
            const float* x = current_dataset(*data->state).sys.atom.x + range_beg;
            const float* y = current_dataset(*data->state).sys.atom.y + range_beg;
            const float* z = current_dataset(*data->state).sys.atom.z + range_beg;

            size_t temp_pos = md_temp_get_pos();
            defer { md_temp_set_pos_back(temp_pos); };
            float* r = (float*)md_temp_push(sizeof(float) * range_len);
            md_atom_extract_radii(r, range_beg, range_len, &current_dataset(*data->state).sys.atom);

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
                        current_dataset(*data->state).trajectory_data.backbone_angles.data + current_dataset(*data->state).trajectory_data.backbone_angles.stride * data->frames[1],
                        current_dataset(*data->state).trajectory_data.backbone_angles.data + current_dataset(*data->state).trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    const md_backbone_angles_t* src_angle = data->t < 0.5f ? src_angles[0] : src_angles[1];
                    MEMCPY(current_dataset(*data->state).sys.protein_backbone.segment.angle, src_angle, current_dataset(*data->state).sys.protein_backbone.segment.count * sizeof(md_backbone_angles_t));
                });

                tasks[num_tasks++] = angle_task;
                break;
            }
            case InterpolationMode::Linear: {
                task_system::ID angle_task = task_system::create_pool_task(STR_LIT("## Compute Backbone Angles"), (uint32_t)sys.protein_backbone.segment.count, [data = &payload](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
                    (void)thread_num;
                    const md_backbone_angles_t* src_angles[2] = {
                        current_dataset(*data->state).trajectory_data.backbone_angles.data + current_dataset(*data->state).trajectory_data.backbone_angles.stride * data->frames[1],
                        current_dataset(*data->state).trajectory_data.backbone_angles.data + current_dataset(*data->state).trajectory_data.backbone_angles.stride * data->frames[2],
                    };
                    md_system_t& mol = current_dataset(*data->state).sys;
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
                        current_dataset(*data->state).trajectory_data.backbone_angles.data + current_dataset(*data->state).trajectory_data.backbone_angles.stride * data->frames[0],
                        current_dataset(*data->state).trajectory_data.backbone_angles.data + current_dataset(*data->state).trajectory_data.backbone_angles.stride * data->frames[1],
                        current_dataset(*data->state).trajectory_data.backbone_angles.data + current_dataset(*data->state).trajectory_data.backbone_angles.stride * data->frames[2],
                        current_dataset(*data->state).trajectory_data.backbone_angles.data + current_dataset(*data->state).trajectory_data.backbone_angles.stride * data->frames[3],
                    };
                    md_system_t& sys = current_dataset(*data->state).sys;
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
        if (md_array_size(current_dataset(*state).interpolated_properties.secondary_structure) != sys.protein_backbone.segment.count) {
			MD_LOG_ERROR("Secondary structure array size does not match the number of segments.");
        }
        size_t num_backbone_segments = sys.protein_backbone.segment.count;
        task_system::ID ss_task = task_system::create_pool_task(STR_LIT("## Interpolate Secondary Structures"), (uint32_t)num_backbone_segments, [data = &payload, mode](uint32_t range_beg, uint32_t range_end, uint32_t thread_num) {
            (void)thread_num;
            const md_secondary_structure_t* src_ss[4] = {
                (md_secondary_structure_t*)current_dataset(*data->state).trajectory_data.secondary_structure.data + current_dataset(*data->state).trajectory_data.secondary_structure.stride * data->frames[0],
                (md_secondary_structure_t*)current_dataset(*data->state).trajectory_data.secondary_structure.data + current_dataset(*data->state).trajectory_data.secondary_structure.stride * data->frames[1],
                (md_secondary_structure_t*)current_dataset(*data->state).trajectory_data.secondary_structure.data + current_dataset(*data->state).trajectory_data.secondary_structure.stride * data->frames[2],
                (md_secondary_structure_t*)current_dataset(*data->state).trajectory_data.secondary_structure.data + current_dataset(*data->state).trajectory_data.secondary_structure.stride * data->frames[3],
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
                    current_dataset(*data->state).sys.protein_backbone.segment.secondary_structure[i] = ss;
                    current_dataset(*data->state).interpolated_properties.secondary_structure[i] = md_gl_secondary_structure_convert(ss);
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
                    current_dataset(*data->state).sys.protein_backbone.segment.secondary_structure[i] = src_ss_nearest[i];
                    current_dataset(*data->state).interpolated_properties.secondary_structure[i] = ss_gl_i;
                }
                break;
            }
            case InterpolationMode::CubicSpline: {
                const md_gl_secondary_structure_t ss_coil = { 0,0 };
                const md_gl_secondary_structure_t ss_sheet = { .sheet = 1.0f };

                for (size_t i = range_beg; i < range_end; ++i) {
                    md_secondary_structure_t ss[4] = { src_ss[0][i], src_ss[1][i], src_ss[2][i], src_ss[3][i] };

                    md_gl_secondary_structure_t ss_gl[4] = {
                        md_gl_secondary_structure_convert(ss[0]),
                        md_gl_secondary_structure_convert(ss[1]),
                        md_gl_secondary_structure_convert(ss[2]),
                        md_gl_secondary_structure_convert(ss[3]),
                    };

                    auto is_eq = [](md_gl_secondary_structure_t a, md_gl_secondary_structure_t b) {
                        return a.helix == b.helix && a.sheet == b.sheet;
                    };

                    // Cleanup isolated coils temporally to reduce noise during transitions.
                    // This is a common issue with secondary structure assignment during transitions, where a segment might flip between coil and helix/sheet rapidly, creating a flickering effect.
                    // We compare indices 0, 1, 2, 3 and if idx 1 or 2 differs from the other 3 (which are the same), we set 1 to be the same as the others, effectively removing isolated coil assignments.
#if 0
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
                    current_dataset(*data->state).sys.protein_backbone.segment.secondary_structure[i] = src_ss_nearest[i];
                    current_dataset(*data->state).interpolated_properties.secondary_structure[i] = ss_gl_i;
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

            md_gl_secondary_structure_t* ss_gl = current_dataset(*data->state).interpolated_properties.secondary_structure;
            for (size_t i = 0; i < current_dataset(*data->state).sys.protein_backbone.range.count; ++i) {
                size_t range_beg = current_dataset(*data->state).sys.protein_backbone.range.offset[i];
                size_t range_end = current_dataset(*data->state).sys.protein_backbone.range.offset[i + 1];
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

        current_dataset(*state).dirty_gpu_buffers |= MolBit_DirtySecondaryStructure;
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
    current_dataset(*state).aabb_min = aabb_min;
    current_dataset(*state).aabb_max = aabb_max;
    sys.unitcell = payload.unitcell;

#if 0
    if (mol.unitcell.flags) {
        vec3_t c = mol.unitcell.basis * vec3_set1(0.5f);
        current_dataset(*state).model_mat = mat4_translate_vec3(-c);
    }
#endif

    current_dataset(*state).dirty_gpu_buffers |= MolBit_DirtyPosition;
}

void reset_view(ApplicationState* state, const md_bitfield_t* target, bool move_camera, bool smooth_transition) {
    ASSERT(state);

    md_allocator_i* temp_arena = state->allocator.frame;
    md_vm_arena_temp_t tmp = md_vm_arena_temp_begin(temp_arena);
    defer { md_vm_arena_temp_end(tmp); };

    Dataset& current = current_dataset(*state);

    if (!current.sys.atom.count) return;
    const auto& sys = current.sys;
    size_t popcount = 0;
    if (target) {
        popcount = md_bitfield_popcount(target);
    }

    int32_t* indices = nullptr;
    if (0 < popcount && popcount < sys.atom.count) {
        indices = (int32_t*)md_vm_arena_push_array(temp_arena, int32_t, popcount);
        size_t len = md_bitfield_iter_extract_indices(indices, popcount, md_bitfield_iter_create(target));
        if (len > popcount || len > sys.atom.count) {
            MD_LOG_DEBUG("Error: Invalid number of indices");
            len = MIN(popcount, sys.atom.count);
        }
    }

    size_t count = popcount ? popcount : sys.atom.count;
    vec3_t com = md_util_com_compute(sys.atom.x, sys.atom.y, sys.atom.z, nullptr, indices, count, &sys.unitcell);
    mat3_t C = mat3_covariance_matrix(sys.atom.x, sys.atom.y, sys.atom.z, nullptr, indices, count, com);
    mat3_eigen_t eigen = mat3_eigen(C);
    mat3_t PCA = mat3_orthonormalize(mat3_extract_rotation(eigen.vectors));
    mat4_t Ri  = mat4_from_mat3(PCA);

    // Compute min and maximum extent along the PCA axes
    vec3_t min_ext = vec3_set1( FLT_MAX);
    vec3_t max_ext = vec3_set1(-FLT_MAX);
    mat3_t basis = mat3_transpose(PCA);

    const float radius = 1.0f;

    // Transform the gto (x,y,z,cutoff) into the PCA frame to find the min and max extend within it
    for (size_t i = 0; i < count; ++i) {
        int32_t idx = indices ? indices[i] : (int32_t)i;
        vec3_t xyz = { sys.atom.x[idx], sys.atom.y[idx], sys.atom.z[idx] };

        vec3_t p = mat4_mul_vec3(Ri, xyz, 1.0f);
        min_ext = vec3_min(min_ext, vec3_sub_f(p, radius));
        max_ext = vec3_max(max_ext, vec3_add_f(p, radius));
    }

    vec3_t optimal_pos;
    quat_t optimal_ori;
    float  optimal_dist;
    camera_compute_optimal_view(&optimal_pos, &optimal_ori, &optimal_dist, basis, min_ext, max_ext);

    if (move_camera) {
        current.view.animation.target_position    = optimal_pos;
        current.view.animation.target_orientation = optimal_ori;
        current.view.animation.target_distance    = optimal_dist;

        if (!smooth_transition) {
            current.view.camera.position       = optimal_pos;
            current.view.camera.orientation    = optimal_ori;
            current.view.camera.focus_distance = optimal_dist;
        }
    }

    const vec3_t cell_ext = mat3_mul_vec3(md_unitcell_basis_mat3(&sys.unitcell), vec3_set1(1.0f));
    const float  max_cell_ext = vec3_reduce_max(cell_ext);
    const float  max_aabb_ext = vec3_reduce_max(vec3_sub(max_ext, min_ext));

    current.view.camera.near_plane = 1.0f;
    current.view.camera.far_plane = 100000.0f;
    current.view.trackball_param.max_distance = MAX(max_cell_ext, max_aabb_ext) * 10.0f;
}
