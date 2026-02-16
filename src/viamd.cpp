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

static void init_all_representations(ApplicationState* state);

void draw_info_window(const ApplicationState& state, uint32_t picking_idx) {
    const auto& sys = state.mold.sys;

    if (picking_idx == INVALID_PICKING_IDX) return;

    md_strb_t sb = md_strb_create(state.allocator.frame);
    defer {
        md_strb_free(&sb);
    };

    if (picking_idx < sys.atom.count) {
        int atom_idx = picking_idx;
        int local_idx = atom_idx;
        const vec3_t pos = md_atom_coord(&sys.atom, atom_idx);
        str_t type = md_atom_name(&sys.atom, atom_idx);
        md_atomic_number_t z = md_atom_atomic_number(&sys.atom, atom_idx);
        str_t elem = z ? md_util_element_name(z)   : str_t{};
        str_t symb = z ? md_util_element_symbol(z) : str_t{};

        int comp_idx = md_comp_find_by_atom_idx(&sys.comp, atom_idx);
        str_t comp_name = {};
        int comp_seq_id = 0;
        if (comp_idx != -1) {
            comp_name   = md_comp_name(&sys.comp, comp_idx);
            comp_seq_id = md_comp_seq_id(&sys.comp, comp_idx);
            md_urange_t range = md_comp_atom_range(&sys.comp, comp_idx);
            local_idx = atom_idx - range.beg;
        }

        int inst_idx = md_system_inst_find_by_atom_idx(&sys, atom_idx);
        str_t chain_id = {};
        if (inst_idx != -1) {
            chain_id = md_inst_id(&sys.inst, inst_idx);
        }

        // External indices begin with 1 not 0
        md_strb_fmt(&sb, "atom[%i][%i]: %.*s %.*s %.*s (%.2f, %.2f, %.2f)\n", atom_idx + 1, local_idx + 1, STR_ARG(type), STR_ARG(elem), STR_ARG(symb), pos.x, pos.y, pos.z);
        if (comp_idx != -1) {
            md_strb_fmt(&sb, "res[%i]: %.*s %i\n", comp_idx + 1, STR_ARG(comp_name), comp_seq_id);
        }
        if (inst_idx != -1) {
            md_strb_fmt(&sb, "chain[%i]: %.*s\n", inst_idx + 1, STR_ARG(chain_id));
        }

        uint32_t flags = sys.atom.flags[atom_idx];
        if (flags) {
            sb += "flags: ";
            if (flags & MD_FLAG_HETERO)             { sb += "HETERO "; }
            if (flags & MD_FLAG_AMINO_ACID)         { sb += "AMINO "; }
            //if (flags & MD_FLAG_SIDE_CHAIN)         { sb += "SIDE_CHAIN "; }
            if (flags & MD_FLAG_NUCLEOTIDE)         { sb += "NUCLEOTIDE "; }
            if (flags & MD_FLAG_NUCLEOBASE)         { sb += "NUCLEOBASE "; }
            //if (flags & MD_FLAG_NUCLEOSIDE)         { sb += "NUCLEOSIDE "; }
            if (flags & MD_FLAG_WATER)              { sb += "WATER "; }
            if (flags & MD_FLAG_ION)                { sb += "ION "; }
            if (flags & MD_FLAG_BACKBONE)           { sb += "BACKBONE "; }
            if (flags & MD_FLAG_SP)                 { sb += "SP "; }
            if (flags & MD_FLAG_SP2)                { sb += "SP2 "; }
            if (flags & MD_FLAG_SP3)                { sb += "SP3 "; }
            if (flags & MD_FLAG_AROMATIC)           { sb += "AROMATIC "; }
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
            md_atom_pair_t pair = sys.bond.pairs[bond_idx];
            md_flags_t    flags = sys.bond.flags[bond_idx];
            char bond_flags_buf[256] = {};
            int  len = 0;

            typedef struct {
                md_flags_t flag;
                const char* label;
            } bond_flag_label_t;

            bond_flag_label_t bond_flag_map[] = {
                {MD_BOND_FLAG_COVALENT,     "COVALENT"},
                {MD_BOND_FLAG_DOUBLE,       "DOUBLE"},
                {MD_BOND_FLAG_TRIPLE,       "TRIPLE"},
                {MD_BOND_FLAG_QUADRUPLE,    "QUADRUPLE"},
                {MD_BOND_FLAG_AROMATIC,     "AROMATIC"},
                {MD_BOND_FLAG_INTER,        "INTER"},
                {MD_BOND_FLAG_COORDINATE,   "COORD"},
				{MD_BOND_FLAG_METAL,        "METAL"},
				{MD_BOND_FLAG_USER_DEFINED, "USER"},
            };

            for (size_t i = 0; i < ARRAY_SIZE(bond_flag_map); ++i) {
                if (flags & bond_flag_map[i].flag) {
                    len += snprintf(bond_flags_buf + len, 256 - len, "%s ", bond_flag_map[i].label);
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

    if (state->mold.traj) {
        load::traj::close(state->mold.traj);
        state->mold.traj = nullptr;
    }
    state->files.trajectory[0] = '\0';

    state->mold.sys.unitcell = {};
    md_array_free(state->timeline.x_values,  state->allocator.persistent);
    md_array_free(state->display_properties, state->allocator.persistent);

    md_array_free(state->trajectory_data.backbone_angles.data,     state->allocator.persistent);
    md_array_free(state->trajectory_data.secondary_structure.data, state->allocator.persistent);
}

void init_trajectory_data(ApplicationState* data) {
    size_t num_frames = md_trajectory_num_frames(data->mold.traj);
    if (num_frames > 0) {
        size_t min_frame = 0;
        size_t max_frame = num_frames - 1;
        md_trajectory_header_t header;
        md_trajectory_get_header(data->mold.traj, &header);

        double min_time = header.frame_times[0];
        double max_time = header.frame_times[num_frames - 1];

        data->timeline.view_range = {min_time, max_time};
        data->timeline.filter.beg_frame = min_frame;
        data->timeline.filter.end_frame = max_frame;

        md_array_resize(data->timeline.x_values, num_frames, data->allocator.persistent);
        for (size_t i = 0; i < num_frames; ++i) {
            data->timeline.x_values[i] = header.frame_times[i];
        }

        data->animation.frame = CLAMP(data->animation.frame, (double)min_frame, (double)max_frame);
        int64_t frame_idx = CLAMP((int64_t)(data->animation.frame + 0.5), 0, (int64_t)max_frame);

        md_trajectory_frame_header_t frame_header;
        md_trajectory_load_frame(data->mold.traj, frame_idx, &frame_header, data->mold.sys.atom.x, data->mold.sys.atom.y, data->mold.sys.atom.z);
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
                md_trajectory_i* traj = load::traj::get_raw_trajectory(data->mold.traj);

                md_allocator_i* temp_arena = md_vm_arena_create(GIGABYTES(1));
                defer { md_vm_arena_destroy(temp_arena); };

                const size_t capacity = ALIGN_TO(sys.atom.count, 16);
                float* x = (float*)md_vm_arena_push(temp_arena, sizeof(float) * capacity);
                float* y = (float*)md_vm_arena_push(temp_arena, sizeof(float) * capacity);
                float* z = (float*)md_vm_arena_push(temp_arena, sizeof(float) * capacity);

                for (uint32_t frame_idx = range_beg; frame_idx < range_end; ++frame_idx) {
                    md_trajectory_frame_header_t frame_header;
                    md_backbone_angles_t* bb_dst = data->trajectory_data.backbone_angles.data + data->trajectory_data.backbone_angles.stride * frame_idx;
                    md_secondary_structure_t* ss_dst = data->trajectory_data.secondary_structure.data + data->trajectory_data.secondary_structure.stride * frame_idx;

                    md_trajectory_load_frame(traj, frame_idx, &frame_header, x, y, z);
                    md_util_backbone_angles_compute(bb_dst, data->trajectory_data.backbone_angles.stride, x, y, z, &frame_header.unitcell, &sys.protein_backbone);
                    md_util_backbone_secondary_structure_infer(ss_dst, data->trajectory_data.secondary_structure.stride, x, y, z, &frame_header.unitcell, &sys.protein_backbone);
                }
                });

            uint64_t time = (uint64_t)md_time_current();
            task_system::ID main_task = task_system::create_main_task(STR_LIT("Update Trajectory Data"), [data, t0 = time]() {
                uint64_t t1 = (uint64_t)md_time_current();
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

bool load_trajectory_data(ApplicationState* data, str_t filename, md_trajectory_loader_i* loader, LoadTrajectoryFlags flags) {
    md_trajectory_i* traj = load::traj::open_file(filename, loader, &data->mold.sys, data->allocator.persistent, flags);
    if (traj) {
        free_trajectory_data(data);
        data->mold.traj = traj;
        str_copy_to_char_buf(data->files.trajectory, sizeof(data->files.trajectory), filename);
        init_trajectory_data(data);
        data->animation.frame = 0;
        return true;
    }

    return false;
}

void init_molecule_data(ApplicationState* data) {
    if (data->mold.sys.atom.count) {

        data->picking.idx = INVALID_PICKING_IDX;
        data->selection.atom_idx.hovered = -1;
        data->selection.atom_idx.right_click = -1;
        data->selection.bond_idx.hovered = -1;
        data->selection.bond_idx.right_click = -1;

        data->mold.gl_mol = md_gl_mol_create(&data->mold.sys);
        if (data->mold.sys.protein_backbone.segment.count > 0) {
            data->interpolated_properties.secondary_structure = md_array_create(md_gl_secondary_structure_t, data->mold.sys.protein_backbone.segment.count, data->mold.sys_alloc);
        }

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

    //md_molecule_free(&data->mold.sys, persistent_alloc);
    md_arena_allocator_reset(data->mold.sys_alloc);
    MEMSET(&data->mold.sys, 0, sizeof(data->mold.sys));

    md_gl_mol_destroy(data->mold.gl_mol);
    MEMSET(data->files.molecule, 0, sizeof(data->files.molecule));

    data->interpolated_properties.secondary_structure = nullptr;

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

    viamd::event_system_broadcast_event(viamd::EventType_ViamdTopologyFree, viamd::EventPayloadType_ApplicationState, data);
}

bool load_dataset_from_file(ApplicationState* data, const LoadParam& param) {
    ASSERT(data);

    str_t path_to_file = md_path_make_canonical(param.file_path, data->allocator.frame);
    if (path_to_file) {
        if (param.sys_loader) {
            interrupt_async_tasks(data);
            free_trajectory_data(data);
            free_molecule_data(data);

            if (!param.sys_loader->init_from_file(&data->mold.sys, path_to_file, param.sys_loader_arg, data->mold.sys_alloc)) {
                VIAMD_LOG_ERROR("Failed to load molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));
                return false;
            }
            VIAMD_LOG_SUCCESS("Successfully loaded molecular data from file '" STR_FMT "'", STR_ARG(path_to_file));

            str_copy_to_char_buf(data->files.molecule, sizeof(data->files.molecule), path_to_file);
            data->files.coarse_grained = param.coarse_grained;
            // @NOTE: If the dataset is coarse-grained, then postprocessing must be aware
            md_util_postprocess_flags_t flags = param.coarse_grained ? 0 : MD_UTIL_POSTPROCESS_ALL;
            md_util_molecule_postprocess(&data->mold.sys, data->mold.sys_alloc, flags);
            init_molecule_data(data);

            // @NOTE: Some files contain both atomic coordinates and trajectory
            if (param.traj_loader) {
                VIAMD_LOG_INFO("File may also contain trajectory, attempting to load trajectory");
            } else {
                return true;
            }
        }

        if (param.traj_loader) {
            if (!data->mold.sys.atom.count) {
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
                } else if (str_eq(ident, STR_LIT("OrbIdx"))) {
                    viamd::extract_int(rep->electronic_structure.mo_idx, arg);
                } else if (str_eq(ident, STR_LIT("OrbMoIdx"))) {
                    viamd::extract_int(rep->electronic_structure.mo_idx, arg);
                } else if (str_eq(ident, STR_LIT("OrbNtoIdx"))) {
                    viamd::extract_int(rep->electronic_structure.nto_idx, arg);
                } else if (str_eq(ident, STR_LIT("OrbRes"))) {
                    int res;
                    viamd::extract_int(res, arg);
                    rep->electronic_structure.resolution = (VolumeResolution)res;
                } else if (str_eq(ident, STR_LIT("OrbType"))) {
                    int type;
                    viamd::extract_int(type, arg);
                    rep->electronic_structure.type = (ElectronicStructureType)type;
                } else if (str_eq(ident, STR_LIT("OrbIso"))) {
                    viamd::extract_flt(rep->electronic_structure.iso_psi.values[0], arg);
                    rep->electronic_structure.iso_psi.values[1] = - rep->electronic_structure.iso_psi.values[0];
                    rep->electronic_structure.iso_den.values[0] = rep->electronic_structure.iso_psi.values[0] * rep->electronic_structure.iso_psi.values[0];
                } else if (str_eq(ident, STR_LIT("OrbColPos"))) {
                    viamd::extract_vec4(rep->electronic_structure.iso_psi.colors[0], arg);
                } else if (str_eq(ident, STR_LIT("OrbColNeg"))) {
                    viamd::extract_vec4(rep->electronic_structure.iso_psi.colors[1], arg);
                } else if (str_eq(ident, STR_LIT("OrbColDen"))) {
                    viamd::extract_vec4(rep->electronic_structure.iso_den.colors[0], arg);
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
    
    str_copy_to_char_buf(data->files.workspace, sizeof(data->files.workspace), filename);
    
    data->files.coarse_grained  = new_coarse_grained;

    load::LoaderState loader_state = {};
    load::init_loader_state(&loader_state, new_molecule_file, data->allocator.frame);

    LoadParam param = {};
    param.file_path = new_molecule_file;
    param.sys_loader = loader_state.sys_loader;
    param.traj_loader = loader_state.traj_loader;
    param.coarse_grained = new_coarse_grained;
    param.sys_loader_arg = loader_state.sys_loader_arg;

    if (new_molecule_file && load_dataset_from_file(data, param)) {
        str_copy_to_char_buf(data->files.molecule, sizeof(data->files.molecule), new_molecule_file);
    } else {
        data->files.molecule[0]   = '\0';
    }

    if (new_trajectory_file) {
        load::init_loader_state(&loader_state, new_trajectory_file, data->allocator.frame);
        param.sys_loader = 0;
        param.traj_loader = loader_state.traj_loader;
        param.file_path = new_trajectory_file;
        if (load_dataset_from_file(data, param)) {
            str_copy_to_char_buf(data->files.trajectory, sizeof(data->files.trajectory), new_trajectory_file);
        }
        data->animation.frame = new_frame;
    } else {
        data->files.trajectory[0] = '\0';
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
    viamd::write_flt(state,  STR_LIT("Distance"), app_state->view.camera.focus_distance);
    viamd::write_int(state,  STR_LIT("Mode"), (int)app_state->view.mode);


    for (size_t i = 0; i < md_array_size(app_state->representation.reps); ++i) {
        const Representation& rep = app_state->representation.reps[i];
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
            viamd::write_flt(state,  STR_LIT("ElectronicStructureIso"),      rep.electronic_structure.iso_psi.values[0]);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColPos"),   rep.electronic_structure.iso_psi.colors[0]);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColNeg"),   rep.electronic_structure.iso_psi.colors[1]);
            viamd::write_vec4(state, STR_LIT("ElectronicStructureColDen"),   rep.electronic_structure.iso_den.colors[0]);
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
        rep->prop.idx = 0;
        rep->prop.range_beg = state->representation.info.atom_properties[0].value_min;
        rep->prop.range_end = state->representation.info.atom_properties[0].value_max;
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
    if (rep.electronic_structure.vol.tex_id != 0) gl::free_texture(&rep.electronic_structure.vol.tex_id);
    if (rep.electronic_structure.dvr.tf_tex != 0) gl::free_texture(&rep.electronic_structure.dvr.tf_tex);
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

static inline bool rep_type_uses_atomic_colors(RepresentationType type) {
    return type != RepresentationType::ElectronicStructure;
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

            if (md_array_size(state->representation.info.atom_properties) > 0) {
                float* values = (float*)md_vm_arena_push(frame_alloc, sizeof(float) * num_atoms);
                EvalAtomProperty eval = {
                    .property_id = state->representation.info.atom_properties[rep->prop.idx].id,
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
                            .mol = &state->mold.sys,
                            .traj = state->mold.traj,
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
        size_t num_mos = state->representation.info.alpha.num_orbitals;
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
                    .type = rep->electronic_structure.type,
                    .major_idx = (int)orb_idx,
                    .minor_idx = (int)sub_idx,
                    .samples_per_angstrom = samples_per_angstrom[(int)rep->electronic_structure.resolution],
                    .dst_volume = &rep->electronic_structure.vol,
                };
                viamd::event_system_broadcast_event(viamd::EventType_ViamdRepresentationEvalElectronicStructure, viamd::EventPayloadType_EvalElectronicStructure, &data);

                if (data.output_written) {
#if !VIAMD_RECOMPUTE_ORBITAL_PER_FRAME
                    rep->electronic_structure.vol_hash = vol_hash;
#endif
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

    if (state->mold.sys.comp.count == 0) {
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

        if (state->mold.sys.inst.count > 1) {
            color = ColorMapping::InstId;
        } else {
            size_t res_count = md_inst_comp_count(&state->mold.sys.inst, 0);
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