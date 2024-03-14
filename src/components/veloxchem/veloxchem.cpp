#define IMGUI_DEFINE_MATH_OPERATORS

#include <event.h>
#include <viamd.h>
#include <core/md_arena_allocator.h>
#include <task_system.h>

#include <core/md_log.h>
#include <md_vlx.h>
#include <md_gto.h>

#include <gfx/volumerender_utils.h>
#include <gfx/gl_utils.h>
#include <gfx/immediate_draw_utils.h>

#include <imgui_widgets.h>
#include <implot_widgets.h>

struct VeloxChem : viamd::EventHandler {
	VeloxChem() { viamd::event_system_register_handler(*this); }

	bool show_volume = false;
	bool show_window = false;
	md_vlx_data_t vlx{};
	int vol_dim = 128;
	mat4_t model_mat = {};

	float* vol_data = 0;
	vec3_t min_box = {};
	vec3_t max_box = {};
	vec3_t step_size = {};

	uint32_t vol_texture = 0;
	uint32_t tf_texture = 0;

	vec3_t clip_min = { 0, 0, 0 };
	vec3_t clip_max = { 1, 1, 1 };

	bool bounding_box_enabled = false;
	vec4_t bounding_box_color = { 0, 0, 0, 1 };

	int mo_idx = 0;
	int homo_idx = 0;
	int lumo_idx = 0;

	task_system::ID compute_volume_task = 0;

	md_gto_data_t pgto = { 0 };

	// double r_values[100] = {0};
	// double orb_values[5][100] = {0};
	// double r_density [3][100] = {0};

	float density_scale = 1.0f;

	struct {
		bool enabled = true;
		size_t count = 2;
		float values[8] = { 0.05f, -0.05 };
		vec4_t colors[8] = { {215.f / 255.f, 25.f / 255.f, 28.f / 255.f, 0.75f}, {44.f / 255.f, 123.f / 255.f, 182.f / 255.f, 0.75f} };
	} iso_surface;

	md_allocator_i* arena = 0;

	size_t num_orbitals() const { return vlx.scf.alpha.orbitals.dim[1]; }

	size_t num_cgtos() const { return vlx.scf.alpha.orbitals.dim[0]; }

	void process_events(const viamd::Event* events, size_t num_events) final {
		for (size_t event_idx = 0; event_idx < num_events; ++event_idx) {
			const viamd::Event& e = events[event_idx];

			switch (e.type) {
			case viamd::EventType_ViamdInitialize: {
				ApplicationState& state = *(ApplicationState*)e.payload;
				arena = md_arena_allocator_create(state.allocator.persistent, MEGABYTES(1));
				break;
			}
			case viamd::EventType_ViamdShutdown:
				md_arena_allocator_destroy(arena);
				break;
			case viamd::EventType_ViamdFrameTick:
				draw_window();
				break;
			case viamd::EventType_ViamdDrawMenu:
				ImGui::Checkbox("VeloxChem", &show_window);
				break;
			case viamd::EventType_ViamdPostRender: {
				ApplicationState& state = *(ApplicationState*)e.payload;
				draw_volume(state);
				break;
			}
			case viamd::EventType_ViamdTopologyInit: {
				ApplicationState& state = *(ApplicationState*)e.payload;
				str_t ext;
				str_t top_file = str_from_cstr(state.files.molecule);
				if (extract_ext(&ext, top_file) && str_eq_ignore_case(ext, STR_LIT("out"))) {
					MD_LOG_INFO("Attempting to load VeloxChem data from file '" STR_FMT "'", STR_ARG(top_file));
					md_vlx_data_free(&vlx);
					if (md_vlx_data_parse_file(&vlx, top_file, arena)) {
						MD_LOG_INFO("Successfully loaded VeloxChem data");
						show_window = true;

						size_t num_rows = vlx.scf.alpha.orbitals.dim[0];
						size_t num_cols = vlx.scf.alpha.orbitals.dim[1];

						for (int i = 0; i < (int)vlx.scf.alpha.occupations.count; ++i) {
							if (vlx.scf.alpha.occupations.data[i] > 0) {
								homo_idx = i;
							}
							else {
								lumo_idx = i;
								break;
							}
						}

						size_t num_pgtos = md_vlx_pgto_count(&vlx);
						size_t stride = ROUND_UP(num_pgtos, 16);
						size_t pgto_elem_bytes = sizeof(float) * 6 + sizeof(int) * 4;
						size_t tot_bytes = num_pgtos * pgto_elem_bytes;
						void* mem = md_arena_allocator_push_aligned(arena, tot_bytes, 64);
						MEMSET(mem, 0, tot_bytes);
						pgto.count = num_pgtos;
						pgto.x = (float*)mem;
						pgto.y = (float*)mem + stride * 1;
						pgto.z = (float*)mem + stride * 2;
						pgto.neg_alpha = (float*)mem + stride * 3;
						pgto.coeff = (float*)mem + stride * 4;
						pgto.cutoff = (float*)mem + stride * 5;
						pgto.i = (int*)mem + stride * 6;
						pgto.j = (int*)mem + stride * 7;
						pgto.k = (int*)mem + stride * 8;
						pgto.l = (int*)mem + stride * 9;

						vol_data = (float*)md_arena_allocator_push(arena, sizeof(float) * vol_dim * vol_dim * vol_dim);

						min_box = vec3_set1(FLT_MAX);
						max_box = vec3_set1(-FLT_MAX);

						for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
							vec3_t coord = vec3_set((float)vlx.geom.coord_x[i], (float)vlx.geom.coord_y[i], (float)vlx.geom.coord_z[i]);
							min_box = vec3_min(min_box, coord);
							max_box = vec3_max(max_box, coord);
						}

						const float pad = 2.0f;
						min_box = vec3_sub_f(min_box, pad);
						max_box = vec3_add_f(max_box, pad);
						step_size = vec3_div_f(vec3_sub(max_box, min_box), (float)vol_dim);

						model_mat = volume::compute_model_to_world_matrix(min_box, max_box);

						update_volume();
						show_volume = true;

						/*
						vec3_t coords[100];
						for (size_t i = 0; i < ARRAY_SIZE(coords); ++i) {
							double z = lerp(0.0, 2.0, i / (double)(ARRAY_SIZE(coords) - 1));
							r_values[i] = z;
							coords[i] = vec3_set(0,0, z * 0.529177210903);
						}
						for (int mo_idx = 0; mo_idx < 5; ++mo_idx) {
							for (size_t i = 0; i < num_mo_coeffs; ++i) {
								mo_coeffs[i] = vlx.scf.alpha.orbitals.data[i * num_cols + mo_idx];
							}
							md_vlx_get_mo(orb_values[mo_idx], coords, ARRAY_SIZE(coords), &vlx.geom, &vlx.basis, mo_coeffs, num_rows);
							const double s = sign(orb_values[mo_idx][10]);
							for (size_t i = 0; i < ARRAY_SIZE(coords); ++i) {
								orb_values[mo_idx][i] *= s;
							}
						}
						for (size_t i = 0; i < ARRAY_SIZE(r_values); ++i) {
							const double r2 = pow(r_values[i], 2);
							r_density[0][i] = 4 * PI * 2 * r2 *  pow(orb_values[0][i], 2);
							r_density[1][i] = 4 * PI * 2 * r2 *  pow(orb_values[1][i], 2);
							r_density[2][i] = 4 * PI * 2 * r2 * (pow(orb_values[2][i], 2) + pow(orb_values[3][i],2) + pow(orb_values[4][i],2));
						}
						*/
					}
					else {
						MD_LOG_INFO("Failed to load VeloxChem data");
						md_arena_allocator_reset(arena);
						vlx = {};
						pgto = {};
						mo_idx = 0;
						show_volume = false;
					}
				}
				break;
			}
			case viamd::EventType_ViamdTopologyFree:
				md_arena_allocator_reset(arena);
				vlx = {};
				pgto = {};
				mo_idx = 0;
				show_volume = false;
				break;
			default:
				break;
			}
		}
	}

	void update_volume() {
		if (task_system::task_is_running(compute_volume_task)) {
			task_system::task_interrupt(compute_volume_task);
		}
		else {
#define BLK_DIM 8

			if (!md_vlx_extract_alpha_mo_pgtos(&pgto, &vlx, mo_idx)) {
				MD_LOG_ERROR("Failed to extract alpha orbital for orbital index: %i", mo_idx);
				return;
			}

			uint32_t num_blocks = (vol_dim / BLK_DIM) * (vol_dim / BLK_DIM) * (vol_dim / BLK_DIM);
			// We evaluate the in parallel over smaller NxNxN blocks

			compute_volume_task = task_system::pool_enqueue(
				STR_LIT("Compute Volume"), 0, num_blocks,
				[](uint32_t range_beg, uint32_t range_end, void* user_data) {
					VeloxChem* data = (VeloxChem*)user_data;

					int num_blk = data->vol_dim / BLK_DIM;

					md_grid_t grid = {
						.data = data->vol_data,
						.dim = {data->vol_dim, data->vol_dim, data->vol_dim},
						.origin = {data->min_box.x + 0.5f * data->step_size.x, data->min_box.y + 0.5f * data->step_size.y,
								   data->min_box.z + 0.5f * data->step_size.z},
						.stepsize = {data->step_size.x, data->step_size.y, data->step_size.z},
					};

					// Conversion from �ngstr�m to Bohr
					const float factor = 1.0 / 0.529177210903;

					grid.origin[0] *= factor;
					grid.origin[1] *= factor;
					grid.origin[2] *= factor;

					grid.stepsize[0] *= factor;
					grid.stepsize[1] *= factor;
					grid.stepsize[2] *= factor;

					for (uint32_t i = range_beg; i < range_end; ++i) {
						// Determine block index
						int blk_z = i & (num_blk - 1);
						int blk_y = (i / num_blk) & (num_blk - 1);
						int blk_x = i / (num_blk * num_blk);

						const int beg_idx[3] = { blk_x * BLK_DIM, blk_y * BLK_DIM, blk_z * BLK_DIM };
						const int end_idx[3] = { blk_x * BLK_DIM + BLK_DIM, blk_y * BLK_DIM + BLK_DIM, blk_z * BLK_DIM + BLK_DIM };

						md_gto_grid_evaluate_sub(&grid, beg_idx, end_idx, &data->pgto);
					}
				},
				this);

#undef BLK_DIM

			task_system::main_enqueue(
				STR_LIT("Update Volume"),
				[](void* user_data) {
					VeloxChem* data = (VeloxChem*)user_data;
					gl::init_texture_3D(&data->vol_texture, data->vol_dim, data->vol_dim, data->vol_dim, GL_R32F);
					gl::set_texture_3D_data(data->vol_texture, data->vol_data, GL_R32F);
				},
				this, compute_volume_task);
		}
	}

	void draw_volume(ApplicationState& state) {
		if (!show_volume) return;

		volume::RenderDesc desc = {
			.render_target =
				{
					.depth = state.gbuffer.deferred.depth,
					.color = state.gbuffer.deferred.post_tonemap,
					.normal = state.gbuffer.deferred.normal,
					.width = state.gbuffer.width,
					.height = state.gbuffer.height,
				},
			.texture =
				{
					.volume = vol_texture,
					.transfer_function = tf_texture,
				},
			.matrix =
				{
					.model = model_mat,
					.view = state.view.param.matrix.current.view,
					.proj = state.view.param.matrix.current.proj_jittered,
				},
			.clip_volume =
				{
					.min = clip_min,
					.max = clip_max,
				},
			.global_scaling =
				{
					.density = density_scale,
				},
			.iso_surface =
				{
					.count = iso_surface.count,
					.values = iso_surface.values,
					.colors = iso_surface.colors,
				},
			.isosurface_enabled = iso_surface.enabled,
			.direct_volume_rendering_enabled = false,
			.voxel_spacing = step_size,
		};

		volume::render_volume(desc);

		if (bounding_box_enabled) {
			if (model_mat != mat4_t{ 0 }) {
				immediate::set_model_view_matrix(mat4_mul(state.view.param.matrix.current.view, model_mat));
				immediate::set_proj_matrix(state.view.param.matrix.current.proj);
				immediate::draw_box_wireframe(vec3_set1(0), vec3_set1(1), bounding_box_color);
				immediate::render();
			}
		}
	}

	// Draws the veloxchem stats
	void draw_veloxchem() {
		ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);  // This sets the window size the first time it is created
		// Begin is what creates the window
		if (ImGui::Begin("VeloxChem", &show_window)) {
			if (vlx.geom.num_atoms) {
				if (ImGui::TreeNode("Geometry")) {
					ImGui::Text("Num Atoms:           %6zu", vlx.geom.num_atoms);
					ImGui::Text("Num Alpha Electrons: %6zu", vlx.geom.num_alpha_electrons);
					ImGui::Text("Num Beta Electrons:  %6zu", vlx.geom.num_beta_electrons);
					ImGui::Text("Molecular Charge:    %6i", vlx.geom.molecular_charge);
					ImGui::Text("Spin Multiplicity:   %6i", vlx.geom.spin_multiplicity);
					ImGui::Spacing();
					ImGui::Text("Atom      Coord X      Coord Y      Coord Z");
					for (size_t i = 0; i < vlx.geom.num_atoms; ++i) {
						ImGui::Text("%4s %12.6f %12.6f %12.6f", vlx.geom.atom_symbol[i].buf, vlx.geom.coord_x[i], vlx.geom.coord_y[i],
							vlx.geom.coord_z[i]);
					}
					ImGui::TreePop();
				}
			}
			ImGui::ColorEdit4("Color Positive", iso_surface.colors[0].elem);
			ImGui::ColorEdit4("Color Negative", iso_surface.colors[1].elem);

			static double iso_val = 0.05f;
			const double iso_min = 1.0e-6;
			const double iso_max = 5.0;
			ImGui::SliderScalar("Iso Value", ImGuiDataType_Double, &iso_val, &iso_min, &iso_max, "%.6f", ImGuiSliderFlags_Logarithmic);
			iso_surface.values[0] = (float)iso_val;
			iso_surface.values[1] = -(float)iso_val;

			/*
			char buf[32];
			snprintf(buf, sizeof(buf), "%i", mo_idx);
			const ImVec2 draw_list_size = ImVec2(-1, 200);
			ImGui::Text("Molecular Orbitals");
			if (ImGui::BeginListBox("##Molecular Orbitals", draw_list_size)) {
				for (int n = 0; n < (int)vlx.scf.alpha.orbitals.dim[0]; ++n) {
					const char* lbl = "";
					if (homo_idx == n) {
						lbl = "HOMO";
					} else if (lumo_idx == n) {
						lbl = "LUMO";
					}
					const double energy = vlx.scf.alpha.energies.data[n];
					snprintf(buf, sizeof(buf), "%i %s", n + 1, lbl);
					if (ImGui::Selectable(buf, mo_idx == n)) {
						mo_idx = n;
						update_volume();
					}

					// Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
					if (mo_idx == n) { ImGui::SetItemDefaultFocus(); }
				}
				ImGui::EndListBox();
			}
			*/

			const float TEXT_BASE_HEIGHT = ImGui::GetTextLineHeightWithSpacing();

			enum {
				Col_Idx,
				Col_Occ,
				Col_Ene,
				Col_Lbl,
			};

			const ImGuiTableFlags flags = ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_RowBg |
				ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV | ImGuiTableFlags_NoBordersInBody |
				ImGuiTableFlags_ScrollY;
			if (ImGui::BeginTable("Molecular Orbitals", 4, flags, ImVec2(0.0f, TEXT_BASE_HEIGHT * 15), 0.0f)) {
				// Declare columns
				// We use the "user_id" parameter of TableSetupColumn() to specify a user id that will be stored in the sort specifications.
				// This is so our sort function can identify a column given our own identifier. We could also identify them based on their index!
				// Demonstrate using a mixture of flags among available sort-related flags:
				// - ImGuiTableColumnFlags_DefaultSort
				// - ImGuiTableColumnFlags_NoSort / ImGuiTableColumnFlags_NoSortAscending / ImGuiTableColumnFlags_NoSortDescending
				// - ImGuiTableColumnFlags_PreferSortAscending / ImGuiTableColumnFlags_PreferSortDescending
				ImGui::TableSetupColumn("Index", ImGuiTableColumnFlags_DefaultSort | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Idx);
				ImGui::TableSetupColumn("Occupation", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Occ);
				ImGui::TableSetupColumn("Energy", ImGuiTableColumnFlags_PreferSortDescending | ImGuiTableColumnFlags_WidthFixed, 0.0f, Col_Ene);
				ImGui::TableSetupColumn("Label", ImGuiTableColumnFlags_NoSort | ImGuiTableColumnFlags_WidthStretch, 0.0f, Col_Lbl);
				ImGui::TableSetupScrollFreeze(0, 1);  // Make row always visible
				ImGui::TableHeadersRow();

				// Demonstrate using clipper for large vertical lists
				// ImGuiListClipper clipper;
				// clipper.Begin(num_orbitals);
				// while (clipper.Step()) {
				for (int n = 0; n < num_orbitals(); n++) {
					// Display a data item
					ImGui::PushID(n + 1);
					ImGui::TableNextRow();
					bool is_selected = mo_idx == n;
					ImGui::TableNextColumn();
					char buf[32];
					snprintf(buf, sizeof(buf), "%i", n + 1);
					ImGuiSelectableFlags selectable_flags = ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap;
					if (ImGui::Selectable(buf, is_selected, selectable_flags)) {
						mo_idx = n;
						update_volume();
					}
					ImGui::TableNextColumn();
					ImGui::Text("%.1f", vlx.scf.alpha.occupations.data[n]);
					ImGui::TableNextColumn();
					ImGui::Text("%.4f", vlx.scf.alpha.energies.data[n]);
					ImGui::TableNextColumn();
					const char* lbl = (n == homo_idx) ? "HOMO" : (n == lumo_idx) ? "LUMO" : "";
					ImGui::TextUnformatted(lbl);
					ImGui::PopID();
				}
				//}
				ImGui::EndTable();
			}

			/*
			if (ImPlot::BeginPlot("MO")) {
				ImPlot::SetupAxes("Distance to nucleus (Bohr)", "Orbital value (a.u.)");
				ImPlot::PlotLine("1s",  r_values, orb_values[0], ARRAY_SIZE(r_values));
				ImPlot::PlotLine("2s",  r_values, orb_values[1], ARRAY_SIZE(r_values));
				ImPlot::PlotLine("2p1", r_values, orb_values[2], ARRAY_SIZE(r_values));
				ImPlot::PlotLine("2p2", r_values, orb_values[3], ARRAY_SIZE(r_values));
				ImPlot::PlotLine("2p3", r_values, orb_values[4], ARRAY_SIZE(r_values));
				ImPlot::EndPlot();
			}
			if (ImPlot::BeginPlot("Density")) {
				ImPlot::SetupAxes("Distance to nucleus (Bohr)", "Radial densities (a.u.)");
				ImPlot::PlotLine("1s", r_values, r_density[0], ARRAY_SIZE(r_values));
				ImPlot::PlotLine("2s", r_values, r_density[1], ARRAY_SIZE(r_values));
				ImPlot::PlotLine("2p", r_values, r_density[2], ARRAY_SIZE(r_values));
				ImPlot::EndPlot();
			}
			*/
		}
		ImGui::End();
	}

	void draw_scf() {
		// We set up iterations as doubles for easier use
		md_array(double) iter = md_array_create(double, vlx.scf.iter.count, arena);
		for (int i = 0; i < vlx.scf.iter.count; ++i) {
			iter[i] = (double)vlx.scf.iter.iteration[i];
		}

		// The actual plot
		ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
		if (ImGui::Begin("SCF", &show_window)) {
			// We draw 2 plots as "Energy total" has values in a different range then the rest of the data
			if (ImPlot::BeginSubplots("##AxisLinking", 2, 1, ImVec2(-1, -1), ImPlotSubplotFlags_LinkCols)) {
				if (ImPlot::BeginPlot("SCF")) {

					ImPlot::SetupAxisLimits(ImAxis_X1, 1.0, (int)vlx.scf.iter.count);
					ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);
					ImPlot::SetupAxes("Iterations", "eV");

					ImPlot::PlotLine("Density Change", iter, vlx.scf.iter.density_change, (int)vlx.scf.iter.count);
					ImPlot::PlotLine("Energy Change", iter, vlx.scf.iter.energy_change, (int)vlx.scf.iter.count);
					ImPlot::PlotLine("Gradient Norm", iter, vlx.scf.iter.gradient_norm, (int)vlx.scf.iter.count);
					ImPlot::PlotLine("Max Gradient", iter, vlx.scf.iter.max_gradient, (int)vlx.scf.iter.count);
				}
				ImPlot::EndPlot();

				if (ImPlot::BeginPlot("SCF")) {

					ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);

					ImPlot::PlotLine("Energy Total", iter, vlx.scf.iter.energy_total, (int)vlx.scf.iter.count);
				}
				ImPlot::EndPlot();
			}
			ImPlot::EndSubplots();
		}
		ImGui::End();
	}

	md_array(double) create_distributed_array(double start, double end, int count) {
		md_array(double) a = md_array_create(double, count, arena);
		double width = end - start;
		double inv_count = 1.0 / (double)count;

		for (int i = 0; i < count; i++) {
			a[i] = start + width * inv_count * (double)i;
		}
		return a;
	}

	md_array(double) gaussian_broadening(double* x_peaks, double* y_peaks, double sigma, double* x_array) {
		int num_x = (int)md_array_size(x_array);
		int num_e = (int)md_array_size(x_peaks);
		double tot = 0.0;
		md_array(double) y_return_array = md_array_create(double, num_x, arena);

		for (int xi = 0; xi < num_x; xi++) {
			tot = 0.0;
			for (int ei = 0; ei < num_e; ei++) {
				tot += y_peaks[ei] * exp(-(pow(x_peaks[ei] - x_array[xi], 2) / (2 * pow(sigma, 2))));
			}
			y_return_array[xi] = tot;
		}
		return y_return_array;
	}

	md_array(double) lorentzian_broadening(double* x_peaks, double* y_peaks, double sigma, double* x_array) {
		int num_x = (int)md_array_size(x_array);
		int num_e = (int)md_array_size(x_peaks);
		md_array(double) y_return_array = md_array_create(double, num_x, arena);
		double tot = 0.0;

		for (int xi = 0; xi < num_x; xi++) {
			//y_return_array[xi] = 0.0;
			tot = 0.0;
			for (int ei = 0; ei < num_e; ei++) {
				//y_return_array[xi] = y_return_array[xi] + y_peaks[ei] * (sigma / 2.0) / (pow(x_array[xi] - x_peaks[ei], 2) + pow(sigma / 2.0, 2));
				tot += y_peaks[ei] / (1 + pow((x_array[xi] - x_peaks[ei]) / sigma, 2));
				//tot += (y_peaks[ei] * sigma) / ((pow((x_array[xi] - x_peaks[ei]), 2) + pow(2 * sigma, 2)) / PI); // I don't 

			}
			y_return_array[xi] = tot;
		}
		return y_return_array;
	}
	/*
	double max(double* array) {
		for (int i = 0; i < )
	}
	*/

	md_array(double) scale_array(double* array, double factor) {
		md_array(double) return_array = md_array_create(double, md_array_size(array), arena);
		for (int i = 0; i < md_array_size(return_array); i++) {
			return_array[i] = array[i] * factor;
		}
		return return_array;
	}

	/*
	md_array(float) to_float(double array) {
		md_array(float) = md_array_create(float, md_array_size(array), arena)
		for (int i = 0; i < md_array_size(array); i++) {

		}
	}
	*/

	void draw_rsp() {

		// TODO: Implement broadening
		// Gaussian https://mdommett.github.io/blog/interpolation-with-gaussian-broadening/

		int count = (int)vlx.rsp.num_excited_states;
		static float sigma = 0.1;
		static int plot_points = 500;

		const char* broadening_options[] = { "Gaussian","Lorentzian" };
		static int broadening_current = 0;

		double con_fac = 1;
		const char* x_lable = "debug";
		const double eV2nm = 1239.84193;
		const char* convert_options[] = { "eV", "nm" };
		static int current_conversion = 1;

		// The actual plot
		ImGui::SetNextWindowSize({ 300, 350 }, ImGuiCond_FirstUseEver);
		if (ImGui::Begin("RSP", &show_window)) {
			// We draw 2 plots as "Energy total" has values in a different range then the rest of the data
			ImGui::DragFloat("Broadening", &sigma, 0.01, 0, 1);
			ImGui::DragInt("Num plot points", &plot_points, 1, 100, 1000);

			/*
			if (ImGui::BeginCombo("Broadening type", "Gaussian")) {
				ImGui::Selectable("Lor", true);
			}
			ImGui::EndCombo();
			*/

			ImGui::Combo("combo", &broadening_current, broadening_options, IM_ARRAYSIZE(broadening_options));
			//ImGui::ShowDemoWindow();

			ImGui::Combo("X lable", &current_conversion, convert_options, IM_ARRAYSIZE(convert_options));
			switch (current_conversion) {
			case 0:
				con_fac = 1;
				x_lable = "eV";
				break;
			case 1:
				con_fac = eV2nm;
				x_lable = "nm";
				break;
			}
			
			if (ImPlot::BeginPlot("Spectra")) {
				// ImPlot::SetupAxisLimits(ImAxis_X1, 1.0, vlx.scf.iter.count);
				ImPlot::SetupLegend(ImPlotLocation_East, ImPlotLegendFlags_Outside);
				ImPlot::SetupAxes(x_lable, "epsilon", ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);

				ImPlot::PlotBars("Exited States", scale_array(vlx.rsp.absorption_ev, con_fac), vlx.rsp.absorption_osc_str, count, 0.01 * con_fac);

				md_array(double) x_array = create_distributed_array(vlx.rsp.absorption_ev[0] - 1, vlx.rsp.absorption_ev[count - 1] + 1, plot_points);
				md_array(double) scaled_array = scale_array(x_array, con_fac);
				double test0 = scaled_array[0];
				double test1 = scaled_array[1];
				if (broadening_current == 0) {
					md_array(double) gaussian_array = gaussian_broadening(vlx.rsp.absorption_ev, vlx.rsp.absorption_osc_str, sigma, x_array);
					ImPlot::PlotLine("Gaussian", scaled_array, gaussian_array, plot_points);
				}
				else if (broadening_current == 1) {
					md_array(double) lorentzian_array = lorentzian_broadening(vlx.rsp.absorption_ev, vlx.rsp.absorption_osc_str, sigma, x_array);
					ImPlot::PlotLine("Lorentzian", scaled_array, lorentzian_array, plot_points);
				}
				//md_array(double) gaussian_arrayln2 = gaussian_broadening(vlx.rsp.absorption_ev, vlx.rsp.absorption_osc_str, sigma, x_array);
				//ImPlot::PlotLine("Gaussian ln2", x_array, gaussian_arrayln2, plot_points, 0, 10);

			}
			ImPlot::EndPlot();
			
		}
		ImGui::End();
	}

	void draw_window() {
		if (!show_window) return;
		draw_veloxchem();
		draw_scf();
		draw_rsp();
	}
};
static VeloxChem instance = {};
