#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_molecule.h>
#include <task_system.h>

namespace ramachandran {
	void initialize();
	void shutdown();
};

// Representation (density, colormap, iso)
struct rama_rep_t {
	uint32_t map_tex[4];
	uint32_t iso_tex[4];
	float    den_sum[4];	// Density sum for each ramachandran type (General, Glycine, Proline and PreProline)
	uint32_t den_tex;		// This uses 4 channels, one for each ramachandran type (General, Glycine, Proline and PreProline)
};

struct rama_data_t {
	rama_rep_t ref;
	rama_rep_t full;
	rama_rep_t filt;
};

struct rama_colormap_t {
	const uint32_t* colors;
	uint32_t count;
	float min_value;
	float max_value;
};

struct rama_isomap_t {
	const float* values;
	const uint32_t* level_colors;	// Optional
	const uint32_t* contour_colors; // Optional
	uint32_t count;
};

bool rama_init(rama_data_t* data);
bool rama_free(rama_data_t* data);

task_system::ID rama_rep_compute_density(rama_rep_t* rep, const md_backbone_angles_t* angles, const uint32_t* rama_type_indices[4], uint32_t frame_beg, uint32_t frame_end, uint32_t frame_stride);

// Computes the iso levels given a set of percentiles, e.g. (0.85) will compute which (density) value best corresponds to that
//bool rama_rep_compute_density_levels(float* out_levels[4], const rama_rep_t* rep, const float* percentiles, int64_t num_percentiles);

// Display resolution is the max dim of the displayed texture resolution, this is used as a hint to get the renderings to look better (anti-aliased lines etc dependent on the effective view resolution)
void rama_rep_render_map(rama_rep_t* rep, const float viewport[4], const rama_colormap_t colormap[4], uint32_t display_res);
void rama_rep_render_iso(rama_rep_t* rep, const float viewport[4], const rama_isomap_t isomap[4], uint32_t display_res);
