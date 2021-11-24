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
	uint32_t den_tex; // This uses 4 channels for each ramachandran type
	uint32_t map_tex[4];
	uint32_t iso_tex[4];
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
	const uint32_t* colors;
	uint32_t count;
};

bool rama_init(rama_data_t* data);
bool rama_free(rama_data_t* data);

task_system::ID rama_rep_compute_density(rama_rep_t* rep, const md_backbone_angles_t* angles, const uint32_t* rama_type_indices[4], uint32_t frame_beg, uint32_t frame_end, uint32_t frame_stride);

void rama_rep_render_map(rama_rep_t* rep, const float viewport[4], const rama_colormap_t colormap[4]);
void rama_rep_render_iso(rama_rep_t* rep, const float viewport[4], const rama_isomap_t isomap[4]);
