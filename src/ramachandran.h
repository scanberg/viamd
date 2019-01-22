#pragma once

#include <core/vector_types.h>
#include <core/array_types.h>
#include <core/math_utils.h>
#include "image.h"
// enum class RamachandranConformationClassification { None, BetaHigh, BetaMid, BetaLow, AlphaHigh, AlphaMid, AlphaLow, LeftAlphaHigh, LeftAlphaMid, LeftAlphaLow, PMid, PLow };

namespace ramachandran {

void initialize();
void shutdown();

// Accumulation texture
void clear_accumulation_texture();
// Radius is given as percentage of normalized texture space coordinates (1.0 = 1% of texture width and height)
void compute_accumulation_texture(Array<const vec2> angles, vec4 color, float radius = 1.f, float outline = 0.f);

enum Region_ { Region_None, Region_AlphaHigh, Region_AlphaMid, Region_BetaHigh, Region_BetaMid, Region_LeftAlphaHigh, Region_LeftAlphaMid, Region_PMid };

constexpr float h_red = 0.0f / 360.0f;
constexpr float h_green = 120.0f / 360.0f;
constexpr float h_blue = 240.0f / 360.0f;
constexpr float h_yellow = 60.0f / 360.0f;
constexpr float c = 0.5f;
constexpr float l = 0.8f;
constexpr float a_high = 0.8f;
constexpr float a_mid = 0.2f;

const vec4 DEF_COL_BACKGROUND = vec4(1, 1, 1, 1);
const vec4 DEF_COL_ALPHA_HIGH = vec4(math::hcl_to_rgb(vec3(h_red, c, l)), a_high);
const vec4 DEF_COL_ALPHA_MID = vec4(math::hcl_to_rgb(vec3(h_red, c, l)), a_mid);
const vec4 DEF_COL_BETA_HIGH = vec4(math::hcl_to_rgb(vec3(h_blue, c, l)), a_high);
const vec4 DEF_COL_BETA_MID = vec4(math::hcl_to_rgb(vec3(h_blue, c, l)), a_mid);
const vec4 DEF_COL_LEFT_ALPHA_HIGH = vec4(math::hcl_to_rgb(vec3(h_green, c, l)), a_high);
const vec4 DEF_COL_LEFT_ALPHA_MID = vec4(math::hcl_to_rgb(vec3(h_green, c, l)), a_mid);
const vec4 DEF_COL_P_MID = vec4(math::hcl_to_rgb(vec3(h_yellow, c, l)), a_mid);

struct ColorMap {
    vec4 region_color[8] = {DEF_COL_BACKGROUND, DEF_COL_ALPHA_HIGH, DEF_COL_ALPHA_MID, DEF_COL_BETA_HIGH, DEF_COL_BETA_MID, DEF_COL_LEFT_ALPHA_HIGH, DEF_COL_LEFT_ALPHA_MID, DEF_COL_P_MID};
};

void init_gui_map(const ColorMap& color_map = {}, int blur_level = 2);
void init_segmentation_map(const ColorMap& color_map = {}, int blur_level = 2);
void init_color_map(const ColorMap& color_map = {}, int blur_level = 2);

const Image& get_gui_image();
const Image& get_segmentation_image();
const Image& get_color_image();

uint32 get_accumulation_texture();
uint32 get_gui_texture();
uint32 get_segmentation_texture();
uint32 get_color_texture();

}  // namespace ramachandran
