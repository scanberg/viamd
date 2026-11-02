#pragma once

#include "nci_surface_model.h"

namespace topology {

void draw_nci_panel(NciState& state, bool critical_points_enabled, std::vector<uint32_t>& critical_point_list_scope,
                    std::vector<uint32_t>& critical_point_environment_selection);

}  // namespace topology
