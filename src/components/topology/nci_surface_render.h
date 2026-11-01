#pragma once

#include "nci_surface_model.h"

struct MeshSurfaceRequest;
struct RenderKey;

namespace topology {

void draw_nci_legends(const NciState& state, int first_slot);
void nci_append_mesh_surfaces(NciState& state, MeshSurfaceRequest& request);
void nci_hash_render_state(const NciState& state, RenderKey& key);

}  // namespace topology
