#pragma once

#include <viamd.h>

namespace openmm_interface {
    // Public interface for other components to trigger energy minimization
    void minimize_energy_if_available(ApplicationState& state);
}