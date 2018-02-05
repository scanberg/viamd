#pragma once

#include <mol/molecule.h>

struct GroStructure : MoleculeStructure {
    vec3 box_vectors[3];
};

// This is only valid for NVT types
