#pragma once

#include <mol/molecule.h>

struct GroStructure : MoleculeStructure {
    mat3 box;
};

// This is only valid for NVT types
