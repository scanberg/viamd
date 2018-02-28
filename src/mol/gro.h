#pragma once

#include <mol/molecule.h>

struct GroStructure : MoleculeInterface {
    mat3 box;
};

struct GroData : MoleculeData {
};

// This is only valid for NVT types
