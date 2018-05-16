#pragma once

#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>

struct MoleculeDynamic {
	MoleculeStructure molecule{};
	MoleculeTrajectory trajectory{};
};