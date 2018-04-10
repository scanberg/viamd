#pragma once

#include <mol/molecule.h>
#include <mol/trajectory.h>

struct MoleculeDynamic {
	MoleculeStructure molecule{};
	Trajectory trajectory{};
};