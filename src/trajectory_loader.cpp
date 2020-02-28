#include "trajectory_loader.h"
#include <core/hash.h>
#include <mol/xtc_utils.h>
#include <mol/pdb_utils.h>
#include <mol/gro_utils.h>

struct MoleculeLoader {
    u64 id;
    bool (*read_molecule)(MoleculeStructure* mol, CStringView filename);
};

static constexpr MoleculeLoader mol_loaders[] {
    {hash::crc64("pdb"), pdb::load_molecule_from_file},
    {hash::crc64("gro"), gro::load_molecule_from_file},
};

struct TrajectoryLoader {
    u64 id;
    i32 (*read_num_frames)(CStringView filename);
    DynamicArray<i64> (*read_frame_offsets)(CStringView filename);
    //bool (*read_data)(Array<u8> dst, i64 offset, i64 size, CStringView filename);
    bool (*extract_frame)(TrajectoryFrame* frame, i32 num_atoms, Array<u8> data);
};

static constexpr TrajectoryLoader traj_loaders[] {
    {hash::crc64("xtc"), xtc::read_num_frames, xtc::read_frame_offsets, xtc::decompress_trajectory_frame},
    {hash::crc64("pdb"), pdb::read_num_frames, pdb::read_frame_offsets, pdb::extract_trajectory_frame},
};

namespace load {

namespace mol {
static MoleculeLoader* get_loader(CStringView extension) {
    StringBuffer<8> ext = extension;
    for (char* c : ext) {
        *c = to_lower(*c);
    }
    u64 id = hash::crc64(ext);
    for (int i = 0; i < ARRAY_SIZE(mol_loaders); i++) {
        if (id == mol_loaders[i].id) {
            return &mol_loaders[i];
        }
    }
    return NULL;
}

bool load_molecule(MoleculeStructure* mol, CStringView filename) {
    MoleculeLoader* loader = get_loader(get_file_extension(filename));
    if (loader) {
        return loader->read_molecule(mol, filename);
    }
    LOG_ERROR("Could not find loader for file '%.*s'", filename.length(), filename.beg());
    return false;
}

task::TaskID load_molecule_async(MoleculeStructure* mol, CStringView filename) {
    MoleculeLoader* loader = get_loader(get_file_extension(filename));
    if (loader) {
        task::TaskID id = task::create_task("Loading Molecule", [mol, CStringBuffer<256> file = filename, loader](task::TaskSetRange, task::TaskData){
            loader->read_molecule(mol, file);
        });
        return id;
    }
    LOG_ERROR("Could not find loader for file '%.*s'", filename.length(), filename.beg());
    return 0;
}

}  // namespace mol

namespace traj {
static TrajectoryLoader* get_loader(CStringView extension) {
    StringBuffer<8> ext = extension;
    for (char* c : ext) {
        *c = to_lower(*c);
    }
    u64 id = hash::crc64(ext);
    for (int i = 0; i < ARRAY_SIZE(traj_loaders); i++) {
        if (id == traj_loaders[i].id) {
            return &traj_loaders[i];
        }
    }
    return NULL;
}

bool load_trajectory(MoleculeTrajectory* traj, CStringView filename) {

}

task::TaskID load_trajectory_async(MoleculeTrajectory* traj, CStringView filename) {

}

}  // namespace traj

}  // namespace task