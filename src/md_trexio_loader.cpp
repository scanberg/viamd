#if VIAMD_TREXIO

#include "md_trexio_loader.h"
#include "trexio_data.h"

#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_str.h>
#include <md_system.h>
#include <md_util.h>

#include <string>

#define TREXIO_BOHR_TO_ANGSTROM 0.529177210903

extern "C" {

bool md_trexio_path_is_trexio(str_t filename) {
    std::string filepath(filename.ptr, filename.len);
    return trexio_data::is_trexio_file(filepath.c_str());
}

static bool trexio_to_system(md_system_t* sys, const trexio_data::TrexioData* data, md_allocator_i* alloc) {
    ASSERT(sys);
    ASSERT(data);
    ASSERT(alloc);

    if (data->atoms.empty()) {
        MD_LOG_ERROR("TREXIO to system: no nuclei in TREXIO data");
        return false;
    }

    MEMSET(sys, 0, sizeof(md_system_t));
    sys->alloc = alloc;

    const size_t num_atoms = data->atoms.size();
    const size_t reserve_size = ALIGN_TO(num_atoms, 16);
    md_array_ensure(sys->atom.x, reserve_size, alloc);
    md_array_ensure(sys->atom.y, reserve_size, alloc);
    md_array_ensure(sys->atom.z, reserve_size, alloc);
    md_array_ensure(sys->atom.type_idx, reserve_size, alloc);
    md_array_ensure(sys->atom.flags, reserve_size, alloc);

    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unknown"), 0, 0.0f, 0.0f, 0, (md_flags_t)0, alloc);

    for (size_t i = 0; i < num_atoms; ++i) {
        const trexio_data::Atom& atom = data->atoms[i];
        md_atomic_number_t z = (md_atomic_number_t)CLAMP(atom.atomic_number, 0, 118);
        str_t element_symbol = md_atomic_number_symbol(z);
        if (!element_symbol.ptr || element_symbol.len == 0) {
            element_symbol = STR_LIT("Unknown");
        }

        md_atom_type_idx_t atom_type_idx = md_atom_type_find_or_add(
            &sys->atom.type,
            element_symbol,
            z,
            md_atomic_number_mass(z),
            md_atomic_number_vdw_radius(z),
            0,
            (md_flags_t)0,
            alloc);

        sys->atom.count += 1;
        md_array_push(sys->atom.x, (float)(atom.x * TREXIO_BOHR_TO_ANGSTROM), alloc);
        md_array_push(sys->atom.y, (float)(atom.y * TREXIO_BOHR_TO_ANGSTROM), alloc);
        md_array_push(sys->atom.z, (float)(atom.z * TREXIO_BOHR_TO_ANGSTROM), alloc);
        md_array_push(sys->atom.flags, (md_flags_t)0, alloc);
        md_array_push(sys->atom.type_idx, atom_type_idx, alloc);
    }

    md_array_push(sys->component.name, make_label(STR_LIT("TREXIO")), alloc);
    md_array_push(sys->component.seq_id, 1, alloc);
    md_array_push(sys->component.atom_offset, 0, alloc);
    md_array_push(sys->component.atom_offset, (uint32_t)num_atoms, alloc);
    md_array_push(sys->component.flags, (md_flags_t)0, alloc);
    sys->component.count = 1;

    md_util_system_infer_covalent_bonds(sys);
    MD_LOG_INFO("Loaded TREXIO file: %zu atoms, %zu bonds", sys->atom.count, sys->bond.count);
    return true;
}

bool md_trexio_system_init_from_file(md_system_t* sys, str_t filename, md_allocator_i* alloc) {
    std::string filepath(filename.ptr, filename.len);
    std::string error;
    trexio_data::TrexioData data = trexio_data::parse_trexio_file(filepath, &error);
    if (data.atoms.empty()) {
        MD_LOG_ERROR("Failed to parse TREXIO file '%.*s': %s", (int)filename.len, filename.ptr, error.c_str());
        return false;
    }
    return trexio_to_system(sys, &data, alloc);
}

} // extern "C"

#endif // VIAMD_TREXIO
