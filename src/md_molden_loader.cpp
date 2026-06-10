/**
 * @file md_molden_loader.cpp
 * @brief Molden file system loader for VIAMD integration
 * 
 * This file implements direct system loader entry points for Molden files,
 * converting MoldenData structures to md_system_t for visualization.
 * 
 * Phase 3 Implementation:
 * - Convert Molden atoms to md_system_t atom data
 * - Handle coordinate unit conversions (AU vs Angstrom)
 * - Infer bonds using VDW radii
 * - Provide loader interface for VIAMD
 */

#include "molden.h"

#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_str.h>
#include <md_system.h>
#include <md_util.h>

#include <cmath>

extern "C" {

/**
 * @brief Convert MoldenData to md_system_t
 * 
 * This function takes parsed Molden data and creates an md_system_t structure
 * suitable for visualization in VIAMD.
 * 
 * @param sys Output md_system_t structure
 * @param data Input MoldenData structure
 * @param alloc Memory allocator
 * @return true on success, false on error
 */
static bool molden_to_system(md_system_t* sys, const molden::MoldenData* data, md_allocator_i* alloc) {
    ASSERT(sys);
    ASSERT(data);
    ASSERT(alloc);
    
    if (data->atoms.empty()) {
        MD_LOG_ERROR("Molden to system: No atoms in Molden data");
        return false;
    }
    
    MEMSET(sys, 0, sizeof(md_system_t));
    sys->alloc = alloc;
    
    const size_t num_atoms = data->atoms.size();
    const size_t reserve_size = ALIGN_TO(num_atoms, 16);
    
    // Allocate arrays for atom data
    md_array_ensure(sys->atom.x, reserve_size, alloc);
    md_array_ensure(sys->atom.y, reserve_size, alloc);
    md_array_ensure(sys->atom.z, reserve_size, alloc);
    md_array_ensure(sys->atom.type_idx, reserve_size, alloc);
    md_array_ensure(sys->atom.flags, reserve_size, alloc);
    
    // Setup unknown atom type as fallback
    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unknown"), 0, 0.0f, 0.0f, 0, (md_flags_t)0, alloc);
    
    // Convert atoms
    for (size_t i = 0; i < num_atoms; ++i) {
        const molden::Atom& atom = data->atoms[i];
        
        // Convert coordinates (handle AU vs Angstrom)
        float x = atom.x;
        float y = atom.y;
        float z = atom.z;
        
        if (data->coord_unit == molden::CoordinateUnit::AtomicUnit) {
            // Convert from Bohr to Angstrom
            x = static_cast<float>(molden::util::au_to_angstrom(x));
            y = static_cast<float>(molden::util::au_to_angstrom(y));
            z = static_cast<float>(molden::util::au_to_angstrom(z));
        }
        
        // Get element symbol
        str_t element_symbol = str_from_cstr(atom.element_symbol.c_str());
        
        // Get atomic properties
        md_atomic_number_t atomic_number = (md_atomic_number_t)atom.atomic_number;
        float mass = md_atomic_number_mass(atomic_number);
        float radius = md_atomic_number_vdw_radius(atomic_number);
        
        // Find or add atom type
        md_atom_type_idx_t atom_type_idx = md_atom_type_find_or_add(
            &sys->atom.type, 
            element_symbol, 
            atomic_number, 
            mass, 
            radius, 
            0,
            (md_flags_t)0,
            alloc
        );
        
        // Add atom data
        sys->atom.count += 1;
        md_array_push(sys->atom.x, x, alloc);
        md_array_push(sys->atom.y, y, alloc);
        md_array_push(sys->atom.z, z, alloc);
        md_array_push(sys->atom.flags, (md_flags_t)0, alloc);
        md_array_push(sys->atom.type_idx, atom_type_idx, alloc);
    }

    // Seed a single component spanning the full molecule so mdlib can infer
    // entity/instance data during postprocess.
    md_array_push(sys->component.name, make_label(STR_LIT("MOL")), alloc);
    md_array_push(sys->component.seq_id, 1, alloc);
    md_array_push(sys->component.atom_offset, 0, alloc);
    md_array_push(sys->component.atom_offset, (uint32_t)num_atoms, alloc);
    md_array_push(sys->component.flags, (md_flags_t)0, alloc);
    sys->component.count = 1;
    
    // Infer bonds using VDW radii
    // This uses the mdlib utility function to compute covalent bonds
    md_util_system_infer_covalent_bonds(sys);
    
    MD_LOG_INFO("Loaded Molden file: %zu atoms, %zu bonds", sys->atom.count, sys->bond.count);
    
    return true;
}

/**
 * @brief Initialize md_system_t from Molden file content (string)
 */
bool md_molden_system_init_from_str(md_system_t* sys, str_t str, md_allocator_i* alloc) {
    std::string content(str.ptr, str.len);
    std::string error;
    
    molden::MoldenData data = molden::parse_molden_string(content, &error);
    
    if (data.atoms.empty()) {
        MD_LOG_ERROR("Failed to parse Molden string: %s", error.c_str());
        return false;
    }
    
    return molden_to_system(sys, &data, alloc);
}

/**
 * @brief Initialize md_system_t from Molden file
 */
bool md_molden_system_init_from_file(md_system_t* sys, str_t filename, md_allocator_i* alloc) {
    std::string filepath(filename.ptr, filename.len);
    std::string error;
    
    molden::MoldenData data = molden::parse_molden_file(filepath, &error);
    
    if (data.atoms.empty()) {
        MD_LOG_ERROR("Failed to parse Molden file '%.*s': %s", (int)filename.len, filename.ptr, error.c_str());
        return false;
    }
    
    return molden_to_system(sys, &data, alloc);
}

} // extern "C"
