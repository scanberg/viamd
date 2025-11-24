/**
 * @file molden.cpp
 * @brief Implementation for Molden file format data structures
 * 
 * This file contains implementation details for the Molden data structures.
 * Phase 1 focuses on structure definitions and validation utilities.
 * 
 * Parsing implementation will be added in Phase 2.
 * 
 * Implementation Notes:
 * 
 * 1. Molden Format Overview:
 *    The Molden format is a text-based format for storing molecular orbital data.
 *    It consists of several sections, each starting with a tag in square brackets.
 * 
 * 2. Key Sections and Parsing Strategy:
 * 
 *    [Atoms] Section:
 *    - Contains atomic coordinates
 *    - Format: element_symbol atom_idx atomic_number x y z
 *    - Must handle both AU and Angs units
 *    - Parsing: Read line by line, split by whitespace, validate 6 fields
 * 
 *    [GTO] Section:
 *    - Contains Gaussian basis set information
 *    - Format is hierarchical:
 *        atom_index  0
 *        shell_type  num_primitives  [scale_factor]
 *        exponent1   coefficient1
 *        ...
 *        blank line (separator)
 *    - Parsing: State machine to track atom -> shell -> primitive levels
 *    - Edge case: SP shells have two coefficients per line
 * 
 *    [MO] Section:
 *    - Contains molecular orbital information
 *    - Format:
 *        Ene= energy
 *        Spin= Alpha|Beta
 *        Occup= occupation
 *        Sym= symmetry (optional)
 *        coefficient_1
 *        coefficient_2
 *        ...
 *    - Parsing: Read key-value pairs, then read coefficient list
 *    - Edge case: Number of coefficients must match total basis functions
 * 
 * 3. Edge Case Handling Strategy:
 * 
 *    a) Missing or Optional Fields:
 *       - Symmetry labels in MO section (use empty string)
 *       - Scale factors in GTO section (default to 1.0)
 *       - Atomic numbers (can infer from element symbol)
 * 
 *    b) SP Shells:
 *       - These require special handling as they have two coefficient columns
 *       - Need to expand into separate S and P shells for processing
 *       - Coefficient ordering: First for S, second for P
 * 
 *    c) Coordinate Units:
 *       - Must be detected from [Atoms] header
 *       - Convert to consistent internal representation
 *       - Default to Angstroms if not specified (common convention)
 * 
 *    d) Basis Function Ordering:
 *       - Cartesian: xx, yy, zz, xy, xz, yz (for d)
 *       - Spherical: Standard spherical harmonics ordering
 *       - Must track total count for MO coefficient validation
 * 
 *    e) Multiple Geometries:
 *       - Some files contain [GEOMETRIES] section with multiple structures
 *       - Phase 1: Store only first geometry
 *       - Future: Support trajectory-like data
 * 
 * 4. Validation Strategy:
 * 
 *    - Verify coefficient count matches basis function count
 *    - Check atom indices are positive and reasonable
 *    - Validate shell type characters
 *    - Ensure occupation numbers are in valid range [0.0, 2.0]
 *    - Check for duplicate atom indices
 *    - Verify exponents are positive
 * 
 * 5. Data Structure Mapping:
 * 
 *    Molden Text           ->  C++ Structure
 *    -----------               --------------
 *    [Atoms] section       ->  std::vector<Atom>
 *    Each atom line        ->  Atom struct
 *    [GTO] per atom        ->  AtomBasisSet
 *    Each shell block      ->  ContractedShell
 *    Each primitive line   ->  PrimitiveGaussian
 *    [MO] per orbital      ->  MolecularOrbital
 *    MO properties         ->  energy, spin, occupation, symmetry
 *    MO coefficients       ->  std::vector<double>
 * 
 * 6. Memory Management:
 * 
 *    - Use std::vector for dynamic arrays (automatic memory management)
 *    - Reserve capacity when size is known to avoid reallocations
 *    - Move semantics for efficient data transfer
 * 
 * 7. Integration Points (Future Phases):
 * 
 *    Phase 2 (Parsing):
 *    - Add md_molden_data_parse_str() function
 *    - Add md_molden_data_parse_file() function
 *    - Implement section parsers for each Molden section
 * 
 *    Phase 3 (Loader Integration):
 *    - Implement md_molden_system_loader() interface
 *    - Convert MoldenData to md_system_t structure
 *    - Register with VIAMD loader system
 *    - Add .molden extension to supported formats
 * 
 *    Phase 4 (Visualization):
 *    - Convert basis sets to md_gto_t format for rendering
 *    - Implement orbital visualization using existing GTO infrastructure
 *    - Add UI controls for orbital selection
 *    - Support isosurface rendering of orbitals
 */

#include "molden.h"
#include <algorithm>
#include <cmath>
#include <cctype>

namespace molden {

namespace util {

/**
 * @brief Validate MoldenData structure consistency
 * 
 * Checks:
 * - Atoms are present
 * - Basis sets reference valid atoms
 * - Orbital coefficient counts match basis function count
 * - Occupation numbers are valid
 * 
 * @param data The MoldenData to validate
 * @return true if valid, false otherwise
 * 
 * Future: Return detailed error messages via optional parameter
 */
bool validate_molden_data(const MoldenData& data) {
    // Check atoms exist
    if (data.atoms.empty()) {
        return false;
    }
    
    // Validate atom indices are positive
    for (const auto& atom : data.atoms) {
        if (atom.atom_index <= 0 || atom.atomic_number <= 0) {
            return false;
        }
    }
    
    // Validate basis sets reference valid atoms
    for (const auto& basis : data.basis_sets) {
        bool found = false;
        for (const auto& atom : data.atoms) {
            if (atom.atom_index == basis.atom_index) {
                found = true;
                break;
            }
        }
        if (!found) {
            return false;
        }
        
        // Validate shell data
        for (const auto& shell : basis.shells) {
            if (shell.primitives.empty()) {
                return false;
            }
            
            // Check exponents are positive
            for (const auto& prim : shell.primitives) {
                if (prim.exponent <= 0.0) {
                    return false;
                }
            }
        }
    }
    
    // Validate molecular orbitals
    for (const auto& mo : data.orbitals) {
        // Check occupation is in valid range [0.0, 2.0]
        if (mo.occupation < 0.0 || mo.occupation > 2.0) {
            return false;
        }
        
        // Check coefficient count matches total basis functions
        if (data.total_basis_functions > 0) {
            if (mo.coefficients.size() != data.total_basis_functions) {
                return false;
            }
        }
    }
    
    return true;
}

/**
 * @brief Calculate total number of basis functions in a MoldenData structure
 * 
 * Counts basis functions from all shells in all atom basis sets.
 * Handles special cases like SP shells and Cartesian vs Spherical formats.
 * 
 * @param data The MoldenData to analyze
 * @return Total number of basis functions
 */
size_t calculate_total_basis_functions(const MoldenData& data) {
    size_t total = 0;
    
    for (const auto& basis : data.basis_sets) {
        for (const auto& shell : basis.shells) {
            total += get_num_basis_functions(shell.shell_type, data.basis_format);
        }
    }
    
    return total;
}

/**
 * @brief Count number of orbitals by spin type
 * 
 * @param data The MoldenData to analyze
 * @param spin The spin type to count (Alpha or Beta)
 * @return Number of orbitals with specified spin
 */
size_t count_orbitals_by_spin(const MoldenData& data, SpinType spin) {
    size_t count = 0;
    
    for (const auto& mo : data.orbitals) {
        if (mo.spin == spin) {
            count++;
        }
    }
    
    return count;
}

/**
 * @brief Convert element symbol to atomic number
 * 
 * This is a fallback for files that omit atomic numbers.
 * Covers common elements found in quantum chemistry.
 * Uses case-insensitive comparison for robustness.
 * 
 * @param symbol Element symbol (e.g., "C", "H", "O")
 * @return Atomic number, or 0 if unknown
 */
int32_t element_symbol_to_atomic_number(const std::string& symbol) {
    // Helper for case-insensitive comparison
    auto to_upper = [](char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); };
    
    // Handle empty or invalid input
    if (symbol.empty() || symbol.length() > 2) {
        return 0;
    }
    
    // Single character elements
    if (symbol.length() == 1) {
        char c = to_upper(symbol[0]);
        switch (c) {
            case 'H': return 1;
            case 'B': return 5;
            case 'C': return 6;
            case 'N': return 7;
            case 'O': return 8;
            case 'F': return 9;
            case 'P': return 15;
            case 'S': return 16;
            case 'K': return 19;
            case 'V': return 23;
            default: return 0;
        }
    }
    
    // Two character elements - use a more efficient lookup
    // Pack two characters into a single uint16_t for fast comparison
    // Example: "He" -> (('H' << 8) | 'E') = 0x4845
    char c1 = to_upper(symbol[0]);
    char c2 = to_upper(symbol[1]);
    uint16_t key = (static_cast<uint16_t>(c1) << 8) | static_cast<uint16_t>(c2);
    
    switch (key) {
        case (('H' << 8) | 'E'): return 2;  // He
        case (('L' << 8) | 'I'): return 3;  // Li
        case (('B' << 8) | 'E'): return 4;  // Be
        case (('N' << 8) | 'E'): return 10; // Ne
        case (('N' << 8) | 'A'): return 11; // Na
        case (('M' << 8) | 'G'): return 12; // Mg
        case (('A' << 8) | 'L'): return 13; // Al
        case (('S' << 8) | 'I'): return 14; // Si
        case (('C' << 8) | 'L'): return 17; // Cl
        case (('A' << 8) | 'R'): return 18; // Ar
        case (('C' << 8) | 'A'): return 20; // Ca
        case (('S' << 8) | 'C'): return 21; // Sc
        case (('T' << 8) | 'I'): return 22; // Ti
        case (('C' << 8) | 'R'): return 24; // Cr
        case (('M' << 8) | 'N'): return 25; // Mn
        case (('F' << 8) | 'E'): return 26; // Fe
        case (('C' << 8) | 'O'): return 27; // Co
        case (('N' << 8) | 'I'): return 28; // Ni
        case (('C' << 8) | 'U'): return 29; // Cu
        case (('Z' << 8) | 'N'): return 30; // Zn
        default: return 0; // Unknown
    }
}

/**
 * @brief Parse shell type from string (handles 'sp' case)
 * 
 * Uses case-insensitive comparison without string allocation.
 * 
 * @param str Shell type string (e.g., "s", "p", "sp")
 * @return ShellType enum value
 */
ShellType parse_shell_type(const std::string& str) {
    if (str.empty() || str.length() > 2) {
        return ShellType::Unknown;
    }
    
    // Helper for case-insensitive character comparison
    auto to_lower = [](char c) { return static_cast<char>(std::tolower(static_cast<unsigned char>(c))); };
    
    // Single character shell types
    if (str.length() == 1) {
        char c = to_lower(str[0]);
        switch (c) {
            case 's': return ShellType::S;
            case 'p': return ShellType::P;
            case 'd': return ShellType::D;
            case 'f': return ShellType::F;
            case 'g': return ShellType::G;
            default: return ShellType::Unknown;
        }
    }
    
    // Two character - check for 'sp'
    if (to_lower(str[0]) == 's' && to_lower(str[1]) == 'p') {
        return ShellType::SP;
    }
    
    return ShellType::Unknown;
}

/**
 * @brief Convert spin string to SpinType enum
 * 
 * Uses case-insensitive comparison without string allocation.
 * 
 * @param str Spin string (e.g., "Alpha", "Beta")
 * @return SpinType enum value
 */
SpinType parse_spin_type(const std::string& str) {
    if (str.empty()) {
        return SpinType::Unknown;
    }
    
    // Helper for case-insensitive character comparison
    auto to_lower = [](char c) { return static_cast<char>(std::tolower(static_cast<unsigned char>(c))); };
    
    // Check first character to quickly differentiate
    char first = to_lower(str[0]);
    
    if (first == 'a' && str.length() >= 5) {
        // Check for "alpha"
        if (to_lower(str[1]) == 'l' && to_lower(str[2]) == 'p' && 
            to_lower(str[3]) == 'h' && to_lower(str[4]) == 'a') {
            return SpinType::Alpha;
        }
    } else if (first == 'b' && str.length() >= 4) {
        // Check for "beta"
        if (to_lower(str[1]) == 'e' && to_lower(str[2]) == 't' && 
            to_lower(str[3]) == 'a') {
            return SpinType::Beta;
        }
    }
    
    return SpinType::Unknown;
}

} // namespace util

} // namespace molden


// Future Phase 2 - Parsing Functions:
// 
// namespace molden {
// 
// /**
//  * @brief Parse Molden data from string buffer
//  * 
//  * Future implementation will:
//  * - Tokenize input by sections ([Atoms], [GTO], [MO])
//  * - Parse each section into corresponding data structures
//  * - Validate data consistency
//  * - Handle edge cases documented above
//  * 
//  * @param str Input string containing Molden file contents
//  * @return MoldenData structure with parsed data
//  */
// MoldenData parse_molden_string(const std::string& str);
// 
// /**
//  * @brief Parse Molden data from file
//  * 
//  * @param filename Path to Molden file
//  * @return MoldenData structure with parsed data
//  */
// MoldenData parse_molden_file(const std::string& filename);
// 
// } // namespace molden


// Future Phase 3 - Loader Integration:
//
// Integration with VIAMD loader system will require:
//
// 1. C wrapper functions for C API compatibility:
//    extern "C" {
//        md_system_loader_i* md_molden_system_loader(void);
//        bool md_molden_system_init(md_system_t* sys, const MoldenData* data, md_allocator_i* alloc);
//    }
//
// 2. Conversion functions:
//    - MoldenData -> md_system_t (atomic structure)
//    - MoldenData basis sets -> md_gto_data_t (for orbital visualization)
//
// 3. Registration in loader.cpp:
//    - Add SYS_LOADER_MOLDEN to sys_loader_t enum
//    - Add entry to sys_loader_name[] array
//    - Add ".molden;.mold" to sys_loader_ext[] array  
//    - Add md_molden_system_loader() to sys_loader[] array
