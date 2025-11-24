/**
 * @file molden.h
 * @brief Data structures for Molden file format representation
 * 
 * This header defines the C++ data structures for representing Molden file data.
 * The Molden format is used to store molecular geometry, basis functions, and 
 * molecular orbital information from quantum chemistry calculations.
 * 
 * Phase 1 Scope:
 * - Define data structures matching Molden format sections
 * - Document edge cases and format variations
 * - Provide clear mapping between Molden syntax and C++ types
 * 
 * Future Phases:
 * - Parsing implementation (Phase 2)
 * - Integration with VIAMD loader system (Phase 3)
 * - Visualization and UI components (Phase 4)
 * 
 * Reference: Molden format is primarily documented through examples from
 * quantum chemistry programs like GAMESS, Gaussian, ORCA, etc.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace molden {

/**
 * @brief Coordinate unit type for atomic positions
 * 
 * Molden files can specify coordinates in either Angstroms or atomic units (Bohr).
 * Format: [Atoms] (AU|Angs)
 */
enum class CoordinateUnit : uint8_t {
    Unknown,
    Angstrom,   // Standard Angstroms
    AtomicUnit  // Bohr (a.u.)
};

/**
 * @brief Spin type for molecular orbitals
 * 
 * Format in [MO] section: Spin= Alpha|Beta
 */
enum class SpinType : uint8_t {
    Unknown,
    Alpha,
    Beta
};

/**
 * @brief Gaussian shell type designation
 * 
 * Shell types in [GTO] section:
 * - s, p, d, f, g: standard angular momentum shells
 * - sp: combined S and P shell (requires special handling)
 * 
 * Edge case: 'sp' shells contain coefficients for both S and P in the same block
 */
enum class ShellType : uint8_t {
    Unknown,
    S,    // l = 0
    P,    // l = 1
    D,    // l = 2
    F,    // l = 3
    G,    // l = 4
    SP    // Combined S and P (special case)
};

/**
 * @brief Format for angular basis functions
 * 
 * Edge case: d, f, g orbitals can be represented in:
 * - Cartesian: xx, xy, xz, yy, yz, zz (6 for d, 10 for f, etc.)
 * - Spherical: using spherical harmonics (5 for d, 7 for f, etc.)
 * 
 * This affects the number of coefficients per orbital.
 */
enum class BasisFormat : uint8_t {
    Unknown,
    Cartesian,  // Cartesian GTOs (default in many programs)
    Spherical   // Spherical harmonic GTOs
};

/**
 * @brief Atomic data from [Atoms] section
 * 
 * Format: element_symbol  atom_index  atomic_number  x  y  z
 * Example: C   1   6   0.000000   0.000000   0.000000
 * 
 * Edge cases:
 * - atom_index may not be consecutive (gaps possible)
 * - Some files omit atomic_number (can be inferred from symbol)
 */
struct Atom {
    std::string element_symbol;  // Element symbol (e.g., "C", "H", "O")
    int32_t     atom_index;      // 1-based index from file
    int32_t     atomic_number;   // Atomic number (Z)
    float       x;               // X coordinate
    float       y;               // Y coordinate  
    float       z;               // Z coordinate
};

/**
 * @brief Primitive Gaussian function within a contracted shell
 * 
 * Each primitive is defined by an exponent (alpha) and contraction coefficient.
 * Format in [GTO] section:
 *   exponent  coefficient [coefficient_sp]
 * 
 * Edge case: SP shells have two coefficients per primitive
 */
struct PrimitiveGaussian {
    double exponent;       // Gaussian exponent (alpha)
    double coefficient;    // Contraction coefficient
    double coefficient_sp; // Second coefficient for SP shells (0.0 if not SP)
};

/**
 * @brief Contracted Gaussian shell
 * 
 * Format in [GTO] section:
 *   shell_type  num_primitives  scale_factor
 *   exponent1  coefficient1
 *   exponent2  coefficient2
 *   ...
 * 
 * A contracted shell is a linear combination of primitive Gaussians.
 */
struct ContractedShell {
    ShellType shell_type;                        // Type of shell (s, p, d, f, g, sp)
    std::vector<PrimitiveGaussian> primitives;   // Primitive Gaussians
    double scale_factor;                         // Optional scaling factor (default 1.0)
};

/**
 * @brief Basis set data for a single atom
 * 
 * Format in [GTO] section:
 *   atom_index  0
 *   [shell data]
 *   [shell data]
 *   blank line (separator)
 * 
 * Edge cases:
 * - The "0" after atom_index is always present but meaning unclear
 * - Multiple atoms may share the same basis set (not indicated in format)
 */
struct AtomBasisSet {
    int32_t atom_index;                      // 1-based atom index this basis applies to
    std::vector<ContractedShell> shells;     // Contracted shells for this atom
};

/**
 * @brief Molecular orbital data from [MO] section
 * 
 * Format:
 *   Ene= energy_value
 *   Spin= Alpha|Beta
 *   Occup= occupation_number
 *   Sym= symmetry_label
 *   coefficient_1
 *   coefficient_2
 *   ...
 * 
 * Edge cases:
 * - Sym field is optional (some programs omit it)
 * - Coefficients correspond to basis functions in order defined in [GTO]
 * - Number of coefficients must match total number of basis functions
 * - For Cartesian d: 6 functions, Spherical d: 5 functions
 */
struct MolecularOrbital {
    double energy;                       // Orbital energy (a.u.)
    SpinType spin;                       // Alpha or Beta
    double occupation;                   // Occupation number (0.0 to 2.0)
    std::string symmetry;                // Symmetry label (optional, may be empty)
    std::vector<double> coefficients;    // MO coefficients
};

/**
 * @brief Complete Molden file data
 * 
 * This structure represents all data that can be extracted from a Molden file.
 * 
 * Sections:
 * - [Molden Format]: Optional header
 * - [Title]: Optional title
 * - [Atoms]: Atomic coordinates (required)
 * - [GTO]: Basis set data (required for MO visualization)
 * - [MO]: Molecular orbital data (optional)
 * - [FREQ]: Vibrational frequencies (not implemented in Phase 1)
 * - [FR-COORD]: Fractional coordinates (not implemented in Phase 1)
 * - [GEOMETRIES]: Multiple geometries (not implemented in Phase 1)
 * 
 * Edge cases:
 * - Not all sections are always present
 * - Section order may vary
 * - Section names are case-insensitive in many parsers
 * - Some programs add custom sections
 */
struct MoldenData {
    // Metadata
    std::string title;                           // Optional title from [Title] section
    CoordinateUnit coord_unit;                   // Coordinate units (AU or Angs)
    BasisFormat basis_format;                    // Cartesian or Spherical (default Cartesian)
    
    // Atomic structure
    std::vector<Atom> atoms;                     // Atomic coordinates from [Atoms]
    
    // Basis set information
    std::vector<AtomBasisSet> basis_sets;        // Basis sets from [GTO]
    
    // Molecular orbitals
    std::vector<MolecularOrbital> orbitals;      // MOs from [MO]
    
    // Statistics for validation
    size_t total_basis_functions;                // Total number of basis functions
    size_t num_alpha_orbitals;                   // Number of alpha orbitals
    size_t num_beta_orbitals;                    // Number of beta orbitals
    
    /**
     * @brief Initialize default values
     */
    MoldenData() 
        : coord_unit(CoordinateUnit::Unknown)
        , basis_format(BasisFormat::Cartesian)
        , total_basis_functions(0)
        , num_alpha_orbitals(0)
        , num_beta_orbitals(0)
    {}
};

/**
 * @brief Helper functions for Molden data structures
 */
namespace util {

/**
 * @brief Convert coordinate from AU to Angstroms
 */
inline double au_to_angstrom(double au) {
    constexpr double BOHR_TO_ANGSTROM = 0.529177210903;
    return au * BOHR_TO_ANGSTROM;
}

/**
 * @brief Convert coordinate from Angstroms to AU
 */
inline double angstrom_to_au(double angstrom) {
    constexpr double ANGSTROM_TO_BOHR = 1.0 / 0.529177210903;
    return angstrom * ANGSTROM_TO_BOHR;
}

/**
 * @brief Get number of basis functions for a shell type
 * 
 * @param shell Shell type
 * @param format Basis format (Cartesian or Spherical)
 * @return Number of basis functions
 * 
 * Edge case handling:
 * - SP shells count as 1 + 3 = 4 functions
 * - Cartesian d has 6 functions (xx, yy, zz, xy, xz, yz)
 * - Spherical d has 5 functions (spherical harmonics)
 */
inline size_t get_num_basis_functions(ShellType shell, BasisFormat format = BasisFormat::Cartesian) {
    switch (shell) {
        case ShellType::S:  return 1;
        case ShellType::P:  return 3;
        case ShellType::SP: return 4;  // 1 (s) + 3 (p)
        case ShellType::D:  return (format == BasisFormat::Cartesian) ? 6 : 5;
        case ShellType::F:  return (format == BasisFormat::Cartesian) ? 10 : 7;
        case ShellType::G:  return (format == BasisFormat::Cartesian) ? 15 : 9;
        default: return 0;
    }
}

/**
 * @brief Convert shell type character to enum
 */
inline ShellType char_to_shell_type(char c) {
    switch (c) {
        case 's': case 'S': return ShellType::S;
        case 'p': case 'P': return ShellType::P;
        case 'd': case 'D': return ShellType::D;
        case 'f': case 'F': return ShellType::F;
        case 'g': case 'G': return ShellType::G;
        default: return ShellType::Unknown;
    }
}

/**
 * @brief Convert shell type enum to string
 */
inline const char* shell_type_to_string(ShellType shell) {
    switch (shell) {
        case ShellType::S:  return "s";
        case ShellType::P:  return "p";
        case ShellType::D:  return "d";
        case ShellType::F:  return "f";
        case ShellType::G:  return "g";
        case ShellType::SP: return "sp";
        default: return "unknown";
    }
}

/**
 * @brief Validate MoldenData structure consistency
 * 
 * Checks for valid atoms, basis sets, and molecular orbitals.
 * @param data The MoldenData to validate
 * @return true if valid, false otherwise
 */
bool validate_molden_data(const MoldenData& data);

/**
 * @brief Calculate total number of basis functions in a MoldenData structure
 * 
 * @param data The MoldenData to analyze
 * @return Total number of basis functions
 */
size_t calculate_total_basis_functions(const MoldenData& data);

/**
 * @brief Count number of orbitals by spin type
 * 
 * @param data The MoldenData to analyze
 * @param spin The spin type to count (Alpha or Beta)
 * @return Number of orbitals with specified spin
 */
size_t count_orbitals_by_spin(const MoldenData& data, SpinType spin);

/**
 * @brief Convert element symbol to atomic number
 * 
 * @param symbol Element symbol (e.g., "C", "H", "O")
 * @return Atomic number, or 0 if unknown
 */
int32_t element_symbol_to_atomic_number(const std::string& symbol);

/**
 * @brief Parse shell type from string (handles 'sp' case)
 * 
 * @param str Shell type string (e.g., "s", "p", "sp")
 * @return ShellType enum value
 */
ShellType parse_shell_type(const std::string& str);

/**
 * @brief Convert spin string to SpinType enum
 * 
 * @param str Spin string (e.g., "Alpha", "Beta")
 * @return SpinType enum value
 */
SpinType parse_spin_type(const std::string& str);

} // namespace util

/**
 * @brief Parse Molden data from string buffer
 * 
 * Parses a complete Molden file from a string buffer.
 * Extracts atoms, basis sets, and molecular orbitals.
 * 
 * @param str Input string containing Molden file contents
 * @param error_msg Optional pointer to string for error messages
 * @return MoldenData structure with parsed data, or empty structure on error
 * 
 * Example usage:
 * @code
 *   std::string error;
 *   MoldenData data = molden::parse_molden_string(file_contents, &error);
 *   if (data.atoms.empty()) {
 *       std::cerr << "Parse error: " << error << std::endl;
 *   }
 * @endcode
 */
MoldenData parse_molden_string(const std::string& str, std::string* error_msg = nullptr);

/**
 * @brief Parse Molden data from file
 * 
 * Reads and parses a Molden file from disk.
 * 
 * @param filename Path to Molden file
 * @param error_msg Optional pointer to string for error messages
 * @return MoldenData structure with parsed data, or empty structure on error
 * 
 * Example usage:
 * @code
 *   std::string error;
 *   MoldenData data = molden::parse_molden_file("molecule.molden", &error);
 *   if (!molden::util::validate_molden_data(data)) {
 *       std::cerr << "Invalid data: " << error << std::endl;
 *   }
 * @endcode
 */
MoldenData parse_molden_file(const std::string& filename, std::string* error_msg = nullptr);

} // namespace molden
