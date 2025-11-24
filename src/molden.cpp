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
#include <fstream>
#include <sstream>
#include <iostream>

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

// =============================================================================
// Phase 2 - Parsing Functions
// =============================================================================

namespace {

/**
 * @brief Trim whitespace from both ends of a string
 */
std::string trim(const std::string& str) {
    const char* whitespace = " \t\r\n";
    size_t start = str.find_first_not_of(whitespace);
    if (start == std::string::npos) {
        return "";
    }
    size_t end = str.find_last_not_of(whitespace);
    return str.substr(start, end - start + 1);
}

/**
 * @brief Split string by whitespace
 */
std::vector<std::string> split_whitespace(const std::string& str) {
    std::vector<std::string> tokens;
    std::string current;
    
    for (char c : str) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            if (!current.empty()) {
                tokens.push_back(current);
                current.clear();
            }
        } else {
            current += c;
        }
    }
    
    if (!current.empty()) {
        tokens.push_back(current);
    }
    
    return tokens;
}

/**
 * @brief Case-insensitive string comparison
 */
bool iequals(const std::string& a, const std::string& b) {
    if (a.length() != b.length()) {
        return false;
    }
    
    for (size_t i = 0; i < a.length(); ++i) {
        if (std::tolower(static_cast<unsigned char>(a[i])) != 
            std::tolower(static_cast<unsigned char>(b[i]))) {
            return false;
        }
    }
    
    return true;
}

/**
 * @brief Check if string starts with prefix (case-insensitive)
 */
bool istarts_with(const std::string& str, const std::string& prefix) {
    if (str.length() < prefix.length()) {
        return false;
    }
    
    for (size_t i = 0; i < prefix.length(); ++i) {
        if (std::tolower(static_cast<unsigned char>(str[i])) != 
            std::tolower(static_cast<unsigned char>(prefix[i]))) {
            return false;
        }
    }
    
    return true;
}

/**
 * @brief Replace all occurrences of 'D' or 'd' with 'E' for Fortran-style scientific notation
 */
std::string normalize_fortran_number(const std::string& str) {
    std::string result = str;
    for (char& c : result) {
        if (c == 'D' || c == 'd') {
            c = 'E';
        }
    }
    return result;
}

/**
 * @brief Parse a double value, handling Fortran-style scientific notation
 */
double parse_double(const std::string& str) {
    std::string normalized = normalize_fortran_number(str);
    return std::stod(normalized);
}

} // anonymous namespace

/**
 * @brief Parse [Atoms] section
 * 
 * Format: element_symbol atom_index atomic_number x y z
 * Header: [Atoms] (AU|Angs)
 */
bool parse_atoms_section(
    const std::vector<std::string>& lines,
    size_t& line_idx,
    MoldenData& data,
    std::string& error_msg
) {
    if (line_idx >= lines.size()) {
        error_msg = "Unexpected end of file in [Atoms] section";
        return false;
    }
    
    // Parse header to get coordinate units
    std::string header_line = trim(lines[line_idx]);
    
    // Check for coordinate unit specification
    if (header_line.find("AU") != std::string::npos || 
        header_line.find("au") != std::string::npos) {
        data.coord_unit = CoordinateUnit::AtomicUnit;
    } else if (header_line.find("Angs") != std::string::npos || 
               header_line.find("angs") != std::string::npos ||
               header_line.find("ANGS") != std::string::npos) {
        data.coord_unit = CoordinateUnit::Angstrom;
    } else {
        // Default to Angstroms if not specified
        data.coord_unit = CoordinateUnit::Angstrom;
    }
    
    line_idx++;
    
    // Parse atom lines
    while (line_idx < lines.size()) {
        std::string line = trim(lines[line_idx]);
        
        // Empty line or new section starts
        if (line.empty() || line[0] == '[') {
            break;
        }
        
        auto tokens = split_whitespace(line);
        if (tokens.size() < 5) {
            error_msg = "Invalid atom line (expected at least 5 fields): " + line;
            return false;
        }
        
        Atom atom;
        atom.element_symbol = tokens[0];
        
        try {
            atom.atom_index = std::stoi(tokens[1]);
            
            // Check if atomic number is present or needs to be inferred
            size_t coord_start = 2;
            if (tokens.size() >= 6) {
                // Assume format: symbol index Z x y z
                atom.atomic_number = std::stoi(tokens[2]);
                coord_start = 3;
            } else {
                // Infer atomic number from symbol
                atom.atomic_number = util::element_symbol_to_atomic_number(atom.element_symbol);
                coord_start = 2;
            }
            
            if (coord_start + 2 >= tokens.size()) {
                error_msg = "Not enough coordinate values in atom line: " + line;
                return false;
            }
            
            atom.x = static_cast<float>(parse_double(tokens[coord_start]));
            atom.y = static_cast<float>(parse_double(tokens[coord_start + 1]));
            atom.z = static_cast<float>(parse_double(tokens[coord_start + 2]));
            
        } catch (const std::exception& e) {
            error_msg = "Error parsing atom line: " + line + " (" + e.what() + ")";
            return false;
        }
        
        data.atoms.push_back(atom);
        line_idx++;
    }
    
    return true;
}

/**
 * @brief Parse [GTO] section
 * 
 * Format:
 *   atom_index 0
 *   shell_type num_primitives [scale_factor]
 *   exponent coefficient [coefficient_sp]
 *   ...
 *   blank line
 */
bool parse_gto_section(
    const std::vector<std::string>& lines,
    size_t& line_idx,
    MoldenData& data,
    std::string& error_msg
) {
    line_idx++; // Skip [GTO] header
    
    while (line_idx < lines.size()) {
        std::string line = trim(lines[line_idx]);
        
        // New section starts
        if (!line.empty() && line[0] == '[') {
            break;
        }
        
        // Skip blank lines
        if (line.empty()) {
            line_idx++;
            continue;
        }
        
        // Parse atom header: "atom_index 0"
        auto tokens = split_whitespace(line);
        if (tokens.size() < 2) {
            error_msg = "Invalid GTO atom header: " + line;
            return false;
        }
        
        AtomBasisSet basis;
        try {
            basis.atom_index = std::stoi(tokens[0]);
        } catch (const std::exception& e) {
            error_msg = "Error parsing atom index in GTO section: " + line;
            return false;
        }
        
        line_idx++;
        
        // Parse shells for this atom
        while (line_idx < lines.size()) {
            line = trim(lines[line_idx]);
            
            // Empty line signals end of this atom's basis
            if (line.empty()) {
                line_idx++;
                break;
            }
            
            // New section or new atom
            if (line[0] == '[' || (tokens = split_whitespace(line)).size() >= 2) {
                // Check if this is a new atom header (two integers)
                bool is_new_atom = false;
                if (tokens.size() >= 2) {
                    try {
                        std::stoi(tokens[0]);
                        std::stoi(tokens[1]);
                        is_new_atom = true;
                    } catch (...) {
                        is_new_atom = false;
                    }
                }
                
                if (line[0] == '[' || is_new_atom) {
                    break;
                }
            }
            
            // Parse shell header: "shell_type num_primitives [scale_factor]"
            tokens = split_whitespace(line);
            if (tokens.size() < 2) {
                error_msg = "Invalid shell header: " + line;
                return false;
            }
            
            ContractedShell shell;
            shell.shell_type = util::parse_shell_type(tokens[0]);
            
            if (shell.shell_type == ShellType::Unknown) {
                error_msg = "Unknown shell type: " + tokens[0];
                return false;
            }
            
            int num_primitives;
            try {
                num_primitives = std::stoi(tokens[1]);
                shell.scale_factor = (tokens.size() >= 3) ? parse_double(tokens[2]) : 1.0;
            } catch (const std::exception& e) {
                error_msg = "Error parsing shell header: " + line;
                return false;
            }
            
            line_idx++;
            
            // Parse primitives
            for (int i = 0; i < num_primitives; ++i) {
                if (line_idx >= lines.size()) {
                    error_msg = "Unexpected end of file while parsing primitives";
                    return false;
                }
                
                line = trim(lines[line_idx]);
                tokens = split_whitespace(line);
                
                if (tokens.size() < 2) {
                    error_msg = "Invalid primitive line: " + line;
                    return false;
                }
                
                PrimitiveGaussian prim;
                try {
                    prim.exponent = parse_double(tokens[0]);
                    prim.coefficient = parse_double(tokens[1]);
                    
                    // SP shells have two coefficients
                    if (shell.shell_type == ShellType::SP && tokens.size() >= 3) {
                        prim.coefficient_sp = parse_double(tokens[2]);
                    } else {
                        prim.coefficient_sp = 0.0;
                    }
                } catch (const std::exception& e) {
                    error_msg = "Error parsing primitive: " + line + " (" + e.what() + ")";
                    return false;
                }
                
                shell.primitives.push_back(prim);
                line_idx++;
            }
            
            basis.shells.push_back(shell);
        }
        
        data.basis_sets.push_back(basis);
    }
    
    return true;
}

/**
 * @brief Parse [MO] section
 * 
 * Format:
 *   Ene= energy
 *   Spin= Alpha|Beta
 *   Occup= occupation
 *   Sym= symmetry (optional)
 *   coefficient_1
 *   coefficient_2
 *   ...
 */
bool parse_mo_section(
    const std::vector<std::string>& lines,
    size_t& line_idx,
    MoldenData& data,
    std::string& error_msg
) {
    line_idx++; // Skip [MO] header
    
    while (line_idx < lines.size()) {
        std::string line = trim(lines[line_idx]);
        
        // New section starts
        if (!line.empty() && line[0] == '[') {
            break;
        }
        
        // Skip blank lines
        if (line.empty()) {
            line_idx++;
            continue;
        }
        
        // Start of a new MO - expect Ene= line
        if (!istarts_with(line, "Ene=") && !istarts_with(line, "ene=")) {
            line_idx++;
            continue;
        }
        
        MolecularOrbital mo;
        
        // Parse Ene= line
        size_t eq_pos = line.find('=');
        if (eq_pos == std::string::npos) {
            error_msg = "Invalid MO energy line: " + line;
            return false;
        }
        
        try {
            mo.energy = parse_double(trim(line.substr(eq_pos + 1)));
        } catch (const std::exception& e) {
            error_msg = "Error parsing MO energy: " + line;
            return false;
        }
        
        line_idx++;
        
        // Parse Spin= line
        if (line_idx >= lines.size()) {
            error_msg = "Unexpected end of file after Ene= line";
            return false;
        }
        
        line = trim(lines[line_idx]);
        if (!istarts_with(line, "Spin=") && !istarts_with(line, "spin=")) {
            error_msg = "Expected Spin= line, got: " + line;
            return false;
        }
        
        eq_pos = line.find('=');
        if (eq_pos != std::string::npos) {
            std::string spin_str = trim(line.substr(eq_pos + 1));
            mo.spin = util::parse_spin_type(spin_str);
            
            if (mo.spin == SpinType::Unknown) {
                error_msg = "Unknown spin type: " + spin_str;
                return false;
            }
        }
        
        line_idx++;
        
        // Parse Occup= line
        if (line_idx >= lines.size()) {
            error_msg = "Unexpected end of file after Spin= line";
            return false;
        }
        
        line = trim(lines[line_idx]);
        if (!istarts_with(line, "Occup=") && !istarts_with(line, "occup=")) {
            error_msg = "Expected Occup= line, got: " + line;
            return false;
        }
        
        eq_pos = line.find('=');
        if (eq_pos != std::string::npos) {
            try {
                mo.occupation = parse_double(trim(line.substr(eq_pos + 1)));
            } catch (const std::exception& e) {
                error_msg = "Error parsing occupation: " + line;
                return false;
            }
        }
        
        line_idx++;
        
        // Parse optional Sym= line
        if (line_idx < lines.size()) {
            line = trim(lines[line_idx]);
            if (istarts_with(line, "Sym=") || istarts_with(line, "sym=")) {
                eq_pos = line.find('=');
                if (eq_pos != std::string::npos) {
                    mo.symmetry = trim(line.substr(eq_pos + 1));
                }
                line_idx++;
            }
        }
        
        // Parse coefficients (continue until we hit empty line, new MO, or new section)
        while (line_idx < lines.size()) {
            line = trim(lines[line_idx]);
            
            // Check for end of this MO
            if (line.empty() || 
                line[0] == '[' || 
                istarts_with(line, "Ene=") || 
                istarts_with(line, "ene=")) {
                break;
            }
            
            // Line can contain multiple coefficients
            auto tokens = split_whitespace(line);
            
            // Skip lines that are just a single integer (orbital index)
            // This is a common variation where the MO index is printed
            if (tokens.size() == 1) {
                try {
                    int idx = std::stoi(tokens[0]);
                    // If it's a small positive integer, it's probably an index
                    if (idx > 0 && idx < 10000 && tokens[0].find('.') == std::string::npos) {
                        line_idx++;
                        continue;
                    }
                } catch (...) {
                    // Not an integer, continue parsing as coefficient
                }
            }
            
            for (const auto& token : tokens) {
                try {
                    double coeff = parse_double(token);
                    mo.coefficients.push_back(coeff);
                } catch (const std::exception& e) {
                    error_msg = "Error parsing MO coefficient: " + token;
                    return false;
                }
            }
            
            line_idx++;
        }
        
        data.orbitals.push_back(mo);
    }
    
    return true;
}

/**
 * @brief Parse Molden data from string buffer
 * 
 * @param str Input string containing Molden file contents
 * @param error_msg Output parameter for error messages
 * @return MoldenData structure with parsed data, or empty structure on error
 */
MoldenData parse_molden_string(const std::string& str, std::string* error_msg) {
    MoldenData data;
    std::string local_error;
    
    // Split into lines
    std::vector<std::string> lines;
    std::string current_line;
    for (char c : str) {
        if (c == '\n' || c == '\r') {
            if (!current_line.empty() || (c == '\n' && !lines.empty())) {
                lines.push_back(current_line);
                current_line.clear();
            }
        } else {
            current_line += c;
        }
    }
    if (!current_line.empty()) {
        lines.push_back(current_line);
    }
    
    // Parse sections
    size_t line_idx = 0;
    
    while (line_idx < lines.size()) {
        std::string line = trim(lines[line_idx]);
        
        // Skip empty lines and comments
        if (line.empty() || line[0] == '#' || line[0] == '!') {
            line_idx++;
            continue;
        }
        
        // Check for section headers
        if (line[0] == '[') {
            // Extract section name
            size_t end_bracket = line.find(']');
            if (end_bracket == std::string::npos) {
                local_error = "Invalid section header: " + line;
                if (error_msg) *error_msg = local_error;
                return MoldenData();
            }
            
            std::string section_name = trim(line.substr(1, end_bracket - 1));
            
            // Handle different sections
            if (iequals(section_name, "Title")) {
                // Parse title (next non-empty line)
                line_idx++;
                if (line_idx < lines.size()) {
                    data.title = trim(lines[line_idx]);
                    line_idx++;
                }
            }
            else if (iequals(section_name, "Atoms")) {
                if (!parse_atoms_section(lines, line_idx, data, local_error)) {
                    if (error_msg) *error_msg = local_error;
                    return MoldenData();
                }
            }
            else if (iequals(section_name, "GTO")) {
                if (!parse_gto_section(lines, line_idx, data, local_error)) {
                    if (error_msg) *error_msg = local_error;
                    return MoldenData();
                }
            }
            else if (iequals(section_name, "MO")) {
                if (!parse_mo_section(lines, line_idx, data, local_error)) {
                    if (error_msg) *error_msg = local_error;
                    return MoldenData();
                }
            }
            else if (iequals(section_name, "5D")) {
                data.basis_format = BasisFormat::Spherical;
                line_idx++;
            }
            else if (iequals(section_name, "5D7F") || iequals(section_name, "5D10F")) {
                data.basis_format = BasisFormat::Spherical;
                line_idx++;
            }
            else {
                // Unknown section - skip it
                line_idx++;
            }
        } else {
            line_idx++;
        }
    }
    
    // Calculate statistics
    data.total_basis_functions = util::calculate_total_basis_functions(data);
    data.num_alpha_orbitals = util::count_orbitals_by_spin(data, SpinType::Alpha);
    data.num_beta_orbitals = util::count_orbitals_by_spin(data, SpinType::Beta);
    
    return data;
}

/**
 * @brief Parse Molden data from file
 * 
 * @param filename Path to Molden file
 * @param error_msg Output parameter for error messages
 * @return MoldenData structure with parsed data, or empty structure on error
 */
MoldenData parse_molden_file(const std::string& filename, std::string* error_msg) {
    // Read file into string
    std::ifstream file(filename);
    if (!file.is_open()) {
        if (error_msg) {
            *error_msg = "Failed to open file: " + filename;
        }
        return MoldenData();
    }
    
    std::string content((std::istreambuf_iterator<char>(file)),
                        std::istreambuf_iterator<char>());
    file.close();
    
    return parse_molden_string(content, error_msg);
}

/**
 * @brief Main loader function for VIAMD integration
 * 
 * This is a convenience wrapper around parse_molden_file that provides
 * a simple C-compatible interface for loading Molden files.
 * 
 * @param filepath Path to Molden file
 * @return MoldenData structure with parsed data, or empty structure on error
 */
MoldenData load_molden_file(const char* filepath) {
    if (!filepath) {
        return MoldenData();
    }
    
    std::string error;
    MoldenData data = parse_molden_file(std::string(filepath), &error);
    
    if (!error.empty()) {
        // Log error to stderr for debugging
        std::cerr << "Molden loader error: " << error << std::endl;
    }
    
    return data;
}

} // namespace molden


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
