#pragma once

#if VIAMD_TREXIO

#include <cstdint>
#include <string>
#include <vector>

namespace trexio_data {

struct Atom {
    std::string label;
    int32_t atomic_number = 0;
    double charge = 0.0;
    double x = 0.0; // Bohr
    double y = 0.0; // Bohr
    double z = 0.0; // Bohr
};

struct PrimitiveGaussian {
    double exponent = 0.0;
    double coefficient = 0.0;
};

struct Shell {
    int32_t atom_index = 0;
    int32_t angular_momentum = 0;
    double factor = 1.0;
    std::vector<PrimitiveGaussian> primitives;
};

struct Orbital {
    double energy = 0.0;
    double occupation = 0.0;
    int32_t spin = 0; // TREXIO convention: 0 alpha/unspecified, 1 beta
    std::string symmetry;
    std::vector<double> coefficients;
};

struct TrexioData {
    std::string path;
    std::string basis_type;
    std::string mo_type;

    int32_t electron_num = 0;
    int32_t electron_up_num = 0;
    int32_t electron_dn_num = 0;
    int32_t ao_num = 0;
    int32_t mo_num = 0;
    int32_t cartesian = 0;

    std::vector<Atom> atoms;
    std::vector<Shell> shells;
    std::vector<int32_t> ao_shell;
    std::vector<Orbital> orbitals;
};

bool is_trexio_file(const char* filepath);
TrexioData parse_trexio_file(const std::string& filepath, std::string* error_msg = nullptr);

} // namespace trexio_data

#endif // VIAMD_TREXIO
