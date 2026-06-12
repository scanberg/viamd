#if VIAMD_TREXIO

#include "trexio_data.h"

#include <trexio.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <memory>

namespace trexio_data {

static const char* trexio_error_string(trexio_exit_code rc) {
    return trexio_string_of_error(rc);
}

static bool trexio_call_ok(trexio_exit_code rc, std::string* error_msg, const char* what) {
    if (rc == TREXIO_SUCCESS) return true;
    if (error_msg) {
        char buf[256];
        snprintf(buf, sizeof(buf), "%s: %s", what, trexio_error_string(rc));
        *error_msg = buf;
    }
    return false;
}

static trexio_t* open_readonly(const char* filepath, std::string* error_msg) {
    trexio_exit_code rc = TREXIO_SUCCESS;
    trexio_t* file = trexio_open(filepath, 'r', TREXIO_AUTO, &rc);
    if (file) return file;

    file = trexio_open(filepath, 'r', TREXIO_HDF5, &rc);
    if (file) return file;

    file = trexio_open(filepath, 'r', TREXIO_TEXT, &rc);
    if (file) return file;

    if (error_msg) {
        char buf[256];
        snprintf(buf, sizeof(buf), "Failed to open TREXIO file: %s", trexio_error_string(rc));
        *error_msg = buf;
    }
    return nullptr;
}

static std::string read_string_attr(trexio_t* file, trexio_exit_code (*has_fn)(trexio_t* const), trexio_exit_code (*read_fn)(trexio_t* const, char* const, const int32_t), int32_t max_len = 128) {
    if (!file || has_fn(file) != TREXIO_SUCCESS) return {};
    std::string str;
    str.resize((size_t)max_len, '\0');
    if (read_fn(file, str.data(), max_len) != TREXIO_SUCCESS) return {};
    str.resize(str.find('\0') == std::string::npos ? str.size() : str.find('\0'));
    return str;
}

static std::vector<std::string> read_string_array(trexio_t* file, int32_t count, trexio_exit_code (*has_fn)(trexio_t* const), trexio_exit_code (*read_fn)(trexio_t* const, char**, const int32_t), int32_t max_len = 64) {
    std::vector<std::string> result;
    if (!file || count <= 0 || has_fn(file) != TREXIO_SUCCESS) return result;

    result.resize((size_t)count);
    std::vector<std::unique_ptr<char[]>> storage((size_t)count);
    std::vector<char*> ptrs((size_t)count);
    for (int32_t i = 0; i < count; ++i) {
        storage[(size_t)i] = std::make_unique<char[]>((size_t)max_len);
        storage[(size_t)i][0] = '\0';
        ptrs[(size_t)i] = storage[(size_t)i].get();
    }

    if (read_fn(file, ptrs.data(), max_len) != TREXIO_SUCCESS) {
        result.clear();
        return result;
    }

    for (int32_t i = 0; i < count; ++i) {
        result[(size_t)i] = ptrs[(size_t)i] ? ptrs[(size_t)i] : "";
    }
    return result;
}

bool is_trexio_file(const char* filepath) {
    if (!filepath || !filepath[0]) return false;
    trexio_exit_code rc = TREXIO_SUCCESS;
    trexio_t* file = trexio_open(filepath, 'r', TREXIO_AUTO, &rc);
    if (!file) file = trexio_open(filepath, 'r', TREXIO_HDF5, &rc);
    if (!file) file = trexio_open(filepath, 'r', TREXIO_TEXT, &rc);
    if (!file) return false;

    const bool ok = trexio_has_nucleus_num(file) == TREXIO_SUCCESS;
    trexio_close(file);
    return ok;
}

TrexioData parse_trexio_file(const std::string& filepath, std::string* error_msg) {
    TrexioData data;
    data.path = filepath;

    trexio_t* file = open_readonly(filepath.c_str(), error_msg);
    if (!file) return data;

    auto close_file = [&]() { trexio_close(file); file = nullptr; };

    int32_t nucleus_num = 0;
    if (!trexio_call_ok(trexio_read_nucleus_num(file, &nucleus_num), error_msg, "read nucleus.num")) {
        close_file();
        return {};
    }
    if (nucleus_num <= 0) {
        if (error_msg) *error_msg = "TREXIO file has no nuclei";
        close_file();
        return {};
    }

    std::vector<double> charges((size_t)nucleus_num, 0.0);
    std::vector<double> coords((size_t)nucleus_num * 3, 0.0);
    if (!trexio_call_ok(trexio_read_nucleus_charge(file, charges.data()), error_msg, "read nucleus.charge") ||
        !trexio_call_ok(trexio_read_nucleus_coord(file, coords.data()), error_msg, "read nucleus.coord")) {
        close_file();
        return {};
    }

    std::vector<std::string> labels = read_string_array(file, nucleus_num, trexio_has_nucleus_label, trexio_read_nucleus_label);
    data.atoms.resize((size_t)nucleus_num);
    for (int32_t i = 0; i < nucleus_num; ++i) {
        Atom& atom = data.atoms[(size_t)i];
        atom.label = (i < (int32_t)labels.size()) ? labels[(size_t)i] : std::string();
        atom.charge = charges[(size_t)i];
        atom.atomic_number = (int32_t)std::lround(atom.charge);
        atom.x = coords[(size_t)i * 3 + 0];
        atom.y = coords[(size_t)i * 3 + 1];
        atom.z = coords[(size_t)i * 3 + 2];
    }

    if (trexio_has_electron_num(file) == TREXIO_SUCCESS) trexio_read_electron_num(file, &data.electron_num);
    if (trexio_has_electron_up_num(file) == TREXIO_SUCCESS) trexio_read_electron_up_num(file, &data.electron_up_num);
    if (trexio_has_electron_dn_num(file) == TREXIO_SUCCESS) trexio_read_electron_dn_num(file, &data.electron_dn_num);

    data.basis_type = read_string_attr(file, trexio_has_basis_type, trexio_read_basis_type);
    data.mo_type = read_string_attr(file, trexio_has_mo_type, trexio_read_mo_type);

    int32_t shell_num = 0;
    int32_t prim_num = 0;
    if (trexio_has_basis_shell_num(file) == TREXIO_SUCCESS &&
        trexio_has_basis_prim_num(file) == TREXIO_SUCCESS &&
        trexio_read_basis_shell_num(file, &shell_num) == TREXIO_SUCCESS &&
        trexio_read_basis_prim_num(file, &prim_num) == TREXIO_SUCCESS &&
        shell_num > 0 && prim_num > 0) {

        std::vector<int32_t> nucleus_index((size_t)shell_num, 0);
        std::vector<int32_t> shell_ang_mom((size_t)shell_num, 0);
        std::vector<double> shell_factor((size_t)shell_num, 1.0);
        std::vector<int32_t> shell_index((size_t)prim_num, 0);
        std::vector<double> exponent((size_t)prim_num, 0.0);
        std::vector<double> coefficient((size_t)prim_num, 0.0);
        std::vector<double> prim_factor((size_t)prim_num, 1.0);

        bool basis_ok = true;
        basis_ok &= trexio_read_basis_nucleus_index(file, nucleus_index.data()) == TREXIO_SUCCESS;
        basis_ok &= trexio_read_basis_shell_ang_mom(file, shell_ang_mom.data()) == TREXIO_SUCCESS;
        basis_ok &= trexio_read_basis_shell_index(file, shell_index.data()) == TREXIO_SUCCESS;
        basis_ok &= trexio_read_basis_exponent(file, exponent.data()) == TREXIO_SUCCESS;
        basis_ok &= trexio_read_basis_coefficient(file, coefficient.data()) == TREXIO_SUCCESS;
        if (trexio_has_basis_shell_factor(file) == TREXIO_SUCCESS) {
            basis_ok &= trexio_read_basis_shell_factor(file, shell_factor.data()) == TREXIO_SUCCESS;
        }
        if (trexio_has_basis_prim_factor(file) == TREXIO_SUCCESS) {
            basis_ok &= trexio_read_basis_prim_factor(file, prim_factor.data()) == TREXIO_SUCCESS;
        }

        if (basis_ok) {
            data.shells.resize((size_t)shell_num);
            for (int32_t s = 0; s < shell_num; ++s) {
                Shell& shell = data.shells[(size_t)s];
                shell.atom_index = nucleus_index[(size_t)s];
                shell.angular_momentum = shell_ang_mom[(size_t)s];
                shell.factor = shell_factor[(size_t)s];
            }
            for (int32_t p = 0; p < prim_num; ++p) {
                const int32_t s = shell_index[(size_t)p];
                if (s < 0 || s >= shell_num) continue;
                PrimitiveGaussian prim = {
                    .exponent = exponent[(size_t)p],
                    .coefficient = coefficient[(size_t)p] * prim_factor[(size_t)p] * shell_factor[(size_t)s],
                };
                data.shells[(size_t)s].primitives.push_back(prim);
            }
        }
    }

    if (trexio_has_ao_cartesian(file) == TREXIO_SUCCESS) trexio_read_ao_cartesian(file, &data.cartesian);
    if (trexio_has_ao_num(file) == TREXIO_SUCCESS) trexio_read_ao_num(file, &data.ao_num);
    if (data.ao_num > 0 && trexio_has_ao_shell(file) == TREXIO_SUCCESS) {
        data.ao_shell.resize((size_t)data.ao_num);
        if (trexio_read_ao_shell(file, data.ao_shell.data()) != TREXIO_SUCCESS) {
            data.ao_shell.clear();
        }
    }

    if (data.ao_num > 0 && trexio_has_mo_num(file) == TREXIO_SUCCESS && trexio_has_mo_coefficient(file) == TREXIO_SUCCESS) {
        trexio_read_mo_num(file, &data.mo_num);
        if (data.mo_num > 0) {
            std::vector<double> coefficients((size_t)data.mo_num * (size_t)data.ao_num, 0.0);
            if (trexio_read_mo_coefficient(file, coefficients.data()) == TREXIO_SUCCESS) {
                data.orbitals.resize((size_t)data.mo_num);
                for (int32_t mo = 0; mo < data.mo_num; ++mo) {
                    data.orbitals[(size_t)mo].coefficients.assign(
                        coefficients.begin() + (ptrdiff_t)((size_t)mo * (size_t)data.ao_num),
                        coefficients.begin() + (ptrdiff_t)((size_t)(mo + 1) * (size_t)data.ao_num));
                }

                std::vector<double> occupation((size_t)data.mo_num, 0.0);
                std::vector<double> energy((size_t)data.mo_num, 0.0);
                std::vector<int32_t> spin((size_t)data.mo_num, 0);
                if (trexio_has_mo_occupation(file) == TREXIO_SUCCESS && trexio_read_mo_occupation(file, occupation.data()) == TREXIO_SUCCESS) {
                    for (int32_t i = 0; i < data.mo_num; ++i) data.orbitals[(size_t)i].occupation = occupation[(size_t)i];
                }
                if (trexio_has_mo_energy(file) == TREXIO_SUCCESS && trexio_read_mo_energy(file, energy.data()) == TREXIO_SUCCESS) {
                    for (int32_t i = 0; i < data.mo_num; ++i) data.orbitals[(size_t)i].energy = energy[(size_t)i];
                }
                if (trexio_has_mo_spin(file) == TREXIO_SUCCESS && trexio_read_mo_spin(file, spin.data()) == TREXIO_SUCCESS) {
                    for (int32_t i = 0; i < data.mo_num; ++i) data.orbitals[(size_t)i].spin = spin[(size_t)i];
                }

                std::vector<std::string> sym = read_string_array(file, data.mo_num, trexio_has_mo_symmetry, trexio_read_mo_symmetry);
                for (int32_t i = 0; i < data.mo_num && i < (int32_t)sym.size(); ++i) {
                    data.orbitals[(size_t)i].symmetry = sym[(size_t)i];
                }
            }
        }
    }

    close_file();
    return data;
}

} // namespace trexio_data

#endif // VIAMD_TREXIO
