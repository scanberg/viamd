// viamd_nci_reference: the NCI category references of a VeloxChem calculation, without the GUI.
//
// The same pipeline the NCI panel runs before "Export NCI samples...", for the category reference rows only:
// the analytic critical point search, the chemical labels (with ring surfaces for tied ring CPs), the CPs grouped
// into the reference categories (hydrogen bonds, rings, pi interactions, other), and one scatter-only extraction over
// every NCI CP with the uniform reference settings (nci_surface::reference_params: s <= 1.0, rho <= 0.07, 2.7 Bohr,
// both signs of lambda2), each sample to its nearest CP. Writes, in the format of the GUI export,
// <prefix>-density.csv (per category) and <prefix>-anchors.csv (per CP: its label at every level of detail, atoms,
// elements, molecules, group charges, contact class, type, rho and Hessian eigenvalues) for nci_reference_regions.py;
// <prefix>-cp-density.csv, the density per CP, and <prefix>-samples.bin, up to 50,000 raw samples with their CP, so
// CPs may be regrouped and the map rebinned later without recomputing; and <prefix>.xyz, the geometry.
//
// The density is evaluated on the CPU (md_topo_gto_density_derivatives), only within reach of the anchors.

#include <md_gto.h>
#include <md_qm.h>
#include <md_system.h>
#include <md_topo.h>
#include <md_util.h>
#include <md_vlx.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <topology/critical_point_distribution_core.h>
#include <topology/nci_surface_core.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

constexpr double angstrom_to_bohr = 1.0 / 0.529177210903;
constexpr double bohr_to_angstrom = 0.529177210903;
// As the NCI panel: the grid spans the atoms padded by 2.5 Bohr, clipped to the anchors' reach
constexpr double atom_margin = 2.5;
constexpr int category_count = 4;
static_assert(category_count == nci_surface::reference_category_count);
const char* category_names[category_count] = {
    nci_surface::reference_category_name(0), nci_surface::reference_category_name(1),
    nci_surface::reference_category_name(2), nci_surface::reference_category_name(3)};

int category_of(nci_surface::Family family) {
    const size_t category = nci_surface::reference_category(family);
    return category < nci_surface::reference_category_count ? (int)category : -1;
}

// One set of extraction settings, the same for every CP; without --settings the reference's
struct SettingSet {
    std::string tag;   // empty: the reference settings, written to PREFIX-density.csv
    double s_cut = nci_surface::reference_scatter_s, rho_max = nci_surface::reference_rho_max;
    double radius = nci_surface::reference_radius;   // Bohr
};

struct Options {
    std::string input, prefix;
    std::vector<SettingSet> sets;
    double spacing = 0.2;     // Bohr, the panel's Medium
    size_t max_samples = 64u << 20;
    double rho_min = 1.0e-4;  // the CP search threshold, as the panel's default
};

void usage() {
    fprintf(stderr,
            "usage: viamd_nci_reference [--spacing BOHR] [--max-samples N] [--rho-min AU] [--settings TAG:S:RHO[:R] ...]\n"
            "                           INPUT.h5 OUTPUT_PREFIX\n"
            "Writes OUTPUT_PREFIX-density.csv and OUTPUT_PREFIX-anchors.csv: the NCI category references of the\n"
            "calculation, as VIAMD's NCI export writes them; OUTPUT_PREFIX-cp-density.csv, the same per CP;\n"
            "OUTPUT_PREFIX-samples.bin, raw samples with their CP; and OUTPUT_PREFIX.xyz. --spacing defaults to 0.2 Bohr (Medium).\n"
            "--settings (repeatable) uses s <= S, rho <= RHO and radius R (Bohr, default %.1f) instead of the\n"
            "reference's (s <= %.1f, rho <= %.2f). Each set writes OUTPUT_PREFIX-TAG-density.csv and\n"
            "OUTPUT_PREFIX-TAG-cp-density.csv and -samples.bin; the critical points are searched once for all sets.\n",
            nci_surface::reference_radius, nci_surface::reference_scatter_s, nci_surface::reference_rho_max);
}

bool parse(Options* o, int argc, char** argv) {
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const auto value = [&](double* out) {
            if (i + 1 >= argc) return false;
            char* end = nullptr;
            *out = std::strtod(argv[++i], &end);
            return end && *end == 0 && std::isfinite(*out) && *out > 0.0;
        };
        double v = 0.0;
        if (a == "--spacing") { if (!value(&o->spacing)) return false; }
        else if (a == "--rho-min") { if (!value(&o->rho_min)) return false; }
        else if (a == "--max-samples") { if (!value(&v)) return false; o->max_samples = (size_t)v; }
        else if (a == "--settings") {
            if (i + 1 >= argc) return false;
            const std::string spec = argv[++i];
            const size_t c1 = spec.find(':'), c2 = c1 == std::string::npos ? c1 : spec.find(':', c1 + 1);
            if (c1 == 0 || c2 == std::string::npos) return false;
            SettingSet set;
            set.tag = spec.substr(0, c1);
            char* end = nullptr;
            set.s_cut = std::strtod(spec.c_str() + c1 + 1, &end);
            if (end != spec.c_str() + c2) return false;
            set.rho_max = std::strtod(spec.c_str() + c2 + 1, &end);
            if (*end == ':') set.radius = std::strtod(end + 1, &end);
            if (*end || !(set.s_cut > 0.0) || !(set.rho_max > 0.0) || !(set.radius > 0.0)) return false;
            o->sets.push_back(set);
        }
        else if (a == "-h" || a == "--help") return false;
        else if (!a.empty() && a[0] == '-') { fprintf(stderr, "unknown option %s\n", a.c_str()); return false; }
        else positional.push_back(a);
    }
    if (positional.size() != 2) return false;
    o->input = positional[0];
    o->prefix = positional[1];
    if (o->sets.empty()) o->sets.push_back({});
    return true;
}

double seconds_since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

std::string number(double value) {
    char text[48];
    snprintf(text, sizeof(text), "%.9g", value);
    return text;
}

void csv_row(std::string& out, const std::vector<std::string>& row) {
    for (size_t i = 0; i < row.size(); ++i) {
        if (i) out += ',';
        const std::string& field = row[i];
        if (field.find_first_of(",\"\n\r") == std::string::npos) {
            out += field;
            continue;
        }
        out += '"';
        for (char c : field) {
            if (c == '"') out += '"';
            out += c;
        }
        out += '"';
    }
    out += '\n';
}

bool write_file(const std::string& path, const std::string& text) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return false;
    const bool ok = fwrite(text.data(), 1, text.size(), f) == text.size();
    return fclose(f) == 0 && ok;
}

// The basis atoms' positions in Bohr: through the QM atom map when the file has one
bool basis_positions(std::vector<float>* xyz, const md_system_t& sys, const md_system_state_t& state, size_t num_basis_atoms) {
    const md_attribute_t* map_attr = md_attributes_find(&sys.attributes, STR_LIT("qm/atom/system_index"));
    const uint32_t* map = map_attr ? (const uint32_t*)md_attribute_view(map_attr, MD_ATTRIBUTE_TYPE_U32, 1, 1) : nullptr;
    if (map && map_attr->format.shape[0] < num_basis_atoms) return false;
    xyz->resize(3 * num_basis_atoms);
    for (size_t i = 0; i < num_basis_atoms; ++i) {
        const size_t idx = map ? (size_t)map[i] : i;
        if (idx >= state.num_atoms) return false;
        (*xyz)[3 * i + 0] = (float)(state.xyz[idx].x * angstrom_to_bohr);
        (*xyz)[3 * i + 1] = (float)(state.xyz[idx].y * angstrom_to_bohr);
        (*xyz)[3 * i + 2] = (float)(state.xyz[idx].z * angstrom_to_bohr);
    }
    return true;
}

// rho at the grid samples within reach of the anchors (0 elsewhere, out of every NCI region)
bool sample_density(std::vector<float>* rho, const nci_surface::Grid& grid, const std::vector<double>& anchors,
                    double reach, const md_topo_gto_desc_t& desc) {
    const int dx = grid.dim[0], dy = grid.dim[1], dz = grid.dim[2];
    const size_t n = (size_t)dx * dy * dz;
    std::vector<uint8_t> wanted(n, 0);
    for (size_t p = 0; p < anchors.size(); p += 3) {
        int lo[3], hi[3];
        for (int a = 0; a < 3; ++a) {
            lo[a] = std::max(0, (int)std::floor((anchors[p + a] - reach - grid.origin[a]) / grid.spacing[a]));
            hi[a] = std::min(grid.dim[a] - 1, (int)std::ceil((anchors[p + a] + reach - grid.origin[a]) / grid.spacing[a]));
        }
        for (int k = lo[2]; k <= hi[2]; ++k) for (int j = lo[1]; j <= hi[1]; ++j) for (int i = lo[0]; i <= hi[0]; ++i) {
            const int ijk[3] = {i, j, k};
            double d2 = 0.0;
            for (int a = 0; a < 3; ++a) {
                const double d = grid.origin[a] + ijk[a] * grid.spacing[a] - anchors[p + a];
                d2 += d * d;
            }
            if (d2 <= reach * reach) wanted[i + (size_t)dx * (j + (size_t)dy * k)] = 1;
        }
    }
    rho->assign(n, 0.0f);
    constexpr size_t chunk = 1u << 20;
    std::vector<double> xyz, value;
    std::vector<size_t> index;
    xyz.reserve(3 * chunk);
    index.reserve(chunk);
    const auto flush = [&]() {
        if (index.empty()) return true;
        value.assign(index.size(), NAN);
        if (!md_topo_gto_density_derivatives(value.data(), nullptr, nullptr, xyz.data(), index.size(), &desc)) return false;
        for (size_t t = 0; t < index.size(); ++t) (*rho)[index[t]] = std::isfinite(value[t]) ? (float)value[t] : 0.0f;
        xyz.clear();
        index.clear();
        return true;
    };
    for (int k = 0; k < dz; ++k) for (int j = 0; j < dy; ++j) for (int i = 0; i < dx; ++i) {
        const size_t v = i + (size_t)dx * (j + (size_t)dy * k);
        if (!wanted[v]) continue;
        xyz.push_back(grid.origin[0] + i * grid.spacing[0]);
        xyz.push_back(grid.origin[1] + j * grid.spacing[1]);
        xyz.push_back(grid.origin[2] + k * grid.spacing[2]);
        index.push_back(v);
        if (index.size() == chunk && !flush()) return false;
    }
    return flush();
}

}  // namespace

int main(int argc, char** argv) {
    Options opt;
    if (!parse(&opt, argc, argv)) {
        usage();
        return 2;
    }
    const auto t_start = std::chrono::steady_clock::now();

    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), 64u << 20);
    md_system_t sys = {};
    md_system_state_t state = {};
    sys.alloc = alloc;
    state.alloc = alloc;
    if (!md_vlx_system_init_from_file(&sys, &state, str_from_cstr(opt.input.c_str()))) {
        fprintf(stderr, "%s: not a VeloxChem file this build can read\n", opt.input.c_str());
        return 1;
    }
    if (!md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL))
        fprintf(stderr, "%s: nothing could be inferred; labels use the bare atoms\n", opt.input.c_str());

    md_gto_basis_t basis = {};
    if (!md_gto_basis_extract_attributes(&basis, &sys.attributes, alloc)) {
        fprintf(stderr, "%s: no GTO basis\n", opt.input.c_str());
        return 1;
    }
    const size_t num_ao = md_gto_basis_num_ao(&basis);
    std::vector<float> atom_xyz;
    if (!basis_positions(&atom_xyz, sys, state, md_gto_basis_num_atoms(&basis))) {
        fprintf(stderr, "%s: no positions for the basis atoms\n", opt.input.c_str());
        return 1;
    }
    // QTAIM is defined on the total density; a file with one spin channel publishes the total as alpha
    const md_attribute_t* density_attr = md_attributes_find(&sys.attributes, STR_LIT("orbital/total/density"));
    if (!density_attr) density_attr = md_attributes_find(&sys.attributes, STR_LIT("orbital/alpha/density"));
    const size_t dim = density_attr ? md_qm_extract_symmetric_f64(nullptr, 0, density_attr, md_attribute_slice_all()) : 0;
    if (dim == 0 || dim != num_ao) {
        fprintf(stderr, "%s: no density matrix matching the basis (%zu AOs)\n", opt.input.c_str(), num_ao);
        return 1;
    }
    std::vector<double> density_matrix(dim * dim);
    if (md_qm_extract_symmetric_f64(density_matrix.data(), density_matrix.size(), density_attr, md_attribute_slice_all()) != dim) {
        fprintf(stderr, "%s: the density matrix could not be read\n", opt.input.c_str());
        return 1;
    }

    md_topo_gto_desc_t desc = {};
    desc.basis = &basis;
    desc.atom_xyz = atom_xyz.data();
    desc.atom_xyz_stride = sizeof(float) * 3;
    desc.density_matrix = density_matrix.data();
    desc.rho_min = opt.rho_min;
    desc.trace_separatrices = true;

    // Critical points, labelled as the panel labels them
    auto t = std::chrono::steady_clock::now();
    md_topo_extremum_graph_t graph = {};
    graph.alloc = md_get_heap_allocator();
    md_topo_gto_info_t info = {};
    if (!md_topo_compute_extremum_graph_gto(&graph, &info, &desc) || !graph.vertices_f64) {
        fprintf(stderr, "%s: the critical point search failed\n", opt.input.c_str());
        return 1;
    }
    const double cp_seconds = seconds_since(t);
    std::vector<critical_point_distribution::Position> atoms(state.num_atoms);
    for (size_t i = 0; i < atoms.size(); ++i) atoms[i] = {state.xyz[i].x, state.xyz[i].y, state.xyz[i].z};
    std::vector<critical_point_distribution::ChemicalLabel> labels = critical_point_distribution::chemical_labels(graph, sys, atoms);
    std::vector<uint32_t> rcps;
    for (uint32_t i = 0; i < (uint32_t)graph.num_vertices; ++i)
        if (graph.types[i] == MD_TOPO_JOIN_SADDLE && critical_point_distribution::ring_surface_may_decide(labels[i])) rcps.push_back(i);
    if (!rcps.empty()) {
        constexpr uint32_t capacity = 32;
        std::vector<uint32_t> maxima(rcps.size() * capacity), counts(rcps.size());
        if (md_topo_gto_ring_surface_maxima(maxima.data(), counts.data(), capacity, &desc, &graph, rcps.data(), rcps.size())) {
            std::vector<std::vector<uint32_t>> surfaces(graph.num_vertices);
            for (size_t k = 0; k < rcps.size(); ++k)
                surfaces[rcps[k]].assign(maxima.begin() + k * capacity, maxima.begin() + k * capacity + std::min(counts[k], capacity));
            labels = critical_point_distribution::chemical_labels(graph, sys, atoms, surfaces);
        }
    }

    // Every NCI CP is an anchor; the CP's type decides its category, as in the export (-1: no reference category)
    std::vector<double> anchors;
    std::vector<int> anchor_category;
    std::vector<uint32_t> anchor_cp;
    std::array<size_t, category_count> category_cps = {};
    const auto element = [&](int atom) -> std::string {
        if (atom < 0 || (size_t)atom >= state.num_atoms) return "";
        const str_t symbol = md_atomic_number_symbol(md_atom_atomic_number(&sys.atom, atom));
        return std::string(symbol.ptr, symbol.len);
    };
    const auto molecule = [&](int atom) -> std::string {
        if (atom < 0 || (size_t)atom >= state.num_atoms) return "";
        return std::to_string((int)md_system_molecule_find_by_atom_idx(&sys, atom));
    };
    const auto charge = [&](int atom) -> std::string {
        if (atom < 0 || (size_t)atom >= state.num_atoms) return "";
        return std::to_string(critical_point_distribution::group_charge(sys, atom));
    };
    const auto cp_type_name = [](md_topo_critical_point_type_t type) -> std::string {
        switch (type) {
            case MD_TOPO_SPLIT_SADDLE: return "BCP";
            case MD_TOPO_JOIN_SADDLE: return "RCP";
            case MD_TOPO_MINIMUM: return "CCP";
            case MD_TOPO_MAXIMUM: return "NCP";
            default: return "";
        }
    };
    std::string anchors_csv;
    csv_row(anchors_csv, {"layer", "layer_id", "surface_id", "parameter_group", "reference_category", "cp_index", "cp_category",
                          "chemical_environment", "cp_x (angstrom)", "cp_y (angstrom)", "cp_z (angstrom)", "kind",
                          "interaction", "group", "detail", "atom_a", "atom_b", "element_a", "element_b",
                          "molecule_a", "molecule_b", "charge_a", "charge_b", "cp_type", "rho (au)", "lambda1 (au)", "lambda2 (au)", "lambda3 (au)",
                          "laplacian (au)", "finer_1", "finer_2", "finer_3"});
    for (uint32_t i = 0; i < (uint32_t)graph.num_vertices; ++i) {
        const nci_surface::Family family = nci_surface::family_of(labels[i]);
        const nci_surface::Family parameter_family = nci_surface::parameter_family_of(labels[i]);
        if (!nci_surface::has_nci(parameter_family)) continue;
        const int category = category_of(family);
        const auto& v = graph.vertices_f64[i];
        anchors.insert(anchors.end(), {v.x, v.y, v.z});
        anchor_category.push_back(category);
        anchor_cp.push_back(i);
        if (category >= 0) ++category_cps[category];
        std::string environment;
        std::vector<std::string> path = critical_point_distribution::environment_path(labels[i], sys);
        path.push_back(critical_point_distribution::display_name(labels[i]));
        for (const std::string& part : path) environment += (environment.empty() ? "" : " > ") + part;
        const auto& label = labels[i];
        const auto finer = [&](size_t k) { return k < label.finer.size() ? label.finer[k] : std::string(); };
        csv_row(anchors_csv, {"Full system", "-1", "0", nci_surface::family_name(parameter_family),
                              category >= 0 ? category_names[category] : "Unassigned",
                              std::to_string(i), nci_surface::family_name(family), environment,
                              number(v.x * bohr_to_angstrom), number(v.y * bohr_to_angstrom), number(v.z * bohr_to_angstrom),
                              critical_point_distribution::kind_name(label.kind), label.interaction, label.group, label.detail,
                              std::to_string(label.atom[0]), std::to_string(label.atom[1]),
                              element(label.atom[0]), element(label.atom[1]), molecule(label.atom[0]), molecule(label.atom[1]),
                              charge(label.atom[0]), charge(label.atom[1]),
                              cp_type_name(graph.types[i]), number(v.value), number(v.hessian_eigenvalues[0]),
                              number(v.hessian_eigenvalues[1]), number(v.hessian_eigenvalues[2]),
                              number(v.hessian_eigenvalues[0] + v.hessian_eigenvalues[1] + v.hessian_eigenvalues[2]),
                              finer(0), finer(1), finer(2)});
    }

    // The geometry, so the CPs can be reclassified later
    std::string xyz = std::to_string(state.num_atoms) + "\n" + opt.input + "\n";
    for (size_t i = 0; i < state.num_atoms; ++i)
        xyz += element((int)i) + " " + number(state.xyz[i].x) + " " + number(state.xyz[i].y) + " " + number(state.xyz[i].z) + "\n";

    double atoms_lo[3] = {INFINITY, INFINITY, INFINITY}, atoms_hi[3] = {-INFINITY, -INFINITY, -INFINITY};
    for (size_t i = 0; i < atom_xyz.size(); i += 3)
        for (int a = 0; a < 3; ++a) {
            atoms_lo[a] = std::min(atoms_lo[a], (double)atom_xyz[i + a] - atom_margin);
            atoms_hi[a] = std::max(atoms_hi[a], (double)atom_xyz[i + a] + atom_margin);
        }

    std::string header;
    csv_row(header, {"layer", "layer_id", "surface_id", "sample_scope", "parameter_group", "category",
                          "extraction_settings", "bin_x", "bin_y", "sign(lambda2)rho_min (au)", "sign(lambda2)rho_max (au)",
                          "s_min (dimensionless)", "s_max (dimensionless)", "volume (bohr^3)", "total_volume (bohr^3)",
                          "rho_cutoff (au)", "s_isovalue (dimensionless)", "grid_spacing_x (bohr)", "grid_spacing_y (bohr)",
                          "grid_spacing_z (bohr)", "resolution", "connected", "radius_cap", "near_radius (bohr)"});
    std::string cp_header;
    csv_row(cp_header, {"cp_index", "reference_category", "bin_x", "bin_y", "volume (bohr^3)"});
    std::vector<std::string> density_csv(opt.sets.size(), header), cp_density_csv(opt.sets.size(), cp_header);
    // Raw samples: "NCISAMP1", uint64 count, uint64 eligible count (the samples are a uniform random subset of it),
    // float64 voxel volume (bohr^3), then per sample float32 sign(lambda2)rho, float32 s, uint32 cp_index
    std::vector<std::string> samples_bin(opt.sets.size());
    printf("%s: %zu atoms, %zu AOs, %zu CPs (%zu NCI) in %.1f s (%s)\n", opt.input.c_str(), (size_t)state.num_atoms, num_ao,
           (size_t)graph.num_vertices, anchors.size() / 3, cp_seconds, info.complete ? "complete" : "incomplete");
    for (int c = 0; c < category_count; ++c) printf("  %-15s %5zu CPs\n", category_names[c], category_cps[c]);

    if (!anchors.empty()) {
        t = std::chrono::steady_clock::now();
        // One grid and one density for all sets: sized for the largest radius
        double reach = 0.0;
        for (const SettingSet& set : opt.sets) reach = std::max(reach, set.radius);
        double lo[3], hi[3];
        for (int a = 0; a < 3; ++a) {
            double alo = INFINITY, ahi = -INFINITY;
            for (size_t p = a; p < anchors.size(); p += 3) {
                alo = std::min(alo, anchors[p]);
                ahi = std::max(ahi, anchors[p]);
            }
            lo[a] = std::max(atoms_lo[a], alo - reach - 1.0);
            hi[a] = std::min(atoms_hi[a], ahi + reach + 1.0);
        }
        if (!(hi[0] > lo[0] && hi[1] > lo[1] && hi[2] > lo[2])) {
            fprintf(stderr, "%s: the CPs are outside the atoms' box\n", opt.input.c_str());
            return 1;
        }
        double h = opt.spacing;
        int gdim[3];
        for (int pass = 0; pass < 2; ++pass) {
            size_t total = 1;
            for (int a = 0; a < 3; ++a) {
                gdim[a] = std::max(8, ((int)std::ceil((hi[a] - lo[a]) / h) + 7) / 8 * 8);
                total *= (size_t)gdim[a];
            }
            if (total <= opt.max_samples) break;
            h *= std::cbrt((double)total / (double)opt.max_samples) * 1.01;
        }
        nci_surface::Grid grid;
        for (int a = 0; a < 3; ++a) {
            grid.dim[a] = gdim[a];
            grid.spacing[a] = (hi[a] - lo[a]) / gdim[a];
            grid.origin[a] = lo[a] + 0.5 * grid.spacing[a];
        }
        const double h_max = std::max({grid.spacing[0], grid.spacing[1], grid.spacing[2]});
        std::vector<float> rho;
        // The flood stays within the radius; the derivatives there reach two samples further
        if (!sample_density(&rho, grid, anchors, reach + 3.0 * h_max, desc)) {
            fprintf(stderr, "%s: the density could not be evaluated\n", opt.input.c_str());
            return 1;
        }
        const double rho_seconds = seconds_since(t);

        constexpr int bins = nci_surface::scatter_density_bins;
        for (size_t k = 0; k < opt.sets.size(); ++k) {
            const SettingSet& set = opt.sets[k];
            const auto t_set = std::chrono::steady_clock::now();
            nci_surface::Params params = nci_surface::reference_params(anchors);
            params.s_iso = params.scatter_s = set.s_cut;
            params.rho_max = set.rho_max;
            params.near_radius = set.radius;
            nci_surface::Mesh mesh;
            nci_surface::Scatter scatter;
            if (!nci_surface::extract(&mesh, rho.data(), grid, params, nullptr, nullptr, &scatter)) {
                fprintf(stderr, "%s: the extraction failed\n", opt.input.c_str());
                return 1;
            }
            {
                std::string& out = samples_bin[k];
                const auto put = [&](const void* p, size_t n) { out.append((const char*)p, n); };
                const uint64_t count = scatter.s.size(), eligible = scatter.eligible_count;
                const double voxel = grid.spacing[0] * grid.spacing[1] * grid.spacing[2];
                out = "NCISAMP1";
                put(&count, 8);
                put(&eligible, 8);
                put(&voxel, 8);
                for (size_t p = 0; p < count; ++p) {
                    const uint32_t cp = p < scatter.anchor.size() && scatter.anchor[p] < anchor_cp.size() ? anchor_cp[scatter.anchor[p]] : UINT32_MAX;
                    put(&scatter.signed_rho[p], 4);
                    put(&scatter.s[p], 4);
                    put(&cp, 4);
                }
            }
            const nci_surface::ScatterDensity& density = scatter.density;
            if (!(density.rho_extent > 0.0f) || density.volume.size() != (size_t)bins * bins) continue;
            std::array<std::vector<double>, category_count> volume;
            std::array<double, category_count> total = {}, positive = {};
            for (auto& v : volume) v.assign((size_t)bins * bins, 0.0);
            for (const nci_surface::AnchorBin& b : density.by_anchor) {
                if (b.anchor >= anchor_category.size() || b.bin >= (size_t)bins * bins) continue;
                const int c = anchor_category[b.anchor];
                csv_row(cp_density_csv[k], {std::to_string(anchor_cp[b.anchor]), c >= 0 ? category_names[c] : "Unassigned",
                                            std::to_string(b.bin % bins), std::to_string(b.bin / bins), number(b.volume)});
                if (c < 0) continue;
                volume[c][b.bin] += b.volume;
                total[c] += b.volume;
                if ((int)(b.bin % bins) >= bins / 2) positive[c] += b.volume;
            }
            const double w = 2.0 * density.rho_extent / bins, bh = nci_surface::scatter_s_max / (double)bins;
            for (int c = 0; c < category_count; ++c) {
                const nci_surface::Family settings = nci_surface::reference_settings_family(c);
                for (int by = 0; by < bins; ++by) for (int bx = 0; bx < bins; ++bx) {
                    const double v = volume[c][bx + (size_t)bins * by];
                    if (!(v > 0.0)) continue;
                    csv_row(density_csv[k], {"Full system", "-1", "0", "category", nci_surface::family_name(settings), category_names[c],
                                             "category_reference_uniform", std::to_string(bx), std::to_string(by),
                                             number(-density.rho_extent + bx * w), number(-density.rho_extent + (bx + 1) * w),
                                             number(by * bh), number((by + 1) * bh), number(v), number(total[c]),
                                             number(params.rho_max), number(params.scatter_s), number(grid.spacing[0]),
                                             number(grid.spacing[1]), number(grid.spacing[2]), "-1", "true", "true",
                                             number(params.near_radius)});
                }
                if (category_cps[c])
                    printf("  %-15s %-8s %.1f Bohr^3 (%.0f%% positive)\n", category_names[c], set.tag.empty() ? "default" : set.tag.c_str(),
                           total[c], total[c] > 0 ? 100.0 * positive[c] / total[c] : 0.0);
            }
            printf("  grid %dx%dx%d at %.3f Bohr, density %.1f s, extraction %.1f s, %.1f Bohr^3 in all\n", gdim[0], gdim[1], gdim[2],
                   h_max, rho_seconds, seconds_since(t_set), density.total_volume);
        }
    }

    const auto write = [&](const std::string& path, const std::string& text) {
        if (write_file(path, text)) return true;
        fprintf(stderr, "could not write %s\n", path.c_str());
        return false;
    };
    if (!write(opt.prefix + "-anchors.csv", anchors_csv) || !write(opt.prefix + ".xyz", xyz)) return 1;
    for (size_t k = 0; k < opt.sets.size(); ++k) {
        const std::string stem = opt.prefix + (opt.sets[k].tag.empty() ? "" : "-" + opt.sets[k].tag);
        if (!write(stem + "-density.csv", density_csv[k]) || !write(stem + "-cp-density.csv", cp_density_csv[k]) ||
            (!samples_bin[k].empty() && !write(stem + "-samples.bin", samples_bin[k])))
            return 1;
        printf("  wrote %s-density.csv\n", stem.c_str());
    }
    printf("  done in %.1f s\n", seconds_since(t_start));
    md_topo_extremum_graph_free(&graph);
    return 0;
}
