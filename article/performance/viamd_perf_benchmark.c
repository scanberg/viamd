#include <md_vlx.h>
#include <md_gto.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_grid.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define DEFAULT_REPEATS 9
#define DEFAULT_GRID_DIM 64
#define GTO_CUTOFF 1.0e-6

#ifndef ARTICLE_PERF_DATA_DIR
#define ARTICLE_PERF_DATA_DIR "ext/mdlib/test_data"
#endif

#define STR2(x) #x
#define STR(x) STR2(x)

typedef struct Dataset {
    const char* label;
    const char* rel_path;
} Dataset;

typedef struct Metrics {
    size_t atoms;
    size_t aos;
    size_t states;
    double file_kb;
} Metrics;

typedef bool (*bench_fn_t)(const Dataset* dataset, int grid_dim);

static int cmp_double(const void* a, const void* b) {
    const double da = *(const double*)a;
    const double db = *(const double*)b;
    return (da > db) - (da < db);
}

static const char* dataset_path(const Dataset* dataset) {
    static char path[4096];
    snprintf(path, sizeof(path), "%s/%s", STR(ARTICLE_PERF_DATA_DIR), dataset->rel_path);
    return path;
}

static double file_size_kb(const char* path) {
    struct stat st;
    if (stat(path, &st) != 0) return 0.0;
    return (double)st.st_size / 1024.0;
}

static Metrics read_metrics(const Dataset* dataset) {
    Metrics m = {0};
    const char* path = dataset_path(dataset);
    m.file_kb = file_size_kb(path);

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
    md_vlx_t* vlx = md_vlx_create(arena);
    if (md_vlx_parse_file(vlx, str_from_cstr(path))) {
        m.atoms = md_vlx_number_of_atoms(vlx);
        m.aos = md_vlx_scf_number_of_atomic_orbitals(vlx);
        m.states = md_vlx_rsp_number_of_excited_states(vlx);
    }
    md_arena_allocator_destroy(arena);
    return m;
}

static bool bench_parse(const Dataset* dataset, int grid_dim) {
    (void)grid_dim;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
    md_vlx_t* vlx = md_vlx_create(arena);
    bool ok = md_vlx_parse_file(vlx, str_from_cstr(dataset_path(dataset)));
    md_arena_allocator_destroy(arena);
    return ok;
}

static bool bench_basis_extract(const Dataset* dataset, int grid_dim) {
    (void)grid_dim;
    bool ok = false;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
    md_vlx_t* vlx = md_vlx_create(arena);
    if (md_vlx_parse_file(vlx, str_from_cstr(dataset_path(dataset)))) {
        md_gto_basis_t basis = {0};
        ok = md_vlx_gto_basis_extract(&basis, vlx, arena);
    }
    md_arena_allocator_destroy(arena);
    return ok;
}

static bool init_mo_grid(md_grid_t* grid, vec3_t** atom_xyz, md_gto_basis_t* basis, const double** coeffs, md_vlx_t* vlx, md_allocator_i* arena, int grid_dim) {
    if (!grid || !atom_xyz || !basis || !coeffs || !vlx) return false;

    const size_t num_atoms = md_vlx_number_of_atoms(vlx);
    const size_t num_aos = md_vlx_scf_number_of_atomic_orbitals(vlx);
    if (num_atoms == 0 || num_aos == 0) return false;

    const dvec3_t* coords = md_vlx_atom_coordinates(vlx);
    if (!coords) return false;

    *atom_xyz = md_arena_allocator_push(arena, sizeof(vec3_t) * num_atoms);
    vec3_t min_box = vec3_set1(FLT_MAX);
    vec3_t max_box = vec3_set1(-FLT_MAX);
    for (size_t i = 0; i < num_atoms; ++i) {
        (*atom_xyz)[i] = vec3_set((float)(coords[i].x * ANGSTROM_TO_BOHR),
                                  (float)(coords[i].y * ANGSTROM_TO_BOHR),
                                  (float)(coords[i].z * ANGSTROM_TO_BOHR));
        min_box = vec3_min(min_box, (*atom_xyz)[i]);
        max_box = vec3_max(max_box, (*atom_xyz)[i]);
    }

    min_box = vec3_sub1(min_box, 6.0f);
    max_box = vec3_add1(max_box, 6.0f);
    vec3_t spacing = vec3_div1(vec3_sub(max_box, min_box), (float)grid_dim);
    vec3_t origin = vec3_add(min_box, vec3_mul1(spacing, 0.5f));

    *grid = (md_grid_t){
        .orientation = mat3_ident(),
        .origin = origin,
        .spacing = spacing,
        .dim = {grid_dim, grid_dim, grid_dim},
    };

    if (!md_vlx_gto_basis_extract(basis, vlx, arena)) return false;
    size_t mo_idx = md_vlx_scf_homo_idx(vlx, MD_VLX_SPIN_ALPHA);
    *coeffs = md_vlx_scf_mo_coefficients(vlx, mo_idx, MD_VLX_SPIN_ALPHA);
    return *coeffs != NULL;
}

static bool bench_mo_volume(const Dataset* dataset, int grid_dim) {
    bool ok = false;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(128));
    md_vlx_t* vlx = md_vlx_create(arena);
    if (!md_vlx_parse_file(vlx, str_from_cstr(dataset_path(dataset)))) {
        md_arena_allocator_destroy(arena);
        return false;
    }

    md_grid_t grid = {0};
    vec3_t* atom_xyz = NULL;
    md_gto_basis_t basis = {0};
    const double* coeffs = NULL;
    if (init_mo_grid(&grid, &atom_xyz, &basis, &coeffs, vlx, arena, grid_dim)) {
        size_t max_gtos = md_gto_pgto_count(&basis);
        md_gto_t* gtos = md_arena_allocator_push(arena, sizeof(md_gto_t) * max_gtos);
        size_t num_gtos = md_gto_expand_with_ao_coeffs(gtos, &basis, (const float*)atom_xyz, sizeof(vec3_t), coeffs, GTO_CUTOFF);
        size_t num_points = md_grid_num_points(&grid);
        float* volume = md_arena_allocator_push_zero(arena, sizeof(float) * num_points);
        md_gto_grid_evaluate(volume, &grid, gtos, num_gtos, MD_GTO_EVAL_MODE_PSI);
        volatile float sink = volume[num_points / 2];
        (void)sink;
        ok = true;
    }

    md_arena_allocator_destroy(arena);
    return ok;
}

static bool bench_transition_density_extract(const Dataset* dataset, int grid_dim) {
    (void)grid_dim;
    bool ok = false;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(128));
    md_vlx_t* vlx = md_vlx_create(arena);
    if (md_vlx_parse_file(vlx, str_from_cstr(dataset_path(dataset)))) {
        if (md_vlx_rsp_number_of_excited_states(vlx) > 0) {
            size_t dim = md_vlx_rsp_transition_density_matrix_size(vlx, 0);
            if (dim > 0) {
                double* matrix = md_arena_allocator_push(arena, sizeof(double) * dim * dim);
                ok = md_vlx_rsp_transition_density_matrix_extract(matrix, vlx, 0, MD_VLX_TRANSITION_DENSITY_DIFFERENCE) == dim;
                volatile double sink = matrix[(dim * dim) / 2];
                (void)sink;
            }
        }
    }
    md_arena_allocator_destroy(arena);
    return ok;
}

static void run_benchmark(const Dataset* dataset, const char* operation, bench_fn_t fn, int grid_dim, int repeats) {
    Metrics m = read_metrics(dataset);
    double* samples = calloc((size_t)repeats, sizeof(double));
    int n = 0;

    // Warm-up.
    fn(dataset, grid_dim);

    for (int i = 0; i < repeats; ++i) {
        md_timestamp_t t0 = md_time_now();
        bool ok = fn(dataset, grid_dim);
        md_timestamp_t t1 = md_time_now();
        if (ok) {
            samples[n++] = md_time_as_milliseconds(t1 - t0);
        }
    }

    if (n == 0) {
        free(samples);
        return;
    }

    qsort(samples, (size_t)n, sizeof(double), cmp_double);
    double sum = 0.0;
    for (int i = 0; i < n; ++i) sum += samples[i];
    double median = samples[n / 2];
    double mean = sum / (double)n;
    double min_v = samples[0];
    double max_v = samples[n - 1];

    printf("%s,%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%zu,%zu,%zu,%.3f\n",
           dataset->label, operation, grid_dim, n, min_v, median, mean, max_v,
           m.atoms, m.aos, m.states, m.file_kb);
    free(samples);
}

int main(int argc, char** argv) {
    int grid_dim = DEFAULT_GRID_DIM;
    int repeats = DEFAULT_REPEATS;
    if (argc > 1) grid_dim = atoi(argv[1]);
    if (argc > 2) repeats = atoi(argv[2]);
    if (grid_dim <= 0) grid_dim = DEFAULT_GRID_DIM;
    if (repeats <= 0) repeats = DEFAULT_REPEATS;

    const Dataset datasets[] = {
        {"H2O", "vlx/h2o.h5"},
        {"Molecule", "vlx/mol.h5"},
        {"TQ SCF", "vlx/tq.scf.results.h5"},
        {"Acrolein RSP", "vlx/acro-rsp.h5"},
    };

    printf("dataset,operation,grid_dim,repeats,min_ms,median_ms,mean_ms,max_ms,atoms,aos,states,file_kb\n");
    for (size_t i = 0; i < sizeof(datasets) / sizeof(datasets[0]); ++i) {
        run_benchmark(&datasets[i], "HDF5 parse", bench_parse, grid_dim, repeats);
        run_benchmark(&datasets[i], "Basis extraction", bench_basis_extract, grid_dim, repeats);
        run_benchmark(&datasets[i], "MO volume", bench_mo_volume, grid_dim, repeats);
        run_benchmark(&datasets[i], "Transition density", bench_transition_density_extract, grid_dim, repeats);
    }
    return 0;
}
