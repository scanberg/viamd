#pragma once

#include <stddef.h>
#include <stdint.h>

#include <core/md_vec_math.h>

// The frame of the QM atoms as a file had them (qm/atom/coordinate, Angstrom), and the rigid motion taking it onto
// where those atoms are: what a file published in that frame (a cavity's tesserae, a grid) follows the atoms with.
//
// x = R x_file + t (column major R). xyz holds the positions now, `scale` of them per Angstrom (1: Angstrom,
// 1 / 0.529...: bohr), R and t act on the file's positions so scaled; idx the entry of each QM atom in xyz (NULL: the
// QM atom's own index). The identity when the atoms have not moved; with fewer than three atoms only a translation.
struct QmPlacement {
    double R[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    double t[3] = {0, 0, 0};
    double scale = 1.0;                 // the units of x per Angstrom
};

bool qm_place(QmPlacement* out, const double* qm_xyz, uint32_t num_qm_atoms, const md_float3_t* xyz, size_t count, const uint32_t* idx, double scale);

// A point (Angstrom, the file's frame) where it is now, in the units of the placement
static inline void qm_place_point(const QmPlacement& p, const double in[3], double out[3]) {
    const double x[3] = {in[0] * p.scale, in[1] * p.scale, in[2] * p.scale};
    for (int k = 0; k < 3; ++k) out[k] = p.R[k] * x[0] + p.R[3 + k] * x[1] + p.R[6 + k] * x[2] + p.t[k];
}

// A direction, rotated
static inline void qm_place_dir(const QmPlacement& p, const double in[3], double out[3]) {
    for (int k = 0; k < 3; ++k) out[k] = p.R[k] * in[0] + p.R[3 + k] * in[1] + p.R[6 + k] * in[2];
}
