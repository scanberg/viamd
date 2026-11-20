#include <qm_placement.h>

#include <core/md_geom.h>

#include <math.h>
#include <stdlib.h>

bool qm_place(QmPlacement* out, const double* qm_xyz, uint32_t m, const md_float3_t* xyz, size_t count, const uint32_t* idx, double scale) {
    *out = QmPlacement{};
    out->scale = scale;
    if (!qm_xyz || !xyz || m == 0) return false;
    for (uint32_t i = 0; i < m; ++i) {
        if ((idx ? idx[i] : i) >= count) return false;
    }
    // The rigid motion taking the QM atoms of the file (scaled) onto where they are (Kabsch)
    double cf[3] = {0, 0, 0}, cw[3] = {0, 0, 0};
    for (uint32_t i = 0; i < m; ++i) {
        const md_float3_t w = xyz[idx ? idx[i] : i];
        for (int k = 0; k < 3; ++k) cf[k] += qm_xyz[3 * i + k] * scale / m;
        cw[0] += w.x / m; cw[1] += w.y / m; cw[2] += w.z / m;
    }
    md_float4_t* ref = (md_float4_t*)malloc(sizeof(md_float4_t) * m);
    md_float4_t* trg = (md_float4_t*)malloc(sizeof(md_float4_t) * m);
    if (!ref || !trg) { free(ref); free(trg); return false; }
    double dev = 0;
    for (uint32_t i = 0; i < m; ++i) {
        const md_float3_t w = xyz[idx ? idx[i] : i];
        ref[i] = {(float)(qm_xyz[3 * i + 0] * scale - cf[0]), (float)(qm_xyz[3 * i + 1] * scale - cf[1]), (float)(qm_xyz[3 * i + 2] * scale - cf[2]), 1.0f};
        trg[i] = {(float)(w.x - cw[0]), (float)(w.y - cw[1]), (float)(w.z - cw[2]), 1.0f};
        dev = fmax(dev, fabs(ref[i].x - trg[i].x) + fabs(ref[i].y - trg[i].y) + fabs(ref[i].z - trg[i].z));
    }
    // R maps the target onto the reference, R (q - cw) = p - cf; the file's frame goes the other way: q = R^T (p - cf) + cw
    if (m >= 3 && dev > 1e-4 * scale) {
        const md_float3x3_t M = md_geom_optimal_rotation_rel(ref, trg, m);
        for (int c = 0; c < 3; ++c)
            for (int r = 0; r < 3; ++r) out->R[3 * c + r] = M.col[r].elem[c];   // the transpose, column major
    }
    free(ref);
    free(trg);
    const double* R = out->R;
    for (int k = 0; k < 3; ++k) out->t[k] = cw[k] - (R[k] * cf[0] + R[3 + k] * cf[1] + R[6 + k] * cf[2]);
    return true;
}
