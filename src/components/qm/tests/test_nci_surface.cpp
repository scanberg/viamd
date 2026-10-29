#include "utest.h"
#include <qm/nci_surface_core.h>
#include <qm/critical_point_distribution_core.h>

#include <array>
#include <atomic>
#include <cmath>
#include <vector>

namespace {
// Exponential atoms a exp(-b r), sampled on a grid
struct Atoms {
    std::vector<std::array<double, 3>> centre;
    double a = 0.5, b = 1.8;
    double rho(const double x[3]) const {
        double r = 0.0;
        for (const auto& c : centre) r += a * std::exp(-b * std::sqrt((x[0] - c[0]) * (x[0] - c[0]) + (x[1] - c[1]) * (x[1] - c[1]) + (x[2] - c[2]) * (x[2] - c[2])));
        return r;
    }
    std::vector<float> sample(const nci_surface::Grid& g) const {
        std::vector<float> v((size_t)g.dim[0] * g.dim[1] * g.dim[2]);
        for (int k = 0; k < g.dim[2]; ++k) for (int j = 0; j < g.dim[1]; ++j) for (int i = 0; i < g.dim[0]; ++i) {
            const double x[3] = {g.origin[0] + i * g.spacing[0], g.origin[1] + j * g.spacing[1], g.origin[2] + k * g.spacing[2]};
            v[i + (size_t)g.dim[0] * (j + (size_t)g.dim[1] * k)] = (float)rho(x);
        }
        return v;
    }
};

nci_surface::Grid box(double half_x, double half_yz, double h) {
    nci_surface::Grid g;
    const double half[3] = {half_x, half_yz, half_yz};
    for (int a = 0; a < 3; ++a) {
        g.dim[a] = (int)std::lround(2.0 * half[a] / h) + 1;
        g.origin[a] = -half[a];
        g.spacing[a] = h;
    }
    return g;
}
}  // namespace

UTEST(nci_surface, finite_differences) {
    Atoms pair;
    pair.centre = {{-3.0, 0.0, 0.0}, {3.0, 0.0, 0.0}};
    const auto g = box(6.0, 4.0, 0.1);
    const auto rho = pair.sample(g);
    // At (0, 1, 0.5): rho, gradient analytic
    const int i = 60, j = 50, k = 45;
    double s, signed_rho;
    nci_surface::sample(&s, &signed_rho, rho.data(), g, i, j, k);
    const double x[3] = {0.0, 1.0, 0.5};
    const double r = pair.rho(x);
    double grad[3] = {0, 0, 0};
    for (const auto& c : pair.centre) {
        const double d[3] = {x[0] - c[0], x[1] - c[1], x[2] - c[2]};
        const double dist = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        for (int a = 0; a < 3; ++a) grad[a] += -pair.b * pair.a * std::exp(-pair.b * dist) * d[a] / dist;
    }
    EXPECT_NEAR(s, critical_point_distribution::reduced_density_gradient(r, grad), 0.02 * s);
    EXPECT_NEAR(std::abs(signed_rho), r, 1e-6);
}

UTEST(nci_surface, pair_and_crowded) {
    Atoms pair;
    pair.centre = {{-3.0, 0.0, 0.0}, {3.0, 0.0, 0.0}};
    const auto g = box(6.0, 5.0, 0.2);
    const auto rho = pair.sample(g);
    nci_surface::Params params;
    nci_surface::Mesh mesh;
    std::atomic<float> progress = 0.0f;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, &progress));
    EXPECT_EQ(progress.load(), 1.0f);
    const size_t nv = mesh.value.size();
    ASSERT_GT(nv, (size_t)50);
    ASSERT_EQ(mesh.position.size(), 3 * nv);
    ASSERT_EQ(mesh.normal.size(), 3 * nv);
    ASSERT_GT(mesh.index.size(), (size_t)0);
    EXPECT_EQ(mesh.index.size() % 3, (size_t)0);
    // A disk across the bond path, between the atoms, attractive (lambda2 < 0), facing out of the disk
    size_t negative = 0, far = 0;
    float max_x = -1e9f;
    size_t at_max = 0;
    for (size_t v = 0; v < nv; ++v) {
        negative += mesh.value[v] < 0.0f;
        far += std::abs(mesh.position[3 * v]) > 2.0f;
        if (mesh.position[3 * v] > max_x) { max_x = mesh.position[3 * v]; at_max = v; }
    }

    EXPECT_EQ(far, (size_t)0);
    EXPECT_GT(negative, nv * 95 / 100);
    EXPECT_GT(mesh.normal[3 * at_max], 0.5f);
    for (uint32_t t : mesh.index) ASSERT_LT(t, (uint32_t)nv);

    // A third atom closing a crowded triangle: the ring's depletion, lambda2 > 0, shows
    Atoms triangle;
    triangle.centre = {{-3.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {0.0, 4.5, 0.0}};
    const auto g3 = box(6.0, 7.0, 0.2);
    const auto rho3 = triangle.sample(g3);
    nci_surface::Mesh crowded;
    nci_surface::Scatter scatter;
    ASSERT_TRUE(nci_surface::extract(&crowded, rho3.data(), g3, params, nullptr, nullptr, &scatter));
    size_t positive = 0;
    for (float v : crowded.value) positive += v > 0.0f;
    EXPECT_GT(positive, (size_t)10);
    size_t attractive_samples = 0, repulsive_samples = 0;
    for (float v : scatter.signed_rho) {
        attractive_samples += v < 0;
        repulsive_samples += v > 0;
    }
    EXPECT_GT(attractive_samples, (size_t)0);
    EXPECT_GT(repulsive_samples, (size_t)0);

    // Only near given points: nothing near a far corner; cancelled: nothing
    const double corner[3] = {-5.5, -4.5, -4.5};
    params.near_xyz = std::span<const double>(corner, 3);
    params.near_radius = 1.0;
    nci_surface::Mesh none;
    ASSERT_TRUE(nci_surface::extract(&none, rho.data(), g, params));
    EXPECT_EQ(none.value.size(), (size_t)0);
    const double centre[3] = {0.0, 0.0, 0.0};
    params.near_xyz = std::span<const double>(centre, 3);
    params.near_radius = 0.5;
    nci_surface::Mesh near;
    ASSERT_TRUE(nci_surface::extract(&near, rho.data(), g, params));
    EXPECT_GT(near.value.size(), (size_t)0);
    EXPECT_LT(near.value.size(), nv);
    std::atomic<int32_t> cancel = 1;
    EXPECT_FALSE(nci_surface::extract(&none, rho.data(), g, {}, &cancel));
    EXPECT_TRUE(none.position.empty());
    EXPECT_TRUE(none.index.empty());
    cancel.store(0);
    ASSERT_TRUE(nci_surface::extract(&none, rho.data(), g, params, &cancel, &progress));
    EXPECT_EQ(none.value.size(), near.value.size());
    EXPECT_EQ(progress.load(), 1.0f);
}

UTEST(nci_surface, scatter_samples_volume_and_respects_scope) {
    Atoms pair;
    pair.centre = {{-3.0, 0.0, 0.0}, {3.0, 0.0, 0.0}};
    const auto g = box(6.0, 5.0, 0.2);
    const auto rho = pair.sample(g);
    nci_surface::Params params;
    const double centre[] = {0, 0, 0};
    params.near_xyz = centre;
    params.near_radius = 1.0;
    nci_surface::Mesh mesh;
    nci_surface::Scatter scatter;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, nullptr, &scatter));
    ASSERT_GT(scatter.s.size(), (size_t)0);
    EXPECT_EQ(scatter.s.size(), scatter.signed_rho.size());
    EXPECT_EQ(scatter.s.size(), scatter.eligible_count);
    size_t expected = 0, away_from_isovalue = 0;
    for (int k = 1; k < g.dim[2] - 1; ++k) for (int j = 1; j < g.dim[1] - 1; ++j) for (int i = 1; i < g.dim[0] - 1; ++i) {
        const double x = g.origin[0] + i * g.spacing[0];
        const double y = g.origin[1] + j * g.spacing[1];
        const double z = g.origin[2] + k * g.spacing[2];
        if (x * x + y * y + z * z > params.near_radius * params.near_radius) continue;
        const float r = rho[i + (size_t)g.dim[0] * (j + (size_t)g.dim[1] * k)];
        if (!(r > 0) || r > params.rho_max) continue;
        double s, signed_rho;
        nci_surface::sample(&s, &signed_rho, rho.data(), g, i, j, k);
        if (s > nci_surface::scatter_s_max) continue;
        ASSERT_LT(expected, scatter.s.size());
        EXPECT_NEAR(scatter.s[expected], s, 1e-6);
        EXPECT_NEAR(scatter.signed_rho[expected], signed_rho, 1e-6);
        away_from_isovalue += std::abs(s - params.s_iso) > 0.1;
        ++expected;
    }
    EXPECT_EQ(expected, scatter.eligible_count);
    EXPECT_GT(away_from_isovalue, (size_t)0);

    params.s_iso = 0.9;
    nci_surface::Scatter changed_iso;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, nullptr, &changed_iso));
    ASSERT_EQ(changed_iso.s.size(), scatter.s.size());
    for (size_t i = 0; i < scatter.s.size(); ++i) {
        EXPECT_EQ(changed_iso.s[i], scatter.s[i]);
        EXPECT_EQ(changed_iso.signed_rho[i], scatter.signed_rho[i]);
    }
    params.rho_max = 1e-8;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, nullptr, &changed_iso));
    EXPECT_TRUE(changed_iso.s.empty());
    EXPECT_EQ(changed_iso.eligible_count, (size_t)0);
}

UTEST(nci_surface, scatter_is_capped_and_deterministic) {
    const nci_surface::Grid g = {{42, 42, 42}, {0, 0, 0}, {0.2, 0.2, 0.2}};
    std::vector<float> rho(42 * 42 * 42);
    for (size_t i = 0; i < rho.size(); ++i) rho[i] = 0.01f + (float)(i % 42) * 1e-5f;
    nci_surface::Mesh mesh;
    nci_surface::Scatter a, b;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, {}, nullptr, nullptr, &a));
    EXPECT_EQ(a.eligible_count, (size_t)(40 * 40 * 40));
    ASSERT_EQ(a.s.size(), nci_surface::scatter_max_points);
    ASSERT_EQ(a.signed_rho.size(), a.s.size());
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, {}, nullptr, nullptr, &b));
    ASSERT_EQ(a.s.size(), b.s.size());
    for (size_t i = 0; i < a.s.size(); ++i) {
        EXPECT_EQ(a.s[i], b.s[i]);
        EXPECT_EQ(a.signed_rho[i], b.signed_rho[i]);
        EXPECT_LE(a.s[i], nci_surface::scatter_s_max);
        EXPECT_LE(std::abs(a.signed_rho[i]), 0.05f);
    }
    std::atomic<int32_t> cancel = 1;
    EXPECT_FALSE(nci_surface::extract(&mesh, rho.data(), g, {}, &cancel, nullptr, &a));
    EXPECT_TRUE(a.s.empty());
    EXPECT_TRUE(a.signed_rho.empty());
    EXPECT_EQ(a.eligible_count, (size_t)0);
}
