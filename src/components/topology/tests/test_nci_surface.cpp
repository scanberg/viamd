#include "utest.h"
#include <topology/nci_surface_core.h>
#include <topology/critical_point_distribution_core.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <random>
#include <utility>
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

UTEST(nci_surface, critical_points_group_without_chemical_labels) {
    const nci_surface::CriticalPointAnchor points[] = {
        {{{1.0, 2.0, 3.0}}, nci_surface::Family::HydrogenBond, true},
        {{{4.0, 5.0, 6.0}}, nci_surface::Family::PiStacking, false},
        {{{7.0, 8.0, 9.0}}, nci_surface::Family::Unassigned, true},
    };

    const nci_surface::CriticalPointGroups groups = nci_surface::group_critical_points(points);
    EXPECT_EQ(groups.all_anchors[(size_t)nci_surface::Family::HydrogenBond].size(), (size_t)3);
    EXPECT_EQ(groups.all_anchors[(size_t)nci_surface::Family::PiStacking].size(), (size_t)3);
    EXPECT_EQ(groups.all_anchors[(size_t)nci_surface::Family::Unassigned].size(), (size_t)3);
    EXPECT_EQ(groups.selected_anchors[(size_t)nci_surface::Family::HydrogenBond].size(), (size_t)3);
    EXPECT_TRUE(groups.selected_anchors[(size_t)nci_surface::Family::PiStacking].empty());
    EXPECT_EQ(groups.selected_anchors[(size_t)nci_surface::Family::Unassigned][0], 7.0);
}

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
    params.scatter_s = 1.0;
    nci_surface::Mesh mesh;
    nci_surface::Scatter scatter;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, nullptr, &scatter));
    ASSERT_GT(scatter.s.size(), (size_t)0);
    EXPECT_EQ(scatter.s.size(), scatter.signed_rho.size());
    ASSERT_EQ(scatter.voxel.size(), scatter.s.size());
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
        if (s > params.scatter_s) continue;
        ASSERT_LT(expected, scatter.s.size());
        EXPECT_NEAR(scatter.s[expected], s, 1e-6);
        EXPECT_NEAR(scatter.signed_rho[expected], signed_rho, 1e-6);
        EXPECT_EQ(scatter.voxel[expected], (uint32_t)(i + g.dim[0] * (j + g.dim[1] * k)));
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
        EXPECT_EQ(a.voxel[i], b.voxel[i]);
        EXPECT_LE(a.s[i], nci_surface::scatter_s_max);
        EXPECT_LE(std::abs(a.signed_rho[i]), 0.05f);
    }
    std::atomic<int32_t> cancel = 1;
    EXPECT_FALSE(nci_surface::extract(&mesh, rho.data(), g, {}, &cancel, nullptr, &a));
    EXPECT_TRUE(a.s.empty());
    EXPECT_TRUE(a.signed_rho.empty());
    EXPECT_EQ(a.eligible_count, (size_t)0);
}

UTEST(nci_surface, connected_scope_selects_whole_components_with_optional_cap) {
    Atoms pairs;
    pairs.centre = {{-3, -6, 0}, {3, -6, 0}, {-3, 6, 0}, {3, 6, 0}};
    const auto g = box(6.0, 11.0, 0.2);
    const auto rho = pairs.sample(g);
    nci_surface::Mesh all, chosen, capped;
    ASSERT_TRUE(nci_surface::extract(&all, rho.data(), g, {}));
    ASSERT_GT(all.value.size(), (size_t)0);
    nci_surface::Params params;
    const double anchor[] = {0.07, -6.03, 0.04};
    params.near_xyz = anchor;
    params.connected = true;
    params.radius_cap = false;
    params.near_radius = 0.5;
    nci_surface::Scatter scatter;
    ASSERT_TRUE(nci_surface::extract(&chosen, rho.data(), g, params, nullptr, nullptr, &scatter));
    ASSERT_GT(chosen.index.size(), (size_t)0);
    EXPECT_LT(chosen.value.size(), all.value.size());
    bool beyond_radius = false;
    for (size_t v = 0; v < chosen.value.size(); ++v) {
        EXPECT_LT(chosen.position[3 * v + 1], 0.0f);
        double d2 = 0;
        for (int a = 0; a < 3; ++a) d2 += std::pow(chosen.position[3 * v + a] - anchor[a], 2);
        beyond_radius |= d2 > params.near_radius * params.near_radius;
    }
    EXPECT_TRUE(beyond_radius);
    ASSERT_GT(scatter.s.size(), (size_t)0);
    ASSERT_EQ(scatter.voxel.size(), scatter.s.size());
    for (size_t v = 0; v < scatter.s.size(); ++v) {
        const uint32_t voxel = scatter.voxel[v];
        const int i = voxel % g.dim[0], j = (voxel / g.dim[0]) % g.dim[1], k = voxel / (g.dim[0] * g.dim[1]);
        EXPECT_LT(g.origin[1] + j * g.spacing[1], 0.0);
        double s, r;
        nci_surface::sample(&s, &r, rho.data(), g, i, j, k);
        EXPECT_NEAR(scatter.s[v], s, 1e-6);
        EXPECT_NEAR(scatter.signed_rho[v], r, 1e-6);
    }
    params.radius_cap = true;
    nci_surface::Scatter cap_scatter;
    ASSERT_TRUE(nci_surface::extract(&capped, rho.data(), g, params, nullptr, nullptr, &cap_scatter));
    EXPECT_GT(capped.index.size(), (size_t)0);
    EXPECT_LT(capped.value.size(), chosen.value.size());
    for (uint32_t voxel : cap_scatter.voxel) {
        const int ijk[3] = {(int)(voxel % g.dim[0]), (int)((voxel / g.dim[0]) % g.dim[1]), (int)(voxel / (g.dim[0] * g.dim[1]))};
        double d2 = 0;
        for (int a = 0; a < 3; ++a) d2 += std::pow(g.origin[a] + ijk[a] * g.spacing[a] - anchor[a], 2);
        EXPECT_LE(d2, params.near_radius * params.near_radius + 1e-10);
    }
    const double both[] = {0, -6, 0, 0, 6, 0};
    params.near_xyz = both;
    params.radius_cap = false;
    nci_surface::Mesh united;
    ASSERT_TRUE(nci_surface::extract(&united, rho.data(), g, params));
    const double upper[] = {0, 6, 0};
    params.near_xyz = upper;
    nci_surface::Mesh other;
    ASSERT_TRUE(nci_surface::extract(&other, rho.data(), g, params));
    EXPECT_EQ(united.index.size(), chosen.index.size() + other.index.size());
    EXPECT_EQ(united.value.size(), chosen.value.size() + other.value.size());
    // The global field also contains an unclaimed low-density component between the pairs.
    EXPECT_LT(united.value.size(), all.value.size());
    const double absent[] = {-5, 0, -9};
    params.near_xyz = absent;
    ASSERT_TRUE(nci_surface::extract(&united, rho.data(), g, params, nullptr, nullptr, &scatter));
    EXPECT_TRUE(united.index.empty());
    EXPECT_TRUE(scatter.s.empty());
}

UTEST(nci_surface, scatter_s_widens_the_scatter_without_changing_the_surface) {
    Atoms pair;
    pair.centre = {{-3, 0, 0}, {3, 0, 0}};
    const auto g = box(6.0, 5.0, 0.2);
    const auto rho = pair.sample(g);
    nci_surface::Params params;
    const double anchor[] = {0, 0, 0};
    params.near_xyz = anchor;
    params.connected = true;
    params.near_radius = 4.0;
    params.s_iso = 0.4;
    nci_surface::Mesh base, wide;
    nci_surface::Scatter base_scatter, wide_scatter, same_scatter;
    ASSERT_TRUE(nci_surface::extract(&base, rho.data(), g, params, nullptr, nullptr, &base_scatter));
    ASSERT_GT(base.index.size(), (size_t)0);
    params.scatter_s = 1.0;
    ASSERT_TRUE(nci_surface::extract(&wide, rho.data(), g, params, nullptr, nullptr, &wide_scatter));
    EXPECT_TRUE(wide.index == base.index);
    EXPECT_TRUE(wide.position == base.position);
    EXPECT_GT(wide_scatter.eligible_count, base_scatter.eligible_count);
    EXPECT_GT(wide_scatter.density.total_volume, base_scatter.density.total_volume);
    float base_max = 0, wide_max = 0;
    for (float s : base_scatter.s) base_max = std::max(base_max, s);
    for (float s : wide_scatter.s) wide_max = std::max(wide_max, s);
    EXPECT_GT(wide_max, base_max);
    // Equal to s_iso it is the default
    params.scatter_s = params.s_iso;
    ASSERT_TRUE(nci_surface::extract(&wide, rho.data(), g, params, nullptr, nullptr, &same_scatter));
    EXPECT_EQ(same_scatter.eligible_count, base_scatter.eligible_count);
    EXPECT_TRUE(same_scatter.voxel == base_scatter.voxel);
    params.scatter_s = -1.0;
    EXPECT_FALSE(nci_surface::extract(&wide, rho.data(), g, params));
}

UTEST(nci_surface, scatter_only_matches_surface_extraction) {
    Atoms pair;
    pair.centre = {{-3, 0, 0}, {3, 0, 0}};
    const auto g = box(6.0, 5.0, 0.2);
    const auto rho = pair.sample(g);
    nci_surface::Params params;
    nci_surface::Mesh mesh;
    nci_surface::Scatter a, b;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, nullptr, &a));
    ASSERT_GT(mesh.index.size(), (size_t)0);
    params.surface = false;
    std::atomic<float> progress = 0;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, &progress, &b));
    EXPECT_TRUE(mesh.position.empty());
    EXPECT_TRUE(mesh.index.empty());
    EXPECT_EQ(progress.load(), 1.0f);
    EXPECT_EQ(a.eligible_count, b.eligible_count);
    ASSERT_EQ(a.s.size(), b.s.size());
    for (size_t i = 0; i < a.s.size(); ++i) {
        EXPECT_EQ(a.s[i], b.s[i]);
        EXPECT_EQ(a.signed_rho[i], b.signed_rho[i]);
        EXPECT_EQ(a.voxel[i], b.voxel[i]);
    }

}

UTEST(nci_surface, surface_isovalue_preserves_reference_scatter) {
    Atoms pair;
    pair.centre = {{-3, 0, 0}, {3, 0, 0}};
    const auto g = box(6.0, 5.0, 0.2);
    const auto rho = pair.sample(g);
    const double anchor[] = {0, 0, 0};
    auto params = nci_surface::reference_params(anchor);
    params.surface = true;
    params.s_iso = 0.4;
    nci_surface::Mesh a, b;
    nci_surface::Scatter sa, sb;
    ASSERT_TRUE(nci_surface::extract(&a, rho.data(), g, params, nullptr, nullptr, &sa));
    params.s_iso = 0.7;
    ASSERT_TRUE(nci_surface::extract(&b, rho.data(), g, params, nullptr, nullptr, &sb));
    EXPECT_FALSE(a.position == b.position);
    EXPECT_TRUE(sa.voxel == sb.voxel);
    EXPECT_TRUE(sa.signed_rho == sb.signed_rho);
    EXPECT_TRUE(sa.s == sb.s);
    EXPECT_EQ(sa.eligible_count, sb.eligible_count);
    EXPECT_NEAR(sa.density.total_volume, sb.density.total_volume, 1e-8);
    for (float s : sb.s) EXPECT_LE(s, (float)nci_surface::reference_scatter_s);
}

UTEST(nci_surface, invalid_scope_parameters) {
    const auto g = box(1, 1, 0.2);
    std::vector<float> rho((size_t)g.dim[0] * g.dim[1] * g.dim[2], 0.01f);
    nci_surface::Params params;
    nci_surface::Mesh mesh;
    const double malformed[] = {0, 0};
    params.near_xyz = malformed;
    EXPECT_FALSE(nci_surface::extract(&mesh, rho.data(), g, params));
    const double nonfinite[] = {NAN, 0, 0};
    params.near_xyz = nonfinite;
    EXPECT_FALSE(nci_surface::extract(&mesh, rho.data(), g, params));
    params.near_xyz = {};
    params.s_iso = INFINITY;
    EXPECT_FALSE(nci_surface::extract(&mesh, rho.data(), g, params));
    params.s_iso = 0.5;
    params.near_radius = -1;
    EXPECT_FALSE(nci_surface::extract(&mesh, rho.data(), g, params));
    params.radius_cap = false;
    EXPECT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params));
    params.rho_max = NAN;
    EXPECT_FALSE(nci_surface::extract(&mesh, rho.data(), g, params));
}

UTEST(nci_surface, contact_classes_and_reference_categories) {
    using nci_surface::Family;
    // The CP contact classes (critical_point_distribution's contact_class) each have a type
    const char* classes[] = {"Hydrogen bond", "Halogen bond", "\xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80",
                             "X-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80", "C-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80",
                             "C-H\xC2\xB7\xC2\xB7\xC2\xB7X", "H\xC2\xB7\xC2\xB7\xC2\xB7H", "van der Waals contact",
                             "Ionic hydrogen bond", "Salt bridge", "\xCF\x80-hole bond", "Triel bond", "Aerogen bond",
                             "Hydride bond", "Agostic", "Coordination", "Ion pair", "Metal\xC2\xB7\xC2\xB7\xC2\xB7metal"};
    for (const char* name : classes) {
        bool listed = false;
        for (const auto& c : nci_surface::contact_classes()) listed |= std::string_view(c.name) == name;
        EXPECT_TRUE(listed);
    }
    for (const auto& c : nci_surface::contact_classes()) {
        ASSERT_TRUE(c.description != nullptr);
        EXPECT_GT(std::string_view(c.description).size(), (size_t)10);
        EXPECT_FALSE(critical_point_distribution::contact_environment_name(c.name).empty());
    }
    EXPECT_EQ(nci_surface::contact_family("unknown"), Family::OtherContact);
    EXPECT_EQ(nci_surface::reference_category(Family::HydrogenBond), 0u);
    EXPECT_EQ(nci_surface::reference_category(Family::Ring), 1u);
    EXPECT_EQ(nci_surface::reference_category(Family::PiStacking), 2u);
    EXPECT_EQ(nci_surface::reference_category(Family::HalogenBond), 3u);
    EXPECT_EQ(nci_surface::reference_category(Family::Cage), 3u);
    EXPECT_EQ(nci_surface::reference_category(Family::Bond), nci_surface::reference_category_count);
    EXPECT_EQ(nci_surface::reference_category(Family::Unassigned), nci_surface::reference_category_count);
    for (size_t c = 0; c < nci_surface::reference_category_count; ++c)
        EXPECT_EQ(nci_surface::reference_category(nci_surface::reference_settings_family(c)), c);
    EXPECT_EQ(nci_surface::reference_settings_family(nci_surface::reference_category_count), Family::Unassigned);
}

UTEST(nci_surface, families_follow_chemical_labels_not_character_or_detail) {
    using critical_point_distribution::ChemicalLabel;
    using critical_point_distribution::Kind;
    using nci_surface::Family;
    ChemicalLabel label;
    EXPECT_EQ(nci_surface::family_of(label), Family::Unassigned);
    label.kind = Kind::Contact;
    label.interaction = "Hydrogen bond";
    label.group = "Hydrogen bond (intermolecular): water";
    label.finer = {"O-H...O in water"};
    EXPECT_EQ(nci_surface::family_of(label), Family::HydrogenBond);
    label.interaction = "\xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80";
    EXPECT_EQ(nci_surface::family_of(label), Family::PiStacking);
    label.interaction = "Halogen bond";
    EXPECT_EQ(nci_surface::family_of(label), Family::HalogenBond);
    // Every pi contact shares the pi-stacking type
    label.interaction = "C-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80";
    EXPECT_EQ(nci_surface::family_of(label), Family::PiStacking);
    label.interaction = "X-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80";
    EXPECT_EQ(nci_surface::family_of(label), Family::PiStacking);
    label.interaction = "C-H\xC2\xB7\xC2\xB7\xC2\xB7X";
    EXPECT_EQ(nci_surface::family_of(label), Family::OtherContact);
    // Rings and cages stay rings and cages; closed by contacts, they take their closing interaction's parameters
    label.kind = Kind::Ring;
    EXPECT_EQ(nci_surface::family_of(label), Family::Ring);
    EXPECT_EQ(nci_surface::parameter_family_of(label), Family::OtherContact);
    label.interaction = "C-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80";
    EXPECT_EQ(nci_surface::parameter_family_of(label), Family::PiStacking);
    label.interaction = "\xCF\x80\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80";
    EXPECT_EQ(nci_surface::parameter_family_of(label), Family::PiStacking);
    label.kind = Kind::Cage;
    EXPECT_EQ(nci_surface::family_of(label), Family::Cage);
    EXPECT_EQ(nci_surface::parameter_family_of(label), Family::PiStacking);
    label.interaction = "Hydrogen bond";
    EXPECT_EQ(nci_surface::parameter_family_of(label), Family::HydrogenBond);
    label.interaction.clear();
    EXPECT_EQ(nci_surface::parameter_family_of(label), Family::Cage);
    label.kind = Kind::Ring;
    EXPECT_EQ(nci_surface::parameter_family_of(label), Family::Ring);
    label.kind = Kind::Contact;
    label.interaction = "Hydrogen bond";
    EXPECT_EQ(nci_surface::parameter_family_of(label), Family::HydrogenBond);
    label.kind = Kind::Bond;
    EXPECT_EQ(nci_surface::family_of(label), Family::Bond);
    label.kind = Kind::BondPath;
    EXPECT_EQ(nci_surface::family_of(label), Family::BondPath);
    label.kind = Kind::Atom;
    EXPECT_EQ(nci_surface::family_of(label), Family::Atom);
    label.kind = Kind::NonNuclear;
    EXPECT_EQ(nci_surface::family_of(label), Family::Atom);
    EXPECT_FALSE(nci_surface::has_nci(Family::Bond));
    EXPECT_FALSE(nci_surface::has_nci(Family::BondPath));
    EXPECT_FALSE(nci_surface::has_nci(Family::Atom));
    EXPECT_TRUE(nci_surface::has_nci(Family::Ring));
    EXPECT_TRUE(nci_surface::has_nci(Family::Cage));
    EXPECT_TRUE(nci_surface::has_nci(Family::HalogenBond));
    EXPECT_TRUE(nci_surface::has_nci(Family::Unassigned));
}

UTEST(nci_surface, uniform_defaults_and_mixed_membership) {
    using nci_surface::Family;
    const auto hb = nci_surface::family_bit(Family::HydrogenBond);
    const auto ring = nci_surface::family_bit(Family::Ring);
    ASSERT_TRUE(nci_surface::single_family(hb));
    ASSERT_TRUE(nci_surface::single_family(ring));
    EXPECT_FALSE(nci_surface::single_family(hb | ring));
    EXPECT_FALSE(nci_surface::single_family(0));
    EXPECT_EQ(nci_surface::reference_scatter_s, 1.0);
    EXPECT_EQ(nci_surface::reference_rho_max, 0.07);
    EXPECT_EQ(nci_surface::reference_radius, 2.7);
    const auto uniform = nci_surface::uniform_defaults();
    EXPECT_EQ(uniform.s_iso, 0.7);
    EXPECT_EQ(uniform.scatter_s, 1.0);
    EXPECT_EQ(uniform.rho_max, 0.07);
    EXPECT_EQ(uniform.radius, 2.7);
    EXPECT_EQ(uniform.color_range, 0.03);
}

UTEST(nci_surface, selection_key_changes_update_default_view_state) {
    using nci_surface::View;
    View view = View::Overview;
    uint64_t previous = 0;
    EXPECT_TRUE(nci_surface::follow_selection(&view, &previous, 10, false));
    EXPECT_EQ(view, View::Overview);
    EXPECT_TRUE(nci_surface::follow_selection(&view, &previous, 20, true));
    EXPECT_EQ(view, View::Preview);
    view = View::Saved;
    EXPECT_FALSE(nci_surface::follow_selection(&view, &previous, 20, true));
    EXPECT_EQ(view, View::Saved);
    EXPECT_TRUE(nci_surface::follow_selection(&view, &previous, 30, true));
    EXPECT_EQ(view, View::Preview);
    view = View::Overview;
    EXPECT_FALSE(nci_surface::follow_selection(&view, &previous, 30, true));
    EXPECT_EQ(view, View::Overview);
    EXPECT_TRUE(nci_surface::follow_selection(&view, &previous, 10, false));
    EXPECT_EQ(view, View::Overview);
}

UTEST(nci_surface, reference_parameters_are_uniform) {
    const double anchors[] = {0, 0, 0};
    const auto p = nci_surface::reference_params(anchors);
    EXPECT_EQ(p.scatter_s, nci_surface::reference_scatter_s);
    EXPECT_EQ(p.rho_max, nci_surface::reference_rho_max);
    EXPECT_EQ(p.near_radius, nci_surface::reference_radius);
    EXPECT_TRUE(p.connected);
    EXPECT_TRUE(p.radius_cap);
    EXPECT_TRUE(p.by_anchor);
}

UTEST(nci_surface, family_previews_extract_independent_geometries) {
    Atoms pair;
    pair.centre = {{-3, 0, 0}, {3, 0, 0}};
    const auto g = box(6.0, 5.0, 0.2);
    const auto rho = pair.sample(g);
    const double anchor[] = {0, 0, 0};
    const auto params = [&](double radius) {
        nci_surface::Params p;
        p.s_iso = nci_surface::reference_scatter_s;
        p.rho_max = nci_surface::reference_rho_max;
        p.near_radius = radius;
        // Up to the plot's top, so the radius rather than s bounds the scatter
        p.scatter_s = nci_surface::scatter_s_max;
        p.near_xyz = anchor;
        p.connected = true;
        return p;
    };
    const auto hb = params(1.7);
    const auto stacking = params(3.7);
    nci_surface::Mesh a, b, repeated;
    nci_surface::Scatter sa, sb;
    ASSERT_TRUE(nci_surface::extract(&a, rho.data(), g, hb, nullptr, nullptr, &sa));
    ASSERT_TRUE(nci_surface::extract(&b, rho.data(), g, stacking, nullptr, nullptr, &sb));
    ASSERT_GT(a.index.size(), (size_t)0);
    ASSERT_GT(b.index.size(), (size_t)0);
    // Same s and rho, different radius: the wider pi-stacking radius reaches samples the hydrogen-bond one does not
    EXPECT_GT(sb.s.size(), sa.s.size());
    for (const auto& samples : {&sa, &sb}) {
        const double radius = samples == &sa ? hb.near_radius : stacking.near_radius;
        for (uint32_t v : samples->voxel) {
            const int ijk[] = {(int)(v % g.dim[0]), (int)((v / g.dim[0]) % g.dim[1]), (int)(v / (g.dim[0] * g.dim[1]))};
            double d2 = 0;
            for (int axis = 0; axis < 3; ++axis) d2 += std::pow(g.origin[axis] + ijk[axis] * g.spacing[axis], 2);
            EXPECT_LE(d2, radius * radius + 1e-10);
        }
    }
    ASSERT_TRUE(nci_surface::extract(&repeated, rho.data(), g, hb));
    ASSERT_EQ(a.position.size(), repeated.position.size());
    for (size_t i = 0; i < a.position.size(); ++i) EXPECT_EQ(a.position[i], repeated.position[i]);
    EXPECT_EQ(hb.s_iso, stacking.s_iso);
    EXPECT_LT(hb.near_radius, stacking.near_radius);
}

UTEST(nci_surface, scatter_density_contour_tracks_the_selected_sample_cloud) {
    nci_surface::Scatter scatter;
    scatter.signed_rho.assign(100, -0.015f);
    scatter.s.assign(100, 0.5f);
    std::vector<std::array<float, 4>> outline, repeated;
    ASSERT_TRUE(nci_surface::scatter_density_outline(&outline, scatter, 0.05f));
    ASSERT_GT(outline.size(), (size_t)4);
    bool diagonal = false;
    for (const auto& line : outline) {
        for (int v : {0, 2}) {
            EXPECT_LT(line[v], -0.009f);
            EXPECT_GT(line[v], -0.021f);
            EXPECT_GT(line[v + 1], 0.35f);
            EXPECT_LT(line[v + 1], 0.65f);
        }
        diagonal |= std::abs(line[0] - line[2]) > 1e-8f && std::abs(line[1] - line[3]) > 1e-8f;
    }
    EXPECT_TRUE(diagonal);
    ASSERT_TRUE(nci_surface::scatter_density_outline(&repeated, scatter, 0.05f));
    ASSERT_EQ(outline.size(), repeated.size());
    for (size_t i = 0; i < outline.size(); ++i)
        for (int a = 0; a < 4; ++a) EXPECT_EQ(outline[i][a], repeated[i][a]);
    scatter = {};
    ASSERT_TRUE(nci_surface::scatter_density_outline(&outline, scatter, 0.05f));
    EXPECT_TRUE(outline.empty());
    scatter.s.push_back(1);
    EXPECT_FALSE(nci_surface::scatter_density_outline(&outline, scatter, 0.05f));
    scatter.signed_rho.push_back(NAN);
    EXPECT_FALSE(nci_surface::scatter_density_outline(&outline, scatter, 0.05f));
    EXPECT_FALSE(nci_surface::scatter_density_outline(&outline, {}, 0.0f));
}

UTEST(nci_surface, union_outline_encloses_each_scatter_outline_with_a_margin) {
    // A dense cloud and a sparse one: each at its own 90%-of-volume level, so the sparse one is not lost
    nci_surface::Scatter dense, sparse;
    dense.signed_rho.assign(1000, -0.02f);
    dense.s.assign(1000, 0.5f);
    sparse.signed_rho.assign(10, 0.02f);
    sparse.s.assign(10, 0.8f);
    const auto bounds = [](const std::vector<std::array<float, 4>>& lines) {
        std::array<float, 4> b = {1e9f, 1e9f, -1e9f, -1e9f};
        for (const auto& l : lines) for (int v : {0, 2}) {
            b[0] = std::min(b[0], l[v]); b[1] = std::min(b[1], l[v + 1]);
            b[2] = std::max(b[2], l[v]); b[3] = std::max(b[3], l[v + 1]);
        }
        return b;
    };
    std::vector<std::array<float, 4>> outline, merged;
    const nci_surface::Scatter* both[] = {&dense, &sparse};
    ASSERT_TRUE(nci_surface::scatter_union_outline(&merged, both, 0.05f, 1));
    const auto u = bounds(merged);
    for (const nci_surface::Scatter* scatter : both) {
        ASSERT_TRUE(nci_surface::scatter_density_outline(&outline, *scatter, 0.05f));
        const auto o = bounds(outline);
        EXPECT_LT(u[0], o[0]);
        EXPECT_LT(u[1], o[1]);
        EXPECT_GT(u[2], o[2]);
        EXPECT_GT(u[3], o[3]);
    }
    const nci_surface::Scatter* none[] = {&dense, nullptr};
    EXPECT_FALSE(nci_surface::scatter_union_outline(&merged, none, 0.05f, 1));
    ASSERT_TRUE(nci_surface::scatter_union_outline(&merged, {}, 0.05f, 1));
    EXPECT_TRUE(merged.empty());
}

UTEST(nci_surface, expected_region_edge_holds_ninety_percent_of_the_samples_and_the_core_half) {
    // A broad Gaussian cloud: the edge holds 90% of its samples whatever its peak height
    nci_surface::Scatter cloud;
    std::mt19937 rng(7);
    std::normal_distribution<float> rho(-0.01f, 0.006f), s(0.8f, 0.2f);
    for (int i = 0; i < 20000; ++i) {
        cloud.signed_rho.push_back(rho(rng));
        cloud.s.push_back(s(rng));
    }
    const nci_surface::Scatter* samples[] = {&cloud};
    const auto region = nci_surface::expected_region(samples, 0.05f, 0);
    ASSERT_TRUE(region.valid);
    EXPECT_GT(region.core_outline.size(), (size_t)4);
    const auto counts = nci_surface::expected_region_counts(cloud, region);
    const double inside = (double)counts.inside / (counts.inside + counts.outside);
    EXPECT_GT(inside, 0.87);
    EXPECT_LT(inside, 0.95);
    // The core lies inside the edge
    for (const auto& line : region.core_outline)
        EXPECT_TRUE(nci_surface::in_expected_region(region, line[0], line[1]));
}

UTEST(nci_surface, expected_regions_are_irregular_and_derived_from_reference_samples) {
    nci_surface::Scatter hydrogen_bonds;
    for (int i = 0; i < 100; ++i) {
        hydrogen_bonds.signed_rho.push_back(-0.025f + (i % 5) * 0.0004f);
        hydrogen_bonds.s.push_back(0.35f + (i / 5) * 0.004f);
    }
    nci_surface::Scatter reference;
    const nci_surface::Scatter* samples[] = {&hydrogen_bonds};
    const auto region = nci_surface::expected_region(samples, 0.05f);
    ASSERT_TRUE(region.valid);
    EXPECT_EQ(region.mask.size(), (size_t)nci_surface::scatter_field_bins * nci_surface::scatter_field_bins);
    EXPECT_EQ(region.density.size(), (size_t)nci_surface::scatter_field_bins * nci_surface::scatter_field_bins);
    EXPECT_GT(region.outline.size(), (size_t)4);
    EXPECT_TRUE(nci_surface::in_expected_region(region, -0.024f, 0.39f));
    EXPECT_FALSE(nci_surface::in_expected_region(region, 0.025f, 1.5f));
    EXPECT_FALSE(nci_surface::in_expected_region(region, -0.024f, 2.1f));

    reference.signed_rho = {-0.024f, -0.024f, 0.025f, -0.024f};
    reference.s = {0.39f, 0.39f, 1.5f, 2.1f};
    const auto counts = nci_surface::expected_region_counts(reference, region);
    EXPECT_EQ(counts.inside, 2u);
    EXPECT_EQ(counts.outside, 1u);  // Samples above scatter_s_max are ineligible
    std::vector<std::array<float, 2>> fill;
    nci_surface::expected_region_fill(&fill, region);
    EXPECT_FALSE(fill.empty());
    EXPECT_EQ(fill.size() % 3, 0u);
    EXPECT_FALSE(nci_surface::expected_region({}, 0.05f).valid);
}

UTEST(nci_surface, union_fill_lies_inside_the_union_outline) {
    nci_surface::Scatter cloud;
    cloud.signed_rho.assign(500, -0.02f);
    cloud.s.assign(500, 0.5f);
    const nci_surface::Scatter* one[] = {&cloud};
    std::vector<std::array<float, 4>> outline;
    std::vector<std::array<float, 2>> fill;
    ASSERT_TRUE(nci_surface::scatter_union_outline(&outline, one, 0.05f, 0));
    ASSERT_TRUE(nci_surface::scatter_union_fill(&fill, one, 0.05f, 0));
    ASSERT_FALSE(fill.empty());
    EXPECT_EQ(fill.size() % 3, 0u);
    float lo_r = 1e9f, lo_s = 1e9f, hi_r = -1e9f, hi_s = -1e9f;
    for (const auto& l : outline) for (int v : {0, 2}) {
        lo_r = std::min(lo_r, l[v]); hi_r = std::max(hi_r, l[v]);
        lo_s = std::min(lo_s, l[v + 1]); hi_s = std::max(hi_s, l[v + 1]);
    }
    double area = 0.0;
    for (size_t i = 0; i < fill.size(); i += 3) {
        for (size_t k = 0; k < 3; ++k) {
            EXPECT_TRUE(fill[i + k][0] >= lo_r - 1e-6f && fill[i + k][0] <= hi_r + 1e-6f);
            EXPECT_TRUE(fill[i + k][1] >= lo_s - 1e-6f && fill[i + k][1] <= hi_s + 1e-6f);
        }
        const auto& a = fill[i]; const auto& b = fill[i + 1]; const auto& c = fill[i + 2];
        area += 0.5 * std::abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]));
    }
    EXPECT_GT(area, 0.0);
}

UTEST(nci_surface, heatmap_peaks_at_one_where_the_samples_are) {
    nci_surface::Scatter a, b;
    a.signed_rho.assign(100, -0.02f);
    a.s.assign(100, 0.5f);
    b.signed_rho.assign(10, 0.03f);
    b.s.assign(10, 1.5f);
    const nci_surface::Scatter* both[] = {&a, &b};
    nci_surface::Heatmap map;
    ASSERT_TRUE(nci_surface::scatter_heatmap(&map, both, 0.05f));
    ASSERT_EQ(map.values.size(), (size_t)map.bins * map.bins);
    const auto bin = [&](float rho, float s) {
        const int x = (int)((rho + 0.05f) / 0.1f * map.bins), y = (int)(s / nci_surface::scatter_s_max * map.bins);
        return map.values[x + (size_t)map.bins * y];
    };
    EXPECT_NEAR(bin(-0.02f, 0.5f), 1.0f, 1e-6f);
    EXPECT_GT(bin(0.03f, 1.5f), 0.0f);
    EXPECT_LT(bin(0.03f, 1.5f), 0.5f);
    EXPECT_EQ(*std::max_element(map.values.begin(), map.values.end()), 1.0f);
}

namespace {
nci_surface::Scatter context_scatter(const Atoms& atoms, double h) {
    const auto g = box(5.5, 3.5, h);
    const auto rho = atoms.sample(g);
    nci_surface::Params params;
    params.rho_max = 0.05;
    params.radius_cap = false;
    params.surface = false;
    nci_surface::Mesh mesh;
    nci_surface::Scatter scatter;
    if (!nci_surface::extract(&mesh, rho.data(), g, params, nullptr, nullptr, &scatter)) return {};
    return scatter;
}
double heatmap_difference(const nci_surface::Scatter& a, const nci_surface::Scatter& b) {
    nci_surface::Heatmap ma, mb;
    const nci_surface::Scatter* sa[] = {&a};
    const nci_surface::Scatter* sb[] = {&b};
    if (!nci_surface::scatter_heatmap(&ma, sa, 0.05f) || !nci_surface::scatter_heatmap(&mb, sb, 0.05f)) return INFINITY;
    double diff = 0.0, sum = 0.0;
    for (size_t i = 0; i < ma.values.size(); ++i) {
        diff += std::abs(ma.values[i] - mb.values[i]);
        sum += ma.values[i] + mb.values[i];
    }
    return sum > 0.0 ? 2.0 * diff / sum : INFINITY;
}
}  // namespace

UTEST(nci_surface, continuous_density_covers_the_eligible_volume_and_converges) {
    Atoms pair;
    pair.centre = {{-2.6, 0.0, 0.0}, {2.6, 0.0, 0.0}};
    const auto coarse = context_scatter(pair, 0.2);
    const auto fine = context_scatter(pair, 0.1);
    ASSERT_EQ(coarse.density.volume.size(), (size_t)nci_surface::scatter_density_bins * nci_surface::scatter_density_bins);
    ASSERT_EQ(fine.density.volume.size(), coarse.density.volume.size());
    EXPECT_GT(fine.density.cells, (size_t)0);
    // The same region as the eligible grid samples, each standing for its voxel
    EXPECT_NEAR(coarse.density.total_volume, coarse.eligible_count * 0.008, 0.1 * coarse.eligible_count * 0.008);
    EXPECT_NEAR(fine.density.total_volume, fine.eligible_count * 0.001, 0.1 * fine.eligible_count * 0.001);
    EXPECT_NEAR(coarse.density.total_volume, fine.density.total_volume, 0.1 * fine.density.total_volume);

    auto coarse_points = coarse, fine_points = fine;
    coarse_points.density = {};
    fine_points.density = {};
    const double continuous = heatmap_difference(coarse, fine);
    const double points = heatmap_difference(coarse_points, fine_points);
    printf("continuous density change h 0.2 -> 0.1: %.3f, point histogram: %.3f\n", continuous, points);
    EXPECT_LT(continuous, points);
}

UTEST(nci_surface, continuous_density_respects_scope_and_cancel) {
    Atoms pair;
    pair.centre = {{-2.6, 0.0, 0.0}, {2.6, 0.0, 0.0}};
    const auto g = box(5.5, 3.5, 0.2);
    const auto rho = pair.sample(g);
    const double anchor[3] = {0.0, 0.0, 0.0};
    nci_surface::Params params;
    params.rho_max = 0.05;
    params.near_xyz = anchor;
    params.near_radius = 1.0;
    params.surface = false;
    nci_surface::Mesh mesh;
    nci_surface::Scatter scatter;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, nullptr, &scatter));
    EXPECT_GT(scatter.density.total_volume, 0.0);
    EXPECT_LT(scatter.density.total_volume, 4.0 / 3.0 * 3.15 * 1.3 * 1.3 * 1.3);
    // A contact at signed rho < 0 only: no volume on the positive side near the cutoff's middle
    double negative = 0.0, positive = 0.0;
    constexpr int bins = nci_surface::scatter_density_bins;
    for (int y = 0; y < bins; ++y) for (int x = 0; x < bins; ++x)
        (x < bins / 2 ? negative : positive) += scatter.density.volume[x + (size_t)bins * y];
    EXPECT_GT(negative, 0.0);

    std::atomic<int32_t> cancel = 1;
    EXPECT_FALSE(nci_surface::extract(&mesh, rho.data(), g, params, &cancel, nullptr, &scatter));
}

UTEST(nci_surface, scatter_keeps_both_signs_below_a_hard_s_cut) {
    // A three-atom ring: bond contacts with lambda2 < 0 on its edges, the ring's lambda2 > 0 region inside
    Atoms ring;
    const double r = 3.0;
    for (int a = 0; a < 3; ++a)
        ring.centre.push_back({r * std::cos(2.0944 * a), r * std::sin(2.0944 * a), 0.0});
    const auto g = box(5.5, 5.5, 0.15);
    const auto rho = ring.sample(g);
    const double bond[3] = {0.5 * (ring.centre[0][0] + ring.centre[1][0]), 0.5 * (ring.centre[0][1] + ring.centre[1][1]), 0.0};
    nci_surface::Mesh mesh;
    nci_surface::Scatter scatter;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, nci_surface::reference_params(bond), nullptr, nullptr, &scatter));
    size_t positive = 0, negative = 0;
    for (size_t i = 0; i < scatter.s.size(); ++i) {
        EXPECT_LE(scatter.s[i], (float)nci_surface::reference_scatter_s);
        (scatter.signed_rho[i] > 0.0f ? positive : negative)++;
    }
    EXPECT_GT(positive, (size_t)0);
    EXPECT_GT(negative, (size_t)0);
    // No density above the cut either
    constexpr int bins = nci_surface::scatter_density_bins;
    const int cut = (int)std::ceil(nci_surface::reference_scatter_s / nci_surface::scatter_s_max * bins);
    double above = 0.0;
    for (int y = cut; y < bins; ++y) for (int x = 0; x < bins; ++x) above += scatter.density.volume[x + (size_t)bins * y];
    EXPECT_EQ(above, 0.0);
    EXPECT_GT(scatter.density.total_volume, 0.0);
}

UTEST(nci_surface, scatter_reaches_the_plot_s_range_beyond_the_reference) {
    // The settings may raise the scatter's s past the reference's to show the tail of the plot up to s = 2
    Atoms ring;
    const double r = 3.0;
    for (int a = 0; a < 3; ++a)
        ring.centre.push_back({r * std::cos(2.0944 * a), r * std::sin(2.0944 * a), 0.0});
    const auto g = box(5.5, 5.5, 0.15);
    const auto rho = ring.sample(g);
    const double bond[3] = {0.5 * (ring.centre[0][0] + ring.centre[1][0]), 0.5 * (ring.centre[0][1] + ring.centre[1][1]), 0.0};
    nci_surface::Params params = nci_surface::reference_params(bond);
    params.scatter_s = nci_surface::scatter_s_max;
    params.surface = false;
    nci_surface::Mesh mesh;
    nci_surface::Scatter scatter;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, nullptr, &scatter));
    size_t tail = 0;
    for (float s : scatter.s) {
        EXPECT_LE(s, nci_surface::scatter_s_max);
        tail += s > (float)nci_surface::reference_scatter_s;
    }
    EXPECT_GT(tail, (size_t)0);
}

UTEST(nci_surface, display_tail_leaves_the_reference_scatter_untouched) {
    // Points above scatter_s up to display_s go to the tail only; the scatter and its density stay the reference's
    Atoms ring;
    const double r = 3.0;
    for (int a = 0; a < 3; ++a)
        ring.centre.push_back({r * std::cos(2.0944 * a), r * std::sin(2.0944 * a), 0.0});
    const auto g = box(5.5, 5.5, 0.15);
    const auto rho = ring.sample(g);
    const double bond[3] = {0.5 * (ring.centre[0][0] + ring.centre[1][0]), 0.5 * (ring.centre[0][1] + ring.centre[1][1]), 0.0};
    nci_surface::Mesh mesh;
    nci_surface::Scatter plain, tailed;
    nci_surface::Params params = nci_surface::reference_params(bond);
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, nullptr, &plain));
    EXPECT_TRUE(plain.tail.s.empty());
    params.display_s = nci_surface::scatter_s_max;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, params, nullptr, nullptr, &tailed));
    EXPECT_TRUE(tailed.s == plain.s);
    EXPECT_TRUE(tailed.signed_rho == plain.signed_rho);
    EXPECT_EQ(tailed.eligible_count, plain.eligible_count);
    EXPECT_TRUE(tailed.density.volume == plain.density.volume);
    ASSERT_FALSE(tailed.tail.s.empty());
    EXPECT_EQ(tailed.tail.s.size(), tailed.tail.signed_rho.size());
    EXPECT_EQ(tailed.tail.s.size(), tailed.tail.anchor.size());
    for (float s : tailed.tail.s) {
        EXPECT_GT(s, (float)nci_surface::reference_scatter_s);
        EXPECT_LE(s, nci_surface::scatter_s_max);
    }
}

UTEST(nci_surface, nearest_anchor_partitions_the_scatter) {
    Atoms ring;
    const double r = 3.0;
    for (int a = 0; a < 3; ++a)
        ring.centre.push_back({r * std::cos(2.0944 * a), r * std::sin(2.0944 * a), 0.0});
    const auto g = box(5.5, 5.5, 0.15);
    const auto rho = ring.sample(g);
    // The three bond midpoints and the ring centre
    std::vector<double> anchors;
    for (int a = 0; a < 3; ++a) {
        const auto& p = ring.centre[a];
        const auto& q = ring.centre[(a + 1) % 3];
        anchors.insert(anchors.end(), {0.5 * (p[0] + q[0]), 0.5 * (p[1] + q[1]), 0.0});
    }
    anchors.insert(anchors.end(), {0.0, 0.0, 0.0});
    nci_surface::Mesh mesh;
    nci_surface::Scatter all;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, nci_surface::reference_params(anchors), nullptr, nullptr, &all));
    ASSERT_EQ(all.anchor.size(), all.s.size());
    ASSERT_FALSE(all.density.by_anchor.empty());
    for (uint32_t a : all.anchor) EXPECT_LT(a, 4u);

    // The subsets of a partition add up to the whole
    nci_surface::Scatter sum;
    for (uint8_t a = 0; a < 4; ++a) {
        std::vector<uint8_t> keep(4, 0);
        keep[a] = 1;
        const auto part = nci_surface::scatter_subset(all, keep);
        for (uint32_t b : part.anchor) EXPECT_EQ(b, (uint32_t)a);
        EXPECT_GT(part.density.total_volume, 0.0);
        nci_surface::scatter_append(&sum, part);
    }
    EXPECT_EQ(sum.s.size(), all.s.size());
    EXPECT_NEAR(sum.density.total_volume, all.density.total_volume, 1e-6 * all.density.total_volume);
    for (size_t b = 0; b < all.density.volume.size(); ++b)
        EXPECT_NEAR(sum.density.volume[b], all.density.volume[b], 1e-3f * (all.density.volume[b] + 1e-6f));

    // With the others as competitors, one anchor's own extraction is its share of the partition
    const std::span<const double> first(anchors.data(), 3), rest(anchors.data() + 3, anchors.size() - 3);
    nci_surface::Scatter own;
    ASSERT_TRUE(nci_surface::extract(&mesh, rho.data(), g, nci_surface::reference_params(first, rest), nullptr, nullptr, &own));
    const std::vector<uint8_t> keep_first = {1, 0, 0, 0};
    const auto share = nci_surface::scatter_subset(all, keep_first);
    EXPECT_GT(own.density.total_volume, 0.0);
    // The connected growth from one seed may miss parts reached only through the others' regions
    EXPECT_LE(own.density.total_volume, share.density.total_volume * (1.0 + 1e-6));
    EXPECT_GT(own.density.total_volume, 0.5 * share.density.total_volume);
}

UTEST(nci_surface, builtin_reference_map_has_nested_levels_per_category) {
    const auto& data = nci_surface::reference_data();
    ASSERT_EQ(data.bins, nci_surface::scatter_field_bins);
    const auto& map = nci_surface::builtin_reference_map();
    ASSERT_TRUE(map.valid);
    EXPECT_EQ(map.rho_extent, data.rho_extent);
    for (size_t c = 0; c < nci_surface::reference_category_count; ++c) {
        const auto& category = map.categories[c];
        ASSERT_TRUE(category.valid);
        EXPECT_GT(category.systems, 0);
        double sum = 0.0;
        for (int i = 0; i < data.bins * data.bins; ++i) sum += data.fields[c][i];
        EXPECT_NEAR(sum, 1.0, 1e-3);
        // The levels increase in what they enclose: 98 % is the outermost, lowest threshold; 90 % the inside one
        ASSERT_EQ(category.thresholds.size(), nci_surface::reference_default_fractions.size());
        EXPECT_GT(category.thresholds[2], 0.0f);
        EXPECT_GE(category.thresholds[0], category.thresholds[1]);
        EXPECT_GE(category.thresholds[1], category.thresholds[2]);
        EXPECT_EQ(category.inside_threshold, category.thresholds[1]);
        for (const auto& isoline : category.isolines) EXPECT_FALSE(isoline.empty());
    }
    // Hydrogen bonds are attractive: most of their volume is at negative signed rho (both signs are sampled, so the
    // repulsive side next to the bond holds the rest), and they hold the largest share of all categories at strongly
    // negative signed rho. Their peak is taken on the attractive side; the far corner is outside every category.
    const auto attractive = [&](size_t c, int columns) {
        double sum = 0.0;
        for (int y = 0; y < data.bins; ++y)
            for (int x = 0; x < columns; ++x) sum += data.fields[c][x + data.bins * y];
        return sum;
    };
    EXPECT_GT(attractive(0, data.bins / 2), 0.5);
    for (size_t c = 1; c < nci_surface::reference_category_count; ++c)
        EXPECT_GT(attractive(0, data.bins / 4), attractive(c, data.bins / 4));
    const auto& hb = data.fields[0];
    int peak = 0;
    for (int i = 1; i < data.bins * data.bins; ++i)
        if (i % data.bins < data.bins / 2 && hb[i] > hb[peak]) peak = i;
    const float peak_rho = -data.rho_extent + 2.0f * data.rho_extent * ((peak % data.bins) + 0.5f) / data.bins;
    const float peak_s = data.s_max * ((peak / data.bins) + 0.5f) / data.bins;
    EXPECT_LT(peak_rho, 0.0f);
    EXPECT_GE(nci_surface::reference_level(map, 0, peak_rho, peak_s), 2);
    for (size_t c = 0; c < nci_surface::reference_category_count; ++c)
        EXPECT_EQ(nci_surface::reference_level(map, c, data.rho_extent * 0.99f, 1.98f), 0);
    EXPECT_EQ(nci_surface::reference_level(map, 0, 1.0f, 0.5f), 0);
    EXPECT_EQ(nci_surface::reference_level(map, 0, NAN, 0.5f), 0);
    // At the hydrogen-bond peak, hydrogen bonds are the most likely category; the shares sum to one, or are all 0 off
    // the map
    const auto shares = nci_surface::reference_shares(map, peak_rho, peak_s);
    float total = 0.0f;
    for (float share : shares) total += share;
    EXPECT_NEAR(total, 1.0f, 1e-5f);
    for (size_t c = 1; c < shares.size(); ++c) EXPECT_GT(shares[0], shares[c]);
    for (float share : nci_surface::reference_shares(map, 1.0f, 0.5f)) EXPECT_EQ(share, 0.0f);
}

UTEST(nci_surface, reference_map_levels_counts_and_dominant_image) {
    // Two synthetic categories: a narrow peak left of centre and a broad one right of it, on the scatter plot bins
    constexpr int bins = nci_surface::scatter_field_bins;
    std::vector<float> narrow(bins * bins), broad(bins * bins), empty(bins * bins, 0.0f);
    for (int y = 0; y < bins; ++y) for (int x = 0; x < bins; ++x) {
        const auto g = [&](float cx, float cy, float w) {
            return std::exp(-((x - cx) * (x - cx) + (y - cy) * (y - cy)) / (2 * w * w));
        };
        narrow[x + bins * y] = g(30, 48, 3);
        broad[x + bins * y] = g(60, 48, 10);
    }
    for (auto* f : {&narrow, &broad}) {
        double sum = 0;
        for (float v : *f) sum += v;
        for (float& v : *f) v = (float)(v / sum);
    }
    nci_surface::ReferenceData data;
    data.bins = bins;
    data.rho_extent = 0.05f;
    data.s_max = nci_surface::scatter_s_max;
    data.names = {"a", "b", "c", "d"};
    data.systems = {1, 2, 0, 0};
    data.fields = {narrow.data(), broad.data(), empty.data(), nullptr};
    nci_surface::ReferenceMap map;
    ASSERT_TRUE(nci_surface::build_reference_map(&map, data, nci_surface::reference_default_fractions));
    EXPECT_TRUE(map.categories[0].valid);
    EXPECT_TRUE(map.categories[1].valid);
    EXPECT_FALSE(map.categories[2].valid);
    EXPECT_FALSE(map.categories[3].valid);
    const auto rho_at = [&](float x) { return -data.rho_extent + 2.0f * data.rho_extent * x / (bins + 1); };
    const auto s_at = [&](float y) { return nci_surface::scatter_s_max * y / (bins + 1); };
    EXPECT_EQ(nci_surface::reference_level(map, 0, rho_at(31), s_at(49)), 3);
    EXPECT_EQ(nci_surface::reference_level(map, 0, rho_at(61), s_at(49)), 0);
    EXPECT_EQ(nci_surface::reference_level(map, 1, rho_at(61), s_at(49)), 3);
    EXPECT_EQ(nci_surface::reference_level(map, 2, rho_at(61), s_at(49)), 0);

    // Counts: inside the 90 % area
    nci_surface::Scatter scatter;
    scatter.signed_rho = {rho_at(31), rho_at(31), rho_at(80), NAN};
    scatter.s = {s_at(49), s_at(50), s_at(10), 0.5f};
    const auto counts = nci_surface::reference_counts(scatter, map, 0);
    EXPECT_EQ(counts.inside, (size_t)2);
    EXPECT_EQ(counts.outside, (size_t)1);
    EXPECT_EQ(nci_surface::reference_counts(scatter, map, 2).inside + nci_surface::reference_counts(scatter, map, 2).outside, (size_t)0);

    // Each category dominates its own peak; the corners are empty
    nci_surface::ReferenceImage image;
    nci_surface::reference_image(&image, map, 97);
    ASSERT_EQ(image.category.size(), (size_t)97 * 97);
    EXPECT_EQ(image.category[31 + 97 * 49], 0);
    EXPECT_EQ(image.level[31 + 97 * 49], 3);
    EXPECT_EQ(image.category[61 + 97 * 49], 1);
    EXPECT_EQ(image.category[0], -1);
    EXPECT_EQ(image.level[0], 0);
    // The combined map: depth 1 in a category's core, 0 outside it, graded from the edge depth in between
    ASSERT_EQ(image.depth.size(), (size_t)97 * 97 * nci_surface::reference_category_count);
    const auto depth = [&](int i, int j, size_t c) { return image.depth[((size_t)i + 97 * j) * nci_surface::reference_category_count + c]; };
    EXPECT_EQ(depth(31, 49, 0), 1.0f);
    EXPECT_EQ(depth(61, 49, 0), 0.0f);
    EXPECT_EQ(depth(0, 0, 0), 0.0f);
    for (float d : image.depth) EXPECT_TRUE(d == 0.0f || (d >= nci_surface::reference_edge_depth && d <= 1.0f));

    // Other levels: their own thresholds and isolines; the inside region stays the 90 % one
    const float one_level[] = {0.7f};
    nci_surface::ReferenceMap single;
    ASSERT_TRUE(nci_surface::build_reference_map(&single, data, one_level));
    ASSERT_EQ(single.categories[0].isolines.size(), (size_t)1);
    EXPECT_EQ(single.categories[0].inside_threshold, map.categories[0].inside_threshold);
    EXPECT_EQ(nci_surface::reference_level(single, 0, rho_at(31), s_at(49)), 1);
    EXPECT_EQ(nci_surface::reference_counts(scatter, single, 0).inside, counts.inside);
    nci_surface::reference_image(&image, single, 97);
    EXPECT_EQ(depth(31, 49, 0), 1.0f);
    const float unsorted[] = {0.9f, 0.5f}, out_of_range[] = {0.5f, 1.0f};
    EXPECT_FALSE(nci_surface::build_reference_map(&single, data, unsorted));
    EXPECT_FALSE(nci_surface::build_reference_map(&single, data, out_of_range));
    EXPECT_FALSE(nci_surface::build_reference_map(&single, data, std::span<const float>()));

    data.bins = bins / 2;
    EXPECT_FALSE(nci_surface::build_reference_map(&map, data, nci_surface::reference_default_fractions));
    EXPECT_FALSE(map.valid);
}

UTEST(nci_surface, reference_subfamilies_nest_in_their_categories) {
    const nci_surface::ReferenceMap& map = nci_surface::builtin_reference_map();
    ASSERT_TRUE(map.valid);
    ASSERT_FALSE(map.subfamilies.empty());
    const nci_surface::ReferenceData& data = nci_surface::reference_data();
    for (size_t i = 0; i < map.subfamilies.size(); ++i) {
        const nci_surface::ReferenceCategory& sub = map.subfamilies[i];
        ASSERT_LT(sub.category, map.categories.size());
        EXPECT_TRUE(sub.valid);
        EXPECT_GE(sub.systems, 5);
        ASSERT_EQ(map.fractions.size(), (size_t)3);
        EXPECT_NEAR(map.fractions[1], 0.9f, 1e-6f);
        ASSERT_EQ(sub.isolines.size(), map.fractions.size());
        EXPECT_FALSE(sub.isolines[1].empty());
        EXPECT_TRUE(sub.overlap > 0.0f && sub.overlap <= 1.0f);
        // Somewhere it is at its deepest level, and there its category's subfamilies share the point
        bool found = false;
        for (int y = 0; y < data.bins && !found; ++y) for (int x = 0; x < data.bins && !found; ++x) {
            const float rho = data.rho_extent * (2.0f * (x + 0.5f) / data.bins - 1.0f);
            const float s = data.s_max * (y + 0.5f) / data.bins;
            if (nci_surface::reference_subfamily_level(map, i, rho, s) < (int)map.fractions.size()) continue;
            found = true;
            const std::vector<float> shares = nci_surface::reference_subfamily_shares(map, sub.category, rho, s);
            ASSERT_EQ(shares.size(), map.subfamilies.size());
            float sum = 0.0f;
            for (size_t j = 0; j < shares.size(); ++j) {
                if (map.subfamilies[j].category != sub.category) EXPECT_EQ(shares[j], 0.0f);
                sum += shares[j];
            }

            EXPECT_NEAR(sum, 1.0f, 1e-4f);
            if (!sub.parent[0]) EXPECT_GT(shares[i], 0.0f);
            else EXPECT_EQ(shares[i], 0.0f);
        }
        EXPECT_TRUE(found);
    }
}

UTEST(nci_surface, reference_refinements_cover_every_category) {
    const auto& map = nci_surface::builtin_reference_map();
    ASSERT_GT(map.subfamilies.size(), (size_t)64);
    std::array<bool, nci_surface::reference_category_count> refined = {};
    bool nh_o = false;
    for (const auto& sub : map.subfamilies) {
        if (!sub.parent[0]) continue;
        bool parent_found = false;
        for (const auto& parent : map.subfamilies)
            parent_found |= parent.category == sub.category && !parent.parent[0] &&
                            std::string_view(parent.name) == sub.parent;
        EXPECT_TRUE(parent_found);
        refined[sub.category] = true;
        nh_o |= std::string_view(sub.name).find("N-H\xC2\xB7\xC2\xB7\xC2\xB7O") != std::string_view::npos;
    }
    for (bool available : refined) EXPECT_TRUE(available);
    EXPECT_TRUE(nh_o);
}
