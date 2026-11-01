#include "utest.h"
#include <topology/nci_surface_core.h>
#include <topology/critical_point_distribution_core.h>

#include <algorithm>
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
        if (s > nci_surface::scatter_s_max) continue;
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
    label.interaction = "C-H\xC2\xB7\xC2\xB7\xC2\xB7\xCF\x80";
    EXPECT_EQ(nci_surface::family_of(label), Family::OtherContact);
    // Rings and cages stay rings and cages; closed by contacts, they take their closing interaction's parameters
    label.kind = Kind::Ring;
    EXPECT_EQ(nci_surface::family_of(label), Family::Ring);
    EXPECT_EQ(nci_surface::parameter_family_of(label), Family::OtherContact);
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

UTEST(nci_surface, calibrated_family_defaults_and_mixed_membership) {
    using nci_surface::Family;
    const auto hb = nci_surface::family_bit(Family::HydrogenBond);
    const auto ring = nci_surface::family_bit(Family::Ring);
    ASSERT_TRUE(nci_surface::single_family(hb));
    ASSERT_TRUE(nci_surface::single_family(ring));
    EXPECT_FALSE(nci_surface::single_family(hb | ring));
    EXPECT_FALSE(nci_surface::single_family(0));
    auto defaults = nci_surface::family_defaults(hb);
    EXPECT_TRUE(defaults.established);
    EXPECT_EQ(defaults.s_iso, 1.0);
    EXPECT_EQ(defaults.rho_max, 0.07);
    EXPECT_EQ(defaults.radius, 1.7);
    defaults = nci_surface::family_defaults(ring);
    EXPECT_TRUE(defaults.established);
    EXPECT_EQ(defaults.s_iso, 1.0);
    EXPECT_EQ(defaults.rho_max, 0.05);
    EXPECT_EQ(defaults.radius, 2.0);
    defaults = nci_surface::family_defaults(nci_surface::family_bit(Family::PiStacking));
    EXPECT_TRUE(defaults.established);
    EXPECT_EQ(defaults.s_iso, 0.85);
    EXPECT_EQ(defaults.rho_max, 0.07);
    EXPECT_EQ(defaults.radius, 3.7);
    defaults = nci_surface::family_defaults(nci_surface::family_bit(Family::Cage));
    EXPECT_TRUE(defaults.established);
    EXPECT_EQ(defaults.s_iso, 1.0);
    EXPECT_EQ(defaults.rho_max, 0.05);
    EXPECT_EQ(defaults.radius, 2.0);
    defaults = nci_surface::family_defaults(nci_surface::family_bit(Family::HalogenBond));
    EXPECT_FALSE(defaults.established);
    EXPECT_EQ(defaults.s_iso, 0.9);
    EXPECT_EQ(defaults.rho_max, 0.07);
    EXPECT_EQ(defaults.radius, 2.2);
    const auto uniform = nci_surface::uniform_defaults();
    EXPECT_EQ(uniform.s_iso, 0.7);
    EXPECT_EQ(uniform.rho_max, 0.07);
    EXPECT_EQ(uniform.color_range, 0.03);

    EXPECT_FALSE(nci_surface::family_defaults(hb | ring).established);
    EXPECT_FALSE(nci_surface::family_defaults(0).established);
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

UTEST(nci_surface, unassigned_family_uses_generic_parameters) {
    const auto unknown = nci_surface::family_defaults(nci_surface::family_bit(nci_surface::Family::Unassigned));
    EXPECT_FALSE(unknown.established);
    EXPECT_EQ(unknown.s_iso, 0.5);
    EXPECT_EQ(unknown.rho_max, 0.05);
    EXPECT_EQ(unknown.radius, 2.7);
}

UTEST(nci_surface, family_previews_extract_independent_geometries) {
    Atoms pair;
    pair.centre = {{-3, 0, 0}, {3, 0, 0}};
    const auto g = box(6.0, 5.0, 0.2);
    const auto rho = pair.sample(g);
    const double anchor[] = {0, 0, 0};
    const auto params = [&](nci_surface::Family family) {
        const auto defaults = nci_surface::family_defaults(nci_surface::family_bit(family));
        nci_surface::Params p;
        p.s_iso = defaults.s_iso;
        p.rho_max = defaults.rho_max;
        p.near_radius = defaults.radius;
        p.near_xyz = anchor;
        p.connected = true;
        return p;
    };
    const auto hb = params(nci_surface::Family::HydrogenBond);
    const auto stacking = params(nci_surface::Family::PiStacking);
    nci_surface::Mesh a, b, repeated;
    nci_surface::Scatter sa, sb;
    ASSERT_TRUE(nci_surface::extract(&a, rho.data(), g, hb, nullptr, nullptr, &sa));
    ASSERT_TRUE(nci_surface::extract(&b, rho.data(), g, stacking, nullptr, nullptr, &sb));
    ASSERT_GT(a.index.size(), (size_t)0);
    ASSERT_GT(b.index.size(), (size_t)0);
    EXPECT_NE(a.position.size(), b.position.size());
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
    EXPECT_EQ(hb.s_iso, 1.0);
    EXPECT_EQ(stacking.s_iso, 0.85);
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
    // A dense cloud and a sparse one: each at its own 10%-of-peak level, so the sparse one is not lost
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

UTEST(nci_surface, expected_regions_match_the_spike_sign_of_each_type) {
    using nci_surface::Family;
    const auto hb = nci_surface::expected_region(Family::HydrogenBond);
    const auto ring = nci_surface::expected_region(Family::Ring);
    const auto pi = nci_surface::expected_region(Family::PiStacking);
    ASSERT_TRUE(hb.valid && ring.valid && pi.valid);
    EXPECT_LT(hb.rho_max, 0.0f);
    EXPECT_GT(ring.rho_min, 0.0f);
    EXPECT_TRUE(pi.rho_min < 0.0f && pi.rho_max > 0.0f);
    EXPECT_FALSE(nci_surface::expected_region(Family::Unassigned).valid);
    EXPECT_FALSE(nci_surface::expected_region(Family::Bond).valid);

    nci_surface::Scatter scatter;
    scatter.signed_rho = {-0.03f, -0.002f, 0.02f, -0.03f};
    scatter.s = {0.3f, 0.4f, 0.5f, 1.8f};
    const auto counts = nci_surface::expected_region_counts(scatter, Family::HydrogenBond);
    EXPECT_EQ(counts.inside, 1u);
    EXPECT_EQ(counts.outside, 2u);   // the sample above s_max is not counted
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
