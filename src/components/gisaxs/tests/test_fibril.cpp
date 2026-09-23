#include "utest.h"

#include "../fibril_core.h"

#include <math.h>
#include <stdio.h>
#include <complex>
#include <vector>

// The fibril sweep turns slices of coarse grained beads into a continuous density along the fibril. The properties
// that matter for scattering: electrons are conserved, the density is continuous along the axis (no residual of the
// bead spacing), and the cross-section follows the orientation (twist) of the slices.

namespace {

constexpr double PI_D = 3.14159265358979323846;

struct Chain {
    std::vector<float> x, y, z, e;
    std::vector<uint32_t> off;
};

// A hexagonal 7 bead slice (center + 6), twisted by 'twist' per slice, along x. Wrapped into a periodic box of 'box'.
Chain make_chain(int n, double h, double r, double twist, double box) {
    Chain c;
    for (int j = 0; j < n; ++j) {
        const double cx = j * h;
        const double th = j * twist;
        auto put = [&](double dy, double dz, double e) {
            c.x.push_back((float)(box > 0 ? fmod(cx, box) : cx));
            c.y.push_back((float)(50.0 + dy));
            c.z.push_back((float)(50.0 + dz));
            c.e.push_back((float)e);
        };
        put(0, 0, 300);
        for (int b = 0; b < 6; ++b) put(r * cos(th + b * PI_D / 3), r * sin(th + b * PI_D / 3), 250);
    }
    c.off = {0, (uint32_t)n};
    return c;
}

FibrilInput make_input(const Chain& c, double box) {
    FibrilInput in;
    in.num_slices = c.off.back();
    in.stride = 7;
    in.x = c.x.data(); in.y = c.y.data(); in.z = c.z.data(); in.e = c.e.data();
    in.chain_offset = c.off.data();
    in.num_chains = c.off.size() - 1;
    in.box[0] = box;
    in.ref_bead = 1;
    return in;
}

double axial_structure_factor(const std::vector<float>& x, const std::vector<float>& w, double q) {
    std::complex<double> s = 0;
    double t = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        s += (double)w[i] * std::exp(std::complex<double>(0, -q * x[i]));
        t += w[i];
    }
    return std::abs(s) / t;
}

}  // namespace

UTEST(fibril, mean_layout_recovers_slice_geometry) {
    const Chain c = make_chain(30, 21.0, 12.0, 10.0 * PI_D / 180.0, 0.0);
    const FibrilInput in = make_input(c, 0.0);
    std::vector<float> u, v, f;
    ASSERT_TRUE(fibril_mean_layout(u, v, f, in));
    ASSERT_EQ(u.size(), (size_t)7);
    EXPECT_NEAR(u[0], 0.0f, 1e-3f);
    EXPECT_NEAR(v[0], 0.0f, 1e-3f);
    for (int b = 1; b < 7; ++b) {
        EXPECT_NEAR(sqrt(u[b] * u[b] + v[b] * v[b]), 12.0, 1e-2);
    }
    EXPECT_NEAR(f[0], 300.0f / 1800.0f, 1e-5f);
}

UTEST(fibril, sweep_conserves_electrons_and_is_continuous) {
    // Crosses the periodic boundary on purpose
    const double h = 21.0;
    const Chain c = make_chain(40, h, 12.0, 10.0 * PI_D / 180.0, 300.0);
    const FibrilInput in = make_input(c, 300.0);
    std::vector<float> u, v, f;
    ASSERT_TRUE(fibril_mean_layout(u, v, f, in));

    const FibrilTemplate beads = fibril_template_beads(u, v, f, 4.0f);
    const FibrilTemplate disk = fibril_template_disk(u, v, 19.1f, 8.0f);
    for (const FibrilTemplate* t : {&beads, &disk}) {
        FibrilOutput out;
        char err[256] = "";
        ASSERT_TRUE(fibril_sweep(&out, in, *t, 0.0, err, sizeof(err)));
        EXPECT_NEAR(out.stats.electrons_out / out.stats.electrons_in, 1.0, 1e-6);
        EXPECT_NEAR(out.stats.mean_spacing, h, 0.2);
        // The bead period must not survive in the swept density (it is exactly 1 for the beads themselves)
        EXPECT_LT(axial_structure_factor(out.x, out.w, 2.0 * PI_D / h), 1e-3);
    }
}

UTEST(fibril, sweep_follows_twist) {
    const double h = 21.0, tw = 25.0 * PI_D / 180.0;
    const Chain c = make_chain(20, h, 12.0, tw, 0.0);
    const FibrilInput in = make_input(c, 0.0);
    std::vector<float> u, v, f;
    ASSERT_TRUE(fibril_mean_layout(u, v, f, in));
    // One sample per segment (at the midpoints), component 1 is the anchor at angle 0
    const FibrilTemplate t = fibril_template_beads(u, v, f, (float)h);
    FibrilOutput out;
    char err[256] = "";
    ASSERT_TRUE(fibril_sweep(&out, in, t, 0.0, err, sizeof(err)));
    const size_t G = t.gauss.size();
    double ref = 0.0, max_dev = 0.0;
    bool first = true;
    for (size_t i = 0; i + G <= out.x.size(); i += G) {
        const double X = out.x[i];
        if (X < 0.0 || X > (20 - 1) * h) continue;     // interior segments only
        const double a = atan2(out.z[i + 1] - 50.0, out.y[i + 1] - 50.0);
        const double d = remainder(a - (X / h) * tw, 2.0 * PI_D);
        if (first) { ref = d; first = false; }
        max_dev = fmax(max_dev, fabs(remainder(d - ref, 2.0 * PI_D)));
    }
    EXPECT_FALSE(first);
    EXPECT_LT(max_dev, 1e-3);
}

UTEST(fibril, template_file_and_handedness) {
    char path[] = "fibril_template_test.txt";
    FILE* fp = fopen(path, "w");
    ASSERT_TRUE(fp != NULL);
    fprintf(fp, "# test\nanchor 0 0\nanchor 10 0\nanchor 0 5\ngauss 0 0 2 3\ngauss 10 0 1 3\ngauss 0 5 1 4\n");
    fclose(fp);

    FibrilTemplate t;
    char err[256] = "";
    ASSERT_TRUE(fibril_template_load(&t, path, err, sizeof(err)));
    remove(path);
    ASSERT_EQ(t.anchor_u.size(), (size_t)3);
    ASSERT_EQ(t.gauss.size(), (size_t)3);
    EXPECT_NEAR(t.gauss[0].weight, 0.5f, 1e-6f);

    // Reference layout: the anchors mirrored and rotated by 90 degrees -> the template has to be mirrored
    std::vector<float> ru, rv;
    for (size_t i = 0; i < t.anchor_u.size(); ++i) {
        const float uu = t.anchor_u[i], vv = -t.anchor_v[i];
        ru.push_back(-vv);
        rv.push_back(uu);
    }
    double rms = 1.0;
    EXPECT_TRUE(fibril_template_match_handedness(t, ru, rv, &rms));
    EXPECT_LT(rms, 1e-4);
    EXPECT_FALSE(fibril_template_match_handedness(t, ru, rv, &rms));
}
