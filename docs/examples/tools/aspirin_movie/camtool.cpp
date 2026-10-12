// Throwaway: evaluates a movie's camera path with VIAMD's own code. Input on stdin, output on stdout.
#include <gfx/camera_utils.h>
#include <core/md_vec_math.h>
#include <stdio.h>
#include <vector>
#include <math.h>

int main() {
    int n = 0;
    if (scanf("%d", &n) != 1) return 1;
    std::vector<CameraKeyframe> keys(n);
    std::vector<int> track(n);
    for (int i = 0; i < n; ++i) {
        double t, f; float fov, dist, px, py, pz, qx, qy, qz, qw, fcx, fcy, fcz, roll;
        int use_frame, spin, axis, cst, ease, follow, tr;
        scanf("%lf %f %f %f %f %f %f %f %f %f %d %lf %d %d %d %d %d %f %f %f %d %f", &t, &fov, &dist, &px, &py, &pz, &qx, &qy, &qz, &qw, &use_frame, &f, &spin, &axis, &cst, &ease, &follow, &fcx, &fcy, &fcz, &tr, &roll);
        CameraKeyframe& k = keys[i];
        k.time = t; k.fov_y = fov; k.transform.distance = dist;
        k.transform.position = {px, py, pz};
        k.transform.orientation = quat_normalize(quat_t{qx, qy, qz, qw});
        k.use_frame = use_frame != 0; k.frame = f; k.spin_turns = spin; k.spin_axis = (SpinAxis)axis; k.spin_constant_speed = cst != 0;
        k.ease = (KeyEase)ease; k.follow = follow != 0; k.follow_center = {fcx, fcy, fcz}; k.follow_atom = tr >= 0 ? tr : -1; k.roll = roll;
        track[i] = tr;
    }
    double tb, fb, te, fe; int upright_on; float ux, uy, uz;
    scanf("%lf %lf %lf %lf %d %f %f %f", &tb, &fb, &te, &fe, &upright_on, &ux, &uy, &uz);
    int ntracks, nframes;
    scanf("%d %d", &ntracks, &nframes);
    std::vector<vec3_t> pos((size_t)ntracks * nframes);
    for (size_t i = 0; i < pos.size(); ++i) scanf("%f %f %f", &pos[i].x, &pos[i].y, &pos[i].z);
    double t0, t1, dt;
    scanf("%lf %lf %lf", &t0, &t1, &dt);
    const vec3_t up = {ux, uy, uz};
    for (double t = t0; t <= t1 + 1e-9; t += dt) {
        const double frame = camera_keyframes_frame_with_anchors(keys.data(), keys.size(), tb, fb, te, fe, t);
        const double fc = fmin(fmax(frame, 0.0), nframes - 1.0);
        const int f0 = (int)floor(fc), f1 = (int)fmin(f0 + 1, nframes - 1);
        const float w = (float)(fc - f0);
        std::vector<vec3_t> now(n);
        for (int i = 0; i < n; ++i) {
            now[i] = keys[i].follow_center;
            if (track[i] >= 0) {
                const vec3_t a = pos[(size_t)track[i] * nframes + f0], b = pos[(size_t)track[i] * nframes + f1];
                now[i] = a + (b - a) * w;
            }
        }
        ViewTransform vt; float fov;
        camera_keyframes_evaluate(&vt, &fov, keys.data(), keys.size(), t, false, nullptr, now.data(), upright_on ? &up : nullptr);
        const vec3_t la = camera_get_look_at(vt);
        printf("%.4f %.3f  %.3f %.3f %.3f  %.3f %.3f %.3f  %.5f %.3f\n", t, frame, la.x, la.y, la.z, vt.position.x, vt.position.y, vt.position.z, fov, vt.distance);
    }
    return 0;
}
