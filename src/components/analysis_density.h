#pragma once

#include <core/md_vec_math.h>
#include <cmath>
#include <vector>

namespace analysis_density {

enum Display { Scatter, Isolines, FilledLevels, Colormap, DisplayCount };
static const char* const display_names[] = {"Scatter", "Isolines", "Filled levels", "Colormap"};

// A probability density on the shape simplex. Each source kernel is normalized over
// the domain, not wrapped or truncated into a neighbouring corner. The color scale
// is relative to uniform density, never to the number of frames or the grid maximum.
struct Triangle {
    static constexpr int resolution = 64;
    static constexpr float height = 0.86602540378f;
    std::vector<float> values;

    static int index(int i, int j) {
        return j * (resolution + 1) - j * (j - 1) / 2 + i;
    }

    static vec2_t position(int i, int j) {
        return {(i + 0.5f * j) / resolution, height * j / resolution};
    }

    void clear() { values.clear(); }

    void update(const vec2_t* points, size_t frames, size_t groups, size_t beg, size_t end, float sigma) {
        const int n = (resolution + 1) * (resolution + 2) / 2;
        std::vector<float> histogram(n, 0.0f);
        values.assign(n, 0.0f);
        size_t count = 0;
        for (size_t g = 0; g < groups; ++g) {
            for (size_t f = beg; f < end; ++f) {
                const vec2_t p = points[g * frames + f];
                if (!std::isfinite(p.x) || !std::isfinite(p.y)) continue;
                const int j = (int)std::round(p.y / height * resolution);
                const int i = (int)std::round(p.x * resolution - 0.5f * j);
                if (j < 0 || j > resolution || i < 0 || i + j > resolution) continue;
                histogram[index(i, j)] += 1.0f;
                ++count;
            }
        }
        if (!count) return;

        sigma = sigma < 0.25f ? 0.25f : sigma;
        const int radius = (int)std::ceil(3.0f * sigma / height);
        const float inv_sigma = 0.5f / (sigma * sigma);
        std::vector<float> kernel((2 * radius + 1) * (2 * radius + 1));
        for (int dj = -radius; dj <= radius; ++dj) {
            for (int di = -radius; di <= radius; ++di) {
                const float x = di + 0.5f * dj;
                const float y = height * dj;
                kernel[(dj + radius) * (2 * radius + 1) + di + radius] = std::exp(-(x*x + y*y) * inv_sigma);
            }
        }
        for (int j = 0; j <= resolution; ++j) {
            for (int i = 0; i + j <= resolution; ++i) {
                const float source = histogram[index(i, j)];
                if (source == 0.0f) continue;
                float norm = 0.0f;
                for (int dj = -radius; dj <= radius; ++dj) {
                    for (int di = -radius; di <= radius; ++di) {
                        const int x = i + di, y = j + dj;
                        if (x < 0 || y < 0 || x + y > resolution) continue;
                        norm += kernel[(dj + radius) * (2 * radius + 1) + di + radius];
                    }
                }
                for (int dj = -radius; dj <= radius; ++dj) {
                    for (int di = -radius; di <= radius; ++di) {
                        const int x = i + di, y = j + dj;
                        if (x < 0 || y < 0 || x + y > resolution) continue;
                        values[index(x, y)] += source / norm *
                            kernel[(dj + radius) * (2 * radius + 1) + di + radius];
                    }
                }
            }
        }
        const float scale = (float)n / (float)count;
        for (float& v : values) v *= scale;
    }

    static float color_coordinate(float density) { return density / (1.0f + density); }

    std::vector<float> texture_data() const {
        const int dim = resolution + 1;
        std::vector<float> data(dim * dim, 0.0f);
        if (values.empty()) return data;
        for (int y = 0; y < dim; ++y) {
            for (int x = 0; x < dim; ++x) {
                // Extend the nearest domain sample into the immediately adjacent
                // texels. The fragment mask clips the actual triangle, while this
                // extension prevents bilinear sampling from fading its boundary.
                const float coord = x - 0.5f * y;
                const int max_i = resolution - y;
                const int i0 = (int)std::floor(coord);
                const int a = i0 < 0 ? 0 : (i0 > max_i ? max_i : i0);
                const int b = i0 + 1 < 0 ? 0 : (i0 + 1 > max_i ? max_i : i0 + 1);
                const float t = coord - std::floor(coord);
                const float v = values[index(a,y)] * (1-t) + values[index(b,y)] * t;
                data[y * dim + x] = color_coordinate(v);
            }
        }
        return data;
    }
};

}
