#include "isosurface.h"

#include <algorithm>

void IsoSurface::add(float v, const vec4& color) { values.push_back({v, color}); }

void IsoSurface::clear() { values.clear(); }

void IsoSurface::sort() {
    std::sort(values.begin(), values.end(), [](auto& a, auto& b) { return a.first < b.first; });
}

std::pair<std::vector<float>, std::vector<vec4>> IsoSurface::getData() const {
    std::sort(values.begin(), values.end(), [](auto& a, auto& b) { return a.first < b.first; });
    std::pair<std::vector<float>, std::vector<vec4>> result;
    result.first.reserve(values.size());
    result.second.reserve(values.size());
    for (auto& v : values) {
        result.first.push_back(v.first);
        result.second.push_back(v.second);
    }

    return result;
}
