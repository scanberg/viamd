#include "isosurface.h"

#include <algorithm>

void IsoSurface::add(float v, const vec4& color) { 
    if (values.size() >= maxCount) return;
    values.push_back({v, color}); 
}

void IsoSurface::clear() { values.clear(); }

void IsoSurface::sort() {
    std::sort(values.begin(), values.end(), [](auto& a, auto& b) { return a.first < b.first; });
}

std::pair<std::vector<float>, std::vector<vec4>> IsoSurface::getData() const {
    auto isosurfaces = values;
    // isovalues need to be sorted in ascending order for the shader
    std::sort(isosurfaces.begin(), isosurfaces.end(), [](auto& a, auto& b) { return a.first < b.first; });
    std::pair<std::vector<float>, std::vector<vec4>> result;
    result.first.reserve(isosurfaces.size());
    result.second.reserve(isosurfaces.size());
    for (auto& v : isosurfaces) {
        result.first.push_back(v.first);
        result.second.push_back(v.second);
    }

    return result;
}
