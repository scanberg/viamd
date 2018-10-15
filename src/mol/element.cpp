#include <mol/element.h>
#include <string.h>
#include <ctype.h>

namespace element {

// arose NGL (https://github.com/arose/ngl)
constexpr const char* elm1_str[] = {"H", "C", "O", "N", "S", "P"};
constexpr Element elm1_val[] = {Element::H, Element::C, Element::O, Element::N, Element::S, Element::P};
constexpr int32 elm1_size = sizeof(elm1_val) / sizeof(elm1_val[0]);

constexpr const char* elm2_str[] = {"NA", "CL", "FE"};
constexpr Element elm2_val[] = {Element::Na, Element::Cl, Element::Fe};
constexpr int32 elm2_size = sizeof(elm2_val) / sizeof(elm2_val[0]);

Element get_from_string(CString cstr, bool ignore_case) {
    if (cstr.count == 0) return Element::Unknown;

    const char* beg = cstr.beg();
    const char* end = cstr.end();
    while (!isalpha(*beg) && beg != end) beg++;
    if (beg == end) return Element::Unknown;
    const char* tmp = beg + 1;
    while (tmp != end && isalpha(*tmp)) tmp++;
    end = tmp;
    if (end - beg > 3) return Element::Unknown;

    cstr = CString(beg, end - beg);

    if (cstr.size() == 2) {
        for (int32 i = 0; i < elm2_size; i++) {
            if (cstr[0] == elm2_str[i][0] && toupper(cstr[1]) == elm2_str[i][1]) return elm2_val[i];
        }
        for (int32 i = 0; i < elm1_size; i++) {
            if (cstr[0] == elm1_str[i][0]) return elm1_val[i];
        }
    }
    for (int i = 0; i < 119; i++) {
        if (compare(cstr, detail::symbols[i], ignore_case)) {
            return (Element)i;
        }
    }
    return Element::Unknown;
}

}  // namespace element
