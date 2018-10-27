#include <mol/element.h>
#include <string.h>
#include <ctype.h>

namespace element {

/*
// arose NGL (https://github.com/arose/ngl)
constexpr const char* elm1_str[] = {"H", "C", "O", "N", "S", "P"};
constexpr Element elm1_val[] = {Element::H, Element::C, Element::O, Element::N, Element::S, Element::P};
constexpr int32 elm1_size = sizeof(elm1_val) / sizeof(elm1_val[0]);

constexpr const char* elm2_str[] = {"NA", "CL", "FE"};
constexpr Element elm2_val[] = {Element::Na, Element::Cl, Element::Fe};
constexpr int32 elm2_size = sizeof(elm2_val) / sizeof(elm2_val[0]);
*/

Element get_from_string(CString cstr) {
    if (cstr.count == 0) return Element::Unknown;

    // Prune
    const char* c = cstr.beg();
    while (c != cstr.end() && !isalpha(*c)) c++;
    if (c == cstr.end()) return Element::Unknown;
    cstr.data = c;
    while (c != cstr.end() && isalpha(*c)) c++;
    cstr.count = c - cstr.beg();

    // Two or more (Try to match first two)
    if (cstr.length() > 1) {
        for (int32 i = 0; i < 119; i++) {
            auto elem = symbol((Element)i);
            if (strlen(elem) == 2 && cstr[0] == elem[0] && cstr[1] == elem[1]) return (Element)i;
        }
    }

    // Try to match against first character
    for (int32 i = 0; i < 119; i++) {
        auto elem = symbol((Element)i);
        if (strlen(elem) == 1 && cstr[0] == elem[0]) return (Element)i;
    }

    // Try to match against two characters again with the latter in lower case
    if (cstr.length() > 1) {
        for (int32 i = 0; i < 119; i++) {
            auto elem = symbol((Element)i);
            if (strlen(elem) == 2 && cstr[0] == elem[0] && tolower(cstr[1]) == elem[1]) return (Element)i;
        }
    }

    return Element::Unknown;
}

}  // namespace element
