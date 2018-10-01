#include <mol/element.h>
#include <string.h>
#include <ctype.h>

namespace element {

Element get_from_string(CString cstr) {
    if (cstr.count == 0) return Element::Unknown;

    const char* beg = cstr.beg();
    const char* end = cstr.end();
    while (!isalpha(*beg) && beg != end) beg++;
    if (beg == end) return Element::Unknown;
    const char* tmp = beg + 1;
    while (tmp != end && isalpha(*tmp)) tmp++;
    end = tmp;
    if (end - beg > 3) return Element::Unknown;

    for (int i = 0; i < 119; i++) {
        if (compare(CString(beg, end - beg), detail::symbols[i], true)) {
            return (Element)i;
        }
    }
    return Element::Unknown;
}

}  // namespace element
