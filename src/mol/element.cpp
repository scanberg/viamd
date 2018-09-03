#include <mol/element.h>
#include <string.h>
#include <ctype.h>

namespace element {

Element get_from_string(CString cstr) {
    if (cstr.count == 0) return Element::Unknown;

    const char* beg = cstr.beg();
    const char* end = cstr.end();
    while (!isupper(*beg) && beg != end) beg++;
    if (beg == end) return Element::Unknown;
    const char* tmp = beg + 1;
    while (tmp != end && islower(*tmp)) tmp++;
    end = tmp;
    if (end - beg > 3) return Element::Unknown;

    for (size_t i = 0; i < detail::symbols.size(); i++) {
        if (strncmp(beg, detail::symbols[i], strlen(detail::symbols[i])) == 0) {
            return (Element)i;
        }
    }
    return Element::Unknown;
}

}  // namespace element
