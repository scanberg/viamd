#include <mol/element.h>
#include <string.h>
#include <ctype.h>

namespace element {

Element get_from_string(const char* cstr, int length) {
    if (length == -1) length = (int)strlen(cstr);

    if (length == 0) return Element::Unknown;
    // length = length < 3 ? length : 3;

    const char* beg = cstr;
    const char* end = cstr + length;
    while (!isupper(*beg) && beg != end) beg++;
    if (beg == end) return Element::Unknown;
    const char* tmp = beg + 1;
    while (tmp != end && islower(*tmp)) tmp++;
    end = tmp;
    if (end - beg > 3) return Element::Unknown;

    for (size_t i = 0; i < detail::symbols.size(); i++) {
        if (strncmp(beg, detail::symbols[i], strlen(detail::symbols[i])) == 0) {
            return static_cast<Element>(i);
        }
    }
    return Element::Unknown;
}

}  // namespace element