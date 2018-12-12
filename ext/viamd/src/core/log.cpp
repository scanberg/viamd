#include <core/log.h>
#include <stdarg.h>
#include <stdio.h>

namespace logging {

struct CBEntry {
    LoggingFunc func = nullptr;
    void* usr_data = nullptr;
};

static DynamicArray<CBEntry> entries;
static constexpr int32 BUF_SIZE = KILOBYTES(64);
static char buf[BUF_SIZE];

void initialize() {}

void shutdown() {}

void register_backend(LoggingFunc func, void* usr_data) { entries.push_back({func, usr_data}); }

void record(Severity severity, const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    int count = vsnprintf(buf, BUF_SIZE, format, ap);
    va_end(ap);
    CString str(buf, count);

    for (const auto& e : entries) {
        e.func(str, severity, e.usr_data);
    }
}

}  // namespace logging