#pragma once

#include <core/string_utils.h>

#define LOG_NOTE(format, ...) logging::record(logging::Note, format, __VA_ARGS__);
#define LOG_WARNING(format, ...) logging::record(logging::Warning, format, __VA_ARGS__);
#define LOG_ERROR(format, ...) logging::record(logging::Error, format, __VA_ARGS__);
#define LOG_FATAL(format, ...) logging::record(logging::Fatal, format, __VA_ARGS__);

namespace logging {
enum Severity { Note, Warning, Error, Fatal };

typedef void (*LoggingFunc)(CString str, Severity severity, void* usr_data);

void initialize();
void shutdown();

void register_backend(LoggingFunc func, void* usr_data = nullptr);

void record(Severity severity, const char* format, ...);
}  // namespace log