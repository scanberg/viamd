#pragma once

#include "string_utils.h"

#define LOG_NOTE(...) logging::record(logging::Note, __VA_ARGS__);
#define LOG_WARNING(...) logging::record(logging::Warning, __VA_ARGS__);
#define LOG_ERROR(...) logging::record(logging::Error, __VA_ARGS__);
#define LOG_FATAL(...) logging::record(logging::Fatal, __VA_ARGS__);

namespace logging {
enum Severity { Note, Warning, Error, Fatal };

typedef void (*LoggingFunc)(CString str, Severity severity, void* usr_data);

void initialize();
void shutdown();

void register_backend(LoggingFunc func, void* usr_data = nullptr);

void record(Severity severity, const char* format, ...);
}  // namespace logging
