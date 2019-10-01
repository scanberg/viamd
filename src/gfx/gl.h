#pragma once

#if defined(_WIN32) && !defined(APIENTRY)
#define APIENTRY __stdcall
#endif
#include <GL/gl3w.h>
