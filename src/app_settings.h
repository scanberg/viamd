#pragma once

#include <stddef.h>

#include <core/md_str.h>

// Application level settings, persisted in the ImGui .ini file under a single
// [VIAMD][Settings] entry.
//
// The value keeps living where it is used. This layer only binds a key to the
// address of an existing variable, so nothing has to move into a settings blob
// in order to become persistent:
//
//     app_settings::bind(STR_LIT("keep_representations"), &state.settings.keep_representations);
//
// Bind everything first, then call initialize(): it reads the .ini, so every
// bound value is live the moment it returns. A bound address has to outlive the
// ImGui context. After editing a bound value from the UI, call mark_dirty() so
// the change is written back out.

namespace app_settings {

// Installs the .ini handler and reads the file. Call once, after the ImGui
// context exists, after every bind()/on_apply(), and before the first frame.
void initialize();

void bind(str_t key, bool*  value);
void bind(str_t key, int*   value);
void bind(str_t key, float* value);
void bind(str_t key, char*  buf, size_t cap);

// Runs from initialize(), after the .ini has been read, for settings that have to
// be pushed somewhere before they take effect. Runs whether or not a .ini existed,
// so the hook also sees the defaults on a first run.
void on_apply(void (*fn)(void* user_data), void* user_data);

// Flags the .ini as dirty. ImGui writes it out shortly after.
void mark_dirty();

}  // namespace app_settings
