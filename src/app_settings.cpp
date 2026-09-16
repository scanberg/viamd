#include <app_settings.h>
#include <serialization_utils.h>

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_str.h>

#include <imgui.h>
#include <imgui_internal.h>

#include <string.h>

namespace app_settings {

// The set of application settings is known at startup and is small, so the
// tables are fixed capacity and this module allocates nothing.
enum {
    MAX_ENTRIES = 64,
    MAX_HOOKS   = 16,
    MAX_KEY_LEN = 63,
};

enum EntryType {
    EntryType_Bool,
    EntryType_Int,
    EntryType_Float,
    EntryType_Str,
};

struct Entry {
    char      key[MAX_KEY_LEN + 1];
    EntryType type;
    void*     ptr;
    size_t    cap;  // EntryType_Str only
};

struct Hook {
    void (*fn)(void* user_data);
    void* user_data;
};

static Entry  s_entries[MAX_ENTRIES];
static size_t s_entry_count = 0;

static Hook   s_hooks[MAX_HOOKS];
static size_t s_hook_count = 0;

static bool s_initialized = false;
static bool s_hooks_ran   = false;

static Entry* find_entry(str_t key) {
    for (size_t i = 0; i < s_entry_count; ++i) {
        if (str_eq_cstr(key, s_entries[i].key)) {
            return &s_entries[i];
        }
    }
    return nullptr;
}

static void add_entry(str_t key, EntryType type, void* ptr, size_t cap) {
    ASSERT(ptr);
    ASSERT(!s_initialized && "app_settings::bind must be called before initialize");

    if (key.len == 0 || key.len > MAX_KEY_LEN) {
        MD_LOG_ERROR("app_settings: invalid key '" STR_FMT "', not bound", STR_ARG(key));
        return;
    }
    if (find_entry(key)) {
        MD_LOG_ERROR("app_settings: key '" STR_FMT "' is already bound, ignoring rebind", STR_ARG(key));
        return;
    }
    if (s_entry_count == MAX_ENTRIES) {
        MD_LOG_ERROR("app_settings: out of entry slots, '" STR_FMT "' will not be persisted", STR_ARG(key));
        return;
    }

    Entry& entry = s_entries[s_entry_count++];
    str_copy_to_char_buf(entry.key, sizeof(entry.key), key);
    entry.type = type;
    entry.ptr  = ptr;
    entry.cap  = cap;
}

void bind(str_t key, bool*  value) { add_entry(key, EntryType_Bool,  value, 0); }
void bind(str_t key, int*   value) { add_entry(key, EntryType_Int,   value, 0); }
void bind(str_t key, float* value) { add_entry(key, EntryType_Float, value, 0); }
void bind(str_t key, char*  buf, size_t cap) {
    ASSERT(cap > 0);
    add_entry(key, EntryType_Str, buf, cap);
}

void on_apply(void (*fn)(void* user_data), void* user_data) {
    ASSERT(fn);
    ASSERT(!s_initialized && "app_settings::on_apply must be called before initialize");
    if (s_hook_count == MAX_HOOKS) {
        MD_LOG_ERROR("app_settings: out of apply hook slots, hook will not run");
        return;
    }
    s_hooks[s_hook_count++] = {fn, user_data};
}

void mark_dirty() {
    if (ImGui::GetCurrentContext()) {
        ImGui::MarkIniSettingsDirty();
    }
}

// A single entry, so the only name we answer to is "Settings".
static void* read_open(ImGuiContext*, ImGuiSettingsHandler*, const char* name) {
    return strcmp(name, "Settings") == 0 ? (void*)1 : nullptr;
}

static void read_line(ImGuiContext*, ImGuiSettingsHandler*, void*, const char* line) {
    const str_t str = str_trim(str_from_cstr(line));

    size_t loc;
    if (!str_find_char(&loc, str, '=')) {
        return;
    }

    const str_t key = str_trim(str_substr(str, 0, loc));
    const str_t val = str_trim(str_substr(str, loc + 1));

    Entry* entry = find_entry(key);
    if (!entry) {
        // A setting written by another version, or one whose binding was removed.
        // Not an error: it is simply dropped on the next write.
        MD_LOG_DEBUG("app_settings: ignoring unknown setting '" STR_FMT "'", STR_ARG(key));
        return;
    }

    bool ok = false;
    switch (entry->type) {
    case EntryType_Bool:  ok = viamd::extract_bool(*(bool*)entry->ptr,  val); break;
    case EntryType_Int:   ok = viamd::extract_int (*(int*)entry->ptr,   val); break;
    case EntryType_Float: ok = viamd::extract_flt (*(float*)entry->ptr, val); break;
    case EntryType_Str:
        // Taken verbatim: an ini line cannot hold a newline, so there is nothing to unescape.
        str_copy_to_char_buf((char*)entry->ptr, entry->cap, val);
        ok = true;
        break;
    default:
        ASSERT(false);
    }

    if (!ok) {
        MD_LOG_ERROR("app_settings: could not parse value '" STR_FMT "' for setting '" STR_FMT "', keeping default", STR_ARG(val), STR_ARG(key));
    }
}

static void run_hooks() {
    s_hooks_ran = true;
    for (size_t i = 0; i < s_hook_count; ++i) {
        s_hooks[i].fn(s_hooks[i].user_data);
    }
}

static void apply_all(ImGuiContext*, ImGuiSettingsHandler*) {
    run_hooks();
}

static void write_all(ImGuiContext*, ImGuiSettingsHandler* handler, ImGuiTextBuffer* out_buf) {
    if (s_entry_count == 0) {
        return;
    }

    out_buf->reserve(out_buf->size() + (int)s_entry_count * 64);
    out_buf->appendf("[%s][Settings]\n", handler->TypeName);

    for (size_t i = 0; i < s_entry_count; ++i) {
        const Entry& entry = s_entries[i];
        switch (entry.type) {
        case EntryType_Bool:  out_buf->appendf("%s=%d\n", entry.key, *(const bool*)entry.ptr ? 1 : 0); break;
        case EntryType_Int:   out_buf->appendf("%s=%d\n", entry.key, *(const int*)entry.ptr);          break;
        case EntryType_Float: out_buf->appendf("%s=%g\n", entry.key, *(const float*)entry.ptr);        break;
        case EntryType_Str:   out_buf->appendf("%s=%s\n", entry.key, (const char*)entry.ptr);          break;
        default:
            ASSERT(false);
        }
    }

    out_buf->append("\n");
}

void initialize() {
    ASSERT(ImGui::GetCurrentContext() && "app_settings::initialize requires an ImGui context");
    ASSERT(!s_initialized && "app_settings::initialize must be called once");
    // Past the first frame the .ini has already been read and the bindings would never be filled in.
    ASSERT(!ImGui::GetCurrentContext()->SettingsLoaded && "app_settings::initialize must run before the first frame");

    ImGuiSettingsHandler handler;
    handler.TypeName   = "VIAMD";
    handler.TypeHash   = ImHashStr("VIAMD");
    handler.ReadOpenFn = read_open;
    handler.ReadLineFn = read_line;
    handler.ApplyAllFn = apply_all;
    handler.WriteAllFn = write_all;
    ImGui::AddSettingsHandler(&handler);

    s_initialized = true;

    // Read the file here instead of leaving it to the first NewFrame, so every bound value is
    // settled the moment this returns. ImGui flags the settings as loaded and skips its own read.
    if (ImGui::GetIO().IniFilename) {
        ImGui::LoadIniSettingsFromDisk(ImGui::GetIO().IniFilename);
    }

    // Nothing to read, so the load ran no hooks -- but they still have to see the defaults.
    if (!s_hooks_ran) {
        run_hooks();
    }
}

}  // namespace app_settings
