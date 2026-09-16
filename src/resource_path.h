#pragma once

#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_str_builder.h>

namespace viamd {

// viamd's resources live beside the executable, so that the directory as a whole can be copied or
// moved and still run. Every lookup therefore resolves against the executable's own path and never
// against the working directory, which belongs to whoever launched it and says nothing about where
// the program was installed. The build puts them there: see create_copy_resource_dir_target in the
// top level CMakeLists, and the macOS bundle case which lands them in Contents/MacOS for the same
// reason - so that this one rule needs no platform branch.
//
// Returns an empty string when the executable path cannot be determined or does not fit, which the
// caller has to treat as "the resource is not there" rather than as a path.
static inline str_t resource_path(md_strb_t* sb, str_t relative) {
    ASSERT(sb);

    char exe[2048];
    const size_t len = md_path_write_exe(exe, sizeof(exe));
    // A truncated path would silently name a different file, so a full buffer is a failure and not
    // a result. md_path_write_exe writes up to buf_cap, so equality is already suspect.
    if (len == 0 || len >= sizeof(exe)) {
        return {};
    }

    str_t dir = {};
    if (!extract_folder_path(&dir, {exe, len})) {
        return {};
    }

    md_strb_reset(sb);
    md_strb_push_str(sb, dir);
    md_strb_push_str(sb, relative);
    return md_strb_to_str(*sb);
}

}  // namespace viamd
