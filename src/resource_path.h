#pragma once

#include <core/md_path.h>
#include <core/md_str.h>
#include <core/md_writer.h>

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
static inline md_str_t resource_path(md_str_writer_t* sb, md_str_t relative) {
    MD_ASSERT(sb);

    char exe[2048];
    const size_t len = md_path_write_exe(exe, sizeof(exe));
    // A truncated path would silently name a different file, so a full buffer is a failure and not
    // a result. md_path_write_exe writes up to buf_cap, so equality is already suspect.
    if (len == 0 || len >= sizeof(exe)) {
        return {};
    }

    const md_str_t dir = md_path_folder({exe, len});
    if (md_str_empty(dir)) {
        return {};
    }

    md_str_writer_reset(sb);
    md_str_writer_push_str(sb, dir);
    md_str_writer_push_str(sb, relative);
    return md_str_writer_to_str(sb);
}

}  // namespace viamd
