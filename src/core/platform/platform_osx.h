#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <core/string_utils.h>

static DynamicArray<CString> curr_filters;

static int filter_func(const struct dirent *dir) {
    switch (dir->d_type) {
        case DT_REG:
        case DT_LNK:
            {
                const char *s = dir->d_name;
                int s_len = strlen(s);
                for (const auto& f : curr_filters) {
                    int len = s_len - f.count - 1;
                    if (len >= 0) {
                        StringBuffer<32> ext = f;
                        StringBuffer<32> buf;
                        snprintf(buf.beg(), 32, ".%s", ext.beg());
                        if (compare(s + len, buf, true)) return 1;
                    }
                }
            }
            break;
        default:
        break; 
    }

    return 0;
}

DynamicArray<DirEntry> list_directory(CString dir_path, CString filter) {
    struct dirent **files;
    int i;
    int n;

    DynamicArray<DirEntry> res;

    if (filter.count > 0) {
        curr_filters = ctokenize(filter, '|');
        n = scandir (dir_path, &files, filter_func, alphasort);        
    }
    else {
        n = scandir (dir_path, &files, NULL, alphasort);
    }

    if (n >= 0) {
        /* Loop through file names */
        for (i = 0; i < n; i++) {
            struct dirent *ent;

            /* Get pointer to file entry */
            ent = files[i];

            DirEntry entry;

            /* Output file name */
            switch (ent->d_type) {
            case DT_REG:
                entry.type = DirEntry::File;
                //printf ("%s\n", ent->d_name);
                break;

            case DT_DIR:
                entry.type = DirEntry::Dir;
                //printf ("%s/\n", ent->d_name);
                break;

            case DT_LNK:
                entry.type = DirEntry::Link;
                //printf ("%s@\n", ent->d_name);
                break;

            default:
                entry.type = DirEntry::Unknown;
                //printf ("%s*\n", ent->d_name);
            }
            strncpy(entry.name.beg(), ent->d_name, 512);

            res.push_back(entry);
        }

        /* Release file names */
        for (i = 0; i < n; i++) {
            free (files[i]);
        }
        free (files);

    } else {
        StringBuffer<512> buf = dir_path;
        printf ("Cannot open directory %s\n", buf.beg());
    }

    return res;
}

CString get_cwd() {
    return {getcwd(NULL,0)};
}

/*
struct DirectoryWatch {
    FSEventStreamRef fs_stream;
    Url path;
    DynamicArray<Url> files;
};

static DynamicArray<DirectoryWatch> s_watched_dirs;

bool AddDirectoryWatch(CString path_to_watch, void (*callback)(FileEvent)) {
    //Define variables and create a CFArray object containing
     CFString objects containing paths to watch.
     
    CFStringRef path = CFStringCreateWithCString(NULL, path_to_watch, CFStringBuiltInEncodings::kCFStringEncodingUTF8);
    CFArrayRef path_array = CFArrayCreate(NULL, (const void **)&path, 1, NULL);
    FSEventStreamRef stream;
    const CFAbsoluteTime latency = 1.0; // Latency in seconds
    
    // Create the stream, passing in a callback
    stream = FSEventStreamCreate(NULL,
                                 &mycallback,
                                 NULL,
                                 path_array,
                                 kFSEventStreamEventIdSinceNow,
                                 latency,
                                 kFSEventStreamCreateFlagNone);
    
    FSEventStreamScheduleWithRunLoop(stream, CFRunLoopGetCurrent(), kCFRunLoopDefaultMode);
    FSEventStreamStart(stream);
    
    return true;
}

bool RemoveWatch(CString path_to_watch) {
    CString dir = get_directory(path_to_watch);
    CString file = get_file(path_to_watch);
    
    FSEventStreamStop(stream);
    FSEventStreamInvalidate(stream);
    FSEventStreamRelease(stream);
}
*/