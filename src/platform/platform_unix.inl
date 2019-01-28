#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

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
                        StringBuffer<32> buf;
                        snprintf(buf.cstr(), 32, ".%.*s", (int32)f.length(), f.cstr());
                        if (compare_ignore_case(s + len, buf)) return 1;
                    }
                }
            }
            break;
        default:
        break; 
    }

    return 0;
}

Timestamp get_time() {
    timespec t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
    return t.tv_sec * 1000000000 + t.tv_nsec;
}

float compute_delta_ms(Timestamp ts0, Timestamp ts1) {
    return (ts1 - ts0) * 1.0e-6f;
}

DynamicArray<DirectoryEntry> list_directory(CString dir_path) {
    struct dirent **files;
    int n = scandir (dir_path, &files, NULL, alphasort);
	DynamicArray<DirectoryEntry> res{};

    if (n >= 0) {
        /* Loop through file names */
        for (int i = 0; i < n; i++) {
            struct dirent *ent;

            /* Get pointer to file entry */
            ent = files[i];

            DirectoryEntry entry;

            /* Output file name */
            switch (ent->d_type) {
            case DT_REG:
                entry.type = DirectoryEntry::File;
                //printf ("%s\n", ent->d_name);
                break;

            case DT_DIR:
                entry.type = DirectoryEntry::Dir;
                //printf ("%s/\n", ent->d_name);
                break;

            case DT_LNK:
                entry.type = DirectoryEntry::Link;
                //printf ("%s@\n", ent->d_name);
                break;

            default:
                entry.type = DirectoryEntry::Unknown;
                //printf ("%s*\n", ent->d_name);
            }
            strncpy(entry.name.cstr(), ent->d_name, 512);

            res.push_back(entry);
        }

        /* Release file names */
        for (int i = 0; i < n; i++) {
            free (files[i]);
        }
        free (files);

    } else {
        printf ("Cannot open directory %.*s\n", (int32)dir_path.length(), dir_path.cstr());
    }

    return res;
}

CString get_cwd() {
    return { getcwd(data.file_system.cwd.cstr(), 512) };
}

void sleep(int32 milliseconds) {
	usleep(milliseconds * 1000);
}