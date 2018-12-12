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

Timestamp get_time() {
    timespec t;
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &t);
    return *(Timestamp*)(&t);
}

float compute_delta_ms(Timestamp t0, Timestamp t1) {
    timespec ts0 = *(timespec*)(&t0);
    timespec ts1 = *(timespec*)(&t1);
    return (ts1.tv_sec - ts0.tv_sec) * 1000 + (ts1.tv_nsec - ts0.tv_nsec) / 1000;
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
            strncpy(entry.name.beg(), ent->d_name, 512);

            res.push_back(entry);
        }

        /* Release file names */
        for (int i = 0; i < n; i++) {
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
    return { getcwd(data.file_system.cwd.beg(), 512) };
}

void sleep(int32 milliseconds) {
	usleep(milliseconds * 1000);
}