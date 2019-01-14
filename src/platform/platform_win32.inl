#include <Windows.h>
#include <direct.h>
#pragma comment(lib, "User32.lib")

Timestamp get_time() {
    LARGE_INTEGER t;
    QueryPerformanceCounter(&t);
    return t.QuadPart;
}

float compute_delta_ms(Timestamp t0, Timestamp t1) {
    LARGE_INTEGER start, stop, elapsed, frequency;
    QueryPerformanceFrequency(&frequency);
    start.QuadPart = t0;
    stop.QuadPart = t1;
    elapsed.QuadPart = stop.QuadPart - start.QuadPart;
    elapsed.QuadPart *= 1000;
    float ms = (float)((double)elapsed.QuadPart / (double)frequency.QuadPart);
    return ms;
}

DynamicArray<DirectoryEntry> list_directory(CString dir_path) {
    DynamicArray<DirectoryEntry> res;

    WIN32_FIND_DATA ffd;
    HANDLE h_find = INVALID_HANDLE_VALUE;

    // Copy and make a zero terminated string
    StringBuffer<MAX_PATH> dir_buf = dir_path;

    char dir[MAX_PATH];
    snprintf(dir, MAX_PATH, "%s\\*", dir_buf.beg());

    h_find = FindFirstFile(dir, &ffd);
    if (h_find == INVALID_HANDLE_VALUE) {
        printf("ERROR! Could not read directory '%s'", dir);
        return {};
    }

    do {
        DirectoryEntry entry;
        if (ffd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
            entry.type = DirectoryEntry::Dir;
        } else {
            entry.type = DirectoryEntry::File;
        }
        strncpy(entry.name.cstr(), ffd.cFileName, MAX_PATH);
        res.push_back(entry);
    } while (FindNextFile(h_find, &ffd) != 0);

    FindClose(h_find);

    return res;
}

CString get_cwd() { return {_getcwd(data.file_system.cwd.cstr(), 512)}; }

void sleep(int32 milliseconds) { Sleep(milliseconds); }
