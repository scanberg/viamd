#include "profiling.h"
#include <core/hash.h>
#include <core/platform.h>
#include <core/gl.h>

namespace cpu_profiling {

typedef unsigned int Key;

struct Section {
    Key key = 0;
    float32 total_ms = 0;
    platform::Timestamp t0 = 0;

    int32 parent_idx = -1;
    StringBuffer<32> name{};
};

constexpr int MAX_SECTIONS = 1024;
static Section sections[1024];
static int section_count = 0;
static int current_section_idx = -1;

void initialize() {}

void shutdown() {}

int32 find_section_idx(Key key) {
    for (int32 i = 0; i < section_count; i++) {
        if (key == sections[i].key) return i;
    }
    return -1;
}

void push_section(CString sec) {
    Key key = hash::crc32(sec);
    int idx = find_section_idx(key);
    if (idx == -1) {
        idx = section_count++;
        sections[idx] = {key, 0, 0, current_section_idx, sec};
    }
    current_section_idx = idx;
    sections[current_section_idx].t0 = platform::get_time();
}

void pop_section(CString sec) {
#ifdef DEBUG
    ASSERT(current_section_idx != -1);
    Key key = hash::crc32(sec);
    int idx = find_section_idx(key);
    ASSERT(current_section_idx == idx);
#else
    (void)sec;
#endif
    auto& section = sections[current_section_idx];
    section.total_ms += platform::compute_delta_ms(section.t0, platform::get_time());
    current_section_idx = sections[current_section_idx].parent_idx;
}

void clear() {
    section_count = 0;
    current_section_idx = -1;
}
void finish() {}

}  // namespace cpu_profiling

namespace gpu_profiling {

typedef unsigned int Key;

struct Section {
    Key key = 0;
    float total_ms = 0;
    unsigned int query_ids[2] = {0, 0};

    int parent_idx = -1;
    StringBuffer<32> name{};
};

constexpr int MAX_SECTIONS = 1024;
static Section sections[MAX_SECTIONS];
static int section_count = 0;
static int current_section_idx = -1;

void initialize() {}
void shutdown() {}

int32 find_section_idx(Key key) {
    for (int32 i = 0; i < section_count; i++) {
        if (key == sections[i].key) return i;
    }
    return -1;
}

void push_section(CString sec) {
    Key key = hash::crc32(sec);
    int idx = find_section_idx(key);
    if (idx == -1) {
        idx = section_count++;
        sections[idx] = {key, 0, {0, 0}, current_section_idx, sec};
        glGenQueries(2, sections[idx].query_ids);
    }
    current_section_idx = idx;
    glQueryCounter(sections[idx].query_ids[0], GL_TIMESTAMP);
}

void pop_section(CString sec) {
#ifdef DEBUG
    ASSERT(current_section_idx != -1);
    Key key = hash::crc32(sec);
    int idx = find_section_idx(key);
    ASSERT(current_section_idx == idx);
#else
    (void)sec;
#endif
    auto& section = sections[current_section_idx];
    glQueryCounter(section.query_ids[1], GL_TIMESTAMP);
    current_section_idx = section.parent_idx;
}

void clear() { section_count = 0; }

void finish() {
    for (int i = 0; i < section_count; i++) {
        int done = 0;
        while (!done) {
            glGetQueryObjectiv(sections[i].query_ids[1], GL_QUERY_RESULT_AVAILABLE, &done);
        }
        GLuint64 t0, t1;
        glGetQueryObjectui64v(sections[i].query_ids[0], GL_QUERY_RESULT, &t0);
        glGetQueryObjectui64v(sections[i].query_ids[1], GL_QUERY_RESULT, &t1);
        sections[i].total_ms += (t1 - t0) / 1000000.f;
    }
}

int get_depth(int section_idx) {
    int count = 0;
    Section* sec = &sections[section_idx];
    while (sec->parent_idx != -1) {
        count++;
        sec = &sections[sec->parent_idx];
    }
    return count;
}

void print() {
    StringBuffer<32> buff;
    for (int i = 0; i < section_count; i++) {
        int count = 0;
        count += snprintf(buff.beg(), 32, "%*c [%s]:", get_depth(i) * 4, ' ', sections[i].name.cstr());
        printf("%s%8.4f\n", buff.cstr(), sections[i].total_ms);
    }
    printf("\n");
}

};  // namespace gpu_profiling
