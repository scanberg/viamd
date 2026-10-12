#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#else
#include <errno.h>
#include <fcntl.h>
#include <pthread.h>
#include <signal.h>
#include <spawn.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
extern char** environ;
#endif

#include "frame_sink.h"

#include <image.h>

#include <algorithm>
#include <condition_variable>
#include <deque>
#include <mutex>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <thread>

namespace frame_sink {

namespace {

struct Process {
#if defined(_WIN32)
    HANDLE in = nullptr;
    HANDLE proc = nullptr;
#else
    int   in = -1;
    pid_t pid = -1;
#endif
};

#if defined(_WIN32)

std::string quote_arg(const std::string& a) {
    if (!a.empty() && a.find_first_of(" \t\"") == std::string::npos) return a;
    std::string r = "\"";
    for (char c : a) {
        if (c == '"') r += "\\\"";
        else          r += c;
    }
    r += "\"";
    return r;
}

bool process_start(Process& p, const std::vector<std::string>& args, std::string& err) {
    SECURITY_ATTRIBUTES sa = { sizeof(sa), nullptr, TRUE };
    HANDLE rd = nullptr, wr = nullptr;
    if (!CreatePipe(&rd, &wr, &sa, 0)) {
        err = "could not create a pipe";
        return false;
    }
    SetHandleInformation(wr, HANDLE_FLAG_INHERIT, 0);

    std::string cmd;
    for (const std::string& a : args) {
        if (!cmd.empty()) cmd += ' ';
        cmd += quote_arg(a);
    }

    STARTUPINFOA si = {};
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdInput  = rd;
    si.hStdOutput = GetStdHandle(STD_OUTPUT_HANDLE);
    si.hStdError  = GetStdHandle(STD_ERROR_HANDLE);

    PROCESS_INFORMATION pi = {};
    const BOOL ok = CreateProcessA(nullptr, cmd.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW, nullptr, nullptr, &si, &pi);
    CloseHandle(rd);
    if (!ok) {
        CloseHandle(wr);
        err = "could not start '" + args[0] + "' (error " + std::to_string((unsigned long)GetLastError()) + ")";
        return false;
    }
    CloseHandle(pi.hThread);
    p.in = wr;
    p.proc = pi.hProcess;
    return true;
}

bool process_write(Process& p, const void* data, size_t bytes) {
    const char* ptr = (const char*)data;
    while (bytes > 0) {
        DWORD n = 0;
        if (!WriteFile(p.in, ptr, (DWORD)std::min<size_t>(bytes, 1u << 24), &n, nullptr) || n == 0) return false;
        ptr += n;
        bytes -= n;
    }
    return true;
}

// Closes the input and waits. Returns the exit code, or -1 if it could not be had.
int process_finish(Process& p) {
    int code = -1;
    if (p.in) { CloseHandle(p.in); p.in = nullptr; }
    if (p.proc) {
        WaitForSingleObject(p.proc, INFINITE);
        DWORD c = 0;
        if (GetExitCodeProcess(p.proc, &c)) code = (int)c;
        CloseHandle(p.proc);
        p.proc = nullptr;
    }
    return code;
}

void block_sigpipe() {}

#else

bool process_start(Process& p, const std::vector<std::string>& args, std::string& err) {
    int fds[2];
    if (pipe(fds) != 0) {
        err = "could not create a pipe";
        return false;
    }
    // Neither end may leak into other children; dup2 below gives the child its stdin without the flag
    fcntl(fds[0], F_SETFD, FD_CLOEXEC);
    fcntl(fds[1], F_SETFD, FD_CLOEXEC);

    posix_spawn_file_actions_t fa;
    posix_spawn_file_actions_init(&fa);
    posix_spawn_file_actions_adddup2(&fa, fds[0], STDIN_FILENO);

    std::vector<char*> argv;
    for (const std::string& a : args) argv.push_back(const_cast<char*>(a.c_str()));
    argv.push_back(nullptr);

    pid_t pid = -1;
    const int rc = posix_spawnp(&pid, argv[0], &fa, nullptr, argv.data(), environ);
    posix_spawn_file_actions_destroy(&fa);
    ::close(fds[0]);
    if (rc != 0) {
        ::close(fds[1]);
        err = "could not start '" + args[0] + "': " + strerror(rc);
        return false;
    }
    p.in = fds[1];
    p.pid = pid;
    return true;
}

bool process_write(Process& p, const void* data, size_t bytes) {
    const char* ptr = (const char*)data;
    while (bytes > 0) {
        const ssize_t n = write(p.in, ptr, bytes);
        if (n < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        ptr += n;
        bytes -= (size_t)n;
    }
    return true;
}

// Closes the input and waits. Returns the exit code, or -1 if it did not exit normally.
int process_finish(Process& p) {
    int code = -1;
    if (p.in >= 0) { ::close(p.in); p.in = -1; }
    if (p.pid > 0) {
        int st = 0;
        while (waitpid(p.pid, &st, 0) < 0 && errno == EINTR) {}
        if (WIFEXITED(st)) code = WEXITSTATUS(st);
        p.pid = -1;
    }
    return code;
}

// A write to a pipe whose reader has gone away should fail, not take the application down
void block_sigpipe() {
    sigset_t s;
    sigemptyset(&s);
    sigaddset(&s, SIGPIPE);
    pthread_sigmask(SIG_BLOCK, &s, nullptr);
}

#endif

struct Job {
    uint8_t* data;
    int index;
};

}  // namespace

struct Sink {
    Kind kind = Kind::PngSequence;
    std::string dir;
    std::string prefix;
    int width = 0;
    int height = 0;
    size_t frame_bytes = 0;

    std::mutex m;
    std::condition_variable cv_work;
    std::condition_variable cv_free;

    std::vector<uint8_t*> all;
    std::vector<uint8_t*> free_list;
    std::deque<Job> queue;
    std::vector<std::thread> threads;

    bool closing = false;
    int  workers_alive = 0;
    int  submitted = 0;
    int  written = 0;
    int  failed = 0;
    bool ok = true;
    char message[256] = {};

    Process proc;
};

namespace {

// First message wins, it is the one that explains the rest
void set_message_locked(Sink& s, const std::string& msg) {
    s.ok = false;
    if (s.message[0] == '\0') snprintf(s.message, sizeof(s.message), "%s", msg.c_str());
}

std::string to_std(str_t s) {
    return (s.ptr && s.len) ? std::string(s.ptr, s.len) : std::string();
}

void flip_rows(uint8_t* rgba, int width, int height) {
    const size_t stride = (size_t)width * 4;
    std::vector<uint8_t> tmp(stride);
    for (int y = 0; y < height / 2; ++y) {
        uint8_t* a = rgba + (size_t)y * stride;
        uint8_t* b = rgba + (size_t)(height - 1 - y) * stride;
        memcpy(tmp.data(), a, stride);
        memcpy(a, b, stride);
        memcpy(b, tmp.data(), stride);
    }
}

// Empty on success
std::string write_job(Sink& s, const Job& job) {
    if (s.kind == Kind::Ffmpeg) {
        // The vertical flip is ffmpeg's job (-vf vflip), there is no point in doing it here
        if (!process_write(s.proc, job.data, s.frame_bytes)) return "ffmpeg stopped accepting frames";
        return {};
    }

    char path[2048];
    const int len = snprintf(path, sizeof(path), "%s/%s_%05d.png", s.dir.c_str(), s.prefix.c_str(), job.index);
    if (len <= 0 || (size_t)len >= sizeof(path)) return "the output path is too long";

    flip_rows(job.data, s.width, s.height);
    if (!image_write_png(str_t{path, (size_t)len}, job.data, s.width, s.height)) {
        return std::string("could not write '") + path + "'";
    }
    return {};
}

void worker_main(Sink* s, bool finalizes_process) {
    block_sigpipe();

    for (;;) {
        Job job = {};
        {
            std::unique_lock<std::mutex> lk(s->m);
            s->cv_work.wait(lk, [&] { return !s->queue.empty() || s->closing; });
            if (s->queue.empty()) break;
            job = s->queue.front();
            s->queue.pop_front();
        }

        const std::string err = write_job(*s, job);

        {
            std::lock_guard<std::mutex> lk(s->m);
            if (err.empty()) {
                s->written += 1;
            } else {
                s->failed += 1;
                set_message_locked(*s, err);
            }
            s->free_list.push_back(job.data);
        }
        s->cv_free.notify_one();
    }

    if (finalizes_process) {
        const int code = process_finish(s->proc);
        if (code != 0) {
            std::lock_guard<std::mutex> lk(s->m);
            set_message_locked(*s, code < 0 ? std::string("ffmpeg did not exit normally")
                                            : "ffmpeg exited with code " + std::to_string(code) + ", see its output in the terminal");
        }
    }

    {
        std::lock_guard<std::mutex> lk(s->m);
        s->workers_alive -= 1;
    }
    s->cv_free.notify_all();
}

void set_error(char* err, size_t cap, const std::string& msg) {
    if (err && cap) snprintf(err, cap, "%s", msg.c_str());
}

}  // namespace

const char* file_extension(Codec codec) {
    return codec == Codec::Vp9 ? "webm" : "mp4";
}

void ffmpeg_arguments(const Desc& desc, std::vector<std::string>& out) {
    out.clear();
    const std::string exe = desc.ffmpeg.len ? to_std(desc.ffmpeg) : std::string("ffmpeg");
    const std::string dir = to_std(desc.dir);
    const std::string prefix = to_std(desc.prefix);

    char size[64], rate[64], crf[16];
    snprintf(size, sizeof(size), "%dx%d", desc.width, desc.height);
    snprintf(rate, sizeof(rate), "%g", (double)desc.fps);
    snprintf(crf, sizeof(crf), "%d", desc.crf);

    // The picture is flipped because OpenGL reads from the bottom row. yuv420p needs even dimensions, hence the
    // scale, which also does the colour conversion so that it can be tagged bt709.
    out = {
        exe, "-y", "-hide_banner", "-loglevel", "error",
        "-f", "rawvideo", "-pixel_format", "rgba", "-video_size", size, "-framerate", rate, "-i", "-",
        "-vf", "vflip,scale=trunc(iw/2)*2:trunc(ih/2)*2:out_color_matrix=bt709:out_range=tv,format=yuv420p",
    };
    switch (desc.codec) {
    case Codec::H264:
        out.insert(out.end(), {"-c:v", "libx264", "-crf", crf});
        break;
    case Codec::H265:
        out.insert(out.end(), {"-c:v", "libx265", "-crf", crf, "-tag:v", "hvc1", "-x265-params", "log-level=error"});
        break;
    case Codec::Vp9:
        // -b:v 0 makes the crf the only quality control
        out.insert(out.end(), {"-c:v", "libvpx-vp9", "-crf", crf, "-b:v", "0", "-row-mt", "1"});
        break;
    }
    out.insert(out.end(), {"-pix_fmt", "yuv420p", "-colorspace", "bt709", "-color_primaries", "bt709", "-color_trc", "bt709"});
    if (desc.codec != Codec::Vp9) out.insert(out.end(), {"-movflags", "+faststart"});
    out.push_back(dir + "/" + prefix + "." + file_extension(desc.codec));
}

Sink* create(const Desc& desc, char* err, size_t err_cap) {
    if (desc.width <= 0 || desc.height <= 0) {
        set_error(err, err_cap, "the frame size is empty");
        return nullptr;
    }
    if (desc.dir.len == 0) {
        set_error(err, err_cap, "no output folder");
        return nullptr;
    }

    Sink* s = new Sink();
    s->kind = desc.kind;
    s->dir = to_std(desc.dir);
    s->prefix = to_std(desc.prefix);
    s->width = desc.width;
    s->height = desc.height;
    s->frame_bytes = (size_t)desc.width * (size_t)desc.height * 4;

    size_t num_buffers = (size_t)desc.max_in_flight;
    if (num_buffers == 0) num_buffers = std::clamp<size_t>(((size_t)256 << 20) / s->frame_bytes, 2, 6);

    for (size_t i = 0; i < num_buffers; ++i) {
        uint8_t* buf = (uint8_t*)malloc(s->frame_bytes);
        if (!buf) {
            for (uint8_t* b : s->all) free(b);
            delete s;
            set_error(err, err_cap, "out of memory for the frame buffers");
            return nullptr;
        }
        s->all.push_back(buf);
        s->free_list.push_back(buf);
    }

    int num_workers = 1;
    if (desc.kind == Kind::Ffmpeg) {
        std::vector<std::string> args;
        ffmpeg_arguments(desc, args);
        std::string msg;
        if (!process_start(s->proc, args, msg)) {
            for (uint8_t* b : s->all) free(b);
            delete s;
            set_error(err, err_cap, msg);
            return nullptr;
        }
    } else {
        const unsigned hw = std::thread::hardware_concurrency();
        num_workers = (int)std::clamp<unsigned>(hw > 1 ? hw - 1 : 1, 1, 3);
    }

    s->workers_alive = num_workers;
    for (int i = 0; i < num_workers; ++i) {
        s->threads.emplace_back(worker_main, s, desc.kind == Kind::Ffmpeg);
    }
    return s;
}

uint8_t* acquire(Sink* s, bool wait) {
    if (!s) return nullptr;
    std::unique_lock<std::mutex> lk(s->m);
    if (wait) s->cv_free.wait(lk, [&] { return !s->free_list.empty() || s->workers_alive == 0; });
    if (s->free_list.empty()) return nullptr;
    uint8_t* buf = s->free_list.back();
    s->free_list.pop_back();
    return buf;
}

void submit(Sink* s, uint8_t* frame, int index) {
    if (!s || !frame) return;
    {
        std::lock_guard<std::mutex> lk(s->m);
        s->queue.push_back({frame, index});
        s->submitted += 1;
    }
    s->cv_work.notify_one();
}

void release(Sink* s, uint8_t* frame) {
    if (!s || !frame) return;
    {
        std::lock_guard<std::mutex> lk(s->m);
        s->free_list.push_back(frame);
    }
    s->cv_free.notify_one();
}

void close(Sink* s) {
    if (!s) return;
    {
        std::lock_guard<std::mutex> lk(s->m);
        s->closing = true;
    }
    s->cv_work.notify_all();
}

void cancel(Sink* s) {
    if (!s) return;
    {
        std::lock_guard<std::mutex> lk(s->m);
        for (const Job& j : s->queue) s->free_list.push_back(j.data);
        s->queue.clear();
        s->closing = true;
    }
    s->cv_work.notify_all();
    s->cv_free.notify_all();
}

Status status(Sink* s) {
    Status st = {};
    if (!s) {
        st.done = true;
        return st;
    }
    std::lock_guard<std::mutex> lk(s->m);
    st.submitted = s->submitted;
    st.written = s->written;
    st.failed = s->failed;
    st.queued = s->submitted - s->written - s->failed;
    st.done = s->closing && s->workers_alive == 0;
    st.ok = s->ok;
    snprintf(st.message, sizeof(st.message), "%s", s->message);
    return st;
}

void destroy(Sink* s) {
    if (!s) return;
    bool closed;
    {
        std::lock_guard<std::mutex> lk(s->m);
        closed = s->closing;
    }
    // A sink that was closed is draining what it was given, which is waited for rather than dropped
    if (!closed) cancel(s);
    for (std::thread& t : s->threads) {
        if (t.joinable()) t.join();
    }
    for (uint8_t* b : s->all) free(b);
    delete s;
}

}  // namespace frame_sink
