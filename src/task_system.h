#pragma once

#include <core/md_str.h>
#include <core/md_executor.h>

#include <stdint.h>
#include <stddef.h>
#include <atomic>
#include <functional>
#include <memory>
#include <type_traits>

// The application's CPU threads: one pool, owned here, which everything else borrows.
//
// Three ways to use it, by the shape of the work:
//
//   parallel_for   Blocks the caller, which works on the range too. Fork-join, and the bulk of the uses.
//                  From the main thread it only ever helps with work of its own priority, so a frame never
//                  picks up a piece of a background job while it waits.
//
//   Job            Background work which outlives the call that starts it. Owned by the caller, like any
//                  member: no ids, no table, no cap. Its destructor cancels and waits, so a job never runs
//                  past the life of what it works on. A follow-up on the main thread (on_main) runs in
//                  run_main_queue() and nowhere else. Two shapes of work:
//                    start      a body which runs to its end on a pool thread: it may block, and parallel_for
//                    start_steps  work which mostly waits (on the GPU): a step at a time, each returning at once,
//                               called by poll() until the work is done. It holds no thread while it waits.
//
//   detach         Fire and forget, owned here: never cancelled, and waited for by wait_all and shutdown.
//
// Cancellation reaches the work: inside a job body, and inside every parallel_for it issues on whatever
// thread, cancelled() says whether the job was asked to stop, and a parallel_for stops handing out pieces.
//
// The overlay lists Activities: every job with a label is one, and anything else that runs over time
// (work on the GPU over several frames) can register one as well.

namespace task_system {

void     initialize(uint32_t num_threads);  // including the main thread
void     shutdown();                        // ends the step jobs, waits for everything else still under way

uint32_t num_threads();
uint32_t thread_index();                    // [0, num_threads()), 0 is the main thread

// ---------------------------------------------------------------------------------------------------------------------
// Fork-join

typedef void (*RangeFn)(void* user, uint32_t beg, uint32_t end);

// Calls fn over [0, count) in half open pieces of at least grain (but the last), in any order, on any thread,
// the caller one of them. Returns false if the job it runs in was cancelled, and some pieces were skipped.
bool parallel_for(uint32_t count, uint32_t grain, RangeFn fn, void* user);

template <class F>
bool parallel_for(uint32_t count, uint32_t grain, F&& f) {
    using T = std::remove_reference_t<F>;
    return parallel_for(count, grain, [](void* user, uint32_t beg, uint32_t end) { (*static_cast<T*>(user))(beg, end); },
                        const_cast<void*>(static_cast<const void*>(std::addressof(f))));
}

template <class F>
bool parallel_for(uint32_t count, F&& f) { return parallel_for(count, 1, static_cast<F&&>(f)); }

// The pool as an executor, for what mdlib spreads over threads (md_executor.h): parallel_for underneath, with
// thread_index() as thread_idx. Its loops run at the priority of whoever issues them, and stop with their job. NULL
// before initialize(): mdlib then runs the loops on the calling thread.
const md_executor_t* executor();

// ---------------------------------------------------------------------------------------------------------------------
// Inside a job, and inside any parallel_for it issues (on any thread). Outside of a job: false, and nothing.

bool cancelled();
// The same, as a flag for an API which reads one (mdlib descriptors), to be read atomically (md_atomic_load32).
// NULL outside of a job.
const volatile int32_t* cancel_flag();
// How far the job is, 0 to 1 (a negative value: unknown). start_for reports it by itself.
void set_progress(float fraction);
// The same, as a value for an API which writes one (mdlib descriptors). NULL outside of a job.
volatile float* progress_flag();

// ---------------------------------------------------------------------------------------------------------------------
// Jobs

enum class Status : uint8_t {
    Idle,           // never started
    Running,        // under way, or its on_main has yet to run
    Completed,
    Failed,         // a step job whose step said so: its work could not be done
    Cancelled,      // asked to stop, whether before it started, while it ran or after: its results are not to be trusted
};

struct JobState;

class Job {
public:
    Job() = default;
    ~Job();                 // cancel + wait: a job never outlives its owner. A pending on_main is dropped.
    Job(const Job&) = delete;
    Job& operator=(const Job&) = delete;

    using Body      = std::function<void()>;
    using RangeBody = std::function<void(uint32_t beg, uint32_t end)>;
    using OnMain    = std::function<void(Status)>;
    using Step      = std::function<Status(double budget_ms)>;
    using End       = std::function<void()>;

    // An empty label keeps it out of the overlay. Restarting a job which is still under way cancels and waits for
    // it first, and drops its pending on_main.
    void start    (md_str_t label, Body body, OnMain on_main = {});
    void start_for(md_str_t label, uint32_t count, uint32_t grain, RangeBody body, OnMain on_main = {});

    // Work which mostly waits, a step at a time: step is called by poll() until it returns Completed or Failed, never
    // two calls at once, and returns within budget_ms (< 0: it runs the work to its end). It takes its first step
    // here, so that work which is done at once ends at once. Inside it cancelled() and set_progress() work as in a
    // body. A step touches its own job's state only: the thread it runs on is poll()'s. end releases what the work
    // holds, however it ended (done, failed, cancelled), and never waits. The main thread only, for now: the steps
    // run there (poll), and the work they drive shares its GPU stream.
    void start_steps(md_str_t label, Step step, End end = {}, OnMain on_main = {});

    Status status() const;
    bool   running() const { return status() == Status::Running; }
    float  progress() const;    // 0 to 1, negative if unknown

    void   cancel();            // a step job ends at once (end), never waiting; a body sees it in cancelled()
    void   wait();              // for the work (a step job is run to its end here); an on_main still waits for
                                // run_main_queue
    void   cancel_and_wait() { cancel(); wait(); }

    // To recognize this job among the activities
    uint64_t activity_id() const;

private:
    void release();
    JobState* m_state = nullptr;
};

// Runs whatever happens: cancel_all does not reach it (it is for work which must be done, such as giving memory
// back), and wait_all waits for it.
void detach(md_str_t label, std::function<void()> body);

// Every job on the pool at once: a load stops all work which reads what it is about to replace. Step jobs are not
// reached: they hold device memory their owners free in an order of their own, and end with them.
void cancel_all();
void wait_all();

// The step jobs advanced: once a frame, from the main thread, budget_ms shared by all of them (a job in its first
// frames is given more, so that a short one ends within them).
void poll(double budget_ms);

// The main thread's follow-ups (Job on_main). Once per frame, and the only place they run.
void run_main_queue();

// ---------------------------------------------------------------------------------------------------------------------
// Activities: what the overlay shows

constexpr size_t ACTIVITY_LABEL_SIZE = 64;

struct Activity {
    char                 label[ACTIVITY_LABEL_SIZE] = "";
    std::atomic<float>   progress{-1.0f};   // negative: unknown
    std::atomic<int32_t> cancel{0};         // set when it is asked to stop
    bool                 cancellable = true;
    int64_t              t0 = 0;            // md_tick_t: when it began (activity_begin)
    Activity*            prev = nullptr;    // the list, under its lock
    Activity*            next = nullptr;
    bool                 listed = false;
};

void activity_begin(Activity* a, md_str_t label);
void activity_end(Activity* a);

struct ActivityInfo {
    uint64_t id;
    char     label[ACTIVITY_LABEL_SIZE];
    float    progress;
    bool     cancellable;
    double   seconds;       // since it began: what ends within a few frames is better not shown at all
};

// A copy of what is under way: labels are not read while the work ends
size_t activity_list(ActivityInfo* out, size_t cap);
// Asks the activity to stop, if it is still under way
void   activity_cancel(uint64_t id);

}  // namespace task_system
