#include "task_system.h"

#include <TaskScheduler.h>

#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <deque>
#include <mutex>
#include <vector>

// enkiTS underneath. What is used of it is small on purpose: ITaskSet, its priority and WaitforTask, and nothing of
// its dependencies or pinned tasks. Pinned tasks run whenever the main thread waits for anything, which is not a
// place for the main thread's follow-ups; those have a queue of their own here.

namespace task_system {

// Work issued from the main thread (or from within such work) runs at the frame's priority, jobs below it. A wait
// runs nothing of a lower priority than what it waits for, so a frame never picks up a piece of a background job.
static constexpr enki::TaskPriority PRIORITY_FRAME      = enki::TASK_PRIORITY_HIGH;
static constexpr enki::TaskPriority PRIORITY_BACKGROUND = enki::TASK_PRIORITY_LOW;

static enki::TaskScheduler g_ts;
static bool g_initialized = false;

// What the calling thread works for: the job (if any) and the priority its pieces of work are issued at. Carried
// into every piece of a parallel_for, on whatever thread it runs.
struct Context {
    JobState*          job  = nullptr;
    enki::TaskPriority prio = PRIORITY_FRAME;
};
static thread_local Context tl_ctx;

static_assert(sizeof(std::atomic<int32_t>) == sizeof(int32_t) && std::atomic<int32_t>::is_always_lock_free,
              "cancel_flag() hands out the atomic as a plain flag");
static_assert(sizeof(std::atomic<float>) == sizeof(float) && std::atomic<float>::is_always_lock_free,
              "progress_flag() hands out the atomic as a plain value");

// ---------------------------------------------------------------------------------------------------------------------
// Activities

static std::mutex g_activity_mutex;
static Activity*  g_activity_head = nullptr;
static Activity*  g_activity_tail = nullptr;

void activity_begin(Activity* a, md_str_t label) {
    MD_ASSERT(a);
    md_str_copy_to_char_buf(a->label, sizeof(a->label), label);
    a->t0 = (int64_t)md_tick_now();
    std::lock_guard<std::mutex> lock(g_activity_mutex);
    MD_ASSERT(!a->listed);
    a->prev = g_activity_tail;
    a->next = nullptr;
    if (g_activity_tail) g_activity_tail->next = a;
    else                 g_activity_head = a;
    g_activity_tail = a;
    a->listed = true;
}

void activity_end(Activity* a) {
    MD_ASSERT(a);
    std::lock_guard<std::mutex> lock(g_activity_mutex);
    if (!a->listed) return;
    if (a->prev) a->prev->next = a->next;
    else         g_activity_head = a->next;
    if (a->next) a->next->prev = a->prev;
    else         g_activity_tail = a->prev;
    a->prev = a->next = nullptr;
    a->listed = false;
}

size_t activity_list(ActivityInfo* out, size_t cap) {
    const md_tick_t now = md_tick_now();
    std::lock_guard<std::mutex> lock(g_activity_mutex);
    size_t n = 0;
    for (const Activity* a = g_activity_head; a && n < cap; a = a->next) {
        ActivityInfo& info = out[n++];
        info.id = (uint64_t)(uintptr_t)a;
        md_memcpy(info.label, a->label, sizeof(info.label));
        info.progress = a->progress.load(std::memory_order_relaxed);
        info.cancellable = a->cancellable;
        info.seconds = md_tick_to_seconds(now - (md_tick_t)a->t0);
    }
    return n;
}

void activity_cancel(uint64_t id) {
    std::lock_guard<std::mutex> lock(g_activity_mutex);
    for (Activity* a = g_activity_head; a; a = a->next) {
        if ((uint64_t)(uintptr_t)a == id) {
            if (a->cancellable) a->cancel.store(1, std::memory_order_relaxed);
            return;
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Jobs

struct JobTask : enki::ITaskSet {
    JobState* state = nullptr;
    void ExecuteRange(enki::TaskSetPartition, uint32_t) override;
};

struct JobState {
    JobTask             task;
    Activity            activity;          // its cancel is the job's cancel, listed only with a label
    Job::Body           body;
    Job::OnMain         on_main;
    std::atomic<Status> result{Status::Running};
    bool                main_pending = false;   // under g_main_mutex
    bool                labelled = false;
    bool                detached = false;  // runs whatever happens: not reached by cancel_all
    JobState*           prev = nullptr;    // every job under way, under g_job_mutex
    JobState*           next = nullptr;

    // A step job: run by poll(), never by enki (the main thread only)
    bool                is_step = false;
    bool                stepping = false;  // under way: in the step list
    Job::Step           step;
    Job::End            end;
    int64_t             t0 = 0;            // md_tick_t: when it began
    JobState*           step_prev = nullptr;
    JobState*           step_next = nullptr;
};

static std::mutex             g_job_mutex;
static JobState*              g_jobs = nullptr;    // under way: what cancel_all reaches
static std::vector<JobState*> g_detached;          // owned here, freed once enki is done with them

static std::mutex             g_main_mutex;
static std::deque<JobState*>  g_main_queue;

static inline bool job_cancelled(const JobState* s) {
    return s->activity.cancel.load(std::memory_order_relaxed) != 0;
}

static void jobs_link(JobState* s) {
    std::lock_guard<std::mutex> lock(g_job_mutex);
    s->prev = nullptr;
    s->next = g_jobs;
    if (g_jobs) g_jobs->prev = s;
    g_jobs = s;
}

static void jobs_unlink(JobState* s) {
    std::lock_guard<std::mutex> lock(g_job_mutex);
    if (s->prev) s->prev->next = s->next;
    else         g_jobs = s->next;
    if (s->next) s->next->prev = s->prev;
    s->prev = s->next = nullptr;
}

// Once enki has stored the task's completion it touches it no more (the tasks have no dependents), so a complete
// task can be freed.
static void reap_detached() {
    std::lock_guard<std::mutex> lock(g_job_mutex);
    for (size_t i = 0; i < g_detached.size();) {
        if (g_detached[i]->task.GetIsComplete()) {
            delete g_detached[i];
            g_detached[i] = g_detached.back();
            g_detached.pop_back();
        } else {
            ++i;
        }
    }
}

void JobTask::ExecuteRange(enki::TaskSetPartition, uint32_t) {
    JobState* s = state;
    if (!job_cancelled(s)) {
        const Context saved = tl_ctx;
        tl_ctx = {s, PRIORITY_BACKGROUND};
        s->body();
        tl_ctx = saved;
    }
    // Its captures go now, on this thread, rather than whenever the owner next starts or ends it
    s->body = nullptr;

    if (s->labelled) activity_end(&s->activity);
    if (!s->detached) jobs_unlink(s);

    // Before the task completes: whoever looks at it from now on sees the result, or a pending on_main
    s->result.store(job_cancelled(s) ? Status::Cancelled : Status::Completed, std::memory_order_release);
    if (s->on_main) {
        std::lock_guard<std::mutex> lock(g_main_mutex);
        s->main_pending = true;
        g_main_queue.push_back(s);
    }
}

static JobState* job_create(md_str_t label, Job::Body body, Job::OnMain on_main, bool detached) {
    JobState* s = new JobState();
    s->task.state      = s;
    s->task.m_SetSize  = 1;
    s->task.m_Priority = PRIORITY_BACKGROUND;
    s->body     = std::move(body);
    s->on_main  = std::move(on_main);
    s->labelled = label.len > 0;
    s->detached = detached;
    s->activity.cancellable = !detached;
    if (s->labelled) activity_begin(&s->activity, label);
    if (!detached) jobs_link(s);
    return s;
}

static void job_launch(JobState* s) {
    if (g_initialized) {
        g_ts.AddTaskSetToPipe(&s->task);
    } else {
        // No pool (a tool or a test without one): there and then
        s->task.ExecuteRange({0, 1}, 0);
    }
}

static void job_wait(JobState* s) {
    while (g_initialized && !s->task.GetIsComplete()) {
        g_ts.WaitforTask(&s->task);
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Step jobs: the main thread's, so the list takes no lock

// A job in its first frames is given this much of poll's budget, so that one done within it ends at once rather than a
// frame later for each step it takes; after that, STEP_BUDGET_MS a frame
static constexpr double STEP_YOUNG_MS        = 50.0;
static constexpr double STEP_YOUNG_BUDGET_MS = 4.0;
static constexpr double STEP_BUDGET_MS       = 1.0;

static JobState* g_steps_head = nullptr;     // under way, oldest first
static JobState* g_steps_tail = nullptr;

static void steps_unlink(JobState* s) {
    if (s->step_prev) s->step_prev->step_next = s->step_next;
    else              g_steps_head = s->step_next;
    if (s->step_next) s->step_next->step_prev = s->step_prev;
    else              g_steps_tail = s->step_prev;
    s->step_prev = s->step_next = nullptr;
}

// No longer under way: what it holds released, out of the overlay, its result there and its on_main queued
static void step_finish(JobState* s, Status st) {
    steps_unlink(s);
    s->stepping = false;
    if (job_cancelled(s)) st = Status::Cancelled;
    if (s->end) s->end();
    // its captures go now, rather than whenever the owner next starts or ends it
    s->end  = nullptr;
    s->step = nullptr;
    if (st == Status::Completed) s->activity.progress.store(1.0f, std::memory_order_relaxed);
    if (s->labelled) activity_end(&s->activity);
    s->result.store(st, std::memory_order_release);
    if (s->on_main) {
        std::lock_guard<std::mutex> lock(g_main_mutex);
        s->main_pending = true;
        g_main_queue.push_back(s);
    }
}

static void step_run(JobState* s, double budget_ms) {
    if (job_cancelled(s)) {
        step_finish(s, Status::Cancelled);
        return;
    }
    const Context saved = tl_ctx;
    tl_ctx = {s, PRIORITY_FRAME};
    const Status st = s->step(budget_ms);
    tl_ctx = saved;
    if (st != Status::Running) step_finish(s, st);
}

static double step_budget(const JobState* s, md_tick_t now) {
    return md_tick_to_milliseconds(now - (md_tick_t)s->t0) < STEP_YOUNG_MS ? STEP_YOUNG_BUDGET_MS : STEP_BUDGET_MS;
}

void poll(double budget_ms) {
    const md_tick_t t0 = md_tick_now();
    for (JobState* s = g_steps_head; s;) {
        JobState* next = s->step_next;      // it may finish, and leave the list
        const md_tick_t now = md_tick_now();
        const double left = MD_MAX(budget_ms - md_tick_to_milliseconds(now - t0), 0.0);
        step_run(s, MD_MIN(step_budget(s, now), left));
        s = next;
    }
}

static void steps_end_all() {
    while (g_steps_head) step_finish(g_steps_head, Status::Cancelled);
}

Job::~Job() { release(); }

void Job::release() {
    if (!m_state) return;
    cancel();
    if (!m_state->is_step) job_wait(m_state);
    {
        std::lock_guard<std::mutex> lock(g_main_mutex);
        if (m_state->main_pending) {
            for (auto it = g_main_queue.begin(); it != g_main_queue.end(); ++it) {
                if (*it == m_state) {
                    g_main_queue.erase(it);
                    break;
                }
            }
            m_state->main_pending = false;
        }
    }
    delete m_state;
    m_state = nullptr;
}

void Job::start(md_str_t label, Body body, OnMain on_main) {
    MD_ASSERT(body);
    release();
    m_state = job_create(label, std::move(body), std::move(on_main), false);
    job_launch(m_state);
}

void Job::start_steps(md_str_t label, Step step, End end, OnMain on_main) {
    MD_ASSERT(step);
    release();
    JobState* s = new JobState();
    s->is_step  = true;
    s->stepping = true;
    s->step     = std::move(step);
    s->end      = std::move(end);
    s->on_main  = std::move(on_main);
    s->t0       = (int64_t)md_tick_now();
    s->labelled = label.len > 0;
    s->activity.cancellable = true;
    if (s->labelled) activity_begin(&s->activity, label);
    s->step_prev = g_steps_tail;
    if (g_steps_tail) g_steps_tail->step_next = s;
    else              g_steps_head = s;
    g_steps_tail = s;
    m_state = s;
    step_run(s, STEP_YOUNG_BUDGET_MS);
}

void Job::start_for(md_str_t label, uint32_t count, uint32_t grain, RangeBody body, OnMain on_main) {
    MD_ASSERT(body);
    auto run = [count, grain, body = std::move(body)]() {
        JobState* s = tl_ctx.job;
        MD_ASSERT(s);
        s->activity.progress.store(0.0f, std::memory_order_relaxed);
        std::atomic<uint32_t> done{0};
        parallel_for(count, grain, [&](uint32_t beg, uint32_t end) {
            body(beg, end);
            const uint32_t d = done.fetch_add(end - beg, std::memory_order_relaxed) + (end - beg);
            s->activity.progress.store((float)d / (float)MD_MAX(count, 1u), std::memory_order_relaxed);
        });
    };
    start(label, std::move(run), std::move(on_main));
}

Status Job::status() const {
    if (!m_state) return Status::Idle;
    if (m_state->is_step ? m_state->stepping : !m_state->task.GetIsComplete()) return Status::Running;
    {
        std::lock_guard<std::mutex> lock(g_main_mutex);
        if (m_state->main_pending) return Status::Running;
    }
    return m_state->result.load(std::memory_order_acquire);
}

float Job::progress() const {
    return m_state ? m_state->activity.progress.load(std::memory_order_relaxed) : -1.0f;
}

void Job::cancel() {
    if (!m_state) return;
    m_state->activity.cancel.store(1, std::memory_order_relaxed);
    if (m_state->is_step && m_state->stepping) step_finish(m_state, Status::Cancelled);
}

void Job::wait() {
    if (!m_state) return;
    if (m_state->is_step) {
        while (m_state->stepping) step_run(m_state, -1.0);
    } else {
        job_wait(m_state);
    }
}

uint64_t Job::activity_id() const {
    return m_state ? (uint64_t)(uintptr_t)&m_state->activity : 0;
}

void detach(md_str_t label, std::function<void()> body) {
    MD_ASSERT(body);
    reap_detached();
    JobState* s = job_create(label, std::move(body), {}, true);
    {
        std::lock_guard<std::mutex> lock(g_job_mutex);
        g_detached.push_back(s);
    }
    job_launch(s);
}

void cancel_all() {
    std::lock_guard<std::mutex> lock(g_job_mutex);
    for (JobState* s = g_jobs; s; s = s->next) {
        s->activity.cancel.store(1, std::memory_order_relaxed);
    }
}

void wait_all() {
    if (g_initialized) g_ts.WaitforAll();
    reap_detached();
}

void run_main_queue() {
    // One at a time, each taken off the queue before it runs: a follow-up may end the owner of another one
    size_t n;
    {
        std::lock_guard<std::mutex> lock(g_main_mutex);
        n = g_main_queue.size();
    }
    for (size_t i = 0; i < n; ++i) {
        Job::OnMain fn;
        Status st;
        {
            std::lock_guard<std::mutex> lock(g_main_mutex);
            if (g_main_queue.empty()) break;
            JobState* s = g_main_queue.front();
            g_main_queue.pop_front();
            s->main_pending = false;
            fn = std::move(s->on_main);
            s->on_main = nullptr;
            st = s->result.load(std::memory_order_acquire);
        }
        // s is not touched from here: fn may well end its owner
        fn(st);
    }
    reap_detached();
}

// ---------------------------------------------------------------------------------------------------------------------
// Fork-join

struct ParallelSet : enki::ITaskSet {
    RangeFn           fn = nullptr;
    void*             user = nullptr;
    Context           ctx = {};
    uint32_t          count = 0;
    uint32_t          grain = 1;
    std::atomic<bool> skipped{false};

    // enki's range is in grains: its m_MinRange bounds where it splits, not the remainder each split leaves behind
    // (which can be one), so the pieces of at least grain are the ones it is handed in whole grains
    void ExecuteRange(enki::TaskSetPartition range, uint32_t) override {
        if (ctx.job && job_cancelled(ctx.job)) {
            skipped.store(true, std::memory_order_relaxed);
            return;
        }
        const uint32_t beg = (uint32_t)MD_MIN((uint64_t)range.start * grain, (uint64_t)count);
        const uint32_t end = (uint32_t)MD_MIN((uint64_t)range.end * grain, (uint64_t)count);
        const Context saved = tl_ctx;
        tl_ctx = ctx;
        fn(user, beg, end);
        tl_ctx = saved;
    }
};

bool parallel_for(uint32_t count, uint32_t grain, RangeFn fn, void* user) {
    MD_ASSERT(fn);
    const Context ctx = tl_ctx;
    if (count == 0) return !(ctx.job && job_cancelled(ctx.job));
    grain = MD_MAX(grain, 1u);

    const bool serial = !g_initialized || count <= grain || g_ts.GetNumTaskThreads() < 2 || g_ts.GetThreadNum() == enki::NO_THREAD_NUM;
    if (serial) {
        if (!ctx.job) {
            fn(user, 0, count);
            return true;
        }
        // In pieces still, for a cancel to be seen between them
        const uint32_t step = MD_MAX(grain, count / 64);
        for (uint32_t beg = 0; beg < count; beg += step) {
            if (job_cancelled(ctx.job)) return false;
            fn(user, beg, MD_MIN(beg + step, count));
        }
        return true;
    }

    ParallelSet set;
    set.m_SetSize  = (uint32_t)(((uint64_t)count + grain - 1) / grain);
    set.m_MinRange = 1;
    set.count = count;
    set.grain = grain;
    set.m_Priority = ctx.prio;
    set.fn   = fn;
    set.user = user;
    set.ctx  = ctx;
    g_ts.AddTaskSetToPipe(&set);
    g_ts.WaitforTask(&set, ctx.prio);
    return !set.skipped.load(std::memory_order_relaxed);
}

// A piece of an executor loop is under way on this thread: md_executor.h forbids it to issue or wait for parallel
// work, since a thread which waits may run another piece of the same loop with the same thread_idx
static thread_local bool tl_in_executor_piece = false;

static bool executor_parallel_for(void* user, uint32_t count, uint32_t grain, md_range_fn body, void* arg) {
    (void)user;
    return parallel_for(count, grain, [body, arg](uint32_t beg, uint32_t end) {
        MD_ASSERT(!tl_in_executor_piece && "a piece of an executor loop issued or waited for parallel work (md_executor.h)");
        tl_in_executor_piece = true;
        body(arg, beg, end, thread_index());
        tl_in_executor_piece = false;
    });
}

static md_executor_t g_executor = {};

const md_executor_t* executor() {
    return g_initialized ? &g_executor : nullptr;
}

bool cancelled() {
    return tl_ctx.job && job_cancelled(tl_ctx.job);
}

const volatile int32_t* cancel_flag() {
    return tl_ctx.job ? reinterpret_cast<const volatile int32_t*>(&tl_ctx.job->activity.cancel) : nullptr;
}

void set_progress(float fraction) {
    if (tl_ctx.job) tl_ctx.job->activity.progress.store(fraction, std::memory_order_relaxed);
}

volatile float* progress_flag() {
    return tl_ctx.job ? reinterpret_cast<volatile float*>(&tl_ctx.job->activity.progress) : nullptr;
}

// ---------------------------------------------------------------------------------------------------------------------

void initialize(uint32_t num_threads) {
    MD_ASSERT(!g_initialized);
    g_ts.Initialize(num_threads);
    g_initialized = true;
    g_executor = {executor_parallel_for, nullptr, g_ts.GetNumTaskThreads()};
}

void shutdown() {
    steps_end_all();
    if (!g_initialized) return;
    wait_all();
    g_ts.WaitforAllAndShutdown();
    g_initialized = false;
    // What waited for the main thread is not going to get it
    std::lock_guard<std::mutex> lock(g_main_mutex);
    for (JobState* s : g_main_queue) s->main_pending = false;
    g_main_queue.clear();
}

uint32_t num_threads() {
    return g_initialized ? g_ts.GetNumTaskThreads() : 1;
}

uint32_t thread_index() {
    if (!g_initialized) return 0;
    const uint32_t i = g_ts.GetThreadNum();
    return i == enki::NO_THREAD_NUM ? 0 : i;
}

}  // namespace task_system
