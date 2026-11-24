#include "utest.h"

#include <task_system.h>

#include <core/md_os.h>

#include <algorithm>
#include <atomic>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

// task_system owns a global thread pool with no re-entrant initialize, so it is brought up once on
// first use and left running for the life of the process.

namespace {

void ensure_task_system() {
    static bool initialized = false;
    if (!initialized) {
        if (!task_system::executor()) {     // another test may have brought it up (test_density_cache.cpp)
            task_system::initialize(4);
            atexit([]() { task_system::shutdown(); });
        }
        initialized = true;
    }
}

// Polls without waiting through the scheduler, which would run pool work on this thread meanwhile
void poll_until_stopped(const task_system::Job& job) {
    while (job.running()) {
        task_system::run_main_queue();
        md_thread_sleep(1);
    }
}

}  // namespace

UTEST(viamd_task_system, parallel_for_visits_every_index_once) {
    ensure_task_system();

    const uint32_t N = 10000;
    std::vector<std::atomic<int>> visits(N);
    for (uint32_t i = 0; i < N; ++i) visits[i] = 0;

    const bool ok = task_system::parallel_for(N, 1, [&](uint32_t beg, uint32_t end) {
        for (uint32_t i = beg; i < end; ++i) visits[i] += 1;
    });
    EXPECT_TRUE(ok);

    uint32_t wrong = 0;
    for (uint32_t i = 0; i < N; ++i) {
        if (visits[i].load() != 1) wrong += 1;
    }
    EXPECT_EQ(0u, wrong);
}

UTEST(viamd_task_system, ranges_partition_the_whole_range) {
    // Each invocation gets a half open interval, and together they must tile [0, N) exactly. An
    // overlap would double count and a gap would silently drop work.
    ensure_task_system();

    const uint32_t N = 4096;
    const uint32_t GRAIN = 64;
    std::mutex mutex;
    std::vector<std::pair<uint32_t, uint32_t>> ranges;

    task_system::parallel_for(N, GRAIN, [&](uint32_t beg, uint32_t end) {
        std::lock_guard<std::mutex> lock(mutex);
        ranges.push_back({beg, end});
    });

    std::sort(ranges.begin(), ranges.end());
    ASSERT_TRUE(!ranges.empty());
    EXPECT_EQ(0u, ranges.front().first);
    EXPECT_EQ(N,  ranges.back().second);
    for (size_t i = 1; i < ranges.size(); ++i) {
        EXPECT_EQ(ranges[i - 1].second, ranges[i].first);
    }
    for (size_t i = 0; i < ranges.size(); ++i) {
        EXPECT_LT(ranges[i].first, ranges[i].second);
        // the grain is a minimum, but for the last piece
        if (ranges[i].second != N) EXPECT_GE(ranges[i].second - ranges[i].first, GRAIN);
    }
}

UTEST(viamd_task_system, a_job_completes) {
    ensure_task_system();

    std::atomic<int> ran(0);
    task_system::Job job;
    EXPECT_TRUE(job.status() == task_system::Status::Idle);
    job.start(MD_STR_LIT("done"), [&]() { ran = 1; });
    job.wait();

    EXPECT_EQ(1, ran.load());
    EXPECT_TRUE(job.status() == task_system::Status::Completed);
    EXPECT_FALSE(job.running());
}

UTEST(viamd_task_system, a_job_cancelled_before_it_starts_never_runs) {
    ensure_task_system();

    // Every worker is held, so the job below cannot start before it is cancelled
    std::atomic<bool> release(false);
    const uint32_t workers = task_system::num_threads();
    std::vector<std::unique_ptr<task_system::Job>> blockers;
    std::atomic<uint32_t> blocking(0);
    for (uint32_t i = 0; i < workers; ++i) {
        blockers.push_back(std::make_unique<task_system::Job>());
        blockers.back()->start({}, [&]() {
            blocking += 1;
            while (!release) md_thread_sleep(1);
        });
    }
    while (blocking.load() < workers - 1) md_thread_sleep(1);   // the main thread is not one of them

    std::atomic<int> ran(0);
    task_system::Job job;
    job.start({}, [&]() { ran = 1; });
    job.cancel();
    release = true;
    poll_until_stopped(job);
    for (auto& b : blockers) b->wait();

    EXPECT_EQ(0, ran.load());
    EXPECT_TRUE(job.status() == task_system::Status::Cancelled);
}

UTEST(viamd_task_system, cancel_reaches_the_body_and_its_parallel_for) {
    ensure_task_system();

    std::atomic<bool> started(false);
    std::atomic<bool> saw_cancel(false);
    std::atomic<int>  loop_result(-1);
    std::atomic<uint32_t> pieces(0);

    task_system::Job job;
    job.start(MD_STR_LIT("cancel me"), [&]() {
        started = true;
        while (!task_system::cancelled()) md_thread_sleep(1);
        saw_cancel = task_system::cancel_flag() != nullptr && task_system::cancelled();
        // Issued after the cancel: none of its pieces run
        const bool ok = task_system::parallel_for(1000, 1, [&](uint32_t, uint32_t) { pieces += 1; });
        loop_result = ok ? 1 : 0;
    });
    while (!started) md_thread_sleep(1);
    job.cancel();
    job.wait();

    EXPECT_TRUE(saw_cancel.load());
    EXPECT_EQ(0, loop_result.load());
    EXPECT_EQ(0u, pieces.load());
    EXPECT_TRUE(job.status() == task_system::Status::Cancelled);
}

UTEST(viamd_task_system, on_main_runs_in_the_main_queue_only) {
    ensure_task_system();

    std::atomic<int> calls(0);
    task_system::Status got = task_system::Status::Idle;
    task_system::Job job;
    job.start({}, []() {}, [&](task_system::Status st) { calls += 1; got = st; });
    job.wait();

    // The work is done, the follow-up is not: and no wait on this thread runs it
    EXPECT_TRUE(job.running());
    for (int i = 0; i < 20; ++i) {
        task_system::parallel_for(256, 1, [](uint32_t, uint32_t) {});
    }
    EXPECT_EQ(0, calls.load());

    task_system::run_main_queue();
    EXPECT_EQ(1, calls.load());
    EXPECT_TRUE(got == task_system::Status::Completed);
    EXPECT_TRUE(job.status() == task_system::Status::Completed);
}

UTEST(viamd_task_system, an_owner_may_end_in_its_own_on_main) {
    ensure_task_system();

    struct Owner { task_system::Job job; int value = 0; };
    std::atomic<int> seen(0);
    Owner* owner = new Owner();
    owner->value = 7;
    owner->job.start({}, []() {}, [owner, &seen](task_system::Status) {
        seen = owner->value;
        delete owner;
    });
    owner->job.wait();
    task_system::run_main_queue();
    EXPECT_EQ(7, seen.load());
}

UTEST(viamd_task_system, an_ended_owner_drops_its_pending_on_main) {
    ensure_task_system();

    std::atomic<int> calls(0);
    {
        task_system::Job job;
        job.start({}, []() {}, [&](task_system::Status) { calls += 1; });
        job.wait();
    }
    task_system::run_main_queue();
    EXPECT_EQ(0, calls.load());
}

UTEST(viamd_task_system, the_destructor_cancels_and_waits) {
    ensure_task_system();

    std::atomic<bool> started(false);
    std::atomic<bool> finished(false);
    {
        task_system::Job job;
        job.start({}, [&]() {
            started = true;
            while (!task_system::cancelled()) md_thread_sleep(1);
            finished = true;
        });
        while (!started) md_thread_sleep(1);
    }
    EXPECT_TRUE(finished.load());
}

UTEST(viamd_task_system, the_main_thread_never_runs_background_work_in_a_frame_wait) {
    ensure_task_system();

    // Background jobs spread their work over the pool (as the script's evaluation does), so pieces of it sit in
    // the workers' queues: where a waiting main thread would otherwise steal from
    const int NUM_JOBS = 3;
    std::atomic<int> on_main_thread(0);
    std::vector<std::unique_ptr<task_system::Job>> jobs;
    for (int i = 0; i < NUM_JOBS; ++i) {
        jobs.push_back(std::make_unique<task_system::Job>());
        jobs.back()->start({}, [&]() {
            task_system::parallel_for(64, 1, [&](uint32_t, uint32_t) {
                if (task_system::thread_index() == 0) on_main_thread += 1;
                md_thread_sleep(1);
            });
        });
    }
    // Frame work meanwhile, in pieces the workers take some of: the main thread then waits with its own gone
    for (int i = 0; i < 100; ++i) {
        task_system::parallel_for(16, 1, [](uint32_t, uint32_t) { md_thread_sleep(1); });
    }
    for (auto& j : jobs) poll_until_stopped(*j);
    EXPECT_EQ(0, on_main_thread.load());
}

UTEST(viamd_task_system, start_for_reports_progress) {
    ensure_task_system();

    const uint32_t N = 1000;
    std::atomic<uint32_t> sum(0);
    task_system::Job job;
    job.start_for(MD_STR_LIT("range"), N, 16, [&](uint32_t beg, uint32_t end) { sum += end - beg; });
    job.wait();
    EXPECT_EQ(N, sum.load());
    EXPECT_NEAR(1.0f, job.progress(), 1.0e-6f);
}

UTEST(viamd_task_system, a_labelled_job_is_listed_while_it_runs) {
    ensure_task_system();

    std::atomic<bool> release(false);
    std::atomic<bool> started(false);
    task_system::Job job;
    job.start(MD_STR_LIT("Listed"), [&]() {
        started = true;
        while (!release) md_thread_sleep(1);
    });
    while (!started) md_thread_sleep(1);

    task_system::ActivityInfo info[16];
    size_t n = task_system::activity_list(info, 16);
    bool found = false;
    for (size_t i = 0; i < n; ++i) {
        if (info[i].id == job.activity_id()) {
            found = true;
            EXPECT_STREQ("Listed", info[i].label);
            task_system::activity_cancel(info[i].id);
        }
    }
    EXPECT_TRUE(found);
    release = true;
    job.wait();
    EXPECT_TRUE(job.status() == task_system::Status::Cancelled);   // through the activity

    n = task_system::activity_list(info, 16);
    for (size_t i = 0; i < n; ++i) EXPECT_NE(info[i].id, job.activity_id());
}

UTEST(viamd_task_system, detached_work_runs) {
    ensure_task_system();

    std::atomic<int> ran(0);
    for (int i = 0; i < 8; ++i) task_system::detach({}, [&]() { ran += 1; });
    // Not reached by a cancel of everything: what is detached must be done (memory given back)
    task_system::cancel_all();
    task_system::wait_all();
    EXPECT_EQ(8, ran.load());
}

UTEST(viamd_task_system, an_idle_job_is_answered) {
    ensure_task_system();

    task_system::Job job;
    EXPECT_FALSE(job.running());
    EXPECT_EQ(0ull, (unsigned long long)job.activity_id());
    job.wait();
    job.cancel();
}

UTEST(viamd_task_system, the_executor_keeps_its_contract) {
    ensure_task_system();

    // md_executor.h: every index once, thread_idx below num_threads and never shared by two pieces at once
    const md_executor_t* ex = task_system::executor();
    ASSERT_TRUE(ex != nullptr);
    EXPECT_EQ(task_system::num_threads(), md_executor_num_threads(ex));

    struct Ctx {
        std::vector<std::atomic<int>> hits;
        std::vector<std::atomic<int>> in_use;
        std::atomic<int> bad{0};
        uint32_t num_threads = 0;
    } c;
    const uint32_t N = 5000;
    c.num_threads = md_executor_num_threads(ex);
    c.hits = std::vector<std::atomic<int>>(N);
    c.in_use = std::vector<std::atomic<int>>(c.num_threads);
    const bool ok = md_parallel_for(ex, N, 16, [](void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
        Ctx* c = (Ctx*)arg;
        if (thread_idx >= c->num_threads) { c->bad += 1; return; }
        if (c->in_use[thread_idx].fetch_add(1) != 0) c->bad += 1;
        for (uint32_t i = beg; i < end; ++i) c->hits[i] += 1;
        c->in_use[thread_idx].fetch_sub(1);
    }, &c);
    EXPECT_TRUE(ok);
    EXPECT_EQ(0, c.bad.load());
    uint32_t wrong = 0;
    for (uint32_t i = 0; i < N; ++i) wrong += c.hits[i].load() != 1;
    EXPECT_EQ(0u, wrong);
}

UTEST(viamd_task_system, an_executor_loop_stops_with_its_job) {
    ensure_task_system();

    std::atomic<int> result(-1);
    std::atomic<uint32_t> pieces(0);
    std::atomic<bool> started(false);
    task_system::Job job;
    job.start({}, [&]() {
        started = true;
        while (!task_system::cancelled()) md_thread_sleep(1);
        // An mdlib procedure's loop, issued from a cancelled job: no piece runs, and it says so
        const bool ok = md_parallel_for(task_system::executor(), 100, 1, [](void* arg, uint32_t, uint32_t, uint32_t) {
            ((std::atomic<uint32_t>*)arg)->fetch_add(1);
        }, &pieces);
        result = ok ? 1 : 0;
    });
    while (!started) md_thread_sleep(1);
    job.cancel();
    job.wait();
    EXPECT_EQ(0, result.load());
    EXPECT_EQ(0u, pieces.load());
}

// ---------------------------------------------------------------------------------------------------------------------
// Step jobs: work which mostly waits (on the GPU), a step at a time, advanced by poll()

namespace {

// Work which is done after a number of steps, as a GPU job is once the GPU has passed its pieces
struct Counter {
    int  steps_left = 0;
    int  total = 1;
    int  steps = 0;
    int  ends = 0;
    bool fail = false;
};

void start_counter(task_system::Job* job, Counter* c, md_str_t label, task_system::Job::OnMain on_main = {}) {
    job->start_steps(label,
        [c](double budget_ms) {
            c->steps += 1;
            if (c->fail) return task_system::Status::Failed;
            if (budget_ms < 0.0) c->steps_left = 0;
            else if (c->steps_left > 0) c->steps_left -= 1;
            task_system::set_progress(1.0f - (float)c->steps_left / (float)c->total);
            return c->steps_left == 0 ? task_system::Status::Completed : task_system::Status::Running;
        },
        [c]() { c->ends += 1; }, std::move(on_main));
}

bool listed(uint64_t id, float* progress = nullptr) {
    task_system::ActivityInfo info[64];
    const size_t n = task_system::activity_list(info, 64);
    for (size_t i = 0; i < n; ++i) {
        if (info[i].id == id) {
            if (progress) *progress = info[i].progress;
            return true;
        }
    }
    return false;
}

}  // namespace

UTEST(viamd_task_system, a_step_job_runs_to_its_end_by_poll_and_is_listed_meanwhile) {
    ensure_task_system();

    Counter c;
    c.steps_left = c.total = 4;
    task_system::Job job;
    start_counter(&job, &c, MD_STR_LIT("Stepped"));
    EXPECT_EQ(1, c.steps);              // its first step at once
    EXPECT_TRUE(job.running());
    const uint64_t id = job.activity_id();
    float last = -1.0f;
    for (int i = 0; i < 2; ++i) {
        float p = 0.0f;
        EXPECT_TRUE(listed(id, &p));
        EXPECT_GT(p, last);
        last = p;
        task_system::poll(1.0);
        EXPECT_TRUE(job.running());
    }
    task_system::poll(1.0);
    EXPECT_TRUE(job.status() == task_system::Status::Completed);
    EXPECT_NEAR(1.0f, job.progress(), 1.0e-6f);
    EXPECT_FALSE(listed(id));
    EXPECT_EQ(1, c.ends);               // ended once, done: its results are the owner's
    EXPECT_EQ(4, c.steps);
    task_system::poll(1.0);             // no longer stepped
    EXPECT_EQ(4, c.steps);
}

UTEST(viamd_task_system, a_step_job_done_at_once_ends_at_once) {
    ensure_task_system();

    Counter c;
    c.steps_left = c.total = 1;
    task_system::Job job;
    start_counter(&job, &c, MD_STR_LIT("At once"));
    EXPECT_TRUE(job.status() == task_system::Status::Completed);
    EXPECT_EQ(1, c.ends);
}

UTEST(viamd_task_system, a_step_job_cancelled_from_the_overlay) {
    ensure_task_system();

    Counter c;
    c.steps_left = c.total = 100;
    task_system::Job job;
    start_counter(&job, &c, MD_STR_LIT("Overlay"));
    const uint64_t id = job.activity_id();
    ASSERT_TRUE(listed(id));
    task_system::activity_cancel(id);
    task_system::poll(1.0);
    EXPECT_TRUE(job.status() == task_system::Status::Cancelled);
    EXPECT_FALSE(listed(id));
    EXPECT_EQ(1, c.ends);
    EXPECT_EQ(1, c.steps);              // not stepped once cancelled
}

UTEST(viamd_task_system, a_step_job_dropped_by_its_owner_ends_at_once) {
    ensure_task_system();

    Counter c;
    c.steps_left = c.total = 100;
    task_system::Job job;
    start_counter(&job, &c, {});
    job.cancel();
    EXPECT_TRUE(job.status() == task_system::Status::Cancelled);
    EXPECT_EQ(1, c.ends);

    Counter d;
    d.steps_left = d.total = 100;
    {
        task_system::Job scoped;
        start_counter(&scoped, &d, {});
    }   // the destructor cancels
    EXPECT_EQ(1, d.ends);
    task_system::poll(1.0);             // nothing left that points at it
    EXPECT_EQ(1, d.steps);
}

UTEST(viamd_task_system, a_step_job_that_fails) {
    ensure_task_system();

    Counter c;
    c.fail = true;
    task_system::Job job;
    start_counter(&job, &c, MD_STR_LIT("Fails"));
    EXPECT_TRUE(job.status() == task_system::Status::Failed);
    EXPECT_FALSE(listed(job.activity_id()));
    EXPECT_EQ(1, c.ends);
}

UTEST(viamd_task_system, wait_runs_a_step_job_to_its_end) {
    ensure_task_system();

    Counter c;
    c.steps_left = c.total = 100;
    task_system::Job job;
    start_counter(&job, &c, {});
    job.wait();
    EXPECT_TRUE(job.status() == task_system::Status::Completed);
    EXPECT_EQ(2, c.steps);
    EXPECT_EQ(1, c.ends);
}

UTEST(viamd_task_system, a_step_jobs_on_main_runs_in_the_main_queue) {
    ensure_task_system();

    Counter c;
    c.steps_left = c.total = 2;
    int calls = 0;
    task_system::Status got = task_system::Status::Idle;
    task_system::Job job;
    start_counter(&job, &c, {}, [&](task_system::Status st) { calls += 1; got = st; });
    task_system::poll(1.0);
    EXPECT_TRUE(job.running());         // until its follow-up has run
    EXPECT_EQ(0, calls);
    task_system::run_main_queue();
    EXPECT_EQ(1, calls);
    EXPECT_TRUE(got == task_system::Status::Completed);
    EXPECT_TRUE(job.status() == task_system::Status::Completed);
}

UTEST(viamd_task_system, cancel_all_leaves_step_jobs_to_their_owners) {
    ensure_task_system();

    Counter c;
    c.steps_left = c.total = 100;
    task_system::Job job;
    start_counter(&job, &c, {});
    task_system::cancel_all();
    task_system::poll(1.0);
    EXPECT_TRUE(job.running());
    EXPECT_EQ(2, c.steps);
}

UTEST(viamd_task_system, progress_flag_is_the_jobs_progress) {
    ensure_task_system();

    EXPECT_TRUE(task_system::progress_flag() == nullptr);     // outside of a job
    task_system::Job job;
    job.start({}, []() {
        volatile float* p = task_system::progress_flag();
        if (p) *p = 0.5f;
    });
    job.wait();
    EXPECT_NEAR(0.5f, job.progress(), 1.0e-6f);
}
