#include "utest.h"

#include <gpu_job.h>
#include <task_system.h>

#include <string.h>

// GpuJob (gpu_job.h) over a job of its own kind: a counter that is done after a number of polls. What the
// template promises its owners: advanced by gpu_jobs_poll alone, listed in the overlay while it runs, cancellable
// from there, ended once (done or not), and a job whose results are its own handed back by take().

namespace {

struct FakeJob {
    int polls_left = 0;
    int polls = 0;
    int total = 1;
    bool fail = false;
    int* ended = nullptr;     // counts the ends
};

GpuJob::Status fake_poll(void* p, double) {
    FakeJob* j = (FakeJob*)p;
    j->polls += 1;
    if (j->fail) return GpuJob::Status::Failed;
    if (j->polls_left > 0) j->polls_left -= 1;
    return j->polls_left == 0 ? GpuJob::Status::Done : GpuJob::Status::Running;
}

float fake_progress(const void* p) {
    const FakeJob* j = (const FakeJob*)p;
    return 1.0f - (float)j->polls_left / (float)j->total;
}

void fake_end(void* p) {
    FakeJob* j = (FakeJob*)p;
    if (j->ended) *j->ended += 1;
}

const GpuJob::Ops fake_ops       = { fake_poll, fake_progress, fake_end, false };
const GpuJob::Ops fake_ops_taken = { fake_poll, fake_progress, fake_end, true };

void ensure_task_system() {
    static bool initialized = false;
    if (!initialized) {
        task_system::initialize(2);
        initialized = true;
        atexit([]() { task_system::shutdown(); });
    }
}

bool listed(const char* label, float* progress = nullptr, uint64_t* id = nullptr) {
    task_system::ActivityInfo info[64];
    const size_t n = task_system::activity_list(info, 64);
    for (size_t i = 0; i < n; ++i) {
        if (strcmp(info[i].label, label) == 0) {
            if (progress) *progress = info[i].progress;
            if (id) *id = info[i].id;
            return true;
        }
    }
    return false;
}

}  // namespace

UTEST(viamd_gpu_job, runs_to_done_by_the_frame_poll_and_is_listed_meanwhile) {
    ensure_task_system();
    int ended = 0;
    FakeJob fj;
    fj.polls_left = fj.total = 4;
    fj.ended = &ended;
    GpuJob job;
    job.start(fake_ops, &fj, MD_STR_LIT("fake job a"));
    EXPECT_TRUE(job.running());
    EXPECT_TRUE(listed("fake job a"));
    float last = -1.0f;
    for (int i = 0; i < 3; ++i) {
        gpu_jobs_poll(1.0);
        EXPECT_TRUE(job.running());
        float p = 0.0f;
        EXPECT_TRUE(listed("fake job a", &p));
        EXPECT_GT(p, last);
        last = p;
    }
    gpu_jobs_poll(1.0);
    EXPECT_EQ((int)GpuJob::Status::Done, (int)job.status());
    EXPECT_EQ(1.0f, job.progress());
    EXPECT_FALSE(listed("fake job a"));
    EXPECT_EQ(1, ended);            // ended once, when done: its results are the owner's
    EXPECT_EQ(4, fj.polls);
    gpu_jobs_poll(1.0);             // no longer polled
    EXPECT_EQ(4, fj.polls);
    EXPECT_TRUE(job.take() == nullptr);     // not a job whose results are its own
}

UTEST(viamd_gpu_job, cancelled_from_the_overlay) {
    ensure_task_system();
    int ended = 0;
    FakeJob fj;
    fj.polls_left = fj.total = 100;
    fj.ended = &ended;
    GpuJob job;
    job.start(fake_ops, &fj, MD_STR_LIT("fake job b"));
    gpu_jobs_poll(1.0);
    uint64_t id = 0;
    ASSERT_TRUE(listed("fake job b", nullptr, &id));
    task_system::activity_cancel(id);
    gpu_jobs_poll(1.0);
    EXPECT_EQ((int)GpuJob::Status::Cancelled, (int)job.status());
    EXPECT_FALSE(listed("fake job b"));
    EXPECT_EQ(1, ended);
}

UTEST(viamd_gpu_job, dropped_by_its_owner_and_failed) {
    ensure_task_system();
    int ended = 0;
    {
        FakeJob fj;
        fj.polls_left = fj.total = 100;
        fj.ended = &ended;
        GpuJob job;
        job.start(fake_ops, &fj, {});       // not listed
        gpu_jobs_poll(1.0);
        EXPECT_TRUE(job.running());
    }   // the destructor cancels
    EXPECT_EQ(1, ended);
    gpu_jobs_poll(1.0);                     // nothing left that points at it

    FakeJob ff;
    ff.fail = true;
    ff.ended = &ended;
    GpuJob job;
    job.start(fake_ops, &ff, MD_STR_LIT("fake job c"));
    job.advance(1.0);                       // at once, by its owner
    EXPECT_EQ((int)GpuJob::Status::Failed, (int)job.status());
    EXPECT_FALSE(listed("fake job c"));
    EXPECT_EQ(2, ended);

    GpuJob none;
    none.start(fake_ops, nullptr, MD_STR_LIT("fake job d"));
    EXPECT_EQ((int)GpuJob::Status::Failed, (int)none.status());
    EXPECT_FALSE(listed("fake job d"));
}

UTEST(viamd_gpu_job, a_job_whose_results_are_its_own_is_taken) {
    ensure_task_system();
    int ended = 0;
    FakeJob fj;
    fj.polls_left = fj.total = 2;
    fj.ended = &ended;
    GpuJob job;
    job.start(fake_ops_taken, &fj, MD_STR_LIT("fake job e"));
    gpu_jobs_poll(1.0);
    gpu_jobs_poll(1.0);
    EXPECT_EQ((int)GpuJob::Status::Done, (int)job.status());
    EXPECT_EQ(0, ended);            // the owner ends it
    EXPECT_TRUE(job.take() == &fj);
    EXPECT_EQ((int)GpuJob::Status::Idle, (int)job.status());
    EXPECT_TRUE(job.take() == nullptr);
}
