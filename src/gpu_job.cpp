#include <gpu_job.h>

#include <core/md_common.h>
#include <core/md_os.h>

// The jobs under way, oldest first
static GpuJob* g_head = nullptr;
static GpuJob* g_tail = nullptr;

// A job in its first frames is given this much of the budget: one that is done within it ends at once, rather
// than a frame later for each step it takes
#define GPU_JOB_YOUNG_MS      50.0
#define GPU_JOB_YOUNG_BUDGET  4.0
#define GPU_JOB_BUDGET        1.0

GpuJob::~GpuJob() {
    cancel();
}

void GpuJob::start(const Ops& ops, void* job, md_str_t label, bool cancellable) {
    cancel();
    m_ops = ops;
    m_progress = 0.0f;
    m_t0 = (int64_t)md_tick_now();
    if (!job) {
        m_status = Status::Failed;
        return;
    }
    m_job = job;
    m_status = Status::Running;
    m_prev = g_tail;
    m_next = nullptr;
    if (g_tail) g_tail->m_next = this;
    else        g_head = this;
    g_tail = this;
    if (!md_str_empty(label)) {
        m_activity.cancellable = cancellable;
        m_activity.cancel.store(0, std::memory_order_relaxed);
        m_activity.progress.store(0.0f, std::memory_order_relaxed);
        task_system::activity_begin(&m_activity, label);
        m_listed = true;
    }
}

float GpuJob::progress() const {
    return m_status == Status::Done ? 1.0f : m_progress;
}

double GpuJob::seconds() const {
    return m_t0 ? md_tick_to_seconds(md_tick_now() - (md_tick_t)m_t0) : 0.0;
}

// No longer under way: out of the list and the overlay
void GpuJob::finish(Status status) {
    if (m_prev) m_prev->m_next = m_next;
    else if (g_head == this) g_head = m_next;
    if (m_next) m_next->m_prev = m_prev;
    else if (g_tail == this) g_tail = m_prev;
    m_prev = m_next = nullptr;
    if (m_listed) {
        task_system::activity_end(&m_activity);
        m_listed = false;
    }
    m_status = status;
}

void GpuJob::cancel() {
    if (m_status != Status::Running && !m_job) return;
    if (m_job && m_ops.end) m_ops.end(m_job);   // its pieces in flight complete, nothing more is issued
    m_job = nullptr;
    finish(m_status == Status::Running ? Status::Cancelled : m_status);
}

void* GpuJob::take() {
    if (m_status != Status::Done || !m_job) return nullptr;
    void* job = m_job;
    m_job = nullptr;
    m_status = Status::Idle;
    return job;
}

void GpuJob::advance(double budget_ms) {
    if (m_status != Status::Running) return;
    if (m_listed && m_activity.cancel.load(std::memory_order_relaxed)) {
        cancel();
        return;
    }
    const Status st = m_ops.poll(m_job, budget_ms);
    m_progress = m_ops.progress ? m_ops.progress(m_job) : 0.0f;
    if (m_listed) m_activity.progress.store(m_progress, std::memory_order_relaxed);
    if (st == Status::Running) return;
    if (st != Status::Done || !m_ops.owner_takes) {
        // its results are in the owner's buffers: it has done its part
        if (m_ops.end) m_ops.end(m_job);
        m_job = nullptr;
    }
    finish(st);
}

void gpu_jobs_poll(double budget_ms) {
    const md_tick_t t0 = md_tick_now();
    for (GpuJob* j = g_head; j;) {
        GpuJob* next = j->m_next;     // it may finish, and leave the list
        const double spent = md_tick_to_milliseconds(md_tick_now() - t0);
        const double left = MD_MAX(budget_ms - spent, 0.0);
        const double want = j->seconds() * 1.0e3 < GPU_JOB_YOUNG_MS ? GPU_JOB_YOUNG_BUDGET : GPU_JOB_BUDGET;
        j->advance(MD_MIN(want, left));
        j = next;
    }
}

void gpu_jobs_cancel_all() {
    while (g_head) g_head->cancel();
}
