#pragma once

// GPU WORK OVER FRAMES
//
// The template for anything long on the GPU in viamd: an evaluation that must not hold the display off the GPU,
// must not trip a watchdog, and shows how far it is in the overlay.
//
//   mdlib    the kernels, and a JOB of them (md_gto_iso_gpu_evaluate_job_*, md_gto_int_gpu_*_job_*, ...): the
//            work as pieces of a few ms of GPU time, each its own submission, only a little issued ahead of what
//            the GPU has finished (md_gpu_pieces.h). A job is begun by whoever knows what it computes, and owns
//            nothing of that: its inputs and outputs are the caller's.
//
//   here     GpuJob takes the mdlib job and runs it: advanced once a frame for every job at once
//            (gpu_jobs_poll, from the frame loop, with one budget for all of them), listed in the overlay while it
//            runs (a task_system activity: its label, its progress, a cancel button), and ended when it is done,
//            cancelled from the overlay, or dropped by its owner. The owner looks at status() when it likes -
//            nothing has to remember to poll it - and takes its results (its own buffers, read back) once Done.
//
// A new kind of GPU work then takes three steps: kernels that evaluate a range of their work, a job over
// md_gpu_pieces in mdlib, and a GpuJob around it in viamd. Its progress bar comes with the GpuJob.
//
// The GL thread only (the GPU stream is its own).

#include <task_system.h>
#include <core/md_str.h>

#include <stdint.h>

struct md_gto_iso_gpu_eval_job_t;
struct md_gto_iso_gpu_job_t;
struct md_gto_int_gpu_job_t;

class GpuJob {
public:
    enum class Status : uint8_t {
        Idle,           // never started
        Running,
        Done,           // its results are there (in the owner's buffers, as far as the GPU's stream has come)
        Failed,
        Cancelled,      // from the overlay, or by cancel(): its results are not to be read
    };

    // What GpuJob needs of a kind of mdlib job
    struct Ops {
        Status (*poll)(void* job, double budget_ms);    // advance; Done once every piece has passed
        float  (*progress)(const void* job);            // 0 to 1
        void   (*end)(void* job);                       // release, done or not; never waits
        bool   owner_takes;                             // Done, the job is kept for take(): its results are its own
    };

    GpuJob() = default;
    ~GpuJob();                  // cancels
    GpuJob(const GpuJob&) = delete;
    GpuJob& operator=(const GpuJob&) = delete;

    // Takes the mdlib job (NULL: Failed at once), ending one under way first. label: listed in the overlay while it
    // runs, cancellable there; empty: not listed.
    void start(const Ops& ops, void* job, md_str_t label, bool cancellable = true);
    void start(md_gto_iso_gpu_eval_job_t* job, md_str_t label);
    void start(md_gto_int_gpu_job_t* job, md_str_t label);
    // An isosurface extraction: its mesh is the job's own, taken with take() once Done (md_gto_iso_gpu_job_end)
    void start(md_gto_iso_gpu_job_t* job, md_str_t label);

    Status status() const { return m_status; }
    bool   running() const { return m_status == Status::Running; }
    float  progress() const;
    double seconds() const;     // since it began

    // Ends it where it is: what it has in flight completes, nothing more is issued
    void cancel();

    // The mdlib job, Done, given back to the caller who ends it (for a job whose results are its own, such as an
    // extraction's mesh). NULL unless Done; the GpuJob is Idle after.
    void* take();

    // Advances it now, as gpu_jobs_poll does once a frame: for an owner that has just started it, so that a job done
    // within budget_ms ends in this frame rather than the next
    void advance(double budget_ms);

    // Advances every job under way: once a frame, from the GL thread. The budget is shared; a job in its first
    // frames is given more of it, so that a short evaluation ends within them.
    friend void gpu_jobs_poll(double budget_ms);
    friend void gpu_jobs_cancel_all();

private:
    void finish(Status status);

    Ops      m_ops = {};
    void*    m_job = nullptr;
    Status   m_status = Status::Idle;
    int64_t  m_t0 = 0;          // md_tick_t
    float    m_progress = 0.0f;
    task_system::Activity m_activity;
    bool     m_listed = false;
    GpuJob*  m_prev = nullptr;  // the jobs under way
    GpuJob*  m_next = nullptr;
};

void gpu_jobs_poll(double budget_ms);

// The jobs under way, for the shutdown: every one ended (its pieces in flight complete on the GPU)
void gpu_jobs_cancel_all();
