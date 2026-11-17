#include "utest.h"

#include <task_system.h>

#include <core/md_executor.h>
#include <core/md_linalg.h>
#include <core/md_spatial_acc.h>
#include <core/md_unitcell.h>
#include <core/md_allocator.h>
#include <core/md_os.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <algorithm>
#include <atomic>
#include <vector>

/* mdlib's wide sections (md_executor.h, md_run_lanes) on the application's pool: from the main thread, from a job
 * which is cancelled while they run, and the procedures with a wide version against their fork-join one, to the bit.
 *
 * And, on demand, what the wide versions gain on this pool (enkiTS), against the fork-join ones on the same pool:
 *
 *   VIAMD_BENCH=1 viamd_test --filter=viamd_lanes.bench
 *
 * VIAMD_BENCH_THREADS sets the pool's threads when this is the first test to bring it up (default: every logical
 * core), VIAMD_BENCH_REPS the timed runs (default 7, the median is reported). Every line starts with the version
 * below, bumped with every change to what is measured or how. */

#define LANES_BENCH_VERSION "lanes-bench r1"

namespace {

void ensure_pool(uint32_t num_threads) {
    if (!task_system::executor()) {
        task_system::initialize(num_threads);
        atexit([]() { task_system::shutdown(); });
    }
}

struct Wide {
    uint32_t n = 0;
    std::vector<std::atomic<int>> hits;
    std::atomic<int> lanes{0};
    std::atomic<int> bad{0};
    uint32_t max_lanes = 0;
};

void wide_fn(md_lane_t* lane, void* arg) {
    Wide* w = (Wide*)arg;
    w->lanes += 1;
    MD_LANE_FOR(lane, w->n, 7, r) {
        for (uint32_t i = r.beg; i < r.end; ++i) w->hits[i] += 1;
        md_thread_sleep(0);
    }
    const uint32_t count = md_lane_count(lane);
    if (count == 0 || count > w->max_lanes) w->bad += 1;
    uint64_t v = md_lane_idx(lane) == 0 ? 0xC0FFEEull + count : 0;
    v = md_lane_bcast_u64(lane, v, 0);
    if (v != 0xC0FFEEull + count) w->bad += 1;
}

void sym_random(double* a, size_t n, unsigned seed) {
    unsigned s = seed;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            s = s * 1664525u + 1013904223u;
            const double v = (double)(s >> 8) / (double)(1u << 24) - 0.5;
            a[i * n + j] = v;
            a[j * n + i] = v;
        }
        a[i * n + i] += (double)i * 0.1;
    }
}

}  // namespace

UTEST(viamd_lanes, run_on_the_pool_from_the_main_thread) {
    ensure_pool(4);
    const md_executor_t* ex = task_system::executor();
    Wide w;
    w.n = 10000;
    w.hits = std::vector<std::atomic<int>>(w.n);
    w.max_lanes = md_executor_num_threads(ex);
    EXPECT_TRUE(md_run_lanes(ex, w.max_lanes, wide_fn, &w));
    EXPECT_GE(w.lanes.load(), 1);
    EXPECT_LE(w.lanes.load(), (int)w.max_lanes);
    EXPECT_EQ(0, w.bad.load());
    uint32_t wrong = 0;
    for (uint32_t i = 0; i < w.n; ++i) wrong += w.hits[i].load() != 1;
    EXPECT_EQ(0u, wrong);
}

UTEST(viamd_lanes, stop_with_their_job) {
    ensure_pool(4);
    struct Ctx {
        std::atomic<bool> running{false};
        std::atomic<int> result{-1};
        std::atomic<int> syncs_said_false{0};
        std::atomic<int> lanes{0};
        std::atomic<uint64_t> syncs_at_stop[64];
    } c;
    for (auto& s : c.syncs_at_stop) s = 0;
    task_system::Job job;
    job.start({}, [&]() {
        const bool ok = md_run_lanes(task_system::executor(), 64, [](md_lane_t* lane, void* arg) {
            Ctx* c = (Ctx*)arg;
            c->lanes += 1;
            MD_LANE_FOR(lane, 64, 1, r) { md_thread_sleep(1); }
            c->running = true;
            // Syncs until one says the job was cancelled: every lane at the same one
            uint64_t k = 0;
            while (md_lane_sync(lane)) ++k;
            c->syncs_said_false += 1;
            c->syncs_at_stop[md_lane_idx(lane)] = k;
        }, &c);
        c.result = ok ? 1 : 0;
    });
    while (!c.running) md_thread_sleep(1);
    md_thread_sleep(20);
    job.cancel();
    job.wait();
    EXPECT_EQ(0, c.result.load());
    EXPECT_EQ(c.lanes.load(), c.syncs_said_false.load());
    for (int i = 1; i < c.lanes.load(); ++i) EXPECT_EQ(c.syncs_at_stop[0].load(), c.syncs_at_stop[i].load());
}

UTEST(viamd_lanes, wide_is_bitwise_on_the_pool) {
    ensure_pool(4);
    md_executor_t fj = *task_system::executor();
    fj.fork_join = true;
    const md_executor_t* wide = task_system::executor();

    const size_t n = 900;
    std::vector<double> A(n * n), B(n * n), wa(n), wb(n);
    sym_random(A.data(), n, 3u);
    sym_random(B.data(), n, 3u);
    ASSERT_TRUE(md_linalg_sym_eigen_ex(A.data(), wa.data(), n, &fj));
    ASSERT_TRUE(md_linalg_sym_eigen_ex(B.data(), wb.data(), n, wide));
    EXPECT_EQ(0, memcmp(A.data(), B.data(), sizeof(double) * n * n));
    EXPECT_EQ(0, memcmp(wa.data(), wb.data(), sizeof(double) * n));

    const size_t N = 400000;
    std::vector<float> xyz(3 * N);
    unsigned s = 5;
    for (auto& x : xyz) { s = s * 1664525u + 1013904223u; x = (float)((double)(s >> 8) / (double)(1u << 24) * 160.0); }
    const md_unitcell_t cell = md_unitcell_from_extent(160, 160, 160);
    md_spatial_acc_t a = {}, b = {};
    a.alloc = b.alloc = md_get_heap_allocator();
    md_spatial_acc_desc_t desc = {};
    desc.xyz = xyz.data();
    desc.count = N;
    desc.cutoff = 6.0;
    desc.unitcell = &cell;
    desc.executor = &fj;
    md_spatial_acc_init(&a, &desc);
    desc.executor = wide;
    md_spatial_acc_init(&b, &desc);
    ASSERT_EQ(a.num_elems, b.num_elems);
    ASSERT_EQ(a.num_cells, b.num_cells);
    EXPECT_EQ(0, memcmp(a.elem_x, b.elem_x, a.num_elems * sizeof(float)));
    EXPECT_EQ(0, memcmp(a.elem_idx, b.elem_idx, a.num_elems * sizeof(uint32_t)));
    EXPECT_EQ(0, memcmp(a.cell_off, b.cell_off, (a.num_cells + 1) * sizeof(uint32_t)));
    md_spatial_acc_free(&a);
    md_spatial_acc_free(&b);
}

// ### The benchmark

namespace {

double median_of(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return (n & 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

uint64_t fnv(uint64_t h, const void* p, size_t n) {
    const uint8_t* b = (const uint8_t*)p;
    for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ull; }
    return h;
}

double ms_between(md_tick_t t0, md_tick_t t1) { return md_tick_to_milliseconds(t1 - t0); }

// Runs fn(ex, &hash) with the fork-join and the wide executor in turns, and prints a line
template <class F>
void bench_line(const char* what, int reps, F fn) {
    md_executor_t fj = *task_system::executor();
    fj.fork_join = true;
    md_executor_t wide = *task_system::executor();
    wide.fork_join = false;
    uint64_t h_fj = 0, h_wide = 0, h = 0;
    fn(&fj, &h_fj);         // untimed, one of each
    fn(&wide, &h_wide);
    bool same = h_fj == h_wide;
    std::vector<double> t_fj, t_wide;
    for (int r = 0; r < reps; ++r) {
        t_fj.push_back(fn(&fj, &h));
        same &= h == h_fj;
        t_wide.push_back(fn(&wide, &h));
        same &= h == h_fj;
    }
    const double mf = median_of(t_fj), mw = median_of(t_wide);
    printf("%s | %-40s | fork-join %9.3f  wide %9.3f  fj/wide %5.2f | %016llx%s\n", LANES_BENCH_VERSION, what, mf, mw,
           mf / mw, (unsigned long long)h_fj, same ? "" : " DIFFERS");
}

}  // namespace

UTEST(viamd_lanes, bench) {
    if (!getenv("VIAMD_BENCH")) return;
    md_os_sys_info_t si = {};
    md_os_sys_info_query(&si);
    const uint32_t cores = si.num_virtual_cores > 0 ? (uint32_t)si.num_virtual_cores : 4;
    ensure_pool(getenv("VIAMD_BENCH_THREADS") ? (uint32_t)std::max(1, atoi(getenv("VIAMD_BENCH_THREADS"))) : cores);
    const int reps = getenv("VIAMD_BENCH_REPS") ? std::max(1, atoi(getenv("VIAMD_BENCH_REPS"))) : 7;
    printf("%s | enkiTS pool of %u threads (%u logical cores) | median of %d, ms\n", LANES_BENCH_VERSION,
           task_system::num_threads(), cores, reps);

    // An empty phase: a parallel_for of nothing, against a sync (per phase, in microseconds)
    bench_line("empty phase (us)", reps, [](const md_executor_t* ex, uint64_t* hash) -> double {
        *hash = 0;
        const int phases = 2000;
        if (ex->fork_join) {
            const md_tick_t t0 = md_tick_now();
            for (int k = 0; k < phases; ++k) {
                md_parallel_for(ex, md_executor_num_threads(ex), 1, [](void*, uint32_t, uint32_t, uint32_t) {}, nullptr);
            }
            return ms_between(t0, md_tick_now()) * 1000.0 / phases;
        }
        struct D { double us; } d = {0};
        md_run_lanes(ex, md_executor_num_threads(ex), [](md_lane_t* lane, void* arg) {
            // Work enough for every worker to wake and join before the first sync closes the group, then the syncs
            MD_LANE_FOR(lane, 64, 1, r) {
                const md_tick_t t = md_tick_now();
                while (md_tick_to_nanoseconds(md_tick_now() - t) < 100000.0) {}
            }
            const md_tick_t t0 = md_tick_now();
            for (int k = 0; k < 2000; ++k) md_lane_sync(lane);
            if (md_lane_idx(lane) == 0) ((D*)arg)->us = ms_between(t0, md_tick_now()) * 1000.0 / 2000;
        }, &d);
        return d.us;
    });

    for (size_t n : {500, 1000, 2000}) {
        std::vector<double> M(n * n), V(n * n), w(n);
        sym_random(M.data(), n, 7u + (unsigned)n);
        char what[64];
        snprintf(what, sizeof(what), "md_linalg_sym_eigen_ex n = %zu", n);
        bench_line(what, reps, [&](const md_executor_t* ex, uint64_t* hash) -> double {
            V = M;
            const md_tick_t t0 = md_tick_now();
            const bool ok = md_linalg_sym_eigen_ex(V.data(), w.data(), n, ex);
            const double ms = ms_between(t0, md_tick_now());
            *hash = ok ? fnv(fnv(1469598103934665603ull, w.data(), n * sizeof(double)), V.data(), n * n * sizeof(double)) : 0;
            return ms;
        });
    }

    for (size_t N : {100000, 1000000, 4000000}) {
        const double ext = cbrt((double)N / 0.1);
        const md_unitcell_t cell = md_unitcell_from_extent(ext, ext, ext);
        std::vector<float> xyz(3 * N);
        unsigned s = 11;
        for (auto& x : xyz) { s = s * 1664525u + 1013904223u; x = (float)((double)(s >> 8) / (double)(1u << 24) * ext); }
        char what[64];
        snprintf(what, sizeof(what), "md_spatial_acc_init %zu points", N);
        bench_line(what, reps, [&](const md_executor_t* ex, uint64_t* hash) -> double {
            md_spatial_acc_t acc = {};
            acc.alloc = md_get_heap_allocator();
            md_spatial_acc_desc_t desc = {};
            desc.xyz = xyz.data();
            desc.count = N;
            desc.cutoff = 6.0;
            desc.cell_ext = 6.0;
            desc.unitcell = &cell;
            desc.executor = ex;
            const md_tick_t t0 = md_tick_now();
            md_spatial_acc_init(&acc, &desc);
            const double ms = ms_between(t0, md_tick_now());
            uint64_t h = fnv(1469598103934665603ull, acc.elem_x, acc.num_elems * sizeof(float));
            h = fnv(h, acc.elem_idx, acc.num_elems * sizeof(uint32_t));
            *hash = fnv(h, acc.cell_off, (acc.num_cells + 1) * sizeof(uint32_t));
            md_spatial_acc_free(&acc);
            return ms;
        });
    }
}
