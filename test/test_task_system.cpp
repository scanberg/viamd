#include "utest.h"

#include <task_system.h>

#include <core/md_os.h>

#include <algorithm>
#include <atomic>
#include <mutex>
#include <utility>
#include <vector>

// task_system owns a global thread pool with no re-entrant initialize, so it is brought up once on
// first use and left running for the life of the process.

namespace {

void ensure_task_system() {
    static bool initialized = false;
    if (!initialized) {
        task_system::initialize(4);
        initialized = true;
        atexit([]() { task_system::shutdown(); });
    }
}

}  // namespace

UTEST(viamd_task_system, range_task_visits_every_index_once) {
    ensure_task_system();

    const uint32_t N = 10000;
    std::vector<std::atomic<int>> visits(N);
    for (uint32_t i = 0; i < N; ++i) visits[i] = 0;

    task_system::ID id = task_system::create_pool_task(STR_LIT("visit"), N,
        [&](uint32_t beg, uint32_t end, uint32_t) {
            for (uint32_t i = beg; i < end; ++i) visits[i] += 1;
        }, 1);
    task_system::enqueue_task(id);
    task_system::task_wait_for(id);

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
    std::mutex mutex;
    std::vector<std::pair<uint32_t, uint32_t>> ranges;

    task_system::ID id = task_system::create_pool_task(STR_LIT("partition"), N,
        [&](uint32_t beg, uint32_t end, uint32_t) {
            std::lock_guard<std::mutex> lock(mutex);
            ranges.push_back({beg, end});
        }, 64);
    task_system::enqueue_task(id);
    task_system::task_wait_for(id);

    std::sort(ranges.begin(), ranges.end());
    ASSERT_TRUE(!ranges.empty());
    EXPECT_EQ(0u, ranges.front().first);
    EXPECT_EQ(N,  ranges.back().second);
    for (size_t i = 1; i < ranges.size(); ++i) {
        EXPECT_EQ(ranges[i - 1].second, ranges[i].first);
    }
    for (size_t i = 0; i < ranges.size(); ++i) {
        EXPECT_LT(ranges[i].first, ranges[i].second);
    }
}

UTEST(viamd_task_system, a_dependency_runs_first) {
    ensure_task_system();

    std::atomic<int> stage(0);
    std::atomic<int> observed(-1);

    task_system::ID first = task_system::create_pool_task(STR_LIT("first"), [&]() {
        md_thread_sleep(20);
        stage = 1;
    });
    task_system::ID second = task_system::create_pool_task(STR_LIT("second"), [&]() {
        observed = stage.load();
    });

    task_system::set_task_dependency(second, first);
    task_system::enqueue_task(first);
    task_system::task_wait_for(second);

    EXPECT_EQ(1, observed.load());
}

UTEST(viamd_task_system, a_completed_task_reports_completion) {
    ensure_task_system();

    std::atomic<int> ran(0);
    task_system::ID id = task_system::create_pool_task(STR_LIT("done"), [&]() { ran = 1; });
    task_system::enqueue_task(id);
    task_system::task_wait_for(id);

    EXPECT_EQ(1, ran.load());
    EXPECT_FALSE(task_system::task_is_running(id));

    // An id that was never valid must be answered rather than crashed on
    EXPECT_FALSE(task_system::task_is_running(task_system::INVALID_ID));
    task_system::task_wait_for(task_system::INVALID_ID);
}
