### Atomic Queue, Fixed ring buffer queue for multithreaded programming ###
#project(atomic_queue LANGUAGES CXX)

add_library(atomic_queue INTERFACE)
target_include_directories(atomic_queue INTERFACE .)