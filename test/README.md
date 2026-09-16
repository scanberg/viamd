# viamd tests

Unit tests for viamd's own code, using the same `utest.h` as mdlib's suite so the two behave
alike. Build them with the rest of the project and run everything with one command from the
repository root:

```
cmake .
make -j
ctest --output-on-failure
```

`ctest` picks up both this suite and mdlib's. Pass `-DVIAMD_UNITTEST=OFF` to leave this target out.

## What belongs here

Tests for the application's own units - `src/`, `src/app/`, `src/gfx/`. A unit can be tested only
if it links without an application: no ImGui, no OpenGL, no `ApplicationState`. Today that is

| Unit | Covered by |
|---|---|
| `src/serialization_utils.cpp` | `test_serialization.cpp` — the workspace text format |
| `src/event.cpp` | `test_event.cpp` — delivery, ordering, delays |
| `src/task_system.cpp` | `test_task_system.cpp` — range partitioning, dependencies |
| `src/color_utils.h` | `test_color_utils.cpp` — the colour space conversions, which are header only |
| `src/loader.cpp` | `test_loader.cpp` — the format dispatch tables and what `init` makes of a path |
| `src/gfx/camera_utils.cpp` | `test_camera_utils.cpp` — the view and projection transforms, as inverse pairs |

**Adding a test file:** drop a `.cpp` in this folder. It is globbed, so there is nothing to edit.

**Adding a unit under test:** list it in `UNIT_FILES` in `CMakeLists.txt`. That list stays explicit
on purpose: these files are the application, so which of them can be linked without one is part of
its structure and worth writing down. Components are the opposite case - see below.

## Testing a component

Component tests do **not** live in this folder, and `CMakeLists.txt` here never names a component.
Adding a component to viamd means adding its folder and nothing else; testing one works the same
way. A component declares itself testable by its layout:

```
src/components/foo/foo.cpp            the ImGui drawing and the event handling
src/components/foo/foo_core.{h,cpp}   the computation, linkable without an application
src/components/foo/tests/*.cpp        its tests
```

The `_core` suffix is the declaration of separability - the test build collects `*_core.cpp` from
every component folder, and that component's `tests/` alongside it. A component with no `_core`
unit contributes nothing, which is what keeps a component that is mid-split, or half-landed from
another branch, inert rather than a broken build: its tests simply wait for the unit.

Components as they stand are single translation units mixing computation with drawing and event
handling, so the work is the split itself. Keep the handling and the drawing in `foo.cpp`, and move
the computation into `foo_core.{h,cpp}` taking plain data and returning plain data. The point is
not the file boundary but what it forces: a core that knows about a clearance field or a histogram,
and nothing about viamd, can be checked against cases with closed-form answers rather than by
loading a system and looking at it.

A component's tests include its headers through the component root, as `<foo/foo_core.h>`, so two
components may name their headers alike without colliding.

## Note on global state

The event system and the task system are both global and neither offers a way to unregister. Tests
therefore share them: each event test uses its own event type and its handler ignores everything
else, and the task system is initialized once on first use rather than per test.
