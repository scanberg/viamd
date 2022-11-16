# viamd
Visual Interactive Analysis of Molecular Dynamics

<img src="https://github.com/scanberg/viamd/wiki/img/viamd_environment.png" alt="viamd" width="500"/>

## Status
[![Windows (MSVC 19)](https://github.com/scanberg/viamd/actions/workflows/windows.yml/badge.svg?branch=master)](https://github.com/scanberg/viamd/actions/workflows/windows.yml)
[![Ubuntu 20.04 (GCC 9)](https://github.com/scanberg/viamd/actions/workflows/ubuntu20.yml/badge.svg)](https://github.com/scanberg/viamd/actions/workflows/ubuntu20.yml)
[![MacOS (Clang)](https://github.com/scanberg/viamd/actions/workflows/macos.yml/badge.svg)](https://github.com/scanberg/viamd/actions/workflows/macos.yml)



## Building
### Step 1: Clone the repository:

(Make sure the fetch the submodules as well)

```git clone --recurse-submodules https://github.com/scanberg/viamd.git```

### (Step 1.5: Install dependencies of external libs:) (Only for linux users)
#### glfw
- libx11-dev
- libxrandr-dev
- libxinerama-dev
- libxcursor-dev

#### nativefiledialog
- libgtk-3-dev (optional, fallback is Zenity)
- pkgconf

Example for Ubuntu:
```sudo apt-get install libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev pkgconf```

### Step 2: Configure using CMAKE

https://cmake.org/

### Step 3: Build!
- mkdir build
- cd build/
- cmake ../
- make

### Step 4: Run
- ./viamd 

## Binaries
https://github.com/scanberg/viamd/releases/

## Acknowledgements

https://github.com/ocornut/imgui

https://github.com/BalazsJako/ImGuiColorTextEdit

https://github.com/epezent/implot

https://github.com/glfw/glfw

https://github.com/skaslev/gl3w

https://github.com/dougbinks/enkiTS

https://github.com/max0x7ba/atomic_queue

https://github.com/mlabbe/nativefiledialog

https://github.com/nothings/stb


