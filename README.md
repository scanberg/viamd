# viamd
Visual Interactive Analysis of Molecular Dynamics

<img src="https://github.com/scanberg/viamd/wiki/img/viamd_environment.png" alt="viamd" width="500"/>

## Status
[![build_and_test](https://github.com/scanberg/viamd/workflows/build_and_test/badge.svg?branch=master)](https://github.com/scanberg/viamd/actions)

## Building
### Step 1: Clone the repository:

(Make sure the fetch the submodules as well)

```git clone --recurse-submodules git://github.com/scanberg/viamd.git```

### (Step 1.5: Install dependencies of external libs:) (Only for linux users)
#### glfw
- libx11-dev
- libxrandr-dev
- libxinerama-dev
- libxcursor-dev

#### nativefiledialog
- libgtk-3-dev
- pkgconf

Example for Ubuntu:
#### sudo apt-get install libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libgtk-3-dev pkgconf

### Step 2: Configure using CMAKE

https://cmake.org/

### Step 3: Build!

### Step 4: Run

## Binaries
New version comming soon. The aim is to provide this for Windows x64 and Linux x64 (Tested on Ubuntu).

## Acknowledgements

https://github.com/ocornut/imgui

https://github.com/glfw/glfw

https://github.com/skaslev/gl3w

https://github.com/epezent/implot

https://github.com/BalazsJako/ImGuiColorTextEdit
