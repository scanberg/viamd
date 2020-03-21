# viamd
Visual Interactive Analysis of Molecular Dynamics

## Status
![build_and_test](https://github.com/scanberg/viamd/workflows/build_and_test/badge.svg?branch=master)

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
- libgtk3-dev
- pkgconf

Example for Ubuntu:
#### sudo apt-get install libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libgtk-3-dev pkgconf

### Step 2: Configure using CMAKE

https://cmake.org/

### Step 3: Build!

### Step 4: Profit?

## Binaries

https://github.com/scanberg/viamd/releases/download/v0.1a/viamd_0.1a_osx.zip

https://github.com/scanberg/viamd/releases/download/v0.1a/viamd_0.1a_win64.zip

https://github.com/scanberg/viamd/releases/download/v0.1a/viamd_0.1a_linux.zip

## Acknowledgements

https://github.com/ocornut/imgui

https://github.com/glfw/glfw

https://github.com/skaslev/gl3w

https://github.com/codeplea/tinyexpr
