# VIAMD
Visual Interactive Analysis of Molecular Dynamics

VIAMD is an interactive analysis tool for molecular dynamics (MD) written in C/C++. VIAMD is developed in the Department for Media and Information Technology (LIU, Norrköping). It exposes a rudementary script language that is used to declare operations which are performed over the frames of the trajectory.
The results can then be viewed in the different windows exposed in the application. 
<p align="center">
<img src="https://user-images.githubusercontent.com/38646069/225882518-85309730-47bf-4b3a-9237-cdacbd0f3ecf.png" alt="This is an overview of the viamd software" width="800"/>
</p>


## Status
[![Windows (MSVC 19)](https://github.com/scanberg/viamd/actions/workflows/windows.yml/badge.svg?branch=master)](https://github.com/scanberg/viamd/actions/workflows/windows.yml)
[![Ubuntu 20.04 (GCC 9)](https://github.com/scanberg/viamd/actions/workflows/ubuntu20.yml/badge.svg)](https://github.com/scanberg/viamd/actions/workflows/ubuntu20.yml)
[![MacOS (Clang)](https://github.com/scanberg/viamd/actions/workflows/macos.yml/badge.svg)](https://github.com/scanberg/viamd/actions/workflows/macos.yml)


## Running VIAMD 
### Windows
For windows, we recommend to use the latest binary available on the ![release page](https://github.com/scanberg/viamd/releases/).

### Ubuntu and MacOs
For Linux and MacOs, we recommend you to build VIAMD on your machine following the procedure below.

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
- libgtk-3-dev (recommended for ubuntu,use Zenity for kubuntu)
- pkgconf

Example for Ubuntu:
```sudo apt-get install libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libgtk-3-dev pkgconf```

### Step 2: Configure using CMAKE

https://cmake.org/
(Example using terminal and default build configuration for system)
```
cd viamd
mkdir build
cd build
cmake ..
```

### Step 3: Build!
```cmake --build .```

### Step 4: Run
```
cd bin
./viamd
```

## Binaries
https://github.com/scanberg/viamd/releases/

## Documentation
Documentation about VIAMD is available on the github [wiki](https://github.com/scanberg/viamd/wiki). The two first chapters relate to the [visual](https://github.com/scanberg/viamd/wiki/1.-Visual) and [analysis](https://github.com/scanberg/viamd/wiki/2.-Analysis) features respectively, where we highlight the interactive part of software. The third chapter focus on the VIAMD [language](https://github.com/scanberg/viamd/wiki/3.-Language) used for scripting and the fourth chapter propose a serie of [tutorial](https://github.com/scanberg/viamd/wiki/4.-Tutorials) (under construction). 

## Acknowledgements

https://github.com/glfw/glfw

https://github.com/dougbinks/enkiTS

https://github.com/ocornut/imgui

https://github.com/epezent/implot

https://github.com/BalazsJako/ImGuiColorTextEdit

https://github.com/skaslev/gl3w

https://github.com/max0x7ba/atomic_queue

https://github.com/mlabbe/nativefiledialog

https://github.com/nothings/stb

#
<p align="center">
<img src="https://user-images.githubusercontent.com/38646069/227518757-590d3bfa-ee6b-42eb-b371-25743da1b208.png"  width="150" alt="This is the VIAMD logo" >
</p>


