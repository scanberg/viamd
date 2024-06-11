VIAMD

## TODO ##

    [/] Refactor events into a proper event system (VIAMD | Cleanup)
        [ ] 
    [ ] Implement Integration Tests (VIAMD | Robustness)
    [/] Color representations by Property (VIAMD | Feature)
        [ ] Current implementation is crappy, only one property for all represensations seem supported? (BUG)
        [ ] Expose text-based query for input instead of fixed drop down (Enhancement)

    [ ] Support batched visualization evaluation in script (MDLIB | Enhancement)
    [/] Implement util function to find identical structures based on input structure (MDLIB | Feature)
        [ ] Implement util function to compare if two structures are equivalent (MDLIB | Feature)
    [ ] Implement util function to find maximum common supgraph of two input graphs (MDLIB | Feature)
    [ ] Revise script interface (MDLIB | Cleanup)
        [ ] Revise Property Interface
    [ ] Compress atom type into external table. (MDLIB | Cleanup)
        Referenced by uint8 inside atom fields
        The type then holds: Label, Element (If applicable), Radii
        [ ] Revise interfaces of operations in Utils to handle this change
    [ ] Implement Hashtable. (MDLIB | Feature)
    [ ] DCD reader (MDLIB | Feature)
    [ ] Parallelize Postprocessing (VIAMD / MDLIB)
        [ ] Reformulate operations to support parallel processing (invoked through VIAMD)
    [ ] Implement Hydrogen Bond Determination


### Old design specification ###

Design principles:
- Moderate C++ or (C with some C++ features)
- Example Orthodox C++ [https://gist.github.com/bkaradzic/2e39896bc7d8c34e042b]
- C++17/20 whatever version with minimal dependencies.
- Avoid Templates unless necessary (example: Array might be a candidate for template).
- Maintain interactivity!
- Maintain low compile times! -> Fast iteration

Target Platforms:
- Win32
- OSX 10.9+ (3.3 / 4.1) Ugh
- Linux

License:
- MIT? Very permissive Open source

General:
- Small Math lib (SIMD)
    - vec2
    - vec3
    - vec4
    - mat2
    - mat3
    - mat4
    - dot
    - cross
    - operators + - * /
- Array
- Custom Allocators

Application layer:
- Move all GLFW into here (Window creation)
- OpenFileDialog
- SaveFileDialog
- (EventSystem)
- Provide IMGUI with optional define to enable disable

Rendering:
- ShaderLoading
- Camera (Perspective / Orto)
- CameraController (Trackball, FPS? Convenient swap between)
- VolumeRenderer
- SSAO
- HDR
- Tonemapping
- PBR
- Implicit Rendering?
    - Surface?