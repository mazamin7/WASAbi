<p align="center">
  <img src="logo.png">
</p>

# WASAbi - a Wave-based Acoustic Simulator using ARD

C++ implementation of Adaptive Rectangular Decomposition (ARD) with frequency dependent atmospheric absorption in 2.5D (3D with constant height).

Theory:
> Gerardo Cicalese, Gabriele Ciaramella, Ilario Mazzieri; Addressing atmospheric absorption in adaptive rectangular decomposition. J. Acoust. Soc. Am. 1 October 2024; 156 (4): 2328–2339. [https://doi.org/10.1121/10.0030468](https://doi.org/10.1121/10.0030468)

Extended from ARD-simulator by [@jinnsjj](https://github.com/jinnsjj).
> [https://github.com/jinnsjj/ARD-simulator](https://github.com/jinnsjj/ARD-simulator)

## Input data
`assets/*.txt` records the structure of room on x-y plane. Note that this simulator only supports 2.5D room geometries: z should always be 0 and depth of all partition should be equal.

Input example:

partition:
```
0 0 0 3 3 3  <- partition 0: x, y, z, width, height, depth
3 0 0 3 3 3  <- partition 1: x, y, z, width, height, depth
```

source:
```
1 1 1 <- source 0: x, y, z
```

recorder:
```
1 1 1 <- recorder 0: x, y, z
```

All the values above are in real world scale (meter).

**Don't forget to add an extra blank line at the end of file.**

## Features

- Frequency dependent atmospheric absorption
- Partial absorbing boundaries through PML partitions
- Cross-platform support (Windows & Linux)
- Test cases included

## Building and Running

This project uses **CMake** and provides automated build scripts for Linux and Windows (MinGW), as well as native support for Visual Studio.

### 🟢 Prerequisites
Since this simulator is fully executed on the GPU, you will need:
- An NVIDIA GPU with CUDA compute capability 6.1 or higher (GTX 10-series or newer).
- [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads) (Tested with v12.4).
- Microsoft Visual Studio 2022 (with MSVC Build Tools).
- CMake (bundled with Visual Studio or standalone).

### 🪟 Windows Build
The project uses `build_cuda.bat` to automate the configuration and compilation via `nvcc` and `cl.exe`.

**1. Install Dependencies**
You will need SDL2, SDL2_ttf, and FreeType. Extract them into the project root as expected by `CMakeLists.txt` or configure your library paths manually.

**2. Compile from Source**
Run the automated build script from a standard shell:
```cmd
.\build_cuda.bat
```
This script will:
- Clean any previous `source/build/` directory.
- Initialize the MSVC 64-bit developer environment.
- Call `cmake` pointing to the NVIDIA CUDA compiler (`nvcc`).
- Build the project using `Ninja` with `-arch=sm_61` flags.

**3. Run the Simulator**
Upon a successful build, the executable and all required DLLs/assets will be deployed in the build directory. Run it directly:
```cmd
cd source/build
.\WASAbiApp.exe
```

## Examples

Scene 1:

partition:
```
0 0 0 5 5 5
```
![scene-1.gif](https://i.loli.net/2019/01/25/5c4b06204451f.gif)

Scene 2:

partition:
```
0 0 0 2 2 2
1 2 0 1 1 2
2 1 0 1 3 2
3 2 0 1 2 2
```
![scene-2.gif](https://i.loli.net/2019/01/25/5c4b06215ce95.gif)

Scene 2:

partition:
```
0 0 0 3 3 2
0 3 0 2 1 2
3 0 0 1 2 2
4 0 0 1 1 2
0 4 0 1 1 2
```
![scene-3.gif](https://i.loli.net/2019/01/25/5c4b0622c3267.gif)
