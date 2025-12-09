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

### 🐧 Linux (Recommended for Dev)

**1. Install Prerequisites**
```bash
sudo apt update
sudo apt install build-essential cmake libsdl2-dev libsdl2-ttf-dev libfftw3-dev
```

**2. Build & Run**
Run the provided script to configure, build, and launch:
```bash
chmod +x build_linux.sh
./build_linux.sh

# Or run manually later:
./source/build/WASAbiApp
```

### 🪟 Windows (Option A: MSYS2 / MinGW)
*Best for command-line users who want a Linux-like experience on Windows.*

**1. Install Prerequisites**
1. Download and install [MSYS2](https://www.msys2.org/).
2. Open the **MSYS2 MinGW x64** terminal (Blue Icon).
3. Install dependencies:
```bash
pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake mingw-w64-x86_64-make mingw-w64-x86_64-SDL2 mingw-w64-x86_64-SDL2_ttf mingw-w64-x86_64-fftw
```

**2. Build & Run**
Run this script from the **MSYS2 MinGW x64** terminal:
```bash
./build_windows.sh

# Or run manually later:
./source/build/WASAbiApp.exe
```

### 💜 Windows (Option B: Visual Studio 2022)
*Best for users who prefer a full GUI IDE.*

**1. Install Prerequisites (vcpkg)**
Visual Studio needs **vcpkg** to manage the C++ libraries automatically.
1. Install [vcpkg](https://github.com/microsoft/vcpkg).
2. Install the libraries:
```cmd
vcpkg install sdl2 sdl2-ttf fftw3 --triplet=x64-windows
vcpkg integrate install
```

**2. Build & Run**
1. Open Visual Studio 2022.
2. Select **File > Open > Folder...** and select the project folder.
3. Visual Studio will detect `CMakeLists.txt` and configure automatically.
4. Select **WASAbiApp.exe** from the startup item dropdown (green arrow).
5. Press **F5** to build and run.

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