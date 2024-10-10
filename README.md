<p align="center">
  <img src="logo.png">
</p>

# WASAbi - a Wave-based Acoustic Simulator using ARD

C++ implementation of Adaptive Rectangular Decomposition (ARD) with frequency dependent atmospheric absorption in 2.5D (3D with constant height).

Theory:
> Gerardo Cicalese, Gabriele Ciaramella, Ilario Mazzieri; Addressing atmospheric absorption in adaptive rectangular decomposition. J. Acoust. Soc. Am. 1 October 2024; 156 (4): 2328–2339. https://doi.org/10.1121/10.0030468

Extended from ARD-simulator by [@jinnsjj](https://github.com/jinnsjj).
> https://github.com/jinnsjj/ARD-simulator

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
- Test cases included

## Building and running

Use Visual Studio to build. The solution itself is self-contained, so simply building and running in Visual Studio should work. Win32 mode may cause performance issue, please run under x64 mode.

<!-- ## Note

### FFTW installation note

> <http://www.fftw.org/install/windows.html>

- right click on the project -> properties -> C/C++ -> General -> Additional include Directories.
- right click on the project -> properties -> Linker -> General -> additional library directories.
- right click on the project -> properties -> Linker -> Input -> additional Dependencies.

### SDL installation note

> <https://www.wikihow.com/Set-Up-SDL-with-Visual-Studio-2017>

- right click on the project -> properties -> C/C++ -> General -> Additional include Directories.
- right click on the project -> properties -> Linker -> General -> additional library directories.
- right click on the project -> properties -> Linker -> Input -> additional Dependencies.

### style guide

> <https://google.github.io/styleguide/cppguide.html> -->

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
