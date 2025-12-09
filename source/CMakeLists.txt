cmake_minimum_required(VERSION 3.10)
project(WASAbi_Linux)

set(CMAKE_CXX_STANDARD 17)

# ... (keep your find_package lines here) ...
find_package(SDL2 REQUIRED)
find_package(SDL2_ttf REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(OpenMP REQUIRED)
pkg_check_modules(FFTW3 REQUIRED fftw3)

include_directories(
    ${SDL2_INCLUDE_DIRS}
    ${SDL2_TTF_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIRS}
    .
    ../inih-r58      # <--- NEW: Tells CMake to look in the folder next door
)

file(GLOB SOURCES "*.cpp")

# <--- NEW: We add "../inih-r58/ini.c" so the logic gets compiled too
add_executable(WASAbiApp ${SOURCES} "../inih-r58/ini.c")

target_link_libraries(WASAbiApp
    ${SDL2_LIBRARIES}
    ${SDL2_TTF_LIBRARIES}
    SDL2_ttf
    ${FFTW3_LIBRARIES}
    OpenMP::OpenMP_CXX
    stdc++fs
    m
)