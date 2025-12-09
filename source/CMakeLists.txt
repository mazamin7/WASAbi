cmake_minimum_required(VERSION 3.10)
project(WASAbi_CrossPlatform)

set(CMAKE_CXX_STANDARD 17)

# --- 1. Find Packages ---
find_package(PkgConfig REQUIRED)
pkg_check_modules(SDL2 REQUIRED sdl2)
pkg_check_modules(SDL2_TTF REQUIRED SDL2_ttf)
pkg_check_modules(FFTW3 REQUIRED fftw3)
find_package(OpenMP REQUIRED)

# --- 2. Include Directories ---
include_directories(
    ${SDL2_INCLUDE_DIRS}
    ${SDL2_TTF_INCLUDE_DIRS}
    ${FFTW3_INCLUDE_DIRS}
    .
    ../inih-r58
)

file(GLOB SOURCES "*.cpp")

# --- 3. Create Executable ---
add_executable(WASAbiApp ${SOURCES} "../inih-r58/ini.c")

# --- 4. OS-Specific Logic (Linking & Copying) ---

# Define Asset Paths (Same for both OSs)
set(CURRENT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(SRC_ASSETS  "${CURRENT_DIR}/assets")
set(SRC_CONFIG  "${CURRENT_DIR}/config")
set(SRC_FONT    "${CURRENT_DIR}/font")
set(DEST_DIR    "$<TARGET_FILE_DIR:WASAbiApp>")

if(WIN32)
    # --- WINDOWS SETUP ---
    message(STATUS "Windows detected: Configuring DLLs and Assets...")

    # 1. Link Libraries
    target_link_libraries(WASAbiApp
        mingw32
        SDL2main
        SDL2
        SDL2_ttf
        ${FFTW3_LIBRARIES}
        OpenMP::OpenMP_CXX
        stdc++fs
        m
    )

    # 2. Define DLL Paths
    set(PROJ_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")
    set(MY_SDL2_DLL     "${PROJ_ROOT}/SDL2-2.28.4/lib/x64/SDL2.dll")
    set(MY_TTF_DLL      "${PROJ_ROOT}/SDL2_ttf-2.20.2/lib/x64/SDL2_ttf.dll")
    set(MY_FFTW_DLL     "${PROJ_ROOT}/fftw-3.3.5-dll/x64/libfftw3-3.dll")
    set(MY_FREETYPE_DLL "${PROJ_ROOT}/libfreetype-dll/x64/libfreetype-6.dll")

    get_filename_component(COMPILER_PATH ${CMAKE_CXX_COMPILER} DIRECTORY)
    set(SYS_STDCPP_DLL  "${COMPILER_PATH}/libstdc++-6.dll")
    set(SYS_GCC_DLL     "${COMPILER_PATH}/libgcc_s_seh-1.dll")
    set(SYS_PTHREAD_DLL "${COMPILER_PATH}/libwinpthread-1.dll")

    # 3. Copy Everything (DLLs + Assets)
    add_custom_command(TARGET WASAbiApp POST_BUILD
        # DLLs
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MY_SDL2_DLL}"     "${DEST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MY_TTF_DLL}"      "${DEST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MY_FFTW_DLL}"     "${DEST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MY_FREETYPE_DLL}" "${DEST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SYS_STDCPP_DLL}"  "${DEST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SYS_GCC_DLL}"     "${DEST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${SYS_PTHREAD_DLL}" "${DEST_DIR}"
        # Assets
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_ASSETS}" "${DEST_DIR}/assets"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_CONFIG}" "${DEST_DIR}/config"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_FONT}"   "${DEST_DIR}/font"
        COMMENT "Deploying Windows DLLs and Assets..."
    )

else()
    # --- LINUX SETUP ---
    message(STATUS "Linux detected: Configuring Assets...")

    # 1. Link Libraries
    target_link_libraries(WASAbiApp
        ${SDL2_LIBRARIES}
        ${SDL2_TTF_LIBRARIES}
        ${FFTW3_LIBRARIES}
        OpenMP::OpenMP_CXX
        stdc++fs
        m
    )

    # 2. Copy Assets Only (No DLLs needed on Linux)
    add_custom_command(TARGET WASAbiApp POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_ASSETS}" "${DEST_DIR}/assets"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_CONFIG}" "${DEST_DIR}/config"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_FONT}"   "${DEST_DIR}/font"
        COMMENT "Deploying Linux Assets..."
    )
endif()