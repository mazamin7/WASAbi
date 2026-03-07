cmake_minimum_required(VERSION 3.18)
project(WASAbi_CrossPlatform LANGUAGES C CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# --- 1. Find CUDA and OpenMP ---
find_package(OpenMP REQUIRED)
find_package(CUDAToolkit REQUIRED) # Find CUDA Toolkit for cuFFT and cudart

# --- 2. Include Directories ---
set(PROJ_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/..")

include_directories(
    ${PROJ_ROOT}/SDL2-2.28.4/include
    ${PROJ_ROOT}/SDL2_ttf-2.20.2/include
    .
    ../inih-r58
)

file(GLOB ALL_SOURCES "*.cpp" "*.cu")
set(SOURCES ${ALL_SOURCES})

# --- 3. Create Executable ---
set(INI_SOURCE "D:/Documents/projects/WASAbi/inih-r58/ini.c")
add_executable(WASAbiApp ${SOURCES} ${INI_SOURCE})

# --- 4. OS-Specific Logic (Linking & Copying) ---
set(CURRENT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
set(SRC_ASSETS  "${CURRENT_DIR}/assets")
set(SRC_CONFIG  "${CURRENT_DIR}/config")
set(SRC_FONT    "${CURRENT_DIR}/font")
set(DEST_DIR    "$<TARGET_FILE_DIR:WASAbiApp>")

if(WIN32)
    # --- WINDOWS SETUP ---
    message(STATUS "Windows detected: Configuring NVCC MSVC DLLs and Assets...")

    # 1. Link Libraries Manually using included .lib files
    target_link_libraries(WASAbiApp
        ${PROJ_ROOT}/SDL2-2.28.4/lib/x64/SDL2.lib
        ${PROJ_ROOT}/SDL2_ttf-2.20.2/lib/x64/SDL2_ttf.lib
        OpenMP::OpenMP_CXX
        CUDA::cufft
        CUDA::cudart
    )
    target_link_options(WASAbiApp PRIVATE "/SUBSYSTEM:CONSOLE")

    # 2. Define DLL Paths for Copying
    set(MY_SDL2_DLL     "${PROJ_ROOT}/SDL2-2.28.4/lib/x64/SDL2.dll")
    set(MY_TTF_DLL      "${PROJ_ROOT}/SDL2_ttf-2.20.2/lib/x64/SDL2_ttf.dll")
    set(MY_FREETYPE_DLL "${PROJ_ROOT}/libfreetype-dll/x64/libfreetype-6.dll")
    # CUDA Runtime DLLs (CUDA 12.4 naming convention)
    set(MY_CUDART_DLL   "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.4/bin/cudart64_12.dll")
    set(MY_CUFFT_DLL    "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.4/bin/cufft64_11.dll")

    # 3. Copy Everything (DLLs + Assets)
    add_custom_command(TARGET WASAbiApp POST_BUILD
        # SDL/font DLLs
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MY_SDL2_DLL}"     "${DEST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MY_TTF_DLL}"      "${DEST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MY_FREETYPE_DLL}" "${DEST_DIR}"
        # CUDA DLLs
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MY_CUDART_DLL}"   "${DEST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MY_CUFFT_DLL}"    "${DEST_DIR}"
        # Assets (Disabled copying per user request - use source originals)
        # COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_ASSETS}" "${DEST_DIR}/assets"
        # COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_CONFIG}" "${DEST_DIR}/config"
        # COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_FONT}"   "${DEST_DIR}/font"
        COMMENT "Deploying Windows DLLs and Assets..."
    )
else()
    # LINUX SETUP omitted for brevity since we are compiling on Windows
endif()