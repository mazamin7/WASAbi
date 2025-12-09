#!/bin/bash

# --- Configuration ---
SOURCE_DIR="source"
BUILD_DIR="${SOURCE_DIR}/build"
EXECUTABLE_NAME="WASAbiApp.exe"
CLEAN_BUILD=false

# --- Helper Function for Errors ---
error_exit() {
    echo "❌ Error: $1"
    exit 1
}

# --- 0. Argument Parsing ---
# Run './build_windows.sh clean' to wipe the build folder and start fresh
if [ "$1" == "clean" ]; then
    CLEAN_BUILD=true
fi

# --- 1. Navigate to Project Root ---
# Ensures the script works even if called from a different directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || error_exit "Could not determine script directory."

# Verify source folder exists
if [ ! -d "$SOURCE_DIR" ]; then
    error_exit "Source directory '$SOURCE_DIR' not found in $SCRIPT_DIR."
fi

# --- 2. Clean Build (If requested) ---
if [ "$CLEAN_BUILD" = true ]; then
    echo "🧹 Cleaning build directory..."
    rm -rf "$BUILD_DIR"
fi

# --- 3. Prepare Build Directory ---
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR" || error_exit "Could not enter build directory."

# --- 4. Configure CMake ---
echo "⚙️  Configuring CMake..."
# We explicitly ask for MinGW Makefiles to ensure consistency, 
# unless a configuration already exists.
if [ ! -f "CMakeCache.txt" ]; then
    cmake .. -G "MinGW Makefiles" || error_exit "CMake configuration failed."
else
    cmake .. || error_exit "CMake reconfiguration failed."
fi

# --- 5. Build Project ---
echo "🔨 Building Project..."
# This is the Robust Fix: 'cmake --build' works with BOTH Ninja and Makefiles.
# It automatically calls the correct underlying tool.
cmake --build . --parallel "$(nproc)" || error_exit "Compilation failed."

# --- 6. Success Message ---
echo "--------------------------------------"
echo "✅ Build Success!"
echo "To run the app:"
echo "./$BUILD_DIR/$EXECUTABLE_NAME"
echo "--------------------------------------"