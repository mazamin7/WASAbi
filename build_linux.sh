#!/bin/bash

# 1. Enter source directory
cd source || { echo "Source folder not found"; exit 1; }

# 2. Create build folder if it doesn't exist
mkdir -p build
cd build

# 3. Configure with CMake
echo "Configuring CMake..."
cmake ..

# 4. Build
echo "Building Project..."
make -j$(nproc)  # Uses all CPU cores for faster build

# 5. Success Message
if [ $? -eq 0 ]; then
    echo "--------------------------------------"
	echo "✅ Build Success!"
    echo "--------------------------------------"
else
    echo "Build Failed."
fi