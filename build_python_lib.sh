#!/bin/bash
# build_python_lib.sh - Build shared library for Python interface

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Building VERSE Python interface..."
echo ""

# Check if verse.c exists
if [ ! -f c/verse.c ]; then
    echo "Error: verse.c not found in c/ subdirectory"
    exit 1
fi

# Detect OS and build
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Detected macOS, building libverse.dylib..."
    cc -O3 -fPIC -dynamiclib -o python/libverse.dylib c/verse.c -lm
    
    if [ $? -eq 0 ] && [ -f python/libverse.dylib ]; then
        echo "Successfully built python/libverse.dylib"
        ls -lh python/libverse.dylib
    else
        echo "Failed to build libverse.dylib"
        exit 1
    fi
    
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    echo "Detected Linux, building libverse.so..."
    cc -O3 -fPIC -shared -o python/libverse.so c/verse.c -lm
    
    if [ $? -eq 0 ] && [ -f python/libverse.so ]; then
        echo "Successfully built python/libverse.so"
        ls -lh python/libverse.so
    else
        echo "Failed to build libverse.so"
        exit 1
    fi
    
elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    echo "Detected Windows, building verse.dll..."
    cc -O3 -shared -o python/verse.dll c/verse.c -lm
    
    if [ $? -eq 0 ] && [ -f python/verse.dll ]; then
        echo "Successfully built python/verse.dll"
        ls -lh python/verse.dll
    else
        echo "Failed to build verse.dll"
        exit 1
    fi
    
else
    echo "Unsupported OS: $OSTYPE"
    exit 1
fi

echo ""
echo "=========================================="
echo "VERSE library built successfully!"
echo ""
echo "To install the Python package, run:"
echo "  pip install -e ."
echo "=========================================="
