#!/bin/bash
# build_verse.sh - Build script for verse.c shared library

set -e  # Exit immediately if any command fails

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "Building VERSE library in: $SCRIPT_DIR"

# Check if verse.c exists
if [ ! -f verse.c ]; then
    echo "Error: verse.c not found in current directory"
    exit 1
fi

# Detect OS and build
if [[ "$OSTYPE" == "darwin"* ]]; then
    # macOS
    echo "Detected macOS, building libverse.dylib..."
    cc -O3 -fPIC -dynamiclib -o libverse.dylib verse.c -lm
    
    if [ $? -eq 0 ] && [ -f libverse.dylib ]; then
        echo "Successfully built libverse.dylib"
        ls -lh libverse.dylib
    else
        echo "Failed to build libverse.dylib"
        exit 1
    fi
    
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    # Linux
    echo "Detected Linux, building libverse.so..."
    cc -O3 -fPIC -shared -o libverse.so verse.c -lm
    
    if [ $? -eq 0 ] && [ -f libverse.so ]; then
        echo "Successfully built libverse.so"
        ls -lh libverse.so
    else
        echo "Failed to build libverse.so"
        exit 1
    fi
    
elif [[ "$OSTYPE" == "msys" ]] || [[ "$OSTYPE" == "cygwin" ]]; then
    # Windows
    echo "Detected Windows, building verse.dll..."
    cc -O3 -shared -o verse.dll verse.c -lm
    
    if [ $? -eq 0 ] && [ -f verse.dll ]; then
        echo "Successfully built verse.dll"
        ls -lh verse.dll
    else
        echo "Failed to build verse.dll"
        exit 1
    fi
    
else
    echo "Unsupported OS: $OSTYPE"
    exit 1
fi

echo ""
echo "VERSE library built successfully!"
echo "You can now use it with verse_python.py"