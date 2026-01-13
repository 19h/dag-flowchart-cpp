#!/bin/bash

# Build the graph layout tool

set -e

if command -v cmake &> /dev/null; then
    mkdir -p build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make -j$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
    echo "Build complete: ./build/graph"
else
    clang++ -std=c++17 -O3 -ffast-math -Iinclude -o graph \
        main.cpp \
        src/GraphGridLayout.cpp \
        src/SegmentOffset.cpp \
        src/LinearProgramming.cpp \
        src/DotParser.cpp \
        src/SvgRenderer.cpp
    echo "Build complete: ./graph"
fi
