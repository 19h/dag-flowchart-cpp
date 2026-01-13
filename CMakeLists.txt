cmake_minimum_required(VERSION 3.16)
project(nodite CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

set(SOURCES
    main.cpp
    src/layout/GraphGridLayout.cpp
    src/layout/CoordinateConverter.cpp
    src/placement/BlockPlacer.cpp
    src/graph/GraphTraversal.cpp
    src/routing/SegmentOffset.cpp
    src/routing/EdgeRouter.cpp
    src/optimization/LinearProgramming.cpp
    src/optimization/LayoutOptimizer.cpp
    src/io/DotParser.cpp
    src/io/SvgRenderer.cpp
)

add_executable(graph ${SOURCES})

target_compile_options(graph PRIVATE
    $<$<CONFIG:Release>:-O3 -ffast-math>
    $<$<CONFIG:Debug>:-g -O0>
)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(graph PRIVATE -Wall -Wextra)
endif()
