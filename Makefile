# Makefile for nodite graph layout library

CXX := clang++
CXXFLAGS := -std=c++17 -Wall -Wextra -Iinclude

# Build type: debug or release (default: release)
BUILD ?= release

ifeq ($(BUILD),debug)
    CXXFLAGS += -g -O0 -DDEBUG
    BUILD_DIR := build/debug
else
    CXXFLAGS += -O3 -ffast-math -DNDEBUG
    BUILD_DIR := build/release
endif

# Output binary
TARGET := $(BUILD_DIR)/graph

# Source files
SOURCES := \
    main.cpp \
    src/layout/GraphGridLayout.cpp \
    src/layout/CoordinateConverter.cpp \
    src/placement/BlockPlacer.cpp \
    src/graph/GraphTraversal.cpp \
    src/routing/SegmentOffset.cpp \
    src/routing/EdgeRouter.cpp \
    src/optimization/LinearProgramming.cpp \
    src/optimization/LayoutOptimizer.cpp \
    src/io/DotParser.cpp \
    src/io/SvgRenderer.cpp

# Object files
OBJECTS := $(SOURCES:%.cpp=$(BUILD_DIR)/%.o)

# Header dependencies
DEPS := $(OBJECTS:.o=.d)

# Default target
all: $(TARGET)

# Link
$(TARGET): $(OBJECTS)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile with automatic dependency generation
$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c -o $@ $<

# Include dependency files
-include $(DEPS)

# Run the program
run: $(TARGET)
	$(TARGET)

# Run and verify output matches reference
test: $(TARGET)
	$(TARGET)
	@diff -q graph_reference.svg graph.svg && echo "Test passed: output matches reference"

# Clean build artifacts
clean:
	rm -rf build/debug build/release

# Clean all generated files
distclean: clean
	rm -f graph.svg

# Debug build shortcut
debug:
	$(MAKE) BUILD=debug

# Release build shortcut
release:
	$(MAKE) BUILD=release

# Show help
help:
	@echo "Usage: make [target] [BUILD=debug|release]"
	@echo ""
	@echo "Targets:"
	@echo "  all       Build the project (default)"
	@echo "  run       Build and run the program"
	@echo "  test      Build, run, and verify output matches reference"
	@echo "  clean     Remove build artifacts"
	@echo "  distclean Remove build artifacts and generated files"
	@echo "  debug     Build with debug flags"
	@echo "  release   Build with release flags (default)"
	@echo "  help      Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make                  # Release build"
	@echo "  make BUILD=debug      # Debug build"
	@echo "  make debug            # Debug build (shortcut)"
	@echo "  make test             # Build and verify output"

.PHONY: all run test clean distclean debug release help
