#include "layout/GraphGridLayout.h"
#include "graph/GraphTraversal.h"
#include "placement/BlockPlacer.h"
#include "routing/EdgeRouter.h"
#include "layout/CoordinateConverter.h"
#include "optimization/LayoutOptimizer.h"

#include <algorithm>

// ============================================================================
// Constructor
// ============================================================================

/**
 * Initializes the layout with predefined settings based on layout type.
 *
 * - Narrow: Tight packing, no parent centering
 * - Medium: Balanced, parent centered between children
 * - Wide: Spread out, no LP optimization
 */
GraphGridLayout::GraphGridLayout(LayoutType layoutType)
    : GraphLayout({})
{
    switch (layoutType) {
    case LayoutType::Narrow:
        tightSubtreePlacement = true;
        parentBetweenDirectChild = false;
        useLayoutOptimization = true;
        break;
    case LayoutType::Medium:
        tightSubtreePlacement = false;
        parentBetweenDirectChild = true;
        useLayoutOptimization = true;
        break;
    case LayoutType::Wide:
        tightSubtreePlacement = false;
        parentBetweenDirectChild = true;
        useLayoutOptimization = false;
        break;
    }
}

// ============================================================================
// Main Layout Algorithm
// ============================================================================

/**
 * Computes the complete layout for a graph.
 *
 * The algorithm proceeds in four phases:
 * 1. Block Placement: Assign grid positions using tree-based algorithm
 * 2. Edge Routing: Compute orthogonal edge paths
 * 3. Coordinate Conversion: Transform grid to pixel coordinates
 * 4. Optimization (optional): Minimize edge lengths using LP
 */
void GraphGridLayout::CalculateLayout(Graph& blocks, uint64_t entry,
                                      int& width, int& height) const {
    if (blocks.empty()) return;

    // ========================================================================
    // Initialize Layout State
    // ========================================================================

    LayoutState state;
    state.blocks = &blocks;

    // Validate entry point, default to first block if invalid
    if (blocks.find(entry) == blocks.end()) {
        entry = blocks.begin()->first;
    }

    // Create grid block entries for all blocks
    for (const auto& [id, _] : blocks) {
        GridBlock block;
        block.id = id;
        state.grid_blocks[id] = block;
    }

    // ========================================================================
    // Phase 1: Block Placement
    // ========================================================================
    // Uses tree-based algorithm to assign row/column to each block

    BlockPlacer placer(tightSubtreePlacement, parentBetweenDirectChild);
    auto blockOrder = topoSort(state, entry);
    placer.place(blockOrder, state);

    // ========================================================================
    // Initialize Edge Information
    // ========================================================================

    // Set up edge storage and arrow directions
    for (auto& [blockId, block] : blocks) {
        state.edge[blockId].resize(block.edges.size());
        for (size_t i = 0; i < block.edges.size(); i++) {
            state.edge[blockId][i].dest = block.edges[i].target;
            block.edges[i].arrow = GraphEdge::Down;  // Default to downward arrows
        }
    }

    // Count incoming/outgoing edges for each block (used by edge router)
    for (const auto& [blockId, edges] : state.edge) {
        auto& startBlock = state.grid_blocks[blockId];
        startBlock.outputCount = edges.size();
        for (const auto& edge : edges) {
            state.grid_blocks[edge.dest].inputCount++;
        }
    }

    // ========================================================================
    // Calculate Grid Dimensions
    // ========================================================================

    state.columns = 1;
    state.rows = 1;
    for (const auto& [_, node] : state.grid_blocks) {
        state.rows = std::max(state.rows, size_t(node.row) + 1);
        state.columns = std::max(state.columns, size_t(node.col) + 2);
    }

    // Initialize row heights and column widths from block dimensions
    state.rowHeight.assign(state.rows, 0);
    state.columnWidth.assign(state.columns, 0);
    for (const auto& [id, node] : state.grid_blocks) {
        const auto& block = blocks[id];
        state.rowHeight[node.row] = std::max(block.height, state.rowHeight[node.row]);
        // Each block spans two half-columns
        state.columnWidth[node.col] = std::max(block.width / 2, state.columnWidth[node.col]);
        state.columnWidth[node.col + 1] = std::max(block.width / 2, state.columnWidth[node.col + 1]);
    }

    // ========================================================================
    // Phase 2: Edge Routing
    // ========================================================================
    // Computes orthogonal paths for all edges, assigns offsets to prevent overlaps

    EdgeRouter router(layoutConfig);
    router.route(state);

    // ========================================================================
    // Phase 3: Coordinate Conversion
    // ========================================================================
    // Transforms grid coordinates to pixel coordinates

    CoordinateConverter converter(layoutConfig, verticalBlockAlignmentMiddle);
    converter.convert(state, width, height);

    // ========================================================================
    // Phase 4: Layout Optimization (Optional)
    // ========================================================================
    // Uses linear programming to minimize edge lengths while maintaining constraints

    if (useLayoutOptimization) {
        LayoutOptimizer optimizer(layoutConfig);
        optimizer.optimize(state);
        converter.cropToContent(blocks, width, height);
    }
}
