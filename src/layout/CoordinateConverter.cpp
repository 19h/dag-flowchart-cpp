#include "layout/CoordinateConverter.h"
#include <algorithm>
#include <cassert>

// ============================================================================
// Column/Row Width Adjustment
// ============================================================================

/**
 * Ensures column widths accommodate block widths.
 *
 * Each block spans two half-columns (left and right of its center edge column).
 * This method ensures each half-column is wide enough for the block's half-width.
 */
void CoordinateConverter::adjustColumnWidths(LayoutState& state) const {
    state.rowHeight.assign(state.rows, 0);
    state.columnWidth.assign(state.columns, 0);

    for (const auto& node : state.grid_blocks) {
        const auto& block = (*state.blocks)[node.first];
        const auto& gridPos = node.second;

        // Row height = max height of blocks in that row
        state.rowHeight[gridPos.row] = std::max(block.height, state.rowHeight[gridPos.row]);

        // Block width minus center edge column, split across two half-columns
        int edgeWidth = state.edgeColumnWidth[gridPos.col + 1];
        int columnWidth = (block.width - edgeWidth) / 2;

        // Update both adjacent columns to fit this block
        state.columnWidth[gridPos.col] = std::max(columnWidth, state.columnWidth[gridPos.col]);
        state.columnWidth[gridPos.col + 1] = std::max(columnWidth, state.columnWidth[gridPos.col + 1]);
    }
}

// ============================================================================
// Offset Calculation
// ============================================================================

/**
 * Computes cumulative pixel offsets for interleaved columns.
 *
 * Layout pattern: EdgeCol[0], Col[0], EdgeCol[1], Col[1], ..., EdgeCol[N]
 *
 * Returns the total width (or height when used for rows).
 */
int CoordinateConverter::calculateOffsets(const std::vector<int>& columnWidth,
                                         std::vector<int>& edgeColumnWidth,
                                         std::vector<int>& columnOffset,
                                         std::vector<int>& edgeColumnOffset) {
    assert(edgeColumnWidth.size() == columnWidth.size() + 1);

    int position = 0;
    edgeColumnOffset.resize(edgeColumnWidth.size());
    columnOffset.resize(columnWidth.size());

    // Interleave: edge column, then block column, repeating
    for (size_t i = 0; i < columnWidth.size(); i++) {
        edgeColumnOffset[i] = position;
        position += edgeColumnWidth[i];
        columnOffset[i] = position;
        position += columnWidth[i];
    }

    // Final edge column (right margin)
    edgeColumnOffset.back() = position;
    position += edgeColumnWidth.back();
    return position;
}

// ============================================================================
// Main Conversion Entry Point
// ============================================================================

/**
 * Main conversion: transforms grid coordinates to pixel coordinates.
 *
 * Steps:
 * 1. Calculate column and row offsets from widths/heights
 * 2. Position blocks using calculated offsets
 * 3. Convert edge grid paths to pixel polylines
 * 4. Connect edge endpoints to block boundaries
 */
void CoordinateConverter::convert(LayoutState& state, int& width, int& height) const {
    // Calculate cumulative offsets for columns and rows
    width = calculateOffsets(state.columnWidth, state.edgeColumnWidth,
                            state.columnOffset, state.edgeColumnOffset);
    height = calculateOffsets(state.rowHeight, state.edgeRowHeight,
                             state.rowOffset, state.edgeRowOffset);

    positionBlocks(state);
    convertEdgePolylines(state);
    connectEdgeEnds(*state.blocks);
}

// ============================================================================
// Block Positioning
// ============================================================================

/**
 * Sets pixel (x, y) coordinates for each block.
 *
 * X position: Center of the block's edge column, minus half block width.
 *   This centers the block horizontally over its edge column.
 *
 * Y position: Top of the block's row, optionally centered vertically.
 */
void CoordinateConverter::positionBlocks(LayoutState& state) const {
    for (auto& block : (*state.blocks)) {
        const auto& gridBlock = state.grid_blocks[block.first];
        auto& blockPos = block.second;

        // Center horizontally on edge column
        blockPos.x = state.edgeColumnOffset[gridBlock.col + 1] +
                    state.edgeColumnWidth[gridBlock.col + 1] / 2 -
                    blockPos.width / 2;

        // Position vertically (top-aligned or centered)
        blockPos.y = state.rowOffset[gridBlock.row];
        if (verticalAlignMiddle_) {
            blockPos.y += (state.rowHeight[gridBlock.row] - blockPos.height) / 2;
        }
    }
}

// ============================================================================
// Edge Polyline Conversion
// ============================================================================

/**
 * Converts grid-based edge paths to pixel coordinate polylines.
 *
 * Each edge path consists of alternating vertical and horizontal segments.
 * Grid points are converted as follows:
 * - Odd-indexed points (1, 3, 5, ...): Provide X coordinates (edge column + offset)
 * - Even-indexed points (2, 4, ...): Provide Y coordinates (edge row + offset)
 *
 * The polyline is built incrementally, with each point completing the previous
 * segment and starting a new one.
 */
void CoordinateConverter::convertEdgePolylines(LayoutState& state) const {
    for (auto& [blockId, block] : (*state.blocks)) {
        for (size_t i = 0; i < block.edges.size(); i++) {
            auto& resultEdge = block.edges[i];
            const auto& routingPoints = state.edge[blockId][i];

            resultEdge.polyline.clear();
            // Start point: X will be set by first vertical segment
            resultEdge.polyline.push_back(PointF(0, block.y + block.height));

            for (size_t j = 1; j < routingPoints.points.size(); j++) {
                const auto& point = routingPoints.points[j];

                if (j & 1) {
                    // Odd index: vertical segment, provides X coordinate
                    int x = state.edgeColumnOffset[point.col] + point.offset;
                    resultEdge.polyline.back().setX(x);
                    resultEdge.polyline.push_back(PointF(x, 0));
                } else {
                    // Even index: horizontal segment, provides Y coordinate
                    int y = state.edgeRowOffset[point.row] + point.offset;
                    resultEdge.polyline.back().setY(y);
                    resultEdge.polyline.push_back(PointF(0, y));
                }
            }
        }
    }
}

// ============================================================================
// Edge Endpoint Connection
// ============================================================================

/**
 * Ensures edge endpoints touch their source/target block boundaries.
 *
 * After polyline conversion, the start and end Y coordinates may not exactly
 * match block boundaries (due to row centering, varying block heights, etc.).
 * This method sets:
 * - Start Y = bottom of source block
 * - End Y = top of target block
 */
void CoordinateConverter::connectEdgeEnds(GraphLayout::Graph& graph) {
    for (auto& [id, block] : graph) {
        for (auto& edge : block.edges) {
            const auto& target = graph[edge.target];
            edge.polyline.front().ry() = block.y + block.height;
            edge.polyline.back().ry() = target.y;
        }
    }
}

// ============================================================================
// Content Cropping
// ============================================================================

/**
 * Computes bounding box and shifts all coordinates to minimize canvas size.
 *
 * Algorithm:
 * 1. Find min/max X and Y across all blocks and edge polylines
 * 2. Expand bounds by configured margin spacing
 * 3. Shift all coordinates by -minPos to move content near origin
 * 4. Return final canvas dimensions
 */
void CoordinateConverter::cropToContent(GraphLayout::Graph& graph, int& width, int& height) const {
    // Handle empty graph
    if (graph.empty()) {
        width = std::max(1, config_.edgeHorizontalSpacing);
        height = std::max(1, config_.edgeVerticalSpacing);
        return;
    }

    // Initialize bounds from first block
    const auto& firstBlock = graph.begin()->second;
    int minPos[2] = {firstBlock.x, firstBlock.y};
    int maxPos[2] = {firstBlock.x, firstBlock.y};

    // Helper: expand bounds to include a point
    auto updateBounds = [&](int x, int y) {
        minPos[0] = std::min(minPos[0], x);
        minPos[1] = std::min(minPos[1], y);
        maxPos[0] = std::max(maxPos[0], x);
        maxPos[1] = std::max(maxPos[1], y);
    };

    // ---- Compute bounding box ----
    for (const auto& [id, block] : graph) {
        // Include block corners
        updateBounds(block.x, block.y);
        updateBounds(block.x + block.width, block.y + block.height);

        // Include all edge polyline points
        for (const auto& edge : block.edges) {
            for (const auto& point : edge.polyline) {
                updateBounds(static_cast<int>(point.x), static_cast<int>(point.y));
            }
        }
    }

    // ---- Add margin spacing ----
    minPos[0] -= config_.edgeHorizontalSpacing;
    minPos[1] -= config_.edgeVerticalSpacing;
    maxPos[0] += config_.edgeHorizontalSpacing;
    maxPos[1] += config_.edgeVerticalSpacing;

    // ---- Shift coordinates to origin ----
    for (auto& [id, block] : graph) {
        block.x -= minPos[0];
        block.y -= minPos[1];

        for (auto& edge : block.edges) {
            for (auto& point : edge.polyline) {
                point -= PointF(minPos[0], minPos[1]);
            }
        }
    }

    // ---- Calculate final dimensions ----
    width = maxPos[0] - minPos[0];
    height = maxPos[1] - minPos[1];
}
