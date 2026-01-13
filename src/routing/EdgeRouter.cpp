#include "routing/EdgeRouter.h"
#include "core/BinaryTrees.h"
#include "routing/SegmentOffset.h"
#include "layout/CoordinateConverter.h"

#include <algorithm>
#include <cassert>
#include <cmath>

// ============================================================================
// Public API
// ============================================================================

void EdgeRouter::route(LayoutState& state) const {
    calculateMainColumns(state);
    buildRoughPaths(state);
    refinePlacement(state);
}

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Calculates spacing override when block has too many edges to fit at normal spacing.
 * Returns 0 if default spacing is sufficient, otherwise returns reduced spacing.
 */
int EdgeRouter::getSpacingOverride(int blockWidth, int edgeCount) const {
    if (edgeCount == 0) return 0;
    
    int maxSpacing = blockWidth / edgeCount;
    if (maxSpacing < config_.edgeHorizontalSpacing) {
        return std::max(maxSpacing, 1);
    }
    return 0;
}

// ============================================================================
// Phase 1: Column Selection (Sweep-Line Algorithm)
// ============================================================================

/**
 * Assigns a main column to each edge using a sweep-line algorithm.
 *
 * The algorithm processes the grid row-by-row from top to bottom, tracking
 * which columns are "blocked" by blocks. For each edge, it selects the best
 * available column based on distance minimization and preference rules.
 *
 * Events are processed in order:
 * 1. Block events mark columns as blocked starting at that row
 * 2. Edge events select the best column for the edge
 *
 * Edge columns are interleaved with block columns:
 *   EdgeCol 0 | BlockCol 0 | EdgeCol 1 | BlockCol 1 | EdgeCol 2 | ...
 */
void EdgeRouter::calculateMainColumns(LayoutState& state) const {
    // ---- Event Structure ----
    // Represents either a block occupying a column or an edge needing routing
    struct Event {
        uint64_t blockId;
        size_t edgeId;
        int row;
        enum Type { Edge = 0, Block = 1 } type;
    };

    // ---- Build Event List ----
    std::vector<Event> events;
    events.reserve(state.grid_blocks.size() * 2);

    for (const auto& [blockId, gridBlock] : state.grid_blocks) {
        // Block event: marks this column as blocked from this row onward
        events.push_back({blockId, 0, gridBlock.row, Event::Block});

        const auto& inputBlock = (*state.blocks)[blockId];
        int startRow = gridBlock.row + 1;

        // Initialize edge storage for this block
        auto& gridEdges = state.edge[blockId];
        gridEdges.resize(inputBlock.edges.size());

        // Edge events: each edge needs a column assigned
        for (size_t i = 0; i < inputBlock.edges.size(); i++) {
            auto targetId = inputBlock.edges[i].target;
            gridEdges[i].dest = targetId;
            const auto& targetGridBlock = state.grid_blocks[targetId];
            int endRow = targetGridBlock.row;
            // Edge event at max(startRow, endRow) ensures we process after source
            events.push_back({blockId, i, std::max(startRow, endRow), Event::Edge});
        }
    }

    // ---- Sort Events ----
    // Process by row, with blocks before edges (so columns are marked blocked first)
    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        if (a.row != b.row) return a.row < b.row;
        return static_cast<int>(a.type) < static_cast<int>(b.type);
    });

    // ---- Sweep Line Data Structure ----
    // Tracks the row at which each column becomes blocked
    // Value at column c = row where block starts (or -1 if unblocked)
    PointSetMinTree blockedColumns(state.columns + 1, -1);

    // ---- Process Events ----
    for (const auto& event : events) {
        if (event.type == Event::Block) {
            // Mark this column as blocked starting at this row
            const auto& block = state.grid_blocks[event.blockId];
            blockedColumns.set(block.col + 1, event.row);
            continue;
        }

        // ---- Edge Event: Select Best Column ----
        const auto& block = state.grid_blocks[event.blockId];
        int column = block.col + 1;  // Edge column adjacent to source block
        auto& edge = state.edge[event.blockId][event.edgeId];
        const auto& targetBlock = state.grid_blocks[edge.dest];
        auto topRow = std::min(block.row + 1, targetBlock.row);
        auto targetColumn = targetBlock.col + 1;

        // Priority 1: Use source column if unblocked
        if (blockedColumns.valueAtPoint(column) < topRow) {
            edge.mainColumn = column;
            continue;
        }

        // Priority 2: Use target column if unblocked
        if (blockedColumns.valueAtPoint(targetColumn) < topRow) {
            edge.mainColumn = targetColumn;
            continue;
        }

        // Find nearest unblocked columns to left and right
        auto nearestLeft = blockedColumns.rightMostLessThan(column, topRow);
        auto nearestRight = blockedColumns.leftMostLessThan(column, topRow);
        assert(nearestLeft != -1 && nearestRight != -1);

        // Calculate total path distance for each option
        auto distanceLeft = column - nearestLeft + abs(targetColumn - nearestLeft);
        auto distanceRight = nearestRight - column + abs(targetColumn - nearestRight);

        // Priority 3: For upward edges, prefer adjacent columns for cleaner routing
        if (targetBlock.row < block.row) {
            if (targetColumn < column &&
                blockedColumns.valueAtPoint(column + 1) < topRow &&
                column - targetColumn <= distanceLeft + 2) {
                edge.mainColumn = column + 1;
                continue;
            }
            if (targetColumn > column &&
                blockedColumns.valueAtPoint(column - 1) < topRow &&
                targetColumn - column <= distanceRight + 2) {
                edge.mainColumn = column - 1;
                continue;
            }
        }

        // Priority 4: Choose column with shorter total path distance
        if (distanceLeft != distanceRight) {
            edge.mainColumn = (distanceLeft < distanceRight) ? nearestLeft : nearestRight;
        } else {
            // Tie-breaker: alternate based on edge index for visual balance
            edge.mainColumn = (event.edgeId < state.edge[event.blockId].size() / 2)
                ? nearestLeft : nearestRight;
        }
    }
}

// ============================================================================
// Phase 2: Path Construction
// ============================================================================

/**
 * Builds orthogonal paths for all edges based on assigned main columns.
 *
 * Each path consists of up to 6 points forming an orthogonal route:
 *   [0] Start: Below source block center
 *   [1] Vertical: In source or main column
 *   [2] Horizontal: Connect to main column (if needed)
 *   [3] Vertical: Main column segment
 *   [4] Horizontal: Connect to target column (if needed)
 *   [5] End: Above target block center
 *
 * Points are assigned a "kind" value indicating horizontal alignment:
 *   -2: Far left (wrapping around left side)
 *   -1: Left half
 *    0: Center
 *   +1: Right half
 *   +2: Far right (wrapping around right side)
 */
void EdgeRouter::buildRoughPaths(LayoutState& state) const {
    for (const auto& [blockId, block] : state.grid_blocks) {
        auto& blockEdges = state.edge[blockId];

        for (auto& edge : blockEdges) {
            const auto& start = block;
            const auto& target = state.grid_blocks[edge.dest];

            // ---- Point 0: Start below source ----
            edge.addPoint(start.row + 1, start.col + 1);

            // ---- Points 1-2: Route to main column ----
            if (edge.mainColumn != start.col + 1) {
                // Horizontal offset direction: -1 = left, +1 = right
                edge.addPoint(start.row + 1, start.col + 1,
                    edge.mainColumn < start.col + 1 ? -1 : 1);
                // Vertical direction hint for upward edges
                edge.addPoint(start.row + 1, edge.mainColumn,
                    target.row <= start.row ? -2 : 0);
            }

            // ---- Point 3: Main column vertical segment ----
            // Determine kind based on whether main column is to left/right of both endpoints
            int mainColumnKind = 0;
            if (edge.mainColumn < start.col + 1 && edge.mainColumn < target.col + 1) {
                mainColumnKind = +2;  // Main column is left of both, align right
            } else if (edge.mainColumn > start.col + 1 && edge.mainColumn > target.col + 1) {
                mainColumnKind = -2;  // Main column is right of both, align left
            } else if (edge.mainColumn == start.col + 1 && edge.mainColumn != target.col + 1) {
                mainColumnKind = edge.mainColumn < target.col + 1 ? 1 : -1;
            } else if (edge.mainColumn == target.col + 1 && edge.mainColumn != start.col + 1) {
                mainColumnKind = edge.mainColumn < start.col + 1 ? 1 : -1;
            }

            edge.addPoint(target.row, edge.mainColumn, mainColumnKind);

            // ---- Points 4-5: Route from main column to target ----
            if (target.col + 1 != edge.mainColumn) {
                // Horizontal direction hint for upward edges
                edge.addPoint(target.row, target.col + 1,
                    target.row <= start.row ? 2 : 0);
                // Final horizontal segment
                edge.addPoint(target.row, target.col + 1,
                    target.col + 1 < edge.mainColumn ? 1 : -1);
            }

            // ---- Spacing Overrides ----
            // Reduce spacing for blocks with many edges
            auto startSpacingOverride = getSpacingOverride(
                (*state.blocks)[start.id].width, start.outputCount);
            auto targetSpacingOverride = getSpacingOverride(
                (*state.blocks)[target.id].width, target.inputCount);

            edge.points.front().spacingOverride = startSpacingOverride;
            edge.points.back().spacingOverride = targetSpacingOverride;

            // Propagate tighter spacing to adjacent points
            if (edge.points.size() <= 2) {
                if (startSpacingOverride && startSpacingOverride < targetSpacingOverride) {
                    edge.points.back().spacingOverride = startSpacingOverride;
                }
            } else {
                edge.points[1].spacingOverride = startSpacingOverride;
            }

            // ---- Calculate Path Length for Priority ----
            int length = 0;
            for (size_t i = 1; i < edge.points.size(); i++) {
                length += abs(edge.points[i].row - edge.points[i - 1].row) +
                         abs(edge.points[i].col - edge.points[i - 1].col);
            }
            // Longer paths get higher priority; downward edges break ties
            edge.secondaryPriority = 2 * length + (target.row >= start.row ? 1 : 0);
        }
    }
}

// ============================================================================
// Phase 3: Offset Refinement
// ============================================================================

/**
 * Calculates pixel offsets for edge segments to prevent overlapping.
 *
 * This is done in two passes:
 * 1. Vertical segments: Calculate horizontal offsets within edge columns
 * 2. Horizontal segments: Calculate vertical offsets within edge rows
 *
 * Each pass:
 * - Extracts segments from edge points
 * - Creates obstacle regions from block boundaries
 * - Calls calculateSegmentOffsets() to assign non-overlapping offsets
 * - Updates column/row widths if more space is needed
 */
void EdgeRouter::refinePlacement(LayoutState& state) const {
    int edgeIndex = 0;

    // Helper: Create EdgeSegment from a GridPoint
    auto segmentFromPoint = [&edgeIndex](const GridPoint& point, const GridEdge& edge, 
                                         int y0, int y1, int x) {
        EdgeSegment segment;
        segment.y0 = y0;
        segment.y1 = y1;
        segment.x = x;
        segment.edgeIndex = edgeIndex++;
        segment.kind = point.kind;
        segment.spacingOverride = point.spacingOverride;
        segment.secondaryPriority = edge.secondaryPriority;
        return segment;
    };

    std::vector<EdgeSegment> segments;
    std::vector<NodeSide> rightSides, leftSides;
    std::vector<int> edgeOffsets;

    // ========================================================================
    // Pass 1: Vertical Segments (horizontal offset calculation)
    // ========================================================================

    // ---- Extract vertical segments (odd-indexed points) ----
    for (auto& edgeListIt : state.edge) {
        for (const auto& edge : edgeListIt.second) {
            for (size_t j = 1; j < edge.points.size(); j += 2) {
                segments.push_back(segmentFromPoint(edge.points[j], edge,
                    edge.points[j - 1].row * 2, edge.points[j].row * 2,
                    edge.points[j].col));
            }
        }
    }

    // ---- Create block side obstacles ----
    // Blocks create obstacles that edges must route around
    for (auto& blockIt : state.grid_blocks) {
        auto& node = blockIt.second;
        auto width = (*state.blocks)[blockIt.first].width;
        auto leftWidth = width / 2;
        auto rightWidth = width - leftWidth;
        int row = node.row * 2 + 1;

        // Left side of block affects edge column to its left
        leftSides.push_back({node.col, row, row, leftWidth});
        // Right side of block affects edge column to its right
        rightSides.push_back({node.col + 1, row, row, rightWidth});
    }

    // ---- Initialize edge column widths ----
    state.edgeColumnWidth.assign(state.columns + 1, config_.blockHorizontalSpacing);
    state.edgeColumnWidth[0] = state.edgeColumnWidth.back() = config_.edgeHorizontalSpacing;
    edgeOffsets.resize(edgeIndex);

    // ---- Calculate offsets for vertical segments ----
    calculateSegmentOffsets(segments, edgeOffsets, state.edgeColumnWidth, rightSides, leftSides,
                          state.columnWidth, 2 * state.rows + 1,
                          config_.edgeHorizontalSpacing);

    // Center edges within their allocated space
    centerEdges(edgeOffsets, state.edgeColumnWidth, segments);
    edgeIndex = 0;

    // ---- Helper: Copy segment offsets back to edge points ----
    auto copySegmentsToEdges = [&](bool col) {
        int idx = 0;
        for (auto& edgeListIt : state.edge) {
            for (auto& edge : edgeListIt.second) {
                for (size_t j = col ? 1 : 2; j < edge.points.size(); j += 2) {
                    int offset = edgeOffsets[idx++];
                    if (col) {
                        // Clamp offset to keep edge within block bounds
                        GraphLayout::GraphBlock* block = nullptr;
                        if (j == 1) {
                            block = &(*state.blocks)[edgeListIt.first];
                        } else if (j + 1 == edge.points.size()) {
                            block = &(*state.blocks)[edge.dest];
                        }
                        if (block) {
                            int blockWidth = block->width;
                            int edgeColumnWidth = state.edgeColumnWidth[edge.points[j].col];
                            offset = std::max(-blockWidth / 2 + edgeColumnWidth / 2, offset);
                            offset = std::min(edgeColumnWidth / 2 + std::min(blockWidth, edgeColumnWidth) / 2, offset);
                        }
                    }
                    edge.points[j].offset = offset;
                }
            }
        }
    };

    // ---- Adjust column widths for wrap-around edges ----
    auto oldColumnWidths = state.columnWidth;
    CoordinateConverter converter(config_, false);
    converter.adjustColumnWidths(state);

    // Compensate offsets for edges that wrap around blocks
    for (auto& segment : segments) {
        auto& offset = edgeOffsets[segment.edgeIndex];
        if (segment.kind == -2) {
            // Edge wraps around left: adjust for expanded column
            offset -= (state.edgeColumnWidth[segment.x - 1] / 2 + state.columnWidth[segment.x - 1])
                     - oldColumnWidths[segment.x - 1];
        } else if (segment.kind == 2) {
            // Edge wraps around right: adjust for expanded column
            offset += (state.edgeColumnWidth[segment.x + 1] / 2 + state.columnWidth[segment.x])
                     - oldColumnWidths[segment.x];
        }
    }

    // Calculate final column offsets and copy to edges
    CoordinateConverter::calculateOffsets(state.columnWidth, state.edgeColumnWidth,
                                         state.columnOffset, state.edgeColumnOffset);
    copySegmentsToEdges(true);

    // ========================================================================
    // Pass 2: Horizontal Segments (vertical offset calculation)
    // ========================================================================

    segments.clear();
    leftSides.clear();
    rightSides.clear();
    edgeIndex = 0;

    // ---- Extract horizontal segments (even-indexed points > 0) ----
    for (auto& edgeListIt : state.edge) {
        for (const auto& edge : edgeListIt.second) {
            for (size_t j = 2; j < edge.points.size(); j += 2) {
                // y0, y1 are the x-coordinates of adjacent vertical segments
                int y0 = state.edgeColumnOffset[edge.points[j - 1].col] + edge.points[j - 1].offset;
                int y1 = state.edgeColumnOffset[edge.points[j + 1].col] + edge.points[j + 1].offset;
                segments.push_back(segmentFromPoint(edge.points[j], edge, y0, y1, edge.points[j].row));
            }
        }
    }

    edgeOffsets.resize(edgeIndex);

    // ---- Create block top/bottom obstacles ----
    for (auto& blockIt : state.grid_blocks) {
        auto& node = blockIt.second;
        auto blockWidth = (*state.blocks)[node.id].width;
        // Calculate block's horizontal span in pixel coordinates
        int leftSide = state.edgeColumnOffset[node.col + 1] +
                      state.edgeColumnWidth[node.col + 1] / 2 - blockWidth / 2;
        int rightSide = leftSide + blockWidth;

        int h = (*state.blocks)[blockIt.first].height;
        int topProfile = state.rowHeight[node.row];
        int bottomProfile = h;

        // Top and bottom of block create obstacles for horizontal edges
        leftSides.push_back({node.row, leftSide, rightSide, topProfile});
        rightSides.push_back({node.row, leftSide, rightSide, bottomProfile});
    }

    // ---- Initialize edge row heights ----
    state.edgeRowHeight.assign(state.rows + 1, config_.blockVerticalSpacing);
    state.edgeRowHeight[0] = state.edgeRowHeight.back() = config_.edgeVerticalSpacing;
    edgeOffsets.resize(edgeIndex);

    // ---- Calculate offsets for horizontal segments ----
    auto compressedCoordinates = compressCoordinates(segments, leftSides, rightSides);
    calculateSegmentOffsets(segments, edgeOffsets, state.edgeRowHeight, rightSides, leftSides,
                          state.rowHeight, compressedCoordinates, config_.edgeVerticalSpacing);

    copySegmentsToEdges(false);
}
