#ifndef EDGE_ROUTER_H
#define EDGE_ROUTER_H

#include "core/GridBlock.h"
#include "core/EdgeSegment.h"

/**
 * @class EdgeRouter
 * @brief Routes edges between blocks using a column-based orthogonal routing algorithm.
 *
 * This class implements a three-phase edge routing algorithm that produces
 * orthogonal (right-angle) edge paths between source and target blocks:
 *
 * 1. **Column Selection** (calculateMainColumns): Uses a sweep-line algorithm
 *    to assign each edge a "main column" - the vertical column where most of
 *    the edge's vertical travel will occur. The algorithm processes events
 *    (blocks and edges) by row, tracking which columns are blocked.
 *
 * 2. **Path Construction** (buildRoughPaths): Constructs initial edge paths
 *    as sequences of grid points. Each edge typically has 2-6 points forming
 *    an orthogonal path: start -> horizontal to main column -> vertical ->
 *    horizontal to target -> end.
 *
 * 3. **Offset Refinement** (refinePlacement): Calculates pixel offsets for
 *    each edge segment to prevent overlapping. Uses SegmentOffset utilities
 *    to assign offsets within edge columns while avoiding node boundaries.
 *
 * @par Edge Path Structure:
 * Each edge path consists of alternating horizontal and vertical segments:
 * @verbatim
 *   Point 0: Start (below source block)
 *   Point 1: First vertical segment (in main column or source column)
 *   Point 2: Horizontal segment (connects to main column if different)
 *   Point 3: Main vertical segment
 *   Point 4: Horizontal segment (connects to target column if different)
 *   Point 5: End (above target block)
 * @endverbatim
 *
 * @par Point "kind" Values:
 * Each point has a "kind" that indicates its horizontal alignment within
 * the edge column:
 * - -2: Align to left edge of column (for edges wrapping around left)
 * - -1: Align to left half
 * -  0: Center (default)
 * - +1: Align to right half
 * - +2: Align to right edge of column (for edges wrapping around right)
 *
 * @par Usage:
 * @code
 *   EdgeRouter router(config);
 *   router.route(state);
 *   // state.edge now contains routed paths with offsets
 * @endcode
 *
 * @see GridEdge for edge path representation
 * @see EdgeSegment for segment offset calculation data
 * @see SegmentOffset for offset calculation utilities
 */
class EdgeRouter {
public:
    /**
     * @brief Constructs an EdgeRouter with the specified layout configuration.
     * @param config Layout configuration containing spacing parameters
     */
    explicit EdgeRouter(const GraphLayout::LayoutConfig& config) : config_(config) {}

    /**
     * @brief Routes all edges in the layout state.
     *
     * This is the main entry point that orchestrates the routing pipeline:
     * 1. calculateMainColumns() - Assign vertical columns to edges
     * 2. buildRoughPaths() - Construct initial orthogonal paths
     * 3. refinePlacement() - Calculate offsets to prevent overlaps
     *
     * @param state Layout state containing blocks and edges to route
     *
     * @pre state.grid_blocks must be populated with column positions
     * @post state.edge contains routed paths with point offsets
     * @post state.edgeColumnWidth contains calculated column widths
     * @post state.edgeRowHeight contains calculated row heights
     */
    void route(LayoutState& state) const;

private:
    const GraphLayout::LayoutConfig& config_;

    /**
     * @brief Assigns a main column to each edge using a sweep-line algorithm.
     *
     * The main column is the vertical column where an edge's primary vertical
     * travel occurs. The algorithm:
     *
     * 1. Creates events for all blocks (at their row) and edges (at destination row)
     * 2. Sorts events by row, with blocks processed before edges
     * 3. Sweeps through events, tracking which columns are blocked by blocks
     * 4. For each edge, selects the best unblocked column based on:
     *    - Preference for source or target column if unblocked
     *    - Distance minimization to nearest unblocked column
     *    - Tie-breaking based on edge index for visual consistency
     *
     * @par Data Structure:
     * Uses PointSetMinTree to efficiently query:
     * - Whether a column is blocked at a given row
     * - The nearest unblocked column to the left/right
     *
     * @param state Layout state with grid_blocks; populates state.edge[].mainColumn
     */
    void calculateMainColumns(LayoutState& state) const;

    /**
     * @brief Constructs initial orthogonal paths for all edges.
     *
     * For each edge from source to target:
     * 1. Start point: Just below the source block center
     * 2. If main column differs from source: Add horizontal segment
     * 3. Vertical segment in main column
     * 4. If main column differs from target: Add horizontal segment
     * 5. End point: Just above the target block center
     *
     * Also calculates:
     * - Point "kind" values for horizontal alignment within edge columns
     * - Spacing overrides for blocks with many edges (to prevent crowding)
     * - Secondary priority for edge ordering (based on path length)
     *
     * @param state Layout state with mainColumn set; populates edge points
     */
    void buildRoughPaths(LayoutState& state) const;

    /**
     * @brief Calculates pixel offsets for edge segments to prevent overlapping.
     *
     * This method processes edges in two passes:
     *
     * **Pass 1 - Vertical Segments:**
     * - Extracts all vertical segments (odd-indexed points)
     * - Creates node side obstacles from block boundaries
     * - Calls calculateSegmentOffsets() to assign horizontal offsets
     * - Adjusts column widths if segments need more space
     *
     * **Pass 2 - Horizontal Segments:**
     * - Extracts all horizontal segments (even-indexed points > 0)
     * - Creates node side obstacles from block top/bottom
     * - Calls calculateSegmentOffsets() to assign vertical offsets
     * - Adjusts row heights if segments need more space
     *
     * @param state Layout state with rough paths; populates point offsets
     *
     * @see calculateSegmentOffsets() in SegmentOffset.h
     * @see centerEdges() for final offset centering
     */
    void refinePlacement(LayoutState& state) const;

    /**
     * @brief Calculates spacing override for blocks with many edges.
     *
     * When a block has many outgoing or incoming edges, the default spacing
     * between edge attachment points may be too wide to fit within the block
     * width. This method calculates a reduced spacing.
     *
     * @param blockWidth Width of the block in pixels
     * @param edgeCount Number of edges attached to this side of the block
     * @return Reduced spacing value, or 0 if default spacing is sufficient
     */
    int getSpacingOverride(int blockWidth, int edgeCount) const;
};

#endif // EDGE_ROUTER_H
