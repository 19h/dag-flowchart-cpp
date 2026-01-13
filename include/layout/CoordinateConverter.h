#ifndef COORDINATE_CONVERTER_H
#define COORDINATE_CONVERTER_H

#include "core/GridBlock.h"
#include <vector>

/**
 * @class CoordinateConverter
 * @brief Converts grid coordinates to pixel coordinates for final layout rendering.
 *
 * This class is responsible for the final phase of layout computation, transforming
 * the abstract grid-based representation into concrete pixel coordinates suitable
 * for rendering. It handles:
 *
 * 1. **Offset Calculation**: Computing cumulative offsets for interleaved columns
 *    (block columns and edge columns) and rows (block rows and edge rows).
 *
 * 2. **Block Positioning**: Converting grid cell positions to pixel (x, y) coordinates,
 *    centering blocks within their cells.
 *
 * 3. **Edge Polyline Conversion**: Converting grid-based edge paths to pixel coordinates,
 *    applying calculated offsets for each segment.
 *
 * 4. **Content Cropping**: Computing a bounding box around all content and shifting
 *    coordinates to minimize canvas size with appropriate margins.
 *
 * @par Coordinate System:
 * The grid uses interleaved columns for blocks and edges:
 * @verbatim
 *   | EdgeCol 0 | BlockCol 0 | EdgeCol 1 | BlockCol 1 | EdgeCol 2 |
 *   |-----------|------------|-----------|------------|-----------|
 *   | margin    | block      | edges     | block      | margin    |
 * @endverbatim
 *
 * Each edge column contains vertical edge segments routed between blocks.
 * Similarly, rows are interleaved with edge rows between block rows.
 *
 * @par Usage:
 * @code
 *   CoordinateConverter converter(config, verticalAlignMiddle);
 *   converter.adjustColumnWidths(state);  // Ensure columns fit blocks
 *   int width, height;
 *   converter.convert(state, width, height);  // Main conversion
 *   converter.cropToContent(graph, width, height);  // Minimize canvas
 * @endcode
 *
 * @see EdgeRouter for edge path computation
 * @see LayoutOptimizer for row/column width calculation
 */
class CoordinateConverter {
public:
    /**
     * @brief Constructs a CoordinateConverter with layout configuration.
     * @param config Layout configuration containing spacing parameters
     * @param verticalAlignMiddle If true, center blocks vertically in their rows;
     *        if false, align to top of row
     */
    CoordinateConverter(const GraphLayout::LayoutConfig& config, bool verticalAlignMiddle)
        : config_(config)
        , verticalAlignMiddle_(verticalAlignMiddle) {}

    /**
     * @brief Main entry point: converts all grid coordinates to pixels.
     *
     * Performs three steps:
     * 1. calculateOffsets() for columns and rows
     * 2. positionBlocks() to set block (x, y) coordinates
     * 3. convertEdgePolylines() to create edge path coordinates
     *
     * @param[in,out] state Layout state with grid positions; outputs pixel coordinates
     * @param[out] width Total width of the layout in pixels
     * @param[out] height Total height of the layout in pixels
     *
     * @pre state must have columnWidth, edgeColumnWidth, rowHeight, edgeRowHeight set
     * @post All blocks in state.blocks have valid x, y coordinates
     * @post All edges have populated polyline vectors
     */
    void convert(LayoutState& state, int& width, int& height) const;

    /**
     * @brief Computes bounding box and shifts coordinates to minimize canvas.
     *
     * Finds the minimum and maximum x/y coordinates across all blocks and edges,
     * then shifts all coordinates so the content starts near the origin with
     * appropriate margins (based on config spacing).
     *
     * @param[in,out] graph Graph with positioned blocks and edges; coordinates are shifted
     * @param[out] width Final canvas width after cropping
     * @param[out] height Final canvas height after cropping
     */
    void cropToContent(GraphLayout::Graph& graph, int& width, int& height) const;
    
    /**
     * @brief Connects edge endpoints to their source/target block boundaries.
     *
     * Sets the Y coordinate of edge start points to the bottom of the source block,
     * and end points to the top of the target block. This ensures edges visually
     * connect to blocks even if block heights vary.
     *
     * @param graph Graph with positioned blocks and edges
     *
     * @note This is called after main conversion but before cropping
     */
    static void connectEdgeEnds(GraphLayout::Graph& graph);

    /**
     * @brief Calculates cumulative offsets for interleaved column/edge-column structure.
     *
     * Given widths for block columns and edge columns, computes the starting
     * pixel offset for each column. The interleaving pattern is:
     *   EdgeCol[0], Col[0], EdgeCol[1], Col[1], ..., EdgeCol[N]
     *
     * @param columnWidth Widths of block columns
     * @param edgeColumnWidth Widths of edge columns (should be columnWidth.size() + 1)
     * @param[out] columnOffset Starting offset for each block column
     * @param[out] edgeColumnOffset Starting offset for each edge column
     * @return Total width (sum of all column and edge column widths)
     *
     * @par Example:
     * If columnWidth = [100, 80] and edgeColumnWidth = [20, 30, 20]:
     *   edgeColumnOffset = [0, 120, 230]
     *   columnOffset = [20, 150]
     *   returns 250
     */
    static int calculateOffsets(const std::vector<int>& columnWidth,
                               std::vector<int>& edgeColumnWidth,
                               std::vector<int>& columnOffset,
                               std::vector<int>& edgeColumnOffset);

    /**
     * @brief Adjusts column widths to ensure blocks fit within their cells.
     *
     * For each block, ensures that the column widths on either side are at least
     * large enough to accommodate half the block width (minus the central edge column).
     * Also updates row heights to fit the tallest block in each row.
     *
     * @param[in,out] state Layout state with grid_blocks; updates columnWidth, rowHeight
     */
    void adjustColumnWidths(LayoutState& state) const;

private:
    const GraphLayout::LayoutConfig& config_;
    bool verticalAlignMiddle_;

    /**
     * @brief Converts grid block positions to pixel coordinates.
     *
     * For each block:
     * - X = center of its edge column, minus half block width
     * - Y = top of its row (optionally centered vertically)
     *
     * @param state Layout state with offsets computed
     */
    void positionBlocks(LayoutState& state) const;

    /**
     * @brief Converts grid edge paths to pixel coordinate polylines.
     *
     * For each edge, converts its sequence of grid points into a polyline
     * of pixel coordinates. Odd-indexed points contribute X coordinates,
     * even-indexed points contribute Y coordinates.
     *
     * @param state Layout state with positioned blocks and routed edges
     */
    void convertEdgePolylines(LayoutState& state) const;
};

#endif // COORDINATE_CONVERTER_H
