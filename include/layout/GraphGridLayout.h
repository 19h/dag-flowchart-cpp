#ifndef GRAPH_GRID_LAYOUT_H
#define GRAPH_GRID_LAYOUT_H

#include "core/GraphLayout.h"
#include "core/GridBlock.h"

#include <cstdint>

/**
 * @class GraphGridLayout
 * @brief Main orchestrator for hierarchical graph layout using a grid-based approach.
 *
 * This class implements a layered (Sugiyama-style) graph layout algorithm optimized
 * for control flow graphs and similar hierarchical structures. The algorithm proceeds
 * in four main phases:
 *
 * 1. **Block Placement** (BlockPlacer): Assigns row and column positions to blocks
 *    using a tree-based subtree placement algorithm.
 *
 * 2. **Edge Routing** (EdgeRouter): Computes orthogonal edge paths between blocks,
 *    assigning each edge a "main column" and calculating segment offsets to prevent
 *    overlapping.
 *
 * 3. **Coordinate Conversion** (CoordinateConverter): Transforms abstract grid
 *    positions to concrete pixel coordinates, handling column/row width calculation
 *    and edge polyline generation.
 *
 * 4. **Layout Optimization** (LayoutOptimizer): Optional post-processing that uses
 *    linear programming to minimize edge lengths while maintaining constraints.
 *
 * @par Layout Types:
 * - **Narrow**: Tight placement, suitable for compact displays
 * - **Medium**: Balanced spacing, good for general use
 * - **Wide**: More spread out, no LP optimization for maximum clarity
 *
 * @par Configuration Options:
 * - tightSubtreePlacement: Pack subtrees more closely
 * - parentBetweenDirectChild: Center parent between immediate children
 * - verticalBlockAlignmentMiddle: Center blocks vertically in rows
 * - useLayoutOptimization: Enable LP-based optimization
 *
 * @par Usage:
 * @code
 *   GraphGridLayout layout(GraphGridLayout::LayoutType::Medium);
 *   GraphLayout::Graph blocks = ...;  // Populate blocks and edges
 *   int width, height;
 *   layout.CalculateLayout(blocks, entryBlockId, width, height);
 *   // blocks now contain x, y coordinates and edge polylines
 * @endcode
 *
 * @see BlockPlacer for tree-based block placement
 * @see EdgeRouter for orthogonal edge routing
 * @see CoordinateConverter for grid-to-pixel conversion
 * @see LayoutOptimizer for LP-based optimization
 */
class GraphGridLayout : public GraphLayout {
public:
    /**
     * @brief Predefined layout configurations.
     */
    enum class LayoutType {
        Medium,   ///< Balanced spacing, parent centered between children
        Wide,     ///< More spread out, no LP optimization
        Narrow,   ///< Tighter packing, tight subtree placement
    };

    /**
     * @brief Constructs a GraphGridLayout with the specified layout type.
     * @param layoutType Predefined configuration (defaults to Medium)
     */
    explicit GraphGridLayout(LayoutType layoutType = LayoutType::Medium);

    /**
     * @brief Computes the layout for a graph.
     *
     * This is the main entry point that orchestrates all layout phases.
     *
     * @param[in,out] blocks Map of block IDs to GraphBlock structures.
     *        Input: width, height, edges for each block.
     *        Output: x, y coordinates and edge polylines.
     * @param entry The ID of the entry block (starting point for layout)
     * @param[out] width Total width of the layout canvas
     * @param[out] height Total height of the layout canvas
     *
     * @pre blocks must contain at least one block
     * @pre Each block must have valid width and height
     * @pre Edge targets must reference valid block IDs
     * @post All blocks have valid x, y positions
     * @post All edges have populated polyline vectors
     */
    void CalculateLayout(Graph& blocks, uint64_t entry, 
                        int& width, int& height) const override;

    /**
     * @brief Enables/disables tight subtree placement.
     *
     * When enabled, subtrees are placed as close as possible using
     * silhouette-based collision detection. When disabled, bounding
     * boxes are used for more uniform spacing.
     *
     * @param enabled True to enable tight placement
     */
    void setTightSubtreePlacement(bool enabled) { tightSubtreePlacement = enabled; }

    /**
     * @brief Enables/disables parent centering between direct children.
     *
     * When enabled, parent blocks are centered between their immediate
     * children rather than the full subtree bounding box.
     *
     * @param enabled True to center between direct children
     */
    void setParentBetweenDirectChild(bool enabled) { parentBetweenDirectChild = enabled; }

    /**
     * @brief Enables/disables vertical centering of blocks in rows.
     *
     * When enabled, blocks are centered vertically within their row.
     * When disabled, blocks are aligned to the top of their row.
     *
     * @param enabled True to center blocks vertically
     */
    void setVerticalBlockAlignmentMiddle(bool enabled) { verticalBlockAlignmentMiddle = enabled; }

    /**
     * @brief Enables/disables LP-based layout optimization.
     *
     * When enabled, a linear programming solver is used to minimize
     * edge lengths while maintaining spacing constraints. This produces
     * cleaner layouts but takes additional computation time.
     *
     * @param enabled True to enable optimization
     */
    void setLayoutOptimization(bool enabled) { useLayoutOptimization = enabled; }

private:
    bool tightSubtreePlacement = false;
    bool parentBetweenDirectChild = false;
    bool verticalBlockAlignmentMiddle = false;
    bool useLayoutOptimization = true;
};

#endif // GRAPH_GRID_LAYOUT_H
