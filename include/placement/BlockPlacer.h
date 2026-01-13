#ifndef BLOCK_PLACER_H
#define BLOCK_PLACER_H

#include "core/GridBlock.h"
#include <vector>

/**
 * @class BlockPlacer
 * @brief Computes horizontal positions for blocks in a layered graph layout.
 *
 * This class implements a subtree-based placement algorithm that:
 * 1. Assigns row numbers to blocks based on their depth in the graph
 * 2. Extracts a spanning tree from the DAG structure
 * 3. Identifies merge points where multiple branches converge
 * 4. Places subtrees side-by-side while tracking their silhouettes
 *
 * The algorithm uses linked lists to efficiently track the left and right
 * boundaries (silhouettes) of each subtree, allowing O(n) placement where
 * n is the number of nodes.
 *
 * @par Configuration Options:
 * - tightPlacement: When true, subtrees are placed as close as possible.
 *   When false, bounding boxes are used for shorter subtrees.
 * - parentBetweenChildren: When true, parent is centered between direct
 *   children. When false, parent is centered in subtree bounding box.
 *
 * @par Usage:
 * @code
 *   BlockPlacer placer(tightPlacement, parentBetweenChildren);
 *   auto order = topoSort(state, entryId);
 *   placer.place(order, state);
 * @endcode
 */
class BlockPlacer {
public:
    /**
     * @brief Constructs a BlockPlacer with the specified placement strategy.
     * @param tightPlacement If true, minimize horizontal space between subtrees
     * @param parentBetweenChildren If true, center parent between direct children
     */
    BlockPlacer(bool tightPlacement, bool parentBetweenChildren)
        : tightPlacement_(tightPlacement)
        , parentBetweenChildren_(parentBetweenChildren) {}

    /**
     * @brief Computes column positions for all blocks.
     *
     * This is the main entry point that orchestrates the placement pipeline:
     * 1. assignRows() - Assigns vertical layer to each block
     * 2. selectTree() - Extracts spanning tree edges from DAG
     * 3. findMergePoints() - Identifies diamond-pattern merge points
     * 4. placeSubtrees() - Computes horizontal positions bottom-up
     *
     * @param blockOrder Topologically sorted block IDs (leaves first)
     * @param state Layout state containing grid_blocks to be positioned
     *
     * @post All blocks in state.grid_blocks have valid col values
     */
    void place(const std::vector<uint64_t>& blockOrder, LayoutState& state) const;

private:
    bool tightPlacement_;
    bool parentBetweenChildren_;

    /**
     * @brief Identifies merge points in diamond-shaped control flow patterns.
     *
     * A merge point occurs when multiple paths from a single parent converge
     * to the same grandchild. This information is used to improve alignment
     * of if-then-else style structures.
     *
     * @par Algorithm:
     * For each block with multiple children:
     * 1. Check if all children have edges to a common grandchild
     * 2. If so, mark that grandchild as a merge point
     * 3. Adjust initial column offset for better centering
     *
     * @param state Layout state with tree edges populated
     */
    void findMergePoints(LayoutState& state) const;

    /**
     * @brief Places subtrees using a bottom-up silhouette-tracking algorithm.
     *
     * Processes blocks in reverse topological order (leaves to roots).
     * For each block:
     * 1. If leaf: Initialize with unit width and simple silhouette
     * 2. If internal: Merge child subtrees side-by-side, tracking silhouettes
     *
     * The silhouette of a subtree is represented as a linked list of relative
     * column offsets, one per row. This allows efficient merging of subtrees
     * without recomputing absolute positions.
     *
     * @par Silhouette Representation:
     * Each row stores the horizontal distance from the previous row.
     * Example: [0, 2, -1] means "start at 0, move right 2, move left 1"
     *
     * @param blockOrder Blocks in reverse topological order
     * @param state Layout state to populate with column positions
     */
    void placeSubtrees(const std::vector<uint64_t>& blockOrder, LayoutState& state) const;
};

#endif // BLOCK_PLACER_H
