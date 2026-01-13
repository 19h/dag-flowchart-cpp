#include "placement/BlockPlacer.h"
#include "graph/GraphTraversal.h"
#include "core/LinkedListPool.h"

#include <algorithm>
#include <cassert>
#include <climits>

// ============================================================================
// Public Interface
// ============================================================================

void BlockPlacer::place(const std::vector<uint64_t>& blockOrder, LayoutState& state) const {
    // Step 1: Assign each block to a row based on graph depth
    assignRows(state, blockOrder);

    // Step 2: Extract spanning tree edges from the DAG
    selectTree(state);

    // Step 3: Identify merge points for better alignment of branches
    findMergePoints(state);

    // Step 4: Compute column positions using subtree placement
    placeSubtrees(blockOrder, state);
}

// ============================================================================
// Merge Point Detection
// ============================================================================

void BlockPlacer::findMergePoints(LayoutState& state) const {
    for (auto& [blockId, block] : state.grid_blocks) {
        // Skip blocks without multiple children
        if (block.tree_edge.size() < 2) continue;

        // Find the potential merge point (grandchild of this block)
        GridBlock* mergeBlock = nullptr;
        int grandChildCount = 0;

        for (const auto childId : block.tree_edge) {
            const auto& childBlock = state.grid_blocks[childId];
            if (!childBlock.tree_edge.empty()) {
                mergeBlock = &state.grid_blocks[childBlock.tree_edge[0]];
            }
            grandChildCount += childBlock.tree_edge.size();
        }

        // Valid merge point: exactly one grandchild reachable via tree edges
        if (!mergeBlock || grandChildCount != 1) continue;

        // Count how many children connect to the merge point via DAG edges
        int blocksGoingToMerge = 0;
        int blockWithTreeEdge = 0;

        for (const auto childId : block.tree_edge) {
            const auto& childBlock = state.grid_blocks[childId];

            // Check if this child has a non-tree edge to the merge point
            const bool connectsToMerge = std::any_of(
                childBlock.dag_edge.begin(),
                childBlock.dag_edge.end(),
                [&](uint64_t target) { return target == mergeBlock->id; }
            );

            if (connectsToMerge) {
                // Track which child has the tree edge (vs DAG edge) to merge
                if (childBlock.tree_edge.size() == 1) {
                    blockWithTreeEdge = blocksGoingToMerge;
                }
                blocksGoingToMerge++;
            } else {
                break;
            }
        }

        // Record merge point and adjust column positioning
        if (blocksGoingToMerge > 0) {
            block.mergeBlock = mergeBlock->id;
            // Offset the tree-edge child to center the merge structure
            state.grid_blocks[block.tree_edge[blockWithTreeEdge]].col =
                blockWithTreeEdge * 2 - (blocksGoingToMerge - 1);
        }
    }
}

// ============================================================================
// Subtree Placement Algorithm
// ============================================================================

void BlockPlacer::placeSubtrees(const std::vector<uint64_t>& blockOrder, LayoutState& state) const {
    // Pool for managing silhouette linked lists efficiently
    // Two lists per node: left boundary and right boundary
    LinkedListPool<int> sides(blockOrder.size() * 2);

    // Process blocks bottom-up (leaves first, then parents)
    for (const uint64_t blockId : blockOrder) {
        auto& block = state.grid_blocks[blockId];

        // --------------------------------------------------------------------
        // Case 1: Leaf Node
        // Initialize with unit-width silhouette
        // --------------------------------------------------------------------
        if (block.tree_edge.empty()) {
            block.row_count = 1;
            block.col = 0;
            block.lastRowRight = 2;       // Right edge at column 2
            block.lastRowLeft = 0;        // Left edge at column 0
            block.leftPosition = 0;       // Leftmost extent
            block.rightPosition = 2;      // Rightmost extent
            block.leftSideShape = sides.makeList(0);   // Left silhouette
            block.rightSideShape = sides.makeList(2);  // Right silhouette
            continue;
        }

        // --------------------------------------------------------------------
        // Case 2: Internal Node
        // Merge child subtrees and position this node above them
        // --------------------------------------------------------------------

        // Start with the first child's properties
        auto& firstChild = state.grid_blocks[block.tree_edge[0]];
        auto leftSide = firstChild.leftSideShape;
        auto rightSide = firstChild.rightSideShape;

        block.row_count = firstChild.row_count;
        block.lastRowRight = firstChild.lastRowRight;
        block.lastRowLeft = firstChild.lastRowLeft;
        block.leftPosition = firstChild.leftPosition;
        block.rightPosition = firstChild.rightPosition;

        // Merge remaining children one by one
        for (size_t i = 1; i < block.tree_edge.size(); i++) {
            auto& child = state.grid_blocks[block.tree_edge[i]];

            // Walk both silhouettes to find minimum horizontal offset
            int minOffset = INT_MIN;
            int leftPos = 0, rightPos = 0;
            auto leftIt = sides.head(rightSide);        // Current subtree's right edge
            auto rightIt = sides.head(child.leftSideShape);  // New child's left edge
            int maxLeftWidth = 0;
            int minRightPos = child.col;

            while (leftIt && rightIt) {
                leftPos += *leftIt;
                rightPos += *rightIt;
                minOffset = std::max(minOffset, leftPos - rightPos);
                maxLeftWidth = std::max(maxLeftWidth, leftPos);
                minRightPos = std::min(minRightPos, rightPos);
                ++leftIt;
                ++rightIt;
            }

            // Calculate final offset based on placement strategy
            int rightTreeOffset;
            if (tightPlacement_) {
                // Pack subtrees as close as possible
                rightTreeOffset = minOffset;
            } else {
                // Use bounding box for the shorter subtree
                rightTreeOffset = leftIt
                    ? maxLeftWidth - child.leftPosition
                    : block.rightPosition - minRightPos;
            }

            // Apply offset to child and merge silhouettes
            child.col += rightTreeOffset;

            if (leftIt) {
                // Current subtree is taller - extend with child's right edge
                *leftIt -= (rightTreeOffset + child.lastRowRight - leftPos);
                rightSide = sides.append(child.rightSideShape, sides.splitTail(rightSide, leftIt));
            } else if (rightIt) {
                // Child subtree is taller - extend left edge with child's left edge
                *rightIt += (rightPos + rightTreeOffset - block.lastRowLeft);
                leftSide = sides.append(leftSide, sides.splitTail(child.leftSideShape, rightIt));
                rightSide = child.rightSideShape;
                block.lastRowRight = child.lastRowRight + rightTreeOffset;
                block.lastRowLeft = child.lastRowLeft + rightTreeOffset;
            } else {
                // Same height - just take child's right edge
                rightSide = child.rightSideShape;
            }

            // Update combined silhouette
            *sides.head(rightSide) += rightTreeOffset;
            block.row_count = std::max(block.row_count, child.row_count);
            block.leftPosition = std::min(block.leftPosition, child.leftPosition + rightTreeOffset);
            block.rightPosition = std::max(block.rightPosition, rightTreeOffset + child.rightPosition);
        }

        // --------------------------------------------------------------------
        // Position Parent Node Above Children
        // --------------------------------------------------------------------
        int col;
        if (parentBetweenChildren_) {
            // Center between direct children
            col = 0;
            for (const uint64_t childId : block.tree_edge) {
                col += state.grid_blocks[childId].col;
            }
            col /= block.tree_edge.size();
        } else {
            // Center in subtree bounding box, constrained by outer children
            col = (block.rightPosition + block.leftPosition) / 2 - 1;
            col = std::max(col, state.grid_blocks[block.tree_edge.front()].col - 1);
            col = std::min(col, state.grid_blocks[block.tree_edge.back()].col + 1);
        }

        // Update block position and extend silhouettes upward
        block.col += col;
        block.row_count += 1;
        block.leftPosition = std::min(block.leftPosition, block.col);
        block.rightPosition = std::max(block.rightPosition, block.col + 2);

        *sides.head(leftSide) -= block.col;
        block.leftSideShape = sides.append(sides.makeList(block.col), leftSide);

        *sides.head(rightSide) -= block.col + 2;
        block.rightSideShape = sides.append(sides.makeList(block.col + 2), rightSide);

        // Make children's positions relative to parent
        for (const uint64_t childId : block.tree_edge) {
            state.grid_blocks[childId].col -= block.col;
        }
    }

    // ========================================================================
    // Finalization: Convert Relative to Absolute Positions
    // ========================================================================

    // Place root nodes (row 0) sequentially from left to right
    int nextEmptyColumn = 0;
    for (auto& [_, block] : state.grid_blocks) {
        if (block.row == 0) {
            int offset = -block.leftPosition;
            block.col += nextEmptyColumn + offset;
            nextEmptyColumn = block.rightPosition + offset + nextEmptyColumn;
        }
    }

    // Propagate absolute positions top-down
    for (auto it = blockOrder.rbegin(); it != blockOrder.rend(); ++it) {
        auto& block = state.grid_blocks[*it];
        assert(block.col >= 0);
        for (const uint64_t childId : block.tree_edge) {
            state.grid_blocks[childId].col += block.col;
        }
    }
}
