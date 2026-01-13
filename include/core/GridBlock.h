#ifndef GRID_BLOCK_H
#define GRID_BLOCK_H

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "core/GraphLayout.h"
#include "core/LinkedListPool.h"

/**
 * @brief Internal grid-based representation of a graph block.
 * 
 * Used during layout calculation to track grid positions and tree structure.
 */
struct GridBlock {
    uint64_t id;                          ///< Block identifier
    std::vector<uint64_t> tree_edge;      ///< Outgoing tree edges (spanning tree)
    std::vector<uint64_t> dag_edge;       ///< Outgoing DAG edges (after cycle removal)
    std::size_t has_parent = false;       ///< Whether this node has a parent in the tree
    int inputCount = 0;                   ///< Number of incoming edges
    int outputCount = 0;                  ///< Number of outgoing edges

    int row_count = 0;                    ///< Height of subtree in rows
    int col = 0;                          ///< Column position in grid
    int row = 0;                          ///< Row position in grid

    uint64_t mergeBlock = 0;              ///< ID of merge point block (if any)

    int lastRowLeft;                      ///< Left side of last row of subtree
    int lastRowRight;                     ///< Right side of last row of subtree
    int leftPosition;                     ///< Leftmost position of subtree
    int rightPosition;                    ///< Rightmost position of subtree
    LinkedListPool<int>::List leftSideShape;   ///< Shape of left side of subtree
    LinkedListPool<int>::List rightSideShape;  ///< Shape of right side of subtree
};

/**
 * @brief A point in the edge routing grid.
 */
struct GridPoint {
    int row;
    int col;
    int offset;
    int16_t kind;              ///< Type of point for routing decisions
    int16_t spacingOverride;   ///< Override default spacing (0 = use default)
};

/**
 * @brief Edge representation in the grid coordinate system.
 */
struct GridEdge {
    uint64_t dest;             ///< Destination block ID
    int mainColumn = -1;       ///< Primary vertical channel for this edge
    std::vector<GridPoint> points;  ///< Routing points
    int secondaryPriority;     ///< Priority for ordering during offset calculation

    void addPoint(int row, int col, int16_t kind = 0) {
        this->points.push_back({row, col, 0, kind, 0});
    }
};

/**
 * @brief Complete state for layout calculation.
 */
struct LayoutState {
    std::unordered_map<uint64_t, GridBlock> grid_blocks;
    std::unordered_map<uint64_t, GraphLayout::GraphBlock>* blocks = nullptr;
    std::unordered_map<uint64_t, std::vector<GridEdge>> edge;
    
    size_t rows = -1;
    size_t columns = -1;
    
    std::vector<int> columnWidth;
    std::vector<int> rowHeight;
    std::vector<int> edgeColumnWidth;
    std::vector<int> edgeRowHeight;

    std::vector<int> columnOffset;
    std::vector<int> rowOffset;
    std::vector<int> edgeColumnOffset;
    std::vector<int> edgeRowOffset;
};

/// Map from block ID to GridBlock
using GridBlockMap = std::unordered_map<uint64_t, GridBlock>;

#endif // GRID_BLOCK_H
