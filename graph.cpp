#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <queue>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "BinaryTrees.h"
#include "LinkedListPool.h"

// Our simple replacement for QPointF.
struct PointF {
    double x;
    double y;
    PointF() : x(0), y(0) {}
    PointF(double x, double y) : x(x), y(y) {}
    void setX(double newX) { x = newX; }
    void setY(double newY) { y = newY; }
    // Allow in–place modification like Qt's rx()/ry()
    double& rx() { return x; }
    double& ry() { return y; }
    const double& rx() const { return x; }
    const double& ry() const { return y; }
    PointF& operator-=(const PointF &other) { x -= other.x; y -= other.y; return *this; }
};

// Our replacement for QPolygonF (a vector of PointF).
typedef std::vector<PointF> Polyline;

class GraphLayout
{
public:
    struct GraphEdge
    {
        uint64_t target;
        Polyline polyline;
        enum ArrowDirection { Down, Left, Up, Right, None };
        ArrowDirection arrow = ArrowDirection::Down;

        explicit GraphEdge(uint64_t target) : target(target) {}
    };

    struct GraphBlock
    {
        int x = 0;
        int y = 0;
        int width = 0;
        int height = 0;
        // This is a unique identifier, e.g. offset in the case of rizin blocks
        uint64_t entry;
        // Edges
        std::vector<GraphEdge> edges;
    };
    using Graph = std::unordered_map<uint64_t, GraphBlock>;

    struct LayoutConfig
    {
        int blockVerticalSpacing = 40;
        int blockHorizontalSpacing = 20;
        int edgeVerticalSpacing = 10;
        int edgeHorizontalSpacing = 10;
    };

    GraphLayout(const LayoutConfig &layout_config) : layoutConfig(layout_config) {}
    virtual ~GraphLayout() {}
    virtual void CalculateLayout(Graph &blocks, uint64_t entry, int &width, int &height) const = 0;
    virtual void setLayoutConfig(const LayoutConfig &config) { this->layoutConfig = config; };

protected:
    LayoutConfig layoutConfig;
};


/**
 * @brief Graph layout algorithm on layered graph layout approach. For simplicity all the nodes are
 * placed in a grid.
 */
class GraphGridLayout : public GraphLayout
{
public:
    enum class LayoutType {
        Medium,
        Wide,
        Narrow,
    };

    GraphGridLayout(LayoutType layoutType = LayoutType::Medium);
    virtual void CalculateLayout(Graph &blocks, uint64_t entry, int &width, int &height) const override;
    void setTightSubtreePlacement(bool enabled) { tightSubtreePlacement = enabled; }
    void setParentBetweenDirectChild(bool enabled) { parentBetweenDirectChild = enabled; }
    void setverticalBlockAlignmentMiddle(bool enabled) { verticalBlockAlignmentMiddle = enabled; }
    void setLayoutOptimization(bool enabled) { useLayoutOptimization = enabled; }

private:
    // Configuration flags:
    bool tightSubtreePlacement = false;
    bool parentBetweenDirectChild = false;
    bool verticalBlockAlignmentMiddle = false;
    bool useLayoutOptimization = true;

    // Internal grid–block representation.
    struct GridBlock {
        uint64_t id;
        std::vector<uint64_t> tree_edge;   // outgoing tree–edges
        std::vector<uint64_t> dag_edge;    // outgoing DAG–edges
        std::size_t has_parent = false;
        int inputCount = 0;
        int outputCount = 0;

        int row_count = 0;
        int col = 0;
        int row = 0;

        uint64_t mergeBlock = 0;

        int lastRowLeft;   // left side of last row of subtree
        int lastRowRight;  // right side of last row of subtree
        int leftPosition;
        int rightPosition;
        LinkedListPool<int>::List leftSideShape;
        LinkedListPool<int>::List rightSideShape;
    };

    struct Point {
        int row;
        int col;
        int offset;
        int16_t kind;
        int16_t spacingOverride;
    };

    struct GridEdge {
        uint64_t dest;
        int mainColumn = -1;
        std::vector<Point> points;
        int secondaryPriority;

        void addPoint(int row, int col, int16_t kind = 0) {
            this->points.push_back({ row, col, 0, kind, 0 });
        }
    };

    struct LayoutState {
        std::unordered_map<uint64_t, GridBlock> grid_blocks;
        std::unordered_map<uint64_t, GraphBlock> *blocks = nullptr;
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

    using GridBlockMap = std::unordered_map<uint64_t, GridBlock>;

    // Internal helper methods.
    void findMergePoints(LayoutState &state) const;
    void computeAllBlockPlacement(const std::vector<uint64_t> &blockOrder,
                                  LayoutState &layoutState) const;
    static std::vector<uint64_t> topoSort(LayoutState &state, uint64_t entry);
    static void assignRows(LayoutState &state, const std::vector<uint64_t> &blockOrder);
    static void selectTree(LayoutState &state);
    void routeEdges(LayoutState &state) const;
    void calculateEdgeMainColumn(LayoutState &state) const;
    void roughRouting(LayoutState &state) const;
    void elaborateEdgePlacement(LayoutState &state) const;
    void adjustColumnWidths(LayoutState &state) const;
    static int calculateColumnOffsets(const std::vector<int> &columnWidth,
                                      std::vector<int> &edgeColumnWidth,
                                      std::vector<int> &columnOffset,
                                      std::vector<int> &edgeColumnOffset);
    void convertToPixelCoordinates(LayoutState &state, int &width, int &height) const;
    void cropToContent(Graph &graph, int &width, int &height) const;
    void connectEdgeEnds(Graph &graph) const;
    void optimizeLayout(LayoutState &state) const;
};


/** @class GraphGridLayout

Basic familiarity with graph algorithms is recommended.

# Terms used:
- **Vertex**, **node**, **block** - see the definition of graph. Within this text
vertex/node/block are used interchangeably due to the code being purposed for visualizing basic
block control flow graph.
- **edge** - see the definition of graph.
- **DAG** - directed acyclic graph, a graph using directed edges which doesn't have cycles. A DAG
may contain loops if following them would require going in both directions of edges. Example 1->2
1->3 3->2 is a DAG, 2->1 1->3 3->2 isn't a DAG.
- **DFS** - depth first search, a graph traversal algorithm
- **toposort** - topological sorting, the process of ordering a DAG vertices that results in all
edges going from vertices earlier in the toposort order to vertices later in toposort order. There
are multiple algorithms implementing toposort. A single DAG can have multiple valid topological
orderings, a toposort algorithm can be designed to prioritize a specific one from all valid toposort
orders. Example: for graph 1->4, 2->1, 2->3, 3->4 valid topological orders are [2,1,3,4] and
[2,3,1,4].

# High level algorithm structure
1. Select a subset of edges that form a DAG (remove cycles)
2. Toposort the DAG
3. Choose a subset of edges that form a tree and assign layers
4. Assign node positions within grid using tree structure, child subtrees are placed side by side
with parent on top
5. Perform edge routing
6. Calculate column and row pixel positions based on node sizes and amount edges between the rows
7. [optional] Layout compacting


Contrary to many other layered graph-drawing algorithms this implementation doesn't perform node
reordering to minimize edge crossing. This simplifies the implementation, and preserves the original
control-flow structure for conditional jumps ( true jump on one side, false jump on other). Due to
most of the control flow resulting from structured programming constructs like if/then/else and
loops, the resulting layout is usually readable without node reordering within layers.


# Grid
To simplify the layout algorithm, its initial steps assume that all nodes have the same size and
that edges are zero-width. After nodes placement and edges rounting, the row/column of nodes is
known as well as the amount of edges between each pair of rows. Using this information, positions
are converted from grid cells to pixel coordinates. Routing zero-width edges between rows can also
be interpreted as every second row and column being reserved for edges. The row numbers in code are
using the first interpretation. To allow better centering of nodes one above other, each node is 2
columns wide and 1 row high.

\image html graph_grid.svg

# 1-2 Cycle removal and toposort

Cycle removal and toposort are done in a single DFS traversal. In case the entrypoint
is part of a loop, the DFS starts from the entrypoint. This ensures that the entrypoint is at the
top of resulting layout, if possible. The resulting toposort order is used in many of the following
layout steps that require calculating some property of a vertex based on a child property or the
other way around. Using toposort order, such operations can be implemented by array iteration in
either forward/backward direction. To prevent running out of stack memory when processing large
graphs, DFS is implemented non-recursively.

# Row assignment

Rows are assigned in toposort order from top to bottom, with nodes row being max(predecessor.row)+1.
This ensures that loop back-edges are the only edges going from lower to higher layers.

To further simply node placement, a subset of edges is selected which forms a tree. This turns a DAG
drawing problem into a tree drawing problem. For each node in level n the following nodes with
level exactly n+1 are greedily assigned as child nodes in the tree. If a node already has a parent
assigned then the corresponding edge is not part of the tree.

# Node placement

Since the graph has been reduced to a tree, node placement is more or less putting subtrees side by
side with parent on top. There is some room for interpretation as to what exactly 'side by side'
means and where exactly 'on top' is: drawing the graph either too dense or too sparse may make it
less readable, so there are configuration options which allow choosing these things resulting in
more or less dense layout.

Once the subtrees are placed side by side, the parent node can be placed either in the middle of
the horizontal bounds or in the middle of its direct children. The first option results in narrower
layout and more vertical columns, while the second option results in more spread out layout which
may help seeing where each edge goes.

In compact mode two subtrees are placed side by side accounting for their shape. In wider
mode the bounding box of the shorter subtree is used instead of its exact shape. This gives slightly
sparser layout without being too wide.

\image html graph_parent_placement.svg

# Edge routing
Edge routing can be split into: main column selection, rough routing, and segment offset
calculation.

Transition from source to target row is done using a single vertical segment. This segment is called
the 'main column'.

Main columns are computed using a sweep line: blocks and edges are processed as events top to
bottom based off their row (max(start row, end row) for edges). Blocked columns are tracked in a
tree structure which allows searching nearest column with at least last N rows empty. The column
of the starting block is favored for the main column, otherwise the target block's column is chosen
if it is not blocked. If both the source and target columns are blocked, nearest unblocked column
is chosen. An empty column can always be found, in the worst case there are empty columns at the
sides of drawing. If two columns are equally close, the tie is broken based on whether the edge is a
true or false branch. In case of upward edges it is allowed to choose a column on the outside which
is slightly further than nearest empty to reduce the chance of producing tilted figure 8 shaped
crossing between two blocks.

Due to nodes being placed in a grid, horizontal segments of edges can't intersect with any nodes.
The path for edges is chosen so that it consists of at most 5 segments, typically resulting in
sideways U shape or square Z shape:
- short vertical segment from node to horizontal line
- move to empty column
- vertical segment between starting row and end row
- horizontal segment to target node column
- short vertical segment connecting to target node

There are 3 special cases:
- source and target nodes are in the same column with no nodes between - single vertical segment
- column bellow stating node is empty - segments 1-3 are merged
- column above target node is empty - segments 3-5 are merged

After rough routing segment offsets are calculated relative to their corresponding edge column. This
ensures that two segments don't overlap. Segment offsets within each column are assigned greedily
with some heuristics for assignment order to reduce amount of edge crossings and result in more
visually pleasing output for a typical CFG graph. Each segment gets assigned an offset that is
maximum of previously assigned offsets overlapping with current segment + segment spacing.

Assignment order is based on:
- direction of previous and last segment - helps reducing crossings and place the segments between
nodes
- segment length - reduces crossing when segment endpoints have the same structure as valid
parentheses expression
- edge length - establishes some kind of order when single node is connected to many edges,
typically a block with switch statement or block after switch statement.

# Layout compacting

Doing the layout on a grid limits the minimal spacing to the widest block within a column and
tallest block within a row. One common case is a function-entry block being wider due to the
function name, causing wide horizontal space between branching blocks. Another case is rows in two
parallel columns being aligned.

\image html layout_compacting.svg

Both problems are mitigated by squishing the graph. Compressing in each of the two direction is done
separately. The process is defined as liner program. Each variable represents a position of edge
segment or node in the direction being optimized.

The following constraints are used:
- Keep the order with nearest segments.
- If a node has two outgoing edges, one to the left and one to the right, keep them
on the corresponding side of the node's center.
- Equality constraint to keep relative position between nodes and and segments directly connected to
them.
- For all blocks connected by forward edge, keep the vertical distance at least as big as configured
block vertical spacing. This helps when vertical block-spacing is set bigger than double edge
spacing and an edge shadows relationship between two blocks.
- Equality constraint to keep a node centered when control flow merges.

In the vertical direction the objective function minimizes y positions of nodes and lengths of
vertical segments. In the horizontal direction the objective function minimizes the lengths of
horizontal segments.

In the resulting linear program all constraints besides x_i >= 0 consist of exactly two variables:
either x_i - x_j <= c_k or x_i = x_j + c_k.

Since a perfect solution isn't necessary and to avoid worst case performance, the current
implementation isn't using a general purpose linear solver. Instead, each variable is modified
until a constraint is satisfied and afterwards variables are grouped and modified together.

*/

GraphGridLayout::GraphGridLayout(GraphGridLayout::LayoutType layoutType)
    : GraphLayout({}) // empty config; user may set config later
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

std::vector<uint64_t> GraphGridLayout::topoSort(LayoutState &state, uint64_t entry) 
{
    auto &blocks = *state.blocks;
    std::vector<uint64_t> blockOrder;
    enum class State : uint8_t { NotVisited = 0, InStack, Visited };
    std::unordered_map<uint64_t, State> visited;
    visited.reserve(state.blocks->size());
    std::stack<std::pair<uint64_t, size_t>> stack;

    auto dfsFragment = [&](uint64_t first) {
        visited[first] = State::InStack;
        stack.push({first, 0});

        while (!stack.empty()) {
            auto [v, edge_index] = stack.top();
            const auto &block = blocks[v];

            if (edge_index < block.edges.size()) {
                ++stack.top().second;
                auto target = block.edges[edge_index].target;
                auto &targetState = visited[target];

                if (targetState == State::NotVisited) {
                    targetState = State::InStack;
                    stack.push({target, 0});
                    state.grid_blocks[v].dag_edge.push_back(target);
                } else if (targetState == State::Visited) {
                    state.grid_blocks[v].dag_edge.push_back(target);
                }
            } else {
                stack.pop();
                visited[v] = State::Visited;
                blockOrder.push_back(v);
            }
        }
    };

    dfsFragment(entry);
    for (const auto &[id, _] : blocks) {
        if (visited[id] == State::NotVisited) {
            dfsFragment(id);
        }
    }

    return blockOrder;
}

void GraphGridLayout::assignRows(LayoutState &state, const std::vector<uint64_t> &blockOrder)
{
    // Traverse blocks in reverse topological order to assign row levels
    for (auto it = blockOrder.rbegin(); it != blockOrder.rend(); ++it) {
        const auto &block = state.grid_blocks[*it];
        const int nextLevel = block.row + 1;

        // Ensure all child nodes are at least one row below their parent
        for (const auto target : block.dag_edge) {
            auto &targetBlock = state.grid_blocks[target];
            targetBlock.row = std::max(targetBlock.row, nextLevel);
        }
    }
}

void GraphGridLayout::selectTree(LayoutState &state)
{
    // Select a spanning tree from the DAG edges
    for (auto &[_, block] : state.grid_blocks) {
        for (const auto targetId : block.dag_edge) {
            auto &targetBlock = state.grid_blocks[targetId];
            
            // Only select edges to nodes in the next row that don't already have a parent
            if (!targetBlock.has_parent && targetBlock.row == block.row + 1) {
                // Need non-const reference to modify tree_edge
                auto &mutableBlock = const_cast<GridBlock&>(block);
                mutableBlock.tree_edge.push_back(targetId);
                targetBlock.has_parent = true;
            }
        }
    }
}

void GraphGridLayout::CalculateLayout(GraphLayout::Graph &blocks, uint64_t entry, int &width,
                                    int &height) const
{
    if (blocks.empty()) {
        return;
    }

    LayoutState layoutState;
    layoutState.blocks = &blocks;

    // Use first block as entry if specified entry doesn't exist
    if (blocks.find(entry) == blocks.end()) {
        entry = blocks.begin()->first;
    }

    // Initialize grid blocks
    for (const auto &[id, _] : blocks) {
        GridBlock block;
        block.id = id;
        layoutState.grid_blocks[id] = block;
    }

    // Calculate layout
    auto blockOrder = topoSort(layoutState, entry);
    computeAllBlockPlacement(blockOrder, layoutState);

    // Set up edge information
    for (auto &[blockId, block] : blocks) {
        layoutState.edge[blockId].resize(block.edges.size());
        for (size_t i = 0; i < block.edges.size(); i++) {
            layoutState.edge[blockId][i].dest = block.edges[i].target;
            block.edges[i].arrow = GraphEdge::Down;
        }
    }

    // Calculate input/output counts
    for (const auto &[blockId, edges] : layoutState.edge) {
        auto &startBlock = layoutState.grid_blocks[blockId];
        startBlock.outputCount = edges.size();
        for (const auto &edge : edges) {
            auto &targetBlock = layoutState.grid_blocks[edge.dest];
            targetBlock.inputCount++;
        }
    }

    // Calculate grid dimensions
    layoutState.columns = 1;
    layoutState.rows = 1;
    for (const auto &[_, node] : layoutState.grid_blocks) {
        layoutState.rows = std::max(layoutState.rows, size_t(node.row) + 1);
        layoutState.columns = std::max(layoutState.columns, size_t(node.col) + 2);
    }

    // Calculate row heights and column widths
    layoutState.rowHeight.assign(layoutState.rows, 0);
    layoutState.columnWidth.assign(layoutState.columns, 0);
    for (const auto &[id, node] : layoutState.grid_blocks) {
        const auto &block = blocks[id];
        layoutState.rowHeight[node.row] = 
            std::max(block.height, layoutState.rowHeight[node.row]);
        layoutState.columnWidth[node.col] =
            std::max(block.width / 2, layoutState.columnWidth[node.col]);
        layoutState.columnWidth[node.col + 1] =
            std::max(block.width / 2, layoutState.columnWidth[node.col + 1]);
    }

    // Finalize layout
    routeEdges(layoutState);
    convertToPixelCoordinates(layoutState, width, height);
    
    if (useLayoutOptimization) {
        optimizeLayout(layoutState);
        cropToContent(blocks, width, height);
    }
}

void GraphGridLayout::findMergePoints(GraphGridLayout::LayoutState &state) const 
{
    for (auto &[blockId, block] : state.grid_blocks) {
        // Find merge point by looking at grandchildren
        GridBlock *mergeBlock = nullptr;
        int grandChildCount = 0;
        
        for (const auto edge : block.tree_edge) {
            const auto &targetBlock = state.grid_blocks[edge];
            if (!targetBlock.tree_edge.empty()) {
                mergeBlock = &state.grid_blocks[targetBlock.tree_edge[0]];
            }
            grandChildCount += targetBlock.tree_edge.size();
        }

        // Skip if no valid merge point found
        if (!mergeBlock || grandChildCount != 1) {
            continue;
        }

        // Count blocks that merge to the same point
        int blocksGoingToMerge = 0;
        int blockWithTreeEdge = 0;
        
        for (const auto edge : block.tree_edge) {
            const auto &targetBlock = state.grid_blocks[edge];
            
            // Check if this block connects to merge point
            const bool goesToMerge = std::any_of(
                targetBlock.dag_edge.begin(),
                targetBlock.dag_edge.end(),
                [&](const auto &target) { return target == mergeBlock->id; }
            );

            if (goesToMerge) {
                if (targetBlock.tree_edge.size() == 1) {
                    blockWithTreeEdge = blocksGoingToMerge;
                }
                blocksGoingToMerge++;
            } else {
                break;
            }
        }

        // Update merge block and column position
        if (blocksGoingToMerge > 0) {
            block.mergeBlock = mergeBlock->id;
            state.grid_blocks[block.tree_edge[blockWithTreeEdge]].col = 
                blockWithTreeEdge * 2 - (blocksGoingToMerge - 1);
        }
    }
}

void GraphGridLayout::computeAllBlockPlacement(const std::vector<uint64_t>& blockOrder,
                                             LayoutState& layoutState) const 
{
    assignRows(layoutState, blockOrder);
    selectTree(layoutState);
    findMergePoints(layoutState);

    // Shapes of subtrees are maintained using linked lists. Each value within list is column
    // relative to previous row. This allows moving things around by changing only first value.
    LinkedListPool<int> sides(blockOrder.size() * 2); // Two sides per node

    // Process nodes bottom-to-top to ensure subtrees are processed before parents
    for (auto blockId : blockOrder) {
        auto& block = layoutState.grid_blocks[blockId];

        // Handle leaf nodes
        if (block.tree_edge.empty()) {
            block.row_count = 1;
            block.col = 0;
            block.lastRowRight = 2;
            block.lastRowLeft = 0;
            block.leftPosition = 0;
            block.rightPosition = 2;
            block.leftSideShape = sides.makeList(0);
            block.rightSideShape = sides.makeList(2);
            continue;
        }

        // Handle internal nodes
        auto& firstChild = layoutState.grid_blocks[block.tree_edge[0]];
        auto leftSide = firstChild.leftSideShape;
        auto rightSide = firstChild.rightSideShape;
        
        // Initialize from first child
        block.row_count = firstChild.row_count;
        block.lastRowRight = firstChild.lastRowRight;
        block.lastRowLeft = firstChild.lastRowLeft;
        block.leftPosition = firstChild.leftPosition;
        block.rightPosition = firstChild.rightPosition;

        // Place children subtrees side by side
        for (size_t i = 1; i < block.tree_edge.size(); i++) {
            auto& child = layoutState.grid_blocks[block.tree_edge[i]];
            int minPos = INT_MIN;
            int leftPos = 0;
            int rightPos = 0;
            auto leftIt = sides.head(rightSide);
            auto rightIt = sides.head(child.leftSideShape);
            int maxLeftWidth = 0;
            int minRightPos = child.col;

            // Process overlapping parts of subtrees
            while (leftIt && rightIt) {
                leftPos += *leftIt;
                rightPos += *rightIt;
                minPos = std::max(minPos, leftPos - rightPos);
                maxLeftWidth = std::max(maxLeftWidth, leftPos);
                minRightPos = std::min(minRightPos, rightPos);
                ++leftIt;
                ++rightIt;
            }

            // Calculate offset for right subtree
            int rightTreeOffset;
            if (tightSubtreePlacement) {
                rightTreeOffset = minPos; // Place subtrees as close as possible
            } else {
                // Use bounding box for shortest subtree
                rightTreeOffset = leftIt ? 
                    maxLeftWidth - child.leftPosition :
                    block.rightPosition - minRightPos;
            }

            // Update child position and shapes
            child.col += rightTreeOffset;
            if (leftIt) {
                *leftIt -= (rightTreeOffset + child.lastRowRight - leftPos);
                rightSide = sides.append(child.rightSideShape, 
                                       sides.splitTail(rightSide, leftIt));
            } else if (rightIt) {
                *rightIt += (rightPos + rightTreeOffset - block.lastRowLeft);
                leftSide = sides.append(leftSide, 
                                      sides.splitTail(child.leftSideShape, rightIt));
                rightSide = child.rightSideShape;
                block.lastRowRight = child.lastRowRight + rightTreeOffset;
                block.lastRowLeft = child.lastRowLeft + rightTreeOffset;
            } else {
                rightSide = child.rightSideShape;
            }

            *sides.head(rightSide) += rightTreeOffset;
            block.row_count = std::max(block.row_count, child.row_count);
            block.leftPosition = std::min(block.leftPosition, 
                                        child.leftPosition + rightTreeOffset);
            block.rightPosition = std::max(block.rightPosition,
                                         rightTreeOffset + child.rightPosition);
        }

        // Calculate parent position
        int col;
        if (parentBetweenDirectChild) {
            // Average of children positions
            col = 0;
            for (auto target : block.tree_edge) {
                col += layoutState.grid_blocks[target].col;
            }
            col /= block.tree_edge.size();
        } else {
            // Between leftmost and rightmost children
            col = (block.rightPosition + block.leftPosition) / 2 - 1;
            col = std::max(col, layoutState.grid_blocks[block.tree_edge.front()].col - 1);
            col = std::min(col, layoutState.grid_blocks[block.tree_edge.back()].col + 1);
        }

        // Update block position and shape
        block.col += col;
        block.row_count += 1;
        block.leftPosition = std::min(block.leftPosition, block.col);
        block.rightPosition = std::max(block.rightPosition, block.col + 2);

        *sides.head(leftSide) -= block.col;
        block.leftSideShape = sides.append(sides.makeList(block.col), leftSide);

        *sides.head(rightSide) -= block.col + 2;
        block.rightSideShape = sides.append(sides.makeList(block.col + 2), rightSide);

        // Make children positions relative to parent
        for (auto target : block.tree_edge) {
            layoutState.grid_blocks[target].col -= block.col;
        }
    }

    // Place root nodes and calculate final positions
    int nextEmptyColumn = 0;
    for (auto& [_, block] : layoutState.grid_blocks) {
        if (block.row == 0) {
            auto offset = -block.leftPosition;
            block.col += nextEmptyColumn + offset;
            nextEmptyColumn = block.rightPosition + offset + nextEmptyColumn;
        }
    }

    // Convert relative positions to absolute (top-to-bottom)
    for (auto it = blockOrder.rbegin(); it != blockOrder.rend(); ++it) {
        auto& block = layoutState.grid_blocks[*it];
        assert(block.col >= 0);
        for (auto childId : block.tree_edge) {
            layoutState.grid_blocks[childId].col += block.col;
        }
    }
}
void GraphGridLayout::routeEdges(GraphGridLayout::LayoutState& state) const {
    calculateEdgeMainColumn(state);
    roughRouting(state);
    elaborateEdgePlacement(state);
}

void GraphGridLayout::calculateEdgeMainColumn(GraphGridLayout::LayoutState& state) const {
    struct Event {
        uint64_t blockId;
        size_t edgeId;
        int row;
        enum Type { Edge = 0, Block = 1 } type;
    };

    // Create events for blocks and edges
    std::vector<Event> events;
    events.reserve(state.grid_blocks.size() * 2);

    for (const auto& [blockId, gridBlock] : state.grid_blocks) {
        events.push_back({blockId, 0, gridBlock.row, Event::Block});

        const auto& inputBlock = (*state.blocks)[blockId];
        int startRow = gridBlock.row + 1;

        auto& gridEdges = state.edge[blockId];
        gridEdges.resize(inputBlock.edges.size());

        for (size_t i = 0; i < inputBlock.edges.size(); i++) {
            auto targetId = inputBlock.edges[i].target;
            gridEdges[i].dest = targetId;
            const auto& targetGridBlock = state.grid_blocks[targetId];
            int endRow = targetGridBlock.row;
            events.push_back({blockId, i, std::max(startRow, endRow), Event::Edge});
        }
    }

    // Sort events by row and type
    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        if (a.row != b.row) {
            return a.row < b.row;
        }
        return static_cast<int>(a.type) < static_cast<int>(b.type);
    });

    // Process events to choose main column for each edge
    PointSetMinTree blockedColumns(state.columns + 1, -1);

    for (const auto& event : events) {
        if (event.type == Event::Block) {
            const auto& block = state.grid_blocks[event.blockId];
            blockedColumns.set(block.col + 1, event.row);
            continue;
        }

        const auto& block = state.grid_blocks[event.blockId];
        int column = block.col + 1;
        auto& edge = state.edge[event.blockId][event.edgeId];
        const auto& targetBlock = state.grid_blocks[edge.dest];
        auto topRow = std::min(block.row + 1, targetBlock.row);
        auto targetColumn = targetBlock.col + 1;

        // Try to use source block's column first
        if (blockedColumns.valueAtPoint(column) < topRow) {
            edge.mainColumn = column;
            continue;
        }

        // Try target block's column next
        if (blockedColumns.valueAtPoint(targetColumn) < topRow) {
            edge.mainColumn = targetColumn;
            continue;
        }

        // Find nearest available columns on both sides
        auto nearestLeft = blockedColumns.rightMostLessThan(column, topRow);
        auto nearestRight = blockedColumns.leftMostLessThan(column, topRow);
        assert(nearestLeft != -1 && nearestRight != -1);

        // Calculate distances including path to target
        auto distanceLeft = column - nearestLeft + abs(targetColumn - nearestLeft);
        auto distanceRight = nearestRight - column + abs(targetColumn - nearestRight);

        // Special handling for upward edges to avoid crossings
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

        // Choose side based on total distance or edge index
        if (distanceLeft != distanceRight) {
            edge.mainColumn = (distanceLeft < distanceRight) ? nearestLeft : nearestRight;
        } else {
            edge.mainColumn = (event.edgeId < state.edge[event.blockId].size() / 2) 
                ? nearestLeft 
                : nearestRight;
        }
    }
}

void GraphGridLayout::roughRouting(GraphGridLayout::LayoutState &state) const
{
    // Helper to calculate edge spacing based on block width and edge count
    auto getSpacingOverride = [this](int blockWidth, int edgeCount) {
        if (edgeCount == 0) return 0;
        
        int maxSpacing = blockWidth / edgeCount;
        if (maxSpacing < layoutConfig.edgeHorizontalSpacing) {
            return std::max(maxSpacing, 1);
        }
        return 0;
    };

    // Route each edge from source to target block
    for (const auto &[blockId, block] : state.grid_blocks) {
        auto &blockEdges = state.edge[blockId];
        
        for (auto &edge : blockEdges) {
            const auto &start = block;
            const auto &target = state.grid_blocks[edge.dest];

            // Add initial point at source block
            edge.addPoint(start.row + 1, start.col + 1);

            // Route from source to main column if different
            if (edge.mainColumn != start.col + 1) {
                edge.addPoint(start.row + 1, start.col + 1, 
                    edge.mainColumn < start.col + 1 ? -1 : 1);
                edge.addPoint(start.row + 1, edge.mainColumn,
                    target.row <= start.row ? -2 : 0);
            }

            // Determine main column routing type
            int mainColumnKind = 0;
            if (edge.mainColumn < start.col + 1 && edge.mainColumn < target.col + 1) {
                mainColumnKind = +2;
            } else if (edge.mainColumn > start.col + 1 && edge.mainColumn > target.col + 1) {
                mainColumnKind = -2;
            } else if (edge.mainColumn == start.col + 1 && edge.mainColumn != target.col + 1) {
                mainColumnKind = edge.mainColumn < target.col + 1 ? 1 : -1;
            } else if (edge.mainColumn == target.col + 1 && edge.mainColumn != start.col + 1) {
                mainColumnKind = edge.mainColumn < start.col + 1 ? 1 : -1;
            }

            // Add main column point
            edge.addPoint(target.row, edge.mainColumn, mainColumnKind);

            // Route from main column to target if different
            if (target.col + 1 != edge.mainColumn) {
                edge.addPoint(target.row, target.col + 1,
                    target.row <= start.row ? 2 : 0);
                edge.addPoint(target.row, target.col + 1,
                    target.col + 1 < edge.mainColumn ? 1 : -1);
            }

            // Handle edge spacing overrides for blocks with many edges
            auto startSpacingOverride = getSpacingOverride(
                (*state.blocks)[start.id].width, start.outputCount);
            auto targetSpacingOverride = getSpacingOverride(
                (*state.blocks)[target.id].width, target.inputCount);

            edge.points.front().spacingOverride = startSpacingOverride;
            edge.points.back().spacingOverride = targetSpacingOverride;

            if (edge.points.size() <= 2) {
                if (startSpacingOverride && startSpacingOverride < targetSpacingOverride) {
                    edge.points.back().spacingOverride = startSpacingOverride;
                }
            } else {
                edge.points[1].spacingOverride = startSpacingOverride;
            }

            // Calculate total path length for priority
            int length = 0;
            for (size_t i = 1; i < edge.points.size(); i++) {
                length += abs(edge.points[i].row - edge.points[i - 1].row) +
                         abs(edge.points[i].col - edge.points[i - 1].col);
            }
            edge.secondaryPriority = 2 * length + (target.row >= start.row ? 1 : 0);
        }
    }
}

namespace {
/**
 * @brief Single segment of an edge. An edge can be drawn using multiple horizontal and vertical
 * segments. x y meaning matches vertical segments. For horizontal segments axis are swapped.
 */
struct EdgeSegment
{
    int y0;
    int y1;
    int x;
    int edgeIndex;
    int secondaryPriority;
    int16_t kind;
    int16_t spacingOverride; //< segment spacing override, 0 if default spacing should be used
};
struct NodeSide
{
    int x;
    int y0;
    int y1;
    int size; //< block size in the x axis direction
};
}

/**
 * @brief Calculate segment offsets relative to their column
 *
 * Argument naming uses terms for vertical segments, but the function can be used for horizontal
 * segments as well.
 *
 * @param segments Segments that need to be processed.
 * @param edgeOffsets Output argument for returning segment offsets relative to their columns.
 * @param edgeColumnWidth InOut argument describing how much column with edges take. Initial value
 * used as minimal value. May be increased depending on amount of segments in each column and how
 * tightly they are packed.
 * @param nodeRightSide Right side of nodes. Used to reduce space reserved for edges by placing them
 * between nodes.
 * @param nodeLeftSide Same as right side.
 * @param columnWidth Width of each column
 * @param H All the segment and node coordinates y0 and y1 are expected to be in range [0;H)
 * @param segmentSpacing The expected spacing between two segments in the same column. Actual
 * spacing may be smaller for nodes with many edges.
 */
void calculateSegmentOffsets(std::vector<EdgeSegment>& segments, 
                           std::vector<int>& edgeOffsets,
                           std::vector<int>& edgeColumnWidth,
                           std::vector<NodeSide>& nodeRightSide,
                           std::vector<NodeSide>& nodeLeftSide,
                           const std::vector<int>& columnWidth, 
                           size_t H, 
                           int segmentSpacing) {
    // Ensure y0 <= y1 for all segments
    for (auto& segment : segments) {
        if (segment.y0 > segment.y1) {
            std::swap(segment.y0, segment.y1);
        }
    }

    // Sort segments by x, kind, size and priority
    std::sort(segments.begin(), segments.end(), [](const EdgeSegment& a, const EdgeSegment& b) {
        if (a.x != b.x) return a.x < b.x;
        if (a.kind != b.kind) return a.kind < b.kind;
        
        auto aSize = a.y1 - a.y0;
        auto bSize = b.y1 - b.y0;
        if (aSize != bSize) {
            return (a.kind != 1) ? aSize < bSize : aSize > bSize;
        }
        if (a.kind != 1) {
            return a.secondaryPriority < b.secondaryPriority;
        }
        return a.secondaryPriority > b.secondaryPriority;
    });

    // Sort node sides by x coordinate
    auto compareNode = [](const NodeSide& a, const NodeSide& b) { return a.x < b.x; };
    std::sort(nodeRightSide.begin(), nodeRightSide.end(), compareNode);
    std::sort(nodeLeftSide.begin(), nodeLeftSide.end(), compareNode);

    RangeAssignMaxTree maxSegment(H, INT_MIN);
    auto nextSegmentIt = segments.begin();
    auto rightSideIt = nodeRightSide.begin();
    auto leftSideIt = nodeLeftSide.begin();

    // Process segments column by column
    while (nextSegmentIt != segments.end()) {
        int x = nextSegmentIt->x;

        // Handle left side of column
        int leftColumnWidth = (x > 0) ? columnWidth[x - 1] : 0;
        maxSegment.setRange(0, H, -leftColumnWidth);

        // Process right side nodes from previous column
        while (rightSideIt != nodeRightSide.end() && rightSideIt->x + 1 <= x) {
            if (rightSideIt->x + 1 == x) {
                maxSegment.setRange(rightSideIt->y0, rightSideIt->y1 + 1,
                                  rightSideIt->size - leftColumnWidth);
            }
            rightSideIt++;
        }

        // Process left side segments
        while (nextSegmentIt != segments.end() && nextSegmentIt->x == x && nextSegmentIt->kind <= 1) {
            int y = maxSegment.rangeMaximum(nextSegmentIt->y0, nextSegmentIt->y1 + 1);
            if (nextSegmentIt->kind != -2) {
                y = std::max(y, 0);
            }
            y += nextSegmentIt->spacingOverride ? nextSegmentIt->spacingOverride : segmentSpacing;
            maxSegment.setRange(nextSegmentIt->y0, nextSegmentIt->y1 + 1, y);
            edgeOffsets[nextSegmentIt->edgeIndex] = y;
            nextSegmentIt++;
        }

        // Handle right side of column
        auto firstRightSideSegment = nextSegmentIt;
        auto middleWidth = std::max(maxSegment.rangeMaximum(0, H), 0);
        int rightColumnWidth = (x < static_cast<int>(columnWidth.size())) ? columnWidth[x] : 0;

        maxSegment.setRange(0, H, -rightColumnWidth);

        // Process left side nodes in current column
        while (leftSideIt != nodeLeftSide.end() && leftSideIt->x <= x) {
            if (leftSideIt->x == x) {
                maxSegment.setRange(leftSideIt->y0, leftSideIt->y1 + 1,
                                  leftSideIt->size - rightColumnWidth);
            }
            leftSideIt++;
        }

        // Process right side segments
        while (nextSegmentIt != segments.end() && nextSegmentIt->x == x) {
            int y = maxSegment.rangeMaximum(nextSegmentIt->y0, nextSegmentIt->y1 + 1);
            y += nextSegmentIt->spacingOverride ? nextSegmentIt->spacingOverride : segmentSpacing;
            maxSegment.setRange(nextSegmentIt->y0, nextSegmentIt->y1 + 1, y);
            edgeOffsets[nextSegmentIt->edgeIndex] = y;
            nextSegmentIt++;
        }

        // Adjust right side segment positions
        auto rightSideMiddle = std::max(maxSegment.rangeMaximum(0, H), 0);
        rightSideMiddle = std::max(rightSideMiddle, edgeColumnWidth[x] - middleWidth - segmentSpacing);
        
        for (auto it = firstRightSideSegment; it != nextSegmentIt; ++it) {
            edgeOffsets[it->edgeIndex] = 
                middleWidth + (rightSideMiddle - edgeOffsets[it->edgeIndex]) + segmentSpacing;
        }
        
        edgeColumnWidth[x] = middleWidth + segmentSpacing + rightSideMiddle;
    }
}

/**
 * @brief Center the segments to the middle of edge columns when possible.
 * @param segmentOffsets offsets relative to the left side edge column.
 * @param edgeColumnWidth widths of edge columns
 * @param segments either all horizontal or all vertical edge segments
 */
static void centerEdges(std::vector<int>& segmentOffsets, 
                       const std::vector<int>& edgeColumnWidth,
                       const std::vector<EdgeSegment>& segments) {
    // Split segments in each edge column into non-intersecting chunks.
    // Center each chunk separately by processing segment endpoints sorted by x and y.
    // Track active segment count - when it reaches 0, there is a gap between chunks.
    struct Event {
        int x;
        int y; 
        int index;
        bool start;
    };

    std::vector<Event> events;
    events.reserve(segments.size() * 2);

    // Create events for segment start/end points
    for (const auto& segment : segments) {
        auto offset = segmentOffsets[segment.edgeIndex];
        
        // Skip segments outside edge column bounds since they're hard to adjust safely
        if (offset >= 0 && offset <= edgeColumnWidth[segment.x]) {
            events.push_back({segment.x, segment.y0, segment.edgeIndex, true});
            events.push_back({segment.x, segment.y1, segment.edgeIndex, false}); 
        }
    }

    // Sort by x, y, with starts before ends
    std::sort(events.begin(), events.end(), [](const Event& a, const Event& b) {
        if (a.x != b.x) return a.x < b.x;
        if (a.y != b.y) return a.y < b.y;
        return int(a.start) > int(b.start);
    });

    // Process events to find and center chunks
    auto it = events.begin();
    while (it != events.end()) {
        int left = segmentOffsets[it->index];
        int right = left;
        auto chunkStart = it++;
        int activeSegments = 1;

        // Find chunk bounds
        while (activeSegments > 0) {
            activeSegments += it->start ? 1 : -1;
            int offset = segmentOffsets[it->index];
            left = std::min(left, offset);
            right = std::max(right, offset);
            it++;
        }

        // Center the chunk
        int spacing = (edgeColumnWidth[chunkStart->x] - (right - left)) / 2 - left;
        for (auto segment = chunkStart; segment != it; segment++) {
            if (segment->start) {
                segmentOffsets[segment->index] += spacing;
            }
        }
    }
}

/**
 * @brief Convert segment coordinates from arbitrary range to continuous range starting at 0
 * 
 * Takes segment y-coordinates in an arbitrary range and maps them to a continuous range
 * starting at 0 while preserving their relative ordering.
 *
 * @param segments Edge segments to compress coordinates for
 * @param leftSides Left sides of nodes to compress coordinates for 
 * @param rightSides Right sides of nodes to compress coordinates for
 * @return Size of the compressed coordinate range
 */
static int compressCoordinates(std::vector<EdgeSegment>& segments,
                             std::vector<NodeSide>& leftSides,
                             std::vector<NodeSide>& rightSides) {
    // Collect all unique y-coordinates
    std::vector<int> positions;
    positions.reserve((segments.size() + leftSides.size()) * 2);
    
    for (const auto& segment : segments) {
        positions.push_back(segment.y0);
        positions.push_back(segment.y1);
    }
    
    for (const auto& side : leftSides) {
        positions.push_back(side.y0);
        positions.push_back(side.y1);
    }

    // Sort and remove duplicates
    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());

    // Map old coordinates to new compressed indices
    auto positionToIndex = [&](int position) {
        auto index = std::lower_bound(positions.begin(), positions.end(), position) - positions.begin();
        assert(index < positions.size());
        return static_cast<int>(index);
    };

    // Update segment coordinates
    for (auto& segment : segments) {
        segment.y0 = positionToIndex(segment.y0);
        segment.y1 = positionToIndex(segment.y1);
    }

    // Update node side coordinates
    assert(leftSides.size() == rightSides.size());
    for (size_t i = 0; i < leftSides.size(); i++) {
        auto newY0 = positionToIndex(leftSides[i].y0);
        auto newY1 = positionToIndex(leftSides[i].y1);
        leftSides[i].y0 = rightSides[i].y0 = newY0;
        leftSides[i].y1 = rightSides[i].y1 = newY1;
    }

    return static_cast<int>(positions.size());
}

void GraphGridLayout::elaborateEdgePlacement(GraphGridLayout::LayoutState &state) const
{
    int edgeIndex = 0;

    // Helper to create edge segments from points
    auto segmentFromPoint = [&edgeIndex](const Point &point, const GridEdge &edge, int y0, int y1, int x) {
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
    std::vector<NodeSide> rightSides;
    std::vector<NodeSide> leftSides;
    std::vector<int> edgeOffsets;

    // Create vertical segments
    for (auto &edgeListIt : state.edge) {
        for (const auto &edge : edgeListIt.second) {
            for (size_t j = 1; j < edge.points.size(); j += 2) {
                segments.push_back(
                    segmentFromPoint(edge.points[j], edge,
                                   edge.points[j - 1].row * 2,  // edges in even rows
                                   edge.points[j].row * 2, 
                                   edge.points[j].col)
                );
            }
        }
    }

    // Create node sides for blocks
    for (auto &blockIt : state.grid_blocks) {
        auto &node = blockIt.second;
        auto width = (*state.blocks)[blockIt.first].width;
        auto leftWidth = width / 2;
        auto rightWidth = width - leftWidth;
        int row = node.row * 2 + 1;  // blocks in odd rows
        
        leftSides.push_back({node.col, row, row, leftWidth});
        rightSides.push_back({node.col + 1, row, row, rightWidth});
    }

    // Initialize column widths and calculate offsets
    state.edgeColumnWidth.assign(state.columns + 1, layoutConfig.blockHorizontalSpacing);
    state.edgeColumnWidth[0] = state.edgeColumnWidth.back() = layoutConfig.edgeHorizontalSpacing;
    edgeOffsets.resize(edgeIndex);

    calculateSegmentOffsets(segments, edgeOffsets, state.edgeColumnWidth, rightSides, leftSides,
                          state.columnWidth, 2 * state.rows + 1,
                          layoutConfig.edgeHorizontalSpacing);
    
    centerEdges(edgeOffsets, state.edgeColumnWidth, segments);
    edgeIndex = 0;

    // Helper to copy segment offsets back to edges
    auto copySegmentsToEdges = [&](bool col) {
        int edgeIndex = 0;
        for (auto &edgeListIt : state.edge) {
            for (auto &edge : edgeListIt.second) {
                for (size_t j = col ? 1 : 2; j < edge.points.size(); j += 2) {
                    int offset = edgeOffsets[edgeIndex++];
                    
                    if (col) {
                        GraphBlock *block = nullptr;
                        if (j == 1) {
                            block = &(*state.blocks)[edgeListIt.first];
                        } else if (j + 1 == edge.points.size()) {
                            block = &(*state.blocks)[edge.dest];
                        }
                        
                        if (block) {
                            int blockWidth = block->width;
                            int edgeColumnWidth = state.edgeColumnWidth[edge.points[j].col];
                            offset = std::max(-blockWidth / 2 + edgeColumnWidth / 2, offset);
                            offset = std::min(edgeColumnWidth / 2 + std::min(blockWidth, edgeColumnWidth) / 2,
                                           offset);
                        }
                    }
                    edge.points[j].offset = offset;
                }
            }
        }
    };

    // Adjust column widths and offsets
    auto oldColumnWidths = state.columnWidth;
    adjustColumnWidths(state);

    for (auto &segment : segments) {
        auto &offset = edgeOffsets[segment.edgeIndex];
        if (segment.kind == -2) {
            offset -= (state.edgeColumnWidth[segment.x - 1] / 2 + state.columnWidth[segment.x - 1])
                     - oldColumnWidths[segment.x - 1];
        } else if (segment.kind == 2) {
            offset += (state.edgeColumnWidth[segment.x + 1] / 2 + state.columnWidth[segment.x])
                     - oldColumnWidths[segment.x];
        }
    }

    calculateColumnOffsets(state.columnWidth, state.edgeColumnWidth, state.columnOffset,
                         state.edgeColumnOffset);
    copySegmentsToEdges(true);

    // Create horizontal segments
    segments.clear();
    leftSides.clear();
    rightSides.clear();
    edgeIndex = 0;

    for (auto &edgeListIt : state.edge) {
        for (const auto &edge : edgeListIt.second) {
            for (size_t j = 2; j < edge.points.size(); j += 2) {
                int y0 = state.edgeColumnOffset[edge.points[j - 1].col] + edge.points[j - 1].offset;
                int y1 = state.edgeColumnOffset[edge.points[j + 1].col] + edge.points[j + 1].offset;
                segments.push_back(
                    segmentFromPoint(edge.points[j], edge, y0, y1, edge.points[j].row)
                );
            }
        }
    }

    edgeOffsets.resize(edgeIndex);

    // Create node sides for blocks with exact coordinates
    for (auto &blockIt : state.grid_blocks) {
        auto &node = blockIt.second;
        auto blockWidth = (*state.blocks)[node.id].width;
        int leftSide = state.edgeColumnOffset[node.col + 1] + 
                      state.edgeColumnWidth[node.col + 1] / 2 - blockWidth / 2;
        int rightSide = leftSide + blockWidth;

        int h = (*state.blocks)[blockIt.first].height;
        int freeSpace = state.rowHeight[node.row] - h;
        int topProfile = state.rowHeight[node.row];
        int bottomProfile = h;

        if (verticalBlockAlignmentMiddle) {
            topProfile -= freeSpace / 2;
            bottomProfile += freeSpace / 2;
        }

        leftSides.push_back({node.row, leftSide, rightSide, topProfile});
        rightSides.push_back({node.row, leftSide, rightSide, bottomProfile});
    }

    // Calculate final horizontal segment positions
    state.edgeRowHeight.assign(state.rows + 1, layoutConfig.blockVerticalSpacing);
    state.edgeRowHeight[0] = state.edgeRowHeight.back() = layoutConfig.edgeVerticalSpacing;
    edgeOffsets.resize(edgeIndex);

    auto compressedCoordinates = compressCoordinates(segments, leftSides, rightSides);
    calculateSegmentOffsets(segments, edgeOffsets, state.edgeRowHeight, rightSides, leftSides,
                          state.rowHeight, compressedCoordinates,
                          layoutConfig.edgeVerticalSpacing);
    
    copySegmentsToEdges(false);
}

void GraphGridLayout::adjustColumnWidths(GraphGridLayout::LayoutState &state) const
{
    // Initialize row heights and column widths to 0
    state.rowHeight.assign(state.rows, 0);
    state.columnWidth.assign(state.columns, 0);

    // Process each node to determine maximum row heights and column widths
    for (const auto &node : state.grid_blocks) {
        const auto &block = (*state.blocks)[node.first];
        const auto &gridPos = node.second;

        // Update maximum height for this row
        state.rowHeight[gridPos.row] = std::max(block.height, state.rowHeight[gridPos.row]);

        // Calculate width needed for each column around the node
        int edgeWidth = state.edgeColumnWidth[gridPos.col + 1];
        int columnWidth = (block.width - edgeWidth) / 2;

        // Update maximum widths for the two columns this node spans
        state.columnWidth[gridPos.col] = std::max(columnWidth, state.columnWidth[gridPos.col]);
        state.columnWidth[gridPos.col + 1] = std::max(columnWidth, state.columnWidth[gridPos.col + 1]);
    }
}

int GraphGridLayout::calculateColumnOffsets(const std::vector<int>& columnWidth,
                                          std::vector<int>& edgeColumnWidth,
                                          std::vector<int>& columnOffset,
                                          std::vector<int>& edgeColumnOffset) {
    assert(edgeColumnWidth.size() == columnWidth.size() + 1);

    int position = 0;
    edgeColumnOffset.resize(edgeColumnWidth.size());
    columnOffset.resize(columnWidth.size());

    // Calculate offsets alternating between edge columns and regular columns
    for (size_t i = 0; i < columnWidth.size(); i++) {
        edgeColumnOffset[i] = position;
        position += edgeColumnWidth[i];

        columnOffset[i] = position; 
        position += columnWidth[i];
    }

    // Handle final edge column
    edgeColumnOffset.back() = position;
    position += edgeColumnWidth.back();

    return position;
}

void GraphGridLayout::convertToPixelCoordinates(LayoutState &state, int &width, int &height) const 
{
    // Calculate final pixel offsets for columns and rows
    width = calculateColumnOffsets(state.columnWidth, state.edgeColumnWidth, 
                                 state.columnOffset, state.edgeColumnOffset);
    height = calculateColumnOffsets(state.rowHeight, state.edgeRowHeight,
                                  state.rowOffset, state.edgeRowOffset);

    // Position blocks in pixel coordinates
    for (auto &block : (*state.blocks)) {
        const auto &gridBlock = state.grid_blocks[block.first];
        auto &blockPos = block.second;

        // Center block horizontally in its grid cell
        blockPos.x = state.edgeColumnOffset[gridBlock.col + 1] + 
                    state.edgeColumnWidth[gridBlock.col + 1] / 2 -
                    blockPos.width / 2;

        // Position block vertically, optionally centering in cell
        blockPos.y = state.rowOffset[gridBlock.row];
        if (verticalBlockAlignmentMiddle) {
            blockPos.y += (state.rowHeight[gridBlock.row] - blockPos.height) / 2;
        }
    }

    // Convert edge routing points to pixel coordinates
    for (auto &[blockId, block] : (*state.blocks)) {
        for (size_t i = 0; i < block.edges.size(); i++) {
            auto &resultEdge = block.edges[i];
            const auto &routingPoints = state.edge[blockId][i];

            // Start edge from bottom of block
            resultEdge.polyline.clear();
            resultEdge.polyline.push_back(PointF(0, block.y + block.height));

            // Convert routing points to pixel coordinates
            for (size_t j = 1; j < routingPoints.points.size(); j++) {
                const auto &point = routingPoints.points[j];
                
                if (j & 1) { // Vertical segment
                    int x = state.edgeColumnOffset[point.col] + point.offset;
                    resultEdge.polyline.back().setX(x);
                    resultEdge.polyline.push_back(PointF(x, 0));
                } else { // Horizontal segment  
                    int y = state.edgeRowOffset[point.row] + point.offset;
                    resultEdge.polyline.back().setY(y);
                    resultEdge.polyline.push_back(PointF(0, y));
                }
            }
        }
    }

    connectEdgeEnds(*state.blocks);
}

void GraphGridLayout::cropToContent(Graph &graph, int &width, int &height) const
{
    // Handle empty graph case
    if (graph.empty()) {
        width = std::max(1, layoutConfig.edgeHorizontalSpacing);
        height = std::max(1, layoutConfig.edgeVerticalSpacing);
        return;
    }

    // Initialize bounds from first block
    const auto &firstBlock = graph.begin()->second;
    int minPos[2] = {firstBlock.x, firstBlock.y};
    int maxPos[2] = {firstBlock.x, firstBlock.y};

    // Helper to update bounding box
    auto updateBounds = [&](int x, int y) {
        minPos[0] = std::min(minPos[0], x);
        minPos[1] = std::min(minPos[1], y);
        maxPos[0] = std::max(maxPos[0], x);
        maxPos[1] = std::max(maxPos[1], y);
    };

    // Find bounds of all blocks and edges
    for (const auto &[id, block] : graph) {
        // Include block bounds
        updateBounds(block.x, block.y);
        updateBounds(block.x + block.width, block.y + block.height);

        // Include edge polyline points
        for (const auto &edge : block.edges) {
            for (const auto &point : edge.polyline) {
                updateBounds(static_cast<int>(point.x), static_cast<int>(point.y));
            }
        }
    }

    // Add spacing margins
    minPos[0] -= layoutConfig.edgeHorizontalSpacing;
    minPos[1] -= layoutConfig.edgeVerticalSpacing;
    maxPos[0] += layoutConfig.edgeHorizontalSpacing;
    maxPos[1] += layoutConfig.edgeVerticalSpacing;

    // Shift all coordinates to start at origin
    for (auto &[id, block] : graph) {
        // Shift block position
        block.x -= minPos[0];
        block.y -= minPos[1];

        // Shift edge polylines
        for (auto &edge : block.edges) {
            for (auto &point : edge.polyline) {
                point -= PointF(minPos[0], minPos[1]);
            }
        }
    }

    // Set final dimensions
    width = maxPos[0] - minPos[0];
    height = maxPos[1] - minPos[1];
}

void GraphGridLayout::connectEdgeEnds(Graph& graph) const {
    for (auto& [id, block] : graph) {
        for (auto& edge : block.edges) {
            const auto& target = graph[edge.target];
            edge.polyline.front().ry() = block.y + block.height;
            edge.polyline.back().ry() = target.y;
        }
    }
}

/// Either equality or inequality x_i <= x_j + c
using Constraint = std::pair<std::pair<int, int>, int>;

/**
 * @brief Single pass of linear program optimizer.
 * 
 * Changes variables until a constraint is hit, afterwards the two variables are changed together.
 * 
 * @param n Number of variables
 * @param objectiveFunction Coefficients for function \f$\sum c_i x_i\f$ which needs to be minimized
 * @param inequalities Inequality constraints \f$x_{e_i} - x_{f_i} \leq b_i\f$
 * @param equalities Equality constraints \f$x_{e_i} - x_{f_i} = b_i\f$
 * @param solution Input/output argument, returns results, needs to be initialized with a feasible solution
 * @param stickWhenNotMoving Variable grouping strategy
 */
static void optimizeLinearProgramPass(
    size_t n,
    std::vector<int> objectiveFunction,
    std::vector<Constraint> inequalities,
    std::vector<Constraint> equalities,
    std::vector<int>& solution,
    bool stickWhenNotMoving
) {
    // Initialize data structures
    std::vector<int> group(n);
    std::iota(group.begin(), group.end(), 0); // Initially each variable is in its own group
    
    assert(n == objectiveFunction.size());
    assert(n == solution.size());
    
    std::vector<size_t> edgeCount(n);
    LinkedListPool<size_t> edgePool(inequalities.size() * 2);
    std::vector<decltype(edgePool)::List> edges(n);
    std::vector<uint8_t> processed(n);
    
    // Smallest variable value in the group relative to main one
    // Used to maintain implicit x_i >= 0 constraint
    std::vector<int> groupRelativeMin(n, 0);

    // Helper functions
    auto getGroup = [&](int v) {
        while (group[v] != v) {
            group[v] = group[group[v]];
            v = group[v];
        }
        return v;
    };

    auto joinGroup = [&](int a, int b) { 
        group[getGroup(b)] = getGroup(a); 
    };

    // Build edge lists
    for (auto& constraint : inequalities) {
        int a = constraint.first.first;
        int b = constraint.first.second;
        size_t index = &constraint - &inequalities.front();
        
        edges[a] = edgePool.append(edges[a], edgePool.makeList(index));
        edges[b] = edgePool.append(edges[b], edgePool.makeList(index));
        edgeCount[a]++;
        edgeCount[b]++;
    }

    auto joinSegmentGroups = [&](int a, int b) {
        a = getGroup(a);
        b = getGroup(b);
        joinGroup(a, b);
        
        edgeCount[a] += edgeCount[b];
        objectiveFunction[a] += objectiveFunction[b];
        
        int internalEdgeCount = 0;
        auto writeIt = edgePool.head(edges[b]);
        
        // Update inequalities and remove constraints between grouped variables
        for (auto it = edgePool.head(edges[b]); it; ++it) {
            auto& constraint = inequalities[*it];
            int other = constraint.first.first + constraint.first.second - b;
            
            if (getGroup(other) == a) {
                // Skip inequalities where both variables are now in same group
                internalEdgeCount++;
                continue;
            }
            
            *writeIt++ = *it;
            
            // Modify inequalities for attached group relative to main group
            int diff = solution[a] - solution[b];
            if (b == constraint.first.first) {
                constraint.first.first = a;
                constraint.second += diff;
            } else {
                constraint.first.second = a;
                constraint.second -= diff;
            }
        }
        
        edges[a] = edgePool.append(edges[a], edgePool.splitHead(edges[b], writeIt));
        edgeCount[a] -= internalEdgeCount;
        groupRelativeMin[a] = std::min(groupRelativeMin[a], 
                                     groupRelativeMin[b] + solution[b] - solution[a]);
    };

    // Process equalities
    for (auto& equality : equalities) {
        int a = getGroup(equality.first.first);
        int b = getGroup(equality.first.second);
        
        if (a == b) {
            equality = {{0, 0}, 0};
            continue;
        }
        
        // Always join smallest group to bigger one
        if (edgeCount[a] > edgeCount[b]) {
            std::swap(a, b);
            std::swap(equality.first.first, equality.first.second);
            equality.second = -equality.second;
        }
        
        joinSegmentGroups(b, a);
        equality = {{a, b}, solution[a] - solution[b]};
        processed[a] = 1;
    }

    // Process groups starting with smallest
    std::priority_queue<
        std::pair<int, int>,
        std::vector<std::pair<int, int>>,
        std::greater<std::pair<int, int>>
    > queue;
    
    for (size_t i = 0; i < n; i++) {
        if (!processed[i]) {
            queue.push({edgeCount[i], i});
        }
    }

    while (!queue.empty()) {
        int g = queue.top().second;
        int size = queue.top().first;
        queue.pop();
        
        if ((size_t)size != edgeCount[g] || processed[g]) {
            continue;
        }
        
        int direction = objectiveFunction[g];
        if (direction == 0) {
            continue;
        }

        // Find first constraint hit by changing variable in desired direction
        int limitingGroup = -1;
        int smallestMove = 0;
        
        if (direction < 0) {
            smallestMove = INT_MAX;
            for (auto it = edgePool.head(edges[g]); it; ++it) {
                auto& inequality = inequalities[*it];
                if (g == inequality.first.second) {
                    continue;
                }
                
                int other = inequality.first.second;
                if (getGroup(other) == g) {
                    continue;
                }
                
                int move = solution[other] + inequality.second - solution[g];
                if (move < smallestMove) {
                    smallestMove = move;
                    limitingGroup = other;
                }
            }
        } else {
            smallestMove = -solution[g] - groupRelativeMin[g]; // Keep all variables >= 0
            for (auto it = edgePool.head(edges[g]); it; ++it) {
                auto& inequality = inequalities[*it];
                if (g == inequality.first.first) {
                    continue;
                }
                
                int other = inequality.first.first;
                if (getGroup(other) == g) {
                    continue;
                }
                
                int move = solution[other] - inequality.second - solution[g];
                if (move > smallestMove) {
                    smallestMove = move;
                    limitingGroup = other;
                }
            }
        }

        assert(smallestMove != INT_MAX);
        if (smallestMove == INT_MAX) {
            // Unbound variable - linear program wasn't set up correctly
            // Don't change it instead of stretching graph to infinity
            smallestMove = 0;
        }
        
        solution[g] += smallestMove;
        if (smallestMove == 0 && !stickWhenNotMoving) {
            continue;
        }
        
        processed[g] = 1;
        if (limitingGroup != -1) {
            joinSegmentGroups(limitingGroup, g);
            if (!processed[limitingGroup]) {
                queue.push({edgeCount[limitingGroup], limitingGroup});
            }
            equalities.push_back({{g, limitingGroup}, solution[g] - solution[limitingGroup]});
        }
    }

    // Update solutions based on equalities
    for (auto it = equalities.rbegin(); it != equalities.rend(); ++it) {
        solution[it->first.first] = solution[it->first.second] + it->second;
    }
}

/**
 * @brief Linear programming solver that optimizes a linear objective function subject to constraints
 * 
 * Note: Does not guarantee finding the optimal solution, but finds a feasible solution that improves
 * the objective value.
 *
 * @param n Number of variables in the linear program
 * @param objectiveFunction Coefficients ci for minimizing sum(ci * xi) 
 * @param inequalities Constraints of form x[ei] - x[fi] <= bi
 * @param equalities Constraints of form x[ei] - x[fi] = bi
 * @param solution Input: Initial feasible solution, Output: Optimized solution
 */
static void optimizeLinearProgram(size_t n,
                                const std::vector<int>& objectiveFunction,
                                std::vector<Constraint> inequalities, 
                                const std::vector<Constraint>& equalities,
                                std::vector<int>& solution) {
    // Remove redundant inequalities with same variable pairs
    std::sort(inequalities.begin(), inequalities.end());
    auto uniqueEnd = std::unique(inequalities.begin(), inequalities.end(),
                               [](const Constraint& a, const Constraint& b) {
                                   return a.first == b.first;
                               });
    inequalities.erase(uniqueEnd, inequalities.end());

    // Run single optimization pass
    static const int ITERATIONS = 1;
    for (int i = 0; i < ITERATIONS; i++) {
        optimizeLinearProgramPass(n, objectiveFunction, inequalities, equalities, solution, true);
    }
}

namespace {
struct Segment {
    int x;          // x-coordinate of segment
    int variableId; // ID of variable this segment belongs to
    int y0;         // Start y-coordinate
    int y1;         // End y-coordinate
};
}

/**
 * @brief Creates an inequality constraint between two variables
 *
 * @param a First variable index
 * @param posA Position of first variable
 * @param b Second variable index  
 * @param posB Position of second variable
 * @param minSpacing Minimum required spacing between variables
 * @param positions Vector of current positions
 * @return Constraint representing a - b <= value
 */
static Constraint createInequality(size_t a, int posA, size_t b, int posB, int minSpacing,
                                 const std::vector<int>& positions) {
    minSpacing = std::min(minSpacing, posB - posA);
    return {{a, b}, posB - positions[b] - (posA - positions[a]) - minSpacing};
}

/**
 * @brief Create inequality constraints from segments which preserves their relative order on single axis.
 *
 * @param segments List of edge segments and block sides
 * @param positions Initial element positions before optimization
 * @param blockCount Number of variables representing blocks. Segments with variableId < blockCount 
 *                  represent one side of a block.
 * @param variableGroup Used to check if segments are part of the same edge and spacing can be reduced
 * @param blockSpacing Minimal spacing between blocks
 * @param segmentSpacing Minimal spacing between two edge segments. Spacing may be less if values in
 *                      positions are closer than this.
 * @param inequalities Output vector for resulting inequalities. Initial values are preserved.
 */
static void createInequalitiesFromSegments(std::vector<Segment> segments,
                                         const std::vector<int>& positions,
                                         const std::vector<size_t>& variableGroup,
                                         int blockCount,
                                         int blockSpacing,
                                         int segmentSpacing,
                                         std::vector<Constraint>& inequalities) {
    // Map used as binary search tree: y_position -> segment{variableId, x_position}
    // Used to track which segment was last seen at each y position
    std::map<int, std::pair<int, int>> lastSegments;
    lastSegments[-1] = {-1, -1}; // Sentinel value

    // Sort segments by x coordinate
    std::sort(segments.begin(), segments.end(), 
              [](const Segment& a, const Segment& b) { return a.x < b.x; });

    for (const auto& segment : segments) {
        // Find last segment that starts before or at current segment's y0
        auto startPos = lastSegments.lower_bound(segment.y0);
        --startPos; // Safe because map has sentinel at -1
        
        auto lastSegment = startPos->second;
        auto it = startPos;

        // Process all segments that overlap with current segment's y range
        while (it != lastSegments.end() && it->first <= segment.y1) {
            int prevSegmentVariable = it->second.first;
            int prevSegmentPos = it->second.second;

            if (prevSegmentVariable != -1) {
                int minSpacing = segmentSpacing;

                // Handle special spacing cases
                if (prevSegmentVariable < blockCount && segment.variableId < blockCount) {
                    // Skip inequality between two sides of same block
                    if (prevSegmentVariable == segment.variableId) {
                        ++it;
                        continue;
                    }
                    minSpacing = blockSpacing;
                }
                else if (variableGroup[prevSegmentVariable] == variableGroup[segment.variableId]) {
                    minSpacing = 0; // No spacing needed between segments of same edge
                }

                inequalities.push_back(createInequality(prevSegmentVariable, prevSegmentPos,
                                                      segment.variableId, segment.x, 
                                                      minSpacing, positions));
            }

            lastSegment = it->second;
            ++it;
        }

        // Update lastSegments map
        if (startPos->first < segment.y0) {
            startPos++;
        }
        lastSegments.erase(startPos, it); // Remove segments covered by current one
        
        // Add current segment and remaining part of partially covered segment
        lastSegments[segment.y0] = {segment.variableId, segment.x};
        lastSegments[segment.y1] = lastSegment;
    }
}

void GraphGridLayout::optimizeLayout(GraphGridLayout::LayoutState &state) const
{
    // Map block IDs to sequential indices for optimization
    std::unordered_map<uint64_t, int> blockMapping;
    size_t blockIndex = 0;
    for (auto &blockIt : *state.blocks) {
        blockMapping[blockIt.first] = blockIndex++;
    }

    // Initialize variable groups - initially each variable is in its own group
    std::vector<size_t> variableGroups(blockMapping.size());
    std::iota(variableGroups.begin(), variableGroups.end(), 0);

    // Vectors for optimization problem
    std::vector<int> objectiveFunction;
    std::vector<Constraint> inequalities;
    std::vector<Constraint> equalities;
    std::vector<int> solution;

    // Helper functions for building optimization problem
    auto addObjective = [&](size_t a, int posA, size_t b, int posB) {
        objectiveFunction.resize(std::max(objectiveFunction.size(), std::max(a, b) + 1));
        if (posA < posB) {
            objectiveFunction[b] += 1;
            objectiveFunction[a] -= 1;
        } else {
            objectiveFunction[a] += 1;
            objectiveFunction[b] -= 1;
        }
    };

    auto addInequality = [&](size_t a, int posA, size_t b, int posB, int minSpacing) {
        inequalities.push_back(createInequality(a, posA, b, posB, minSpacing, solution));
    };

    auto addBlockSegmentEquality = [&](uint64_t blockId, int edgeVariable, int edgeVariablePos) {
        int blockPos = (*state.blocks)[blockId].x;
        int blockVariable = blockMapping[blockId];
        equalities.push_back({{blockVariable, edgeVariable}, blockPos - edgeVariablePos});
    };

    auto setFeasibleSolution = [&](size_t variable, int value) {
        solution.resize(std::max(solution.size(), variable + 1));
        solution[variable] = value;
    };

    auto copyVariablesToPositions = [&](const std::vector<int> &solution, bool horizontal = false) {
#ifndef NDEBUG
        for (auto v : solution) {
            assert(v >= 0);
        }
#endif
        size_t variableIndex = blockMapping.size();
        for (auto &blockIt : *state.blocks) {
            auto &block = blockIt.second;
            for (auto &edge : block.edges) {
                for (int i = 1 + int(horizontal); i < edge.polyline.size(); i += 2) {
                    int x = solution[variableIndex++];
                    if (horizontal) {
                        edge.polyline[i].ry() = x;
                        edge.polyline[i - 1].ry() = x;
                    } else {
                        edge.polyline[i].rx() = x;
                        edge.polyline[i - 1].rx() = x;
                    }
                }
            }
            int blockVariable = blockMapping[blockIt.first];
            (horizontal ? block.y : block.x) = solution[blockVariable];
        }
    };

    // First optimize vertical positions
    std::vector<Segment> segments;
    segments.reserve(state.blocks->size() * 4);
    size_t variableIndex = blockMapping.size();
    size_t edgeIndex = 0;

    objectiveFunction.assign(blockMapping.size(), 1);

    // Process blocks and edges for vertical optimization
    for (auto &blockIt : *state.blocks) {
        auto &block = blockIt.second;
        int blockVariable = blockMapping[blockIt.first];

        // Add vertical spacing constraints between connected blocks
        for (auto &edge : block.edges) {
            auto &targetBlock = (*state.blocks)[edge.target];
            if (block.y < targetBlock.y) {
                int spacing = block.height + layoutConfig.blockVerticalSpacing;
                inequalities.push_back({{blockVariable, blockMapping[edge.target]}, -spacing});
            }

            // Process edge segments
            if (edge.polyline.size() >= 3) {
                for (int i = 2; i < edge.polyline.size(); i += 2) {
                    int y0 = edge.polyline[i - 1].x;
                    int y1 = edge.polyline[i].x;
                    if (y0 > y1) std::swap(y0, y1);
                    
                    int x = edge.polyline[i].y;
                    segments.push_back({x, int(variableIndex), y0, y1});
                    variableGroups.push_back(blockMapping.size() + edgeIndex);
                    setFeasibleSolution(variableIndex, x);
                    
                    if (i > 2) {
                        int prevX = edge.polyline[i - 2].y;
                        addObjective(variableIndex, x, variableIndex - 1, prevX);
                    }
                    variableIndex++;
                }
                edgeIndex++;
            }
        }

        // Add block segments
        segments.push_back({block.y, blockVariable, block.x, block.x + block.width});
        segments.push_back({block.y + block.height, blockVariable, block.x, block.x + block.width});
        setFeasibleSolution(blockVariable, block.y);
    }

    // Create inequalities from segments and optimize
    createInequalitiesFromSegments(std::move(segments), solution, variableGroups,
                                 blockMapping.size(), layoutConfig.blockVerticalSpacing,
                                 layoutConfig.edgeVerticalSpacing, inequalities);

    objectiveFunction.resize(solution.size());
    optimizeLinearProgram(solution.size(), objectiveFunction, inequalities, equalities, solution);
    copyVariablesToPositions(solution, true);
    connectEdgeEnds(*state.blocks);

    // Reset for horizontal optimization
    variableGroups.resize(blockMapping.size());
    solution.clear();
    equalities.clear();
    inequalities.clear();
    objectiveFunction.clear();
    segments.clear();
    variableIndex = blockMapping.size();
    edgeIndex = 0;

    // Process blocks and edges for horizontal optimization
    for (auto &blockIt : *state.blocks) {
        auto &block = blockIt.second;
        
        // Process edges
        for (auto &edge : block.edges) {
            if (edge.polyline.size() < 2) continue;

            size_t firstEdgeVariable = variableIndex;
            for (int i = 1; i < edge.polyline.size(); i += 2) {
                int y0 = edge.polyline[i - 1].y;
                int y1 = edge.polyline[i].y;
                if (y0 > y1) std::swap(y0, y1);
                
                int x = edge.polyline[i].x;
                segments.push_back({x, int(variableIndex), y0, y1});
                variableGroups.push_back(blockMapping.size() + edgeIndex);
                setFeasibleSolution(variableIndex, x);
                
                if (i > 2) {
                    int prevX = edge.polyline[i - 2].x;
                    addObjective(variableIndex, x, variableIndex - 1, prevX);
                }
                variableIndex++;
            }

            // Add block-segment equalities
            size_t lastEdgeVariableIndex = variableIndex - 1;
            addBlockSegmentEquality(blockIt.first, firstEdgeVariable, edge.polyline[1].x);
            addBlockSegmentEquality(edge.target, lastEdgeVariableIndex, segments.back().x);
            edgeIndex++;
        }

        // Add block segments
        int blockVariable = blockMapping[blockIt.first];
        segments.push_back({block.x, blockVariable, block.y, block.y + block.height});
        segments.push_back({block.x + block.width, blockVariable, block.y, block.y + block.height});
        setFeasibleSolution(blockVariable, block.x);
    }

    // Create inequalities and optimize horizontal positions
    createInequalitiesFromSegments(std::move(segments), solution, variableGroups,
                                 blockMapping.size(), layoutConfig.blockHorizontalSpacing,
                                 layoutConfig.edgeHorizontalSpacing, inequalities);

    objectiveFunction.resize(solution.size());

    // Add horizontal centering constraints
    for (auto &blockIt : *state.blocks) {
        auto &block = blockIt.second;
        int blockVariable = blockMapping[blockIt.first];
        
        if (block.edges.size() == 2) {
            auto &blockLeft = (*state.blocks)[block.edges[0].target];
            auto &blockRight = (*state.blocks)[block.edges[1].target];
            auto middle = block.x + block.width / 2;
            
            if (blockLeft.x + blockLeft.width < middle && blockRight.x > middle) {
                addInequality(blockMapping[block.edges[0].target], blockLeft.x + blockLeft.width,
                            blockVariable, middle, layoutConfig.blockHorizontalSpacing / 2);
                addInequality(blockVariable, middle, blockMapping[block.edges[1].target],
                            blockRight.x, layoutConfig.blockHorizontalSpacing / 2);
                
                auto &gridBlock = state.grid_blocks[blockIt.first];
                if (gridBlock.mergeBlock) {
                    auto &mergeBlock = (*state.blocks)[gridBlock.mergeBlock];
                    if (mergeBlock.x + mergeBlock.width / 2 == middle) {
                        equalities.push_back(
                            {{blockVariable, blockMapping[gridBlock.mergeBlock]},
                             block.x - mergeBlock.x});
                    }
                }
            }
        }
    }

    optimizeLinearProgram(solution.size(), objectiveFunction, inequalities, equalities, solution);
    copyVariablesToPositions(solution);
}

/**
 * @brief Adapter that converts vertical graph layout to horizontal by swapping coordinates
 */
class GraphHorizontalAdapter : public GraphLayout {
public:
    explicit GraphHorizontalAdapter(std::unique_ptr<GraphLayout> layout)
        : GraphLayout({}), layout(std::move(layout))
    {
        swapLayoutConfigDirection();
    }

    void CalculateLayout(Graph &blocks, uint64_t entry, int &width, int &height) const override
    {
        // Swap width/height of all blocks
        for (auto &block : blocks) {
            std::swap(block.second.width, block.second.height);
        }

        // Calculate layout with swapped dimensions
        layout->CalculateLayout(blocks, entry, height, width);

        // Convert back to horizontal layout
        for (auto &block : blocks) {
            std::swap(block.second.width, block.second.height);
            std::swap(block.second.x, block.second.y);
            for (auto &edge : block.second.edges) {
                for (auto &point : edge.polyline) {
                    std::swap(point.rx(), point.ry());
                }
                switch (edge.arrow) {
                case GraphEdge::Down:
                    edge.arrow = GraphEdge::Right;
                    break;
                case GraphEdge::Left:
                    edge.arrow = GraphEdge::Up;
                    break;
                case GraphEdge::Up:
                    edge.arrow = GraphEdge::Left;
                    break;
                case GraphEdge::Right:
                    edge.arrow = GraphEdge::Down;
                    break;
                case GraphEdge::None:
                    edge.arrow = GraphEdge::None;
                    break;
                }
            }
        }
    }

    void setLayoutConfig(const LayoutConfig &config) override
    {
        GraphLayout::setLayoutConfig(config);
        swapLayoutConfigDirection();
        layout->setLayoutConfig(config);
    }

private:
    void swapLayoutConfigDirection()
    {
        std::swap(layoutConfig.edgeVerticalSpacing, layoutConfig.edgeHorizontalSpacing);
        std::swap(layoutConfig.blockVerticalSpacing, layoutConfig.blockHorizontalSpacing);
    }

    std::unique_ptr<GraphLayout> layout;
};

struct GraphBlockWithText : public GraphLayout::GraphBlock {
    std::string text;
};

// Graph type mapping IDs to blocks with text
using GraphWithText = std::unordered_map<uint64_t, GraphBlockWithText>;

// Layout subclass that uses text blocks but keeps same layout algorithm
class GraphGridLayoutWithText : public GraphGridLayout {
public:
    using GraphGridLayout::GraphGridLayout;
};

int main() {
    GraphWithText graph;

    // Parse DOT file
    std::ifstream dotFile("graph.dot");
    if (!dotFile.is_open()) {
        std::cerr << "Error: Could not open graph.dot\n";
        return 1;
    }

    std::string line;
    std::getline(dotFile, line); // Skip "digraph {" line

    uint64_t nodeCounter = 0;
    std::unordered_map<std::string, uint64_t> nodeIds;

    // Parse nodes and edges
    while (std::getline(dotFile, line)) {
        if (line.empty() || line == "}") {
            continue;
        }

        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t"));

        // Skip comments and graph attributes
        if (line[0] == '#' || 
            line.find("subgraph") != std::string::npos ||
            line.find("rankdir") != std::string::npos ||
            line.find("node") != std::string::npos ||
            line.find("compound") != std::string::npos ||
            line.find("bgcolor") != std::string::npos ||
            line.find("label=") == 0) {
            continue;
        }

        if (line.find("->") == std::string::npos) {
            // Node definition
            size_t labelStart = line.find("[label=\"");
            if (labelStart == std::string::npos) {
                continue;
            }

            // Extract node ID and label
            std::string nodeId = line.substr(0, labelStart);
            nodeId.erase(0, nodeId.find_first_not_of(" \t"));
            nodeId.erase(nodeId.find_last_not_of(" \t") + 1);

            size_t labelEnd = line.find("\"", labelStart + 8);
            if (labelEnd == std::string::npos) {
                continue;
            }
            std::string label = line.substr(labelStart + 8, labelEnd - (labelStart + 8));

            // Remove quotes from nodeId if present
            if (nodeId.front() == '"' && nodeId.back() == '"') {
                nodeId = nodeId.substr(1, nodeId.size() - 2);
            }

            // Create node if new
            if (nodeIds.find(nodeId) == nodeIds.end()) {
                nodeIds[nodeId] = ++nodeCounter;
            }

            // Estimate text width (assuming each character is 8 pixels wide)
            int textWidth = std::max(100, static_cast<int>(label.length() * 8));

            // Add node to graph
            GraphBlockWithText node;
            node.entry = nodeIds[nodeId];
            node.width = textWidth;
            node.height = 50;
            node.text = label;
            graph[node.entry] = node;

        } else {
            // Edge definition
            size_t arrowPos = line.find("->");
            if (arrowPos == std::string::npos) {
                continue;
            }

            // Extract source and target nodes
            std::string fromStr = line.substr(0, arrowPos);
            std::string toStr = line.substr(arrowPos + 2);

            // Clean up node strings
            fromStr.erase(0, fromStr.find_first_not_of(" \t"));
            fromStr.erase(fromStr.find_last_not_of(" \t;[") + 1);
            
            toStr.erase(0, toStr.find_first_not_of(" \t"));
            size_t bracketPos = toStr.find("[");
            if (bracketPos != std::string::npos) {
                toStr = toStr.substr(0, bracketPos);
            }
            toStr.erase(toStr.find_last_not_of(" \t;[") + 1);

            // Remove quotes from node strings if present
            if (fromStr.front() == '"' && fromStr.back() == '"') {
                fromStr = fromStr.substr(1, fromStr.size() - 2);
            }
            if (toStr.front() == '"' && toStr.back() == '"') {
                toStr = toStr.substr(1, toStr.size() - 2);
            }

            // Create nodes if needed
            for (const auto& nodeStr : {fromStr, toStr}) {
                if (nodeIds.find(nodeStr) == nodeIds.end()) {
                    // Estimate text width (assuming each character is 8 pixels wide)
                    int textWidth = std::max(100, static_cast<int>(nodeStr.length() * 8));

                    nodeIds[nodeStr] = ++nodeCounter;
                    GraphBlockWithText node;
                    node.entry = nodeIds[nodeStr];
                    node.width = textWidth;
                    node.height = 50;
                    node.text = nodeStr;
                    graph[node.entry] = node;
                }
            }

            // Add edge
            graph[nodeIds[fromStr]].edges.push_back(GraphLayout::GraphEdge(nodeIds[toStr]));
        }
    }

    auto use_horizontal_layout = false;

    // Calculate layout
    //GraphGridLayoutWithText layout(GraphGridLayoutWithText::LayoutType::Medium);
    //int layoutWidth, layoutHeight;
    //layout.CalculateLayout(reinterpret_cast<GraphLayout::Graph&>(graph), 1, layoutWidth, layoutHeight);

    std::unique_ptr<GraphLayout> baseLayout = 
        std::make_unique<GraphGridLayoutWithText>(GraphGridLayoutWithText::LayoutType::Medium);
    
    // Optionally wrap with horizontal adapter
    std::unique_ptr<GraphLayout> layout = use_horizontal_layout ? 
        std::make_unique<GraphHorizontalAdapter>(std::move(baseLayout)) :
        std::move(baseLayout);

    int layoutWidth, layoutHeight;
    layout->CalculateLayout(reinterpret_cast<GraphLayout::Graph&>(graph), 1, 
                          layoutWidth, layoutHeight);

    // Generate SVG output
    std::ofstream svg("graph.svg");
    if (!svg.is_open()) {
        std::cerr << "Error: Could not open graph.svg for writing\n";
        return 1;
    }

    svg << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
        "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
        // Overall background: very dark gray
        << "<svg width=\"" << layoutWidth << "\" height=\"" << layoutHeight
        << "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" "
        << "style=\"background-color:#1e1f22;\">\n";

    // Draw edges
    for (const auto& [_, block] : graph) {
        for (const auto& edge : block.edges) {
            // Use a subtle gray stroke to match lines in the screenshot
            svg << "<polyline points=\"";
            for (const auto& pt : edge.polyline) {
                svg << pt.x << "," << pt.y << " ";
            }
            svg << "\" style=\"fill:none;stroke:#bae67e;stroke-width:1\"/>\n";

            // Arrow in the same subtle gray
            if (edge.polyline.size() >= 2) {
                const auto& p1 = edge.polyline[edge.polyline.size() - 2];
                const auto& p2 = edge.polyline[edge.polyline.size() - 1];
                double dx = p2.x - p1.x;
                double dy = p2.y - p1.y;
                double len = std::sqrt(dx*dx + dy*dy);
                if (len > 0) {
                    dx /= len;
                    dy /= len;
                    double arrowSize = 5;
                    double x1 = p2.x - arrowSize * (dx + dy/2);
                    double y1 = p2.y - arrowSize * (dy - dx/2);
                    double x2 = p2.x - arrowSize * (dx - dy/2);
                    double y2 = p2.y - arrowSize * (dy + dx/2);

                    svg << "<polygon points=\""
                        << p2.x << "," << p2.y << " "
                        << x1 << "," << y1 << " "
                        << x2 << "," << y2
                        << "\" style=\"fill:#bae67e\"/>\n";
                }
            }
        }
    }

    // Draw nodes
    for (const auto& [_, block] : graph) {
        // Block fill: dark gray, stroke: medium gray
        svg << "<rect x=\"" << block.x << "\" y=\"" << block.y
            << "\" width=\"" << block.width << "\" height=\"" << block.height
            << "\" style=\"fill:#1c1f24;stroke:#646464;stroke-width:1\"/>\n"
            // Text: light gray
            << "<text x=\"" << (block.x + block.width / 2)
            << "\" y=\"" << (block.y + block.height / 2)
            << "\" text-anchor=\"middle\" alignment-baseline=\"middle\" "
            << "font-size=\"12\" fill=\"#ed9366\" font-family=\"Jetbrains Mono\" font-weight=\"bold\">";
        // Replace lone & with &amp; but don't touch existing &amp;
        std::string text = block.text;
        std::transform(text.begin(), text.end(), text.begin(), ::toupper);
        // First convert & to &amp;
        size_t pos = 0;
        while ((pos = text.find('&', pos)) != std::string::npos) {
            if (pos + 4 >= text.length() || text.substr(pos, 5) != "&amp;") {
                text.replace(pos, 1, "&amp;");
                pos += 5;
            } else {
                pos++;
            }
        }
        svg << text << "</text>\n";
    }

    svg << "</svg>\n";
    svg.close();

    std::cout << "Graph layout rendered to graph.svg (" 
              << layoutWidth << "x" << layoutHeight << ")\n";
    return 0;
}