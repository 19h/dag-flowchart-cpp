#ifndef GRAPH_LAYOUT_H
#define GRAPH_LAYOUT_H

#include <cstdint>
#include <unordered_map>
#include <vector>

#include "core/PointF.h"

/**
 * @brief Abstract base class for graph layout algorithms.
 * 
 * Provides the interface and common data structures for laying out
 * directed graphs with nodes (blocks) and edges.
 */
class GraphLayout {
public:
    /**
     * @brief Represents a directed edge in the graph.
     */
    struct GraphEdge {
        uint64_t target;      ///< ID of the target node
        Polyline polyline;    ///< Routed path from source to target
        
        enum ArrowDirection { Down, Left, Up, Right, None };
        ArrowDirection arrow = ArrowDirection::Down;

        explicit GraphEdge(uint64_t target) : target(target) {}
    };

    /**
     * @brief Represents a node (block) in the graph.
     */
    struct GraphBlock {
        int x = 0;            ///< X position (set after layout)
        int y = 0;            ///< Y position (set after layout)
        int width = 0;        ///< Block width in pixels
        int height = 0;       ///< Block height in pixels
        uint64_t entry;       ///< Unique identifier for this block
        std::vector<GraphEdge> edges;  ///< Outgoing edges
    };

    /// Graph is a map from block IDs to block data
    using Graph = std::unordered_map<uint64_t, GraphBlock>;

    /**
     * @brief Configuration for layout spacing.
     */
    struct LayoutConfig {
        int blockVerticalSpacing = 40;
        int blockHorizontalSpacing = 20;
        int edgeVerticalSpacing = 10;
        int edgeHorizontalSpacing = 10;
    };

    explicit GraphLayout(const LayoutConfig& layout_config) 
        : layoutConfig(layout_config) {}
    
    virtual ~GraphLayout() {}

    /**
     * @brief Calculate positions for all blocks and route all edges.
     * 
     * @param blocks The graph to lay out (modified in place)
     * @param entry ID of the entry/root node
     * @param width Output: total width of the layout
     * @param height Output: total height of the layout
     */
    virtual void CalculateLayout(Graph& blocks, uint64_t entry, 
                                 int& width, int& height) const = 0;

    virtual void setLayoutConfig(const LayoutConfig& config) {
        this->layoutConfig = config;
    }

protected:
    LayoutConfig layoutConfig;
};

#endif // GRAPH_LAYOUT_H
