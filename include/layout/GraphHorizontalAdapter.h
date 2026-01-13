#ifndef GRAPH_HORIZONTAL_ADAPTER_H
#define GRAPH_HORIZONTAL_ADAPTER_H

#include "core/GraphLayout.h"

#include <memory>
#include <utility>

/**
 * @brief Adapter that converts vertical graph layout to horizontal by swapping coordinates.
 * 
 * Wraps any GraphLayout implementation and transposes the result to produce
 * a left-to-right layout instead of top-to-bottom.
 */
class GraphHorizontalAdapter : public GraphLayout {
public:
    explicit GraphHorizontalAdapter(std::unique_ptr<GraphLayout> layout)
        : GraphLayout({}), layout(std::move(layout))
    {
        swapLayoutConfigDirection();
    }

    void CalculateLayout(Graph& blocks, uint64_t entry, int& width, int& height) const override
    {
        // Swap width/height of all blocks
        for (auto& block : blocks) {
            std::swap(block.second.width, block.second.height);
        }

        // Calculate layout with swapped dimensions
        layout->CalculateLayout(blocks, entry, height, width);

        // Convert back to horizontal layout
        for (auto& block : blocks) {
            std::swap(block.second.width, block.second.height);
            std::swap(block.second.x, block.second.y);
            for (auto& edge : block.second.edges) {
                for (auto& point : edge.polyline) {
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

    void setLayoutConfig(const LayoutConfig& config) override
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

#endif // GRAPH_HORIZONTAL_ADAPTER_H
