#include "io/DotParser.h"
#include "layout/GraphGridLayout.h"
#include "layout/GraphHorizontalAdapter.h"
#include "io/SvgRenderer.h"

#include <iostream>
#include <memory>

/**
 * @brief Layout subclass that uses text blocks but keeps same layout algorithm.
 */
class GraphGridLayoutWithText : public GraphGridLayout {
public:
    using GraphGridLayout::GraphGridLayout;
};

int main() {
    GraphWithText graph;

    // Parse DOT file
    if (!parseDotFile("graph.dot", graph)) {
        return 1;
    }

    auto useHorizontalLayout = false;

    // Create layout engine
    std::unique_ptr<GraphLayout> baseLayout =
        std::make_unique<GraphGridLayoutWithText>(GraphGridLayoutWithText::LayoutType::Medium);

    // Optionally wrap with horizontal adapter
    std::unique_ptr<GraphLayout> layout = useHorizontalLayout ?
        std::make_unique<GraphHorizontalAdapter>(std::move(baseLayout)) :
        std::move(baseLayout);

    // Calculate layout
    int layoutWidth, layoutHeight;
    layout->CalculateLayout(reinterpret_cast<GraphLayout::Graph&>(graph), 1,
                          layoutWidth, layoutHeight);

    // Generate SVG output
    if (!renderSvg("graph.svg", graph, layoutWidth, layoutHeight)) {
        return 1;
    }

    std::cout << "Graph layout rendered to graph.svg ("
              << layoutWidth << "x" << layoutHeight << ")\n";

    return 0;
}
