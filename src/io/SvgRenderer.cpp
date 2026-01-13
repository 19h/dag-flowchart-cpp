#include "io/SvgRenderer.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

bool renderSvg(const std::string& filename,
               const GraphWithText& graph,
               int width,
               int height) {
    std::ofstream svg(filename);
    if (!svg.is_open()) {
        std::cerr << "Error: Could not open " << filename << " for writing\n";
        return false;
    }

    svg << "<?xml version=\"1.0\" standalone=\"no\"?>\n"
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
        "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
        // Overall background: very dark gray
        << "<svg width=\"" << width << "\" height=\"" << height
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
        
        // Escape special characters and convert to uppercase
        std::string text = block.text;
        std::transform(text.begin(), text.end(), text.begin(), ::toupper);
        
        // Replace & with &amp;
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

    return true;
}
