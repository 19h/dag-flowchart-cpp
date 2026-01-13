#ifndef SVG_RENDERER_H
#define SVG_RENDERER_H

#include "io/DotParser.h"

#include <string>

/**
 * @brief Render a graph to an SVG file.
 * 
 * Generates an SVG with:
 * - Dark theme background (#1e1f22)
 * - Nodes as rectangles with labels
 * - Edges as polylines with arrowheads
 * 
 * @param filename Output SVG file path
 * @param graph The graph to render (must have layout already calculated)
 * @param width Total width of the layout
 * @param height Total height of the layout
 * @return true if rendering succeeded, false on error
 */
bool renderSvg(const std::string& filename, 
               const GraphWithText& graph,
               int width, 
               int height);

#endif // SVG_RENDERER_H
