#ifndef DOT_PARSER_H
#define DOT_PARSER_H

#include "core/GraphLayout.h"

#include <string>
#include <unordered_map>

/**
 * @brief Extended GraphBlock that includes text label.
 */
struct GraphBlockWithText : public GraphLayout::GraphBlock {
    std::string text;
};

/// Graph type mapping IDs to blocks with text
using GraphWithText = std::unordered_map<uint64_t, GraphBlockWithText>;

/**
 * @brief Parse a DOT format graph file.
 * 
 * Supports basic DOT format with nodes and edges.
 * Node labels are extracted from [label="..."] attributes.
 * 
 * @param filename Path to the DOT file
 * @param graph Output graph (modified in place)
 * @return true if parsing succeeded, false on error
 */
bool parseDotFile(const std::string& filename, GraphWithText& graph);

#endif // DOT_PARSER_H
