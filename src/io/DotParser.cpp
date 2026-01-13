#include "io/DotParser.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <string>

bool parseDotFile(const std::string& filename, GraphWithText& graph) {
    std::ifstream dotFile(filename);
    if (!dotFile.is_open()) {
        std::cerr << "Error: Could not open " << filename << "\n";
        return false;
    }

    std::string line;
    std::getline(dotFile, line); // Skip "digraph {" line

    uint64_t nodeCounter = 0;
    std::unordered_map<std::string, uint64_t> nodeIds;

    // Helper to extract node name from string
    auto extractNodeName = [](std::string s) -> std::string {
        std::string trimmed = s;
        // Trim leading whitespace
        trimmed.erase(0, trimmed.find_first_not_of(" \t"));
        // If the node is quoted, extract the content inside quotes.
        if (!trimmed.empty() && trimmed.front() == '"') {
            size_t closingQuote = trimmed.find('"', 1);
            if (closingQuote != std::string::npos) {
                return trimmed.substr(1, closingQuote - 1);
            }
        }
        // Fallback for unquoted nodes: remove any attribute part starting with " ["
        size_t bracketPos = trimmed.find(" [");
        if (bracketPos != std::string::npos) {
            trimmed = trimmed.substr(0, bracketPos);
        }
        // Trim trailing whitespace, semicolons, etc.
        trimmed.erase(trimmed.find_last_not_of(" \t;") + 1);
        return trimmed;
    };

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
            ((line.size() >= 4) && line.substr(0, 4) == "node" &&
            (line.size() == 4 || isspace(line[4]) || line[4]=='[')) ||
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

            fromStr = extractNodeName(fromStr);
            toStr = extractNodeName(toStr);

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

    return true;
}
