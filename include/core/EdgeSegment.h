#ifndef EDGE_SEGMENT_H
#define EDGE_SEGMENT_H

#include <cstdint>

/**
 * @brief Single segment of an edge path.
 * 
 * An edge can be drawn using multiple horizontal and vertical segments.
 * The x/y meaning matches vertical segments; for horizontal segments 
 * the axes are swapped.
 */
struct EdgeSegment {
    int y0;                 ///< Start coordinate
    int y1;                 ///< End coordinate
    int x;                  ///< Position on perpendicular axis
    int edgeIndex;          ///< Index for storing offset result
    int secondaryPriority;  ///< Priority for ordering
    int16_t kind;           ///< Type of segment for routing decisions
    int16_t spacingOverride; ///< Segment spacing override (0 = use default)
};

/**
 * @brief One side of a node for collision detection.
 */
struct NodeSide {
    int x;      ///< Position on primary axis
    int y0;     ///< Start on secondary axis
    int y1;     ///< End on secondary axis
    int size;   ///< Block size in the x axis direction
};

/**
 * @brief Segment representation for linear programming optimization.
 */
struct Segment {
    int x;              ///< Position coordinate
    int variableId;     ///< ID of LP variable this segment belongs to
    int y0;             ///< Start coordinate
    int y1;             ///< End coordinate
};

#endif // EDGE_SEGMENT_H
