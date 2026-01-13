#ifndef SEGMENT_OFFSET_H
#define SEGMENT_OFFSET_H

#include <vector>
#include "core/EdgeSegment.h"

/**
 * @brief Calculate segment offsets relative to their column.
 *
 * Argument naming uses terms for vertical segments, but the function can be 
 * used for horizontal segments as well.
 *
 * @param segments Segments that need to be processed.
 * @param edgeOffsets Output argument for returning segment offsets relative to their columns.
 * @param edgeColumnWidth InOut argument describing how much column with edges take. 
 *        Initial value used as minimal value. May be increased depending on amount 
 *        of segments in each column and how tightly they are packed.
 * @param nodeRightSide Right side of nodes. Used to reduce space reserved for edges 
 *        by placing them between nodes.
 * @param nodeLeftSide Same as right side.
 * @param columnWidth Width of each column
 * @param H All the segment and node coordinates y0 and y1 are expected to be in range [0;H)
 * @param segmentSpacing The expected spacing between two segments in the same column. 
 *        Actual spacing may be smaller for nodes with many edges.
 */
void calculateSegmentOffsets(std::vector<EdgeSegment>& segments,
                            std::vector<int>& edgeOffsets,
                            std::vector<int>& edgeColumnWidth,
                            std::vector<NodeSide>& nodeRightSide,
                            std::vector<NodeSide>& nodeLeftSide,
                            const std::vector<int>& columnWidth,
                            size_t H,
                            int segmentSpacing);

/**
 * @brief Center the segments to the middle of edge columns when possible.
 * 
 * @param segmentOffsets offsets relative to the left side edge column.
 * @param edgeColumnWidth widths of edge columns
 * @param segments either all horizontal or all vertical edge segments
 */
void centerEdges(std::vector<int>& segmentOffsets,
                const std::vector<int>& edgeColumnWidth,
                const std::vector<EdgeSegment>& segments);

/**
 * @brief Convert segment coordinates from arbitrary range to continuous range starting at 0.
 * 
 * Takes segment y-coordinates in an arbitrary range and maps them to a continuous range
 * starting at 0 while preserving their relative ordering.
 *
 * @param segments Edge segments to compress coordinates for
 * @param leftSides Left sides of nodes to compress coordinates for
 * @param rightSides Right sides of nodes to compress coordinates for
 * @return Size of the compressed coordinate range
 */
int compressCoordinates(std::vector<EdgeSegment>& segments,
                       std::vector<NodeSide>& leftSides,
                       std::vector<NodeSide>& rightSides);

#endif // SEGMENT_OFFSET_H
