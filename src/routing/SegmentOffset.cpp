#include "routing/SegmentOffset.h"
#include "core/BinaryTrees.h"

#include <algorithm>
#include <cassert>
#include <climits>
#include <vector>

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

void centerEdges(std::vector<int>& segmentOffsets,
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

int compressCoordinates(std::vector<EdgeSegment>& segments,
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
