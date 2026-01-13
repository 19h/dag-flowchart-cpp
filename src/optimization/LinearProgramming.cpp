#include "optimization/LinearProgramming.h"
#include "core/LinkedListPool.h"

#include <algorithm>
#include <cassert>
#include <climits>
#include <map>
#include <numeric>
#include <queue>
#include <vector>

/**
 * @brief Single pass of linear program optimizer.
 *
 * Changes variables until a constraint is hit, afterwards the two variables 
 * are changed together.
 */
static void optimizeLinearProgramPass(
    size_t n,
    std::vector<int> objectiveFunction,
    std::vector<Constraint> inequalities,
    std::vector<Constraint> equalities,
    std::vector<int>& solution,
    bool stickWhenNotMoving
) {
    // Initialize data structures
    std::vector<int> group(n);
    std::iota(group.begin(), group.end(), 0); // Initially each variable is in its own group

    assert(n == objectiveFunction.size());
    assert(n == solution.size());

    std::vector<size_t> edgeCount(n);
    LinkedListPool<size_t> edgePool(inequalities.size() * 2);
    std::vector<decltype(edgePool)::List> edges(n);
    std::vector<uint8_t> processed(n);

    // Smallest variable value in the group relative to main one
    // Used to maintain implicit x_i >= 0 constraint
    std::vector<int> groupRelativeMin(n, 0);

    // Helper functions
    auto getGroup = [&](int v) {
        while (group[v] != v) {
            group[v] = group[group[v]];
            v = group[v];
        }
        return v;
    };

    auto joinGroup = [&](int a, int b) {
        group[getGroup(b)] = getGroup(a);
    };

    // Build edge lists
    for (auto& constraint : inequalities) {
        int a = constraint.first.first;
        int b = constraint.first.second;
        size_t index = &constraint - &inequalities.front();

        edges[a] = edgePool.append(edges[a], edgePool.makeList(index));
        edges[b] = edgePool.append(edges[b], edgePool.makeList(index));
        edgeCount[a]++;
        edgeCount[b]++;
    }

    auto joinSegmentGroups = [&](int a, int b) {
        a = getGroup(a);
        b = getGroup(b);
        joinGroup(a, b);

        edgeCount[a] += edgeCount[b];
        objectiveFunction[a] += objectiveFunction[b];

        int internalEdgeCount = 0;
        auto writeIt = edgePool.head(edges[b]);

        // Update inequalities and remove constraints between grouped variables
        for (auto it = edgePool.head(edges[b]); it; ++it) {
            auto& constraint = inequalities[*it];
            int other = constraint.first.first + constraint.first.second - b;

            if (getGroup(other) == a) {
                // Skip inequalities where both variables are now in same group
                internalEdgeCount++;
                continue;
            }

            *writeIt++ = *it;

            // Modify inequalities for attached group relative to main group
            int diff = solution[a] - solution[b];
            if (b == constraint.first.first) {
                constraint.first.first = a;
                constraint.second += diff;
            } else {
                constraint.first.second = a;
                constraint.second -= diff;
            }
        }

        edges[a] = edgePool.append(edges[a], edgePool.splitHead(edges[b], writeIt));
        edgeCount[a] -= internalEdgeCount;
        groupRelativeMin[a] = std::min(groupRelativeMin[a],
                                     groupRelativeMin[b] + solution[b] - solution[a]);
    };

    // Process equalities
    for (auto& equality : equalities) {
        int a = getGroup(equality.first.first);
        int b = getGroup(equality.first.second);

        if (a == b) {
            equality = {{0, 0}, 0};
            continue;
        }

        // Always join smallest group to bigger one
        if (edgeCount[a] > edgeCount[b]) {
            std::swap(a, b);
            std::swap(equality.first.first, equality.first.second);
            equality.second = -equality.second;
        }

        joinSegmentGroups(b, a);
        equality = {{a, b}, solution[a] - solution[b]};
        processed[a] = 1;
    }

    // Process groups starting with smallest
    std::priority_queue<
        std::pair<int, int>,
        std::vector<std::pair<int, int>>,
        std::greater<std::pair<int, int>>
    > queue;

    for (size_t i = 0; i < n; i++) {
        if (!processed[i]) {
            queue.push({edgeCount[i], i});
        }
    }

    while (!queue.empty()) {
        int g = queue.top().second;
        int size = queue.top().first;
        queue.pop();

        if ((size_t)size != edgeCount[g] || processed[g]) {
            continue;
        }

        int direction = objectiveFunction[g];
        if (direction == 0) {
            continue;
        }

        // Find first constraint hit by changing variable in desired direction
        int limitingGroup = -1;
        int smallestMove = 0;

        if (direction < 0) {
            smallestMove = INT_MAX;
            for (auto it = edgePool.head(edges[g]); it; ++it) {
                auto& inequality = inequalities[*it];
                if (g == inequality.first.second) {
                    continue;
                }

                int other = inequality.first.second;
                if (getGroup(other) == g) {
                    continue;
                }

                int move = solution[other] + inequality.second - solution[g];
                if (move < smallestMove) {
                    smallestMove = move;
                    limitingGroup = other;
                }
            }
        } else {
            smallestMove = -solution[g] - groupRelativeMin[g]; // Keep all variables >= 0
            for (auto it = edgePool.head(edges[g]); it; ++it) {
                auto& inequality = inequalities[*it];
                if (g == inequality.first.first) {
                    continue;
                }

                int other = inequality.first.first;
                if (getGroup(other) == g) {
                    continue;
                }

                int move = solution[other] - inequality.second - solution[g];
                if (move > smallestMove) {
                    smallestMove = move;
                    limitingGroup = other;
                }
            }
        }

        assert(smallestMove != INT_MAX);
        if (smallestMove == INT_MAX) {
            // Unbound variable - linear program wasn't set up correctly
            // Don't change it instead of stretching graph to infinity
            smallestMove = 0;
        }

        solution[g] += smallestMove;
        if (smallestMove == 0 && !stickWhenNotMoving) {
            continue;
        }

        processed[g] = 1;
        if (limitingGroup != -1) {
            joinSegmentGroups(limitingGroup, g);
            if (!processed[limitingGroup]) {
                queue.push({edgeCount[limitingGroup], limitingGroup});
            }
            equalities.push_back({{g, limitingGroup}, solution[g] - solution[limitingGroup]});
        }
    }

    // Update solutions based on equalities
    for (auto it = equalities.rbegin(); it != equalities.rend(); ++it) {
        solution[it->first.first] = solution[it->first.second] + it->second;
    }
}

void optimizeLinearProgram(size_t n,
                          const std::vector<int>& objectiveFunction,
                          std::vector<Constraint> inequalities,
                          const std::vector<Constraint>& equalities,
                          std::vector<int>& solution) {
    // Remove redundant inequalities with same variable pairs
    std::sort(inequalities.begin(), inequalities.end());
    auto uniqueEnd = std::unique(inequalities.begin(), inequalities.end(),
                               [](const Constraint& a, const Constraint& b) {
                                   return a.first == b.first;
                               });
    inequalities.erase(uniqueEnd, inequalities.end());

    // Run single optimization pass
    static const int ITERATIONS = 1;
    for (int i = 0; i < ITERATIONS; i++) {
        optimizeLinearProgramPass(n, objectiveFunction, inequalities, equalities, solution, true);
    }
}

Constraint createInequality(size_t a, int posA, size_t b, int posB,
                           int minSpacing, const std::vector<int>& positions) {
    minSpacing = std::min(minSpacing, posB - posA);
    return {{a, b}, posB - positions[b] - (posA - positions[a]) - minSpacing};
}

void createInequalitiesFromSegments(std::vector<Segment> segments,
                                   const std::vector<int>& positions,
                                   const std::vector<size_t>& variableGroup,
                                   int blockCount,
                                   int blockSpacing,
                                   int segmentSpacing,
                                   std::vector<Constraint>& inequalities) {
    // Map used as binary search tree: y_position -> segment{variableId, x_position}
    // Used to track which segment was last seen at each y position
    std::map<int, std::pair<int, int>> lastSegments;
    lastSegments[-1] = {-1, -1}; // Sentinel value

    // Sort segments by x coordinate
    std::sort(segments.begin(), segments.end(),
              [](const Segment& a, const Segment& b) { return a.x < b.x; });

    for (const auto& segment : segments) {
        // Find last segment that starts before or at current segment's y0
        auto startPos = lastSegments.lower_bound(segment.y0);
        --startPos; // Safe because map has sentinel at -1

        auto lastSegment = startPos->second;
        auto it = startPos;

        // Process all segments that overlap with current segment's y range
        while (it != lastSegments.end() && it->first <= segment.y1) {
            int prevSegmentVariable = it->second.first;
            int prevSegmentPos = it->second.second;

            if (prevSegmentVariable != -1) {
                int minSpacing = segmentSpacing;

                // Handle special spacing cases
                if (prevSegmentVariable < blockCount && segment.variableId < blockCount) {
                    // Skip inequality between two sides of same block
                    if (prevSegmentVariable == segment.variableId) {
                        ++it;
                        continue;
                    }
                    minSpacing = blockSpacing;
                }
                else if (variableGroup[prevSegmentVariable] == variableGroup[segment.variableId]) {
                    minSpacing = 0; // No spacing needed between segments of same edge
                }

                inequalities.push_back(createInequality(prevSegmentVariable, prevSegmentPos,
                                                      segment.variableId, segment.x,
                                                      minSpacing, positions));
            }

            lastSegment = it->second;
            ++it;
        }

        // Update lastSegments map
        if (startPos->first < segment.y0) {
            startPos++;
        }
        lastSegments.erase(startPos, it); // Remove segments covered by current one

        // Add current segment and remaining part of partially covered segment
        lastSegments[segment.y0] = {segment.variableId, segment.x};
        lastSegments[segment.y1] = lastSegment;
    }
}
