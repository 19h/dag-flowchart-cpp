#ifndef LINEAR_PROGRAMMING_H
#define LINEAR_PROGRAMMING_H

#include <cstddef>
#include <utility>
#include <vector>

#include "core/EdgeSegment.h"

/// Constraint: either equality or inequality x_i <= x_j + c
using Constraint = std::pair<std::pair<int, int>, int>;

/**
 * @brief Linear programming solver for layout optimization.
 * 
 * Optimizes a linear objective function subject to constraints.
 * Note: Does not guarantee finding the optimal solution, but finds a 
 * feasible solution that improves the objective value.
 *
 * @param n Number of variables in the linear program
 * @param objectiveFunction Coefficients ci for minimizing sum(ci * xi)
 * @param inequalities Constraints of form x[ei] - x[fi] <= bi
 * @param equalities Constraints of form x[ei] - x[fi] = bi
 * @param solution Input: Initial feasible solution, Output: Optimized solution
 */
void optimizeLinearProgram(size_t n,
                          const std::vector<int>& objectiveFunction,
                          std::vector<Constraint> inequalities,
                          const std::vector<Constraint>& equalities,
                          std::vector<int>& solution);

/**
 * @brief Creates an inequality constraint between two variables.
 *
 * @param a First variable index
 * @param posA Position of first variable
 * @param b Second variable index
 * @param posB Position of second variable
 * @param minSpacing Minimum required spacing between variables
 * @param positions Vector of current positions
 * @return Constraint representing a - b <= value
 */
Constraint createInequality(size_t a, int posA, size_t b, int posB, 
                           int minSpacing, const std::vector<int>& positions);

/**
 * @brief Create inequality constraints from segments which preserves their relative order.
 *
 * @param segments List of edge segments and block sides
 * @param positions Initial element positions before optimization
 * @param blockCount Number of variables representing blocks
 * @param variableGroup Used to check if segments are part of the same edge
 * @param blockSpacing Minimal spacing between blocks
 * @param segmentSpacing Minimal spacing between two edge segments
 * @param inequalities Output vector for resulting inequalities
 */
void createInequalitiesFromSegments(std::vector<Segment> segments,
                                   const std::vector<int>& positions,
                                   const std::vector<size_t>& variableGroup,
                                   int blockCount,
                                   int blockSpacing,
                                   int segmentSpacing,
                                   std::vector<Constraint>& inequalities);

#endif // LINEAR_PROGRAMMING_H
