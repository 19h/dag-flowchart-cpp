#ifndef LAYOUT_OPTIMIZER_H
#define LAYOUT_OPTIMIZER_H

#include "core/GridBlock.h"
#include "optimization/LinearProgramming.h"
#include <unordered_map>
#include <vector>

/**
 * @class LayoutOptimizer
 * @brief Optimizes graph layout using linear programming to minimize edge lengths.
 *
 * This class performs post-processing optimization on an already-positioned graph
 * to improve visual quality by:
 * - Minimizing total edge length (straightening edges where possible)
 * - Maintaining minimum spacing between blocks and edges
 * - Centering parent blocks between their children
 * - Aligning merge points with their parent blocks
 *
 * The optimization is performed in two phases:
 * 1. **Vertical (Y-axis)**: Optimize block Y positions and horizontal edge segment Y positions
 * 2. **Horizontal (X-axis)**: Optimize block X positions and vertical edge segment X positions
 *
 * @par Linear Programming Formulation:
 * Each block and edge segment becomes a variable. The problem is:
 * - **Minimize**: Sum of edge segment offsets from their neighbors (reduces edge length)
 * - **Subject to**:
 *   - Spacing constraints: Adjacent elements must maintain minimum spacing
 *   - Edge-block constraints: Edge endpoints must stay attached to their blocks
 *   - Centering constraints: Parents centered between children (soft constraint)
 *
 * @par Segment-Based Collision Detection:
 * To prevent overlaps, blocks and edges are represented as vertical/horizontal
 * segments in a sweep-line algorithm. The algorithm generates inequality constraints
 * for all pairs of segments that could potentially overlap.
 *
 * @par Usage:
 * @code
 *   LayoutOptimizer optimizer(config);
 *   optimizer.optimize(state);
 *   // state.blocks now contains optimized positions
 * @endcode
 *
 * @see LinearProgramming.h for the LP solver implementation
 * @see createInequalitiesFromSegments() for segment-based constraint generation
 */
class LayoutOptimizer {
public:
    /**
     * @brief Constructs a LayoutOptimizer with layout configuration.
     * @param config Layout configuration containing spacing parameters
     */
    explicit LayoutOptimizer(const GraphLayout::LayoutConfig& config) : config_(config) {}

    /**
     * @brief Main entry point: optimizes block and edge positions.
     *
     * Performs two-phase optimization:
     * 1. optimizeVertical() - Optimize Y coordinates
     * 2. optimizeHorizontal() - Optimize X coordinates
     *
     * @param[in,out] state Layout state with positioned blocks and edges;
     *        positions are modified in-place
     *
     * @pre Blocks must have valid x, y, width, height
     * @pre Edges must have valid polyline coordinates
     * @post Block and edge positions are optimized while maintaining constraints
     */
    void optimize(LayoutState& state) const;

private:
    const GraphLayout::LayoutConfig& config_;

    /**
     * @brief Internal context for tracking optimization state.
     *
     * This struct encapsulates all temporary data needed during optimization,
     * including variable mappings, constraints, and the LP solution.
     */
    struct OptimizationContext {
        /// Maps block IDs to variable indices (0..N-1 for N blocks)
        std::unordered_map<uint64_t, int> blockMapping;
        
        /// Group ID for each variable (blocks share group 0..N-1, edges get unique groups)
        std::vector<size_t> variableGroups;
        
        /// Objective function coefficients (minimize sum of objective[i] * solution[i])
        std::vector<int> objective;
        
        /// Inequality constraints: a - b >= constant
        std::vector<Constraint> inequalities;
        
        /// Equality constraints: a - b = constant
        std::vector<Constraint> equalities;
        
        /// Current solution vector (feasible values during construction, optimal after solving)
        std::vector<int> solution;
        
        /// Segments for sweep-line collision detection
        std::vector<Segment> segments;
        
        /// Next available variable index (beyond block variables)
        size_t variableIndex = 0;
        
        /// Next available edge group index
        size_t edgeIndex = 0;

        /**
         * @brief Resets context for a new optimization phase.
         * @param numBlocks Number of blocks (determines initial variable count)
         */
        void reset(size_t numBlocks);

        /**
         * @brief Sets the initial feasible value for a variable.
         * @param variable Variable index
         * @param value Initial value (used as starting point for LP solver)
         */
        void setFeasibleValue(size_t variable, int value);

        /**
         * @brief Adds a term to the objective to minimize distance between two variables.
         *
         * If posA < posB, adds coefficient +1 to b and -1 to a (minimizes b - a).
         * Otherwise, adds +1 to a and -1 to b (minimizes a - b).
         *
         * @param a First variable index
         * @param posA Current position of first variable
         * @param b Second variable index
         * @param posB Current position of second variable
         */
        void addToObjective(size_t a, int posA, size_t b, int posB);
    };

    /**
     * @brief Creates block-to-variable mapping.
     * @param state Layout state with blocks
     * @param ctx Context to populate with mapping
     */
    void initializeBlockMapping(LayoutState& state, OptimizationContext& ctx) const;

    /**
     * @brief Optimizes vertical (Y) coordinates of blocks and horizontal edge segments.
     *
     * Creates variables for:
     * - Block Y positions (variables 0..N-1)
     * - Horizontal edge segment Y positions (variables N..)
     *
     * Constraints:
     * - Parent blocks must be above children (with spacing)
     * - Segments must not overlap (generated by sweep-line)
     *
     * @param state Layout state to optimize
     * @param ctx Optimization context
     */
    void optimizeVertical(LayoutState& state, OptimizationContext& ctx) const;

    /**
     * @brief Optimizes horizontal (X) coordinates of blocks and vertical edge segments.
     *
     * Creates variables for:
     * - Block X positions (variables 0..N-1)
     * - Vertical edge segment X positions (variables N..)
     *
     * Constraints:
     * - Edge endpoints must stay attached to their blocks (equality)
     * - Segments must not overlap (generated by sweep-line)
     * - Parents with 2 children: centered between children (if possible)
     *
     * @param state Layout state to optimize
     * @param ctx Optimization context
     */
    void optimizeHorizontal(LayoutState& state, OptimizationContext& ctx) const;

    /**
     * @brief Applies optimized positions from solution to state.
     *
     * Copies solution values back to block positions and edge polylines.
     *
     * @param state Layout state to update
     * @param ctx Context containing solution
     * @param horizontal If true, applying Y positions; if false, applying X positions
     */
    void applyPositions(LayoutState& state, const OptimizationContext& ctx, bool horizontal) const;
};

#endif // LAYOUT_OPTIMIZER_H
