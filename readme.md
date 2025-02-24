<p align="center">
  <img src="https://github.com/19h/dag-flowchart-cpp/raw/master/graph.svg" width="600px"/>
</p>

<hr/>

<h5 align="center">
A lightweight, dependency-free graph layout algorithm for C++ projects.
<br/>
Drop-in solution for visualizing directed graphs with optimized node placement and edge routing.
<br/>
<br/>
Creates clean, layered visualizations using advanced layout techniques
<br/>
to produce readable diagrams of control flow and call graphs.
</h5>

<hr/>

```bash
$ clang++ -std=c++17 -o graph -Ofast graph.cpp
```

### tl;dr

- Advanced graph layout algorithm producing layered visualizations
- Zero dependencies - ready to use in any C++ project
- Uses DFS, topological sorting and tree placement for layouts
- Platform-independent from embedded systems to servers
- Handles node positioning and edge routing with minimal memory ✨

### tl;dr

A fast, focused graph layout engine that transforms complex graphs into clear, optimized visualizations. Built with proven algorithms and zero dependencies.

## Usage

This library implements a grid‐based graph layout engine for advanced visualization of directed graphs (e.g. control flow, call graphs, state machines, dependency graphs). It consists of several core components:

- **BinaryTrees.h**  
  Contains auxiliary balanced tree structures used internally for efficient event and range queries during layout.

- **LinkedListPool.h**  
  Provides a custom memory pool for managing linked lists, which are employed to maintain subtree shape profiles and dynamic constraint lists during layout computation.

- **GraphLayout (in graph.cpp)**  
  Defines the abstract base class `GraphLayout` and concrete implementations such as `GraphGridLayout` (with configurable modes: Narrow, Medium, or Wide) and the adapter `GraphHorizontalAdapter` (for swapping axis orientations). These classes encapsulate algorithms for cycle removal, topological sorting, subtree placement, multi-segment edge routing, and constraint-based layout optimization.

> **Note:** The example `main()` function in *graph.cpp* is provided solely as a demonstration. **Remove or exclude this function** when integrating the library into your project.

### Integration Steps

1. **Include Required Files**

   Add the following files to your project:
   - `BinaryTrees.h`
   - `LinkedListPool.h`
   - `graph.cpp` (with the demonstration `main()` removed)

2. **Compile the Library**

   Compile these files together with your project. The library makes extensive use of STL containers and modern C++ idioms. Configure your build system (e.g., CMake, Makefiles) accordingly.

### Populating the Graph Structure

The library operates on a graph defined as an unordered map from unique identifiers (typically `uint64_t`) to graph nodes (instances of `GraphLayout::GraphBlock` or a user-extended type such as `GraphBlockWithText`). Each node (or block) encapsulates:
  
- **Geometry:** Initial dimensions (`width` and `height`), later refined to pixel coordinates.
- **Connectivity:** A unique identifier (`entry`) and a vector of outgoing edges (`GraphEdge`), each containing routing information (polylines and arrow directions).
- **Optional Metadata:** In extended versions (e.g., `GraphBlockWithText`), text labels or other user-defined data may be included.

#### Creative Approaches to Graph Population

Since the library imposes no specific input format, you have complete flexibility in how you create and populate the graph:

- **Programmatic Construction:**  
  Build your graph directly in C++ by instantiating nodes and connecting them with edges. For example, if you are analyzing a control flow graph from a compiler or debugger, create a node for each basic block and add edges for jumps or function calls.

  ```cpp
  #include "BinaryTrees.h"
  #include "LinkedListPool.h"
  #include "graph.cpp"  // Ensure the main() function is removed

  // Define your graph type (you can extend GraphBlock to include custom metadata)
  using Graph = std::unordered_map<uint64_t, GraphLayout::GraphBlock>;

  Graph graph;

  // Create nodes programmatically.
  graph[1] = GraphLayout::GraphBlock{0, 0, 120, 50, 1};
  graph[2] = GraphLayout::GraphBlock{0, 0, 100, 50, 2};
  // Connect nodes with edges.
  graph[1].edges.push_back(GraphLayout::GraphEdge(2));
  // Continue adding nodes and edges as needed...
  ```

- **Dynamic Graph Generation:**  
  Integrate the layout engine with runtime data. For instance, in a GUI-based debugger or performance analyzer, dynamically build graphs representing execution flows or inter-module dependencies. After processing the underlying data, convert it into the required graph structure and invoke the layout engine to compute node positions and edge routes.

- **Data-Driven Layout:**  
  If your application already represents graphs using another format (e.g., JSON, XML, or a custom binary format), write an adapter that converts your data into the library’s graph structure. This decouples the layout logic from input parsing, allowing you to integrate the layout engine regardless of your source data format.

- **Hybrid Visualizations:**  
  Use the library as part of a larger visualization framework. For example, combine the computed layout with custom rendering code (e.g., generating SVG or using a GUI toolkit) to overlay additional information, animate transitions, or support interactive exploration.

### Choosing and Configuring a Layout

The library supports various layout modes and orientations. For vertical layouts with balanced subtree spacing, instantiate `GraphGridLayoutWithText` (a subclass of `GraphGridLayout`). To obtain a horizontal layout, wrap your instance in a `GraphHorizontalAdapter`:

```cpp
// Instantiate the base layout with your preferred mode.
std::unique_ptr<GraphLayout> baseLayout =
    std::make_unique<GraphGridLayoutWithText>(
        GraphGridLayoutWithText::LayoutType::Medium,
    );

// Optionally switch to horizontal orientation.
bool useHorizontalLayout = false;  // Set true for horizontal rendering.
std::unique_ptr<GraphLayout> layout =
    useHorizontalLayout
        ? std::make_unique<GraphHorizontalAdapter>(std::move(baseLayout))
        : std::move(baseLayout);
```

Customize layout parameters (spacing, edge routing preferences, etc.) via the `LayoutConfig` structure:

```cpp
GraphLayout::LayoutConfig config;

config.blockVerticalSpacing = 40;
config.blockHorizontalSpacing = 20;

config.edgeVerticalSpacing = 10;
config.edgeHorizontalSpacing = 10;

layout->setLayoutConfig(config);
```

### Computing and Using the Layout

Invoke the `CalculateLayout` method to compute pixel coordinates for nodes and routing paths for edges:

```cpp
int layoutWidth, layoutHeight;

layout->CalculateLayout(
    reinterpret_cast<GraphLayout::Graph&>(graph),
    /* entry id */ 1,
    layoutWidth,
    layoutHeight,
);
```

After layout computation:
- Each node’s `x` and `y` fields are updated to its position.
- Each edge’s `polyline` vector holds the routed points for drawing smooth connections.
- The overall dimensions (`layoutWidth` and `layoutHeight`) can be used to size your rendering canvas.

### Rendering and Further Processing

With the computed layout:
- **Rendering:**  
  Render nodes (as rectangles, circles, etc.) and edges (as polylines or spline curves) in your GUI or output format. For instance, you can generate SVG elements or use a graphics library to draw directly on a canvas.

- **Interactive Visualization:**  
  Integrate the layout engine into interactive tools where users can explore, zoom, and pan large graphs. The algorithm’s flexibility in node and edge spacing makes it well-suited for dynamic updates.

- **Custom Extensions:**  
  Extend `GraphBlock` or create your own subclass to add domain-specific metadata (e.g., profiling data, annotations, or debug information) that can be rendered alongside the layout.

# Graph Grid Layout Algorithm

This algorithm produces a layered (grid) layout for directed graphs (for example, control‐flow or call graphs). It takes as input a graph (provided in DOT format) and assigns grid coordinates to nodes (“blocks”) along with carefully routed edge paths. The final layout is output as an SVG image.

The implementation is broken down into several components:

• A preprocessing phase that “linearizes” the graph using DFS and topological sorting so that nodes can be assigned rows.  
• A node placement phase in which a spanning tree is extracted and subtrees are merged side‐by side using a linked–list data structure to track subtree “shapes.”  
• An edge routing phase that computes “main columns” (vertical channels) for each edge using range–minimum queries and a sweeping algorithm.  
• A layout compaction and optimization phase that runs a custom linear–programming solver to adjust spacing between nodes and edge segments.  
• Finally, conversion from grid coordinates to pixel positions yields the final SVG output.

The source code is organized over several files. For example, “BinaryTrees.h” contains segment trees used to query and update ranges quickly, “LinkedListPool.h” implements pooled linked lists to merge subtree outlines, and “graph.cpp” contains the bulk of the layout algorithm (together with parsing input and producing SVG output).

---

## 1. Preprocessing and Topological Sorting

The algorithm begins by converting the input graph (which may contain cycles) into a directed acyclic graph (DAG) in order to assign layers:

• A non‐recursive DFS is performed (see the function `topoSort`) starting from the entry node.  
• Each node is “visited” and the DFS tracks both:
 – Tree edges (used later for spanning tree extraction), and  
 – Back or cross edges (which remain as additional “DAG” edges once cycles are removed).

The DFS produces a topological ordering (an array of node IDs) that is later used to assign each node a row index:
  row[node] = max { row[predecessor] } + 1.

---

## 2. Row Assignment and Spanning Tree Selection

Using the reverse topological order the algorithm assigns each node a layer (row):

• For every node in reverse order, its row is set to one greater than the maximum row among its direct predecessors.  
• Next, from the available DAG edges the algorithm “chooses” a subset to form a spanning tree.  
 – For each node (from top to bottom) the algorithm greedily selects—as a child—a candidate of the next row (if that child is not already assigned a parent).  
 – In doing so the resulting structure is a rooted tree where the parent is above (i.e. in a lower–numbered row) its children.

This spanning tree makes the node placement problem similar to a tree–drawing problem since subtrees will be placed side by side.

---

## 3. Subtree Placement and Node Positioning

Once the rows are fixed, the algorithm computes the horizontal (column) positions of nodes on the grid. This is done in a bottom–up manner:

• **Leaf nodes:**  
 – For each leaf node, a “shape” is created (represented as a linked list stored in a memory pool via `LinkedListPool`) to capture its horizontal extent.  
 – Their initial column positions, left/right bounds, and row counts are set.

• **Internal nodes:**  
 – For each internal node, the algorithm “merges” the shapes of its child subtrees. During the merge, it computes the relative offsets needed to pack the subtrees as tightly as possible.  
 – Two configuration flags control the placement:
  – One flag decides whether the parent should be centered between its leftmost and rightmost child (yielding a wider layout) or be placed at the average position of the children.  
  – Another flag (tight subtree placement) allows subtrees to be packed more closely by using their exact profiles rather than bounding boxes.

• After the merging procedure, the algorithm makes one final pass to convert all relative column positions into absolute pixel coordinates. It calculates the maximum width required per column and the maximum height per row so that the grid conversion is consistent.

---

## 4. Edge Routing

Once nodes are laid out in the grid, the next task is to route the edges between them so that curves do not overlap nodes and the paths are easy to follow.

### a. Main Column Assignment

• For each edge—defined by its source and target nodes—a vertical “channel” (the “main column”) is chosen for the edge’s primary vertical segment.  
• The algorithm uses a sweep–line approach along with a segment tree (implemented by `PointSetMinTree`) to determine the closest available columns that are not “blocked” by a node.
 – First, it tries the source node’s column; if blocked it tries the target, and if both are blocked it finds the nearest available column using range–queries.

### b. Rough Edge Routing

• Once the main column is chosen, several routing points are inserted to form the edge’s path.  
• Typically, an edge path consists of up to five segments:  
 1. A short vertical segment from the source node;  
 2. A horizontal segment moving toward the main column;  
 3. A long vertical segment along the main column between source and target rows;  
 4. A horizontal segment moving toward the target node; and  
 5. A short vertical segment connecting finally to the target.
• Additional data (such as segment “kind” and spacing overrides) is stored with each routing point to guide later fine–tuning.

### c. Edge Offset Calculation

• Multiple edges passing through the same column may “collide.” To avoid this the algorithm computes offsets for each edge segment.  
• It does so by collecting all vertical (and horizontal) segments, sorting them by their positions and priorities and then assigning offsets in a greedy fashion.
• Finally, these offsets are “centered” within the available column width so that the edges look evenly spaced.

---

## 5. Layout Optimization via Linear Programming

After initial placements are computed, the algorithm performs a layout “compaction” step. This is achieved through a lightweight linear programming (LP) solver that adjusts node and edge positions while satisfying spacing constraints.

• **Constraints and Objectives:**  
 – Inequality constraints ensure that two adjacent nodes or edge segments are kept at a minimum spacing apart.  
 – Equality constraints are added in cases where a node must remain centered (for example, when a branch merges) or when a block edge must be tied to a node’s boundary.
 – An objective function (a weighted sum of variable positions) is minimized to “pull” nodes as close together as possible without violating spacing constraints.

• **Solver Implementation:**  
 – The LP solver is implemented in the functions `optimizeLinearProgramPass` and `optimizeLinearProgram`.  
 – It groups variables (each representing a node or an edge segment position) and iteratively adjusts their positions until the constraints are met.  
 – Although the solver does not guarantee a mathematically optimal solution, it produces a “good‐enough” result in practical terms.

• The refinement is performed separately for vertical and horizontal positions. After each optimization pass the computed positions are copied back to the graph data structure.

---

## 6. Conversion to Final Pixel Coordinates and Output

With all node and edge positions computed and optionally compressed via the LP solver:

• Cumulative offsets for columns and rows are calculated to obtain absolute (pixel) positions.  
• A final “cropping” step trims any surplus white space around the graph.
• Edge endpoints are “connected” so that the bottom of a source block meets the start of the corresponding edge and the top of the target block receives the edge.
• Finally, the algorithm writes out the graph as an SVG file. Nodes are drawn as rectangles with labels, and edges are rendered as polyline curves with arrowheads.

---

## 7. Optional Horizontal Layout

Although the default is a vertical (top–to–bottom) layout, a horizontal layout is supported via a simple adapter:

• The `GraphHorizontalAdapter` class wraps any layout object and simply swaps x– and y–coordinates before and after layout calculation.  
• This gives the effect of a left-to–right arranged graph with minimal changes to the core algorithm.

---

## Summary

To summarize, the algorithm works as follows:

1. **Preprocessing:**  
 – Run a DFS to topologically sort and remove cycles.  
 – Assign row numbers based on predecessor rows.

2. **Node Placement:**  
 – Extract a spanning tree and merge child subtrees via linked–list structures.  
 – Compute each node’s horizontal (column) position and convert grid positions to pixels.

3. **Edge Routing:**  
 – For each edge, choose a “main column” using range–minimum queries.  
 – Insert routing points to form a polyline path and compute spacing offsets to avoid collisions.

4. **Layout Optimization:**  
 – Set up a linear programming system with inequality and equality constraints that assures minimal spacing.  
 – Run a series of optimization passes to compact the layout.

5. **Final Conversion and SVG Output:**  
 – Calculate cumulative offsets for columns and rows.  
 – Adjust node and edge positions and write the final layout to an SVG file.

The entire process is implemented in C++ and uses custom memory pools, segment trees, and a greedy LP–solver to achieve a visually pleasing layout.


---

## License

Hark! Behold, another decree of software licensure doth appear before thee!

Verily, thou shalt spend but little time perusing these scrolls, for thy hours are better employed unraveling the mystic riddle of thy once-perfect template metaprogram—now accursed with a staggering forty-seven thousand, two hundred and ninety-three lines of compiler lamentations.

**By rare providence, this code hath been cast into the public commons.**

Wherefore?

Forsooth, life is fleeting, and no good knight hath the leisure for fifty pages of legal decrees nor to tarry awaiting the arbitrations of counsel for the mere use of std::string (as though any soul were ejected for wielding char*).

Copy it, alter it, vend it, or compile it into a maze of bewildering template errors—we care not one jot! Whether thou choosest static linking or dynamic linking, or even that fearsome SFINAE beast that maketh boost::spirit seem as tame as a minstrel’s tune, venture forth and be unrestrained! We vow not to hurl exceptions in its stead (unlike the vexations of thy production code).

Feel at liberty to std::move() this code to any realm thou desirest—yet, in light of thine infamous skirmishes with dangling references, perchance remain true to the art of copying. And nay, enshrouding all within shared_ptr shall not remedy the ills of thy design.

**To ye realms most besotted with copyrights: We hereby renounce all rights to this code, casting it into the public domain with a zeal rivaling that of C++ artisans unraveling fresh incantations to invoke std::function.**

static_assert(true, "Yes, we really mean it"); // Mayst thou discover this gem scarcely within thy build errors!

Lo, the software is bestowed "as is"—a stark and unadorned surprise! No warranties, no pledges of guarantee, not even a modest static_assert attend it.

Should it, by some uncanny alchemy, gain sentience and begin crafting code mightier than thine own, that grievous burden shall be thine alone.

Press onward with thy templates as a battle-hardened artificer well-versed in the trials of code. May thy builds compile without flaw, and thy linker errors be but murmurs lost in the annals of thy log.