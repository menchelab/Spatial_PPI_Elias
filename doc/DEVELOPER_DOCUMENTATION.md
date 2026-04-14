# Cell Layout - Developer Documentation

This document explains the code of my cell_layout program.
This developer documentation focuses on:

- how modules connect,
- what the main objects do,
- how configuration is parsed and validated,
- how constraints and the layout are implemented,
- how to safely extend the system.

---

## 1. module overview

Core modules:

- cli.py - entry point -> argument parsing
- config.py - config schema parsing and constraints
- io_helpers.py - data loading, output writer and helper utilities
- model.py - CellModel and Location objects and the Plotly visualization
- constraints.py - constraints and constrained 3D spring layout implementation

Data files:

- config.json - example configuration with default params based on literature research of 13 compartments
- positions_per_location.tsv - output file produced
- PPI.tsv - protein-protein interaction network
- Localization_table.tsv - Locationnames in config file and localization data must match!!
---

## 2. General Flow of the program

1. cli.run()

  - parses CLI args (--ppi, --localizations, --config, --outdir, optional --only/--exclude, --seed/--iterations, plot flags)
  - creates outdir
  - loads config via config.load_config(config_path)
  - resolves seed and iterations (CLI overrides config defaults)
  - parses --only / --exclude filters
  - loads PPI graph via io_helpers.load_ppi_graph(ppi_path, only, exclude)
  - loads localizations via io_helpers.load_localizations(localizations_path)
  - initializes node attribute via io_helpers.init_locations_attribute(G) (sets G.nodes[n]["locations"] = {})
  - builds CellModel(G)
  - for each LocationSpec in config
  - validates location exists in localization table
  - applies only/exclude filters to localization-derived node list
  - collects node IDs (ignores NaNs)
   -builds constraint via constraints.build_constraint(LocationSpec, ...)
  - creates model.Location(...) and adds it to the CellModel
  - checks if any locations were selected
    - if none: exits early (return 0)
    - if yes: calls cell.assign_all_positions(...)

  - writes positions via io_helpers.write_positions(cell, outdir) (produces positions_per_location.tsv)
  - prints summary via cell.summary()
  - if plotting requested: calls cell.plot_all_locations_3d(...) (Plotly; optional edges/boundaries toggles)

2. CellModel.assign_all_positions()

  - iterates through all Location objects in the model
  - calls Location.assign_positions_to_nodes(...) for each location

3. Location.assign_positions_to_nodes()

  - builds a location-specific subgraph H via Location.make_subgraph(G, node_ids)
  - filters nodes by min_degree (drops low-degree nodes before layout)
  - runs constrained layout via model.spring_layout_3d_constrained(H, constraint, iterations, seed, center, repulsion_strength)
  - stores resulting pos3d in the Location
  - writes back per-node coordinates into G.nodes[node]["locations"][location_name]

4. io_helpers.write_positions(cell, outdir)

  - iterates all nodes and all stored per-location positions
  - writes a tab-separated table positions_per_location.tsv into outdir

5. CellModel.plot_all_locations_3d()

  - builds a Plotly 3D visualization using stored positions
  - optionally renders compartment boundaries and/or graph edges based on CLI plot flags

---

## 3. Data model

### 3.1 Global PPI graph (whole_cell_graph)
Type: networkx.Graph

- Nodes are protein identifiers (strings).
- Each node receives a locations attribute (dict) via init_locations_attribute():
  - G.nodes[node]["locations"] is a mapping: location_name -> (x,y,z)
- This supports proteins assigned to multiple compartments without duplicating nodes.

### 3.2 Location object (model.Location)
A location encapsulates:

- name: compartment label
- node_ids: list of proteins assigned to this location (from localization table)
- min_degree: used for subgraph filtering
- center: translation vector in global coordinates
- constraint: instance of a constraint class (see section 5)
- repulsion_strength, seed, iterations: layout parameters
- pos3d: dict node_id -> (x,y,z) for this compartment

Key methods:

- make_subgraph(whole_graph)
- assign_positions_to_nodes(whole_graph)

Degree filtering:
make_subgraph() removes nodes with degree ≤ min_degree:

So with min_degree=1 you keep nodes with degree ≥ 2.

### 3.3 Cell model (model.CellModel)
A lightweight coordinator:

- holds whole_cell_graph
- holds locations: Dict[str, Location]

Key methods:

- add_location(location)
- assign_all_positions()
- summary()
- plot_all_locations_3d()

---

## 4. Configuration handling (config.py)

### 4.1 Schema objects
- LocationSpec: parsed representation of a single location config
- AppConfig: parsed representation of the whole config

load_config(path):
- reads JSON
- validates each constraint.type is known (key in CONSTRAINTS)
- returns AppConfig(default_params=..., locations=[LocationSpec,...])

### 4.2 Constraints config
CONSTRAINTS maps type names in JSON to Python classes:

CONSTRAINTS = {
    "SphereConstraint": SphereConstraint,
    "ShellConstraint": ShellConstraint,
    "CylinderConstraint": CylinderConstraint,
    "EllipsoidConstraint": EllipsoidConstraint,
    "EllipsoidShellConstraint": EllipsoidShellConstraint,
}

constraint.params must match __init__ parameter names of the class!!

---

## 5. Constraints (constraints.py)

Every constraint class implements a consistent interface:

- sample(n, rng) -> np.ndarray shape (n,3)
  - used for initial positions
- forces(q) -> np.ndarray shape (n,3)
  - computes soft-wall restoring forces
- boundary_traces() -> List[plotly.graph_objects.Trace] (optional)
  - used for visualization boundaries

### 5.1 Soft-wall
Constraints do not do hard projection. Instead they add forces that grow with boundary violation:

- Sphere/Ellipsoid: push inward when outside
- Shell constraints: push inward when outside outer boundary; push outward when inside inner boundary

### 5.2 Sampling
A few sampling choices are implementation-specific and worth remembering:

- SphereConstraint.sample uses rng.random(n) ** (1/3) to approximate uniform volume density.
- EllipsoidConstraint.sample scales random directions by a uniform radius in [0,1] (not cube-root), which biases samples toward the center.
- CylinderConstraint.sample initially places nodes on the radial boundary (r = radius) and uses z ~ U(0, height). Because the layout recenters each iteration (default), the cloud becomes centered around z=0 over time.

If you care about strict uniform sampling, adjust these samplers accordingly.

---

## 6. Constrained 3D spring layout (spring_layout_3d_constrained)

This project’s core numeric routine lives in constraints.py:

- it runs a **Fruchterman-Reingold-style** force-directed layout in **3D**
- it adds **soft constraint forces** from a shape object (SphereConstraint, ShellConstraint, …)
- it optionally **recenters** the point cloud each iteration to avoid drift
- it finally **translates** the layout by a per-compartment center so compartments can be placed side-by-side

### 6.1 Signature and return type

spring_layout_3d_constrained(
    G,
    constraint,
    iterations: int,
    seed: int,
    k: float = None,
    threshold: float = 1e-4,
    recenter_each_iter: bool = True,
    center=(0.0, 0.0, 0.0),
    repulsion_strength: float = 1.0,
) -> Dict[node_id, np.ndarray(3,)]


Return value: a dict mapping each node in G to a length-3 NumPy vector (x, y, z).

### 6.2 Inputs and invariants

Graph (G)
- Expected to be the per-location subgraph H produced by Location.make_subgraph().
- Node ordering is nodes = list(G). This order affects determinism (see 6.9).
- Edge weights are ignored (weight=None), so adjacency is effectively binary.

Constraint object (constraint)
Must implement:
- sample(n, rng) -> (n,3) float array  
  Provides initial 3D positions (in the constraint’s local coordinate frame).
- forces(q) -> (n,3) float array  
  Computes additional force vectors (soft wall behavior) given current positions.

Optional:
- boundary_traces() for visualization; not used by the layout itself.

### 6.3 Dense internal representation

Internally the function converts the graph to dense arrays:

A = nx.to_numpy_array(G, nodelist=nodes, weight=None, dtype=float)  # shape (n,n)
pos = constraint.sample(n, rng)                                     # shape (n,3)

Each iteration then constructs pairwise matrices:

- delta: shape (n, n, 3) (vector from j → i)
- dist:  shape (n, n)    (pairwise distances)

This makes the layout:
- time: O(n² * iterations)
- memory: O(n²) (delta dominates)

### 6.4 Initial spacing parameter k auto-estimation

k acts like the preferred spacing scale of Fruchterman-Reingold.

If k is not provided, it is estimated from the bounding box volume of the sampled initial positions:

extent = np.ptp(pos, axis=0)              # (Lx, Ly, Lz)
bbox_vol = extent[0] * extent[1] * extent[2]
k = (bbox_vol / n) ** (1/3)

- bbox_vol / n approximates the volume available per node
- taking the cube-root gives an average length scale

This ties the force magnitudes to the constraints sampling scale. If a constraints sample() returns points on a very large or very small scale, k will follow.

### 6.5 Force model with repulsion, attraction and constraint

Per iteration, the function computes pairwise offsets:

delta = pos[:, None, :] - pos[None, :, :]         # (n,n,3)
dist  = np.linalg.norm(delta, axis=-1)            # (n,n)
np.clip(dist, 0.01, None, out=dist)               # avoids divide-by-zero

It then forms a scalar coefficient matrix:

coeff = (repulsion_strength * (k*k) / (dist*dist)) - (A * dist / k)


Important: coeff is multiplied by delta (a vector of length dist).
So the actual force magnitudes become:

- repulsion: ||delta|| * (k² / dist²) → k² / dist  
  (matches the classic FR repulsive magnitude)
- attraction (edges only): ||delta|| * (dist / k) → dist² / k  
  (matches the classic FR attractive magnitude)

Then forces are accumulated per node:

disp = np.einsum("ijk,ij->ik", delta, coeff)       # (n,3)
disp += constraint.forces(pos)                     # (n,3)


So the final displacement vector per node is:

- sum of FR-style pairwise forces
- plus a constraint-specific soft-wall correction

#### The role of repulsion_strength
repulsion_strength scales only the repulsive term. Increasing it:
- expands the layout / reduces crowding
- can make constraint wall forces more important (because nodes are pushed outward)

### 6.6 Step limiting via temperature

The routine uses a simple linear cooling schedule:

bounding_box = np.ptp(pos, axis=0)
t  = bounding_box.max() * 0.1
dt = t / (iterations + 1)


Each node’s displacement is then normalized and scaled by t:

disp_length = np.linalg.norm(disp, axis=1)
disp_length = np.clip(disp_length, 1e-12, None)
step = disp * (t / disp_length)[:, None]
pos += step

- each node moves by approximately t each iteration (direction changes, magnitude capped)
- early iterations allow big moves; later iterations shrink moves as t decreases
- very large force magnitudes do not create huge steps (normalization caps them)

This is a good stabilization strategy, but it also means:
- a node with tiny disp still gets scaled up (because of normalization) unless the direction becomes numerically unstable; the 1e-12 clip prevents division-by-zero but can still yield a near-arbitrary direction when disp is extremely small.
- the early-stopping check (6.8) is important to avoid thrashing when the system is near equilibrium.

### 6.7 Recentering and global translation

To avoid the entire point cloud drifting away due to asymmetric forces, the algorithm optionally recenters every iteration:

if recenter_each_iter:
    pos -= pos.mean(axis=0)


This keeps the mean at (0,0,0) in the constraints local frame.

After convergence, it applies the compartment’s translation vector:

pos += np.asarray(center, dtype=float)


So constraints and FR forces operate around the origin, and the compartment is placed into the global cell coordinate system only at the end.

### 6.8 Early stopping condition

The loop exits early if the average overall step becomes small:

if (np.linalg.norm(step) / n) < threshold:
    break

- np.linalg.norm(step) is the Frobenius norm over the (n,3) matrix
- dividing by n yields a rough “per-node” step scale (not exactly the mean step length, but close enough as a heuristic)

### 6.9 Determinism

The random seed  is used only for initial positions:

rng = np.random.default_rng(seed)
pos = constraint.sample(n, rng)


Everything else is deterministic given:
- the same node ordering (nodes = list(G))
- the same NumPy usually stable on the same machine

If you need strict reproducibility across runs:
- ensure graphs are built deterministically
- or explicitly sort nodes: nodes = sorted(G)

### 6.10 Limit cases

Empty graphs
- If n == 0, the k estimation divides by n and will fail.
- Recommended fix:
  - in Location.assign_positions_to_nodes(): skip if H.number_of_nodes() == 0
  or
  - inside spring_layout_3d_constrained(): return {} early when n == 0

Very small graphs
- For n == 1, the algorithm runs but forces are mostly zero -> constraints and recentering dominate. Returning the sampled point (translated by center) is typically fine.

Large graphs
- memory pressure from (n,n,3) arrays is the limiting factor.
- if a compartment has thousands of nodes, consider increasing min_degree to shrink H

### 6.11 How constraints interact with the layout

Constraint classes implement soft walls:
- they apply forces only when points violate a boundary or shell thickness
- force magnitude typically grows linearly with penetration depth (excess)
- the wall parameter inside constraint classes controls stiffness (default 5000.0)

Because the integrator normalizes displacement to t, the wall stiffness influences:
- the direction of motion strongly when boundary violations happen
- but not unbounded step sizes (steps are capped by temperature)

When creating new constraints, aim for:
- smooth, continuous forces and avoid discontinuities
- forces expressed in the same coordinate scale as sample() output
- reasonable wall defaults so boundary correction is effective without dominating everything

## 7. Plotting (CellModel.plot_all_locations_3d)

Plotly visualization combines:

- optional boundary traces from constraint.boundary_traces(), translated by center
- optional edges drawn as line segments within each compartment subgraph
- nodes drawn as markers, one trace per compartment

Translation is handled for both Surface and Scatter3d traces.

The figure hides axes and uses aspectmode="data" for proportional scaling.

Plotting uses loc.pos3d

If plotting is enabled, loc.pos3d must exist and contain all nodes referenced in H.edges(). This is currently true because the layout is computed on H and returned for exactly those nodes!!

---

## 8. I/O (io_helpers.py)

### 8.1 PPI loading
load_ppi_graph(ppi_path, sep="\t"):

- reads two columns named u, v
- removes self-loops (u == v)
- builds an undirected graph

Attention: the sep parameter exists but the function uses a fixed sep="\t". If your files flollowss another format, update this.

### 8.2 Localization loading
load_localizations(loc_path, sep="\t", protein_col="uniprot_id", location_col="location")

Returns:
- loc_to_proteins: Dict[location, Set[protein]]
- protein_to_locations: Dict[protein, List[location]]
- loc_df: pd.DataFrame

### 8.3 Color encoding
rgb_from_location(loc_name):
- derives an MD5 hash of the location name
- maps to a hue in HSV space with fixed saturation/value
- converts to integer RGB (0-255)
- ensures compartment colors are stable across runs

### 8.4 Writing positions
write_positions(cell, outdir, filename="positions_per_location.tsv") creates a row per (node, location):

- r,g,b from location color
- a from closeness-to-centroid within that location:
  - clamped to [100, 255] for visibility

Columns:

text
["Node ID","x","y","z","r","g","b","a","compartment"]


Return value: Path to output file.

---

## 9. CLI behavior (cli.py)

Arguments:

- required:
  - --ppi
  - --localizations
  - --config
  - --outdir
- selection:
  - --only comma-separated
  - --exclude comma-separated
- layout params:
  - --seed
  - --iterations
- plotting:
  - --plot
  - --plot-title
  - --no-boundaries (store_false)
  - --no-edges (store_false)

Filtering logic:

- skip locations missing from localization table (warning)
- apply --only include filter
- apply --exclude filter
- build one Location per remaining spec

Potential bugs:

If a selected location ends up with zero node_ids, the current code still constructs a Location. That can lead to an empty subgraph and failure inside the layout. Consider adding:

if not nodes:
    continue

to be continued :D

---

## 10. Extending the system

### 10.1 Add a new constraint type
1. Implement a new class in constraints.py with:
   - __init__(...)
   - sample(n, rng)
   - forces(q)
   - optional boundary_traces()
2. Add it to CONSTRAINTS in config.py.
3. Use it in config.json:


json filed -> "constraint": {"type": "MyNewConstraint", "params": {...}}


### 10.2 Add new per-location parameters
If you add new parameters to LocationSpec and Location:

- parse JSON in config.load_config()
- set the field in LocationSpec
- pass into Location() construction in cli.run()
- use inside Location.assign_positions_to_nodes() or the layout function


### 10.3 Performance improvements
Because the layout is O(n²):

- Replace dense A with sparse matrices for attraction-only computation.
- Reduce pairwise computations by chunking, especially for GPU acceleration.
- Use NetworkX layout only for initialization, then refine with constraints.

### 10.4 Better compartment isolation
Right now compartments are laid out independently and then translated.
If you want cross-compartment edge influence, you could:
- compute a global layout for all nodes,
- then add constraint forces and per-compartment pinning,
- or add inter-compartment edges with reduced weight.

---

## 11. Testing and reproducibility

With a fixed seed and fixed inputs, results should be deterministic because:
- the initial sampling uses NumPy default_rng(seed)
- the force updates are purely numeric

---
