# Cell Layout - User Documentation

This document explains how to run and use my cell_layout command-line program that generates constrained 3D node coordinates for protein–protein interaction (PPI) subnetworks per subcellular compartment.

It intentionally focuses on practical usage (installation, inputs/outputs, configuration, troubleshooting).  
For the scientific background and the method write-up, see the GitHub repository.

https://github.com/menchelab/SpatialPPI/


---

## 1. What the program does

Given:

- a PPI edge list:        protein A interacts with protein B,
- a localization table:   protein assigned to one or more compartments,
- a JSON config:          defines each compartments geometry and layout parameters

1. builds the global PPI graph,
2. for each configured compartment:
   - selects the proteins localized there,
   - extracts the induced subgraph,
   - filters out weakly connected nodes (degree filter),
   - runs a constrained 3D spring layout inside the compartment geometry,
3. writes a single output table containing (x, y, z) for each (protein, compartment) pair,
4. optionally opens an interactive Plotly 3D visualization** of nodes, edges, and compartment boundaries.

Multi-localized proteins can appear in multiple compartments and will be exported once per compartment.

-> for a greater detail view the Developer Documentation.

---

## 2. Requirements and installation

### 2.1 Python version
The codebase is written for Python 3.12.2.

### 2.2 Python packages
Install these essentiall packages firstly:

- numpy
- pandas
- networkx
- plotly

---

## 3. Input files

### 3.1 PPI edge list
Argument flag: --ppi
Format: 2-column, undirected edge list.

Example:

text
Q9NWU5	O75431
Q9NWU5	P55809
...


Important notes
- The loader reads with a tab separator.
- If your file has a header row, it will be treated like data unless you remove it.
- Follow the Input file format, as otherwise the programm will not work correctly!
- No header row is expected.
- In the current implementation, the delimiter is a tab.

### 3.2 Localization table
Argument flag: --localizations
Format: tabular TSV with at least these columns:

- uniprot_id - protein identifier (must match IDs used in the PPI file)
- location - compartment label

Example:

text
uniprot_id	location
Q9NWU5	mitochondria
Q9NWU5	cytosol
O75431	mitochondria
...


Important notes
- Name matching matters: location labels must match the compartment names in the config exactly.

### 3.3 Configuration JSON
Argument flag: --config
Format: JSON object with

- a default block (global defaults)
- a locations list (one object per compartment)

---

## 4. Configuration reference

### 4.1 Top-level structure

json
{
  "default": {
    "seed": 7,
    "iterations": 350
  },
  "locations": [
    {
      "name": "mitochondria",
      "min_degree": 1,
      "center": [-5500.0, 2000.0, 0.0],
      "repulsion_strength": 1.0,
      "constraint": {
        "type": "EllipsoidConstraint",
        "params": {"axes": [325.0, 325.0, 650.0], "wall": 5000}
      }
    }
  ]
}


### 4.2 Fields per location

- name (string):              Compartment label. Must match localization location values.
- min_degree (int):           Nodes with degree ≤ min_degree are removed from that compartment subgraph before layout.
- center ([x,y,z]):           Compartment translation applied after optimization. Units are the same as your geometry values.
- repulsion_strength (float): Scales the repulsive term in the spring layout for this compartment.
- constraint.type (string):   Constraint class name (see below).
- constraint.params (object): Parameters passed to the constraint class constructor.

### 4.3 Supported constraint types and parameters

All constraints accept a wall parameter that controls how strongly out-of-bounds nodes are pushed back (a soft wall).

#### SphereConstraint
- params: radius (float), wall (float)

Defines a solid sphere centered at the local origin.

#### ShellConstraint
- Parameters: inner_radius (float), outer_radius (float), wall (float)

Defines a spherical shell. Nodes are kept between inner and outer radii.

#### EllipsoidConstraint
- Parameters: axes = [a,b,c], wall

Defines a solid ellipsoid with semi-axes a,b,c.

#### EllipsoidShellConstraint
- Parameters: axes = [a,b,c], outer, wall

Defines an ellipsoidal shell between:
- inner ellipsoid with axes [a,b,c]
- outer ellipsoid with axes [outer*a, outer*b, outer*c]

#### CylinderConstraint
- Parameters: radius (float), height (float), wall (float)

Defines a cylinder aligned with the z-axis. In practice, the layout recenters during optimization; treat height as total extent.

---

## 5. Running the program

### 5.1 simple run

From the directory containing cli.py:

python cli.py \
  --ppi path/to/ppi.tsv \
  --localizations path/to/localizations.tsv \
  --config path/to/config.json \
  --outdir results/


This writes:

- results/positions_per_location.tsv

and prints a short summary to stdout.

### 5.2 Restricting compartments

Run only two locations:

python cli.py ... --only mitochondria,nucleoplasm


Exclude a location:


python cli.py ... --exclude cytosol


### 5.3 Reproducibility controls

Override the random seed and iterations (CLI overrides config defaults):

bash
python cli.py ... --seed 123 --iterations 500


### 5.4 Plotly visualization

Open an interactive 3D plot:

bash
python cli.py ... --plot --plot-title "My 3D cell"


Turn off boundaries or edges:

bash
python cli.py ... --plot --no-boundaries
python cli.py ... --plot --no-edges


- For large subnetworks, edges can dominate rendering time; try --no-edges first!

---

## 6. Outputs

### 6.1 Main output: positions_per_location.tsv

The program writes a single TSV with columns:

- Node ID:      protein identifier (string)
- x, y, z:      3D coordinates (float) 
- r, g, b:      RGB color derived deterministically from the compartment name (0–255) 
- a:            alpha/opacity channel (100–255) 
- compartment:  compartment name (string) 

Alpha (a) is computed from the closeness of the node to the centroid within its compartment. Nodes near the centroid are more opaque.

Example:

text
Node ID	  x	             y	     z	     r	 g	 b	 a	compartment
Q9NWU5	-5287.2718	1965.9563	488.0782	20	72	242	100	mitochondria
O75431	-5809.2447	1988.6209	204.7252	19	55	55	110	nucleoplasm
P55809	-5186.3750	2004.5774	-165.1358	20	72	242	116	mitochondria
...


### 6.2 Multi-localized proteins

If a protein appears in multiple locations in the localization table, it may produce multiple rows-one for each compartment where it was placed.

---

## 7. Interpreting results and typical workflows

### 7.1 run checklist
1. Confirm your PPI IDs and uniprot_id use the same identifier space.
2. Confirm localization location strings match config name strings.
3. Start with one compartment (e.g., mitochondria) to validate the pipeline:
   
   python cli.py ... --only mitochondria --plot
   
4. Scale up to multiple compartments.

### 7.2 Tuning for readability
- Increase iterations for more relaxation -> slower but cleaner
- Increase a compartments repulsion_strength if nodes are too dense.
- Increase wall if nodes appear to leak outside boundaries

### 7.3 Handling sparse compartments
Because nodes with degree ≤ min_degree are removed, small or sparse compartments can end up with very few nodes.
- Lower min_degree (e.g., 0) if you want to retain more nodes.
- Or accept stronger filtering for more stable layouts.

---

## 8. Troubleshooting

### 8.1 No locations selected
You’ll see:
- No locations selected. Check --only/--exclude and your config.

Common causes:
- --only does not match any configured location names.
- All configured names are excluded with --exclude.
- Config names exist but those names do not appear in the localization table.

### 8.2 Unknown constraint type
Your config uses a constraint.type not supported by the code. Use one of:

- SphereConstraint
- ShellConstraint
- CylinderConstraint
- EllipsoidConstraint
- EllipsoidShellConstraint

### 8.3 Nothing placed in a compartment
Likely reasons:
- Most proteins in that location are not in the PPI file.
- Degree filtering removed many nodes (min_degree too strict).
- Location labels don’t match config names exactly.

### 8.4 Plot window does not appear
Plotly needs a browser. Click on the link provided by the program.

### 8.5 Very slow / memory errors
This is expected when a single compartment subgraph becomes large, because the layout uses dense O(n²) pairwise distance matrices.

How to reduce needed ressources:
- Run fewer locations at a time (--only).
- Increase min_degree to shrink the subgraph.
- Disable edges in plotting (--no-edges).
- Consider developer-side optimizations (see developer documentation).

---

## 9. Practical notes and limitations

- This is a visualization-oriented embedding. Coordinates are not physical simulations.
- Coordinate units are the same as your config geometry; keep them consistent.
- Layout results can differ by seed; use --seed to reproduce.
- The PPI reader expects a headerless TSV with two columns.

---