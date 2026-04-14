"""Cell model: locations, constrained layout, and Plotly visualization."""

import numpy as np
import networkx as nx
import plotly.graph_objects as go

from typing import Dict, List, Tuple

from constraints import spring_layout_3d_constrained

class Location:
    """
    Location (Organell) Object.

    Attributes:
    - name:             "Nucleus", "Cytosol" ..
    - nodes:            List of Proteins (Node-IDs) in this location
    - min_degree:       minimum number of degree per node
    - center:           3D coordinates of center (x, y, z)
    - constraint:       3D layout constrained by a certain shape
    """

    def __init__(
        self,
        name: str,
        node_ids: List[str],
        min_degree: int,
        center: Tuple[float, float, float],
        constraint: object,
        
        seed: int, 
        iterations: int,
        repulsion_strength: float,
    ):
        self.name = name
        self.node_ids = node_ids
        self.min_degree = min_degree
        self.center = center
        self.constraint = constraint
        
        self.repulsion_strength = repulsion_strength
        self.seed = seed
        self.iterations = iterations

        self.pos3d = {}

    #Subgraphs

    def make_subgraph(self, whole_graph: nx.Graph):
        """
        Subgraph from big PPI graph for this location with nodes with certain minimum of degrees (min_degree)
        """
        
        nodes_here = [n for n in self.node_ids if n in whole_graph]
        sg = whole_graph.subgraph(nodes_here).copy()
        sg.remove_nodes_from([n for n, d in sg.degree() if d <= self.min_degree]) # <- adjust degree here
        return sg

    # Simple metrics on Subgraphs

    def num_nodes(self, G: nx.Graph):
        H = self.make_subgraph(G)
        return H.number_of_nodes()

    def num_edges(self, G: nx.Graph):
        H = self.make_subgraph(G)
        return H.number_of_edges()
    

    def assign_positions_to_nodes(self, whole_graph: nx.Graph):
        """
        -Assigns all nodes of this location a (x,y,z) position (pos3D)
        -Saves them as Node Attributes in big whole-cell graph.
 
        """

        H = self.make_subgraph(whole_graph)
        
        self.pos3d = spring_layout_3d_constrained(H, self.constraint, seed=self.seed, iterations=self.iterations, center=self.center, repulsion_strength = self.repulsion_strength)
        print(self.name, " layout complete")

        #set Node attributes in big Graph
        for node in self.pos3d.keys():
            #set positions
            whole_graph.nodes[node]["locations"][self.name] = self.pos3d[node]

        

class CellModel:
    """
    holds the whole_cell_graph and a ditionary of all location objects.
    Everything should come together here.
    """

    def __init__(self, whole_cell_graph: nx.Graph):
        self.whole_cell_graph = whole_cell_graph
        self.locations: Dict[str, Location] = {}

    def add_location(self, location: Location):
        self.locations[location.name] = location

    def assign_all_positions(self):
        """
        calls the placement function for all locations
        """
        for loc in self.locations.values():
            loc.assign_positions_to_nodes(self.whole_cell_graph)

    def summary(self):
        print(f"Whole-cell graph: {self.whole_cell_graph.number_of_nodes()} nodes, "
              f"{self.whole_cell_graph.number_of_edges()} edges")
        for name, loc in self.locations.items():
            H = loc.make_subgraph(self.whole_cell_graph)
            print(f"- {name}: {H.number_of_nodes()} nodes, {H.number_of_edges()} edges")

    def plot_all_locations_3d(self, title="3D cell - all locations", show_boundaries=True, show_edges= True):
        fig = go.Figure()

        def _translate_trace(trace, center):
            dx, dy, dz = center

            # Surface (2D-Gitter)
            if getattr(trace, "type", None) == "surface":
                trace.x = (np.asarray(trace.x) + dx)
                trace.y = (np.asarray(trace.y) + dy)
                trace.z = (np.asarray(trace.z) + dz)
                return trace

            # Scatter3d (1D-Listen, ggf. mit None als Trenner)
            def shift_1d(arr, d):
                return [None if v is None else v + d for v in arr]

            trace.x = shift_1d(trace.x, dx)
            trace.y = shift_1d(trace.y, dy)
            trace.z = shift_1d(trace.z, dz)
            return trace

        for loc_name, loc in self.locations.items():
            H = loc.make_subgraph(self.whole_cell_graph)

            if show_boundaries and hasattr(loc.constraint, "boundary_traces"):
                for tr in loc.constraint.boundary_traces():
                    tr = _translate_trace(tr, loc.center)
                    tr.name = f"{loc_name}: {getattr(tr, 'name', 'boundary')}"
                    tr.showlegend = False
                    fig.add_trace(tr)

            # Edges
            if H.number_of_edges() > 0 and show_edges:
                edge_x, edge_y, edge_z = [], [], []
                for u, v in H.edges():
                    x0, y0, z0 = loc.pos3d[u]
                    x1, y1, z1 = loc.pos3d[v]
                    edge_x += [x0, x1, None]
                    edge_y += [y0, y1, None]
                    edge_z += [z0, z1, None]

                fig.add_trace(go.Scatter3d(
                    x=edge_x, y=edge_y, z=edge_z,
                    mode="lines", line=dict(width=1),
                    hoverinfo="none", opacity=0.35,
                    showlegend=False
                ))

            # Nodes
            xs, ys, zs, texts = [], [], [], []
            for node_id, (x, y, z) in loc.pos3d.items():
                xs.append(x); ys.append(y); zs.append(z)
                texts.append(f"{node_id} ({loc_name})")

            fig.add_trace(go.Scatter3d(
                x=xs, y=ys, z=zs,
                mode="markers", name=loc_name,
                hovertext=texts, hoverinfo="text",
                marker=dict(size=2, opacity=0.35)
            ))

        fig.update_layout(
            title=title,
            showlegend=True,
            margin=dict(l=0, r=0, b=0, t=40),
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                aspectmode="data"
            ),
        )
        fig.show()