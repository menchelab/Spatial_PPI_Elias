"""I/O helpers: load PPI network + localization table."""

import pandas as pd
import networkx as nx
from pathlib import Path
import hashlib
import colorsys
import numpy as np

def load_ppi_graph(ppi_path: str, sep: str = "\t") -> nx.Graph:
    """Load an undirected PPI graph from a two-column edge list file."""

    edges = pd.read_csv(ppi_path, sep="\t", names=["u", "v"], dtype=str)
    
    # remove self-loops
    edges = edges[edges["u"] != edges["v"]]
    
    G = nx.from_pandas_edgelist(edges, "u", "v")
    return G

def load_localizations(loc_path: str, sep: str = "\t", protein_col: str = "uniprot_id", location_col: str = "location",):
    """Load localization table and return mapping dicts.
        Returns
        -------
        loc_to_proteins:
            location -> set(protein IDs)
        protein_to_locations:
            protein ID -> sorted list of unique locations
    """

    loc_df = pd.read_csv(loc_path, sep=sep, dtype=str)

    #Mapping dictionaries
    
    loc_to_proteins = (         # Proteins in a single location
        loc_df
        .groupby(location_col)[protein_col]
        .apply(set)
        .to_dict()
    )

    protein_to_locations = (    # Locations of a single protein
    loc_df
    .groupby(protein_col)[location_col]
    .apply(lambda s: sorted(set(s)))
    .to_dict()
    )

    return loc_to_proteins, protein_to_locations, loc_df

def init_locations_attribute(G: nx.Graph) -> None:
    """Ensure every node has a locations dict attribute."""
    for node in G.nodes:
        G.nodes[node]["locations"] = {} #should be a dictionary

def get_nodes_without_location(G: nx.Graph, loc_df):
    all_nodes = set(G.nodes())
    return all_nodes.difference(loc_df["uniprot_id"])

def output_nodes_without_location(G: nx.Graph, loc_df, out_path: str):
    """creates .tsv file with all nodes, where no localiazation data is available"""
    all_nodes = set(G.nodes())
    nodes_without_loc = all_nodes.difference(loc_df["uniprot_id"])
    Path(out_path).write_text("\n".join(sorted(nodes_without_loc)))

def get_nodes_with_location(G: nx.Graph, loc_df):
    all_nodes = set(G.nodes())
    return all_nodes.intersection(loc_df["uniprot_id"])

def print_network(G):
    """Print-Loop for all nodes and its attributes of the Network"""
    for node, attrs in G.nodes(data=True):
        print(node, attrs)


def rgb_from_location(loc_name: str) -> tuple[int, int, int]:
    # stable integer from hash
    h = int.from_bytes(hashlib.md5(loc_name.encode("utf-8")).digest()[:4], "big")
    hue = (h % 360) / 360.0
    r, g, b = colorsys.hsv_to_rgb(hue, 0.7, 0.95)  # sat, val
    return int(r * 255), int(g * 255), int(b * 255)

def write_positions(cell, outdir: str | Path, filename: str = "positions_per_location.tsv"):
    """
    Write a TSV with columns:
    | Node ID | x | y | z | r | g | b | a |

    - RGB is derived from the location name via rgb_from_location()
    - Alpha a is based on closeness to the centroid of the location subnetwork:
        a = 255 for nodes at the centroid
        a -> 0 for nodes farthest from centroid within that location

    - Multi-localized proteins will appear as multiple rows (one per location).
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / filename

    rows = []

    for loc_name, loc in cell.locations.items():
        pos3d = getattr(loc, "pos3d", None) or {}
        if not pos3d:
            continue

        # location color
        r, g, b = rgb_from_location(loc_name)

        #centroid of this location points
        nodes = list(pos3d.keys())
        P = np.array([pos3d[n] for n in nodes], dtype=float)  # shape (N,3)
        C = P.mean(axis=0)                                    # (3,)

        #distances to centroid for alpha scaling
        d = np.linalg.norm(P - C[None, :], axis=1)
        dmax = float(d.max())
        if dmax < 1e-12:
            dmax = 0.0

        # write rows
        for n, (x, y, z) in pos3d.items():
            x = float(x); y = float(y); z = float(z)
            compartment = loc_name
            
            if dmax == 0.0:
                a = 255
            else:
                dist = float(np.linalg.norm(np.array([x, y, z], dtype=float) - C))
                closeness = 1.0 - (dist / dmax)      # 1 at center, 0 at farthest
                a = int(round(closeness * 255))
                if a < 100: a = 100
                if a > 255: a = 255

            rows.append({
                "Node ID": n,
                "x": x, "y": y, "z": z,
                "r": r, "g": g, "b": b, "a": a, "compartment": compartment
            })

    df = pd.DataFrame(rows, columns=["Node ID", "x", "y", "z", "r", "g", "b", "a", "compartment"])
    df.to_csv(out_path, sep="\t", index=False)
    return out_path