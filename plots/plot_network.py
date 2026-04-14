import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import re
from pathlib import Path

# read
df = pd.read_csv("ppi.tsv", sep="\t", header=None, comment="#", engine="python")

# Interpret the first two columns as u, v
df = df.iloc[:, :2].copy()
df.columns = ["u", "v"]

# remove self-loops
df = df[df["u"] != df["v"]]

G = nx.Graph()

edges = list(set(map(tuple, df[["u", "v"]].itertuples(index=False, name=None))))

G.add_edges_from(edges)
N = G.number_of_nodes()
M = G.number_of_edges()

#Select a subgraph for readable plotting

components = sorted(nx.connected_components(G), key=len, reverse=True)
print(components)
#print(len(components))

#Print component sizes
'''
sizes = [len(c) for c in components]
N = G.number_of_nodes()
for i, s in enumerate(sizes, start=1):
    print(f"Komponente {i}: {s} Knoten ({s/N:.1%} der Knoten)")
'''

largest = G.subgraph(components[0]).copy() #take largest component

#limit the plot size
MAX_NODES = 50
H = largest

# Select top nodes by degree and plot the induced subgraph
deg = dict(H.degree())
top_nodes = [n for n, _ in sorted(deg.items(), key=lambda kv: kv[1], reverse=True)[:MAX_NODES]]
H = H.subgraph(top_nodes).copy()

#Plot
plt.figure(figsize=(10, 8))

# Reproducible layout
pos = nx.spring_layout(H, seed=42)

# Node sizes ~ sqrt(degree) for better visibility
degH = dict(H.degree())
sizes = np.array([max(d, 1) for d in degH.values()], dtype=float)
sizes = 80 * np.sqrt(sizes)  # scaled moderately

#Show labels only for smaller graphs
with_labels = H.number_of_nodes() <= 60

# Draw – do not set colors explicitly (use Matplotlib/NetworkX defaults)
nx.draw_networkx(
    H, pos=pos,
    with_labels=with_labels,
    node_size=sizes,
    font_size=8 if with_labels else 0,
    width=0.6
)

plt.axis("off")
plt.title(f"PPI network (plotted subgraph): "
          f"{H.number_of_nodes()} nodes, {H.number_of_edges()} edges\n"
          f"Total: {N} nodes, {M} edges")

out_png = Path("ppi_network.png")
plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()

# Print summary
summary = {
    "All nodes": N,
    "All edges": M,
    "Plotted nodes": H.number_of_nodes(),
    "Plotted edges": H.number_of_edges(),
    "PNG": str(out_png),
}
summary
