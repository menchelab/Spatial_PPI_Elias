import networkx as nx
from pathlib import Path

path = "ppi.tsv"  # zwei Spalten: u v (optional weitere Spalten werden ignoriert)

# Als gerichteten Graphen einlesen
G = nx.read_edgelist(path, create_using=nx.DiGraph(), nodetype=str)

# Bidirektionalität/Reziproziät
reciprocity = nx.reciprocity(G)              # Anteil reziproker Kanten (0..1)
is_bidirectional = reciprocity == 1.0

# Basiskennzahlen
n = G.number_of_nodes()
m = G.number_of_edges()
weak_cc = nx.number_weakly_connected_components(G)
strong_cc = nx.number_strongly_connected_components(G)
density = nx.density(G.to_undirected())
avg_deg = sum(dict(G.degree()).values()) / n
clustering = nx.average_clustering(G.to_undirected())

print(f"Nodes: {n}, Edges: {m}")
print(f"Weakly CC: {weak_cc}, Strongly CC: {strong_cc}")
print(f"Density: {density:.6f}, Avg degree: {avg_deg:.3f}")
print(f"Average clustering: {clustering:.3f}")
print(f"Reciprocity: {reciprocity:.3f}  -> Bidirektional: {is_bidirectional}")
