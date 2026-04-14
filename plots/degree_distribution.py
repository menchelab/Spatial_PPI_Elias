import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ppi_path = "SpatialPPI/data/consensus_ppi_bioplex_biogrid_intact_huri_edgelist.tsv"

edges = pd.read_csv(
    ppi_path,
    sep="\t",
    header=None,
    usecols=[0, 1],
    names=["u", "v"],
    dtype=str,
    low_memory=False
)

# remove self-loops
edges = edges[edges["u"] != edges["v"]]

degrees = pd.concat([edges["u"], edges["v"]]).value_counts().sort_values(ascending=False)

# Plot
fig, ax = plt.subplots(figsize=(6, 4))

ax.hist(
    degrees.values,
    bins=100,
    edgecolor="black",
    color="skyblue",
    linewidth=0.7,
    alpha=0.8
)

ax.set_xlabel("degree (log scale)")
ax.set_ylabel("number of proteins")
ax.set_title("degree distribution of PPI network (log scale)")

ax.set_xlim(1, 600)
ax.grid(True, which="both", axis="y", linestyle="--", alpha=0.5)

plt.tight_layout()
fig.savefig("ppi_degree_hist_logx.png", dpi=180, bbox_inches="tight")

plt.show()
