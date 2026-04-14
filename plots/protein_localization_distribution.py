import pandas as pd
import matplotlib.pyplot as plt

loc = pd.read_csv("SpatialPPI/data/protein_location_HPA_GO.tsv", sep="\t", dtype=str)

counts = (
    loc["location"]
    .value_counts()
    .reset_index(name="count")
)

print(counts)

plt.figure(figsize=(7, 4))
plt.bar(
    counts["location"],
    counts["count"],
    color="skyblue",      
    edgecolor="black",    
    linewidth=0.8,        
)
plt.xticks(rotation=75, ha="right")
plt.xlabel("Location")
plt.ylabel("Count")
plt.title("Top protein locations")
plt.tight_layout()
plt.savefig("location_top.png", dpi=180, bbox_inches="tight")
plt.ylim(0, 8000)
plt.grid(True, which="both", axis="y", linestyle="--", alpha=0.5)
plt.show()
