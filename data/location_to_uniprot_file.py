import pandas as pd
import requests
import time
from io import StringIO

loc_df = pd.read_csv("data/protein_location_HPA_GO.tsv", sep="\t", dtype=str)

map_df = pd.read_csv("data/d_uniprot_to_symbol_20251120.tsv", sep="\t", dtype=str)


map_df = (map_df.drop_duplicates(subset=["symbol"]).set_index("symbol")["uniprotid"])
print(map_df.head())

#Add new column "uniprot_id" to localization data
loc_mapped = loc_df.copy() #make new df
loc_mapped["uniprot_id"] = loc_mapped["protein"].map(map_df) #create column with mapped IDs

print(loc_mapped.head())

#Safe Failed Mapping to file
loc_mapped[loc_mapped["uniprot_id"].isna()].to_csv("location_with_uniprot_fromfile_failed.tsv", sep="\t", index=False)

#Safe mapped localization data to file (includes failed mapping as missing numbers in uniprot_id)
loc_mapped.to_csv("location_with_uniprot_fromfile.tsv", sep="\t", index=False)