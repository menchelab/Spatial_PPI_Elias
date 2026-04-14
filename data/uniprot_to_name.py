import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import re
from pathlib import Path
from io import StringIO

import time
import requests

API_URL = "https://rest.uniprot.org"
POLLING_INTERVAL = 3  # seconds

def submit_id_mapping(from_db: str, to_db: str, ids) -> str:
    """
    Starts a ID-Mapping-Job on UniProt and returns jobId. 
    """
    # UniProt limit: max. 100000 IDs per Job
    ids = list(dict.fromkeys(ids))  # remove duplicates, keep order
    data = {
        "from": from_db,
        "to": to_db,
        "ids": ",".join(ids),
    }
    r = requests.post(f"{API_URL}/idmapping/run", data=data)
    r.raise_for_status()
    return r.json()["jobId"]


def wait_for_results(job_id: str, poll_interval: int = POLLING_INTERVAL) -> None:
    """
    polls /idmapping/status/{jobId}, until results are ready
    """
    status_url = f"{API_URL}/idmapping/status/{job_id}"
    while True:
        r = requests.get(status_url)
        r.raise_for_status()
        j = r.json()

        #UniProt-Doku: wenn "jobStatus" nicht mehr vorhanden ist, enthält die Antwort bereits Ergebnisse/failedIds.
        if "jobStatus" not in j:
            return

        if j["jobStatus"] == "RUNNING":
            print(f"Job {job_id} läuft noch – retry in {poll_interval}s")
            time.sleep(poll_interval)
            continue
        else:
            raise RuntimeError(f"Job {job_id} finished with status {j['jobStatus']}")


def get_results_url(job_id: str) -> str:
    """
    gets redirectURL through /idmapping/details/{jobId}
    """
    details_url = f"{API_URL}/idmapping/details/{job_id}"
    r = requests.get(details_url)
    r.raise_for_status()
    return r.json()["redirectURL"]


def download_tsv_results(results_url: str, fields) -> str:
    """
    gets all Mapping-Results as TSV

    change to /stream/ Endpoint, like documented in UniProt
    "format=tsv" + fields 
    """
    # switch to Stream-Endpoint (all results, not paged)
    if "/stream/" not in results_url:
        results_url = results_url.replace("/results/", "/results/stream/")

    params = {
        "format": "tsv",
        "fields": ",".join(fields),
    }

    r = requests.get(results_url, params=params)
    r.raise_for_status()
    return r.text

def download_json_results(results_url: str):
    r = requests.get(results_url)
    r.raise_for_status()

    return r.json()

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

# Unique UniProt-IDs from v an u as list
all_ids = pd.unique(edges[["u", "v"]].values.ravel("K"))

print(f"{len(all_ids)} unique UniProt-IDs found")

#submit job
job_id = submit_id_mapping(
    from_db="UniProtKB_AC-ID",  # UniProt AC/ID database = Uniprot IDs
    to_db="UniProtKB",          # get UniProtKB-entries + Metadaten
    ids=all_ids
)
print("jobId:", job_id)

# wait until job is finished
wait_for_results(job_id)

# get url for results
results_url = get_results_url(job_id)
print("results_url:", results_url)

# retrieve TSV with fields
#   accession      – UniProt Accession
#   id             – Entry Name (z.B. TP53_HUMAN)
#   gene_names     – Genname
#   protein_name   – Proteinbezeichnung
tsv_text = download_tsv_results(
    results_url,
    fields=["id", "gene_names", "protein_name"],
)

#Get Mapping--------------------------------------------------------------

#TSV to DataFrame
map_df = pd.read_csv(StringIO(tsv_text), sep="\t")

#Remove _HUMAN suffix
map_df["Entry Name"] = map_df["Entry Name"].str.removesuffix("_HUMAN")

#Rename columns
map_df = map_df.rename(columns={
    "From": "uniprot_id",
    "Entry Name": "protein_name",            
    "Gene Names": "gene_names",
    "Protein names": "protein_full_name",
})

print(map_df.head())
print(map_df.columns)

out_file = "uniprot_mapping.tsv"

map_df.to_csv(out_file, sep="\t", index=False)

#Get failed Mapping--------------------------------------------------------
json_data = download_json_results(results_url)

failed_ids = json_data["failedIds"]

failed_df = pd.DataFrame({"From": failed_ids})
print(failed_df.head())
print(f"{len(failed_df)} failed IDs")
failed_df.to_csv("uniprot_failed_ids.tsv", sep="\t", index=False)

#Map PPI------------------------------------------------------------------------

'''
edges_u = edges.merge(
    map_df.rename(columns={"uniprot_id": "u", "protein_name": "u_name"}),
    on="u",
    how="inner"      # nur gemappte u behalten
)

# 2) Merge für v  (alle möglichen Mappings)
edges_uv = edges_u.merge(
    map_df.rename(columns={"uniprot_id": "v", "protein_name": "v_name"}),
    on="v",
    how="inner"      # nur gemappte v behalten
)

# 3) Nur noch die gemappten Namen als Kanten behalten
edges_mapped = (
    edges_uv[["u_name", "v_name"]]
    .rename(columns={"u_name": "u", "v_name": "v"})
    .drop_duplicates()   # optional: identische Kanten nur einmal behalten
)

print(edges_mapped.head())
print(edges_mapped.shape)
edges_mapped.to_csv("mapped_PPI.tsv", sep="\t", index=False)
'''