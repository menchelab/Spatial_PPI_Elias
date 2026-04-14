import pandas as pd
import requests
import time

# 1) Localization-TSV einlesen
loc_path = "SpatialPPI/data/protein_location_HPA_GO.tsv"   # ggf. Pfad anpassen
loc_df = pd.read_csv(loc_path, sep="\t")

# Erwartet Spalten: "protein", "location", "source(s)"
print(loc_df.head())

# 2) UniProt-Search-API vorbereiten
BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
HEADERS = {"accept": "application/json"}

def fetch_uniprot_id(gene_symbol: str) -> str | None:
    """
    Hole die UniProt-Accession für ein Gen-Symbol.
    - nur reviewed (Swiss-Prot)
    - nur Human (organism_id:9606)
    nimmt den ersten Treffer.
    """
    if not gene_symbol or pd.isna(gene_symbol):
        return None

    query = f"gene_exact:{gene_symbol} AND reviewed:true AND organism_id:9606"

    params = {
        "query": query,
        "fields": "accession,protein_name",
        "size": 1,      # nur bester Treffer
    }

    r = requests.get(BASE_URL, headers=HEADERS, params=params)
    r.raise_for_status()
    data = r.json()

    results = data.get("results", [])
    if not results:
        return None

    # laut UniProt-JSON-Struktur
    # https://rest.uniprot.org/uniprotkb/search ... accept: application/json 
    return results[0]["primaryAccession"]

# 3) Für alle einzigartigen Proteine IDs holen
genes = loc_df["protein"].dropna().unique()
print(f"{len(genes)} unterschiedliche Proteine in der Datei")

gene2acc = {}

for g in genes:
    acc = fetch_uniprot_id(g)
    gene2acc[g] = acc
    print(f"{g} -> {acc}")
    #time.sleep(0.2)  # UniProt nicht zu hart spammen (~3 Requests/s empfohlen)

# 4) Mapping zurück ins DataFrame
loc_df["uniprot_id"] = loc_df["protein"].map(gene2acc)

print(loc_df.head())

# Optional: nur Zeilen behalten, wo eine UniProt-ID gefunden wurde
loc_df_mapped = loc_df.dropna(subset=["uniprot_id"]).copy()
print(loc_df_mapped.head())
print(loc_df_mapped.shape)

# Optional speichern
loc_df_mapped.to_csv("location_with_uniprot.tsv",
                     sep="\t", index=False)
