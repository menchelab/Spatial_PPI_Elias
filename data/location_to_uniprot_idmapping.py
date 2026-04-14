import pandas as pd
import requests
import time
from io import StringIO

API_URL = "https://rest.uniprot.org"
POLL_INTERVAL = 3  #seconds

#ID-Mapping functions

def submit_id_mapping(from_db, to_db, ids, taxId):
    """
    send IDs per ID-Mapping job to UniProt
    """
    ids = list(dict.fromkeys(ids))  #remove duplicates

    data = {
        "from": from_db,
        "to": to_db,
        "ids": ",".join(ids),
        "taxId": taxId
    }

    r = requests.post(f"{API_URL}/idmapping/run", data=data)
    r.raise_for_status() #raise Error if error occurs
    
    return r.json()["jobId"]

def wait_for_job(job_id):
    status_url = f"{API_URL}/idmapping/status/{job_id}"
    while True:
        r = requests.get(status_url)
        r.raise_for_status()
        j = r.json()
        if "jobStatus" not in j:   #UniProt documentation: results are ready
            return
        if j["jobStatus"] == "RUNNING":
            print(f"Job {job_id} running …")
            time.sleep(POLL_INTERVAL)
            continue
        raise RuntimeError(f"Job {job_id} finished with status {j['jobStatus']}")

def get_results_url(job_id):
    r = requests.get(f"{API_URL}/idmapping/details/{job_id}")
    r.raise_for_status()
    return r.json()["redirectURL"]


def download_tsv_results(results_url, fields):
    #switch to stream endpoint and load TSV
    if "/stream/" not in results_url:
        results_url = results_url.replace("/results/", "/results/stream/")
    params = {
        "format": "tsv",
        "fields": ",".join(fields),
    }
    r = requests.get(results_url, params=params)
    r.raise_for_status()
    return pd.read_csv(StringIO(r.text), sep="\t", dtype=str) #save http response as dataframe


#read Localization Data
loc_df = pd.read_csv("protein_location_HPA_GO.tsv", sep="\t", dtype=str)

#get all unique proteins to list
unique_proteins = loc_df["protein"].unique().tolist()
print(f"{len(unique_proteins)} different Proteins in this file.")

#Submit Job to UniProt
job_id = submit_id_mapping(
    from_db="Gene_Name",            #Map from Gene-Names
    to_db="UniProtKB-Swiss-Prot",   #To UniProt IDs
    ids=unique_proteins,            #Input: all unique proteins from localization Data
    taxId=9606                      #only Human
)
print("jobId:", job_id)

wait_for_job(job_id)
results_url = get_results_url(job_id)
print("results_url:", results_url)

map_df = download_tsv_results(
    results_url,
    fields=["accession"] #get UniProt ID column called "accession"
)

map_df = (map_df.drop_duplicates(subset=["From"]).set_index("From")["Entry"])

print(map_df.head())

#Add new column "uniprot_id" to localization data
loc_mapped = loc_df.copy() #make new df
loc_mapped["uniprot_id"] = loc_mapped["protein"].map(map_df) #create column with mapped IDs

#normalize names without spaces and lower-case
loc_mapped["location"] = loc_mapped["location"].astype(str).str.strip().str.lower().str.replace(r"\s+", "", regex=True)

print(loc_mapped.head())

#Safe Failed Mapping to file
loc_mapped[loc_mapped["uniprot_id"].isna()].to_csv("location_with_uniprot_failed.tsv", sep="\t", index=False)

#Safe mapped localization data to file (includes failed mapping as missing numbers in uniprot_id)
loc_mapped.to_csv("location_with_uniprot.tsv", sep="\t", index=False)