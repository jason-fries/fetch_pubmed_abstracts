""" Simple example script for bulk downloading PubMed abstacts using BioPython and eFetch.
Entrez reqires setting and email before use and (optionally) and api_key for heavy querying.
"""
from typing import List
import os
import time
import json
from tqdm import tqdm
from Bio import Entrez


Entrez.email = "<YOUR_EMAIL_HERE>"
# see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
# Entrez.api_key = "XXXXXXXXXXXXXXXX"


def fetch_pubmed_abstracts(
    pmids: List,
    outdir: str,
    batch_size: int = 100,
    delay: float = 0.333,
    overwrite: bool = False,
    verbose: bool = False,
):

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if not isinstance(pmids, list):
        pmids = list(pmids)

    n_chunks = int(len(pmids) / batch_size)
    for i in tqdm(
        range(0, n_chunks + (1 if n_chunks * batch_size < len(pmids) else 0))
    ):
        start, end = i * batch_size, (i + 1) * batch_size
        outfname = f"{outdir}/pubmed.{i}.json"
        if os.path.exists(outfname) and not overwrite:
            continue

        query = ",".join(pmids[start:end])
        handle = Entrez.efetch(
            db="pubmed", id=query, rettype="gb", retmode="xml", retmax=batch_size
        )
        record = Entrez.read(handle)
        if len(record["PubmedArticle"]) != len(pmids[start:end]) and verbose:
            print(
                f"Queried {len(pmids[start:end])}, returned {len(record['PubmedArticle'])}"
            )

        time.sleep(delay)
        # dump to JSON
        with open(outfname, "wt") as file:
            file.write(json.dumps(record, indent=2))


# sample PMIDs
document_ids = [
    "18435798",
    "27027316",
    "18668432",
    "22298808",
    "19509253",
    "24571714",
    "15852660",
    "24209620",
    "27585851",
    "16188502",
]
fetch_pubmed_abstracts(document_ids, ".pubmed_data/")
