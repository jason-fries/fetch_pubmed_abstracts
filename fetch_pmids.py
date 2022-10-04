""" Simple example script for bulk downloading PubMed abstacts using BioPython and eFetch.
Entrez reqires setting and email before use and (optionally) and api_key for heavy querying.

Example usage

python fetch_pmids.py \
--input pmids.txt \
--outdir documents/ \
--email "XXXX@XXXXX" \
--api_key "XXXXXXXXXXX"

"""
from typing import List
import os
import time
import json
from tqdm import tqdm
import pandas as pd
from Bio import Entrez
import argparse


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


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, help="PMID list", required=True)
    parser.add_argument("--outdir", type=str, help="output dir", required=True)
    parser.add_argument("--email", type=str, help="email address", required=True)
    parser.add_argument("--api_key", type=str, help="API key")
    parser.add_argument("--batch_size", type=int, default=500, help="batch size")
    args = parser.parse_args()

    # see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
    Entrez.email = args.email
    Entrez.api_key = args.api_key

    assert os.path.exists(args.input)
    if not os.path.exists(args.outdir):
        print(f"Creating {args.outdir}")
        os.mkdir(args.outdir)

    df = pd.read_csv(args.input, header=None, names=["pmid"], dtype={"pmid": str})
    pmids = df.pmid.to_list()
    print(f"File contains {len(pmids)} PMIDs")
    # with/without API key: 3 reqs per/second vs 10 reqs per/second
    delay = 0.15 if args.api_key is not None else 0.33

    fetch_pubmed_abstracts(pmids, args.outdir, args.batch_size, delay)

# # sample PMIDs
# document_ids = [
#     "18435798",
#     "27027316",
#     "18668432",
#     "22298808",
#     "19509253",
#     "24571714",
#     "15852660",
#     "24209620",
#     "27585851",
#     "16188502",
# ]
