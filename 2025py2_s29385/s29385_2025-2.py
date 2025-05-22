#!/usr/bin/env python3
from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt

class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = 'BioScriptEx10'

    def search(self, taxid):
        print(f"Szukam rekordów dla taxID: {taxid}")
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        organism = Entrez.read(handle)[0]["ScientificName"]
        print(f"Organizm: {organism}")
        term = f"txid{taxid}[Organism]"
        handle = Entrez.esearch(db="nucleotide", term=term, usehistory="y")
        result = Entrez.read(handle)
        self.webenv = result["WebEnv"]
        self.query_key = result["QueryKey"]
        return int(result["Count"])

    def fetch(self, retstart=0, retmax=500):
        handle = Entrez.efetch(
            db="nucleotide",
            rettype="gb",
            retmode="text",
            retstart=retstart,
            retmax=retmax,
            webenv=self.webenv,
            query_key=self.query_key
        )
        return list(SeqIO.parse(handle, "genbank"))

def filter_records(records, min_len, max_len):
    return [r for r in records if min_len <= len(r.seq) <= max_len]

def export_csv(records, filename):
    data = [{
        "accession": r.id,
        "length": len(r.seq),
        "description": r.description
    } for r in records]
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    return df

def plot_lengths(df, filename):
    df_sorted = df.sort_values(by="length", ascending=False)
    plt.figure(figsize=(12, 6))
    plt.plot(df_sorted["accession"], df_sorted["length"], marker='o')
    plt.xticks(rotation=90, fontsize=6)
    plt.ylabel("Długość sekwencji")
    plt.xlabel("Numer akcesyjny")
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

def main():
    email = input("Podaj swój adres email do NCBI: ")
    api_key =  input("Podaj swój klucz API do NCBI: ")
    taxid = input("Podaj taxID organizmu: ")
    min_len = int(input("Minimalna długość sekwencji: "))
    max_len = int(input("Maksymalna długość sekwencji: "))

    retriever = NCBIRetriever(email, api_key)
    total = retriever.search(taxid)
    print(f"Znaleziono {total} rekordów. Pobieram...")

    all_records = []
    for start in range(0, total, 500):
        all_records += retriever.fetch(retstart=start, retmax=500)
        if len(all_records) > 1000:
            break  # ogranicz do 1000 rekordów maksymalnie

    filtered = filter_records(all_records, min_len, max_len)
    print(f"Po filtrze długości: {len(filtered)} rekordów")

    df = export_csv(filtered, f"taxid_{taxid}_filtered.csv")
    plot_lengths(df, f"taxid_{taxid}_plot.png")
    print("Zapisano CSV i wykres.")

if __name__ == "__main__":
    main()
