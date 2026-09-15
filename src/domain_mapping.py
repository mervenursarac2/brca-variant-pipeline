import requests
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

UNIPROT_IDS = {
    "BRCA1": "P38398",
    "BRCA2": "P51587",
}


def fetch_domains(uniprot_id):
    """UniProt REST API'sinden domain/bölge bilgilerini çeker"""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    response = requests.get(url)
    response.raise_for_status()
    data = response.json()

    domains = []
    for feature in data.get("features", []):
        if feature["type"] in ("Domain", "Zinc finger", "Repeat", "Region"):
            domains.append({
                "domain_name": feature.get("description", "Unnamed"),
                "feature_type": feature["type"],
                "start": feature["location"]["start"]["value"],
                "end": feature["location"]["end"]["value"],
            })
    return pd.DataFrame(domains)


if __name__ == "__main__":
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)

    all_domains = []
    for gene, uniprot_id in UNIPROT_IDS.items():
        print(f"\n--- {gene} ({uniprot_id}) domain bilgisi çekiliyor ---")
        df = fetch_domains(uniprot_id)
        df["gene_symbol"] = gene
        all_domains.append(df)
        print(df)

    domains_df = pd.concat(all_domains, ignore_index=True)
    domains_df.to_csv(DATA_PROCESSED / "protein_domains.csv", index=False)
    print(f"\nKaydedildi: {DATA_PROCESSED / 'protein_domains.csv'}")