import requests
import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
DATA_PROCESSED.mkdir(parents=True, exist_ok=True)

ENSEMBL_BASE = "https://rest.ensembl.org"


def get_gene_with_transcripts(gene_symbol, species="homo_sapiens"):
    """Gen sembolünden, genişletilmiş (transcript'ler dahil) bilgi çeker"""
    url = f"{ENSEMBL_BASE}/lookup/symbol/{species}/{gene_symbol}"
    params = {"expand": 1}
    response = requests.get(url, headers={"Content-Type": "application/json"}, params=params)
    response.raise_for_status()
    return response.json()


def build_gene_transcript_map(gene_symbol):
    """Bir gen için gene-transcript-protein ilişki tablosu oluşturur"""
    gene_info = get_gene_with_transcripts(gene_symbol)
    gene_id = gene_info["id"]
    transcripts = gene_info.get("Transcript", [])

    records = []
    for t in transcripts:
        translation = t.get("Translation")
        protein_id = translation.get("id") if translation else None

        records.append({
            "gene_symbol": gene_symbol,
            "gene_id": gene_id,
            "transcript_id": t["id"],
            "protein_id": protein_id,
            "biotype": t.get("biotype"),
            "transcript_length": t.get("end", 0) - t.get("start", 0),
            "is_canonical": t.get("is_canonical", 0)
        })

    return pd.DataFrame(records)

def get_canonical_transcripts(df):
    """Sadece canonical (resmi referans) transcript'leri filtreler"""
    canonical_df = df[df["is_canonical"] == 1].reset_index(drop=True)
    return canonical_df


def main():
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)

    genes = ["BRCA1", "BRCA2"]
    all_records = []

    for gene in genes:
        print(f"\n--- {gene} için gene-transcript-protein haritası çekiliyor ---")
        df = build_gene_transcript_map(gene)
        all_records.append(df)
        print(df)

    final_df = pd.concat(all_records, ignore_index=True)
    final_df.to_csv(DATA_PROCESSED / "gene_transcript_protein_map.csv", index=False)
    print(f"\nToplam {len(final_df)} transcript kaydedildi: {DATA_PROCESSED / 'gene_transcript_protein_map.csv'}")

    canonical_df = get_canonical_transcripts(final_df)
    canonical_df.to_csv(DATA_PROCESSED / "canonical_transcripts.csv", index=False)
    print("\n--- Canonical Transcript'ler ---")
    print(canonical_df)

if __name__ == "__main__":
    main()