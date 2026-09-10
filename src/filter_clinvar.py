import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
INPUT_PATH = PROJECT_ROOT / "data" / "raw" / "variant_summary.txt.gz"
OUTPUT_PATH = PROJECT_ROOT / "data" / "processed" / "brca_clinvar_variants.csv"
OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)


def filter_brca_variants(input_path, output_path, chunksize=100_000):
    genes_of_interest = {"BRCA1", "BRCA2"}
    matched_chunks = []

    for i, chunk in enumerate(pd.read_csv(input_path, sep="\t", chunksize=chunksize, low_memory=False)):
        filtered = chunk[chunk["GeneSymbol"].isin(genes_of_interest)]
        if len(filtered) > 0:
            matched_chunks.append(filtered)
        print(f"\rİşlenen satır: {(i+1)*chunksize:,}", end="")

    result = pd.concat(matched_chunks, ignore_index=True)
    result.to_csv(output_path, index=False)
    print(f"\nBulunan BRCA1/2 varyantı: {len(result)}")
    print(f"Kaydedildi: {output_path}")
    return result


if __name__ == "__main__":
    filter_brca_variants(INPUT_PATH, OUTPUT_PATH)