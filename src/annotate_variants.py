import pandas as pd
import re
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


def extract_protein_change(name):
    """ClinVar 'Name' kolonundan p. notasyonunu ayrıştırır"""
    match = re.search(r'\(p\.([A-Za-z]{3})(\d+)([A-Za-z]{3}|=|\*|fs\S*)\)', str(name))
    if not match:
        return None, None, None
    ref_aa, position, change = match.groups()
    return ref_aa, int(position), change


def classify_effect(ref_aa, change):
    """Değişiklik türünü sınıflandırır"""
    if ref_aa is None:
        return "Unknown/Non-coding"
    if "fs" in change:
        return "Frameshift"
    elif change in ("del", "delins"):
        return "In-frame Deletion"
    elif change == "dup":
        return "In-frame Duplication"
    elif change in ("*", "Ter"):
        return "Nonsense"
    elif change == "=":
        return "Synonymous"
    elif re.match(r'^[A-Za-z]{3}$', change):
        return "Missense"
    else:
        return "Other"


def annotate_variants(df):
    """Her varyant için protein değişikliği bilgisini çıkarır"""
    df = df.copy()
    extracted = df["Name"].apply(extract_protein_change)
    df["ref_amino_acid"] = extracted.apply(lambda x: x[0])
    df["protein_position"] = extracted.apply(lambda x: x[1])
    df["amino_acid_change"] = extracted.apply(lambda x: x[2])
    df["effect_type"] = df.apply(
        lambda row: classify_effect(row["ref_amino_acid"], row["amino_acid_change"]), axis=1
    )
    return df


def main():
    df = pd.read_csv(DATA_PROCESSED / "brca_variants_categorized.csv")
    df = annotate_variants(df)

    print("Etki türü dağılımı:")
    print(df["effect_type"].value_counts())

    print("\nEtki türü x Klinik kategori çapraz tablo:")
    print(pd.crosstab(df["effect_type"], df["clinical_category"]))

    print("\nÖrnek satırlar (protein değişikliği bulunan):")
    print(df[df["protein_position"].notna()][["GeneSymbol", "Name", "ref_amino_acid", "protein_position", "amino_acid_change", "effect_type"]].head(10))

    output_path = DATA_PROCESSED / "brca_variants_annotated.csv"
    df.to_csv(output_path, index=False)
    print(f"\nKaydedildi: {output_path}")

if __name__ == "__main__":
    main()