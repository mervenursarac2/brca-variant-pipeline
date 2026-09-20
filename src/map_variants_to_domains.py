import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


def find_domain(gene, position, domains_df):
    """Bir pozisyonun hangi domain(ler)e düştüğünü bulur"""
    if pd.isna(position):
        return None

    matches = domains_df[
        (domains_df["gene_symbol"] == gene) &
        (domains_df["start"] <= position) &
        (domains_df["end"] >= position)
    ]

    if len(matches) == 0:
        return "Outside known domain"

    # Birden fazla domain çakışabilir (örn. Region + Repeat aynı yerde), hepsini birleştir
    return "; ".join(matches["domain_name"].tolist())


def main():
    variants_df = pd.read_csv(DATA_PROCESSED / "brca_variants_annotated.csv")
    domains_df = pd.read_csv(DATA_PROCESSED / "protein_domains.csv")

    variants_df["domain"] = variants_df.apply(
        lambda row: find_domain(row["GeneSymbol"], row["protein_position"], domains_df),
        axis=1
    )

    print("Domain dağılımı (ilk 15):")
    print(variants_df["domain"].value_counts().head(15))

    print("\nSadece Missense varyantlar için domain x klinik kategori:")
    missense_df = variants_df[variants_df["effect_type"] == "Missense"]
    print(pd.crosstab(missense_df["domain"], missense_df["clinical_category"]))

    output_path = DATA_PROCESSED / "brca_variants_with_domains.csv"
        # BRCT domain'lerindeki "Other" kategorisini araştır
    other_brct = variants_df[
        (variants_df["domain"].isin(["BRCT 1", "BRCT 2"])) &
        (variants_df["clinical_category"] == "Other")
    ]
    print("\nBRCT domain'lerindeki 'Other' kategorisinin gerçek etiketleri:")
    print(other_brct["ClinicalSignificance"].value_counts())
    variants_df.to_csv(output_path, index=False)
    print(f"\nKaydedildi: {output_path}")

if __name__ == "__main__":
    main()