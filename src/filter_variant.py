import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"


def load_clean_variants():
    """Temizlenmiş (deduplike edilmiş) ClinVar verisini okur"""
    return pd.read_csv(DATA_PROCESSED / "brca_clinvar_variants_clean.csv")


def categorize_clinical_significance(df):
    """ClinicalSignificance kolonunu sadeleştirilmiş kategorilere ayırır"""
    def categorize(value):
        value = str(value).lower()
        if "pathogenic" in value and "conflicting" not in value:
            return "Pathogenic"
        elif "benign" in value and "conflicting" not in value:
            return "Benign"
        elif "conflicting" in value:
            return "Conflicting"
        elif "uncertain" in value or "vus" in value:
            return "VUS"
        else:
            return "Other"

    df = df.copy()
    df["clinical_category"] = df["ClinicalSignificance"].apply(categorize)
    return df


def filter_low_quality(df):
    """Düşük kaliteli/tek submitter'lı kayıtları işaretler (elemez, işaretler)"""
    df = df.copy()
    # ReviewStatus'a göre güvenilirlik seviyesi ekleyelim
    df["is_high_confidence"] = df["ReviewStatus"].str.contains(
        "expert panel|multiple submitters", case=False, na=False
    )
    return df


def main():
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)

    df = load_clean_variants()
    df = categorize_clinical_significance(df)
    df = filter_low_quality(df)

    print("Kategori dağılımı:")
    print(df["clinical_category"].value_counts())

    print("\nYüksek güvenilirlikli (expert panel / multiple submitters) oranı:")
    print(df["is_high_confidence"].value_counts())

    print("\nKategori x Güvenilirlik çapraz tablo:")
    print(pd.crosstab(df["clinical_category"], df["is_high_confidence"]))

    output_path = DATA_PROCESSED / "brca_variants_categorized.csv"
    df.to_csv(output_path, index=False)
    print(f"\nKaydedildi: {output_path}")

    crosstab_df = pd.crosstab(df["clinical_category"], df["is_high_confidence"])
    crosstab_df.to_csv(DATA_PROCESSED / "clinical_category_confidence_crosstab.csv")

if __name__ == "__main__":
    main()