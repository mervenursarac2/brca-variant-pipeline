import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


def clean_clinvar_data(df):
    """GRCh38'e sabitler, gereksiz kolonları filtreler"""
    df_clean = df[df["Assembly"] == "GRCh38"].copy()
    df_clean = df_clean.reset_index(drop=True)
    print(f"GRCh37 tekrarları temizlendikten sonra: {len(df_clean)} benzersiz varyant")
    return df_clean


def main():
    df = pd.read_csv(DATA_PROCESSED / "brca_clinvar_variants.csv")

    print("Toplam satır:", len(df))
    print("\nKolonlar:")
    print(df.columns.tolist())

    print("\nGen dağılımı:")
    print(df["GeneSymbol"].value_counts())

    print("\nKlinik sınıflandırma dağılımı:")
    print(df["ClinicalSignificance"].value_counts())

    print("\nİlk 5 satır:")
    print(df.head())

    df_clean = clean_clinvar_data(df)
    df_clean.to_csv(DATA_PROCESSED / "brca_clinvar_variants_clean.csv", index=False)


if __name__ == "__main__":
    main()