import pandas as pd
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"

pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


def gene_level_summary(df):
    """Her gen için klinik kategori dağılımını yüzde olarak hesaplar"""
    summary = df.groupby(["GeneSymbol", "clinical_category"]).size().unstack(fill_value=0)
    summary_pct = summary.div(summary.sum(axis=1), axis=0) * 100
    summary_pct["total_variants"] = summary.sum(axis=1)
    return summary_pct.round(1)


def domain_density(variants_df, domains_df):
    """Her domain için, uzunluğa göre normalize edilmiş varyant yoğunluğunu hesaplar"""
    records = []
    for _, domain_row in domains_df.iterrows():
        gene = domain_row["gene_symbol"]
        start, end = domain_row["start"], domain_row["end"]
        length = end - start + 1

        in_domain = variants_df[
            (variants_df["GeneSymbol"] == gene) &
            (variants_df["protein_position"] >= start) &
            (variants_df["protein_position"] <= end)
        ]

        pathogenic_count = (in_domain["clinical_category"] == "Pathogenic").sum()
        total_count = len(in_domain)
        pathogenic_pct = (pathogenic_count / total_count * 100) if total_count > 0 else 0

        records.append({
            "gene_symbol": gene,
            "domain_name": domain_row["domain_name"],
            "length_aa": length,
            "variant_count": total_count,
            "variants_per_100aa": round(total_count / length * 100, 2),
            "pathogenic_count": pathogenic_count,
            "pathogenic_pct": round(pathogenic_pct, 1),
        })

    return pd.DataFrame(records).sort_values("variants_per_100aa", ascending=False)


def effect_type_summary(df):
    """Etki türü bazlı klinik kategori dağılımı (tüm etki türleri için)"""
    summary = df.groupby(["effect_type", "clinical_category"]).size().unstack(fill_value=0)
    summary_pct = summary.div(summary.sum(axis=1), axis=0) * 100
    return summary_pct.round(1)


if __name__ == "__main__":
    df = pd.read_csv(DATA_PROCESSED / "brca_variants_with_domains.csv")
    domains_df = pd.read_csv(DATA_PROCESSED / "protein_domains.csv")

    print("=== Gen Bazlı Özet (%) ===")
    gene_summary = gene_level_summary(df)
    print(gene_summary)
    gene_summary.to_csv(DATA_PROCESSED / "stats_gene_summary.csv")

    print("\n=== Domain Yoğunluk Tablosu (En Yoğun 15) ===")
    density_df = domain_density(df, domains_df)
    print(density_df.head(15))
    density_df.to_csv(DATA_PROCESSED / "stats_domain_density.csv", index=False)

    print("\n=== Etki Türü x Klinik Kategori (%) ===")
    effect_summary = effect_type_summary(df)
    print(effect_summary)
    effect_summary.to_csv(DATA_PROCESSED / "stats_effect_type_summary.csv")

    print("\nTüm istatistik dosyaları data/processed/ içine kaydedildi.")