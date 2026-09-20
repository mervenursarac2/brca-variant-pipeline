import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
RESULTS = PROJECT_ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

CATEGORY_COLORS = {
    "Pathogenic": "#d62728",
    "Benign": "#2ca02c",
    "VUS": "#7f7f7f",
    "Conflicting": "#ff7f0e",
    "Other": "#c7c7c7",
}


def plot_lollipop(variants_df, domains_df, gene_symbol, protein_length, output_name):
    """Protein üzerinde domain'leri ve varyantları, kategori bazlı ayrı şeritlerde gösteren lolipop plot"""
    gene_variants = variants_df[
        (variants_df["GeneSymbol"] == gene_symbol) &
        (variants_df["protein_position"].notna())
    ]
    gene_domains = domains_df[domains_df["gene_symbol"] == gene_symbol]

    categories = ["Pathogenic", "Conflicting", "VUS", "Benign", "Other"]
    y_positions = {cat: i for i, cat in enumerate(categories)}

    fig, ax = plt.subplots(figsize=(16, 7))

    # Domain'leri arka planda dikey şeritler olarak göster (tüm sıraları kapsayacak şekilde)
    for i, row in gene_domains.iterrows():
        ax.axvspan(row["start"], row["end"], color="steelblue", alpha=0.15, zorder=1)

    # Her kategori için kendi satırında dikey lolipop çizgileri
    for category in categories:
        y = y_positions[category]
        cat_variants = gene_variants[gene_variants["clinical_category"] == category]
        positions = cat_variants["protein_position"].values
        color = CATEGORY_COLORS.get(category, "#cccccc")

        ax.vlines(positions, y, y + 0.7, color=color, linewidth=0.6, alpha=0.6, zorder=2)

    ax.set_xlim(0, protein_length)
    ax.set_ylim(-0.3, len(categories))
    ax.set_yticks([y + 0.35 for y in y_positions.values()])
    ax.set_yticklabels(categories)
    ax.set_xlabel("Amino Acid Position")
    ax.set_title(f"{gene_symbol} — Variant Distribution Across Protein Domains (by Clinical Category)")

    plt.tight_layout()
    plt.savefig(RESULTS / output_name, dpi=150)
    plt.close()
    print(f"Kaydedildi: {output_name}")

def plot_domain_pathogenic_rates(density_df, output_name):
    """Domain bazlı Pathogenic oranı bar chart"""
    top_domains = density_df.nlargest(12, "variant_count")
    top_domains = top_domains.sort_values("pathogenic_pct", ascending=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(top_domains["domain_name"] + " (" + top_domains["gene_symbol"] + ")",
            top_domains["pathogenic_pct"], color="crimson")
    ax.set_xlabel("Pathogenic Oranı (%)")
    ax.set_title("Domain Bazlı Pathogenic Varyant Oranı (En Çok Varyant Barındıran 12 Bölge)")

    plt.tight_layout()
    plt.savefig(RESULTS / output_name, dpi=150)
    plt.close()
    print(f"Kaydedildi: {output_name}")


def plot_clinical_category_distribution(df, output_name):
    """Genel klinik kategori dağılımı"""
    counts = df["clinical_category"].value_counts()
    colors = [CATEGORY_COLORS.get(cat, "#cccccc") for cat in counts.index]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(counts.index, counts.values, color=colors)
    ax.set_ylabel("Varyant Sayısı")
    ax.set_title("Genel Klinik Sınıflandırma Dağılımı (BRCA1 + BRCA2)")

    for i, v in enumerate(counts.values):
        ax.text(i, v + 200, str(v), ha="center", fontsize=9)

    plt.tight_layout()
    plt.savefig(RESULTS / output_name, dpi=150)
    plt.close()
    print(f"Kaydedildi: {output_name}")


def plot_effect_type_stacked(effect_summary_df, output_name):
    """Etki türü x klinik kategori stacked bar"""
    plot_df = effect_summary_df[["Pathogenic", "Benign", "VUS", "Conflicting"]]

    fig, ax = plt.subplots(figsize=(10, 6))
    plot_df.plot(kind="bar", stacked=True, ax=ax,
                 color=[CATEGORY_COLORS[c] for c in plot_df.columns])
    ax.set_ylabel("Yüzde (%)")
    ax.set_title("Etki Türüne Göre Klinik Sınıflandırma Dağılımı")
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

    plt.tight_layout()
    plt.savefig(RESULTS / output_name, dpi=150)
    plt.close()
    print(f"Kaydedildi: {output_name}")


if __name__ == "__main__":
    variants_df = pd.read_csv(DATA_PROCESSED / "brca_variants_with_domains.csv")
    domains_df = pd.read_csv(DATA_PROCESSED / "protein_domains.csv")
    density_df = pd.read_csv(DATA_PROCESSED / "stats_domain_density.csv")
    effect_summary_df = pd.read_csv(DATA_PROCESSED / "stats_effect_type_summary.csv", index_col=0)

    plot_lollipop(variants_df, domains_df, "BRCA1", 1863, "lollipop_brca1.png")
    plot_lollipop(variants_df, domains_df, "BRCA2", 3418, "lollipop_brca2.png")
    plot_domain_pathogenic_rates(density_df, "domain_pathogenic_rates.png")
    plot_clinical_category_distribution(variants_df, "clinical_category_distribution.png")
    plot_effect_type_stacked(effect_summary_df, "effect_type_stacked.png")

    print("\nTüm görseller results/ klasörüne kaydedildi.")