import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
SRC = PROJECT_ROOT / "src"
sys.path.insert(0, str(SRC))

DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
RESULTS = PROJECT_ROOT / "results"


def run_step(step_name, output_check_path, step_function, force=False):
    """Bir pipeline adımını çalıştırır; çıktısı zaten varsa atlar (force=True değilse)"""
    print(f"\n{'='*60}")
    print(f"ADIM: {step_name}")
    print('='*60)

    if output_check_path and output_check_path.exists() and not force:
        print(f"Atlanıyor: çıktı zaten mevcut ({output_check_path})")
        return

    try:
        step_function()
        print(f"✓ Tamamlandı: {step_name}")
    except Exception as e:
        print(f"✗ HATA ({step_name}): {e}")
        raise


def main(force=False):
    import sequence_analysis
    import blast_analysis
    import ensembl_mapping
    import download_clinvar
    import filter_clinvar
    import explore_clinvar
    import filter_variant
    import annotate_variants
    import domain_mapping
    import map_variants_to_domains
    import statistics_summary
    import visualize_results

    run_step("1. Sequence Analysis (FASTA indirme + GC + çeviri)",
              DATA_RAW / "brca1_protein.fasta", sequence_analysis.main, force)

    print("\n[NOT] BLAST adımı, önceden manuel indirilmiş XML dosyalarına bağımlı.")
    print("results/brca1_blast_results.xml ve brca2_blast_results.xml mevcut değilse,")
    print("önce web arayüzünden (blast.ncbi.nlm.nih.gov) manuel indirmen gerekiyor.")
    run_step("2. BLAST Sonuçlarını İşleme",
              RESULTS / "brca1_blast_summary.csv", blast_analysis.main, force)

    run_step("3. Ensembl Gene-Transcript-Protein Haritalama",
              DATA_PROCESSED / "canonical_transcripts.csv", ensembl_mapping.main, force)

    print("\n[NOT] ClinVar indirme adımı ~442 MB, ilk çalıştırmada uzun sürebilir.")
    run_step("4. ClinVar Verisi İndirme",
              DATA_RAW / "variant_summary.txt.gz", download_clinvar.main, force)

    run_step("5. BRCA1/2 Varyantlarını Filtreleme",
              DATA_PROCESSED / "brca_clinvar_variants.csv", filter_clinvar.main, force)

    run_step("6. Deduplication (GRCh37/38 Temizleme)",
              DATA_PROCESSED / "brca_clinvar_variants_clean.csv", explore_clinvar.main, force)

    run_step("7. Klinik Kategorileme",
              DATA_PROCESSED / "brca_variants_categorized.csv", filter_variant.main, force)

    run_step("8. Protein-Level Annotation",
              DATA_PROCESSED / "brca_variants_annotated.csv", annotate_variants.main, force)

    run_step("9. UniProt Domain Haritalama",
              DATA_PROCESSED / "protein_domains.csv", domain_mapping.main, force)

    run_step("10. Varyant-Domain Eşleştirme",
              DATA_PROCESSED / "brca_variants_with_domains.csv", map_variants_to_domains.main, force)

    run_step("11. İstatistik Özetleri",
              DATA_PROCESSED / "stats_gene_summary.csv", statistics_summary.main, force)

    run_step("12. Görselleştirme",
              RESULTS / "lollipop_brca1.png", visualize_results.main, force)

    print(f"\n{'='*60}")
    print("PIPELINE TAMAMLANDI")
    print('='*60)


if __name__ == "__main__":
    force_rerun = "--force" in sys.argv
    main(force=force_rerun)