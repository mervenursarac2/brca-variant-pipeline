import os
from pathlib import Path
from dotenv import load_dotenv
from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd

load_dotenv()
Entrez.email = os.getenv("NCBI_EMAIL")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED = PROJECT_ROOT / "data" / "processed"
RESULTS = PROJECT_ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)


def load_protein_sequence(filepath):
    """Kaydedilmiş protein FASTA dosyasını okur"""
    record = SeqIO.read(filepath, "fasta")
    return record


def run_blast(protein_record, output_path, hitlist_size=10):
    """NCBI'da blastp çalıştırır, XML sonucunu kaydeder"""
    print(f"BLAST çalışıyor: {protein_record.id} ... (birkaç dakika sürebilir)")
    result_handle = NCBIWWW.qblast(
        program="blastp",
        database="swissprot",
        sequence=protein_record.seq,
        hitlist_size=hitlist_size
    )
    with open(output_path, "w") as f:
        f.write(result_handle.read())
    print(f"Sonuç kaydedildi: {output_path}")


def parse_blast_results(xml_path):
    """XML sonucunu okuyup DataFrame'e çevirir"""
    with open(xml_path) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    records = []
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            records.append({
                "hit_title": alignment.title,
                "e_value": hsp.expect,
                "identity_pct": round((hsp.identities / hsp.align_length) * 100, 2),
                "align_length": hsp.align_length,
                "query_coverage": round((hsp.align_length / blast_record.query_length) * 100, 2)
            })
    return pd.DataFrame(records)


if __name__ == "__main__":
    # BRCA1 için BLAST
    brca1_protein = load_protein_sequence(DATA_RAW / "brca1_protein.fasta")
    brca1_xml = RESULTS / "brca1_blast_results.xml"
    run_blast(brca1_protein, brca1_xml)

    brca1_df = parse_blast_results(brca1_xml)
    brca1_df.to_csv(RESULTS / "brca1_blast_summary.csv", index=False)
    print(brca1_df.head())