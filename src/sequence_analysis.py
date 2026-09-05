import os
from pathlib import Path

# Script'in bulunduğu yerden proje kökünü bul
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_RAW = PROJECT_ROOT / "data" / "raw"
DATA_RAW.mkdir(parents=True, exist_ok=True)  # klasör yoksa otomatik oluşturur

from dotenv import load_dotenv
from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction

load_dotenv()  # .env dosyasını okur
Entrez.email = os.getenv("NCBI_EMAIL")


def fetch_and_save(accession, db, filepath):
    """NCBI'dan bir kaydı çekip FASTA olarak kaydeder"""
    handle = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    SeqIO.write(record, filepath, "fasta")
    return record


def load_sequence(filepath):
    """Kaydedilmiş FASTA dosyasını okur"""
    return SeqIO.read(filepath, "fasta")


def calculate_gc_content(seq):
    """GC içeriğini yüzde olarak hesaplar"""
    return gc_fraction(seq) * 100


def translate_sequence(seq):
    """mRNA dizisini proteine çevirir"""
    return seq.translate()

def verify_translation_v2(accession, protein_record):
    """CDS'ten çevrilen proteinin gerçek protein ile eşleşip eşleşmediğini kontrol eder"""
    cds_seq, translated_protein = fetch_cds_and_translate(accession)
    translated = str(translated_protein)
    real_protein = str(protein_record.seq)
    
    print(f"Çeviri uzunluğu: {len(translated)}")
    print(f"Gerçek protein uzunluğu: {len(real_protein)}")
    print(f"İlk 30 amino asit (çeviri): {translated[:30]}")
    print(f"İlk 30 amino asit (gerçek): {real_protein[:30]}")
    print(f"Tam eşleşiyor mu: {translated == real_protein}")

def fetch_cds_and_translate(accession):
    """GenBank kaydından CDS (coding sequence) bölgesini bulup çevirir"""
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    
    for feature in record.features:
        if feature.type == "CDS":
            cds_seq = feature.extract(record.seq)
            protein = cds_seq.translate(to_stop=True)  # stop codonda durur, * eklemez
            return cds_seq, protein
    
    raise ValueError(f"{accession} için CDS bulunamadı")


if __name__ == "__main__":
    brca1_mrna = fetch_and_save("NM_007294", "nucleotide", str(DATA_RAW / "brca1_mrna.fasta"))
    brca2_mrna = fetch_and_save("NM_000059", "nucleotide", str(DATA_RAW / "brca2_mrna.fasta"))
    brca1_protein = fetch_and_save("NP_009225", "protein", str(DATA_RAW / "brca1_protein.fasta"))
    brca2_protein = fetch_and_save("NP_000050", "protein", str(DATA_RAW / "brca2_protein.fasta"))

    print("BRCA1 GC content:", calculate_gc_content(brca1_mrna.seq))
    print("BRCA2 GC content:", calculate_gc_content(brca2_mrna.seq))

    print("\n--- BRCA1 Doğrulama (CDS ile) ---")
    verify_translation_v2("NM_007294", brca1_protein)

    print("\n--- BRCA2 Doğrulama (CDS ile) ---")
    verify_translation_v2("NM_000059", brca2_protein)