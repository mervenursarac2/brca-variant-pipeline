# BRCA Variant Analysis & Functional Annotation Pipeline

A computational biology pipeline designed to analyze genomic variants in *BRCA1* and *BRCA2* genes, map mutations to critical functional domains (RING, BRCT, BRC repeats), and assess their potential impact on Homologous Recombination (HR) DNA repair.

## Project Structure
- `src/parsers`: Sequence and VCF data ingestion.
- `src/analysis`: Codon translation, domain mapping, and functional scoring.
- `data/`: Reference sequences and sample variant datasets.