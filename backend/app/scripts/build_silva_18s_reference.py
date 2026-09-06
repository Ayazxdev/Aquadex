"""
build_silva_18s_reference.py
Extract 18S Eukaryotic SSU sequences from SILVA 138 SSURef NR99 FASTA database
and prepare marker-specific reference files for Kraken2 and FAISS vector embeddings.
"""
import os
import sys
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("silva_18s_builder")

def build_18s_fasta(silva_fasta_path: Path, output_dir: Path):
    if not silva_fasta_path.exists():
        logger.error("SILVA FASTA not found at: %s", silva_fasta_path)
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    out_18s_fasta = output_dir / "silva_18s_eukaryotes.fasta"
    out_16s_fasta = output_dir / "silva_16s_prokaryotes.fasta"

    logger.info("Parsing SILVA SSU database from %s...", silva_fasta_path)
    count_18s = 0
    count_16s = 0

    with open(silva_fasta_path, "r", encoding="utf-8", errors="replace") as fin, \
         open(out_18s_fasta, "w", encoding="utf-8") as f_18s, \
         open(out_16s_fasta, "w", encoding="utf-8") as f_16s:

        current_fh = None
        for line in fin:
            if line.startswith(">"):
                header = line.strip()
                if "Eukaryota;" in header:
                    current_fh = f_18s
                    count_18s += 1
                else:
                    current_fh = f_16s
                    count_16s += 1
                current_fh.write(line)
            elif current_fh:
                current_fh.write(line)

    logger.info("✅ Extracted %d Eukaryotic 18S sequences -> %s", count_18s, out_18s_fasta)
    logger.info("✅ Extracted %d Prokaryotic 16S sequences -> %s", count_16s, out_16s_fasta)
    return out_18s_fasta, out_16s_fasta

if __name__ == "__main__":
    base_dir = Path("D:/Aquadex/AQUADEX-main/.edna_data/db/silva")
    fasta = base_dir / "SILVA_138_SSURef_NR99_tax_silva.fasta"
    build_18s_fasta(fasta, base_dir)
