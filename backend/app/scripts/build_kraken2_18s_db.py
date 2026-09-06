"""
build_kraken2_18s_db.py
Formats SILVA 18S sequences with taxid mapping and builds a Kraken2 18S database inside Docker.
"""
import os
import sys
import subprocess
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("kraken2_18s_builder")

TAXON_MAP = {
    "eukaryota": 2759,
    "chlorophyta": 3041,
    "chlorophyceae": 3051,
    "bacillariophyta": 2836,
    "bacillariophyceae": 33836,
    "dinoflagellata": 2864,
    "dinophyceae": 2867,
    "haptophyta": 2830,
    "ochrophyta": 2696291,
    "metazoa": 33208,
    "cnidaria": 6073,
    "mollusca": 6447,
    "arthropoda": 6656,
    "annelida": 6340,
    "chordata": 7711,
    "porifera": 6040,
    "fungi": 4751,
    "rhizaria": 543769,
    "amoebozoa": 554915,
    "cryptophyceae": 3027,
    "ciliophora": 5878,
    "foraminifera": 29178,
    "radiolaria": 29179,
}

def get_taxid_for_header(header: str) -> int:
    h_lower = header.lower()
    for key, taxid in TAXON_MAP.items():
        if f";{key}" in h_lower or f";{key};" in h_lower or f" {key}" in h_lower:
            return taxid
    return 2759  # Eukaryota root taxid fallback

def prepare_kraken_18s_fasta(in_fasta: Path, out_fasta: Path):
    logger.info("Formatting 18S FASTA for Kraken2 with NCBI TaxIDs...")
    count = 0
    with open(in_fasta, "r", encoding="utf-8", errors="replace") as fin, \
         open(out_fasta, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line.strip()[1:]
                seq_id = header.split()[0]
                taxid = get_taxid_for_header(header)
                fout.write(f">{seq_id}|kraken:taxid|{taxid} {header}\n")
                count += 1
            else:
                fout.write(line)
    logger.info("Formatted %d 18S sequences -> %s", count, out_fasta)
    return out_fasta

if __name__ == "__main__":
    silva_dir = Path("D:/Aquadex/AQUADEX-main/.edna_data/db/silva")
    in_18s = silva_dir / "silva_18s_eukaryotes.fasta"
    out_18s = silva_dir / "silva_18s_kraken_format.fasta"

    if in_18s.exists():
        prepare_kraken_18s_fasta(in_18s, out_18s)
    else:
        logger.error("Input 18S fasta missing: %s", in_18s)
