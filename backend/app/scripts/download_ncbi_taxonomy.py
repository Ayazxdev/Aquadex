"""
download_ncbi_taxonomy.py
Downloads and extracts clean NCBI taxdump files (names.dmp, nodes.dmp)
for building custom Kraken2 databases (18S, COI, 16S).
"""
import os
import tarfile
import urllib.request
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("taxonomy_downloader")

NCBI_TAXDUMP_URL = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

def download_and_extract_taxdump(target_dir: Path):
    target_dir.mkdir(parents=True, exist_ok=True)
    tar_path = target_dir / "taxdump.tar.gz"
    
    if not (target_dir / "names.dmp").exists() or not (target_dir / "nodes.dmp").exists():
        logger.info("Downloading NCBI taxdump.tar.gz from %s...", NCBI_TAXDUMP_URL)
        urllib.request.urlretrieve(NCBI_TAXDUMP_URL, tar_path)
        logger.info("Extracting taxonomy files into %s...", target_dir)
        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(path=target_dir)
        logger.info("✅ NCBI taxdump downloaded and extracted successfully.")
    else:
        logger.info("✅ NCBI taxonomy files (names.dmp, nodes.dmp) already present in %s", target_dir)

if __name__ == "__main__":
    db_tax_dir = Path("C:/Users/tabas/.gemini/antigravity-ide/brain/a8f1b0f0-0102-4ca0-9d39-a7a004e06ff7/scratch/edna_data/db/kraken2/18S_SILVA/taxonomy")
    download_and_extract_taxdump(db_tax_dir)
