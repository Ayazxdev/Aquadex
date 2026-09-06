"""
build_kraken2_index.py
Builds the Kraken2 database index (.k2d files) inside Docker once taxonomy and library files are ready.
"""
import subprocess
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("kraken2_indexer")

def build_18s_kraken_db():
    logger.info("Starting Kraken2 build for 18S_SILVA database...")
    cmd = [
        "docker", "exec", "aquadex_backend",
        "kraken2-build", "--build",
        "--db", "/data/db/kraken2/18S_SILVA",
        "--threads", "4",
        "--kmer-len", "35",
        "--minimizer-len", "31",
        "--fast-build"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        logger.info("✅ Kraken2 18S database built successfully!")
        logger.info(result.stdout[:500])
    else:
        logger.error("Kraken2 build failed: %s", result.stderr)

if __name__ == "__main__":
    build_18s_kraken_db()
