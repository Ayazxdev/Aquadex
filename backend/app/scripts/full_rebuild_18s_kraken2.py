"""
full_rebuild_18s_kraken2.py
FULL clean rebuild of the 18S_SILVA Kraken2 database.
Forces regeneration of seqid2taxid map from the fine-grained TaxID FASTA.
"""
import subprocess
import shutil
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("full_rebuild")

OUT_DB_DIR = Path("C:/Users/tabas/.gemini/antigravity-ide/brain/a8f1b0f0-0102-4ca0-9d39-a7a004e06ff7/scratch/edna_data/db/kraken2/18S_SILVA")
LIB_DNA = OUT_DB_DIR / "library_dna.fasta"

def run(cmd, desc=""):
    logger.info("Running: %s", " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.stdout.strip():
        logger.info(r.stdout[:400])
    if r.returncode != 0:
        logger.error("FAILED %s: %s", desc, r.stderr[:400])
        return False
    return True

def full_rebuild():
    # Step 1: Clean everything except library_dna.fasta and taxonomy/
    logger.info("Cleaning old index files...")
    for item in OUT_DB_DIR.iterdir():
        if item.name in ("library_dna.fasta", "taxonomy", "library.fasta"):
            continue
        if item.is_dir():
            shutil.rmtree(item)
            logger.info("Removed dir: %s", item.name)
        elif item.suffix in (".k2d", ".map", ".fasta"):
            item.unlink()
            logger.info("Removed: %s", item.name)
    
    # Copy DNA fasta to library.fasta
    shutil.copy2(LIB_DNA, OUT_DB_DIR / "library.fasta")
    logger.info("Copied library_dna.fasta → library.fasta")

    # Step 2: add-to-library (generates prelim_map with correct TaxIDs)
    if not run(
        ["docker", "exec", "aquadex_backend",
         "kraken2-build", "--add-to-library", "/data/db/kraken2/18S_SILVA/library.fasta",
         "--db", "/data/db/kraken2/18S_SILVA", "--no-masking"],
        "add-to-library"
    ): return False

    # Step 3: build (seqid2taxid is generated fresh from prelim_map in library/added/)
    if not run(
        ["docker", "exec", "aquadex_backend",
         "kraken2-build", "--build",
         "--db", "/data/db/kraken2/18S_SILVA",
         "--threads", "4", "--fast-build"],
        "build"
    ): return False

    logger.info("✅ Full rebuild complete!")
    
    # Verify seqid2taxid has correct taxids now
    seqid_file = OUT_DB_DIR / "seqid2taxid.map"
    if seqid_file.exists():
        lines = seqid_file.read_text(encoding="utf-8").strip().splitlines()[:5]
        logger.info("seqid2taxid sample:\n%s", "\n".join(lines))
    return True

if __name__ == "__main__":
    full_rebuild()
