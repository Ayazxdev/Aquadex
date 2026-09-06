"""
convert_rna_to_dna.py
SILVA SSU reference sequences are stored as RNA (U instead of T).
Illumina sequencing reads are DNA (T not U).
This script converts the 18S SILVA library from RNA to DNA alphabet,
then rebuilds the Kraken2 index so that k-mer matching will work.
"""
import subprocess
import shutil
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("rna_to_dna")

OUT_DB_DIR = Path("C:/Users/tabas/.gemini/antigravity-ide/brain/a8f1b0f0-0102-4ca0-9d39-a7a004e06ff7/scratch/edna_data/db/kraken2/18S_SILVA")
LIB_FASTA = OUT_DB_DIR / "library.fasta"
LIB_DNA = OUT_DB_DIR / "library_dna.fasta"

def convert_rna_to_dna(in_fasta: Path, out_fasta: Path):
    logger.info("Converting RNA (U) → DNA (T) in %s ...", in_fasta)
    count = 0
    with open(in_fasta, "r", encoding="utf-8", errors="replace") as fin, \
         open(out_fasta, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(">"):
                fout.write(line)
                count += 1
            else:
                # Convert RNA to DNA: U→T, u→t; also lowercase to uppercase
                fout.write(line.replace("U", "T").replace("u", "t").upper())
    logger.info("✅ Converted %d sequences to DNA alphabet → %s", count, out_fasta)

def rebuild_kraken2_db():
    # Replace library.fasta with DNA version
    shutil.copy2(LIB_DNA, LIB_FASTA)
    logger.info("Replaced library.fasta with DNA version.")

    # Remove old library directory so kraken re-indexes
    lib_dir = OUT_DB_DIR / "library"
    if lib_dir.exists():
        shutil.rmtree(lib_dir)
    for k2d in OUT_DB_DIR.glob("*.k2d"):
        k2d.unlink()

    logger.info("Running kraken2-build --add-to-library ...")
    r1 = subprocess.run(
        ["docker", "exec", "aquadex_backend",
         "kraken2-build", "--add-to-library", "/data/db/kraken2/18S_SILVA/library.fasta",
         "--db", "/data/db/kraken2/18S_SILVA", "--no-masking"],
        capture_output=True, text=True
    )
    if r1.returncode != 0:
        logger.error("add-to-library failed: %s", r1.stderr)
        return False
    logger.info("Library staged: %s", r1.stdout.strip()[:200])

    logger.info("Running kraken2-build --build ...")
    r2 = subprocess.run(
        ["docker", "exec", "aquadex_backend",
         "kraken2-build", "--build",
         "--db", "/data/db/kraken2/18S_SILVA",
         "--threads", "4", "--fast-build"],
        capture_output=True, text=True
    )
    if r2.returncode != 0:
        logger.error("kraken2-build failed: %s", r2.stderr)
        return False
    logger.info("✅ Kraken2 18S DNA index built!\n%s", r2.stdout[-600:])
    return True

if __name__ == "__main__":
    convert_rna_to_dna(LIB_FASTA, LIB_DNA)
    rebuild_kraken2_db()
