"""
build_silva_18s_with_fine_taxids.py
Rebuilds the 18S SILVA database with fine-grained NCBI TaxIDs mapped per
taxonomic lineage from SILVA headers — enabling Phylum/Class-level resolution.

Run this once. It:
 1. Reads NCBI taxonomy (names.dmp → name→taxid mapping)
 2. Processes SILVA 18S FASTA headers to match taxid from the lineage string
 3. Writes Kraken2-formatted FASTA with proper "|kraken:taxid|NNN" headers
 4. Triggers docker kraken2-build with the new library
"""
import os
import re
import subprocess
import logging
import shutil
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger("silva_18s_fine_builder")

TAXONOMY_DIR = Path("C:/Users/tabas/.gemini/antigravity-ide/brain/a8f1b0f0-0102-4ca0-9d39-a7a004e06ff7/scratch/edna_data/db/kraken2/18S_SILVA/taxonomy")
SILVA_18S_FASTA = Path("D:/Aquadex/AQUADEX-main/.edna_data/db/silva/silva_18s_eukaryotes.fasta")
OUT_DB_DIR = Path("C:/Users/tabas/.gemini/antigravity-ide/brain/a8f1b0f0-0102-4ca0-9d39-a7a004e06ff7/scratch/edna_data/db/kraken2/18S_SILVA")
OUT_FASTA = OUT_DB_DIR / "library_fine.fasta"

def load_ncbi_names(names_dmp: Path) -> dict:
    """Build name → taxid lookup from NCBI names.dmp (scientific + common names)."""
    logger.info("Loading NCBI name->taxid from %s ...", names_dmp)
    name2taxid = {}
    with open(names_dmp, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 4:
                continue
            taxid = int(parts[0])
            name = parts[1].lower()
            name_class = parts[3]
            # Prioritise scientific names
            if name_class in ("scientific name", "synonym", "common name", "includes"):
                if name not in name2taxid:
                    name2taxid[name] = taxid
    logger.info("Loaded %d name→taxid entries.", len(name2taxid))
    return name2taxid

def resolve_taxid_from_silva_lineage(lineage: str, name2taxid: dict) -> int:
    """
    SILVA headers look like:
    >AY846380.1.2583 Eukaryota;Archaeplastida;Chloroplastida;Chlorophyta;Chlorophyceae;Monoraphidium minutum
    
    Walk lineage right-to-left to find the deepest available NCBI taxid.
    """
    parts = [p.strip() for p in lineage.split(";")]
    # Walk from most specific (right) to least specific (left)
    for taxon in reversed(parts):
        name_lower = taxon.lower()
        if name_lower in name2taxid:
            return name2taxid[name_lower]
        # Try without trailing descriptor words ("class", "subphylum", etc.)
        cleaned = re.sub(r"\b(class|subphylum|subclass|order|family|genus|species|sp\.|bacteroidetes)\b", "", name_lower).strip()
        if cleaned in name2taxid:
            return name2taxid[cleaned]
    return 2759  # Fallback: Eukaryota root

def build_fine_kraken_fasta(name2taxid: dict) -> Path:
    logger.info("Building fine-grained TaxID FASTA from %s ...", SILVA_18S_FASTA)
    stats = {"resolved_deep": 0, "fallback": 0}
    
    with open(SILVA_18S_FASTA, "r", encoding="utf-8", errors="replace") as fin, \
         open(OUT_FASTA, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line.strip()[1:]
                # Header format: "AY846380.1.2583 Eukaryota;...;Genus species"
                space_idx = header.find(" ")
                if space_idx < 0:
                    seq_id = header
                    lineage = header
                else:
                    seq_id = header[:space_idx]
                    lineage = header[space_idx + 1:]
                
                taxid = resolve_taxid_from_silva_lineage(lineage, name2taxid)
                if taxid != 2759:
                    stats["resolved_deep"] += 1
                else:
                    stats["fallback"] += 1
                
                fout.write(f">{seq_id}|kraken:taxid|{taxid} {header}\n")
            else:
                fout.write(line)
    
    total = stats["resolved_deep"] + stats["fallback"]
    logger.info("✅ Fine-grained FASTA written: %d / %d resolved below Eukaryota root, %d fallback.",
                stats["resolved_deep"], total, stats["fallback"])
    return OUT_FASTA

def rebuild_kraken2_db():
    """Replace library, re-run kraken2-build inside Docker."""
    # Replace library.fasta
    dest = OUT_DB_DIR / "library.fasta"
    shutil.copy2(OUT_FASTA, dest)
    logger.info("Replaced library.fasta with fine-grained version.")
    
    # Remove old library folder (re-index)
    lib_dir = OUT_DB_DIR / "library"
    if lib_dir.exists():
        shutil.rmtree(lib_dir)
    
    # Remove old .k2d files
    for k2d in OUT_DB_DIR.glob("*.k2d"):
        k2d.unlink()
    
    logger.info("Removed old Kraken2 index files. Running kraken2-build...")
    
    # Step 1: add-to-library
    r1 = subprocess.run(
        ["docker", "exec", "aquadex_backend",
         "kraken2-build", "--add-to-library", "/data/db/kraken2/18S_SILVA/library.fasta",
         "--db", "/data/db/kraken2/18S_SILVA", "--no-masking"],
        capture_output=True, text=True
    )
    if r1.returncode != 0:
        logger.error("add-to-library failed: %s", r1.stderr)
        return False
    logger.info("Library added: %s", r1.stdout.strip())
    
    # Step 2: build index
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
    logger.info("✅ Kraken2 18S index rebuilt!\n%s", r2.stdout[-500:])
    return True

if __name__ == "__main__":
    names_dmp = TAXONOMY_DIR / "names.dmp"
    if not names_dmp.exists():
        logger.error("names.dmp not found at %s. Run download_ncbi_taxonomy.py first.", names_dmp)
        exit(1)
    
    name2taxid = load_ncbi_names(names_dmp)
    build_fine_kraken_fasta(name2taxid)
    rebuild_kraken2_db()
