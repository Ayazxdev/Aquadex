from pathlib import Path
import subprocess
import logging
from typing import List, Tuple, Optional
from utils.io_helpers import compute_sha256

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


class Denoiser:
    """
    Perform dereplication with vsearch and denoising with DADA2 (Rscript).
    """

    def __init__(
        self,
        workdir: Path,
        vsearch_path: str = "vsearch",
        rscript_path: str = "Rscript",
        dada2_script: str = "scripts/run_dada2.R"
    ):
        self.workdir = Path(workdir) / "denoise"
        self.workdir.mkdir(parents=True, exist_ok=True)

        self.vsearch_path = vsearch_path
        self.rscript_path = rscript_path
        self.dada2_script = Path(dada2_script)

    def dereplicate(self, input_files: List[Path]) -> List[Path]:
        derep_files = []
        for infile in input_files:
            outfile = self.workdir / f"derep_{infile.stem}.fasta"
            cmd = [
                self.vsearch_path,
                "--derep_fulllength", str(infile),
                "--output", str(outfile),
                "--sizeout",
                "--minuniquesize", "2"
            ]
            logging.info(f"Running vsearch dereplication on {infile}")
            subprocess.run(cmd, check=True)
            derep_files.append(outfile)
        return derep_files

    def run_dada2(self, derep_files: List[Path]) -> Tuple[Path, Path]:
        asv_fasta = self.workdir / "asvs.fasta"
        asv_table = self.workdir / "asv_table.tsv"

        for derep_file in derep_files:
            cmd = [
                self.rscript_path,
                str(self.dada2_script),
                str(derep_file),
                str(asv_fasta),
                str(asv_table)
            ]
            logging.info(f"Running DADA2 on {derep_file}")
            result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode != 0:
                raise RuntimeError(f"DADA2 failed: {result.stderr.decode()}")

        # Validate output
        if not asv_fasta.exists() or asv_fasta.stat().st_size == 0:
            raise RuntimeError("DADA2 produced empty ASV FASTA")
        if not asv_table.exists() or asv_table.stat().st_size == 0:
            raise RuntimeError("DADA2 produced empty ASV table")

        logging.info(f"ASV fasta checksum: {compute_sha256(asv_fasta)}")
        logging.info(f"ASV table checksum: {compute_sha256(asv_table)}")

        return asv_fasta, asv_table

    def run(self, r1: Path, r2: Optional[Path] = None, update_fn=None) -> Tuple[Path, Path]:
        if update_fn:
            update_fn("denoise", 0.1, "Starting dereplication")

        input_files = [r1] if r2 is None else [r1, r2]
        derep_files = self.dereplicate(input_files)

        if update_fn:
            update_fn("denoise", 0.5, "Running DADA2")

        asv_fasta, asv_table = self.run_dada2(derep_files)

        if update_fn:
            update_fn("denoise", 1.0, "Denoising complete")

        return asv_fasta, asv_table
