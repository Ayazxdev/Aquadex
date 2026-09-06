# pipeline/step_contam.py
from pathlib import Path
import subprocess
import logging
import pandas as pd
from typing import List, Dict
from utils.io_helpers import compute_sha256

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


class ContaminationScreen:
    """
    Perform rapid contamination screening using sourmash and BBMap/BBduk mapping.
    Implements blank subtraction logic.
    """

    def __init__(
        self,
        workdir: Path,
        sourmash_path: str = "sourmash",
        bbmap_path: str = "bbmap.sh",
        blank_threshold: int = 10
    ):
        self.workdir = workdir
        self.sourmash_path = sourmash_path
        self.bbmap_path = bbmap_path
        self.blank_threshold = blank_threshold
        self.contam_dir = workdir / "contamination"
        self.contam_dir.mkdir(parents=True, exist_ok=True)

    def run_sourmash(self, asv_fasta: Path, signature_db: Path) -> Path:
        """
        Compute MinHash signatures and compare against contamination DB.
        """
        report = self.contam_dir / "sourmash_report.json"
        cmd = [
            self.sourmash_path,
            "compute",
            str(asv_fasta),
            "--ksizes", "21,31,51",
            "--out", str(report)
        ]
        logging.info(f"Running sourmash compute: {cmd}")
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            logging.error(f"sourmash failed: {e.stderr.decode()}")
            raise RuntimeError(f"sourmash failed: {e.stderr.decode()}")
        return report

    def map_contaminants(self, asv_fasta: Path, ref_fasta_list: List[Path]) -> Path:
        """
        Map ASVs to contaminant references using BBMap.
        """
        report_tsv = self.contam_dir / "contamination_report.tsv"
        with open(report_tsv, "w") as out:
            for ref in ref_fasta_list:
                cmd = [
                    self.bbmap_path,
                    f"ref={ref}",
                    f"in={asv_fasta}",
                    f"outm={report_tsv}",
                    "nodisk",
                    "maxindel=3"
                ]
                logging.info(f"Running BBMap mapping: {cmd}")
                try:
                    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                except subprocess.CalledProcessError as e:
                    logging.error(f"BBMap failed: {e.stderr.decode()}")
                    raise RuntimeError(f"BBMap failed: {e.stderr.decode()}")
        return report_tsv

    def blank_subtraction(
        self,
        asv_table: Path,
        blank_sample_ids: List[str]
    ) -> Path:
        """
        Subtract ASV counts present in blanks and create cleaned table.
        """
        df = pd.read_csv(asv_table, sep="\t")
        blank_cols = [col for col in blank_sample_ids if col in df.columns]

        def mark_contaminant(row):
            blank_count = row[blank_cols].sum() if blank_cols else 0
            return blank_count >= self.blank_threshold

        df["is_contaminant"] = df.apply(mark_contaminant, axis=1)
        # Subtract blank counts from all samples
        if blank_cols:
            df.loc[:, df.columns.difference(["ASV_ID", "Sequence", "is_contaminant"])] -= df[blank_cols].sum(axis=1)

        blank_sub_file = self.contam_dir / "blank_subtraction_table.csv"
        df.to_csv(blank_sub_file, sep="\t", index=False)

        checksum = compute_sha256(blank_sub_file)
        logging.info(f"Blank-subtracted table saved: {blank_sub_file}, checksum={checksum}")
        return blank_sub_file