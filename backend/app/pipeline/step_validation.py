# backend/app/pipeline/step_validation.py
from pathlib import Path
from typing import List, Dict, Any
import hashlib
import yaml
import json
from utils.subprocess_wrappers import run_cmd
from utils.io_helpers import write_json, write_yaml
from core.config import settings
from Bio import SeqIO
import logging

logger = logging.getLogger("edna_pipeline")

class ValidationStep:
    """
    Handles raw file recording, checksum calculation, basic FASTA validation and writing raw_manifest.
    """

    def __init__(self, run_id: str, raw_dir: Path, uploaded_paths: List[str], metadata: Dict[str, Any]):
        self.run_id = run_id
        self.raw_dir = Path(raw_dir)
        self.uploaded_paths = [Path(p) for p in uploaded_paths]
        self.metadata = metadata or {}
        self.manifest = {
            "run_id": run_id,
            "started_at": None,
            "raw_files": [],
            "metadata": self.metadata,
        }

    @staticmethod
    def compute_checksum(path: Path, algo: str = "sha256") -> str:
        h = hashlib.new(algo)
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(8192), b""):
                h.update(chunk)
        return h.hexdigest()

    def basic_fasta_parse(self, path: Path) -> Dict[str, Any]:
        """
        Parse FASTA/FASTQ with Biopython to detect structural errors and return basic stats.
        """
        stats = {"seq_count": 0, "errors": []}
        try:
            fmt = "fasta"
            # try FASTA first; if fails, try FASTQ
            try:
                recs = list(SeqIO.parse(str(path), "fasta"))
                if not recs:
                    recs = list(SeqIO.parse(str(path), "fastq"))
                    fmt = "fastq"
            except Exception:
                recs = list(SeqIO.parse(str(path), "fastq"))
                fmt = "fastq"
            stats["format"] = fmt
            stats["seq_count"] = len(recs)
        except Exception as e:
            stats["errors"].append(str(e))
            logger.exception("basic_fasta_parse failed for %s", path)
        return stats

    def normalize_headers(self, path: Path, out_path: Path) -> int:
        """
        Normalize sequence headers to SAMPLE|READID style (simple transformation).
        Returns number of records written.
        """
        count = 0
        with open(out_path, "w") as outfh:
            for idx, rec in enumerate(SeqIO.parse(str(path), "fasta")):
                norm = rec.id.replace(" ", "_")
                # ensure uniqueness
                norm = f"{norm}_{idx}"
                outfh.write(f">{norm}\n{str(rec.seq)}\n")
                count += 1
        return count

    def seqkit_stats(self, path: Path) -> str:
        """
        Run seqkit stats and return stdout (or empty string if tool not available).
        """
        try:
            out = run_cmd(["seqkit", "stats", str(path)])
            return out
        except Exception:
            logger.warning("seqkit not available or failed for %s", path)
            return ""

    def write_raw_manifest(self) -> Path:
        """
        Compute checksums and write raw_manifest.yaml
        """
        raw_entries = []
        for p in self.uploaded_paths:
            chk = self.compute_checksum(p)
            stats = self.basic_fasta_parse(p)
            seqkit_out = self.seqkit_stats(p)
            entry = {
                "filename": p.name,
                "path": str(p),
                "checksum_sha256": chk,
                "stats": stats,
                "seqkit": seqkit_out,
            }
            raw_entries.append(entry)
            self.manifest["raw_files"].append(entry)

        manifest_path = Path(settings.WORKDIR) / "raw" / self.run_id / "raw_manifest.yaml"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        write_yaml(self.manifest, manifest_path)
        return manifest_path

    def run_all(self) -> Dict[str, Any]:
        logger.info("ValidationStep starting for run_id=%s", self.run_id)
        self.manifest["started_at"] = None
        try:
            manifest_path = self.write_raw_manifest()
            # Normalize headers for each uploaded fasta and produce normalized files
            normalized_paths = []
            for p in self.uploaded_paths:
                if p.suffix.lower() not in [".fasta", ".fa", ".fastq", ".fq"]:
                    logger.warning("Skipping non-sequence file: %s", p)
                    continue
                out_norm = Path(settings.WORKDIR) / "raw" / self.run_id / f"normalized_{p.name}"
                # attempt normalization; if fails, copy original
                try:
                    count = self.normalize_headers(p, out_norm)
                    normalized_paths.append(str(out_norm))
                except Exception:
                    logger.exception("Header normalization failed for %s; copying original", p)
                    out_norm.write_bytes(p.read_bytes())
                    normalized_paths.append(str(out_norm))

            result = {
                "run_id": self.run_id,
                "manifest_path": str(manifest_path),
                "normalized": normalized_paths,
            }
            return result
        except Exception as e:
            logger.exception("ValidationStep.run_all failed: %s", e)
            raise