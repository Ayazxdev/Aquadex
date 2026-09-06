# pipeline/step_taxonomy.py
from pathlib import Path
import subprocess
import logging
import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from utils.io_helpers import compute_sha256
from app.core.config import settings
from app.services.storage import out_dir
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# ✅ Required files for Kraken2 DB validation
REQUIRED_KRAKEN2_FILES = ["taxonomy", "hash.k2d", "opts.k2d", "taxo.k2d"]

def validate_kraken2_db(db_path: Path):
    if not db_path.exists():
        raise FileNotFoundError(f"Kraken2 DB folder not found: {db_path}")
    for f in REQUIRED_KRAKEN2_FILES:
        if not (db_path / f).exists():
            raise FileNotFoundError(f"Kraken2 DB missing required file: {db_path / f}")

class LayeredTaxonomy:
    def __init__(self, run_id: str, kraken2_path: str = "kraken2",
                 transformer_model_path: str = "/app/models/DNABERT-S",
                 db_path: Path | None = None):
        self.run_id = run_id
        self.workdir = out_dir(run_id) / "taxonomy"
        self.workdir.mkdir(parents=True, exist_ok=True)

        self.kraken2_path = kraken2_path
        self.transformer_model_path = transformer_model_path
        self.db_path = db_path or Path(settings.KRAKEN2_DB)
        validate_kraken2_db(self.db_path)

        self.report_file = None

    def run_kraken2(self, asv_fasta: Path, db_path: Path | None = None) -> Path:
        """
        Run Kraken2 on ASV sequences to assign taxonomy.
        """
        db_to_use = db_path or self.db_path
        validate_kraken2_db(db_to_use)

        output_tsv = self.taxonomy_dir / "kraken2.tsv"
        report_file = self.taxonomy_dir / "kraken2.report"

        cmd = [
            self.kraken2_path,
            "--db", str(db_to_use),
            "--output", str(output_tsv),
            "--report", str(report_file),
            str(asv_fasta)
        ]
        logging.info(f"Running Kraken2: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            err_msg = e.stderr.decode()
            logging.error(f"Kraken2 failed: {err_msg}")
            raise RuntimeError(f"Kraken2 failed: {err_msg}")

        self.report_file = report_file
        return output_tsv

    def run_transformer(self, asv_fasta: Path) -> Path:
        """
        Run transformer-based classifier to produce per-rank probabilities.
        """
        tokenizer = AutoTokenizer.from_pretrained(self.transformer_model_path, local_files_only=True)
        model = AutoModelForSequenceClassification.from_pretrained(self.transformer_model_path, local_files_only=True)
        output_tsv = self.taxonomy_dir / "transformer_taxonomy.tsv"

        seqs = {}
        with open(asv_fasta) as f:
            current_id = None
            current_seq = []
            for line in f:
                if line.startswith(">"):
                    if current_id:
                        seqs[current_id] = "".join(current_seq)
                    current_id = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_id:
                seqs[current_id] = "".join(current_seq)

        records = []
        model.eval()
        with torch.no_grad():
            for asv_id, seq in seqs.items():
                inputs = tokenizer(seq, return_tensors="pt", truncation=True, padding=True)
                logits = model(**inputs).logits.squeeze().softmax(dim=0)
                records.append({
                    "ASV_ID": asv_id,
                    "phylum_prob": float(logits[0]),
                    "class_prob": float(logits[1]),
                    "order_prob": float(logits[2]),
                    "family_prob": float(logits[3]),
                    "genus_prob": float(logits[4])
                })

        df = pd.DataFrame(records)
        df.to_csv(output_tsv, sep="\t", index=False)
        checksum = compute_sha256(output_tsv)
        logging.info(f"Transformer taxonomy saved: {output_tsv}, checksum={checksum}")
        return output_tsv

    def combine_consensus(self, kraken_tsv: Path, transformer_tsv: Path) -> Path:
        """
        Combine Kraken2 and transformer outputs into consensus taxonomy.
        """
        kraken_df = pd.read_csv(kraken_tsv, sep="\t", names=["ASV_ID", "kraken_taxon", "score"])
        transformer_df = pd.read_csv(transformer_tsv, sep="\t")

        consensus_records = []
        for _, row in transformer_df.iterrows():
            asv_id = row["ASV_ID"]
            kraken_tax = kraken_df.loc[kraken_df["ASV_ID"] == asv_id, "kraken_taxon"].values
            kraken_taxon = kraken_tax[0] if len(kraken_tax) > 0 else "Unassigned"

            family_match = (row["family_prob"] >= 0.85)
            genus_match = (row["genus_prob"] >= 0.85)
            consensus_label = kraken_taxon if family_match and genus_match else "Unassigned"

            consensus_records.append({
                "ASV_ID": asv_id,
                "consensus_taxon": consensus_label,
                "family_prob": row["family_prob"],
                "genus_prob": row["genus_prob"]
            })

        consensus_df = pd.DataFrame(consensus_records)
        output_consensus = self.taxonomy_dir / "consensus.tsv"
        consensus_df.to_csv(output_consensus, sep="\t", index=False)
        checksum = compute_sha256(output_consensus)
        logging.info(f"Consensus taxonomy saved: {output_consensus}, checksum={checksum}")
        return output_consensus

    def run(self, asv_fasta: Path, db_path: Path | None = None, update_fn=None) -> Path:
        """
        Run full layered taxonomy: Kraken2 + transformer + consensus.
        """
        if update_fn:
            update_fn("taxonomy", 0.22, "Running Kraken2 taxonomy...")
        kraken_tsv = self.run_kraken2(asv_fasta, db_path)

        if update_fn:
            update_fn("taxonomy", 0.25, "Running transformer model...")
        transformer_tsv = self.run_transformer(asv_fasta)

        if update_fn:
            update_fn("taxonomy", 0.27, "Combining consensus taxonomy...")
        consensus_tsv = self.combine_consensus(kraken_tsv, transformer_tsv)

        return consensus_tsv
