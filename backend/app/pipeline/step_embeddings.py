from pathlib import Path
import logging
import json
import h5py
import torch
from torch import nn
from transformers import AutoTokenizer, AutoModel
import numpy as np
from app.utils.io_helpers import compute_sha256
from typing import Dict

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class DNABERTSEmbeddings:
    def __init__(
        self,
        workdir: Path,
        model_path: str = "/app/models/dnabert",
        batch_size: int = None,
        device: str = None
    ):
        self.workdir = Path(workdir)
        self.model_path = Path(model_path)
        # Determine device with safe GPU initialization test
        if device:
            self.device = device
        else:
            self.device = "cpu"
            if torch.cuda.is_available():
                try:
                    # Test GPU memory allocation and CUDA context initialization
                    test_tensor = torch.zeros(1, device="cuda")
                    del test_tensor
                    self.device = "cuda"
                except Exception as e:
                    logger.warning("CUDA is reported available by PyTorch, but CUDA context initialization failed (%s). Falling back to CPU.", e)
                    self.device = "cpu"

        if batch_size is not None:
            self.batch_size = int(batch_size)
        else:
            self.batch_size = 128 if self.device == "cuda" else 32

        if self.workdir.name == "embeddings":
            self.output_dir = self.workdir
        else:
            self.output_dir = self.workdir / "embeddings"
        self.output_dir.mkdir(parents=True, exist_ok=True)

        logger.info("Loading tokenizer and model from %s on device: %s (batch_size=%d)", self.model_path, self.device, self.batch_size)
        self.tokenizer = AutoTokenizer.from_pretrained(
            pretrained_model_name_or_path=str(self.model_path.resolve()),
            local_files_only=True,
            trust_remote_code=True
        )
        loaded_model = AutoModel.from_pretrained(
            pretrained_model_name_or_path=str(self.model_path.resolve()),
            local_files_only=True,
            trust_remote_code=True
        )
        try:
            self.model = loaded_model.to(self.device)
        except Exception as e:
            logger.warning("Failed to move model to %s (%s). Falling back to CPU execution.", self.device, e)
            self.device = "cpu"
            self.batch_size = 32
            self.model = loaded_model.to("cpu")

        self.model.eval()

        self.max_len = getattr(getattr(self.model, "config", None), "max_position_embeddings", 512)
        logger.info("Using max token length = %d", self.max_len)

    def load_sequences(self, fasta_file: Path) -> Dict[str, str]:
        sequences = {}
        with open(fasta_file) as f:
            current_id = None
            current_seq = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if current_id:
                        sequences[current_id] = "".join(current_seq)
                    current_id = line[1:].strip()
                    current_seq = []
                else:
                    current_seq.append(line)
            if current_id:
                sequences[current_id] = "".join(current_seq)
        logger.info("Loaded %d sequences from %s", len(sequences), fasta_file)
        return sequences

    def batch_infer_embeddings(self, sequences: Dict[str, str], update_fn=None) -> Dict[str, np.ndarray]:
        asv_ids = list(sequences.keys())
        if self.device == "cpu" and len(asv_ids) > 100:
            logger.info("CPU execution detected: Capping embedding inference from %d to top 100 ASVs for performance.", len(asv_ids))
            asv_ids = asv_ids[:100]
        total_asvs = len(asv_ids)
        total_batches = (total_asvs + self.batch_size - 1) // self.batch_size
        embeddings: Dict[str, np.ndarray] = {}

        for batch_idx, i in enumerate(range(0, total_asvs, self.batch_size)):
            batch_ids = asv_ids[i:i + self.batch_size]
            batch_seqs = [sequences[aid] for aid in batch_ids]
            inputs = self.tokenizer(
                batch_seqs,
                return_tensors="pt",
                padding=True,
                truncation=True,
                max_length=self.max_len
            )
            inputs = {k: v.to(self.device) for k, v in inputs.items()}
            with torch.inference_mode():
                outputs = self.model(**inputs)
                last_hidden = outputs.last_hidden_state if hasattr(outputs, "last_hidden_state") else outputs[0]
                emb = last_hidden.mean(dim=1)
                emb = nn.functional.normalize(emb, p=2, dim=1)
                emb_cpu = emb.float().cpu().numpy()
                for idx, aid in enumerate(batch_ids):
                    embeddings[aid] = emb_cpu[idx]

            if update_fn and (batch_idx % max(1, total_batches // 10) == 0 or batch_idx == total_batches - 1):
                progress = 0.30 + 0.10 * ((batch_idx + 1) / total_batches)
                update_fn("embeddings", progress, f"DNABERT-S embedding progress: {min(i + self.batch_size, total_asvs)}/{total_asvs} ASVs ({batch_idx+1}/{total_batches} batches)")

        return embeddings

    def save_embeddings(self, embeddings: Dict[str, np.ndarray]) -> Path:
        h5_path = self.output_dir / "ASV_embeddings.h5"
        json_path = self.output_dir / "embedding_stats.json"

        with h5py.File(h5_path, "w") as hf:
            for asv_id, vec in embeddings.items():
                hf.create_dataset(asv_id, data=vec)

        emb_array = np.stack(list(embeddings.values()))
        stats = {
            "num_asvs": len(embeddings),
            "embedding_dim": int(emb_array.shape[1]),
            "mean": emb_array.mean(axis=0).tolist(),
            "std": emb_array.std(axis=0).tolist()
        }

        with open(json_path, "w") as f:
            json.dump(stats, f, indent=2)

        sha = compute_sha256(h5_path)
        logger.info("Saved %s with SHA-256 %s", h5_path, sha)
        return h5_path

    def run(self, fasta_file: Path, update_fn=None) -> Path:
        """
        Main entry point for embedding any FASTA file.
        Pass the pipeline-generated dereplicated FASTA here.
        """
        logger.info("Starting embedding pipeline on %s", fasta_file)
        sequences = self.load_sequences(fasta_file)
        if not sequences:
            raise ValueError(f"No sequences found in {fasta_file}")
        embeddings = self.batch_infer_embeddings(sequences, update_fn=update_fn)
        return self.save_embeddings(embeddings)


def embed_fasta(fasta_path: str, workdir_path: str, model_path: str, batch_size: int = 16):
    """
    Wrapper function to integrate with your pipeline.
    fasta_path: path to cleaned/dereplicated FASTA from frontend
    workdir_path: run-specific embeddings folder
    model_path: path to local DNABERT model folder
    """
    embedder = DNABERTSEmbeddings(
        workdir=Path(workdir_path),
        model_path=model_path,
        batch_size=batch_size
    )
    return embedder.run(Path(fasta_path))