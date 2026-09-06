from pathlib import Path
import logging
import h5py
import faiss
import numpy as np
import pandas as pd
from typing import List, Dict ,Tuple
from app.utils.io_helpers import compute_sha256

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

class FAISSAnchoring:
    """
    Anchor ASV embeddings to reference embeddings using FAISS and compute anchor confidence scores.
    """
    def __init__(self, workdir: Path, ref_h5_path: Path = None, top_k: int = 5, device: str = None):
        self.workdir = workdir
        self.ref_h5_path = Path(ref_h5_path) if ref_h5_path else Path("/app/models/ref_embeddings.h5")
        self.top_k = top_k
        self.device = device or "cpu"
        self.output_dir = workdir / "anchor"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.index = None
        self.ref_ids = []
        self.ref_embeddings = None

    def load_reference_embeddings(self):
        if not self.ref_h5_path.exists():
            raise FileNotFoundError(f"Reference embeddings file not found: {self.ref_h5_path}")
        with h5py.File(self.ref_h5_path, "r") as hf:
            self.ref_ids = list(hf.keys())
            if not self.ref_ids:
                raise ValueError(f"No datasets found in reference HDF5: {self.ref_h5_path}")
            self.ref_embeddings = np.stack([hf[rid][:] for rid in self.ref_ids]).astype("float32")
        dim = self.ref_embeddings.shape[1]
        self.index = faiss.IndexFlatL2(dim)
        self.index.add(self.ref_embeddings)
        logger.info(f"FAISS index built with {len(self.ref_ids)} reference embeddings.")

    def load_asv_embeddings(self, asv_h5_path: Path) -> Dict[str, np.ndarray]:
        embeddings = {}
        with h5py.File(asv_h5_path, "r") as hf:
            def visit(name, obj):
                if isinstance(obj, h5py.Dataset):
                    safe_name = name.replace("/", "_").replace(";", "_").replace(" ", "_")
                    embeddings[safe_name] = obj[:].astype("float32")
            hf.visititems(visit)
        if not embeddings:
            raise ValueError(f"No ASV embeddings found in {asv_h5_path}")
        return embeddings

    def query_index(self, vec: np.ndarray) -> Tuple[List[str], List[float]]:

        D, I = self.index.search(vec.reshape(1, -1), self.top_k)
        neighbor_ids = [self.ref_ids[i] for i in I[0]]
        distances = D[0].tolist()
        return neighbor_ids, distances

    def compute_anchor_confidence(self, distances: List[float]) -> float:
        return float(np.exp(-np.mean(distances)))

    def run(self, asv_h5_path: Path) -> Path:
        logger.info(f"Using reference embeddings from: {self.ref_h5_path}")
        self.load_reference_embeddings()
        asv_embeddings = self.load_asv_embeddings(asv_h5_path)
        rows = []
        for asv_id, vec in asv_embeddings.items():
            top_refs, distances = self.query_index(vec)
            conf_score = self.compute_anchor_confidence(distances)
            rows.append({
                "ASV_ID": asv_id,
                "top_refs": ";".join(top_refs),
                "distances": ";".join(map(str, distances)),
                "anchor_confidence": conf_score
            })
        df = pd.DataFrame(rows)
        out_path = self.output_dir / "anchor_matches.tsv"
        df.to_csv(out_path, sep="\t", index=False)
        checksum = compute_sha256(out_path)
        logger.info(f"Anchor matches saved: {out_path} (sha256={checksum})")
        return out_path