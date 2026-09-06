# From pipeline/step_clustering.py
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import logging
from typing import Dict

logger = logging.getLogger(__name__)

MAX_UMAP_SAMPLES = 5000  # cap UMAP input for tractable runtime

# Try to import umap; we'll handle runtime failures gracefully
try:
    import umap
except Exception:
    umap = None

import hdbscan
from sklearn.decomposition import PCA

class ClusterRunner:
    """
    Load ASV embeddings, run UMAP (fallback to PCA) + HDBSCAN clustering, and produce a clusters table.
    """

    def __init__(self, workdir: Path, n_neighbors: int = 15, min_dist: float = 0.1,
                 n_components: int = 2, min_cluster_size: int = 2, metric: str = "euclidean"):
        self.workdir = workdir / "clustering"
        self.workdir.mkdir(parents=True, exist_ok=True)

        self.n_neighbors = n_neighbors
        self.min_dist = min_dist
        self.n_components = n_components
        self.min_cluster_size = min_cluster_size
        self.metric = metric

    def load_asv_embeddings(self, asv_h5_path: Path) -> Dict[str, np.ndarray]:
        embeddings: Dict[str, np.ndarray] = {}
        with h5py.File(asv_h5_path, "r") as hf:
            def visit(name, obj):
                if isinstance(obj, h5py.Dataset):
                    safe_name = name.replace("/", "_").replace(";", "_").replace(" ", "_")
                    arr = np.asarray(obj[()], dtype="float32")
                    if arr.ndim == 2 and arr.shape[0] == 1:
                        arr = arr.reshape(arr.shape[1],)
                    embeddings[safe_name] = arr
            hf.visititems(visit)
        return embeddings

    def _reduce(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, dtype="float32")
        if X.ndim != 2:
            raise ValueError("X must be 2D array")
        if umap is not None:
            try:
                reducer = umap.UMAP(
                    n_neighbors=self.n_neighbors,
                    min_dist=self.min_dist,
                    n_components=self.n_components,
                    metric=self.metric,
                    random_state=42,
                    n_jobs=-1,   # use all CPU cores
                )
                emb = reducer.fit_transform(X)
                return emb
            except TypeError as e:
                logger.warning("UMAP TypeError (compatibility). Falling back to PCA. Error: %s", e)
            except Exception as e:
                logger.warning("UMAP failed; falling back to PCA. Error: %s", e)
        else:
            logger.warning("UMAP not available; using PCA fallback.")
        pca = PCA(n_components=self.n_components, random_state=42)
        return pca.fit_transform(X)

    def run(self, asv_h5_path: Path, novelty_report_path: Path) -> Path:
        output_path = self.workdir / "clusters.tsv"
        emb_map = self.load_asv_embeddings(asv_h5_path)
        if not emb_map:
            raise RuntimeError(f"No embeddings found in {asv_h5_path}")
        nov_df = pd.read_csv(novelty_report_path, sep="\t")
        if "ASV_ID" not in nov_df.columns:
            raise ValueError("novelty_report must contain 'ASV_ID' column")
        rows = []
        missing = []
        for asv in nov_df["ASV_ID"].tolist():
            vec = emb_map.get(asv)
            if vec is None:
                missing.append(asv)
                continue
            rows.append((asv, vec))
        if not rows:
            raise RuntimeError("No ASV embeddings matched the novelty report ASV_IDs")
        if missing:
            logger.warning("%d ASV IDs from novelty report missing embeddings; they will be skipped.", len(missing))
        ids = [r[0] for r in rows]
        X = np.stack([r[1] for r in rows], axis=0)

        # Sub-sample if too large (UMAP scales poorly beyond ~5k points)
        n_total = len(ids)
        if n_total > MAX_UMAP_SAMPLES:
            logger.warning(
                "Clustering: %d ASVs exceeds MAX_UMAP_SAMPLES=%d; sub-sampling.",
                n_total, MAX_UMAP_SAMPLES
            )
            rng = np.random.default_rng(42)
            chosen = rng.choice(n_total, MAX_UMAP_SAMPLES, replace=False)
            chosen.sort()
            ids = [ids[i] for i in chosen]
            X = X[chosen]

        reduced = self._reduce(X)
        clusterer = hdbscan.HDBSCAN(min_cluster_size=self.min_cluster_size, metric="euclidean")
        labels = clusterer.fit_predict(reduced)
        out_df = pd.DataFrame({
            "ASV_ID": ids,
            "cluster_id": labels,
        })
        for i in range(reduced.shape[1]):
            out_df[f"dim_{i+1}"] = reduced[:, i]
        nov_map = nov_df.set_index("ASV_ID")["novelty_score"].to_dict() if "novelty_score" in nov_df.columns else {}
        out_df["novelty_score"] = out_df["ASV_ID"].map(nov_map)
        out_df["cluster_size"] = out_df.groupby("cluster_id")["ASV_ID"].transform("count")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(output_path, sep="\t", index=False)
        return output_path