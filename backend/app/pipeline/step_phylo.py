from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage
from typing import Dict, List, Tuple

class PhylogenyRunner:
    """
    Build an embedding-based hierarchical tree (Newick) and a summary table.
    Writes:
      - <run_root>/phylogeny/tree.nwk
      - <run_root>/phylogeny/phylogeny_table.tsv
    """
    def __init__(self,
                 workdir: Path = Path("/app/workdir"),
                 clusters_tsv: Path = None,
                 novelty_tsv: Path = None):
        self.workdir = Path(workdir)
        self.asv_h5 = None  # will be set dynamically in load_embeddings()

        # allow either clusters.tsv or clusters.csv (tolerate pipeline name differences)
        self.clusters_tsv = Path(clusters_tsv) if clusters_tsv else self.workdir / "clustering" / "clusters.tsv"
        self.novelty_tsv = Path(novelty_tsv) if novelty_tsv else self.workdir / "novelty" / "novelty_report.tsv"
        self.out_dir = self.workdir / "phylogeny"
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.tree_path = self.out_dir / "tree.nwk"
        self.table_path = self.out_dir / "phylogeny_table.tsv"

    def _safe_name(self, name: str) -> str:
        return name.replace("/", "_").replace(";", "_").replace(" ", "_")

    def load_embeddings(self) -> Dict[str, np.ndarray]:
        """
        Robustly locate and load ASV embeddings HDF5. Tries multiple candidate paths,
        then does a recursive search under the run folder if needed.
        Returns a dict mapping safe ASV_ID -> numpy array embedding.
        """
        candidates = [
            self.workdir / "embeddings" / "ASV_embeddings.h5",
            self.workdir / "embeddings" / "embeddings" / "ASV_embeddings.h5",
            self.workdir / "embeddings" / "ASV_embeddings.hdf5",
            self.workdir / "embeddings" / "embeddings" / "ASV_embeddings.hdf5",
        ]

        found = None
        for c in candidates:
            if c.exists():
                found = c
                break

        if not found:
            # recursive search under the run folder for any plausible filename
            cand_list = list(self.workdir.rglob("ASV_embeddings*.h5")) + list(self.workdir.rglob("ASV_embeddings*.hdf5"))
            if cand_list:
                found = cand_list[0]

        if not found:
            # As a last-ditch attempt, check top-level /app/workdir (keeps previous behavior visible)
            alt = Path("/app/workdir") / "embeddings" / "ASV_embeddings.h5"
            if alt.exists():
                found = alt

        if not found:
            raise FileNotFoundError(f"ASV embeddings HDF5 not found in run folder {self.workdir}")

        self.asv_h5 = found
        print(f"[step_phylogeny] Loading embeddings from: {self.asv_h5}")

        emb: Dict[str, np.ndarray] = {}
        with h5py.File(self.asv_h5, "r") as hf:
            def visit(name, obj):
                if isinstance(obj, h5py.Dataset):
                    safe = self._safe_name(name)
                    arr = np.array(obj[()]).astype("float32")
                    emb[safe] = arr
            hf.visititems(visit)

        if not emb:
            raise RuntimeError(f"No embeddings found in HDF5: {self.asv_h5}")

        return emb

    def choose_asvs(self, emb_keys: List[str]) -> List[str]:
        ids = []

        # Try both TSV and CSV cluster files (tolerate naming differences)
        possible_cluster_files = [self.clusters_tsv, self.clusters_tsv.with_suffix(".csv")]
        for cf in possible_cluster_files:
            if cf.exists():
                try:
                    dfc = pd.read_csv(cf, sep="\t" if cf.suffix == ".tsv" else ",", dtype=str)
                    if "ASV_ID" in dfc.columns:
                        ids = dfc["ASV_ID"].dropna().astype(str).tolist()
                        break
                except Exception:
                    ids = []

        if not ids and self.novelty_tsv.exists():
            try:
                dfn = pd.read_csv(self.novelty_tsv, sep="\t", dtype=str)
                if "ASV_ID" in dfn.columns:
                    ids = dfn["ASV_ID"].dropna().astype(str).tolist()
            except Exception:
                ids = []

        if not ids:
            ids = list(emb_keys)

        ids_sanitized = [self._safe_name(i) for i in ids]
        seen = set()
        out = []
        for i in ids_sanitized:
            if i in seen:
                continue
            seen.add(i)
            out.append(i)
        return out

    def _linkage_to_newick(self, Z: np.ndarray, labels: List[str]) -> str:
        n = len(labels)
        if n == 1:
            return f"{labels[0]};"
        m = Z.shape[0]
        tree = {}
        heights = {i: 0.0 for i in range(n)}
        for i in range(m):
            left = int(Z[i,0])
            right = int(Z[i,1])
            node_id = n + i
            h = float(Z[i,2])
            tree[node_id] = (left, right, h)
            heights[node_id] = h
        root = n + m - 1
        def build(node_id: int) -> str:
            if node_id < n:
                return labels[node_id]
            left, right, h = tree[node_id]
            left_str = build(left)
            right_str = build(right)
            bl_left = h - heights[left]
            bl_right = h - heights[right]
            left_part = f"{left_str}:{bl_left:.6f}"
            right_part = f"{right_str}:{bl_right:.6f}"
            return f"({left_part},{right_part})"
        newick = build(root) + ";"
        return newick

    def run(self) -> Tuple[Path, Path]:
        emb = self.load_embeddings()
        if not emb:
            raise RuntimeError("No embeddings found in HDF5.")
        selected_ids = self.choose_asvs(list(emb.keys()))
        vectors = []
        ids_present = []
        missing = []
        for asv in selected_ids:
            v = emb.get(asv)
            if v is None:
                missing.append(asv)
                continue
            ids_present.append(asv)
            vectors.append(v)
        if missing:
            print(f"Warning: {len(missing)} ASV IDs referenced by inputs missing embeddings; they were skipped.")
        if not ids_present:
            raise RuntimeError("After skipping missing embeddings, no ASVs remain to build phylogeny.")
        X = np.stack(vectors, axis=0)
        n = X.shape[0]
        if n == 1:
            newick = f"{ids_present[0]};"
        else:
            dists = pdist(X, metric="euclidean")
            Z = linkage(dists, method="average")
            newick = self._linkage_to_newick(Z, ids_present)
        with open(self.tree_path, "w") as fh:
            fh.write(newick + "\n")
        table = pd.DataFrame({"ASV_ID": ids_present})
        # Load cluster info if available (try both tsv/csv)
        possible_cluster_files = [self.clusters_tsv, self.clusters_tsv.with_suffix(".csv")]
        cluster_df = None
        for cf in possible_cluster_files:
            if cf.exists():
                try:
                    sep = "\t" if cf.suffix == ".tsv" else ","
                    cluster_df = pd.read_csv(cf, sep=sep, dtype=str).set_index("ASV_ID")
                    break
                except Exception:
                    cluster_df = None
        if cluster_df is not None:
            table["cluster_id"] = table["ASV_ID"].map(lambda x: cluster_df.loc[x]["cluster_id"] if x in cluster_df.index and "cluster_id" in cluster_df.columns else pd.NA)
            if "umap_1" in cluster_df.columns and "umap_2" in cluster_df.columns:
                table["umap_1"] = table["ASV_ID"].map(lambda x: cluster_df.loc[x]["umap_1"] if x in cluster_df.index else pd.NA)
                table["umap_2"] = table["ASV_ID"].map(lambda x: cluster_df.loc[x]["umap_2"] if x in cluster_df.index else pd.NA)
        else:
            table["cluster_id"] = pd.NA

        if self.novelty_tsv.exists():
            try:
                ndf = pd.read_csv(self.novelty_tsv, sep="\t", dtype=str).set_index("ASV_ID")
                table["novelty_score"] = table["ASV_ID"].map(lambda x: ndf.loc[x]["novelty_score"] if x in ndf.index and "novelty_score" in ndf.columns else pd.NA)
            except Exception:
                table["novelty_score"] = pd.NA
        else:
            table["novelty_score"] = pd.NA

        table["embedding_dim"] = table["ASV_ID"].map(lambda x: int(len(emb[x])) if x in emb else pd.NA)
        table.to_csv(self.table_path, sep="\t", index=False)
        return (self.tree_path, self.table_path)