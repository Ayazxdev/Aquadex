import h5py
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict

class NoveltyScorer:
    def __init__(self, ref_h5_path: Path):
        self.ref_h5_path = ref_h5_path
        # Load reference embeddings once
        with h5py.File(self.ref_h5_path, "r") as hf:
            self.ref_embeddings = {k: np.array(hf[k]) for k in hf.keys()}

    def run(self, anchor_matches_path: Path) -> pd.DataFrame:
        # Load anchor matches
        df = pd.read_csv(anchor_matches_path, sep="\t")
        novelty_scores = []

        for _, row in df.iterrows():
            distances = np.array([float(x) for x in row["distances"].split(";")])
            # Compute novelty as 1 - max similarity (closest ref)
            novelty = 1.0 - np.exp(-distances.min())
            novelty_scores.append(novelty)

        df["novelty_score"] = novelty_scores
        return df