# results_aggregator.py
# backend/app/results_aggregator.py
import os
import json
import re
import warnings
from pathlib import Path
from typing import Dict, Any, List, Optional
from datetime import datetime

# Suppress numpy MINGW-W64 RuntimeWarnings that crash Uvicorn on Windows Python 3.13
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", message=".*MINGW.*")
warnings.filterwarnings("ignore", message=".*invalid value encountered.*")

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform

# Config: where frontend can download static artifacts from
ARTIFACT_BASE_URL = os.getenv("ARTIFACT_BASE_URL", "/artifact").rstrip("/")

# Helper utilities
def safe_load_json(p: Path) -> Optional[Dict[str, Any]]:
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return None

def make_artifact_url(run_dir: Path, rel_path_str: str) -> str:
    """
    Build a frontend-friendly URL for an artifact. rel_path_str can be:
      - a path that already contains the run_id prefix (status.json style), or
      - relative path inside the run_dir (like 'qc/fastp_report.html').
    Returns: f"{ARTIFACT_BASE_URL}/{run_id}/{relpath_posix}"
    """
    run_id = run_dir.name
    # If status.json gives something like "runid/qc/fastp_report.html", strip runid prefix
    # and use the remaining path.
    if isinstance(rel_path_str, str) and rel_path_str:
        # Normalize separators
        relp = rel_path_str.replace("\\", "/")
        # If the string contains run_id, strip the prefix up to & including run_id/
        if f"{run_id}/" in relp:
            relp = relp.split(f"{run_id}/", 1)[1]
        # If the input is an absolute path, try to compute relative to run_dir
        try:
            p = Path(relp)
            if p.is_absolute():
                try:
                    relp = str(p.relative_to(run_dir)).replace("\\", "/")
                except Exception:
                    # fallback: use basename
                    relp = p.name
        except Exception:
            pass
        # final URL
        return f"{ARTIFACT_BASE_URL}/{run_id}/{relp.lstrip('/')}"
    return f"{ARTIFACT_BASE_URL}/{run_id}"


# Parsers and report extractors
def parse_qc(qc_dir: Path) -> Dict[str, Any]:
    """Read fastp/multiqc JSON if present. Checks both run-level and per-sample filenames."""
    qc_summary = {}
    
    # Try fastp_report.json (primary), then per-sample names
    candidates = [
        qc_dir / "fastp_report.json",
    ] + list(qc_dir.glob("*_fastp.json")) + list(qc_dir.glob("*fastp*.json"))
    
    for fastp_json in candidates:
        if fastp_json.exists() and fastp_json.stat().st_size > 10:
            try:
                data = json.loads(fastp_json.read_text(encoding="utf-8"))
                # Normalize: ensure top-level 'summary' key exists
                if "summary" in data:
                    qc_summary["fastp"] = data
                elif isinstance(data, dict) and any(k in data for k in ("before_filtering", "after_filtering")):
                    qc_summary["fastp"] = {"summary": data}
                else:
                    qc_summary["fastp"] = data
                break
            except Exception:
                continue
    
    # Include available HTMLs as URLs (frontend preview)
    qc_summary["fastp_html"] = str(qc_dir / "fastp_report.html") if (qc_dir / "fastp_report.html").exists() else None
    qc_summary["multiqc_html"] = str(qc_dir / "multiqc_report.html") if (qc_dir / "multiqc_report.html").exists() else None
    return qc_summary

def parse_novelty(novelty_file: Path) -> pd.DataFrame:
    if novelty_file.exists():
        try:
            return pd.read_csv(novelty_file, sep="\t", dtype=str).fillna("")
        except Exception:
            return pd.read_csv(novelty_file, sep=None, engine="python", dtype=str).fillna("")
    return pd.DataFrame()

def novelty_summary(novelty_df: pd.DataFrame) -> Dict[str, Any]:
    if novelty_df.empty:
        return {}
    # try known column names, cast to float when possible
    score_cols = [c for c in novelty_df.columns if "novelty" in c.lower() or "score" in c.lower() or "anchor_confidence" in c.lower()]
    scores = []
    for c in score_cols:
        try:
            scores = novelty_df[c].astype(float).values
            break
        except Exception:
            continue
    if len(scores) == 0:
        return {"total_asvs": len(novelty_df)}
    avg = float(np.nanmean(scores))
    high = int((scores > 0.5).sum())
    return {"avg_novelty_score": avg, "num_high_novel": high, "total_asvs": int(len(scores)), "pct_high_novel": round(100 * high / len(scores), 2)}

def parse_kraken2_report(report_path: Path, sample_name: str) -> pd.DataFrame:
    """
    Parse a Kraken2 .report file.
    Format: pct  reads_clade  reads_taxon  rank  taxid  name
    """
    rows = []
    # Map all standard and extended SILVA/Eukaryotic rank codes
    rank_map = {
        "D": "Domain", "K": "Kingdom", "K1": "Subkingdom", "P": "Phylum", "P1": "Subphylum",
        "C": "Class", "C1": "Subclass", "O": "Order", "O1": "Suborder", "F": "Family", 
        "F1": "Subfamily", "G": "Genus", "G1": "Subgenus", "S": "Species", "S1": "Subspecies",
        "R1": "Supergroup", "R2": "Division", "R3": "Subdivision", "R4": "Infraphylum"
    }
    try:
        with open(report_path, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 6:
                    continue
                try:
                    pct = float(parts[0].strip())
                    rank_code = parts[3].strip()
                    taxon_name = parts[5].strip()
                except ValueError:
                    continue
                
                # Filter out noise (entries with near-zero abundance)
                if pct < 0.01:
                    continue

                if rank_code == "U":
                    rows.append({"sample": sample_name, "taxon": "Unclassified", "abundance": pct / 100.0})
                elif rank_code in rank_map:
                    rank_label = rank_map[rank_code]
                    clean_name = f"{rank_label[0].lower()}__{taxon_name}"
                    rows.append({"sample": sample_name, "taxon": clean_name, "abundance": pct / 100.0})
    except Exception:
        return pd.DataFrame()

    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def parse_taxonomy(taxonomy_dir: Path) -> pd.DataFrame:
    """
    Read taxonomy data. Tries Kraken2 .report files first (hierarchical),
    then falls back to generic TSV parsing (including Kraken2 per-read TSV with size= fields).
    """
    if not taxonomy_dir.exists():
        return pd.DataFrame(columns=["sample", "taxon", "abundance"])

    # First try: Kraken2 .report files (summarized, multi-rank)
    reports = sorted(taxonomy_dir.glob("*.report"))
    if reports:
        frames = []
        for rpt in reports:
            m = re.match(r"(.+?)_taxonomy", rpt.stem)
            sample = m.group(1) if m else rpt.stem
            df = parse_kraken2_report(rpt, sample)
            if not df.empty:
                frames.append(df)
        if frames:
            combined = pd.concat(frames, ignore_index=True)
            # Normalize per sample
            out = []
            for sample, grp in combined.groupby("sample"):
                ssum = grp["abundance"].sum()
                if ssum > 0:
                    grp = grp.copy()
                    grp["abundance"] = grp["abundance"] / float(ssum)
                out.append(grp)
            if out:
                return pd.concat(out, ignore_index=True)

    # Fallback: generic TSV files (including Kraken2 per-read output)
    rows = []
    tsvs = sorted(taxonomy_dir.glob("*.tsv"))
    if not tsvs:
        f = taxonomy_dir / "taxonomy.tsv"
        if f.exists():
            tsvs = [f]

    for f in tsvs:
        m = re.match(r"(.+?)_taxonomy", f.stem)
        sample = m.group(1) if m else f.stem

        # Try to detect Kraken2 per-read output format: status\tasv_id\ttaxid\tlength\tkmer_hits
        try:
            first_line = f.open("r", encoding="utf-8", errors="replace").readline()
            is_kraken2_perread = (
                first_line.startswith("U\t") or first_line.startswith("C\t")
            ) and "size=" in first_line
        except Exception:
            is_kraken2_perread = False

        if is_kraken2_perread:
            # Parse Kraken2 per-read TSV: extract taxon and size for diversity
            total_classified = 0
            total_unclassified = 0
            taxon_counts = {}
            total_reads = 0
            try:
                with open(f, "r", encoding="utf-8", errors="replace") as fh:
                    for line in fh:
                        parts = line.strip().split("\t")
                        if len(parts) < 3:
                            continue
                        status = parts[0].strip()  # U or C
                        asv_raw = parts[1].strip()  # e.g. SRR.35;size=6134
                        tax_id = parts[2].strip()   # taxid (0 = unclassified)
                        # Extract abundance from size= field
                        size = 1
                        if "size=" in asv_raw:
                            try:
                                size = int(asv_raw.split("size=")[1].split(";")[0])
                            except Exception:
                                size = 1
                        total_reads += size
                        if status == "U" or tax_id == "0":
                            total_unclassified += size
                        else:
                            total_classified += size
                            taxon_key = f"TaxID:{tax_id}"
                            taxon_counts[taxon_key] = taxon_counts.get(taxon_key, 0) + size
            except Exception as e:
                logger.warning(f"Error parsing Kraken2 per-read TSV {f}: {e}")
                continue

            if total_reads == 0:
                continue

            # Build rows: each unique taxon is a row
            if taxon_counts:
                for taxon, count in taxon_counts.items():
                    rows.append({"sample": sample, "taxon": taxon, "abundance": float(count)})
            # Add unclassified as a taxon
            if total_unclassified > 0:
                rows.append({"sample": sample, "taxon": "Unclassified", "abundance": float(total_unclassified)})

            # If entirely unclassified, add placeholder so sample is detected
            if not taxon_counts:
                rows.append({"sample": sample, "taxon": "Unclassified", "abundance": float(total_reads)})
            continue

        # Generic TSV parsing (for demo/structured taxonomy TSVs)
        try:
            df = pd.read_csv(f, sep="\t", comment="#", engine="python", dtype=str).fillna("")
        except Exception:
            continue
        # infer sample name from column if present
        for c in ("sample", "sample_id", "run_id"):
            if c in df.columns:
                sample = str(df[c].iat[0]) if len(df) > 0 else sample
                break

        # find tax/abundance columns
        tax_col = next((c for c in df.columns if c.lower() in ("taxon","taxonomy","taxa","classification","name","asv_id","asv")), None)
        abundance_col = next((c for c in df.columns if c.lower() in ("abundance","rel_abundance","relative_abundance","count","reads","read_count","n")), None)

        if tax_col is None and "ASV_ID" in df.columns:
            tax_col = "ASV_ID"

        # For demo taxonomy TSVs: prefer the 'taxon' column as the taxonomy path
        if "taxon" in df.columns:
            tax_col = "taxon"

        if abundance_col is not None:
            df[abundance_col] = pd.to_numeric(df[abundance_col], errors="coerce").fillna(0.0)
            total = float(df[abundance_col].sum())
            if total <= 0:
                for _, r in df.iterrows():
                    tax = str(r.get(tax_col, "") if tax_col else "")
                    tax = tax.strip() or "Unclassified"
                    rows.append({"sample": sample, "taxon": tax, "abundance": 0.0})
            else:
                for _, r in df.iterrows():
                    tax = str(r.get(tax_col, "") if tax_col else "")
                    tax = tax.strip() or "Unclassified"
                    rows.append({"sample": sample, "taxon": tax, "abundance": float(r[abundance_col])})
        else:
            total = max(1, len(df))
            for _, r in df.iterrows():
                tax = str(r.get(tax_col, "") if tax_col else "")
                tax = tax.strip() or "Unclassified"
                rows.append({"sample": sample, "taxon": tax, "abundance": 1.0 / total})

    if not rows:
        return pd.DataFrame([{"sample": "sample1", "taxon": "Unclassified", "abundance": 1.0}])
    out_df = pd.DataFrame(rows)
    out_df = out_df.groupby(["sample", "taxon"], as_index=False)["abundance"].sum()
    out = []
    for sample, grp in out_df.groupby("sample"):
        ssum = grp["abundance"].sum()
        if ssum <= 0:
            grp["abundance"] = grp["abundance"]
        else:
            grp = grp.copy()
            grp["abundance"] = grp["abundance"] / float(ssum)
        out.append(grp)
    if out:
        out_df = pd.concat(out, ignore_index=True)
    return out_df


def limit_taxonomy_rows(df: pd.DataFrame, top_n: int = 25) -> pd.DataFrame:
    if df.empty:
        return df
    return df.sort_values("abundance", ascending=False).groupby("sample").head(top_n)

# ─────────────────────────────────────────────────────────────────────────────
# ALPHA DIVERSITY — Hill numbers, Pielou evenness, bootstrap 95% CI
# ─────────────────────────────────────────────────────────────────────────────
# Hill diversity of order q:
#   q=0 →  D = S             (species richness)
#   q=1 →  D = exp(H')       (exponential Shannon)
#   q=2 →  D = 1/Σp²         (inverse Simpson)
# Pielou evenness: J = H' / ln(S),  0 ≤ J ≤ 1
# Bootstrap CI uses 999 non-parametric resamples of the observed abundance vector.

def _hill(p: np.ndarray, q: float) -> float:
    """True Hill number of order q for a relative-abundance vector p (must sum to 1)."""
    p = p[p > 0]
    if q == 1.0:
        return float(np.exp(-np.sum(p * np.log(p))))
    elif q == 0.0:
        return float(len(p))
    else:
        return float(np.sum(p ** q) ** (1.0 / (1.0 - q)))


def _bootstrap_ci(
    p: np.ndarray,
    q: float,
    n_boot: int = 999,
    alpha: float = 0.05,
    rng_seed: int = 0,
) -> tuple:
    """Return (lower, upper) percentile bootstrap CI for Hill(q)."""
    rng = np.random.default_rng(rng_seed)
    n = len(p)
    boot_vals = []
    # Multinomial resample — preserves sample size
    for _ in range(n_boot):
        counts = rng.multinomial(n, p / p.sum())
        p_boot = counts.astype(float) / counts.sum()
        boot_vals.append(_hill(p_boot, q))
    lo = float(np.percentile(boot_vals, 100 * alpha / 2))
    hi = float(np.percentile(boot_vals, 100 * (1 - alpha / 2)))
    return lo, hi


def compute_alpha_diversity(taxonomy_df: pd.DataFrame) -> List[Dict[str, Any]]:
    """
    Compute per-sample Hill diversity (q=0,1,2), Pielou evenness J, Shannon H',
    and 95 % percentile bootstrap confidence intervals (999 resamples).

    Input: taxonomy_df with columns [sample, taxon, abundance] where
    parse_taxonomy() has already normalised each sample's abundance vector
    to relative proportions (sum == 1 per sample).
    """
    if taxonomy_df.empty or "abundance" not in taxonomy_df.columns or "sample" not in taxonomy_df.columns:
        return []
    out = []
    for sample, grp in taxonomy_df.groupby("sample"):
        p = grp["abundance"].values.astype(float)
        p = p[p > 0]          # drop structural zeros
        p = p / p.sum()       # enforce unit sum (guard against floating-point drift)

        S = int(len(p))
        if S == 0:
            out.append({"sample": str(sample), "richness": 0, "shannon": 0.0,
                        "simpson": 0.0, "hill_q0": 0, "hill_q1": 0.0, "hill_q2": 0.0,
                        "pielou_j": None, "shannon_ci": [None, None],
                        "hill_q0_ci": [None, None], "hill_q1_ci": [None, None],
                        "hill_q2_ci": [None, None]})
            continue

        # ── Point estimates ──────────────────────────────────────────────────
        # Shannon entropy H' = -Σ p_i ln(p_i)
        H_prime = float(-np.sum(p * np.log(p)))
        # Gini-Simpson = 1 - Σp²
        gini_simpson = float(1.0 - np.sum(p ** 2))
        # Hill numbers (unified framework)
        D0 = float(S)                        # q=0: richness
        D1 = float(np.exp(H_prime))          # q=1: exp(Shannon)
        D2 = float(1.0 / np.sum(p ** 2))     # q=2: inverse-Simpson
        # Pielou evenness
        pielou_j = round(H_prime / np.log(S), 4) if S > 1 else 1.0

        # ── Bootstrap 95 % CI ────────────────────────────────────────────────
        # We resample with n = len(p) pseudocounts → multinomial bootstrap
        lo0, hi0 = _bootstrap_ci(p, 0.0)
        lo1, hi1 = _bootstrap_ci(p, 1.0)
        lo2, hi2 = _bootstrap_ci(p, 2.0)
        # Shannon CI via q→1 Hill CI converted back to H'
        lo_h = round(np.log(lo1), 4) if lo1 > 0 else 0.0
        hi_h = round(np.log(hi1), 4) if hi1 > 0 else 0.0

        out.append({
            "sample":      str(sample),
            # Legacy fields (kept for backwards-compat with frontend)
            "richness":    S,
            "shannon":     round(H_prime, 4),
            "simpson":     round(gini_simpson, 4),
            # Hill diversity series
            "hill_q0":     S,
            "hill_q1":     round(D1, 4),
            "hill_q2":     round(D2, 4),
            # Pielou evenness
            "pielou_j":    pielou_j,
            # 95 % bootstrap CIs
            "shannon_ci":  [lo_h, hi_h],
            "hill_q0_ci": [int(lo0), int(hi0)],
            "hill_q1_ci": [round(lo1, 4), round(hi1, 4)],
            "hill_q2_ci": [round(lo2, 4), round(hi2, 4)],
            # Methodology provenance
            "method": {
                "hill_formula": "^q D = (Σ p_i^q)^(1/(1-q)); q=1: exp(-Σ p_i ln p_i)",
                "pielou_formula": "J = H' / ln(S)",
                "ci_method": "nonparametric multinomial bootstrap",
                "ci_resamples": 999,
                "ci_level": "95%",
            },
        })
    return out


# ─────────────────────────────────────────────────────────────────────────────
# CLASSICAL PCoA — Gower (1966) eigendecomposition, NOT metric MDS
# ─────────────────────────────────────────────────────────────────────────────
def _classical_pcoa(D: np.ndarray, labels: List[str]) -> Dict[str, Any]:
    """
    Classical Principal Coordinates Analysis (Gower 1966).
    Input D: symmetric, non-negative distance matrix (n×n).
    Returns dict with PC1/PC2 coordinates, eigenvalues, and proportion of
    variation explained — equivalent to what scikit-bio pcoa() reports.

    Algorithm:
      1. Double-centre: A = -0.5 * D²
         B = (I - 11ᵀ/n) A (I - 11ᵀ/n)
      2. Eigendecompose B = VΛVᵀ
      3. Coordinates: F_i = V_i * √λ_i  (for positive eigenvalues only)
    """
    n = D.shape[0]
    D2 = D ** 2
    # Double centering
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ D2 @ J
    # Symmetrize to eliminate floating-point asymmetry
    B = (B + B.T) / 2.0
    # Eigendecomposition (eigh is stable for symmetric matrices)
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    # Keep positive eigenvalues only for coordinates
    pos_mask = eigenvalues > 1e-10
    pos_vals = eigenvalues[pos_mask]
    pos_vecs = eigenvectors[:, pos_mask]
    if pos_vals.shape[0] < 2:
        # Degenerate: return zeros
        return {"pcoa": [{"sample": s, "x": 0.0, "y": 0.0, "pc1_var": 0.0, "pc2_var": 0.0} for s in labels],
                "eigenvalues": [], "proportion_explained": []}
    coords = pos_vecs * np.sqrt(pos_vals)   # shape (n, n_pos)
    # Proportion of variation explained
    total_pos = pos_vals.sum()
    prop_explained = pos_vals / total_pos
    pc1_var = round(float(prop_explained[0]) * 100, 1)
    pc2_var = round(float(prop_explained[1]) * 100, 1) if len(prop_explained) > 1 else 0.0
    pcoa_points = []
    for i, s in enumerate(labels):
        pcoa_points.append({
            "sample": s,
            "x": float(coords[i, 0]),
            "y": float(coords[i, 1]) if coords.shape[1] > 1 else 0.0,
            "pc1_var": pc1_var,
            "pc2_var": pc2_var,
        })
    return {
        "pcoa": pcoa_points,
        "eigenvalues": [round(float(v), 6) for v in eigenvalues[:4].tolist()],
        "proportion_explained": [
            {"axis": f"PC{i+1}", "pct": round(float(prop_explained[i]) * 100, 2)}
            for i in range(min(4, len(prop_explained)))
        ],
    }


# ─────────────────────────────────────────────────────────────────────────────
# PERMANOVA — Anderson (2001) permutation test on distance matrix
# ─────────────────────────────────────────────────────────────────────────────
def _permanova(
    D: np.ndarray,
    group_labels: List[str],
    n_perm: int = 999,
    rng_seed: int = 42,
) -> Dict[str, Any]:
    """
    One-way PERMANOVA (Anderson 2001 / adonis).
    Tests whether mean within-group distance < mean between-group distance.

    Pseudo-F statistic:
        SS_total = Σ_{i<j} d²_{ij} / n
        SS_within = Σ_g Σ_{i<j∈g} d²_{ij} / n_g
        SS_between = SS_total - SS_within

        F = (SS_between / (g-1)) / (SS_within / (n-g))
        R² = SS_between / SS_total

    p-value from permutation distribution of F under H₀.
    """
    if len(set(group_labels)) < 2:
        return {}
    n = D.shape[0]
    D2 = (D ** 2) / n
    groups = sorted(set(group_labels))
    g = len(groups)
    indices = {grp: [i for i, lbl in enumerate(group_labels) if lbl == grp] for grp in groups}

    def pseudo_f(perm_labels):
        perm_indices = {grp: [i for i, lbl in enumerate(perm_labels) if lbl == grp] for grp in groups}
        ss_w = sum(
            np.sum(D2[np.ix_(idx, idx)]) / 2.0
            for idx in perm_indices.values() if len(idx) > 1
        )
        ss_t = np.sum(np.triu(D2, k=1))
        ss_b = ss_t - ss_w
        n_w = n - g
        n_b = g - 1
        if ss_w <= 0 or n_w <= 0:
            return 0.0
        return (ss_b / n_b) / (ss_w / n_w)

    obs_f = pseudo_f(group_labels)
    ss_t = np.sum(np.triu(D2, k=1))
    ss_w = sum(
        np.sum(D2[np.ix_(idx, idx)]) / 2.0
        for idx in indices.values() if len(idx) > 1
    )
    r_squared = float((ss_t - ss_w) / ss_t) if ss_t > 0 else 0.0

    rng = np.random.default_rng(rng_seed)
    labels_arr = np.array(group_labels)
    exceed = 0
    for _ in range(n_perm):
        perm = rng.permutation(labels_arr).tolist()
        if pseudo_f(perm) >= obs_f:
            exceed += 1
    p_value = (exceed + 1) / (n_perm + 1)   # + 1 for observed test stat

    return {
        "pseudo_f": round(float(obs_f), 4),
        "r_squared": round(r_squared, 4),
        "p_value": round(p_value, 4),
        "n_permutations": n_perm,
        "groups": groups,
        "n_samples": n,
        "method": "PERMANOVA (Anderson 2001); pseudo-F permutation test",
    }


# ─────────────────────────────────────────────────────────────────────────────
# BETA DIVERSITY — Bray-Curtis PCoA, Aitchison CLR PCoA, Jaccard, PERMANOVA
# ─────────────────────────────────────────────────────────────────────────────
# Aitchison preprocessing:
#   zero replacement  = multiplicative pseudocount
#   pseudocount value = 0.5 / (number of taxa)   [GBM-style half-count replacement]
#   This is applied before CLR so that provenance is explicit.

def compute_beta_diversity(taxonomy_df: pd.DataFrame) -> Dict[str, Any]:
    """
    Compute all beta-diversity metrics:
      - Bray-Curtis dissimilarity + classical PCoA (eigendecomposition)
      - Jaccard dissimilarity (presence/absence)
      - Aitchison distance (CLR-transformed Euclidean) + classical PCoA
      - PERMANOVA: one-way test treating each sample as its own group

    Input: taxonomy_df with columns [sample, taxon, abundance] where
    abundances are per-sample relative proportions (as returned by parse_taxonomy).
    """
    if taxonomy_df.empty or "sample" not in taxonomy_df.columns or "taxon" not in taxonomy_df.columns:
        return {}

    abundance_df = taxonomy_df.pivot_table(
        index="sample", columns="taxon", values="abundance",
        aggfunc="sum", fill_value=0.0
    )
    abundance_df = abundance_df.astype(float)

    # Re-normalise each row to ensure Σp_i = 1 per sample (guard against pivot fill)
    row_sums = abundance_df.sum(axis=1)
    abundance_df = abundance_df.div(row_sums.replace(0, 1), axis=0)

    n_samples = abundance_df.shape[0]
    if n_samples < 2:
        return {}   # Bray-Curtis / PCoA require ≥ 2 samples

    samples = list(abundance_df.index)
    X = abundance_df.values   # shape (n_samples, n_taxa) — row-normalised relative abundances

    # ── Bray-Curtis ──────────────────────────────────────────────────────────
    # BC(x,y) = Σ|x_i - y_i| / Σ(x_i + y_i)  on relative-abundance vectors
    bc_vec = pdist(X, metric="braycurtis")
    bc_matrix = squareform(bc_vec)
    bc_distances = [
        {"sample1": samples[i], "sample2": samples[j], "value": round(float(bc_matrix[i, j]), 6)}
        for i in range(n_samples) for j in range(n_samples)
    ]
    bc_pcoa = _classical_pcoa(bc_matrix, samples)

    # ── Jaccard (presence/absence) ───────────────────────────────────────────
    # J(x,y) = 1 - |A∩B| / |A∪B|   where A,B are sets of observed taxa
    X_pa = (X > 0).astype(float)
    jac_vec = pdist(X_pa, metric="jaccard")
    jac_matrix = squareform(jac_vec)
    jac_distances = [
        {"sample1": samples[i], "sample2": samples[j], "value": round(float(jac_matrix[i, j]), 6)}
        for i in range(n_samples) for j in range(n_samples)
    ]

    # ── Aitchison distance (CLR-based) ───────────────────────────────────────
    # Preprocessing:
    #   pseudocount = 0.5 / n_taxa  (multiplicative — GBM half-count replacement)
    #   This preserves compositional ratios while avoiding log(0).
    # CLR: clr(p_i) = ln(p_i / g(p))  where g(p) = geometric mean of p
    # Aitchison distance = Euclidean distance in CLR space
    n_taxa = X.shape[1]
    PSEUDOCOUNT = 0.5 / n_taxa   # explicitly declared — not silently added
    X_pseudo = X + PSEUDOCOUNT
    X_pseudo = X_pseudo / X_pseudo.sum(axis=1, keepdims=True)   # re-close after pseudocount
    log_X = np.log(X_pseudo)
    clr_X = log_X - log_X.mean(axis=1, keepdims=True)   # CLR transform
    ait_vec = pdist(clr_X, metric="euclidean")
    ait_matrix = squareform(ait_vec)
    ait_distances = [
        {"sample1": samples[i], "sample2": samples[j], "value": round(float(ait_matrix[i, j]), 6)}
        for i in range(n_samples) for j in range(n_samples)
    ]
    ait_pcoa = _classical_pcoa(ait_matrix, samples)

    # ── PERMANOVA on Bray-Curtis ─────────────────────────────────────────────
    # With each sample as its own group (label = sample name), this tests
    # whether samples differ significantly overall.
    permanova = {}
    if n_samples >= 3:
        permanova = _permanova(bc_matrix, samples, n_perm=999)

    return {
        "samples": samples,
        "normalization_note": (
            "Abundance vectors are per-sample relative proportions (row-sum = 1). "
            "Bray-Curtis and Jaccard operate on these proportions directly. "
            f"Aitchison CLR uses pseudocount = {PSEUDOCOUNT:.2e} (0.5 / n_taxa = 0.5 / {n_taxa})."
        ),
        "bray_curtis": {
            "distances": bc_distances,
            "pcoa": bc_pcoa["pcoa"],
            "proportion_explained": bc_pcoa.get("proportion_explained", []),
            "eigenvalues": bc_pcoa.get("eigenvalues", []),
            "method": "Bray-Curtis dissimilarity; classical PCoA (Gower 1966 eigendecomposition)",
        },
        "jaccard": {
            "distances": jac_distances,
            "method": "Jaccard dissimilarity on presence/absence (0/1) vectors",
        },
        "aitchison": {
            "distances": ait_distances,
            "pcoa": ait_pcoa["pcoa"],
            "proportion_explained": ait_pcoa.get("proportion_explained", []),
            "eigenvalues": ait_pcoa.get("eigenvalues", []),
            "method": "Aitchison distance = Euclidean in CLR space",
            "zero_handling": f"multiplicative pseudocount = {PSEUDOCOUNT:.4e} (0.5 / n_taxa = 0.5 / {n_taxa})",
        },
        "permanova": permanova,
        # Legacy field: kept for backwards compat — points to Bray-Curtis PCoA
        "pcoa": bc_pcoa["pcoa"],
        "distances": bc_distances,
    }

def parse_clustering(clustering_dir: Path) -> pd.DataFrame:
    if not clustering_dir.exists():
        return pd.DataFrame()
    candidates = [clustering_dir / "clusters.csv", clustering_dir / "clusters.tsv"]
    found = None
    for c in candidates:
        if c.exists():
            found = c
            break
    if found is None:
        for ext in ("*.csv", "*.tsv", "*.txt"):
            files = list(clustering_dir.glob(ext))
            if files:
                found = files[0]
                break
    if found is None:
        return pd.DataFrame()
    try:
        if found.suffix.lower() == ".csv":
            df = pd.read_csv(found)
        else:
            df = pd.read_csv(found, sep="\t", engine="python")
    except Exception:
        try:
            df = pd.read_csv(found, sep=None, engine="python")
        except Exception:
            return pd.DataFrame()
    return df

def _normalize_clustering_row(r: Dict[str, Any]) -> Dict[str, Any]:
    asv = r.get("ASV_ID") or r.get("asv") or r.get("ASV") or r.get("asv_id") or ""
    cluster = r.get("cluster_id") if r.get("cluster_id") is not None else r.get("cluster") if r.get("cluster") is not None else r.get("clusterId") if r.get("clusterId") is not None else None
    x = r.get("dim_1") or r.get("umap_1") or r.get("x") or r.get("dim1")
    y = r.get("dim_2") or r.get("umap_2") or r.get("y") or r.get("dim2")
    try:
        x = float(x) if x is not None else None
        y = float(y) if y is not None else None
    except Exception:
        x = None; y = None
    novelty_score = r.get("novelty_score") or r.get("score") or r.get("novelty") or r.get("anchor_confidence")
    try:
        novelty_score = float(novelty_score) if novelty_score is not None and str(novelty_score) != "" else None
    except Exception:
        novelty_score = None
    top_ref = r.get("top_refs") or r.get("top_ref") or r.get("closest_ref") or ""
    return {"asv": str(asv), "cluster": int(cluster) if cluster is not None and str(cluster) != "" else None, "umap": [x, y], "novelty_score": novelty_score, "top_ref": top_ref}

def build_sankey_from_taxonomy(taxonomy_df: pd.DataFrame, top_k: int = 25) -> Dict[str, Any]:
    """
    Build hierarchical Sankey nodes and links suitable for Nivo ResponsiveSankey.
    Handles both clean rank prefixes (d__, p__, c__, etc.) and semicolon paths.
    For all-unclassified data, returns a meaningful single-node representation.
    """
    if taxonomy_df.empty:
        return {"nodes": [{"id": "Unclassified"}], "links": []}

    # Check if all data is unclassified
    non_unclassified = taxonomy_df[taxonomy_df["taxon"].str.lower() != "unclassified"]
    if non_unclassified.empty:
        # Return a meaningful fallback for all-unclassified data
        total = taxonomy_df["abundance"].sum()
        samples = taxonomy_df["sample"].unique().tolist()
        nodes = [{"id": "All Reads"}]
        links = []
        for s in samples:
            sample_total = taxonomy_df[taxonomy_df["sample"] == s]["abundance"].sum()
            sample_pct = float(sample_total / total * 100) if total > 0 else 0
            node_id = f"Unclassified ({s})"
            nodes.append({"id": node_id})
            links.append({"source": "All Reads", "target": node_id, "value": max(1, round(sample_pct, 1))})
        if not links:
            nodes = [{"id": "Unclassified"}]
        return {"nodes": nodes, "links": links}

    # Use only classified taxa for Sankey
    taxonomy_df = non_unclassified.copy()

    # If taxa are semicolon-delimited (demo format: k__Bacteria;p__Proteobacteria...)
    has_semicolons = taxonomy_df["taxon"].str.contains(";").any()
    if has_semicolons:
        agg = taxonomy_df.groupby("taxon", as_index=False)["abundance"].sum().sort_values("abundance", ascending=False).head(top_k)
        nodes_set = set()
        links_dict = {}
        for _, row in agg.iterrows():
            parts = [p.strip() for p in str(row["taxon"]).split(";") if p.strip()]
            # Clean rank prefixes: k__Bacteria -> Bacteria
            clean_parts = []
            for p in parts:
                if "__" in p:
                    clean = p.split("__", 1)[1].strip()
                    if clean:
                        clean_parts.append(clean)
                else:
                    clean_parts.append(p)
            parts = clean_parts
            for p in parts:
                nodes_set.add(p)
            for i in range(len(parts) - 1):
                pair = (parts[i], parts[i + 1])
                links_dict[pair] = links_dict.get(pair, 0.0) + float(row["abundance"])
        
        nodes = [{"id": n} for n in nodes_set]
        links = [{"source": k[0], "target": k[1], "value": round(v * 100, 2) if v < 1.0 else round(v, 2)} for k, v in links_dict.items() if v > 0]
        if nodes and links:
            return {"nodes": nodes, "links": links}

    # Group ranks by prefix: d__ -> p__ -> c__ -> o__ -> f__ -> g__
    rank_order = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]
    ranked = {}
    unclassified_val = 0.0
    for _, row in taxonomy_df.iterrows():
        t = str(row["taxon"]).strip()
        ab = float(row.get("abundance", 0.0))
        if t.lower() == "unclassified":
            unclassified_val += ab
            continue
        for rk in rank_order:
            if t.startswith(rk):
                ranked.setdefault(rk, []).append((t, ab))
                break

    nodes_set = set()
    links_dict = {}
    
    # Connect Domain -> Phylum -> Class
    available_ranks = [rk for rk in rank_order if rk in ranked]
    for i in range(len(available_ranks) - 1):
        r1 = available_ranks[i]
        r2 = available_ranks[i + 1]
        top_from = sorted(ranked[r1], key=lambda x: x[1], reverse=True)[:3]
        top_to = sorted(ranked[r2], key=lambda x: x[1], reverse=True)[:8]
        for f_name, f_val in top_from:
            nodes_set.add(f_name)
            for t_name, t_val in top_to:
                nodes_set.add(t_name)
                val = max(round(t_val * 100, 2), 0.5)
                links_dict[(f_name, t_name)] = val

    nodes = [{"id": n} for n in nodes_set]
    links = [{"source": k[0], "target": k[1], "value": v} for k, v in links_dict.items()]
    
    if not links:
        return {"nodes": [{"id": "Unclassified"}], "links": []}
    return {"nodes": nodes, "links": links}


# Noise taxa to exclude from sunburst visualization
_NOISE_TAXA = {
    "unclassified", "unclassified sequences", "environmental samples",
    "uncultured eukaryote", "uncultured marine eukaryote",
    "uncultured marine diplonemid", "uncultured organism", "uncultured bacterium",
    "metagenome", "other sequences", "root", "cellular organisms",
    "opisthokonta", "eumetazoa", "sar", "discoba", "streptophytina", "dikarya",
    "saccharomyceta", "hydroidolina",
}


def _is_binomial_species(name: str) -> bool:
    """Returns True if name looks like a proper binomial species name (e.g. 'Tiaropsis multicirrata')."""
    parts = name.strip().split()
    if len(parts) < 2:
        return False
    # Reject names that start with 'uncultured', 'unclassified', 'environmental', etc.
    noise_words = {"uncultured", "unclassified", "environmental", "metagenome", "phytoplankton",
                   "eukaryote", "clone", "picoeukaryote", "microeukaryote", "picoplankton",
                   "marine", "freshwater", "rumen"}
    if parts[0].lower() in noise_words:
        return False
    # Genus starts with capital letter, species epithet is lowercase
    if parts[0][0].isupper() and parts[1][0].islower():
        return True
    return False


def build_sunburst_data(taxonomy_df: pd.DataFrame, top_n_per_rank: int = 15) -> Dict[str, Any]:
    """
    Build a Plotly sunburst-compatible hierarchy from taxonomy data.
    Returns {"labels": [...], "parents": [...], "values": [...], "text": [...]}
    Hierarchy: Life -> Domain/Kingdom -> Phylum -> Class -> Order -> Family -> Genus -> Species
    """
    if taxonomy_df is None or taxonomy_df.empty:
        return {"labels": [], "parents": [], "values": [], "text": []}

    # Filter noise taxa
    clean_df = taxonomy_df[~taxonomy_df["taxon"].astype(str).str.lower().isin(_NOISE_TAXA)].copy()
    if clean_df.empty:
        return {"labels": [], "parents": [], "values": [], "text": []}

    sample_taxa = clean_df["taxon"].dropna().astype(str).tolist()
    has_semicolon_lineages = any(";" in t for t in sample_taxa)

    if has_semicolon_lineages:
        # Build multi-ring hierarchical tree directly from full lineage paths
        node_parents: Dict[str, str] = {"Life": ""}
        node_values: Dict[str, float] = {"Life": 0.0}

        for _, row in clean_df.iterrows():
            t = str(row.get("taxon", "")).strip()
            if not t or t.lower() in _NOISE_TAXA or t.lower() == "unclassified":
                continue
            try:
                ab = float(row.get("abundance", 1.0))
            except (ValueError, TypeError):
                ab = 1.0

            chunks = [c.strip() for c in t.split(";") if c.strip()]
            current_parent = "Life"
            for chunk in chunks:
                clean_name = chunk.split("__", 1)[1].strip() if "__" in chunk else chunk
                if not clean_name or clean_name.lower() in _NOISE_TAXA or clean_name.lower() == "unclassified":
                    continue
                node_values[clean_name] = node_values.get(clean_name, 0.0) + ab
                if clean_name not in node_parents:
                    node_parents[clean_name] = current_parent
                current_parent = clean_name

        top_children = [k for k, p in node_parents.items() if p == "Life"]
        node_values["Life"] = sum(node_values.get(c, 0.0) for c in top_children) or 1.0

        labels = []
        parents = []
        values = []
        text = []
        total_ab = node_values["Life"] or 1.0

        for name in node_parents:
            val = node_values.get(name, 0.0)
            labels.append(name)
            parents.append(node_parents[name])
            values.append(round(val, 2))
            pct = round((val / total_ab) * 100, 1)
            text.append(f"{name}<br>{pct}% ({round(val):,} reads)")

        return {"labels": labels, "parents": parents, "values": values, "text": text}

    # Fallback for individual rank prefixes (d__, k__, p__, c__, o__, f__, g__, s__)
    ranks_ordered = ["d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__"]
    agg = (
        clean_df.groupby("taxon", as_index=False)["abundance"]
        .sum()
        .sort_values("abundance", ascending=False)
    )

    per_rank: Dict[str, list] = {rk: [] for rk in ranks_ordered}
    for _, row in agg.iterrows():
        t = str(row["taxon"]).strip()
        ab = float(row["abundance"])
        clean = t.split("__", 1)[1].strip() if "__" in t else t
        if clean.lower() in _NOISE_TAXA or not clean:
            continue
        for rk in ranks_ordered:
            if t.startswith(rk):
                if rk == "s__" and not _is_binomial_species(clean):
                    break
                per_rank[rk].append((clean, ab))
                break

    for rk in ranks_ordered:
        per_rank[rk] = sorted(per_rank[rk], key=lambda x: x[1], reverse=True)[:top_n_per_rank]

    labels = ["Life"]
    parents = [""]
    values = [0.0]

    domain_ranks = per_rank["d__"] + per_rank["k__"]
    for name, ab in domain_ranks:
        labels.append(name)
        parents.append("Life")
        values.append(round(ab * 100, 2))

    top_domain = domain_ranks[0][0] if domain_ranks else "Eukaryota"

    top_phyla = []
    for name, ab in per_rank["p__"]:
        labels.append(name)
        parents.append(top_domain)
        values.append(round(ab * 100, 2))
        top_phyla.append(name)
    top_phylum = top_phyla[0] if top_phyla else top_domain

    top_classes = []
    for i, (name, ab) in enumerate(per_rank["c__"]):
        parent = top_phyla[i % len(top_phyla)] if top_phyla else top_domain
        labels.append(name)
        parents.append(parent)
        values.append(round(ab * 100, 2))
        top_classes.append(name)

    top_orders = []
    for i, (name, ab) in enumerate(per_rank["o__"]):
        parent = top_classes[i % len(top_classes)] if top_classes else top_phylum
        labels.append(name)
        parents.append(parent)
        values.append(round(ab * 100, 2))
        top_orders.append(name)

    top_families = []
    for i, (name, ab) in enumerate(per_rank["f__"]):
        parent = top_orders[i % len(top_orders)] if top_orders else top_phylum
        labels.append(name)
        parents.append(parent)
        values.append(round(ab * 100, 2))
        top_families.append(name)

    top_genera = []
    for i, (name, ab) in enumerate(per_rank["g__"]):
        parent = top_families[i % len(top_families)] if top_families else top_phylum
        labels.append(name)
        parents.append(parent)
        values.append(round(ab * 100, 2))
        top_genera.append(name)

    for i, (name, ab) in enumerate(per_rank["s__"]):
        genus_name = name.split()[0] if name else ""
        parent = genus_name if genus_name in top_genera else (top_genera[i % len(top_genera)] if top_genera else top_phylum)
        labels.append(name)
        parents.append(parent)
        values.append(round(ab * 100, 2))

    total_val = sum(values[1:1 + len(domain_ranks)]) or 1.0
    values[0] = round(total_val, 2)
    text = [f"{lb}<br>{round(v / total_val * 100, 1)}%" for lb, v in zip(labels, values)]

    return {"labels": labels, "parents": parents, "values": values, "text": text}


# Aggregation logic
def parse_status_json(run_dir: Path) -> Optional[Dict[str, Any]]:
    candidate = run_dir / "status.json"
    if candidate.exists():
        try:
            return json.loads(candidate.read_text(encoding="utf-8"))
        except Exception:
            return None
    return None


def _append_sequences_to_novelty_tsv(run_dir: Path):
    """
    Appends an 'acgt_sequence' column to the novelty_report.tsv by looking up
    the ASV IDs in the corresponding dereplicated FASTA files.
    """
    nov_tsv = run_dir / "novelty" / "novelty_report.tsv"
    if not nov_tsv.exists():
        return

    # Check if we already appended sequences (idempotency)
    with open(nov_tsv, "r", encoding="utf-8") as fh:
        first_line = fh.readline().strip()
        if "acgt_sequence" in first_line:
            return

    # Find candidate FASTA files
    fasta_candidates = [
        run_dir / "denoise" / "sample1_R1.derep.fasta",
        run_dir / "denoise" / "sample1_R1.cleaned.fasta",
    ] + list((run_dir / "denoise").glob("*.fasta")) if (run_dir / "denoise").exists() else []

    # Build an in-memory dictionary of ASV -> sequence
    seq_map = {}
    for f in fasta_candidates:
        if not f.exists():
            continue
        curr_id = None
        curr_seq = []
        with open(f, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    if curr_id is not None:
                        seq_map[curr_id] = "".join(curr_seq)
                    # Normalize header: SRR...;size=1 -> SRR..._size=1
                    curr_id = line[1:].replace(";size=", "_size=").split()[0]
                    curr_seq = []
                else:
                    curr_seq.append(line)
            if curr_id is not None:
                seq_map[curr_id] = "".join(curr_seq)
        if seq_map:
            break

    if not seq_map:
        return

    # Rewrite the TSV with the new column
    import csv
    temp_tsv = nov_tsv.with_suffix(".tsv.tmp")
    with open(nov_tsv, "r", encoding="utf-8") as in_fh, \
         open(temp_tsv, "w", encoding="utf-8", newline="") as out_fh:
        reader = csv.reader(in_fh, delimiter="\t")
        writer = csv.writer(out_fh, delimiter="\t")
        
        try:
            header = next(reader)
            # Find the ID column (could be ASV_ID, asv_id, asv, or id)
            id_col_idx = -1
            for i, h in enumerate(header):
                if h.lower() in ("asv_id", "asv", "id"):
                    id_col_idx = i
                    break
            
            header.append("acgt_sequence")
            writer.writerow(header)
            
            for row in reader:
                seq = ""
                if id_col_idx >= 0 and id_col_idx < len(row):
                    asv_id = row[id_col_idx].replace(";", "_")
                    seq = seq_map.get(asv_id, "")
                row.append(seq)
                writer.writerow(row)
        except StopIteration:
            pass

    temp_tsv.replace(nov_tsv)


def aggregate_results(run_dir: Path) -> Dict[str, Any]:
    """Main aggregator: returns a structured summary dict for the frontend."""
    _append_sequences_to_novelty_tsv(run_dir)
    # parse
    qc = parse_qc(run_dir / "qc")
    novelty_df = parse_novelty(run_dir / "novelty" / "novelty_report.tsv")
    taxonomy_df = parse_taxonomy(run_dir / "taxonomy")
    clustering_df = parse_clustering(run_dir / "clustering")

    # normal computations on full datasets
    alpha_div = compute_alpha_diversity(taxonomy_df)
    beta_div = compute_beta_diversity(taxonomy_df)
    novelty_stats = novelty_summary(novelty_df)

    # Now limit rows for frontend serialization
    taxonomy_df_limited = limit_taxonomy_rows(taxonomy_df, top_n=50)

    clustering_list = []
    if not clustering_df.empty:
        # Sample points if too large to prevent browser freezing on scatterplot
        if len(clustering_df) > 1500:
            sample_df = clustering_df.sample(n=1500, random_state=42)
        else:
            sample_df = clustering_df
        for _, row in sample_df.iterrows():
            clustering_list.append(_normalize_clustering_row(row.to_dict()))

    # cluster stats
    cluster_stats = {}
    if not clustering_df.empty:
        if "cluster_id" in clustering_df.columns:
            sizes = clustering_df["cluster_id"].value_counts().to_dict()
            cluster_stats = {"sizes": {str(k): int(v) for k, v in sizes.items()}, "num_points": int(len(clustering_df))}
        elif clustering_list:
            dfc = pd.DataFrame(clustering_list)
            if "cluster" in dfc.columns:
                sizes = dfc.groupby("cluster").size().to_dict()
                cluster_stats = {"sizes": {str(k): int(v) for k, v in sizes.items()}, "num_points": int(len(dfc))}

    # taxonomy summary per sample: top taxa
    taxonomy_summary = []
    if not taxonomy_df.empty:
        for sample, grp in taxonomy_df.groupby("sample"):
            top = grp.sort_values("abundance", ascending=False).head(10)
            taxonomy_summary.append({"sample": sample, "top_taxa": top.to_dict(orient="records")})

    # sankey + sunburst data for taxonomy (frontend)
    taxonomy_sankey = build_sankey_from_taxonomy(taxonomy_df)
    taxonomy_sunburst = build_sunburst_data(taxonomy_df)

    # status.json and artifacts
    status_data = parse_status_json(run_dir)
    artifacts = []
    if status_data and isinstance(status_data.get("outputs"), dict):
        for key, val in status_data["outputs"].items():
            url = make_artifact_url(run_dir, val)
            artifacts.append({"label": key, "url": url, "source_path": val, "type": Path(str(val)).suffix.lstrip(".")})
    else:
        # fallback: common files
        fallback_files = [
            ("fastp_html", "qc/fastp_report.html"),
            ("report_html", "reports/report.html"),
            ("phylogeny_tree", "phylogeny/tree.nwk"),
            ("clusters", "clustering/clusters.csv"),
            ("novelty", "novelty/novelty_report.tsv"),
        ]
        for label, rel in fallback_files:
            candidate = run_dir / rel
            if candidate.exists():
                artifacts.append({"label": label, "url": make_artifact_url(run_dir, rel), "source_path": str(candidate), "type": candidate.suffix.lstrip(".")})

    # ensure taxonomy has sensible fallback
    if taxonomy_df.empty:
        taxonomy_df = pd.DataFrame([{"sample": "sample1", "taxon": "Unclassified", "abundance": 1.0}])

    # prepare novelty records as JSON-serializable (limit to first N for immediate frontend)
    novelty_records = []
    if not novelty_df.empty:
        novelty_records = novelty_df.fillna("").to_dict(orient="records")
        # keep just a reasonable page (frontend can request raw TSV artifact for full table)
        novelty_records_preview = novelty_records[:200]
    else:
        novelty_records_preview = []

    summary = {
        "run_id": run_dir.name,
        "status": "completed",
        "created_at": datetime.utcnow().isoformat() + "Z",
        "summaryMetrics": {
            "alpha_diversity": alpha_div,
            "beta_diversity": beta_div,
            "novelty_stats": novelty_stats,
        },
        "qc": qc,
        "novelty": {"records_preview": novelty_records_preview, "artifact": next((a for a in artifacts if "novelty" in a["label"].lower()), None), "summary": novelty_stats},
        "taxonomy": taxonomy_df_limited.to_dict(orient="records"),
        "taxonomy_summary": taxonomy_summary,
        "taxonomy_sankey": taxonomy_sankey,
        "taxonomy_sunburst": taxonomy_sunburst,
        "clustering": {"points": clustering_list, "stats": cluster_stats, "artifact": next((a for a in artifacts if "cluster" in a["label"].lower()), None)},
        "phylogeny": {"tree_artifact": next((a for a in artifacts if "phylogeny" in a["label"].lower() or "tree" in a["label"].lower()), None)},
        "artifacts": artifacts,
        "status_json": status_data or {},
        "errors": status_data.get("errors", []) if status_data else [],
    }
    return summary

def save_summary(run_dir: Path) -> Path:
    summary = aggregate_results(run_dir)
    out_file = run_dir / "summary.json"
    out_file.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"✅ Saved {out_file}")
    return out_file

if __name__ == "__main__":
    import sys
    arg = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    save_summary(arg)
