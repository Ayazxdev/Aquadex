from pathlib import Path
import pandas as pd
import logging
from typing import Tuple, Optional

class ReportsRunner:
    """
    Build a combined TSV summary and a simple HTML report.
    """

    def __init__(self, workdir: Path):
        self.workdir = Path(workdir)
        self.out_dir = self.workdir / "reports"
        self.out_dir.mkdir(parents=True, exist_ok=True)

        self.qc_html = self.workdir / "qc" / "fastp_report.html"
        self.qc_json = self.workdir / "qc" / "fastp_report.json"

    def run(
        self,
        taxonomy_file: Path,
        novelty_file: Path,
        clustering_file: Path,
        tree_file: Path,
    ) -> Tuple[Path, Path]:
        """
        Main method to run report generation. Accepts all required file paths.
        Always writes a TSV + HTML, even if inputs are missing or empty.
        """
        logging.info("Generating reports")

        # Load inputs
        anchor_df = self._safe_read_csv(self.workdir / "anchor" / "anchor_matches.tsv", sep="\t")
        taxonomy_df = self._safe_read_csv(taxonomy_file, sep="\t")
        novelty_df = self._safe_read_csv(novelty_file, sep="\t")
        cluster_df = self._safe_read_csv(clustering_file, sep="\t")
        phylo_df = self._safe_read_csv(self.workdir / "phylogeny" / "phylogeny_table.tsv", sep="\t")

        # Build summary
        summary_df = pd.DataFrame()

        # Example: if novelty_df is present, add it to summary
        if novelty_df is not None and not novelty_df.empty:
            summary_df = novelty_df.copy()

        # Always save a TSV
        self.tsv_out = self.out_dir / "reports_summary.tsv"
        if not summary_df.empty:
            summary_df.to_csv(self.tsv_out, sep="\t", index=False)
        else:
            # Write a placeholder TSV with just a message
            pd.DataFrame({"message": ["No results available"]}).to_csv(
                self.tsv_out, sep="\t", index=False
            )

        # Build HTML report
        html_parts = [
            "<html><head><title>Pipeline Report</title></head><body>",
            "<h1>Pipeline Report</h1>",
        ]

        if taxonomy_df is None or taxonomy_df.empty:
            html_parts.append("<p><b>Taxonomy:</b> No results available.</p>")
        else:
            html_parts.append("<h2>Taxonomy</h2>")
            html_parts.append(taxonomy_df.head(20).to_html(index=False))

        if novelty_df is None or novelty_df.empty:
            html_parts.append("<p><b>Novelty:</b> No results available.</p>")
        else:
            html_parts.append("<h2>Novelty</h2>")
            html_parts.append(novelty_df.head(20).to_html(index=False))

        if cluster_df is None or cluster_df.empty:
            html_parts.append("<p><b>Clustering:</b> No results available.</p>")
        else:
            html_parts.append("<h2>Clustering</h2>")
            html_parts.append(cluster_df.head(20).to_html(index=False))

        if phylo_df is None or phylo_df.empty:
            html_parts.append("<p><b>Phylogeny:</b> No results available.</p>")
        else:
            html_parts.append("<h2>Phylogeny</h2>")
            html_parts.append(phylo_df.head(20).to_html(index=False))

        html_parts.append("</body></html>")

        self.html_out = self.out_dir / "report.html"
        with open(self.html_out, "w", encoding="utf-8") as fh:
            fh.write("\n".join(html_parts))

        return (self.tsv_out, self.html_out)

    def _safe_read_csv(self, p: Path, **kwargs) -> Optional[pd.DataFrame]:
        if not p.exists():
            return None
        try:
            df = pd.read_csv(p, **kwargs)
            if df.empty:
                return None
            return df
        except Exception as e:
            logging.warning(f"Could not read {p}: {e}")
            return None
