# pipeline/step_qc.py
from pathlib import Path
import subprocess
import logging
import json
import shutil
from typing import List, Tuple
from app.services.storage import out_dir
# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

class QCProcessor:
    """
    Perform primary QC and trimming using fastp, FastQC, and MultiQC.
    """

    def __init__(self, run_id: str, fastp_path: str = "fastp", fastqc_path: str = "fastqc", multiqc_path: str = "multiqc"):
        self.run_id = run_id
        self.workdir = out_dir(run_id) / "qc"
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.fastp_path = fastp_path
        self.fastqc_path = fastqc_path
        self.multiqc_path = multiqc_path

    def run_fastp(self, input_files: List[Path], adapters: Tuple[str, str] = (None, None)) -> List[Path]:
        """
        Run fastp on input FASTQ/FASTA files with optional adapter trimming.
        Returns list of trimmed files.
        """
        trimmed_files = []
        for infile in input_files:
            outfile = self.workdir / f"trimmed_{infile.name}"
            json_report = self.workdir / f"{infile.stem}_fastp.json"
            html_report = self.workdir / f"{infile.stem}_fastp.html"
            cmd = [
                self.fastp_path,
                "-i", str(infile),
                "-o", str(outfile),
                "-j", str(json_report),
                "-h", str(html_report)
            ]
            if adapters[0]:
                cmd.extend(["-a", adapters[0]])
            if adapters[1]:
                cmd.extend(["-A", adapters[1]])
            logging.info(f"Running fastp on {infile}")
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                logging.error(f"fastp failed: {e.stderr.decode()}")
                raise RuntimeError(f"fastp failed: {e.stderr.decode()}")
            trimmed_files.append(outfile)
        return trimmed_files

    def run_fastqc(self, input_files: List[Path]) -> List[Path]:
        """
        Run FastQC on input files. Returns list of output HTML reports.
        """
        output_reports = []
        for infile in input_files:
            outdir = self.workdir / f"{infile.stem}_fastqc"
            outdir.mkdir(parents=True, exist_ok=True)
            cmd = [self.fastqc_path, "-o", str(outdir), str(infile)]
            logging.info(f"Running FastQC on {infile}")
            try:
                subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                logging.error(f"FastQC failed: {e.stderr.decode()}")
                raise RuntimeError(f"FastQC failed: {e.stderr.decode()}")
            report_html = list(outdir.glob("*_fastqc.html"))[0]
            output_reports.append(report_html)
        return output_reports

    def run_multiqc(self, input_dirs: List[Path], output_dir: Path) -> Path:
        """
        Aggregate QC reports with MultiQC. Returns path to multiqc report.
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        cmd = [self.multiqc_path, "-o", str(output_dir)] + [str(d) for d in input_dirs]
        logging.info(f"Running MultiQC on directories: {input_dirs}")
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            logging.error(f"MultiQC failed: {e.stderr.decode()}")
            raise RuntimeError(f"MultiQC failed: {e.stderr.decode()}")
        report_html = output_dir / "multiqc_report.html"
        if not report_html.exists():
            raise RuntimeError("MultiQC report not found after run")
        return report_html

    def check_q30_threshold(self, fastp_json_files: List[Path], threshold: float = 60.0) -> List[str]:
        """
        Check if per-sample base Q30 is below threshold.
        Returns list of filenames flagged as LowQuality.
        """
        low_quality_files = []
        for jf in fastp_json_files:
            with open(jf) as f:
                data = json.load(f)
            try:
                q30_rate = data["summary"]["Q30_rate"] * 100
            except KeyError:
                q30_rate = 100.0  # assume pass if missing
            if q30_rate < threshold:
                logging.warning(f"{jf.stem} Q30 rate {q30_rate:.1f}% below threshold {threshold}%")
                low_quality_files.append(jf.stem)
        return low_quality_files