# pipeline/step_packaging.py
from pathlib import Path
import json
import logging
import shutil
import datetime
import uuid
import subprocess
from typing import Dict, List, Optional, Any
from utils.io_helpers import compute_sha256

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def _run_cmd_capture(cmd: List[str]) -> str:
    """
    Run a command and return stdout. Raise RuntimeError on failure.
    """
    logging.info(f"Running command for version detection: {' '.join(cmd)}")
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True, text=True)
        return res.stdout.strip().splitlines()[0] if res.stdout else ""
    except subprocess.CalledProcessError as e:
        logging.warning(f"Version command failed: {' '.join(cmd)} -> {e}")
        return "unknown"


class Packager:
    """
    Create run manifest, compute checksums, and build submission bundle zip.
    """

    def __init__(self, base_workdir: Path):
        """
        base_workdir: directory where packaging will place run artifacts
        """
        self.base_workdir = base_workdir
        self.base_workdir.mkdir(parents=True, exist_ok=True)

    def _collect_checksums(self, files: Dict[str, Path]) -> Dict[str, str]:
        checks = {}
        for k, p in files.items():
            if p and Path(p).exists():
                checks[k] = compute_sha256(Path(p))
            else:
                checks[k] = None
        return checks

    def _detect_software_versions(self, tools: List[str]) -> Dict[str, str]:
        versions = {}
        for t in tools:
            try:
                out = _run_cmd_capture([t, "--version"])
                versions[t] = out or "unknown"
            except Exception:
                versions[t] = "unknown"
        return versions

    def create_manifest(
        self,
        outputs: Dict[str, Path],
        model_id: Optional[str],
        db_versions: Optional[Dict[str, str]],
        random_seed: Optional[int],
        tools_to_probe: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        Build a run manifest dictionary containing metadata and checksums.
        """
        run_id = f"run_{datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%SZ')}_{uuid.uuid4().hex[:8]}"
        start_time = datetime.datetime.utcnow().isoformat()
        manifest_dir = self.base_workdir / run_id
        manifest_dir.mkdir(parents=True, exist_ok=True)

        # copy outputs into run directory (do not modify originals)
        saved_outputs = {}
        for key, path in outputs.items():
            if path is None:
                saved_outputs[key] = None
                continue
            src = Path(path)
            if not src.exists():
                saved_outputs[key] = None
                continue
            dst = manifest_dir / src.name
            shutil.copy2(src, dst)
            saved_outputs[key] = dst

        # detect software versions
        tools = tools_to_probe or ["fastp", "cutadapt", "vsearch", "Rscript", "diamond", "epa-ng"]
        software_versions = self._detect_software_versions(tools)

        checksums = self._collect_checksums(saved_outputs)

        manifest = {
            "run_id": run_id,
            "start_time": start_time,
            "end_time": datetime.datetime.utcnow().isoformat(),
            "software_versions": software_versions,
            "model_id": model_id or "not_provided",
            "reference_db_versions": db_versions or {},
            "random_seed": random_seed,
            "outputs": {k: (str(v) if v is not None else None) for k, v in saved_outputs.items()},
            "checksums": checksums,
        }

        manifest_path = manifest_dir / "run_manifest.json"
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2)

        return {"manifest": manifest, "manifest_path": manifest_path, "manifest_dir": manifest_dir}

    def create_submission_zip(self, manifest_dir: Path, zip_name: Optional[str] = None) -> Path:
        """
        Zip the manifest directory into submission bundle.
        """
        if zip_name:
            zip_path = manifest_dir.parent / zip_name
        else:
            zip_path = manifest_dir.parent / f"{manifest_dir.name}_submission_bundle.zip"
        logging.info(f"Creating submission bundle: {zip_path}")
        shutil.make_archive(str(zip_path.with_suffix('')), 'zip', root_dir=str(manifest_dir))
        return zip_path

    def package_run(
        self,
        outputs: Dict[str, Path],
        model_id: Optional[str] = None,
        db_versions: Optional[Dict[str, str]] = None,
        random_seed: Optional[int] = None,
        tools_to_probe: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        """
        High-level: create manifest, zip bundle, and return paths and manifest content.
        """
        created = self.create_manifest(outputs, model_id, db_versions, random_seed, tools_to_probe)
        manifest = created["manifest"]
        manifest_dir = created["manifest_dir"]
        manifest_path = created["manifest_path"]
        zip_path = self.create_submission_zip(manifest_dir)
        # compute final zip checksum
        zip_checksum = compute_sha256(zip_path)
        manifest["submission_bundle"] = str(zip_path)
        manifest["submission_bundle_checksum"] = zip_checksum
        # overwrite manifest with updated info
        with open(manifest_path, "w") as fh:
            json.dump(manifest, fh, indent=2)
        return {
            "manifest_path": manifest_path,
            "submission_bundle": zip_path,
            "manifest": manifest
        }