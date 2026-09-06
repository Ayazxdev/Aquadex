from pathlib import Path
from typing import Iterable
from ..core.config import settings

def runs_dir() -> Path:
    return Path(settings.DATA_IN)

def out_dir(run_id: str) -> Path:
    return Path(settings.DATA_OUT) / run_id

def in_dir(run_id: str) -> Path:
    return Path(settings.DATA_IN) / run_id

def ensure_run_dirs(run_id: str) -> tuple[Path, Path]:
    inp = in_dir(run_id)
    outp = out_dir(run_id)
    inp.mkdir(parents=True, exist_ok=True)
    outp.mkdir(parents=True, exist_ok=True)
    return inp, outp
