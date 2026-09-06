import os
from pathlib import Path

def get_base_workdir():
    return os.environ.get("WORKDIR", "/app/workdir")

def ensure_run_dir(run_id):
    if run_id is None:
        raise ValueError("run_id must be provided")
    base = get_base_workdir()
    run_dir = os.path.join(base, str(run_id))
    Path(run_dir).mkdir(parents=True, exist_ok=True)
    return run_dir

def ensure_step_dir(run_id, step_name):
    run_dir = ensure_run_dir(run_id)
    step_dir = os.path.join(run_dir, step_name)
    Path(step_dir).mkdir(parents=True, exist_ok=True)
    return step_dir
