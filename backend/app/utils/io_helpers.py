import json
import yaml
import hashlib
from pathlib import Path

def write_json(obj, path, indent: int = 2):
    """Write Python object as JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=indent)

def write_yaml(obj, path):
    """Write Python object as YAML."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(obj, f, sort_keys=False)

def read_json(path):
    """Read JSON file into Python object."""
    path = Path(path)
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def read_yaml(path):
    """Read YAML file into Python object."""
    path = Path(path)
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def compute_sha256(path, block_size: int = 65536) -> str:
    """Compute SHA256 checksum of a file."""
    sha = hashlib.sha256()
    path = Path(path)
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(block_size), b""):
            sha.update(block)
    return sha.hexdigest()