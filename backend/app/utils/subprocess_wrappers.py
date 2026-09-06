import subprocess
import logging

logger = logging.getLogger("edna_pipeline")

def run_cmd(cmd, cwd=None, shell=False):
    """
    Run a command in subprocess and stream output.
    Returns (stdout, stderr, returncode).
    """
    logger.info(f"Running command: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    process = subprocess.Popen(
        cmd,
        cwd=cwd,
        shell=shell,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        logger.error(f"Command failed [{process.returncode}]: {stderr}")
        raise RuntimeError(f"Command failed: {stderr}")
    return stdout, stderr, process.returncode