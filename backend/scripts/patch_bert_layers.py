"""
Runs at container startup to patch bert_layers.py files.
Replaces the broken Triton flash-attn import block with a simple None assignment
so DNABERT-S falls back to standard PyTorch CUDA attention.
"""
import glob
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def patch_bert_layers(fpath: str) -> None:
    with open(fpath, "r", encoding="utf-8") as f:
        lines = f.readlines()

    new_lines = []
    skip = False
    already_patched = any("flash_attn_qkvpacked_func = None" in l and "try:" not in l for l in lines[:40])
    if already_patched:
        logger.info("Already patched: %s", fpath)
        return

    for line in lines:
        stripped = line.strip()
        if stripped == "try:":
            skip = True
            new_lines.append("flash_attn_qkvpacked_func = None  # Triton disabled at startup\n")
            continue
        if skip:
            # Skip until we hit the logger line which comes after the try/except block
            if stripped.startswith("logger = logging.getLogger") or stripped.startswith("class "):
                skip = False
                new_lines.append(line)
            continue
        new_lines.append(line)

    with open(fpath, "w", encoding="utf-8") as f:
        f.writelines(new_lines)
    logger.info("Patched: %s", fpath)


def main():
    targets = (
        glob.glob("/app/models/**/bert_layers.py", recursive=True)
        + glob.glob("/root/.cache/huggingface/modules/**/bert_layers.py", recursive=True)
    )
    if not targets:
        logger.info("No bert_layers.py files found yet — will patch on first model load.")
        return
    for t in set(targets):
        try:
            patch_bert_layers(t)
        except Exception as exc:
            logger.warning("Could not patch %s: %s", t, exc)


if __name__ == "__main__":
    main()
