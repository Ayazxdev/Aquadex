from fastapi import APIRouter, UploadFile, File, Form, HTTPException
from fastapi.responses import FileResponse
from typing import List, Optional
from ..models.schemas import UploadResponse, RunRequest, RunResponse, StatusResponse, ResultsResponse
from ..services.storage import ensure_run_dirs, out_dir
from ..services.pipeline import start_pipeline, get_status
from ..services import results as results_svc

router = APIRouter(prefix="/api")


@router.post("/upload", response_model=UploadResponse)
async def upload_files(files: List[UploadFile] = File(...), run_id: Optional[str] = Form(None)):
    if not run_id:
        from ..utils.run_id import new_run_id
        run_id = new_run_id()
    inp, _ = ensure_run_dirs(run_id)
    saved = []
    for f in files:
        dest = inp / f.filename
        with dest.open("wb") as w:
            w.write(await f.read())
        saved.append(f.filename)
    return UploadResponse(run_id=run_id, saved_files=saved)


@router.post("/run", response_model=RunResponse)
async def run_pipeline(req: RunRequest):
    ok = start_pipeline(req.run_id, req.marker, req.read_type, req.options)
    if not ok:
        raise HTTPException(status_code=400, detail="Failed to start pipeline")
    return RunResponse(run_id=req.run_id, started=True, message="Pipeline launched")


@router.get("/status/{run_id}", response_model=StatusResponse)
async def status(run_id: str):
    s = get_status(run_id)
    return StatusResponse(
        run_id=run_id,
        status=s.get("status", "unknown"),
        progress=float(s.get("progress", 0.0)),
        message=s.get("message"),
    )


@router.get("/results/{run_id}")
async def results(run_id: str):
    import json
    from ..services.storage import out_dir
    from ..services import results as results_svc
    from fastapi.responses import JSONResponse
    
    headers = {
        "Cache-Control": "no-cache, no-store, must-revalidate",
        "Pragma": "no-cache",
        "Expires": "0"
    }
    
    summary_file = out_dir(run_id) / "summary.json"
    if summary_file.exists():
        try:
            data = json.loads(summary_file.read_text(encoding="utf-8"))
            return JSONResponse(content=data, headers=headers)
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Error reading summary.json: {e}")
            
    # Fallback if summary.json doesn't exist yet
    data = results_svc.compose_dashboard(run_id)
    return JSONResponse(content=data, headers=headers)


@router.get("/taxonomy/{run_id}")
async def taxonomy(run_id: str):
    return {"run_id": run_id, "taxonomyTable": results_svc.read_taxonomy(run_id)}


@router.get("/clustering/{run_id}")
async def clustering(run_id: str):
    rows = results_svc.read_clustering(run_id)
    return {
        "run_id": run_id,
        "clusteringTable": rows,
        "umap_coordinates": results_svc.clustering_to_umap(rows),
    }


@router.get("/artifacts/{run_id}/{path:path}")
async def artifact(run_id: str, path: str):
    base = out_dir(run_id)
    file_path = (base / path).resolve()
    if not str(file_path).startswith(str(base.resolve())):
        raise HTTPException(status_code=400, detail="Invalid path")
    if not file_path.exists():
        raise HTTPException(status_code=404, detail="Not found")
    return FileResponse(str(file_path), filename=file_path.name)


@router.get("/download_all/{run_id}")
async def download_all(run_id: str):
    import zipfile
    from fastapi.responses import FileResponse
    
    base = out_dir(run_id)
    if not base.exists():
        raise HTTPException(status_code=404, detail="Run results not found")
        
    zip_path = base / f"aquadex_results_{run_id}.zip"
    if not zip_path.exists():
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED, compresslevel=1) as zip_file:
            for file_path in base.rglob("*"):
                if file_path.is_file() and file_path.name != zip_path.name:
                    arcname = file_path.relative_to(base)
                    zip_file.write(file_path, arcname=str(arcname))
                    
    return FileResponse(
        str(zip_path),
        media_type="application/zip",
        filename=f"aquadex_results_{run_id}.zip"
    )


@router.get("/asv_sequence/{run_id}/{asv_id:path}")
async def asv_sequence(run_id: str, asv_id: str):
    """
    Returns the actual ACGT nucleotide sequence for a given ASV ID by looking it up
    in the dereplicated FASTA file produced by the denoise step.
    ASV IDs are like: SRR36837148.44_size=500 or SRR36837148.44;size=500
    """
    from pathlib import Path
    from fastapi.responses import JSONResponse

    base = out_dir(run_id)
    fasta_candidates = [
        base / "denoise" / "sample1_R1.derep.fasta",
        base / "denoise" / "sample1_R1.cleaned.fasta",
    ] + list((base / "denoise").glob("*.fasta")) if (base / "denoise").exists() else []

    # Normalize the query ASV ID: convert ; to _ for matching
    query_norm = asv_id.replace(";", "_")

    seq_found = None
    length = 0
    gc_content = 0.0

    for fasta_path in fasta_candidates:
        if not fasta_path.exists():
            continue
        curr_header = None
        curr_seq_parts = []
        with open(fasta_path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    if curr_header is not None:
                        # Check if previous header matches
                        header_norm = curr_header.replace(";", "_")
                        if header_norm == query_norm or header_norm.startswith(query_norm):
                            seq_found = "".join(curr_seq_parts)
                            break
                    curr_header = line[1:].replace(";size=", "_size=")
                    curr_seq_parts = []
                else:
                    curr_seq_parts.append(line)
            # Check last sequence
            if seq_found is None and curr_header is not None:
                header_norm = curr_header.replace(";", "_")
                if header_norm == query_norm or header_norm.startswith(query_norm):
                    seq_found = "".join(curr_seq_parts)

        if seq_found:
            break

    if not seq_found:
        raise HTTPException(status_code=404, detail=f"Sequence not found for ASV: {asv_id}")

    length = len(seq_found)
    g_count = seq_found.upper().count("G")
    c_count = seq_found.upper().count("C")
    gc_content = round((g_count + c_count) / length * 100, 1) if length > 0 else 0.0

    return JSONResponse(content={
        "asv_id": asv_id,
        "sequence": seq_found,
        "length": length,
        "gc_content": gc_content,
    })
