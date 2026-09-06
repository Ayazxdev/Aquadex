from pydantic import BaseModel, Field
from typing import List, Optional, Literal, Dict, Any


class UploadResponse(BaseModel):
    run_id: str
    saved_files: List[str]


class RunRequest(BaseModel):
    run_id: str
    marker: Literal["16S", "18S", "COI"]
    read_type: Literal["short", "long"] = "short"
    options: Dict[str, Any] = {}


class RunResponse(BaseModel):
    run_id: str
    started: bool
    message: str


class StatusResponse(BaseModel):
    run_id: str
    status: Literal["queued", "running", "completed", "failed", "unknown"]
    progress: float = 0.0
    message: str | None = None


class SummaryMetric(BaseModel):
    label: str
    value: str


class NoveltyRow(BaseModel):
    id: str
    noveltyScore: str = ""
    vaeloss: str = ""
    faissDist: str = ""
    epaAnnotation: str = ""
    diamondHit: str = ""
    abundance: str = ""
    depth: str = ""
    location: str = ""
    lat: str = ""
    lon: str = ""


class TaxonomyRow(BaseModel):
    ASV_ID: str = ""
    taxon: str = ""
    abundance: str = ""
    sample: str = ""


class ClusteringRow(BaseModel):
    ASV_ID: str = ""
    cluster_id: str = ""
    dim_1: str = ""
    dim_2: str = ""
    novelty_score: str = ""
    cluster_size: str = ""


class Artifact(BaseModel):
    filename: str
    url: str


class ResultsResponse(BaseModel):
    run_id: str
    status: str = "unknown"
    summaryMetrics: List[SummaryMetric] = []
    noveltyTable: List[NoveltyRow] = []
    taxonomyTable: List[Dict[str, Any]] = []
    clusteringTable: List[Dict[str, Any]] = []
    umap_coordinates: Dict[str, Any] = {}
    qc: Dict[str, Any] = {}
    taxonomySankey: Dict[str, Any] = {}
    taxonomyStackedBar: List[Dict[str, Any]] = []
    artifacts: List[Artifact] = []
