# Backend Service

FastAPI service for pipeline orchestration, classification, and metrics computation.

## Setup

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate | Unix: source .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env
uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
```

## Endpoints

- `POST /api/upload`: Receive raw FASTA/FASTQ sequence files
- `POST /api/run`: Initialize pipeline execution
- `GET /api/status/{run_id}`: Poll run status and progress
- `GET /api/results/{run_id}`: Retrieve summarized run metrics
- `GET /api/artifacts/{run_id}/{path}`: Retrieve generated output files
