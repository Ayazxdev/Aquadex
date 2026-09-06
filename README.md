# AQUADEX

Environmental DNA (eDNA) biodiversity analysis platform with taxonomic classification and novelty scoring.

## Architecture

- **backend**: FastAPI service wrapping Kraken2, FAISS, and downstream reporting.
- **frontend**: React (Vite) interface for file upload, pipeline tracking, and result exploration.

## Quick Start

### Docker (Recommended)

```bash
docker compose up --build
```
Backend runs at `http://localhost:8001`.

### Local Development

#### Backend
```bash
cd backend
python -m venv .venv
# Windows: .venv\Scripts\activate | Unix: source .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env
uvicorn app.main:app --reload --port 8001
```

#### Frontend
```bash
cd frontend
npm install
npm run dev
```
Frontend runs at `http://localhost:5173` and proxies `/api` requests to port 8001.

## Key Endpoints

- `POST /api/upload`: Upload input sequencing files
- `POST /api/run`: Trigger analysis pipeline
- `GET /api/status/{run_id}`: Poll run execution status
- `GET /api/results/{run_id}`: Fetch aggregated run metrics
- `GET /api/artifacts/{run_id}/{path}`: Download run outputs and reports
