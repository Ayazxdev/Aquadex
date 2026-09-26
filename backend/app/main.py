from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from .api.routes import router as api_router

app = FastAPI(title="eDNA Backend", version="1.0.0")

# CORS: allow all origins and Vercel domains
app.add_middleware(
    CORSMiddleware,
    allow_origin_regex=r"https?://.*",
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Support both /api/* and root /* so any URL configuration works seamlessly
app.include_router(api_router, prefix="/api")
app.include_router(api_router)

@app.get("/")
def root():
    return {"ok": True, "service": "AQUADEX eDNA Backend", "health": "/health", "api": "/api"}

@app.get("/health")
def health():
    return {"ok": True}
