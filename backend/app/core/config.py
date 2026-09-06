from pydantic_settings import BaseSettings
from pydantic import Field


class Settings(BaseSettings):
    DATA_IN: str = Field("/data/in", validation_alias="DATA_IN")
    DATA_OUT: str = Field("/data/out", validation_alias="DATA_OUT")
    DATA_DB: str = Field("/data/db", validation_alias="DATA_DB")
    DATA_TMP: str = Field("/data/tmp", validation_alias="DATA_TMP")

    PIPELINE_CMD: str | None = None
    ARTIFACT_BASE_URL: str = "/api/artifacts"

    # Demo mode: synthetic outputs without bioinformatics binaries (default on Windows)
    DEMO_MODE: bool = True

    # Real pipeline paths (used when DEMO_MODE=false, typically inside Docker)
    KRAKEN2_DB: str = "/data/db/kraken2"
    DNABERT_MODEL: str = "/app/models/dnabert-s"

    class Config:
        env_file = ".env"


settings = Settings()
