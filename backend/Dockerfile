# AQUADEX backend — full bioinformatics stack for real pipeline mode
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV DEMO_MODE=false

RUN apt-get update && apt-get install -y \
    git wget curl unzip build-essential cmake \
    python3 python3-pip python3-venv \
    openjdk-11-jre-headless \
    r-base r-base-dev \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev \
    fastp fastqc multiqc seqkit seqtk cutadapt vsearch \
    bbmap ncbi-blast+ kraken2 \
    && rm -rf /var/lib/apt/lists/*

RUN wget http://github.com/bbuchfink/diamond/releases/download/v2.1.13/diamond-linux64.tar.gz \
    && tar xzf diamond-linux64.tar.gz \
    && mv diamond /usr/local/bin/ \
    && rm diamond-linux64.tar.gz

RUN pip3 install --upgrade pip

WORKDIR /app
COPY requirements.txt /app/
RUN pip3 install --no-cache-dir --ignore-installed sympy -r requirements.txt

# Extra runtime deps that the pipeline needs
RUN pip3 install --no-cache-dir h5py>=3.10.0 einops>=0.6.0 "Pillow>=10.0.0" triton

RUN mkdir -p /data/in /data/out /data/db /data/tmp /app/models

COPY app /app/app

# Patch bert_layers.py inside the model folder: disable broken Triton flash-attn,
# fall back to the PyTorch CUDA attention path which works correctly.
COPY scripts/patch_bert_layers.py /app/scripts/patch_bert_layers.py

ENV PORT=10000
EXPOSE 10000
CMD ["/bin/sh", "-c", "python3 /app/scripts/patch_bert_layers.py && uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-10000}"]

