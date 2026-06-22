# EASMS data-curation pipeline (Step 0 QC + Steps 1-9).
# Mirrors the sibling ../Dockerfile: slim Python plus the X libs the RDKit
# wheels need at runtime.
FROM python:3.12-slim

# RDKit wheels link against a couple of X libraries even for headless use.
RUN apt-get update && apt-get install -y --no-install-recommends \
      libxrender1 libxext6 && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install deps first so the layer is cached across code changes.
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# App code + reference files (ASMS Meta Data.csv, MasterLists_sample/, ...).
# Providers.csv is gitignored and NOT baked in -- mount it at runtime (see below).
COPY . .

# The src/ modules use flat imports (`from io_utils import ...`) and rely on the
# script's own directory being on sys.path[0], so Main.py must be invoked
# directly (never `python -m`). Pipeline flags are appended after the image name.
ENTRYPOINT ["python", "src/Main.py"]

# --- Usage ----------------------------------------------------------------
# PRODUCTION (Cloud Run Job) -- keyless. A dedicated service account is attached
# to the job, so ADC comes from the metadata server: NO credentials are baked
# into the image or mounted. Config (Providers.csv) and data live in GCS and are
# passed as gs:// args, so nothing secret is in the image. Use the Makefile for
# the full build + service-account IAM + deploy flow, step by step:
#   make repo build sa iam deploy PROJECT_ID=<PROJECT_ID>   (or: make all ...)
#   make run PROJECT_ID=<PROJECT_ID>                        (execute the job)
#
# Cap fingerprint parallelism to the job's allocated vCPUs (os.cpu_count can
# over-report under cgroup limits): set FP_N_JOBS, e.g. -e FP_N_JOBS=8.
#
# LOCAL DEV ONLY -- mount your *user* ADC (gcloud auth application-default
# login) read-only; never use personal creds in production:
#   docker build -t easms-curation .
#   docker run --rm \
#     -v "$HOME/.config/gcloud:/root/.config/gcloud:ro" \
#     -e GOOGLE_CLOUD_PROJECT=<PROJECT_ID> \
#     easms-curation \
#     --input-file       gs://<bucket>/asms/run.csv \
#     --output-dir       gs://<bucket>/out \
#     --masterlists-dir  gs://<bucket>/library/EASMS12kV1_lib.csv \
#     --providers-csv    gs://<bucket>/config/Providers.csv
#
# Fully local (no GCS): mount data + Providers.csv (QC needs it) to /app:
#   docker run --rm -v "$PWD/data:/data" \
#     -v "$PWD/Providers.csv:/app/Providers.csv:ro" \
#     easms-curation --input-file /data/run.csv --output-dir /data/out
