#!/usr/bin/env bash

set -euo pipefail

# ---- Config (override via env) -------------------------------------------
PROJECT_ID="${PROJECT_ID:?set PROJECT_ID}"
REGION="${REGION:-us-central1}"
AR_REPO="${AR_REPO:-easms}"
JOB_NAME="${JOB_NAME:-easms-curation}"
TAG="${TAG:-latest}"
IMAGE="${IMAGE:-${REGION}-docker.pkg.dev/${PROJECT_ID}/${AR_REPO}/${JOB_NAME}:${TAG}}"


platform_args=()
[ -n "${PLATFORM:-}" ] && platform_args=(--platform "$PLATFORM")

INPUT_FILE="${INPUT_FILE:-gs://asms_sgc-toronto/data/asms_sgcto_20_EASMS12kV1_20260610.csv}"
OUTPUT_DIR="${OUTPUT_DIR:-gs://test-aircheck/asms-test-run-docker}"
MASTERLISTS="${MASTERLISTS:-gs://asms_sgc-toronto/library/EASMS12kV1_lib.csv}"

PROVIDERS_CSV="${PROVIDERS_CSV:-$PWD/Providers.csv}"



# ---- 1. Authenticate Docker to Artifact Registry via gcloud --------------
gcloud auth configure-docker "${REGION}-docker.pkg.dev" --quiet

# ---- 2. Pull the image ---------------------------------------------------
docker pull "${platform_args[@]}" "$IMAGE"

# ---- 3. Run locally ------------------------------------------------------

run_args=(
  --rm
  "${platform_args[@]}"
  -v "$HOME/.config/gcloud:/root/.config/gcloud:ro"
  -e GOOGLE_CLOUD_PROJECT="$PROJECT_ID"

  -e GRPC_VERBOSITY="${GRPC_VERBOSITY:-ERROR}"
)
# Worker count for fingerprinting. Unset -> the container uses os.cpu_count()
# (all cores it can see). Set FP_N_JOBS to cap it.
[ -n "${FP_N_JOBS:-}" ] && run_args+=(-e FP_N_JOBS="$FP_N_JOBS")
[ -n "${DATA_DIR:-}" ]  && run_args+=(-v "${DATA_DIR}:/data")


providers_args=()
if [ -n "${PROVIDERS_CSV:-}" ]; then
  case "$PROVIDERS_CSV" in
    *://*)
      providers_args=(--providers-csv "$PROVIDERS_CSV") ;;
    *)
      [ -f "$PROVIDERS_CSV" ] || { echo "ERROR: PROVIDERS_CSV='$PROVIDERS_CSV' not found on host." >&2; exit 1; }
      providers_abs="$(cd "$(dirname "$PROVIDERS_CSV")" && pwd)/$(basename "$PROVIDERS_CSV")"
      run_args+=(-v "${providers_abs}:/app/Providers.csv:ro")
      providers_args=(--providers-csv /app/Providers.csv) ;;
  esac
fi

# Assemble pipeline args from the variables above; any extra flags passed on the
# command line ("$@") are appended (e.g. --start-from 1 --end-at 7). A flag is
# omitted when its variable is empty.
pipeline_args=()
[ -n "$INPUT_FILE" ]    && pipeline_args+=(--input-file "$INPUT_FILE")
[ -n "$OUTPUT_DIR" ]    && pipeline_args+=(--output-dir "$OUTPUT_DIR")
[ -n "$MASTERLISTS" ]   && pipeline_args+=(--masterlists-dir "$MASTERLISTS")
[ "${#providers_args[@]}" -gt 0 ] && pipeline_args+=("${providers_args[@]}")
pipeline_args+=("$@")

# Nothing to run -> show the pipeline help.
if [ "${#pipeline_args[@]}" -eq 0 ]; then
  pipeline_args=(--help)
fi

exec docker run "${run_args[@]}" "$IMAGE" "${pipeline_args[@]}"
