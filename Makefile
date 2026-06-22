
-include .env

# ---- Config (override as VAR=value) --------------------------------------
PROJECT_ID      ?=
REGION          ?= us-central1
AR_REPO         ?= easms
JOB_NAME        ?= easms-curation
RUNTIME_SA_NAME ?= easms-curation-job
INPUT_BUCKET    ?= asms_sgc-toronto
OUTPUT_BUCKET   ?= test-aircheck
CPU             ?= 8
MEMORY          ?= 16Gi
TASK_TIMEOUT    ?= 3600
MAX_RETRIES     ?= 1
PLATFORM        ?= linux/amd64

INPUT_FILE  ?= gs://$(INPUT_BUCKET)/data/asms_sgcto_1_DEMOLIB_20260101.csv
OUTPUT_DIR  ?= gs://$(OUTPUT_BUCKET)/asms-out
MASTERLISTS ?= gs://$(INPUT_BUCKET)/library/EASMS12kV1_lib.csv
PROVIDERS   ?= gs://$(INPUT_BUCKET)/config/Providers.csv

# Normalize: strip stray whitespace (e.g. trailing spaces or alignment padding
# in .env values), which make otherwise keeps verbatim and which would break
# unquoted command lines -- a space inside $(IMAGE) splits --tag into extra args.
# Buckets are stripped before the paths that interpolate them.
PROJECT_ID      := $(strip $(PROJECT_ID))
REGION          := $(strip $(REGION))
AR_REPO         := $(strip $(AR_REPO))
JOB_NAME        := $(strip $(JOB_NAME))
RUNTIME_SA_NAME := $(strip $(RUNTIME_SA_NAME))
INPUT_BUCKET    := $(strip $(INPUT_BUCKET))
OUTPUT_BUCKET   := $(strip $(OUTPUT_BUCKET))
CPU             := $(strip $(CPU))
MEMORY          := $(strip $(MEMORY))
TASK_TIMEOUT    := $(strip $(TASK_TIMEOUT))
MAX_RETRIES     := $(strip $(MAX_RETRIES))
PLATFORM        := $(strip $(PLATFORM))
INPUT_FILE      := $(strip $(INPUT_FILE))
OUTPUT_DIR      := $(strip $(OUTPUT_DIR))
MASTERLISTS     := $(strip $(MASTERLISTS))
PROVIDERS       := $(strip $(PROVIDERS))

IMAGE    := $(REGION)-docker.pkg.dev/$(PROJECT_ID)/$(AR_REPO)/$(JOB_NAME):latest
SA_EMAIL := $(RUNTIME_SA_NAME)@$(PROJECT_ID).iam.gserviceaccount.com
ARGS     := --input-file,$(INPUT_FILE),--output-dir,$(OUTPUT_DIR),--masterlists-dir,$(MASTERLISTS),--providers-csv,$(PROVIDERS)

.DEFAULT_GOAL := help
.PHONY: help all repo build sa iam deploy run logs pull run-local _require-project

help:
	@echo "EASMS curation -- Cloud Run Job (step by step). Set PROJECT_ID (required)."
	@echo ""
	@echo "  make repo       create the Artifact Registry repo"
	@echo "  make build      build & push the image (Cloud Build)"
	@echo "  make sa         create the runtime service account"
	@echo "  make iam        grant least-privilege bucket roles to the SA"
	@echo "  make deploy     create/update the keyless Cloud Run Job"
	@echo "  make run        execute the job once (baked-in args)"
	@echo "  make logs       list the latest execution"
	@echo "  make all        repo -> build -> sa -> iam -> deploy"
	@echo ""
	@echo "  make pull       pull the image to this machine"
	@echo "  make run-local  pull + run the image locally (uses your user ADC)"
	@echo ""
	@echo "  Image:      $(IMAGE)"
	@echo "  Runtime SA: $(SA_EMAIL)"

_require-project:
	@test -n "$(PROJECT_ID)" || { echo "ERROR: set PROJECT_ID, e.g. make deploy PROJECT_ID=my-proj"; exit 1; }

# 1. Artifact Registry repo (idempotent)
repo: _require-project
	gcloud artifacts repositories describe $(AR_REPO) --location=$(REGION) --project=$(PROJECT_ID) >/dev/null 2>&1 || gcloud artifacts repositories create $(AR_REPO) --repository-format=docker --location=$(REGION) --project=$(PROJECT_ID)

# 2. Build & push the image (Cloud Build; no local Docker needed)
build: _require-project
	gcloud builds submit --tag $(IMAGE) --project=$(PROJECT_ID) .

# 3. Dedicated runtime service account (idempotent)
sa: _require-project
	gcloud iam service-accounts describe $(SA_EMAIL) --project=$(PROJECT_ID) >/dev/null 2>&1 || gcloud iam service-accounts create $(RUNTIME_SA_NAME) --display-name="EASMS curation Cloud Run Job (runtime)" --project=$(PROJECT_ID)

# 4. Least-privilege bucket IAM (read upload/config, read+write output)
iam: _require-project
	gcloud storage buckets add-iam-policy-binding gs://$(INPUT_BUCKET) --member="serviceAccount:$(SA_EMAIL)" --role="roles/storage.objectViewer"
	gcloud storage buckets add-iam-policy-binding gs://$(OUTPUT_BUCKET) --member="serviceAccount:$(SA_EMAIL)" --role="roles/storage.objectAdmin"

# 5. Deploy the Cloud Run Job (create-or-update), keyless
deploy: _require-project
	gcloud run jobs deploy $(JOB_NAME) --image=$(IMAGE) --region=$(REGION) --project=$(PROJECT_ID) --service-account=$(SA_EMAIL) --cpu=$(CPU) --memory=$(MEMORY) --task-timeout=$(TASK_TIMEOUT) --max-retries=$(MAX_RETRIES) --set-env-vars=GOOGLE_CLOUD_PROJECT=$(PROJECT_ID),FP_N_JOBS=$(CPU),GRPC_VERBOSITY=ERROR --args="$(ARGS)"

# Execute the deployed job once with its baked-in args.
run: _require-project
	gcloud run jobs execute $(JOB_NAME) --region=$(REGION) --project=$(PROJECT_ID)

logs: _require-project
	gcloud run jobs executions list --job=$(JOB_NAME) --region=$(REGION) --project=$(PROJECT_ID) --limit=1

all: repo build sa iam deploy
	@echo "Done. Execute with: make run PROJECT_ID=$(PROJECT_ID)"

# --- Local execution (pull the AR image and run on this machine) ----------
# Authenticate Docker to AR and pull. On Apple Silicon the image is amd64, so
# PLATFORM defaults to linux/amd64 (runs under emulation).
pull: _require-project
	gcloud auth configure-docker $(REGION)-docker.pkg.dev --quiet
	docker pull --platform $(PLATFORM) $(IMAGE)

# Pull + run locally. Mounts your *user* ADC (local dev only, not production).
# Inputs are passed as env vars (the script reads them), so there are no
# duplicate flags. Extra ad-hoc flags can still be added inside the script.
run-local: _require-project
	PROJECT_ID=$(PROJECT_ID) REGION=$(REGION) AR_REPO=$(AR_REPO) JOB_NAME=$(JOB_NAME) PLATFORM=$(PLATFORM) FP_N_JOBS=$(CPU) INPUT_FILE=$(INPUT_FILE) OUTPUT_DIR=$(OUTPUT_DIR) MASTERLISTS=$(MASTERLISTS) PROVIDERS_CSV="Providers.csv" ./run_local_from_ar.sh
