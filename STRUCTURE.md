# Project Structure (code)

A quick map of the EASMS data-processing codebase. The pipeline is
**step-sequential with disk checkpoints**: each numbered module writes its own
`StepN_*/` folder, so any step can resume from a saved checkpoint via
`--start-from` / `--end-at`. Steps 1–6 write CSV; Steps 7–8 switch to Parquet
(fingerprints make CSVs huge).

```
EASMS-data-processing/
│
├── src/
│   ├── Main.py                       # Entry point — orchestrates Steps 0–8 per input CSV,
│   │                                 #   step gating (--start-from/--end-at), resume-from-disk, CLI args
│   ├── io_utils.py                   # Filesystem helpers that work on local disk AND gs:// buckets
│   │                                 #   (fsspec/gcsfs); pjoin/listdir/makedirs wrappers
│   │
│   ├── quality_check.py              # Step 0 — pre-processing QC on raw input; per-check funcs,
│   │                                 #   logs grouped by section (+ .log / color-coded .xlsx)
│   ├── separate_protein_files.py     # Step 1 — split raw CSV into one file per TARGET_ID
│   ├── add_scores.py                 # Step 2 — TARGET_VALUE, SELECTIVE/NTC_VALUE, ENRICHMENT,
│   │                                 #   EASMS_ENRICHMENT, MEAN_NONTARGET_VALUES, Welch's PVALUE
│   ├── anomaly_selection.py          # Step 3 — resolve duplicate-SMILES rows w/ conflicting enrichment
│   ├── isomer_handling.py            # Step 4 — split multi-isomer SMILES (';') into one row each
│   ├── add_negatives.py              # Step 5 — pull negatives from the library masterlist (by LIBRARY_NAME)
│   ├── produce_ml_labels.py          # Step 6 — assign AIRCHECK_LABEL (−2..4) from enrichment/pvalue/isomers
│   ├── fingerprint_extraction.py     # Step 7 — RDKit fingerprints/descriptors + renames + binary LABEL
│   ├── fingerprints.py               # FP builders: ECFP4/6, FCFP4/6, MACCS, RDK, AVALON, TOPTOR, ATOMPAIR
│   ├── utils.py                      #   RDKit exception handling + MW/ALOGP descriptor helpers
│   ├── column_selection.py           # Step 8 — split columns into data (Parquet) vs metadata (CSV)
│   │                                 #   driven by ColumnActions.xlsx tags
│   ├── post_quality_check.py         # Post-pipeline QC on Step 7 output (23 checks, ranges/labels/FP lengths)
│   ├── batch_protein_info.py         # Per-batch protein summary — proteins_in_batch.csv (1 row/TARGET_ID)
│   └── fingerprint_extraction_original.py   # Archived earlier version of the FP extractor
│
├── multibatchcodes/                  # Ad-hoc multi-batch analysis & plotting (not part of the pipeline)
│   ├── multi_batch_analysis.py       #   cross-batch aggregation vs SPR gold labels
│   ├── compute_enrichment*.py        #   enrichment variants (median / percentile / normalized)
│   ├── compute_bead_ratio.py         #   bead-ratio + replicate consistency/dropout diagnostics
│   └── plot_regime_comparison*.py    #   regime-comparison / hit-matrix plots
│
├── MasterLists/                      # Library files (negative-sample compound sources)
├── ColumnActions.xlsx                # Tag file: which columns -> data vs metadata in Step 8
├── Providers.csv / RawDataColumns.csv  # Canonical provider + raw-column references (QC uses these)
│
├── PIPELINE.md                       # Step-by-step spec of the 8 pipeline steps
├── QUALITY_CHECKS.md / POST_QC.md    # QC check catalogs (pre + post)
└── USAGE.md / Readme.md              # Run instructions & flags
```

## Pipeline flow at a glance

| Step | Module | Output |
|------|--------|--------|
| 0 | `quality_check.py` | QC log + `.xlsx` (gate; failures skip the file) |
| 1 | `separate_protein_files.py` | `Step1_Separated/` — one CSV per `TARGET_ID` |
| 2 | `add_scores.py` | `Step2_WithScores/` — enrichment + p-value columns |
| 3 | `anomaly_selection.py` | `Step3_AnomalyFiltered/` |
| 4 | `isomer_handling.py` | `Step4_IsomerHandled/` |
| 5 | `add_negatives.py` | `Step5_WithNegatives/` |
| 6 | `produce_ml_labels.py` | `Step6_MLReady/` — `AIRCHECK_LABEL` (last CSV step) |
| 7 | `fingerprint_extraction.py` | `Step7_WithFingerprints/` — Parquet, RDKit FPs + `LABEL` |
| 8 | `column_selection.py` | `Step8_FinalData/` (Parquet) + `Step8_Metadata/` (CSV) |
| post | `post_quality_check.py` | `PostQClog_*` — validates Step 7 output |

See [PIPELINE.md](PIPELINE.md) for the full per-step specification.
