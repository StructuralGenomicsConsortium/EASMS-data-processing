# multibatchcodes

Analysis scripts for the multi-batch E-ASMS data. Only the code is tracked here —
all data, `*_modified/` folders and Excel output stay local (see `.gitignore`).

## The two-step workflow

Everything else is exploratory. For a new set of batches you only need these two,
in this order:

```bash
# 1. compute every analysis column   (raw CSVs  ->  MultiBatch_<name>_modified/)
python multibatchcodes/run_multibatch_pipeline.py --input-dir MultiBatch_<name> \
       --skip-heatmap --skip-matrix

# 2. build the decision report       (_modified  ->  MultiBatch_AcrossBatchDecisions/)
python multibatchcodes/acrossbatch_report.py --dir MultiBatch_<name>_modified
```

Step 2 needs the columns from step 1 and exits with a clear message if they are
missing. Drop `--skip-heatmap --skip-matrix` in step 1 if you also want the hit
matrix PNG and the full matrix CSV/XLSX.

Before running: put one CSV per protein in `MultiBatch_<name>/` (needs
`COMPOUND_ID`, `ASMS_BATCH_NAME`, `POS_INT_REP1/2/3`), keep
`MultiBatchResults/beads_clean.xlsx` in place for the bead columns, and close any
output file that is open in Excel — locked files are skipped and reported.

Timing is dominated by step 1: roughly 10 min for 24 proteins, ~50 min for 136
(the across-batch Mann-Whitney grows with proteins²).

## What step 1 adds

Beyond the older `enrichment_*` / `_norm` / `_pctnorm` families, it appends the
`ns1_*` chain and the across-batch 2×2 as the last columns:

| column | scaling | across-batch background |
|---|---|---|
| `enrich_across_medscale_medbg` | per-batch median | median |
| `enrich_across_medscale_meanbg` | per-batch median | mean |
| `enrich_across_quantscale_medbg` | quantile | median |
| `enrich_across_quantscale_meanbg` | quantile | mean |
| `Normalized score1` | — | robust Z on `log2(E_within)` |

`Normalized score1` = scale all values (3000s included) → median of 3 replicates →
`E = target / mean(non-target medians)` within batch → `log2` → robust Z against
the other proteins carrying the same compound.

Two things worth knowing:

- **The MAD floor is re-derived per dataset**, not hard-coded: the 5th percentile
  of the leave-one-out MAD among compounds detected on ≥75 % of their proteins
  (0.1195 on 20to22, 0.2133 on 20to36 — it really does move). Without it the score
  is undefined for ~44 % of rows and blank for *every* perfectly selective binder.
  Override with `--mad-min`.
- **Per-batch scaling cancels inside a within-batch ratio**, so it changes
  `ns1_enrichment_within` and `Normalized score1` not at all. It only matters for
  the across-batch columns.

## What step 2 produces

```
MultiBatch_AcrossBatchDecisions/
    20to22/
        MultiBatch_20to22_acrossbatch_output.xlsx   combined, one sheet per protein
        per_protein/BRD1.xlsx, BRD9.xlsx, …         the same sheets, split out
    20to36/
        …
```

Per protein: top 200 by each of the five metrics, merged to unique compounds
(`selected_by` / `n_metrics` record which metrics found each row), then split into
three blocks written in this order:

| block | fill | tag |
|---|---|---|
| A new candidate | none | replicate profile, e.g. `3 moderate`, `2 moderate + 1 weak` — sorted by `Normalized score1` |
| B floor / bead | light blue | `has floor` / `bead binder` / `has floor and bead binder` |
| C already selected | gray | `Selected before based on within batch fold number` (`BINARY_LABEL = 1`), kept at the bottom |

Thresholds match the pipeline: strong ≥ 1e6, moderate ≥ 1e4, weak > 3000.
Useful flags: `--top` (default 200), `--no-split`, `--min-detected 1` to stop
all-floored rows competing for a slot.

## Everything else

One-off exploration, kept for reference — the floor/censoring investigation
(`floor_leak_theory.py`, `compare_floor_definitions.py`,
`compare_floor_estimators.py`, `fit_censored_floor.py`,
`plot_intensity_histograms_floor_peak.py`), batch-effect diagnostics
(`protein_median_by_batch.py`, `plot_protein_median_by_batch.py`), scoring
analysis (`analyze_normalized_score1.py`), and the earlier per-step scripts that
`run_multibatch_pipeline.py` now supersedes (`compute_enrichment*.py`,
`hit_matrix*.py`, `export_matrix*.py`).
