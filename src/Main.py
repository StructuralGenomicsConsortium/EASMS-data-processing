# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 09:24:34 2025

@author: shay
test
"""



import os
import argparse
from datetime import datetime
import pandas as pd
from io_utils import pjoin, listdir, makedirs, dirname, isdir
from separate_protein_files import split_protein_data
from add_scores import compute_and_add_scores
from anomaly_selection import filter_anomalous_data
from isomer_handling import handle_isomers
from produce_ml_labels import generate_ml_labels
from add_negatives import add_negative_samples_from_masterlist
from fingerprint_extraction import extract_fingerprints
from column_selection import select_final_columns, load_column_tags
from quality_check import run_quality_checks
from batch_protein_info import write_batch_protein_info
from post_quality_check import run_post_quality_checks

STEP_FOLDERS = {
    1: "Step1_Separated",
    2: "Step2_WithScores",
    3: "Step3_AnomalyFiltered",
    4: "Step4_IsomerHandled",
    5: "Step5_WithNegatives",
    6: "Step6_MLReady",
    7: "Step7_WithFingerprints",
    8: "Step8_FinalData",
}
# Step 8 also writes a per-target metadata CSV (the 'metadata'-tagged columns)
# into this sibling folder.
METADATA_FOLDER = "Step8_Metadata"

# Steps 1-6 save CSV; steps 7-8 save Parquet (fingerprints make CSVs unwieldy).
PARQUET_FROM_STEP = 7


def _step_input_files(start_from, step_dirs, scored_files):
    """Return the list of files to feed into the per-target loop based on start_from.

    Step 8 branches from Step 7's output, so when resuming at 8 we read from
    Step7, not the immediately preceding step.
    """
    if start_from <= 3:
        return scored_files
    load_step = 7 if start_from == 8 else (start_from - 1)
    load_dir = step_dirs[load_step]
    ext = ".parquet" if load_step >= PARQUET_FROM_STEP else ".csv"
    return sorted(
        pjoin(load_dir, f) for f in listdir(load_dir) if f.endswith(ext)
    )


def _unique_protein_id_name(df, fallback):
    """Filename stem for the Step 7/8 outputs, taken from the UNIQUE_PROTEIN_ID
    column (one protein/target per file). Falls back to `fallback` if the column
    is missing or has no value. Characters that are invalid in filenames are
    sanitized — notably the '|' used in protein-complex IDs ('|' -> '+')."""
    if "UNIQUE_PROTEIN_ID" not in df.columns:
        return fallback
    vals = df["UNIQUE_PROTEIN_ID"].dropna()
    if vals.empty:
        return fallback
    name = str(vals.iloc[0]).strip()
    if not name:
        return fallback
    name = name.replace("|", "+")
    for ch in '\\/:*?"<>':
        name = name.replace(ch, "_")
    return name


def process_csv_files(input_files, masterlist_path, output_dir, providers_csv,
                      data_columns, metadata_columns,
                      start_from=0, end_at=8, meta_csv=None, fp_format="array"):
    """Processes the given CSV files through data curation steps 1..8.

    `input_files` is a list of raw CSV paths (local or gs://). For each file, a
    `ProcessedData_<name>/` folder is created under `output_dir`.

    start_from / end_at gate which steps execute. Steps not run are loaded from
    their saved output on disk, so the pipeline can resume from any checkpoint.
    Step 0 (QC) only runs when start_from <= 0.
    """
    for file_path in input_files:
        file_name = os.path.basename(file_path)
        if not file_name.endswith(".csv"):
            continue

        print(f"Processing: {file_name}")

        csv_basename = os.path.splitext(file_name)[0]
        processed_data_dir = pjoin(output_dir, f"ProcessedData_{csv_basename}")
        makedirs(processed_data_dir, exist_ok=True)

        step_dirs = {n: pjoin(processed_data_dir, name) for n, name in STEP_FOLDERS.items()}

        # Per-batch protein summary (one row per TARGET_ID), written next to the
        # step folders. Independent of the step gating below.
        write_batch_protein_info(file_path, processed_data_dir)

        # Step 0: QC (only when start-from == 0)
        if start_from <= 0:
            qc_passed, qc_failures = run_quality_checks(
                file_path,
                processed_data_dir,
                providers_csv=providers_csv,
                masterlist_dir=masterlist_path,
                meta_csv=meta_csv,
            )
            if not qc_passed:
                log_path = pjoin(processed_data_dir, f"QCaircheck{datetime.now().strftime('%Y%m%d')}_{csv_basename}.log")
                bar = "=" * 70
                print()
                print(bar)
                print(f"  Quality checks FAILED for: {file_name}")
                print(bar)
                print(f"  {len(qc_failures)} check(s) did not pass:")
                print()
                for desc, msg in qc_failures:
                    print(f"    - {desc}")
                    print(f"        {msg}")
                print()
                print(f"  Full log: {log_path}")
                print(f"  Skipping {file_name} -- no downstream steps will run.")
                print(bar)
                print()
                continue
            print(f"  Quality checks PASSED for {file_name}.")

        # If end_at < 1 the user only wants QC; nothing else needs Step 1 output.
        if end_at < 1:
            continue

        # ---- Step 1: split by target ----
        if start_from <= 1:
            makedirs(step_dirs[1], exist_ok=True)
            separated_files = split_protein_data(file_path, step_dirs[1])
        else:
            separated_files = sorted(
                pjoin(step_dirs[1], f) for f in listdir(step_dirs[1]) if f.endswith(".csv")
            )
        if end_at < 2:
            continue

        # ---- Step 2: compute scores (batch over all separated files) ----
        if start_from <= 2:
            makedirs(step_dirs[2], exist_ok=True)
            print("\nComputing and Adding Scores to All Separated Files...\n")
            scored_files = compute_and_add_scores(separated_files, output_dir=step_dirs[2])
        else:
            scored_files = sorted(
                pjoin(step_dirs[2], f) for f in listdir(step_dirs[2]) if f.endswith(".csv")
            )
        if end_at < 3:
            continue

        # ---- Steps 3..9: per-target loop ----
        entry_files = _step_input_files(start_from, step_dirs, scored_files)

        for input_file in entry_files:
            base_name = os.path.splitext(os.path.basename(input_file))[0]
            print(f"  Processing separated file: {base_name}")

            if input_file.endswith(".parquet"):
                df = pd.read_parquet(input_file)
            else:
                df = pd.read_csv(input_file)

            # Steps 1-6 are named by TARGET_ID (base_name); Steps 7-8 are named
            # by UNIQUE_PROTEIN_ID. Derived from the real rows here, before Step 5
            # adds negative samples (which have no UNIQUE_PROTEIN_ID).
            out_name = _unique_protein_id_name(df, base_name)

            # Step 3: filter anomalies
            if start_from <= 3 and end_at >= 3:
                print(f"    Step 3: filtering anomalies — {base_name}", flush=True)
                makedirs(step_dirs[3], exist_ok=True)
                df = filter_anomalous_data(df, f"{base_name}.csv")
                df.to_csv(pjoin(step_dirs[3], f"{base_name}.csv"), index=False)

            # Step 4: handle isomers
            if start_from <= 4 and end_at >= 4:
                print(f"    Step 4: handling isomers — {base_name}", flush=True)
                makedirs(step_dirs[4], exist_ok=True)
                df = handle_isomers(df, f"{base_name}.csv")
                df.to_csv(pjoin(step_dirs[4], f"{base_name}.csv"), index=False)

            # Step 5: add negative samples
            if start_from <= 5 and end_at >= 5:
                print(f"    Step 5: adding negative samples — {base_name}", flush=True)
                makedirs(step_dirs[5], exist_ok=True)
                df = add_negative_samples_from_masterlist(df, file_name, masterlist_path)
                df.to_csv(pjoin(step_dirs[5], f"{base_name}.csv"), index=False)

            # Step 6: generate ML labels (last CSV step)
            if start_from <= 6 and end_at >= 6:
                print(f"    Step 6: generating ML labels — {base_name}", flush=True)
                makedirs(step_dirs[6], exist_ok=True)
                df = generate_ml_labels(df)
                df.to_csv(pjoin(step_dirs[6], f"{base_name}.csv"), index=False)

            # Step 7: extract fingerprints + rename + add binary LABEL (first Parquet step)
            if start_from <= 7 and end_at >= 7:
                print(f"    Step 7: extracting fingerprints — {base_name} ({len(df)} rows)...", flush=True)
                makedirs(step_dirs[7], exist_ok=True)
                df = extract_fingerprints(df, fp_format=fp_format)
                df = df.rename(columns={
                    "TARGET_VALUE": "TARGET_INTENSITY_VALUE",
                    "MEAN_NONTARGET_VALUES": "NONTARGET_INTENSITY_VALUE",
                })
                df["LABEL"] = df["BINARY_LABEL"]
                df.to_parquet(pjoin(step_dirs[7], f"{out_name}.parquet"), index=False)

            # Step 8: split the Step 7 dataframe into the final data file
            # ('data'-tagged columns, Parquet) and a metadata sidecar
            # ('metadata'-tagged columns, CSV). select_final_columns returns a
            # new DataFrame, so `df` stays intact between the two calls.
            if start_from <= 8 and end_at >= 8:
                print(f"    Step 8: writing final data + metadata — {base_name}", flush=True)
                makedirs(step_dirs[8], exist_ok=True)
                df_data = select_final_columns(df, data_columns)
                df_data.to_parquet(pjoin(step_dirs[8], f"{out_name}.parquet"), index=False)

                metadata_dir = pjoin(processed_data_dir, METADATA_FOLDER)
                makedirs(metadata_dir, exist_ok=True)
                df_meta = select_final_columns(df, metadata_columns)
                df_meta.to_csv(pjoin(metadata_dir, f"{out_name}.csv"), index=False)

        # ---- Post-pipeline QC ----
        # Runs once per input CSV, against the full Step 7 per-target Parquet
        # output (which still carries every column). Validates value sets,
        # ranges, non-negativity, and fingerprint lengths.
        if end_at >= 7:
            run_post_quality_checks(
                parquet_dir=step_dirs[7],
                log_dir=processed_data_dir,
                csv_basename=csv_basename,
            )


def main(input_files, masterlist_path, output_dir, providers_csv,
         data_columns, metadata_columns,
         start_from=0, end_at=8, meta_csv=None, fp_format="array"):
    """Main function to execute the full data curation pipeline."""
    process_csv_files(input_files, masterlist_path, output_dir, providers_csv,
                      data_columns, metadata_columns,
                      start_from=start_from, end_at=end_at, meta_csv=meta_csv,
                      fp_format=fp_format)

if __name__ == "__main__":
    # Define paths (Modify as needed)
    parser = argparse.ArgumentParser(description="Run EASMS data processing pipeline.")
    # Exactly one of --input-file / --input-dir selects the raw data. Both local
    # paths and gs:// URLs are accepted.
    parser.add_argument(
        "--input-file",
        default=None,
        help="A single raw ASMS CSV to process (local path or gs:// URL). Its ProcessedData_<name>/ folder is created in the file's own folder unless --output-dir is given.",
    )
    parser.add_argument(
        "--input-dir",
        default=None,
        help="A folder of raw ASMS CSVs to process (every *.csv in it, non-recursive). Alternative to --input-file.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Where ProcessedData_*/ folders are created. Default: the folder of --input-file, or --input-dir itself.",
    )
    # Config / reference inputs. Default to the local repo (this checkout) so
    # they don't have to be copied into every data location; override to point
    # at a shared folder or a gs:// path.
    parser.add_argument(
        "--masterlists-dir",
        default=None,
        help="Folder of master-list (library) files. Default: <repo>/MasterLists.",
    )
    parser.add_argument(
        "--providers-csv",
        default=None,
        help="Providers.csv path. Default: <repo>/Providers.csv.",
    )
    parser.add_argument(
        "--meta-csv",
        default=None,
        help="RawDataColumns.csv (canonical column reference) path. Default: <repo>/RawDataColumns.csv.",
    )
    parser.add_argument(
        "--column-actions",
        default=None,
        help="ColumnActions.xlsx tag file (columns 'Column name'/'Action' with values data/metadata/-). "
             "'data' columns -> Step8 final-data Parquet; 'metadata' columns -> Step8 metadata CSV. "
             "Default: <repo>/ColumnActions.xlsx.",
    )
    parser.add_argument(
        "--start-from",
        type=int,
        default=0,
        choices=range(0, 9),
        metavar="N",
        help="Step number to start running from (0-8). 0 = Quality Check. Earlier steps are loaded from their saved output. Default: 0.",
    )
    parser.add_argument(
        "--end-at",
        type=int,
        default=8,
        choices=range(0, 9),
        metavar="N",
        help="Step number to stop after (0-8). 0 = run only Quality Check, then stop. Later steps are skipped. Default: 8.",
    )
    args = parser.parse_args()
    if args.start_from > args.end_at:
        parser.error(f"--start-from ({args.start_from}) cannot be greater than --end-at ({args.end_at}).")
    if bool(args.input_file) == bool(args.input_dir):
        parser.error("provide exactly one of --input-file or --input-dir.")

    # Resolve the raw input file list and the default output location.
    if args.input_file:
        input_files = [args.input_file]
        default_output_dir = dirname(args.input_file)
    else:
        if not isdir(args.input_dir):
            parser.error(
                f"--input-dir is not a folder: {args.input_dir}. "
                "To process a single file, use --input-file instead."
            )
        input_files = sorted(
            pjoin(args.input_dir, f) for f in listdir(args.input_dir) if f.endswith(".csv")
        )
        if not input_files:
            parser.error(f"no .csv files found in --input-dir: {args.input_dir}")
        default_output_dir = args.input_dir
    output_dir = args.output_dir if args.output_dir else default_output_dir

    # Config/reference inputs default to the local repo (parent of this src/ dir).
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    masterlist_path = args.masterlists_dir if args.masterlists_dir else pjoin(repo_root, "MasterLists")
    providers_csv = args.providers_csv if args.providers_csv else pjoin(repo_root, "Providers.csv")
    meta_csv = args.meta_csv if args.meta_csv else pjoin(repo_root, "RawDataColumns.csv")
    column_actions = args.column_actions if args.column_actions else pjoin(repo_root, "ColumnActions.xlsx")

    # How fingerprint values are stored in the Step 7+ output columns:
    #   "array"  (default) -> numpy float32 arrays, ready to use directly.
    #   "string"            -> comma-separated strings (legacy format from earlier
    #                          pipeline versions; consumers must np.fromstring back
    #                          to arrays before use).
    TypeOfFp = "array"

    # Step 8 column selection is driven by the ColumnActions tag file:
    #   'data'     -> final-data Parquet (Step8_FinalData)
    #   'metadata' -> metadata CSV       (Step8_Metadata)
    #   anything else (e.g. '-') is dropped.
    data_columns, metadata_columns = load_column_tags(column_actions)
    print(f"Loaded column tags from {column_actions}: "
          f"{len(data_columns)} data, {len(metadata_columns)} metadata columns.")

    main(input_files, masterlist_path, output_dir, providers_csv,
         data_columns, metadata_columns,
         start_from=args.start_from, end_at=args.end_at, meta_csv=meta_csv,
         fp_format=TypeOfFp)
