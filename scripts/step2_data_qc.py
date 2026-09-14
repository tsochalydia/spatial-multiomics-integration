# %% STEP 2 — QC of both modalities + merge into one AnnData.
#
# Snakemake-adapted version of 2_data_qc.py: paths from the command line, QC
# parameters from the config YAML. The QC logic is unchanged.
#
# INPUT   --input-h5ad   concatenated Xenium AnnData from step 1a
#         --parquet-dir   ID_<slide>_intensity.parquet files from step 1b
# OUTPUT  --output-dir    final_adata.h5ad (RNA in .X/.layers["counts"], protein in
#                         .obsm["protein_expression"]) plus filtered/tracking files

import argparse
import re
import resource
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

# elapsed time + peak memory, flushed so it shows up in SLURM logs immediately
_T0 = time.time()
def log(msg):
    gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    print(f"[{time.time() - _T0:7.1f}s | peak {gb:5.1f} GB] {msg}", flush=True)


# %% ############################## ARGUMENTS #################################

def parse_args():
    p = argparse.ArgumentParser(description="Step 2: QC + modality merge")
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--input-h5ad", required=True, help="concatenated Xenium .h5ad from step 1a")
    p.add_argument("--parquet-dir", required=True, help="CODEX intensity parquets from step 1b")
    p.add_argument("--output-dir", required=True, help="output dir")
    return p.parse_args()

args = parse_args()
cfg = yaml.safe_load(open(args.config_yaml))
qc  = cfg["qc"]

INPUT_H5AD  = Path(args.input_h5ad)
PARQUET_DIR = Path(args.parquet_dir)
OUTPUT_DIR  = Path(args.output_dir)

CELL_UID  = "cell_uid"    # join key, identical in both modalities
CELL_ID   = "cell_id"     # original per-slide id, used for the Xenium Explorer CSVs
SLIDE_KEY = "slide_ID"

# Xenium cell-level thresholds
MIN_COUNTS = qc["min_counts"]
MIN_GENES  = qc["min_genes"]

# Protein markers. MARKERS = None -> every channel except EXCLUDE_CHANNELS.
EXCLUDE_CHANNELS = qc["exclude_channels"]
MARKERS          = qc["markers"]

# Per-slide background correction: subtract the Nth percentile, clip at 0.
BACKGROUND_MARKER     = qc["background_marker"]
BACKGROUND_PERCENTILE = qc["background_percentile"]

# CODEX intensity quantile bounds, computed per slide
LOWER_Q = qc["lower_q"]
UPPER_Q = qc["upper_q"]
MAX_MARKERS_ABOVE_UPPER_Q = qc["max_markers_above_upper_q"]

MARKER_SUFFIX = qc["marker_suffix"]

# Columns a parquet may carry that are not markers (schemas differ between 1b versions)
PARQUET_META_COLS = {CELL_UID, CELL_ID, "cell_labels", SLIDE_KEY, "area_px", "sample_id", "run"}
PARQUET_RE = re.compile(r"^ID_(?P<slide_ID>\d+)_intensity\.parquet$")

# OUTPUT FILE PATHS
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FILTERED_ADATA_OUT           = OUTPUT_DIR / "filtered_adata.h5ad"
CODEX_RAW_OUT                = OUTPUT_DIR / "codex_all_cells_all_markers.parquet"
CODEX_INTENSITY_FILTERED_OUT = OUTPUT_DIR / "codex_all_cells_intensity_filtered.parquet"
MERGED_ADATA_OUT             = OUTPUT_DIR / "final_adata.h5ad"


# %% ################################ HELPERS ##################################

def record_counts(tracking, step_name, data, slide_key=SLIDE_KEY):
    """Record total + per-slide cell counts at one QC step."""
    obs = data.obs if isinstance(data, ad.AnnData) else data
    n_total = data.n_obs if isinstance(data, ad.AnnData) else len(data)

    counts = {"total": n_total}
    if slide_key in obs.columns:
        counts.update(obs[slide_key].value_counts().to_dict())
    tracking[step_name] = counts
    return tracking


def write_cell_group_csvs(cell_group_df, output_dir):
    """One CSV per slide: original cell_id -> group. Feeds Xenium Explorer later."""
    out_folder = Path(output_dir) / "cell_id_x_group"
    out_folder.mkdir(parents=True, exist_ok=True)

    for slide_id, grp in cell_group_df.groupby(SLIDE_KEY, observed=True):
        out = grp[[CELL_ID, "group"]]
        n_dupes = int(out[CELL_ID].duplicated().sum())
        if n_dupes:
            print(f"  [WARNING] {slide_id}: {n_dupes} duplicate cell_ids")
        fname = f"cell_id_x_group_{slide_id}.csv"
        out.to_csv(out_folder / fname, index=False)
        print(f"  {slide_id}: {out['group'].value_counts().to_dict()}  ->  {fname}")

    print(f"[write_cell_group_csvs] Done — files written to {out_folder}")


# ----------------------------------------------------------------------
# Step 1: Xenium QC
# ----------------------------------------------------------------------
def load_xenium(input_path):
    adata = ad.read_h5ad(input_path)

    if adata.obs_names.duplicated().any():
        raise ValueError("Duplicated obs_names in the concatenated Xenium object.")
    adata.obs[CELL_UID] = adata.obs_names.astype(str)

    for col in (CELL_ID, SLIDE_KEY):
        if col not in adata.obs:
            raise KeyError(f"Missing '{col}' in .obs; columns: {adata.obs.columns.tolist()}")
    adata.obs[SLIDE_KEY] = adata.obs[SLIDE_KEY].astype(str)
    return adata


def filter_cells(adata, min_counts, min_genes, output_path=None):
    sc.pp.filter_cells(adata, min_counts=min_counts)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if output_path:
        adata.write_h5ad(output_path)
    return adata


# ----------------------------------------------------------------------
# Step 2: CODEX parquet loading + marker cleanup
# ----------------------------------------------------------------------
def correct_marker_background(df, marker, slide_key=SLIDE_KEY, percentile=10, apply=True):
    """Per-slide background subtraction: value - Nth percentile, clipped at 0."""
    if not apply or marker is None:
        return df, marker

    if marker not in df.columns:
        raise KeyError(
            f"correct_marker_background: marker '{marker}' not found. "
            f"Available: {[c for c in df.columns if c not in PARQUET_META_COLS]}"
        )

    df        = df.copy()
    raw       = df[marker].to_numpy(dtype=float)
    corrected = raw.copy()

    print(f"[correct_marker_background] Per-slide correction: '{marker}' (percentile={percentile})")

    for slide, grp_idx in df.groupby(slide_key, observed=True).groups.items():
        pos        = df.index.get_indexer(grp_idx)
        slide_vals = raw[pos]
        valid_vals = slide_vals[~np.isnan(slide_vals)]

        if len(valid_vals) == 0:
            print(f"  slide {slide}: no valid values — skipping")
            continue

        bg              = np.percentile(valid_vals, percentile)
        slide_corrected = np.clip(slide_vals - bg, 0, None)
        corrected[pos]  = slide_corrected

        n_clipped = int((slide_corrected == 0).sum())
        print(
            f"  slide {slide}: bg ({percentile}th pct) = {bg:.4f} | "
            f"n_cells = {len(slide_vals)} | "
            f"clipped to 0: {n_clipped} ({100 * n_clipped / len(slide_vals):.1f}%)"
        )

    corrected_col = marker + "_cor"
    insert_pos    = df.columns.get_loc(marker)
    df.insert(insert_pos, corrected_col, corrected)
    df = df.drop(columns=[marker])

    print(f"[correct_marker_background] Done — '{marker}' replaced by '{corrected_col}'")
    return df, corrected_col


def load_codex_parquets(folder, markers=None, exclude_channels=(), raw_output_path=None):
    """
    Read every ID_<slide_ID>_intensity.parquet, stack them, and return
    (codex_df, marker_list). slide_ID comes from the filename.
    """
    folder = Path(folder)

    frames = []
    for f in sorted(folder.glob("*.parquet")):
        m = PARQUET_RE.match(f.name)
        if m is None:
            print(f"  [WARNING] skipping unrecognised parquet name: {f.name}")
            continue

        slide_id = m.group("slide_ID")
        df = pd.read_parquet(f)

        if CELL_UID not in df.columns:
            raise KeyError(f"{f.name}: no '{CELL_UID}' column; got {df.columns.tolist()}")

        df[SLIDE_KEY] = slide_id     # authoritative, overrides any stale column
        frames.append(df)
        log(f"  {f.name}: {len(df)} cells")

    if not frames:
        raise ValueError(f"No ID_<slide>_intensity.parquet found in {folder}")

    # join="outer": a channel missing from one slide becomes a NaN column there
    codex_df = pd.concat(frames, ignore_index=True, join="outer")

    if codex_df[CELL_UID].duplicated().any():
        raise ValueError(f"{int(codex_df[CELL_UID].duplicated().sum())} duplicated {CELL_UID}")

    channels = [c for c in codex_df.columns if c not in PARQUET_META_COLS]
    if markers is None:
        markers = [c for c in channels if c not in set(exclude_channels)]
    else:
        missing = [m for m in markers if m not in channels]
        if missing:
            raise KeyError(f"Requested markers not in the parquets: {missing}")

    print(f"\n{len(channels)} channels found, {len(markers)} kept as markers")
    print(f"  excluded: {sorted(set(channels) - set(markers))}")

    keep_meta = [c for c in (CELL_UID, CELL_ID, "cell_labels", SLIDE_KEY) if c in codex_df.columns]
    codex_df  = codex_df[keep_meta + markers]

    if raw_output_path:
        codex_df.to_parquet(raw_output_path, index=False)

    return codex_df, markers


# ----------------------------------------------------------------------
# Step 2b: Per-slide protein-intensity filtering
# ----------------------------------------------------------------------
def filter_protein_intensities(codex_df, markers, upper_q=1.0, lower_q=0.0,
                               min_markers_exceeding=0, slide_key=SLIDE_KEY,
                               tracking=None, output_path=None):
    missing = [m for m in markers if m not in codex_df.columns]
    if missing:
        raise KeyError(
            f"filter_protein_intensities: missing marker columns {missing}. "
            "Pass the marker list returned by load_codex_parquets so corrected "
            "names (e.g. 'CD4_cor') are used."
        )

    codex_df = codex_df.copy()

    if tracking is not None:
        record_counts(tracking, "codex_before_intensity_filter", codex_df, slide_key)

    # Low end: squash noise to 0
    lower_thresh_df = codex_df.groupby(slide_key, observed=True)[markers].quantile(lower_q)
    for m in markers:
        row_lower = codex_df[slide_key].map(lower_thresh_df[m])
        codex_df.loc[codex_df[m] < row_lower, m] = 0

    # High end: drop cells that are extreme in many markers at once
    upper_thresh_df = codex_df.groupby(slide_key, observed=True)[markers].quantile(upper_q)
    exceed = pd.DataFrame(index=codex_df.index)
    for m in markers:
        row_upper = codex_df[slide_key].map(upper_thresh_df[m])
        exceed[m] = codex_df[m] > row_upper

    n_exceeding = exceed.sum(axis=1)
    keep_mask   = n_exceeding < min_markers_exceeding
    codex_df    = codex_df.loc[keep_mask].copy()

    print(f"[filter_protein_intensities] removed {int((~keep_mask).sum())} cells "
          f"exceeding q{upper_q} in >= {min_markers_exceeding} markers")

    if tracking is not None:
        record_counts(tracking, "codex_after_intensity_filter", codex_df, slide_key)

    if output_path:
        codex_df.to_parquet(output_path, index=False)

    thresholds = {
        "lower_thresholds": lower_thresh_df.to_dict(orient="index"),
        "upper_thresholds": upper_thresh_df.to_dict(orient="index"),
    }
    return codex_df, thresholds


# ----------------------------------------------------------------------
# Step 3: Integrate both modalities into one AnnData
# ----------------------------------------------------------------------
def add_protein_expression(adata, codex_df, markers, marker_suffix=MARKER_SUFFIX,
                           obsm_key="protein_expression", drop_unmatched=True,
                           tracking=None):
    codex_indexed = codex_df.set_index(CELL_UID)

    matched = adata.obs_names.isin(codex_indexed.index)
    print(
        f"Matched cells between Xenium and CODEX after filtering: "
        f"{int(matched.sum())} / {adata.n_obs} ({matched.mean() * 100:.1f}%)"
    )

    if tracking is not None:
        record_counts(tracking, "xenium_before_codex_match", adata)

    if drop_unmatched:
        adata = adata[matched].copy()

    aligned = codex_indexed.reindex(adata.obs_names)

    obs_marker_names            = [m + marker_suffix for m in markers]
    adata.obs[obs_marker_names] = aligned[markers].to_numpy()
    adata.obsm[obsm_key]        = aligned[markers].to_numpy(dtype=np.float32)
    adata.uns["protein_markers"] = list(obs_marker_names)

    if tracking is not None:
        record_counts(tracking, "final_after_codex_match", adata)

    return adata


# %% ############################## ORCHESTRATION ##############################

if __name__ == "__main__":
    tracking = {}

    # Load the concatenated Xenium object
    log("loading Xenium")
    adata = load_xenium(INPUT_H5AD)
    log(f"loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    record_counts(tracking, "xenium_before_cell_filter", adata)

    # Snapshot every cell before any filtering, for the group tracking below
    all_cells = adata.obs[[CELL_UID, CELL_ID, SLIDE_KEY]].copy().reset_index(drop=True)

    # Xenium QC filter
    adata = filter_cells(adata, MIN_COUNTS, MIN_GENES, output_path=FILTERED_ADATA_OUT)
    log(f"after Xenium QC: {adata.n_obs} cells")
    record_counts(tracking, "xenium_after_cell_filter", adata)

    xenium_surviving = set(adata.obs_names.astype(str))

    # Load CODEX
    log("loading CODEX parquets")
    codex_df, markers = load_codex_parquets(
        PARQUET_DIR,
        markers=MARKERS,
        exclude_channels=EXCLUDE_CHANNELS,
        raw_output_path=CODEX_RAW_OUT,
    )
    log(f"CODEX: {len(codex_df)} cells x {len(markers)} markers")

    # Step 2b: per-slide background correction
    codex_df, corrected_marker = correct_marker_background(
        codex_df,
        marker=BACKGROUND_MARKER,
        percentile=BACKGROUND_PERCENTILE,
        apply=BACKGROUND_MARKER is not None,
    )
    if corrected_marker != BACKGROUND_MARKER:
        markers = [corrected_marker if m == BACKGROUND_MARKER else m for m in markers]

    # Step 2c: intensity filtering
    codex_df_filtered, intensity_thresholds = filter_protein_intensities(
        codex_df,
        markers,
        upper_q=UPPER_Q,
        lower_q=LOWER_Q,
        min_markers_exceeding=MAX_MARKERS_ABOVE_UPPER_Q,
        tracking=tracking,
        output_path=CODEX_INTENSITY_FILTERED_OUT,
    )
    log(f"CODEX after intensity filter: {len(codex_df_filtered)} cells")

    # Step 3: merge the two modalities
    adata_merged = add_protein_expression(
        adata,
        codex_df_filtered,
        markers,
        marker_suffix=MARKER_SUFFIX,
        drop_unmatched=True,
        tracking=tracking,
    )
    log(f"writing {MERGED_ADATA_OUT}")
    adata_merged.write_h5ad(MERGED_ADATA_OUT)
    log(f"saved -> {MERGED_ADATA_OUT} ({MERGED_ADATA_OUT.stat().st_size / 1e9:.2f} GB)")

    # Step 4: three mutually exclusive groups covering every original cell
    final_ids = set(adata_merged.obs_names.astype(str))
    all_cells["group"] = "xenium_qc_removed"
    all_cells.loc[all_cells[CELL_UID].isin(xenium_surviving), "group"] = "codex_qc_removed"
    all_cells.loc[all_cells[CELL_UID].isin(final_ids),        "group"] = "kept"
    write_cell_group_csvs(all_cells, OUTPUT_DIR)

    # Step 5: count tracking summary
    tracking_df = pd.DataFrame(tracking).T
    print("\n", tracking_df)
    tracking_df.to_csv(OUTPUT_DIR / "tracking_counts.csv")

    # %% checks
    print("\nfinal object:")
    print(adata_merged)
    print("\nprotein markers:", adata_merged.uns["protein_markers"])
    print("\ncells per slide:")
    print(adata_merged.obs[SLIDE_KEY].value_counts().sort_index().to_string())
    print("\nNaN fraction per marker:")
    prot = pd.DataFrame(adata_merged.obsm["protein_expression"],
                        columns=adata_merged.uns["protein_markers"])
    print(prot.isna().mean().to_string())
    print("\ncounts layer present:", "counts" in adata_merged.layers)
    log("done")
