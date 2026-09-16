# %% STEP 1b — Per-cell CODEX intensities from registered .tif files.
#
# Run by workflow/Snakefile: paths come from the command line, parameters from
# config/config.yaml (`slides`, `codex`). The zarr folder is the one step 1a
# wrote, so both steps always use the same zarr stores. For every cell, each
# CODEX channel is averaged over the cell's pixels in the label image.
#
# INPUT   --xenium-zarr-base  per-slide *.zarr from step 1a (cell labels + label images)
#         --codex-base        one folder per slide of *_<CH>_REGISTERED_*.tiff
# OUTPUT  --output-dir        one ID_<slide>_intensity.parquet per slide,
#                             key column "cell_uid" = "<cell_id>_<slide_ID>"

import argparse
import re
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd
import spatialdata as sd
import tifffile
import yaml

# elapsed time + peak memory, flushed so it shows up in SLURM logs immediately
_T0 = time.time()
def log(msg):
    gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    print(f"[{time.time() - _T0:7.1f}s | peak {gb:5.1f} GB] {msg}", flush=True)


# %% ############################## ARGUMENTS #################################

def parse_args():
    p = argparse.ArgumentParser(description="Step 1b: CODEX intensity retrieval")
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--xenium-zarr-base", required=True, help="dir of per-slide zarr from step 1a")
    p.add_argument("--codex-base", required=True, help="dir with one CODEX tif folder per slide")
    p.add_argument("--output-dir", required=True, help="output dir for intensity parquets")
    return p.parse_args()

args = parse_args()
cfg = yaml.safe_load(open(args.config_yaml))

XENIUM_ZARR_BASE = Path(args.xenium_zarr_base)
CODEX_BASE       = Path(args.codex_base)
OUTPUT_BASE      = Path(args.output_dir)

SLIDES      = cfg.get("slides") or None   # None = every slide with both a zarr and a CODEX folder
COMPARTMENT = cfg["codex"]["compartment"]  # "cell" or "nucleus"

# slide_ID = digits before the region part; matches with or without "__Region_N"
ZARR_SLIDE_RE  = re.compile(r"__(?P<slide_ID>\d+)(?=__|\.|$)")   # output-XETG00404__0056777__...
CODEX_SLIDE_RE = re.compile(r"^(?:ID_)?(?P<slide_ID>\d+)")       # ID_0056777__Region_1_scale0_tif

ID_COLS = ["cell_uid", "cell_id", "cell_labels"]

OUTPUT_BASE.mkdir(parents=True, exist_ok=True)


# %% ############################### HELPERS ##################################

def slide_of(name, pattern):
    m = pattern.search(str(name))
    return m.group("slide_ID") if m else None


def channel_of(path):
    """'..._<CHANNEL>_REGISTERED_...tiff' -> CHANNEL, else None."""
    parts = Path(path).name.split("_")
    if "REGISTERED" not in parts:
        return None                      # skips unregistered files
    i = parts.index("REGISTERED")
    return parts[i - 1] if i > 0 else None


def channel_files(folder):
    """{channel: path} for one folder."""
    out = {}
    for f in sorted(Path(folder).glob("*.tif*")):
        ch = channel_of(f)
        if ch is None:
            continue
        if ch in out:
            raise ValueError(f"Two files for channel {ch!r} in {folder}")
        out[ch] = f
    return out


def label_bincount(label_img, weights, n_labels, rows=2048):
    """
    Per-label sum over the image, accumulated in row chunks so no full-size
    int64/float64 copy of a multi-gigapixel array is ever made.
    weights=None -> pixel counts per label.
    """
    out = np.zeros(n_labels + 1, dtype=np.float64)
    for i in range(0, label_img.shape[0], rows):
        lab = label_img[i:i + rows].ravel()
        if weights is None:
            out += np.bincount(lab, minlength=n_labels + 1)
        else:
            w = weights[i:i + rows].ravel().astype(np.float64)
            out += np.bincount(lab, weights=w, minlength=n_labels + 1)
    return out


def load_labels(sdata, compartment):
    want = {"cell": "cell_labels", "nucleus": "nucleus_labels"}[compartment]
    if want not in sdata.labels:
        raise KeyError(f"No {want!r} label image; available: {list(sdata.labels)}")
    return sdata.labels[want].scale0["image"].compute().values


# %% ####################### PAIR ZARR <-> CODEX BY slide_ID ###################

zarr_by_slide = {}
for p in sorted(XENIUM_ZARR_BASE.glob("*.zarr")):
    sid = slide_of(p.name, ZARR_SLIDE_RE)
    if sid:
        zarr_by_slide.setdefault(sid, []).append(p)

codex_by_slide = {}
for p in sorted(x for x in CODEX_BASE.iterdir() if x.is_dir()):
    sid = slide_of(p.name, CODEX_SLIDE_RE)
    if sid is None:
        continue
    if sid in codex_by_slide:
        raise ValueError(f"Two CODEX folders for slide {sid}")
    codex_by_slide[sid] = p

slides = sorted(set(zarr_by_slide) & set(codex_by_slide))
if SLIDES:
    slides = [s for s in slides if s in SLIDES]
if not slides:
    raise ValueError("No slide has both a zarr and a CODEX folder.")

orphan = sorted(set(codex_by_slide) - set(zarr_by_slide))
if orphan:
    raise ValueError(f"CODEX folder with no Xenium zarr: {orphan}")
skipped = sorted(set(zarr_by_slide) - set(codex_by_slide))
if skipped:
    print(f"skipping (no CODEX folder): {skipped}")

for s in slides:                      # whole slides: exactly one zarr each
    if len(zarr_by_slide[s]) > 1:
        raise ValueError(f"Slide {s} matches several zarrs: "
                         f"{[p.name for p in zarr_by_slide[s]]}")

print(f"{len(slides)} slide(s):")
for s in slides:
    print(f"  {s}  {zarr_by_slide[s][0].name}  <->  {codex_by_slide[s].name}")


# %% ##################### DISCOVER THE CHANNEL UNION ##########################
# Filenames only, no images read.

files_by_slide = {s: channel_files(codex_by_slide[s]) for s in slides}
for s, fs in files_by_slide.items():
    if not fs:
        raise ValueError(f"[{s}] no *_REGISTERED_*.tif in {codex_by_slide[s]}")

CHANNELS = []
for s in slides:
    for ch in files_by_slide[s]:
        if ch not in CHANNELS:
            CHANNELS.append(ch)

presence = pd.DataFrame(
    {s: [ch in files_by_slide[s] for ch in CHANNELS] for s in slides}, index=CHANNELS
)
print(f"\n{len(CHANNELS)} channels in the union (0 = NaN column for that slide):")
print(presence.astype(int).to_string())
print("\nCHANNELS =", CHANNELS)


# %% ############################ MAIN PROCESSING ##############################

def process_slide(slide_id):
    sdata = sd.read_zarr(zarr_by_slide[slide_id][0])
    obs = sdata.tables["table"].obs
    if "cell_labels" not in obs:
        raise KeyError(f"[{slide_id}] no 'cell_labels' in table.obs; "
                       f"columns: {obs.columns.tolist()}")

    cell_id = (obs["cell_id"] if "cell_id" in obs else obs.index.to_series()).astype(str)
    labels_idx = obs["cell_labels"].to_numpy(dtype=np.int64)

    df = pd.DataFrame({
        "cell_uid":    cell_id.to_numpy() + f"_{slide_id}",
        "cell_id":     cell_id.to_numpy(),
        "cell_labels": labels_idx,
    })
    if df["cell_uid"].duplicated().any():
        raise ValueError(f"[{slide_id}] duplicated cell_uid")

    label_img = load_labels(sdata, COMPARTMENT)
    log(f"[{slide_id}] labels {label_img.shape} {label_img.dtype} "
        f"({label_img.nbytes / 1e9:.1f} GB)")

    # pixel count per label: the denominator for every channel mean
    n_labels = int(max(labels_idx.max(), label_img.max()))
    counts = label_bincount(label_img, None, n_labels)

    files = files_by_slide[slide_id]
    missing = [c for c in CHANNELS if c not in files]
    log(f"[{slide_id}] {len(files)}/{len(CHANNELS)} channels"
        + (f" | NaN: {missing}" if missing else ""))

    for n, ch in enumerate(CHANNELS, 1):      # same order for every slide
        if ch not in files:
            df[ch] = np.float32("nan")
            continue
        img = tifffile.imread(files[ch])
        if img.shape != label_img.shape:
            raise ValueError(f"[{slide_id}] shape mismatch for {ch}: "
                             f"{img.shape} vs {label_img.shape}")
        sums = label_bincount(label_img, img, n_labels)
        with np.errstate(invalid="ignore", divide="ignore"):
            means = sums / counts             # labels absent from the image -> NaN
        df[ch] = means[labels_idx].astype(np.float32)
        log(f"[{slide_id}] {n}/{len(CHANNELS)} {ch} ({img.dtype})")
        del img, sums, means

    del label_img
    return df[ID_COLS + CHANNELS]


results = {}
for i, s in enumerate(slides, 1):
    log(f"--- slide {i}/{len(slides)}: {s}")
    df = process_slide(s)
    path = OUTPUT_BASE / f"ID_{s}_intensity.parquet"
    df.to_parquet(path, index=False)
    results[s] = df
    log(f"[{s}] {len(df)} cells x {len(CHANNELS)} channels -> {path}")

log("all slides done")


# %% ################################ CHECKS ###################################

combined = pd.concat(results, names=["slide_ID"]).reset_index(level=0)
print("cells:", len(combined), "| unique cell_uid:", combined["cell_uid"].nunique())
assert combined["cell_uid"].is_unique, "cell_uid collision"

print("\ncells per slide:")
print(combined.groupby("slide_ID").size().to_string())

print("\nNaN fraction per channel:")
print(combined[CHANNELS].isna().mean().to_string())

print("\nmedian intensity per slide (look for scale differences between slides):")
print(combined.groupby("slide_ID")[CHANNELS].median().T.to_string())

cols = [list(pd.read_parquet(OUTPUT_BASE / f"ID_{s}_intensity.parquet").columns) for s in slides]
print("\nall parquets share the same column order:", all(c == cols[0] for c in cols))
log("done")
