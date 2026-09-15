# %% STEP 1a — Concatenate the Xenium cell x gene tables of all TMA slides.
#
# Snakemake-adapted version of 1a_xenium_concatenation_test.py: paths come from
# the command line, tunable parameters from the workflow config YAML. The core
# logic (per-slide zarr build, "<cell_id>_<slide_ID>" indexing, concat) is
# unchanged.
#
# INPUT   a directory of raw Xenium runs (each contains experiment.xenium)
# OUTPUT  --zarr-dir     one <run_name>.zarr per slide (read by step 1b)
#         --output-h5ad  concatenated AnnData, raw counts in .X and .layers["counts"]

import argparse
import re
import resource
import shutil
import time
from collections import defaultdict
from pathlib import Path

import anndata as ad
import pandas as pd
import spatialdata as sd
import yaml
from spatialdata_io import xenium

print("anndata", ad.__version__)
print("spatialdata", sd.__version__)

# elapsed time + peak memory, flushed so it shows up in SLURM logs immediately
_T0 = time.time()
def log(msg):
    gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    print(f"[{time.time() - _T0:7.1f}s | peak {gb:5.1f} GB] {msg}", flush=True)


# %% ############################## ARGUMENTS #################################

def parse_args():
    p = argparse.ArgumentParser(description="Step 1a: Xenium concatenation")
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--xenium-raw", required=True, help="dir of raw Xenium runs")
    p.add_argument("--zarr-dir", required=True, help="output dir for per-slide zarr stores")
    p.add_argument("--output-h5ad", required=True, help="output concatenated .h5ad")
    return p.parse_args()

args = parse_args()
cfg = yaml.safe_load(open(args.config_yaml))

DATA_DIR    = Path(args.xenium_raw)
ZARR_DIR    = Path(args.zarr_dir)
OUTPUT_H5AD = Path(args.output_h5ad)

# None = every slide found under DATA_DIR
ONLY_SLIDES      = cfg.get("slides") or None
MAX_FOLDERS      = None
READ_TRANSCRIPTS = cfg["xenium"]["read_transcripts"]
WRITE_ZARR       = cfg["xenium"]["write_zarr"]
OVERWRITE_ZARR   = cfg["xenium"]["overwrite_zarr"]

XENIUM_MARKER     = "experiment.xenium"
DEFAULT_TABLE_KEY = "table"
xen_kwargs = {} if READ_TRANSCRIPTS else {"transcripts": False}

# slide_ID = digits before the region part; identical pattern to the CODEX step,
# so the two modalities can never derive a different key from the same run name.
SLIDE_RE  = re.compile(r"__(?P<slide_id>\d+)(?=__|$)")
REGION_RE = re.compile(r"__(?P<region>Region_\d+)")

assert DATA_DIR.is_dir(), f"DATA_DIR does not exist: {DATA_DIR}"
ZARR_DIR.mkdir(parents=True, exist_ok=True)
OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)


# %% ########################### DISCOVER RUNS ################################

def find_xenium_dirs(parent):
    parent = Path(parent)
    return sorted(p.parent for p in parent.rglob(XENIUM_MARKER))

folders = find_xenium_dirs(DATA_DIR)
if MAX_FOLDERS:
    folders = folders[:MAX_FOLDERS]
assert folders, f"No folders containing {XENIUM_MARKER} under {DATA_DIR}"
print(f"Using {len(folders)} run(s):")
for f in folders:
    print("  ", f.name)


def parse_slide_metadata(run_name):
    m = SLIDE_RE.search(run_name)
    if m is None:
        raise ValueError(f"Cannot parse slide_ID from {run_name!r}")
    r = REGION_RE.search(run_name)
    return {"run_name": run_name,
            "slide_id": m.group("slide_id"),
            "region": r.group("region") if r else None}

metas = {f: parse_slide_metadata(f.name) for f in folders}

if ONLY_SLIDES:
    folders = [f for f in folders if metas[f]["slide_id"] in ONLY_SLIDES]
    metas = {f: metas[f] for f in folders}
    found = {metas[f]["slide_id"] for f in folders}
    assert folders, f"None of ONLY_SLIDES present in {DATA_DIR}: {ONLY_SLIDES}"
    missing = sorted(set(ONLY_SLIDES) - found)
    assert not missing, f"Requested slide(s) not found in {DATA_DIR}: {missing}"
    print(f"filtered to {len(folders)} run(s) for slides {sorted(found)}")

by_slide = defaultdict(list)
for f, m in metas.items():
    by_slide[m["slide_id"]].append(f.name)

summary = pd.DataFrame(list(metas.values()))
summary["n_folders_for_slide"] = summary["slide_id"].map(lambda s: len(by_slide[s]))
print(summary)
print("\nslides:", {s: len(v) for s, v in by_slide.items()})

# Whole slides only: one run per slide, so <cell_id>_<slide_ID> is unique.
multi = {s: v for s, v in by_slide.items() if len(v) > 1}
if multi:
    raise ValueError(f"These slides have more than one run: {multi}. "
                     "This pipeline expects one run per slide.")


# %% ############################### HELPERS ##################################

def extract_table(sdata, table_key=None):
    tables = getattr(sdata, "tables", None)
    if tables is not None and len(tables):
        if table_key and table_key in tables:
            return tables[table_key].copy()
        if DEFAULT_TABLE_KEY in tables:
            return tables[DEFAULT_TABLE_KEY].copy()
        return next(iter(tables.values())).copy()
    return sdata.table.copy()


# sanity read one folder
sdata0 = xenium(str(folders[0]), **xen_kwargs)
print(sdata0)
print(extract_table(sdata0))
del sdata0


def load_sdata(folder, zarr_path):
    # reuse an existing store only if not overwriting AND it reads cleanly
    if WRITE_ZARR and zarr_path.exists() and not OVERWRITE_ZARR:
        try:
            print(f"[{folder.name}] read existing zarr")
            return sd.read_zarr(zarr_path)
        except Exception as e:
            print(f"[{folder.name}] existing zarr unreadable ({type(e).__name__}: {e}) -> rebuilding")
            shutil.rmtree(zarr_path)

    log(f"[{folder.name}] read raw")
    sdata = xenium(str(folder), **xen_kwargs)
    log(f"[{folder.name}] raw read done")
    if WRITE_ZARR:
        ZARR_DIR.mkdir(parents=True, exist_ok=True)
        if zarr_path.exists():
            shutil.rmtree(zarr_path)
        log(f"[{folder.name}] write zarr {zarr_path}")
        try:
            sdata.write(zarr_path, overwrite=True)
        except TypeError:
            sdata.write(zarr_path)
        log(f"[{folder.name}] zarr written")
    return sdata


def process_run(folder):
    meta = metas[folder]
    slide_id = meta["slide_id"]
    zarr_path = ZARR_DIR / f"{folder.name}.zarr"

    sdata = load_sdata(folder, zarr_path)

    adata = extract_table(sdata)
    adata.layers["counts"] = adata.X.copy()
    adata.obs["cell_id_original"] = adata.obs_names.astype(str)
    adata.obs_names = [f"{cid}_{slide_id}" for cid in adata.obs_names]
    adata.obs["slide_ID"] = slide_id
    if meta["region"] is not None:
        adata.obs["region"] = meta["region"]
    adata.obs["xenium_run"] = folder.name
    adata.uns.pop("spatialdata_attrs", None)
    log(f"[{folder.name}] table ready: {adata.n_obs} cells x {adata.n_vars} genes (suffix _{slide_id})")
    return adata


# %% ############################# PROCESS + CONCAT ###########################

adatas = []
for _i, _f in enumerate(folders, 1):
    log(f"--- run {_i}/{len(folders)}: {_f.name}")
    adatas.append(process_run(_f))
log(f"all runs done: {len(adatas)} tables")

assert all("counts" in a.layers for a in adatas), "'counts' missing in some table"
log("concatenating")
combined = ad.concat(
    adatas, axis=0, join="outer", merge="same",
    uns_merge=None, index_unique=None, fill_value=0,
)
log(f"concatenated: {combined.n_obs} cells x {combined.n_vars} genes")

# No renaming: every index entry stays plain <cell_id>_<slide_ID>, or we stop.
n_dup = int(combined.obs_names.duplicated().sum())
if n_dup:
    raise ValueError(f"{n_dup} duplicated cell ids with the _slide_ID suffix.")
print("All cell ids unique with _slide_ID suffix alone.")

for col in ["slide_ID", "region", "xenium_run"]:
    if col in combined.obs:
        combined.obs[col] = combined.obs[col].astype("category")
print(combined)


# %% ################################ CHECKS ##################################

print("cells:", combined.n_obs, "| genes:", combined.n_vars)
print("\ncells per slide_ID:")
print(combined.obs["slide_ID"].value_counts().sort_index())
if "region" in combined.obs:
    print("\ncells per (slide_ID, region):")
    print(combined.obs.groupby(["slide_ID", "region"], observed=True).size())
print("\ncounts layer present:", "counts" in combined.layers)
diff = combined.X - combined.layers["counts"]
print("max |X - counts| (should be 0):", abs(diff).max())

expected = combined.obs["cell_id_original"].astype(str) + "_" + combined.obs["slide_ID"].astype(str)
print("index is plain <cell_id>_<slide_ID>:", bool((combined.obs_names == expected.to_numpy()).all()))
print(combined.obs.head())


# %% ################################# SAVE ###################################

log(f"writing {OUTPUT_H5AD}")
combined.write_h5ad(OUTPUT_H5AD, compression="gzip")
log(f"saved -> {OUTPUT_H5AD} ({OUTPUT_H5AD.stat().st_size / 1e9:.2f} GB)")

# verify a fresh zarr store round-trips (matters for step 1b)
test_store = ZARR_DIR / f"{folders[0].name}.zarr"
sd.read_zarr(test_store)
print("fresh store reads back OK:", test_store)

check = ad.read_h5ad(OUTPUT_H5AD)
print(check)
print("counts layer after reload:", "counts" in check.layers)
log("done")
