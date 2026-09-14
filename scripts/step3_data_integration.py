# %% STEP 3 — totalVI integration of the Xenium RNA + CODEX protein modalities.
#
# Snakemake-adapted version of 3_data_integration.py: paths from the command line,
# training/clustering parameters from the config YAML. Runs on GPU (scvi_env).
#
# INPUT   --input-h5ad  final_adata.h5ad from step 2
# OUTPUT  --save-dir     trained totalVI model, joint latent, denoised RNA/protein,
#                        Leiden clusters, UMAPs, DE table + plots

import argparse
import resource
import time
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import mudata as md
import muon
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import scvi
import seaborn as sns
import torch
import yaml

mpl.rcParams["savefig.dpi"] = 300
scvi.settings.seed = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")

print("Last run with scvi-tools version:", scvi.__version__)

# elapsed time + peak memory, flushed so it shows up in SLURM logs immediately
_T0 = time.time()
def log(msg):
    gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
    print(f"[{time.time() - _T0:7.1f}s | peak {gb:5.1f} GB] {msg}", flush=True)


# %% ############################## ARGUMENTS #################################

def parse_args():
    p = argparse.ArgumentParser(description="Step 3: totalVI integration")
    p.add_argument("--config-yaml", required=True)
    p.add_argument("--input-h5ad", required=True, help="final_adata.h5ad from step 2")
    p.add_argument("--save-dir", required=True, help="output dir")
    return p.parse_args()

args = parse_args()
cfg  = yaml.safe_load(open(args.config_yaml))
ic   = cfg["integration"]

INPUT_H5AD = Path(args.input_h5ad)
SAVE_DIR   = Path(args.save_dir)
SAVE_DIR.mkdir(parents=True, exist_ok=True)

BATCH_KEY = ic["batch_key"]
RNA_LAYER = "counts"

# Duplicate suffixed cell ids: "raise" = stop with a diagnostic; "make_unique" = renumber.
DEDUP = ic["dedup"]

TOTALVI_LATENT_KEY   = "X_totalVI"
TOTALVI_CLUSTERS_KEY = "leiden_totalVI"
RESOLUTION = ic["resolution"]

# RNA normalisation target (visualisation only — totalVI trains on raw counts)
TARGET_SUM = ic["target_sum"]

# Training
MAX_EPOCHS              = ic["max_epochs"]
TRAIN_SIZE              = ic["train_size"]
EARLY_STOPPING_PATIENCE = ic["early_stopping_patience"]

# DE analysis: max features exported per cluster and log2(BF) thresholds
TOP_RNA = ic["top_rna"]
TOP_PRO = ic["top_pro"]
BF_rna     = ic["bf_rna"]
BF_protein = ic["bf_protein"]

# Genes must also be detected in more than this fraction of a cluster's cells (dotplot)
MIN_NONZERO_PROP = ic["min_nonzero_prop"]

MARKER_SUFFIX = cfg["qc"]["marker_suffix"]   # set by step 2, stripped for readable plot labels


# %% ########################### LOAD AND BUILD MUDATA #########################

log("loading adata")
adata = ad.read_h5ad(INPUT_H5AD)

def check_unique_index(adata, batch_key, dedup="raise"):
    dup = adata.obs_names.duplicated(keep=False)
    if not dup.any():
        print(f"[index] all {adata.n_obs} suffixed cell ids are unique")
        return adata

    lines = [f"{int(adata.obs_names.duplicated().sum())} duplicated suffixed cell ids "
             f"({int(dup.sum())} rows involved) in final_adata.h5ad"]
    if batch_key in adata.obs:
        per_id = adata.obs.loc[dup].groupby(adata.obs_names[dup])[batch_key].nunique()
        lines.append(f"  within one slide : {int((per_id == 1).sum())}  "
                     f"(suffix can't fix — check base cell_id uniqueness in 1a/1b)")
        lines.append(f"  spanning >1 slide: {int((per_id > 1).sum())}  "
                     f"(same slide_ID assigned to two runs?)")
    lines.append("  examples: " + ", ".join(map(str, adata.obs_names[dup][:5])))
    report = "\n".join(lines)

    if dedup == "raise":
        raise ValueError(report + "\n  -> set integration.dedup='make_unique' to renumber and continue.")
    print("[index][WARNING]\n" + report)
    adata.obs["obs_name_predup"] = adata.obs_names.astype(str)
    adata.obs_names = ad.utils.make_index_unique(adata.obs_names)
    print(f"[index] renumbered -> {int(adata.obs_names.duplicated().sum())} duplicates remain")
    return adata

adata = check_unique_index(adata, BATCH_KEY, dedup=DEDUP)
if RNA_LAYER not in adata.layers:
    raise KeyError(f"No '{RNA_LAYER}' layer; totalVI must train on raw counts.")
if BATCH_KEY not in adata.obs:
    raise KeyError(f"No '{BATCH_KEY}' in .obs; columns: {adata.obs.columns.tolist()}")
if "protein_markers" not in adata.uns:
    raise KeyError("No .uns['protein_markers']; run step 2 first.")

adata.obs[BATCH_KEY] = adata.obs[BATCH_KEY].astype(str).astype("category")
protein_markers = list(adata.uns["protein_markers"])

log(f"loaded: {adata.n_obs} cells x {adata.n_vars} genes x {len(protein_markers)} proteins")
print("batches:", adata.obs[BATCH_KEY].value_counts().sort_index().to_dict())
print("markers:", protein_markers)

# Normalise + log-transform the RNA modality for visualisation.
sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
sc.pp.log1p(adata)
log("normalized and log-transformed")

# totalVI expects CITE-seq-like ADT counts, so round the CODEX intensities.
X_codex = np.asarray(adata.obsm["protein_expression"], dtype=np.float64)
n_nan = int(np.isnan(X_codex).sum())
X_codex = np.round(np.nan_to_num(X_codex, nan=0.0)).astype(int)
print(f"protein matrix: {X_codex.shape} | NaN set to 0: {n_nan} "
      f"({100 * n_nan / X_codex.size:.2f}% of entries)")

protein_adata = ad.AnnData(
    X=X_codex,
    obs=adata.obs.copy(),
    var=pd.DataFrame(index=protein_markers),
)
protein_adata.obs_names = adata.obs_names
del adata.obsm["protein_expression"]

mdata = md.MuData({"rna": adata, "protein": protein_adata})

# totalVI needs a dense RNA matrix
if sp.issparse(mdata["rna"].X):
    mdata["rna"].X = mdata["rna"].X.toarray()

mdata.write(SAVE_DIR / f"mu_norm{TARGET_SUM}_log1p.h5mu")
log("mdata saved")


# %% ################################ TRAINING #################################

log("setting up totalVI")
scvi.model.TOTALVI.setup_mudata(
    mdata,
    rna_layer=RNA_LAYER,
    protein_layer=None,
    batch_key=BATCH_KEY,
    modalities={
        "rna_layer": "rna",
        "protein_layer": "protein",
        "batch_key": "rna",
    },
)

model = scvi.model.TOTALVI(mdata)
model.train(
    train_size=TRAIN_SIZE,
    max_epochs=MAX_EPOCHS,
    accelerator="gpu",
    devices=1,
    early_stopping=True,
    early_stopping_patience=EARLY_STOPPING_PATIENCE,
    early_stopping_monitor="elbo_validation",
)
log("model trained")
model.save(SAVE_DIR, overwrite=True)
log("model saved")

# Training curve — quick check that the model actually converged
fig, ax = plt.subplots(1, 1)
model.history["elbo_train"].plot(ax=ax, label="train")
model.history["elbo_validation"].plot(ax=ax, label="validation")
ax.set(title="Negative ELBO over training epochs")
ax.legend()
plt.savefig(SAVE_DIR / "elbo_training.png", bbox_inches="tight", dpi=300)
plt.close()

rna     = mdata.mod["rna"]
protein = mdata.mod["protein"]

# Latent representation is stored on the rna modality by convention
rna.obsm[TOTALVI_LATENT_KEY] = model.get_latent_representation()
mdata.write(SAVE_DIR / "mdata_latent.h5mu")
log("latent representation saved")

rna_denoised, protein_denoised = model.get_normalized_expression(n_samples=25, return_mean=True)
rna.layers["denoised_rna"]         = rna_denoised
protein.layers["denoised_protein"] = protein_denoised
protein.layers["protein_foreground_prob"] = 100 * model.get_protein_foreground_probability(
    n_samples=25, return_mean=True
)

# Strip the "_CDX" suffix for readable plot labels (keeps e.g. "CD4_cor" intact)
protein.var["clean_names"] = [
    p[: -len(MARKER_SUFFIX)] if p.endswith(MARKER_SUFFIX) else p for p in protein.var_names
]

mdata.update()
mdata.write(SAVE_DIR / "mdata_denoised_forprob.h5mu")
log("denoised mdata saved")


# %% ############################### CLUSTERING ################################

sc.pp.neighbors(rna, use_rep=TOTALVI_LATENT_KEY)
sc.tl.umap(rna)
sc.tl.leiden(rna, key_added=TOTALVI_CLUSTERS_KEY, resolution=RESOLUTION,
             flavor="igraph", n_iterations=2)
log(f"clustering finished: {rna.obs[TOTALVI_CLUSTERS_KEY].nunique()} clusters")

mdata.update()
mdata.write(SAVE_DIR / "mdata_leiden.h5mu")

# Mirror clusters + latent onto the protein modality so scanpy can build its dendrogram there too
sc.tl.dendrogram(rna, groupby=TOTALVI_CLUSTERS_KEY, use_rep=TOTALVI_LATENT_KEY)
protein.obs[TOTALVI_CLUSTERS_KEY] = rna.obs[TOTALVI_CLUSTERS_KEY]
protein.obsm[TOTALVI_LATENT_KEY]  = rna.obsm[TOTALVI_LATENT_KEY]
sc.tl.dendrogram(protein, groupby=TOTALVI_CLUSTERS_KEY, use_rep=TOTALVI_LATENT_KEY)

mdata.update()
mdata.write(SAVE_DIR / "mdata_leiden_dendrogram.h5mu")
log("dendrogram done")


# %% ################################## UMAPS ##################################

n_clusters = rna.obs[TOTALVI_CLUSTERS_KEY].nunique()
colors  = list(plt.cm.tab20.colors) + list(plt.cm.tab20b.colors) + list(plt.cm.tab20c.colors)
palette = colors[:n_clusters]

muon.pl.embedding(
    mdata, basis="rna:X_umap", color=[f"rna:{TOTALVI_CLUSTERS_KEY}"],
    frameon=False, ncols=1, palette=palette, show=False,
)
plt.savefig(SAVE_DIR / "umap_clusters.png", bbox_inches="tight", dpi=300)
plt.close()

muon.pl.embedding(
    mdata, basis="rna:X_umap", color=[f"rna:{TOTALVI_CLUSTERS_KEY}"],
    frameon=False, ncols=1, legend_loc="on data", palette=palette, show=False,
)
plt.savefig(SAVE_DIR / "umap_clusters_on_data.png", bbox_inches="tight", dpi=300)
plt.close()

muon.pl.embedding(
    mdata, basis="rna:X_umap", layer="protein_foreground_prob",
    color=protein.var_names, frameon=False, ncols=6,
    vmax="p99", wspace=0.1, color_map="cividis", show=False,
)
plt.savefig(SAVE_DIR / "forprob.png", bbox_inches="tight", dpi=300)
plt.close()

muon.pl.embedding(
    mdata, basis="rna:X_umap", layer="denoised_protein",
    color=protein.var_names, frameon=False, ncols=6, vmax="p99", wspace=0.1, show=False,
)
plt.savefig(SAVE_DIR / "denoised_protein.png", bbox_inches="tight", dpi=300)
plt.close()

log("all UMAPs saved")


# %% ############################### DE ANALYSIS ###############################

log("starting DE analysis")
de_df = model.differential_expression(
    groupby=f"rna:{TOTALVI_CLUSTERS_KEY}", delta=0.5, batch_correction=True
)
de_df.to_csv(SAVE_DIR / "differential_expression.csv")

# scvi keeps the protein var_names in the DE index, appending "_protein" only on collisions.
protein_names = set(protein.var_names)

def split_rna_protein(df):
    is_protein = df.index.isin(protein_names) | df.index.str.endswith("_protein")
    return df[~is_protein], df[is_protein]


filtered_table = []
filtered_rna, filtered_pro = {}, {}
cats = rna.obs[TOTALVI_CLUSTERS_KEY].cat.categories

for c in cats:
    cell_type_df = de_df.loc[de_df.comparison == f"{c} vs Rest"].copy()
    cell_type_df = cell_type_df[cell_type_df["lfc_median"] > 0]

    data_rna, data_pro = split_rna_protein(cell_type_df)

    data_rna = data_rna[data_rna["bayes_factor"] > BF_rna]
    data_pro = data_pro[data_pro["bayes_factor"] > BF_protein]   # lower: proteins are noisier

    data_rna = data_rna.sort_values("lfc_median", ascending=False)
    data_pro = data_pro.sort_values("lfc_median", ascending=False)

    top_rna = data_rna.head(TOP_RNA)
    top_pro = data_pro.head(TOP_PRO)

    filtered_table.append({
        "cluster": c,
        "n_cells": int((rna.obs[TOTALVI_CLUSTERS_KEY] == c).sum()),
        "rna_markers": ", ".join(top_rna.index),
        "protein_markers": ", ".join(top_pro.index),
        "rna_bayes_factor": ", ".join(top_rna["bayes_factor"].round(2).astype(str)),
        "protein_bayes_factor": ", ".join(top_pro["bayes_factor"].round(2).astype(str)),
        "n_rna_after_filter": len(data_rna),
        "n_pro_after_filter": len(data_pro),
    })

    filtered_rna[c] = data_rna[
        data_rna["non_zeros_proportion1"] > MIN_NONZERO_PROP
    ].index.tolist()[:2]
    filtered_pro[c] = data_pro.index.tolist()[:3]

marker_table = pd.DataFrame(filtered_table)
marker_table.to_csv(
    SAVE_DIR / f"DE_table_top_{TOP_RNA}_rna_top_{TOP_PRO}_proteins_"
               f"log2BF_rna_{BF_rna}_protein_{BF_protein}.csv",
    index=False,
)
print(marker_table[["cluster", "n_cells", "n_rna_after_filter", "n_pro_after_filter"]].to_string())
log("marker table saved")


# %% ########################### DOTPLOT / MATRIXPLOT ##########################

## Dotplot creation
filtered_pro = {}
filtered_rna = {}
cats = rna.obs[TOTALVI_CLUSTERS_KEY].cat.categories
for c in cats:
    cid = f"{c} vs Rest"
    cell_type_df = de_df.loc[de_df.comparison == cid]
    cell_type_df = cell_type_df.sort_values("lfc_median", ascending=False)

    cell_type_df = cell_type_df[cell_type_df.lfc_median > 0]

    pro_rows = cell_type_df.index.str.contains("protein")
    data_pro = cell_type_df.iloc[pro_rows]
    data_pro = data_pro[data_pro["bayes_factor"] > BF_protein]

    data_rna = cell_type_df.iloc[~pro_rows]
    data_rna = data_rna[data_rna["bayes_factor"] > BF_rna]
    data_rna = data_rna[data_rna["non_zeros_proportion1"] > 0.1]

    filtered_pro[c] = data_pro.index.tolist()[:3]
    filtered_rna[c] = data_rna.index.tolist()[:2]

sc.pl.dotplot(
    rna,
    filtered_rna,
    groupby=TOTALVI_CLUSTERS_KEY,
    dendrogram=True,
    standard_scale="var",
    swap_axes=True,
)

plt.savefig(SAVE_DIR / f'dotplot_top_{TOP_RNA}_rna_top_{TOP_PRO}_proteins_log2BF_rna_{BF_rna}_protein_{BF_protein}.png', bbox_inches="tight", dpi=300)
plt.close()

sc.pl.matrixplot(
    protein, protein.var["clean_names"].tolist(),
    groupby=TOTALVI_CLUSTERS_KEY, gene_symbols="clean_names", dendrogram=True,
    layer="denoised_protein", cmap="Greens", standard_scale="var",
    swap_axes=True, show=False,
)
plt.savefig(SAVE_DIR / f"matrixplot_protein_log2BF_{BF_protein}.png",
            bbox_inches="tight", dpi=300)
plt.close()

log("plots saved — done")
