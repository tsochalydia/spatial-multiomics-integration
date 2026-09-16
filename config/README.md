# Configuration

Two files control a run:

| File | Describes | Start from |
|------|-----------|------------|
| `config/config.yaml` | **the run**: input data, results folder, analysis parameters, job sizes | `config/config.example.yaml` |
| `profiles/slurm/config.yaml` | **the cluster**: SLURM partitions, account, retry behaviour | `profiles/slurm/config.yaml.example` |

Both are git-ignored, because they contain your paths.

`./run.sh` reads `config/config.yaml`. To use another run config, pass it explicitly:

```bash
./run.sh dry --configfile config/another_run.yaml
```

Only that file is read; nothing is inherited from `config/config.yaml`. So every config
file must be **complete**: copy the example and change what you need.

Paths may be absolute, or relative to the repository root.

After changing an analysis parameter, the affected step is **not** redone automatically.
Force it, and everything after it, e.g. `./run.sh submit --forcerun step2_data_qc`. See
"Re-running" in the main README.

---

## `slides`

```yaml
slides: ["0056764", "0056777"]
```

The slide IDs to process, as they appear in the Xenium run folder names
(`output-XETG00404__0056777__Region_1__…`) and CODEX folder names (`ID_0056777__…`).

- **Always put the IDs in quotes.** Without quotes, YAML reads `0056764` as an octal
  number and turns it into `24052`, and no slide is found.
- Steps 1a and 1b process only these slides. Step 2 reads every parquet that step 1b
  wrote, so use a fresh `results_dir` when you change this list.

---

## `paths`

| Key | Meaning |
|-----|---------|
| `xenium_raw` | Folder containing one raw Xenium output folder per slide (each with `experiment.xenium`). Exactly one run per slide. |
| `codex_base` | Folder containing one sub-folder per slide with the registered CODEX images `…_<CHANNEL>_REGISTERED_…tif(f)`, same size as the Xenium image. It must not contain folders for slides without Xenium data in this run: step 1b stops with an error. |
| `results_dir` | Where everything is written. The pipeline creates `step1a/`, `step1b/`, `step2/`, `step3/` and `logs/` inside. |

---

## `container`

| Key | Meaning |
|-----|---------|
| `sif` | Path to `integration.sif`, the Apptainer image all steps run in |
| `extra_binds` | *(optional)* list of extra host folders to make visible inside the container |

The container automatically sees the two input folders, `results_dir`, the scripts and
the folder of the config file. Add `extra_binds` only if the data contains symlinks
pointing somewhere else.

---

## `xenium`: step 1a

| Key | Default | Meaning |
|-----|---------|---------|
| `read_transcripts` | `false` | Also load the transcript points. Not needed for the cell × gene table; slow and memory-hungry. |
| `write_zarr` | `true` | Write one SpatialData `.zarr` per slide. **Must stay `true`**: step 1b reads the cell label images from these stores. |
| `overwrite_zarr` | `false` | `false`: reuse an existing zarr if it can be read (saves ~10 min per slide). `true`: always rebuild it. |

---

## `codex`: step 1b

| Key | Default | Meaning |
|-----|---------|---------|
| `compartment` | `cell` | Which segmentation to average CODEX pixels over: `cell` or `nucleus`. |

---

## `qc`: step 2

**Xenium (RNA) cell filter**

| Key | Default | Meaning |
|-----|---------|---------|
| `min_counts` | `5` | Remove cells with fewer total transcripts |
| `min_genes` | `2` | Remove cells expressing fewer genes |

**Which CODEX channels become protein markers**

| Key | Default | Meaning |
|-----|---------|---------|
| `exclude_channels` | `["DAPI"]` | Channels that are never markers |
| `markers` | `null` | `null`: every channel except `exclude_channels`. Or an explicit list, e.g. `["CD4", "CD8", …]`, to fix the panel. Names must match the channel names in the image file names. |

**Per-slide CODEX clean-up** (thresholds are computed separately for each slide)

| Key | Default | Meaning |
|-----|---------|---------|
| `background_marker` | `CD4` | A marker whose background is subtracted per slide; it is renamed `<marker>_cor`. Must be one of the markers. `null` switches this off. |
| `background_percentile` | `10` | The background subtracted is this percentile of that marker; negative values are set to 0 |
| `lower_q` | `0.001` | Per marker, values below this quantile are set to 0 |
| `upper_q` | `0.999` | Quantile used to flag extreme values |
| `max_markers_above_upper_q` | `4` | Remove a cell if it is above `upper_q` in **this many markers or more**. Must be at least 1; `0` removes every cell. |

**Naming**

| Key | Default | Meaning |
|-----|---------|---------|
| `marker_suffix` | `"_CDX"` | Added to protein names in the merged AnnData (e.g. `CD8_CDX`); step 3 strips it again for plot labels |

`step2/tracking_counts.csv` shows how many cells each filter removed.

---

## `integration`: step 3 (GPU)

**Model and clustering**

| Key | Default | Meaning |
|-----|---------|---------|
| `batch_key` | `slide_ID` | Column in the cell table used as the totalVI batch (here: one batch per slide) |
| `resolution` | `1` | Leiden resolution: higher gives more, smaller clusters |
| `target_sum` | `100` | RNA normalisation for visualisation only. totalVI itself trains on raw counts. |
| `max_epochs` | `400` | Upper limit on training epochs |
| `train_size` | `0.9` | Fraction of cells used for training; the rest is used for validation |
| `early_stopping_patience` | `20` | Stop when the validation loss has not improved for this many epochs |
| `dedup` | `raise` | If two cells end up with the same ID: `raise` stops with a report; `make_unique` renumbers them and continues |

**Marker tables and plots**

| Key | Default | Meaning |
|-----|---------|---------|
| `top_rna` | `10` | Max. RNA markers listed per cluster in the DE table |
| `top_pro` | `5` | Max. protein markers listed per cluster in the DE table |
| `bf_rna` | `1.2` | Minimum `bayes_factor` for an RNA marker. In scvi-tools' differential expression, `bayes_factor` is ln(Bayes factor), the natural logarithm. |
| `bf_protein` | `0.45` | Same for proteins; lower because protein signal is noisier |
| `min_nonzero_prop` | `0.1` | Intended: a gene must be detected in more than this fraction of a cluster's cells to appear in the dotplot. **Currently has no effect**: the plotting code uses a fixed `0.1`. |

`top_rna`, `top_pro`, `bf_rna` and `bf_protein` are part of the output file names, e.g.
`DE_table_top_10_rna_top_5_proteins_lnBF_rna_1.2_protein_0.45.csv`. Changing them
produces new files instead of overwriting the old ones.

> **Note:** after changing `top_rna`, `top_pro` or `bf_rna`, step 3 is **not** redone by
> itself: the workflow only checks for the matrixplot, whose name contains only
> `bf_protein`. Force it with `./run.sh submit --forcerun step3_data_integration`.

---

## `resources`: job sizes

```yaml
resources:
  step1a: { mem_mb_per_cpu: 16000, runtime: 240, threads: 4 }
```

| Key | Meaning |
|-----|---------|
| `mem_mb_per_cpu` | Memory per CPU in MB. Total memory = `mem_mb_per_cpu × threads`. (ETH Euler requires memory per CPU.) |
| `runtime` | Time limit in minutes |
| `threads` | CPUs for the job |

The example values fit the full two-slide dataset (~1.5 million cells). Peak memory in
that run: step 1a ~60 GB, step 1b ~43 GB, step 2 ~5 GB, step 3 ~123 GB. Step 3 took
about 9.5 h. Smaller datasets need much less.

---

## Cluster settings: `profiles/slurm/config.yaml`

| Key | Default in the example | Meaning |
|-----|------------------------|---------|
| `executor` | `slurm` | Submit each step as its own SLURM job |
| `jobs` | `8` | Max. number of jobs at once |
| `default-resources.slurm_partition` | `normal.24h` | Partition for the CPU steps (1a, 1b, 2) |
| `default-resources.slurm_account` | *(commented out)* | Set it if your cluster requires `-A <account>` |
| `default-resources.runtime`, `mem_mb_per_cpu` | `120`, `4000` | Fallbacks for anything without its own value |
| `set-resources.step3_data_integration.slurm_partition` | `gpu.24h` | GPU partition for step 3. Use a CUDA 13 partition (Euler: `cuda13pr.24h`) if torch reports a driver/CUDA error. |
| `keep-going` | `True` | Keep running independent steps if one fails |
| `rerun-incomplete` | `True` | Redo a step whose outputs were left incomplete |
| `restart-times` | `1` | Retry a failed job once |
| `latency-wait` | `60` | Seconds to wait for files to appear on the shared filesystem |

Partition names are specific to ETH Euler; change them for another cluster.
