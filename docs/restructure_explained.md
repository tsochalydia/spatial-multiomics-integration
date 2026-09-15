# The project tidy-up, explained

What changed in this project folder, why, and how we know the results are unaffected.

## In short

- **The analysis did not change.** The four step scripts were moved to a new folder, but
  not a single line in them was edited. The real run's results (in `test/`, despite the
  name) and the old logs in `logs/` were not touched.
- **The packaging changed.** No personal file paths inside the files, all settings in one
  folder, labelled folders, one launch script. Before, the project worked but was hard for
  anyone else to use.
- **Checked in two ways, both passed:** the planned work was compared with the original
  plan after every single change, and the whole workflow ran start to finish on a tiny
  practice dataset.
- **Still to do:** run the original and the tidied workflow on the same data and confirm
  they give the same results. Then put all changes up for review on GitHub (a *pull
  request*) before they become the official version.

## Words used below

| Word | Meaning |
|------|---------|
| **Steps 1a–3** | **1a** collect the gene counts of all cells from all slides into one table · **1b** measure each protein's brightness per cell in the CODEX images · **2** remove bad-quality cells, join genes and proteins · **3** train the model, cluster the cells, draw the plots (needs a graphics card, **GPU**) |
| **Snakemake** | The program that runs the steps in order and skips steps that are already finished. The steps and their order are written in the `Snakefile`; together with the step scripts this is the **workflow**. |
| **Cluster** | The university's shared computers. Jobs wait in a queue (the queue system is **SLURM**) until a machine is free. |
| **Container** | `integration.sif`: one file with all the software inside, so every computer runs exactly the same versions. It can only see the cluster folders it is given access to. |
| **Dry run** | Snakemake prints exactly what it *would* do (steps, files, commands) without doing anything. Takes seconds. |
| **git, commit** | git records the project's history. Each saved change is a **commit** and can be inspected or undone on its own. |

## What we changed, and why

Three phases, one small change at a time, each saved as its own commit.

### Phase 1: remove personal paths

*Problem:* many files contained exact paths on Lydia's account
(`/cluster/project/moor/lydia/...`), so anyone else would have had to find and change them
one by one.

| # | Change | Why |
|---|--------|-----|
| 1.1 | Snakemake reads the settings file you actually give it. | A practice run uses its own settings, not the real run's. |
| 1.2 | One setting, `results_dir`, replaces five separate output paths. The step logs now go inside it, into `logs/` (for the real run: `test/logs/`). The old logs stay in `logs/`. Sub-folder names (`step1a`, `step1b`, …) are unchanged. | Fewer things to fill in. Logs used to go to the same folder whatever the data, so a practice run would have overwritten the real run's logs. The unchanged names keep the finished results recognised. |
| 1.3 | Short paths like `results` are allowed, meaning "inside the project folder". | Works for anyone, whatever their account. |
| 1.4 | The container settings shrank to one line: Snakemake works out by itself which folders the container may see. | Less to fill in, fewer mistakes. |
| 1.5 | Only `run.sh` remains. It finds the software by itself, and its new `submit` option replaces `submit.sh`. | No personal paths. And the two scripts were almost identical, so every fix had to be made twice. |
| 1.6 | Personal paths removed from the software lists; a wrong name in one of them fixed. | Anyone can use the lists. |

### Phase 2: all settings in one place

*Problem:* settings were spread over `config.yaml`, a cluster settings file and the launch
scripts.

| # | Change | Why |
|---|--------|-----|
| 2.1 | All settings files moved into `config/`: real settings, practice settings, and an example to copy. | One place to look. |
| 2.2 | Queue names (e.g. `normal.24h`, `gpu.24h`) moved to the cluster settings in `profiles/slurm/`. Memory and time per step stay in `config/`. The number of GPUs is no longer a setting; the workflow fixes it at 1. | On another cluster, only the cluster file changes. Step 3 trains on exactly one GPU, so asking for more would leave the extra ones unused. |

### Phase 3: labelled folders

*Problem:* the workflow, software lists, container build instructions and settings all sat
together in one folder, so a newcomer couldn't tell what is what.

- **3.1** Files moved into the standard folder layout for Snakemake projects (table below).
- **3.2** The do-not-upload list (`.gitignore`) now also covers results, backup copies and a
  few temporary files, so personal settings and big files can't be uploaded by accident.

| Now | Contains | Before |
|-----|----------|--------|
| `README.md`, `run.sh` | How to start: the only script you run | `README.md`, `run.sh`, `submit.sh` |
| `config/` | **Your settings:** which data, where results go, how big the jobs are | `config.yaml`, `config.yaml.example` |
| `profiles/slurm/` | **Cluster settings:** which queues to use | `profiles/` |
| `workflow/` | `Snakefile` + the four step scripts | `Snakefile`, `scripts/` |
| `envs/` | Software lists with exact versions | `spatial_env.yml`, `scvi_env.yml`, `int_retrieval_env.yml`, `snakemake_env.yml` |
| `container/` | How `integration.sif` is built | `micromamba.dockerfile`, `entrypoint.sh` |
| `docs/` | Background, like this document | – |
| `test_data/` | Tiny practice dataset (kept on the cluster, not uploaded) | – |
| `integration.sif` (13 GB), `integration.tar` (36 GB), `test/`, `logs/` | Unchanged | Same place; not moved, not uploaded |

**Nothing was lost.** The only deleted file is `submit.sh`, and it is still in the git
history. Before the bigger edits of the personal settings files, a backup copy was saved
next to them (the files ending in `.backup…`). A bookmark in the history, the tag
`pre-reorg`, marks the moment before the tidy-up, so that state can be restored at any time.

## How we know nothing broke

### Check 1: the plan, after every change

Before any change, we saved the dry-run plan as the reference. After every change, we made
a new plan and compared it with the reference line by line. The only differences were the
intended ones, new file locations. Inputs, outputs, software, memory, time and every
setting passed to the scripts stayed the same.

For the real data, the dry run also kept saying **"Nothing to be done"**: Snakemake still
sees the finished results as complete and won't redo hours of work.

### Check 2: a full run on tiny practice data

The full data needs more than 10 hours of computing, plus waiting in the queue for a GPU.
So we cut a small square (about 0.4 mm wide) out of each of the two slides:

- **8,697 cells**, about 0.6% of the full data, in exactly the same file layout;
- only whole cells kept, none cut in half at the border;
- the cut is clean: every cell's protein brightness (3 channels checked) is identical to
  the same cell in the original slide.

The practice data lives in `test_data/`, and everything a practice run produces goes into
`test_data/results/`, never into the real results. Because the input is different, its
results are **never compared** with the real run. The run is judged on its own:

1. All four steps finished without errors.
2. Snakemake confirms that every expected output exists.
3. No step log contains an error.
4. Nothing was written outside `test_data/results/`.
5. The steps agree with each other, e.g. the number of cells passed from one step to the
   next adds up.

**All five passed.** The whole workflow, GPU step included, ran in about 18 minutes.

### Hiccups along the way

- One commit accidentally left out `run.sh`. It was added to that same commit before
  anything was uploaded.
- The practice-data script still used a setting name that had been renamed. It was updated.
- A test launched from the wrong folder created an empty helper folder there. It was
  removed.

## How to use it now

```bash
./run.sh dry       # dry run: show the plan
./run.sh submit    # run everything on the cluster
```

**Practice data:** add `--configfile config/config.test.yaml` to either command. This works
only in the cluster copy, because the practice settings and `test_data/` are not uploaded.

**A new user** only has to:

1. install Snakemake once: `micromamba create -n snakemake_env -f envs/snakemake_env.yml`
2. copy `config/config.example.yaml` to `config/config.yaml` and fill in which slides to
   use, where the data is, where results should go, and where `integration.sif` is
3. copy `profiles/slurm/config.yaml.example` to `profiles/slurm/config.yaml`
4. run `./run.sh dry`, then `./run.sh submit`

More detail: [README.md](../README.md) (how to run it, what input data it needs),
[config/README.md](../config/README.md) (every setting),
[container/README.md](../container/README.md) (what is in the container, how to rebuild
it).

## Noticed, but deliberately not touched

Fixing these would change results or behaviour, so they are left for a separate decision.

| Step | Issue |
|------|-------|
| 1a | Searches the raw data folders more thoroughly than needed: slow but harmless. |
| 1b | The setting `compartment: nucleus` would give protein values to the wrong cells for this Xenium version. The runs use `cell`, which is correct. |
| 2 | Reads *every* protein file in its input folder, so leftover files from an older run could sneak in. |
| 3 | Runs only on a machine with a GPU. |
| 3 | The plots (dotplot, matrixplot) choose marker genes with slightly different rules than the marker table. |
| 3 | The setting `min_nonzero_prop` has no effect: the dotplot uses a fixed value of 0.1. |
| 3 | Changing `top_rna`, `top_pro` or `bf_rna` doesn't make Snakemake redo step 3; the rerun has to be forced. |
