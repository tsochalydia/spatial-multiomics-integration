# What happened when we tidied up this project

*A plain-language explanation. No programming knowledge needed.*

---

## 1. What this project does

The project takes two kinds of microscope data from the same tissue slides:

- **Xenium**: which genes (RNA) are active in each cell
- **CODEX**: how much of 24 proteins each cell has

It combines them and sorts the cells into groups (clusters) that can then be named, for
example "T cells" or "macrophages".

It does this in **four steps**, always in the same order:

| Step | In everyday words |
|------|-------------------|
| 1a | Collect the gene counts of all cells from all slides into one table |
| 1b | Measure, for every cell, how bright each protein is in the CODEX images |
| 2  | Remove bad-quality cells and put genes and proteins together |
| 3  | Train a model on both, then group similar cells and draw the plots (needs a graphics card, GPU) |

Think of it as a **recipe**: four cooking steps, a list of ingredients, and some
settings like oven temperature and cooking time.

A program called **Snakemake** is the "kitchen manager". It reads the recipe, starts
each step when the previous one is done, and **skips steps that are already finished**.
The heavy work runs on the university's computer cluster, which has a queue system
(**SLURM**): you ask for a machine, wait your turn, and your job runs.

All the software the steps need is packed into one big file, `integration.sif`. You
can picture it as a **sealed lunchbox**: whatever computer opens it gets exactly the
same software versions.

---

## 2. Why it needed tidying

The recipe worked. The results for the thesis were already produced. But the project
folder was hard for anyone else to use:

1. **Personal addresses were written inside the files.** Many files contained exact
   folder locations on *Lydia's* account, like `/cluster/project/moor/lydia/...`. Anyone
   else would have had to find and change them one by one.
2. **Settings were spread over several places.** Some were in `config.yaml`, some in a
   cluster settings file, some inside the launch scripts.
3. **Everything was in one drawer.** The recipe, the software lists, the lunchbox
   building instructions and the settings all sat together, so a newcomer couldn't tell
   what is what.
4. **There were two almost identical "start" buttons** (`run.sh` and `submit.sh`), so
   every fix had to be made twice.
5. **A test run would have overwritten the real run's notes.** The recipe always wrote
   its log files into the same folder, whatever data it was run on.

**The goal:** tidy all of this **without changing what the analysis calculates**.

---

## 3. What we did, step by step

We made **one small change at a time**. After each change we checked that the recipe
still did exactly the same thing before moving on. Every change was saved separately in
**git**, the project's history book, so any single change can be looked at or undone.

### Phase 1: remove personal addresses

| # | What changed | Why | Does it change the results? |
|---|--------------|-----|-----------------------------|
| 1.1 | The recipe reads the settings file you actually give it, and the log folder can be chosen | So a test run can use its own settings and doesn't overwrite the real run's logs | No |
| 1.2 | One setting, `results_dir` ("put all results here"), replaces five separate output addresses | Fewer things to fill in. The sub-folder names (`step1a`, `step1b`, …) stayed exactly the same, so the finished results are still recognised. | No |
| 1.3 | Short addresses like `results` are allowed (meaning "inside the project folder") | Works for anyone who downloads the project, whatever their account | No |
| 1.4 | The lunchbox settings shrank to one line. The recipe works out by itself which folders the lunchbox may look into. | Less to fill in, and fewer places to make mistakes | No |
| 1.5 | One start button, `run.sh`, finds the software by itself and has a `submit` option that replaces `submit.sh` | No personal addresses, and fixes happen in one place | No |
| 1.6 | Removed personal addresses from the software lists and fixed a wrong name in one of them | Anyone can use the lists | No |

### Phase 2: all settings in one place

| # | What changed | Why | Does it change the results? |
|---|--------------|-----|-----------------------------|
| 2.1 | All settings files moved into a `config/` folder: the real settings, the test settings, and an example to copy | One place to look | No |
| 2.2 | Cluster-specific names (which queue to use, e.g. `normal.24h` or `gpu.24h`) moved to the cluster settings file (`profiles/slurm/`). How much memory and time each step needs stayed in `config/`. | On another cluster, only the cluster file changes | No |

### Phase 3: clear folders

| # | What changed | Why | Does it change the results? |
|---|--------------|-----|-----------------------------|
| 3.1 | Files moved into labelled folders (see the map below) | A newcomer sees at once what is what. It's the standard layout for Snakemake projects. | No |
| 3.2 | The "do not upload" list (`.gitignore`) now also covers results, backup copies and a few temporary files | Personal settings and big files can't be uploaded by accident | No |

### The folder map, before and after

**Before**, everything in one drawer:

```
Snakefile  scripts/  config.yaml  config.yaml.example  run.sh  submit.sh
spatial_env.yml  scvi_env.yml  int_retrieval_env.yml  snakemake_env.yml
micromamba.dockerfile  entrypoint.sh  profiles/  README.md
integration.sif  integration.tar  test/  logs/
```

**After**, labelled drawers:

```
README.md, run.sh      how to start (the only button you press)
config/                YOUR SETTINGS: which data, where results go, how big the jobs are
profiles/slurm/        CLUSTER SETTINGS: which queues to use
workflow/              THE RECIPE: Snakefile + the four step scripts
envs/                  SOFTWARE LISTS: exact versions inside the lunchbox
container/             LUNCHBOX INSTRUCTIONS: how integration.sif is built
test_data/             a tiny practice dataset (kept on the cluster, not uploaded)

integration.sif, integration.tar, test/, logs/   unchanged, not moved, not uploaded
```

---

## 4. What did NOT change

- **The science.** The four step scripts, which do the actual analysis, were only moved
  to a new folder. Not a single line in them was changed.
- **The existing results.** The full results in `test/` and the old logs in `logs/` were
  not touched, moved or deleted.
- **The big files.** `integration.sif` (13 GB) and `integration.tar` (36 GB) stay where
  they were.
- **Nothing was lost.** The only file removed is `submit.sh`, whose job is now done by
  `run.sh submit`, and it's still in the git history. Before the bigger edits of the
  personal settings files, a backup copy was saved next to them (the files ending in
  `.backup…`).

---

## 5. How we made sure nothing broke

### The rehearsal: "dry run"

Snakemake can do a **dry run**: it reads the recipe and prints exactly what it *would*
do (which steps, which files, which commands) without doing anything. It takes seconds
and costs nothing.

- **Before** any change, we saved the dry-run plan as the **reference**.
- **After every** change, we made a new plan and compared it with the reference, line by
  line.
- The only differences allowed were the ones we intended: new file locations. All
  inputs, outputs, software, memory, time and every setting passed to the scripts stayed
  the same.
- The dry run also kept saying **"Nothing to be done"** for the real data. That means
  Snakemake still sees the finished results as complete and would not redo hours of
  work.

### A tiny practice dataset

Running the full data takes more than 10 hours of computing, plus waiting in the queue
for a graphics card.
So we made a **tiny practice dataset** by cutting a small square (about 0.4 mm wide) out
of each of the two slides:

- It holds **8,697 cells**, about 0.6% of the full data, in exactly the same file layout
  as the real data.
- Only whole cells were kept, never cells cut in half at the border.
- We checked that the cut was clean: for every cell in the small piece, the protein
  brightness (in the three protein channels we checked) is identical to the same cell in
  the original slide.

The practice dataset lets the whole recipe run on the cluster in minutes instead of
hours. It lives in `test_data/`, and **everything a practice run produces goes into
`test_data/results/`**, never into the real results.

### The practice run: what we check

The practice run is checked **only on its own terms**. Results from different input
data are never compared.

1. All four steps finished without errors.
2. Snakemake confirms that every expected output exists.
3. No step's notes (log) contain an error.
4. Nothing was written outside `test_data/results/`.
5. The steps agree with each other: for example, the number of cells passed from one
   step to the next adds up.

**Result: all five checks passed.** The tidied recipe ran from start to finish on the
cluster in about 18 minutes, including the graphics-card step. Everything it produced
went into `test_data/results/`, and nothing else was changed.

### Small hiccups along the way (and how they were fixed)

- One saved change accidentally left out a file (`run.sh`). This was noticed and added to
  that same change before anything was uploaded.
- The practice-data script still used a setting name that had been renamed. It was
  updated.
- A test launched from the wrong folder created an empty helper folder there. It was
  removed; it was empty.

---

## 6. How to use the project now

```bash
./run.sh dry                                          # rehearsal: show the plan
./run.sh submit                                       # run everything on the cluster
./run.sh dry --configfile config/config.test.yaml     # rehearsal on the practice data
./run.sh submit --configfile config/config.test.yaml  # run on the practice data
```

A **new person** only has to:

1. install the kitchen manager once:
   `micromamba create -n snakemake_env -f envs/snakemake_env.yml`,
2. copy `config/config.example.yaml` to `config/config.yaml` and fill in where their
   data is, where results should go, and where the lunchbox file (`integration.sif`) is,
3. copy `profiles/slurm/config.yaml.example` to `profiles/slurm/config.yaml`,
4. press `./run.sh dry`, then `./run.sh submit`.

---

## 7. Still to do

- Rewrite the main `README.md` for newcomers, and add short guides for the settings and
  for rebuilding the lunchbox.
- Final check: run the **original, untidied recipe** and the **tidied recipe** on **the
  same data** and confirm they give the same results.
- Offer all changes for review on GitHub (a "pull request") before they become the
  official version.

## Things we noticed but deliberately did not touch

These would change results or behaviour, so they are left for a separate decision:

- **Step 3 marker rules:** the plots (dotplot, matrixplot) choose marker genes with
  slightly different rules than the marker table.
- **Step 2 leftover files:** it reads *every* protein file in its input folder, so
  leftover files from an older run could sneak in.
- **Step 3 needs a GPU:** it can only run on a machine with a graphics card.
- **Step 1a searches slowly:** it searches the raw data folders more thoroughly than
  needed, which is slow but harmless.
- **Step 1b nucleus option is broken:** the setting `compartment: nucleus` would give
  protein values to the wrong cells for this Xenium version. The runs use `cell`, which
  is correct.
- **Step 3 ignores a setting:** `min_nonzero_prop` has no effect; the dotplot uses a
  fixed value of 0.1.
- **Some step 3 changes aren't noticed:** changing `top_rna`, `top_pro` or `bf_rna`
  does not make the recipe redo step 3 by itself (it has to be forced).

---

## Small dictionary

| Word | Meaning here |
|------|--------------|
| **Path / address** | Where a file lives, like `/cluster/project/moor/lydia/...` |
| **Hard-coded path** | An address typed directly into a file instead of being a setting |
| **Config file** | The settings sheet of the recipe |
| **Snakemake** | The kitchen manager that runs the steps in order and skips finished ones |
| **SLURM / cluster** | The shared university computers and their waiting queue |
| **Container (`.sif`)** | The sealed lunchbox with all the software inside |
| **Bind** | Letting the lunchbox look into a folder on the cluster |
| **Dry run** | A rehearsal that prints the plan but runs nothing |
| **git / commit** | The project's history book / one saved change in it |
| **Tag `pre-reorg`** | A bookmark in the history at the moment before the tidy-up, to go back to anytime |
| **`.gitignore`** | The "do not upload" list |
