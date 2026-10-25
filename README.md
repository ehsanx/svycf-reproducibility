# Reproduction materials: design-based central limit theorem for cross-fitted influence functions

Code, recorded results and checks for

> Karim, M. E. A design-based central limit theorem for cross-fitted influence functions under stratified multistage
> sampling. *(Journal details will be added on publication.)*

M. Ehsan Karim, University of British Columbia (ORCID [0000-0002-0346-2871](https://orcid.org/0000-0002-0346-2871)).
This is version v1.0.0 of the reproduction materials.

The paper has two numerical parts:

- a **simulation study** on synthetic finite populations (Section 4; Supplement S2 and S3);
- a **resampling study** on three Service Provision Assessment (SPA) surveys (Section 5; Supplement S4).

This repository contains the implementations and reproduction materials for both. The supplied summaries reproduce
the paper's tables, figure and reported numbers without refitting models, and the sections below explain how to
rerun each study from its required inputs. The per-replicate simulation output is in two archives attached to the
release ([Recorded results](#recorded-results)). The repository contains no survey records.

## Contents

| Path | Contents |
|---|---|
| `src/simulation/` | The simulation engine: an R project with its own estimators, sampler and variance code (`R/`), job entry points (`runs/`), the locked environment (`renv.lock`, `renv/`, `env/`), quadrature rules (`grids/`) and unit tests (`tests/`) |
| `src/application/R/` | The SPA analysis: unpacking, count audit, resampling and estimation |
| `src/displays/` | The script that builds every simulation and SPA display and number of the paper |
| `scripts/simulation/`, `scripts/application/` | Setup, Slurm submission and job scripts; assembly of the recorded per-replicate output |
| `config/` | Settings and machine-readable records: cell tables, frozen inputs, SPA settings and package pins, data files, cluster settings template, Python versions, identities, provenance and archive checksums |
| `results/simulation/`, `results/application/` | Recorded summaries of both studies |
| `results/displays/` | The paper's figure, tables and number macros, as built from `results/` |
| `results/index-cells.csv`, `results/index-displays.csv` | Index from simulation cells and SPA surveys to result files and to the paper's tables and figure |
| `tests/` | Checks: display rebuild, provenance, audit of the per-replicate output, task-7 recovery, independent recomputation of key estimates |
| `vendor/svytmle_0.1.2.tar.gz` | The exact `svytmle` 0.1.2 build used by the SPA study (md5 `c976da847d9d28069123835b4bb71281`) |

Run outputs go to `outputs/` and restricted survey data to `data/`; both are ignored by git.

## Requirements and setup

Run every command from the repository root.

- **Displays:** Python 3.10 with the versions in `config/python-requirements.txt` (pandas 2.3.3, numpy 2.2.6,
  matplotlib 3.10.9).
- **Simulation:** R with a C/C++/Fortran toolchain and a Slurm cluster. The recorded runs used R 4.4.0 on Linux
  (UBC ARC Sockeye); `src/simulation/renv.lock` records R 4.5.1, the version used for the inputs built locally, and
  setup accepts either. The 69 locked packages are restored by `scripts/simulation/setup.sh` (Route 2).
- **SPA study:** approved access to the DHS Program data, R 4.4.0 (the recorded version) with a toolchain, and a Slurm
  cluster. The pinned packages are installed by `scripts/application/setup.sh` (Route 3).

Simulation settings and input checksums are recorded in `config/simulation/frozen-design.json` and the files it
indexes; SPA settings and package versions are in `config/application/`. Use these settings with the commands below
to reproduce the reported analyses.

## Paper-to-code map

Every display of the paper is built by `src/displays/make_sim_displays.py` (with `src/displays/app_displays.py` for
the SPA displays) from the files in `results/simulation/` and `results/application/`:

| Output | In the paper |
|---|---|
| `results/displays/main/fig1_simulation.pdf` | Figure 1 |
| `results/displays/main/tab1_estimates.tex`, `tab1_coverage.tex` | Table 1, panels A and B |
| `results/displays/main/tab2_application.tex` | Table 2 |
| `results/displays/main/sim_numbers.tex` | every number quoted in the text, as LaTeX macros (`\Sim...`, `\App...`) |
| `results/displays/supplement/tabS_*.tex` | Supplement tables |
| `results/displays/extended/tabS_*.tex`, `rec_index.tex` | Tables not printed in the paper: oracle cells at fixed populations, the unequal-allocation pair, the deterministic array, per-population summaries, preparation and dry-run checks, and every estimate with its classification against its benchmark region |

`results/index-displays.csv` gives, for each display file, its label and number in the paper and its inputs;
`results/index-cells.csv` gives, for each simulation cell and SPA survey, where it is reported and its result files.

## Three routes

| Route | Needs | Recorded cost |
|---|---|---|
| [1. Rebuild the displays](#1-rebuild-the-displays) from the recorded results | Python 3.10 with pandas, numpy and matplotlib | under a minute |
| [2. Rerun the simulation](#2-rerun-the-simulation) | R 4.4 or 4.5, a C/C++/Fortran toolchain, a Slurm cluster | 807 allocated core-hours |
| [3. Rerun the SPA study](#3-rerun-the-spa-study) | Approved access to DHS Program data, R, a C/C++/Fortran toolchain, a Slurm cluster | 287 allocated core-hours |

Allocated core-hours are elapsed time times allocated cores, summed over the recorded jobs on UBC ARC Sockeye.

### 1. Rebuild the displays

No model is fitted and no survey data are needed.

```
pip install -r config/python-requirements.txt
python tests/check_displays.py
```

`check_displays.py` runs `src/displays/make_sim_displays.py` into a temporary folder and compares every file, byte for
byte, with `results/displays/`; it prints `38 of 38 display files identical` and exits 0. To write the files in
place instead, run `python src/displays/make_sim_displays.py` (set `DISPLAY_OUT` to write elsewhere). The script
reads only `results/simulation/` and `results/application/` (inputs listed at the top of the script); the outputs
are listed in the [Paper-to-code map](#paper-to-code-map).

### 2. Rerun the simulation

The study has 41 cells (`config/simulation/production-cells.tsv`): many-strata and fixed-strata trajectories, oracle
cells, sampling-fraction and strictness cells, census-only cells, comparison and control cells, and a deterministic
array. Each cell is a Slurm array of one-core tasks; every task draws its populations and samples from keyed
L'Ecuyer-CMRG streams, so a task reproduces the same draws whenever and wherever it runs.

1. **Cluster settings.** `cp config/slurm.env.example config/slurm.env` and fill in your account, partition and
   modules. The recorded runs used `SIM_MODULES="gcc/9.4.0 r/4.4.0"`.
2. **Environment** (once, on a node with internet access; 30 to 60 minutes, packages are compiled):
   ```
   bash scripts/simulation/setup.sh
   ```
   It restores the 69 locked packages of `src/simulation/renv.lock` into `src/simulation/renv/library/`, with
   `svytmle` from `vendor/`; it never touches your own R library. It checks every package against the lock, notes an
   R version other than the two used for the recorded runs (4.4.0 on the cluster, 4.5.1 locally), and writes
   `PACKAGE-COMMIT.txt`, the package identity every run records ([Provenance](#provenance)). The engine's unit
   tests then run with
   `bash -c 'source scripts/simulation/env.sh && Rscript src/simulation/env/run.R src/simulation/tests/testthat.R'`.
3. **Preparation** (6 cells, 21 one-core tasks, at most 13 core-hours): census quantities of every production
   population, the frozen κ_AE record used by the control cell, and the preparation checks.
   ```
   bash scripts/simulation/submit.sh list prep
   bash scripts/simulation/submit.sh prep
   ```
   The recorded κ_AE record (`results/simulation/preparation/prep-panel/prep-panel-20261012T081105Z/kappa_ae.json`)
   is installed first; the panel job then verifies it against freshly computed census components instead of
   deriving it again. `--fresh-kappa` derives it afresh instead.
4. **Production** (879 one-core tasks; the recorded time limits, `config/simulation/production-limits.csv`, total
   2,257 core-hours, and the recorded run used 807):
   ```
   bash scripts/simulation/submit.sh production --then-grade
   ```
   `--then-grade` also queues step 5, which starts only after every production job has ended. `--only P-ms-H25,P-fs-m50`
   submits selected cells, and `--print` shows every `sbatch` command without submitting. Recorded time and memory
   per cell: `results/simulation/production_resources.csv`. One task without Slurm:
   `bash scripts/simulation/submit.sh local production P-ms-H25 1` (task 1 of 40; about an hour).
5. **Summaries and Monte Carlo intervals** (one core, about 3.5 hours, 1 GB). Queued by `--then-grade`, or
   separately with
   ```
   bash scripts/simulation/submit.sh grade --after-production
   ```
   which waits for every job of the latest production submission (`--after ID:ID:...` names jobs explicitly). It
   writes `outputs/simulation/grading/grade-<UTC>/` with the same files as `results/simulation/grading/`: every
   statistic with its 95% population-block bootstrap interval (B = 10,000 keyed resamples), completeness against the
   design, implementation checks, fitting failures and provenance checks. It reads the preparation and dry-run
   records in `results/simulation/` (`--preparation outputs` uses your own preparation outputs). Exit status 0 means
   every cell is complete, no statistic is withheld and every provenance check holds; 3, a cell is incomplete or a
   statistic withheld; 4, a provenance check fails; 5, an error in some cell.

To regrade the recorded output without simulating, assemble it (below) and run
`bash scripts/simulation/submit.sh grade --in <target>/production`.

To build the displays from your own summaries, copy the files of your `grade-<UTC>/` folder over
`results/simulation/grading/` in a working copy and run `python src/displays/make_sim_displays.py`.

### 3. Rerun the SPA study

The SPA microdata (Malawi 2013–14, Haiti 2013, Haiti 2017–18) are available to registered users from
[The DHS Program](https://dhsprogram.com) and may not be redistributed.

1. **Data.** Request access to the Haiti and Malawi SPA surveys and download the 24 Stata zip files listed in
   `config/application/data-files.csv` (name, survey, file type, size, md5). Put them in `data/spa/zips/`, or in
   `$SPA_RAW_DIR/zips/` for a private folder elsewhere. Keep the folder private: the jobs refuse to run when it allows
   group or other access.
2. **Packages** (once, with internet access; packages are compiled): `bash scripts/application/setup.sh`
   (`--dry-run` lists what would be installed). It installs the 62 CRAN packages of `config/application/packages.csv`
   at their pinned versions from the CRAN snapshot of 2026-09-24, and the vendored `svytmle` build, into their own
   library `src/application/library/` (or `SPA_R_LIBRARY`). Your own R library is neither read nor written, and every
   SPA job uses only that library and R's own.
3. **Preparation** (one job): checks the md5 of every zip, unpacks them into `data/spa/`, runs the per-facility count
   audit (`outputs/application/count-audit_*.csv`) and the synthetic self-tests.
   ```
   bash scripts/application/submit.sh prep
   ```
4. **Resampling** (per survey: census overlap, 8 chunks of 64 resamples on 32 cores, 20 reference splits on 20 cores,
   then collection and summaries):
   ```
   bash scripts/application/submit.sh run MW_SPA_2013-14 HT_SPA_2017-18 HT_SPA_2013
   ```
   Walltimes default to about three times the recorded run times (`SPA_RUN_TIME`, `SPA_REF_TIME` in
   `config/slurm.env`). Every resample and split sets its own seed, so results do not depend on chunking or core count.
5. **Frame sizes** for the displays:
   ```
   python scripts/application/extract_frames.py outputs/application/RESAMPLE-RUNS.md > outputs/application/frames.csv
   ```

The outputs in `outputs/application/` have the names and columns of `results/application/`. The resampling refuses to
open any data file unless the settings in `config/application/` are complete and the installed `svytmle` carries the
pinned build stamp, and it checks the md5 of its own code first ([Provenance](#provenance)). Record-level frames stay
in memory; only aggregates are written.

## Recorded results

### Simulation (`results/simulation/`)

| Path | Contents |
|---|---|
| `grading/grades.csv` | Every statistic of every cell and arm: estimate, 95% interval, benchmark region, classification, basis, failed and unresolved replicates |
| `grading/claims.csv`, `population_indexed_summary.csv` | One classification per cell and statistic; population-indexed counts for design-conditional cells |
| `grading/completeness.csv`, `required_checks.csv`, `failed_arms.csv`, `blocked.csv` | Completeness against the design (every cell complete), implementation checks run in every replicate (13,067,483 required, 1 failed), fitting failures (none), statistics withheld by a failed check |
| `grading/provenance.csv`, `provenance_checks.csv`, `round3_inputs.csv`, `grading_manifest.json` | Run folders read and their identities; checksums of the preparation and dry-run records read |
| `grading/displays/*.csv` | The tables behind Figure 1 and the simulation tables |
| `preparation/` | Census quantities of the production populations, the frozen κ_AE record and the preparation checks |
| `execution-check/` | Summary of the dry run of every cell (implementation and running time), with the learner-rate and control-calibration checks the summaries read |
| `unresolved_t2_by_cell_arm.csv` | Replicates whose remainder T₂ could not be resolved numerically (80; recorded as missing, estimates kept) |
| `production_resources.csv` | Per cell: array tasks, longest task, time limit, peak memory |
| `task7_timed_out_printed.csv` | Values printed by a run of task 7 of P-e3-ctrl-m400 that timed out; `tests/verify_task7.py` compares them with the archived task |

The benchmark region of a statistic (`region_lo`, `region_hi` in `grades.csv`) is what the paper calls its benchmark
(Table S5). The `disposition` column compares the statistic's 95% Monte Carlo interval with that region; the paper
reports the estimates and intervals. The values are:

| Disposition (leading text) | Meaning |
|---|---|
| `PASS` | The whole interval lies inside the benchmark region |
| `RED-FLAG` | The whole interval lies outside the benchmark region |
| `UNGRADABLE (MC precision)` | The interval overlaps a boundary of the region, so the comparison is unresolved at the simulated precision |
| `UNGRADABLE (premise: ...)` | A condition the comparison needs was not established; the text in parentheses names it, and the estimate is descriptive |
| `BOUNDARY` | The throttled learner at the boundary rate $\alpha=1/4$, reported descriptively |
| `FAILS AS PREDICTED` | A control condition left its region, as it was designed to |
| `DIAGNOSTIC` | Descriptive output that is not compared with a region |

`basis` records whether a row was compared with its region (`graded`), is descriptive, or is derived from other rows.
In `blocked.csv`, a statistic is withheld because an implementation check of its cell failed; `failed_arms.csv` lists
fitting failures (none). Missing values are `NA` in CSV files and `---` in tables.

One implementation check failed. In the with-replacement code-check cell `P-wr-calib`, 1 of its 112,016 required
checks failed, in replicate 4611 of 8,000 (population 1): the total of the linearized score computed by the `survey`
package and the direct cluster-total formula disagreed beyond the tolerance
`1e-12 * max(abs(total_direct), 1e-300) + 1e-18`, about 6.5e-18 at the saved total 5.5062447824952296e-06. The
comparison total was not saved, so the size of the discrepancy is not known. The variance comparison in the same
check passed. Floating-point cancellation is a possible explanation, not a confirmed one. The cell is a code check
and enters no result reported in the paper; `required_checks.csv` keeps the failure and `blocked.csv` the withheld
statistic.

### Per-replicate simulation output (two archives)

| File | Size | Files | Contents |
|---|---|---|---|
| `simulation-production-outputs.tar.gz` | 280 MB | 6,099 | Every run folder of the 41 production cells except task 7 of P-e3-ctrl-m400 |
| `simulation-task7-recovery.tar.gz` | 77 kB | 9 | Task 7 of P-e3-ctrl-m400 |

Assemble both for the complete output. SHA-256 and md5: `config/assets.csv`.
Download them from this repository's v1.0.0 release. Each run folder holds
`results.csv` (one row per replicate and arm: estimate, variance estimate, interval, remainders, status), `census.csv`
and `census_strata.csv` (census targets), `checks.csv` (implementation checks), `timing.csv`, `t2_refinement.csv` and
`manifest.json` (settings, identities, environment; the cluster path prefix is replaced by `<run-root>`). Run
manifests retain the settings, software versions and labels recorded during execution.

```
python scripts/simulation/assemble_production.py <folder with the two archives> outputs/recorded
python tests/audit_production.py outputs/recorded/production config/simulation/production-cells.tsv
python tests/verify_task7.py outputs/recorded/production/P-e3-ctrl-m400 results/simulation/task7_timed_out_printed.csv
python tests/recompute_grades.py outputs/recorded/production outputs/recorded/production/P-e3-ctrl-m400 results/simulation/grading/grades.csv
python tests/check_provenance.py outputs/recorded/production
```

The assembly checks both archives and every extracted file, and reports 41 cells, 879 run folders and 933,200 result
rows. The audit counts rows, required checks, fitting failures and unresolved T₂ per cell; `recompute_grades.py`
recomputes key estimates with plain statistics (and approximate bootstrap intervals) for comparison with `grades.csv`.

### SPA study (`results/application/`)

Per survey (`MW_SPA_2013-14`, `HT_SPA_2017-18`, `HT_SPA_2013`): `resample_<survey>_f250_replicates.csv` (500
resamples × 8 methods: estimate, standard error, degrees of freedom, interval, status),
`_f250_calibration.csv` (variance ratio with its Monte Carlo standard error, coverage of the whole-frame reference
value), `_reference.csv` and `_reference_splits.csv` (whole-frame reference values: each estimator family run on the
whole frame over 20 fold splits), `_census_overlap.csv`, `_f250_diagnostics_summary.csv`
and `_f250_replicate_diagnostics.csv` (weights, sample sizes, truncation), `_calibration_all_f.csv`,
`_environment.csv` (R and package versions); `count-audit_<survey>_*.csv` (frame counts); `frames.csv` (frame sizes).

The `method` values in these files correspond to the paper's labels as follows:

| `method` in the CSV files | Label in the paper |
|---|---|
| `Fully-Aware` | TMLE, single fit, design-aware |
| `Partially-Aware` | TMLE, single fit, weights only |
| `Non-Aware` | TMLE, single fit, no design |
| `Fully-Aware-CV` | TMLE, internal CV, design-aware |
| `Fully-Aware-CF` | cross-fitted TMLE |
| `AIPW-CF` | cross-fitted AIPW |
| `Fully-Aware-CF [unw+design C]` | cross-fitted TMLE, sensitivity |
| `AIPW-CF [unw+design C]` | cross-fitted AIPW, sensitivity |

### Indexes

`results/index-cells.csv` gives, for each simulation cell and SPA survey, its coordinates, question, Monte Carlo design,
where it is reported in the paper and its result files. `results/index-displays.csv` gives, for each display file, its
label and number in the paper and its inputs.

## Provenance

The supplied results were generated using the implementations and software versions identified in `config/`
(`execution-identity.json`, the locked and pinned package versions). This repository includes adaptations for
portable paths and execution outside the original computing environment; its code comments, and some diagnostic
messages, were rewritten for readers. In the engine files whose messages were rewritten, every other token is
unchanged; `config/execution-identity.json` lists them as adapted. File checksums, and the relationship between each
distributed file and its original implementation (identical, comments edited only, or adapted), are recorded in
`config/provenance.csv`. `python tests/check_provenance.py`
checks every file's integrity and that every file marked "comments edited" has the executable code of its original,
and `python tests/check_displays.py` reproduces the published displays. Runs made with this repository record its
identifier from `config/execution-identity.json` as their package identity.

## Implementation notes and limitations

- The simulation engine is an R project written for this paper, not an installable package. It uses its own
  estimator, fold and variance code and does not call `svytmle`.
- Production, preparation and summary runs refuse to start if any frozen input differs from
  `config/simulation/frozen-design.json`; development runs (`--mode dev`) do not check, for experimentation.
- The summaries read the recorded dry-run summary in `results/simulation/execution-check/`; the wrappers do not rerun
  the dry run.
- Numerical results can differ in the last digits on another platform or R version; the recorded runs used R 4.4.0
  on Linux with the locked package versions.
- The SPA package pins follow the CRAN snapshot of 2026-09-24. The recorded runs logged the versions of `svytmle` and
  six key packages (SuperLearner, tmle, survey, earth, glmnet, ranger), all equal to the pins; the remaining
  dependencies, including recommended packages such as Matrix and survival, are pinned to the snapshot's versions.
- The `svytmle` build stamp `Packaged: 2026-09-24 19:40:00 UTC; wilds` includes the user name of the machine that built
  the package. It is part of the unchanged vendored file and is the pin the code enforces, so it is kept, in the
  tarball and where the runs recorded it.
- Remainders T₁ and T₂ of the R-a TMLE arm are kept in the per-replicate output but are not used in any summary or
  display.

## Data access and privacy

No survey records are included. The SPA results are aggregates over 500 resamples or over the census frame; frame
counts are reported for whole surveys and strata. Use of the SPA microdata is governed by The DHS Program's terms of
use. The simulation output is synthetic.

## Licence and citation

The code is released under the MIT licence (`LICENSE`). The recorded results (the files in `results/` and the two
archives of per-replicate output) are released under the Creative Commons Attribution 4.0 International licence
([CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)). Two files come from other projects and keep their own
licences: `vendor/svytmle_0.1.2.tar.gz` (GPL >= 3) and `src/simulation/renv/activate.R`, the bootstrap script of renv
(MIT). The R and Python packages that the setup scripts install are not distributed here and have their own licences.

To cite these materials, cite M. Ehsan Karim (2026), Reproduction materials: design-based central limit theorem for
cross-fitted influence functions, version v1.0.0, GitHub,
https://github.com/ehsanx/svycf-reproducibility/releases/tag/v1.0.0. The release identifies the exact source commit.
The paper citation above will be updated on publication. Citation metadata: `CITATION.cff`.
