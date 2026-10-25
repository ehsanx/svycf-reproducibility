# THE SIMULATION RUNNER for the fitted, oracle and calibration cells (E1, the sampled E2 cells and E3): one cell,
# a block of populations and their replicate samples, every arm the cell defines, the census targets, the
# remainders and every implementation check. Every quadrature rule comes from the saved, sha256-verified artifact
# grids/quad-rules-dev-v1.tsv (required store).
#
#   Rscript src/simulation/env/run.R src/simulation/runs/run_e2e_dev.R \
#       [--regime many] [--H 25] [--B 2] [--R 10] [--u9-bf 4] [--out-root <dir>]
#   Rscript src/simulation/env/run.R src/simulation/runs/run_e2e_dev.R \
#       --regime fixed --m 200 [--B 2] [--R 5] [--u9-bf 4] [--out-root <dir>]
# The production cells and their arguments are in config/simulation/production-cells.tsv.
#
# Cells:
#   --regime many (default): E1 many strata, m_h = 2, J = 20, H in {25, 50, 100, 200, 400, 800}, m = 2H
#     (prefix views of the 800 x 20 master "MS"). Arms: oracle, Ra-AIPW, Ra-TMLE, Rb-AIPW, Rb-TMLE
#     (m = 1600: oracle and Ra-AIPW only).
#   --regime fixed --m m: E1 fixed strata, H = 5, m_h = m/5, J = 10 m_h (f = 0.10), m in {50, ..., 3200}
#     (prefix views of the 5 x 6400 master "FS"). Arms: oracle, Ra-AIPW, Ra-TMLE (to m = 800) and the R-d
#     AIPW rungs Rd-AIPW-a0xx (R/learner_rd.R, from m = 200; m = 3200 is R-d only, with the internal
#     psi_oracle). No R-b in fixed strata.
#   --regime many --e2c fpc --J j: the E2 sampling-fraction cells: H = 50, m_h = 2, J in {100, 20, 8, 4}
#     (f = 0.02, 0.10, 0.25, 0.50), prefixes of the master "FPC", n0 = 20 (G2: n0/M_k >= 0.40). Arms: oracle
#     at every J; R-a and R-b PSU-CF AIPW at J = 4 only; no TMLE. Writes the census components per stratum
#     (census_strata.csv), the FPC, the G1 quantities and the per-sample UC and two-stage WOR SEs.
#   --regime many --e2c strict --H 100: the E2 strict fitted cell (m = 200, J = 4, f = 0.5, master "MS"):
#     oracle, Ra-AIPW, Rb-AIPW. Its census limb is runs/run_e2_census.R.
#   --regime many --e2c wr --oracle-only 1 [--calib 1]: the E2 with-replacement code check: H = 50, m_h = 2,
#     J = 4, a prefix of the FPC master; WR first stage with selection-keyed stage 2 (draw_sample_wr); every WR
#     variance clusters on the selection. --calib 1 checks, over redraws, the census-centred linearised score
#     with known N against the literal WR target sigma2_m,WR and U2's WR form; the Hajek WR estimator is
#     reported beside it, unchecked.
#   --oracle-only 1 --calib 1: the calibration cell -- U2 and U8's three limbs over redraws, |z| <= 3, plus
#     C3a's both-levels quantities reported without a grade. U8's limbs are checked on the census-centred
#     linearised score with known N (lin_var_designs, R/variance.R), whose design expectations equal the
#     census targets exactly; the Hajek estimators are reported beside them, unchecked.
#   --regime fixed --m 400 --u10 1 --oracle-only 1: the U10 cell (G8; cell label e2c = "u10"): the oracle
#     with the stage-2 weights omitted (R/controls.R C-U10-1), checked over redraws against its exact linear
#     identities E[psi_U10] = Psi_U10 and E[se2_U10] = sigma2_m,U10; its Z keyed to Psi_N and the G6 limb (a)
#     boundary distance are reported.
#   --regime fixed --e3 panel --m 400: the E3 panel: oracle, Ra-AIPW (PSU-CF), Ra-AIPW-nocf, R-c (fit_rc,
#     R/learner_rb.R) as Rc-AIPW-psucf, -nocf, -indcf, and the coarsened individual-CF arm.
#   --regime fixed --e3 vsens --m 200 --V 2|5: the E3 V-sensitivity cell: oracle and Rb-AIPW.
#   --regime fixed --e3 ctrl --m 400: the E3 power-control panel: oracle, the host Ra-AIPW (PSU-CF) and the
#     engineered aligned-error arm AE-AIPW (R/controls.R C-AE-1..4), which reads the frozen kappa_AE record
#     (runs/frozen/kappa_ae.json; no runner solves kappa). Per replicate on the host: Z keyed to
#     Psi_N + 0.30 sigma_m,wor and to the g-formula target Psi_N^gf, and the three wrong-clustering se2
#     variants (each checked against survey's wrong design). Per population: census predictions and the
#     G6 limb (a) boundary distances.
#   --regime fixed --e3 coarse --m 400: the individual-CF leakage-neutral sensitivity: oracle,
#     Rc-AIPW-indcf and Rc-AIPW-indcf-coarse (C_psu coarsened to quartile bins for the learner only) on the
#     same unit folds and level-4 seeds (C-RC-4/5). Coarsening also removes predictive information.
#
# Options:
#   --mode dev|dry|prod (R/run_mode.R): dev = SYNTHETIC root (tests); dry = production root, production
#     population IDs, rep IDs 1000001..; prod = production, rep IDs 1..R. dry and prod refuse to run unless
#     every frozen input has its recorded SHA-256 (config/simulation/frozen-design.json). --pops a:b names the
#     population IDs; --reps a:b the rep IDs; --task i --block k an array block.
#   --rd-alphas a,b: the R-d rungs to run (default all four; production runs 0.25,0.30, because 0.20 and 0.35
#     stopped after their learner-rate check G5). Each rung's level-4 stream is keyed by its own arm, so the
#     rungs that run are unchanged by the ones that do not. --no-rd 1 drops the R-d rungs.
#   --u6 1: U6's limbs for Ra-AIPW: a global fold relabel leaves psi_CF and se2 bit-identical; a fresh
#     within-stratum partition (fold structure "psu-s2") changes psi_CF. Required.
#
# Behaviour:
#   An R-a fit that fails after one Nelder-Mead refit raises ra_fit_failed; its arm (and paired TMLE, or the
#     R-d rung) is recorded FAILED for the replicate -- counted, never dropped, never a failed required check.
#   r_hat = mean over folds of max(||Q_hat - Q0|| (max over a), ||g_hat - g0||) on the census, post-clip, for
#     every fitted arm (production: forests on each population's first replicate only).
#   U7's tolerance is 64 eps max(|psi_oracle|, 1). Populations: the cell's (H, J) block built directly, nesting
#     checked against a direct build of a sub-block (entity addressing; required).
#   T2 of Ra-AIPW, Ra-TMLE and every R-d rung by R/t2_refine.R: a fold fit that does not trigger keeps the
#     primary grid's value bit for bit; a triggered one climbs GL(8,3) -> GL(10,4) -> GL(12,6) until two rungs
#     agree to 5e-7; an unresolved one leaves T2 and T1 NA with T2_status "unresolved" (every rung in
#     t2_refinement.csv). U9's brute-force limb on a selected fold adds an identical-node code check (1e-12)
#     and, for a refined fold, uses CC(13,7). The AE arm keeps its closed form (a deterministic construction
#     inside the clip bound). No forest T1/T2 is computed.
#
# Default cell: E1 many-strata, H = 25, m_h = 2, J = 20, f = 0.10, m = 2H,
# V = 2, n0 = 20. Populations: prefix views (H, 20)
# of the many-strata master (800 x 20), master_id "MS",
# population ids 1..B; R sample draws each. Arms: oracle, Ra-AIPW (PSU-CF, V = 2) and
# Ra-TMLE (one pooled 2-column logistic fluctuation of Ra-AIPW's clipped held-out fits).
#
# Writes runs/out/dev-e2e-m<m>-<UTC stamp>/ (or under --out-root): results.csv, census.csv,
# checks.csv, timing.csv, manifest.json, SUMMARY.md. Exits nonzero if any REQUIRED
# check fails or errors: U1, U3, U4/U5 (folds), U7, U8, finite outputs, and U9
# (P0-grid limb and the brute-force P0 f_v limb on the selected fits, kept
# for R-a), and for Ra-TMLE the conformance rows (initial
# clipped, single pooled 2-column tilt, converged, Q* = independent refit to 1e-8),
# its score equations and U8. The Gram floor lambda_min >= c_H is a PREMISE,
# recorded every replicate with required = FALSE. Informational checks
# are recorded with required = FALSE.
# as REQUIRED checks: the P0 grid hashes equal grids/MANIFEST.json for
# this H; the ambient RNG is unchanged by each replicate and by the whole run;
# R-b's served predictions are identical() to the seed = NULL predictions.
# Timing is recorded by component with no overlap (timing.csv, manifest, SUMMARY).
# The master-prefix nesting check is REQUIRED. In the
# fixed-strata path, also REQUIRED: G3(a)/(a-i)/(c) throttle arithmetic per alpha (q_m equals
# the frozen table; q_m <= pool; integrality), G3(b) master sufficiency; per R-d arm and
# fold: drawn whole-PSU count = q_m, every drawn PSU in the training fold, whole PSUs, no drawn
# PSU scored (disjoint), fewer units than the full pool (no fallback to R-a), pool = m(V-1)/V;
# the arm shares Ra-AIPW's sample, fold vector, training pools and scored units; and R-d's
# U3, U4/U5, U7, U8, U9 brute-force limb (same selected folds as R-a) and finite outputs.
# G3(a-ii)'s design margin and "subsamples differ across alpha" are informational.

via_wrapper <- endsWith(getOption("sim.engine.script", ""), "/runs/run_e2e_dev.R")
if (!via_wrapper)
  stop("run through the locked wrapper: Rscript simulation/engine/env/run.R ",
       "simulation/engine/runs/run_e2e_dev.R [--regime many|fixed] [--H h | --m m | --e2c fpc|strict|wr ... | --e3 panel|vsens|ctrl|coarse ...] [--B n] [--R n] [--u9-bf k] [--u6 1] [--u10 1] [--oracle-only 1 --calib 1] [--out-root dir]")
now <- function() proc.time()[["elapsed"]]
t_start <- now()
utc_start <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
stamp <- format(Sys.time(), "%Y%m%dT%H%M%SZ", tz = "UTC")
# Ambient RNG at start: every stream below is scoped, so it must be untouched at the end
# (a REQUIRED check; also checked around every replicate).
genv <- globalenv()
ambient_rng <- function() {
  kind <- RNGkind()
  list(kind = kind, seed = if (exists(".Random.seed", envir = genv, inherits = FALSE)) get(".Random.seed", envir = genv) else NULL)
}
rng_start <- ambient_rng()
rng_kind_start <- rng_start$kind
# Run time by component, with no overlap: each second of the run is counted once.
COMPONENTS <- c("p0_context", "u9_grid", "populations", "sampling", "ra_fits", "rb_forests", "rc_forests", "rd_fits",
                "remainders", "u9_bruteforce", "targeting", "rhat")
TIME <- new.env(parent = emptyenv())
for (k in COMPONENTS) TIME[[k]] <- 0
timed <- function(component, expr) {
  t0 <- now()
  on.exit(TIME[[component]] <- TIME[[component]] + (now() - t0))
  expr
}
# Peak resident memory of this R process so far (Linux /proc VmHWM, kB; NA elsewhere).
peak_rss_kb_now <- function() {
  f <- "/proc/self/status"
  l <- if (file.exists(f)) grep("^VmHWM:", readLines(f, warn = FALSE), value = TRUE) else character(0)
  if (length(l)) as.numeric(gsub("[^0-9]", "", l[1])) else NA_real_
}

engine_dir <- getOption("sim.engine.dir")
repo_root <- normalizePath(file.path(engine_dir, "..", ".."), winslash = "/", mustWork = TRUE)
args <- getOption("sim.engine.args")
arg_of <- function(name, default) {
  i <- which(args == name)
  if (length(i) && length(args) > i[1]) args[i[1] + 1L] else default
}
arg_int <- function(name, default) {
  x <- suppressWarnings(as.integer(arg_of(name, default)))
  if (length(x) != 1L || is.na(x) || x < 0L) stop(name, " must be a non-negative integer", call. = FALSE)
  x
}

ENGINE_FILES <- c("quadrature.R", "nuisance_truth.R", "p0_grid.R", "remainders.R", "rng_streams.R",
                  "folds.R", "dgp.R", "census.R", "sampling.R", "variance.R", "learner_ra.R",
                  "estimators.R", "tmle.R", "quad_rules_io.R", "learner_rb.R", "learner_rd.R", "controls.R", "gates.R",
                  "run_mode.R", "t2_refine.R")
for (f in ENGINE_FILES) source(file.path(engine_dir, "R", f), local = TRUE)

# ---- configuration --------------------------------------------------------------
A_CPSU <- 1.25                                      # S11, frozen ; explicit
NO_RD <- arg_int("--no-rd", 0L) == 1L                # drop the R-d rungs
# The R-d rungs to run (production: 0.25 and 0.30; 0.20 and 0.35 stopped after G5). Default all four.
RD_RUN <- {
  x <- suppressWarnings(as.numeric(strsplit(arg_of("--rd-alphas", paste(RD_ALPHAS, collapse = ",")), ",", fixed = TRUE)[[1]]))
  if (!length(x) || anyNA(x) || !all(vapply(x, function(a) any(abs(a - RD_ALPHAS) < 1e-12), TRUE)) || anyDuplicated(x))
    stop("--rd-alphas must be a comma-separated subset of ", paste(RD_ALPHAS, collapse = ","), call. = FALSE)
  RD_ALPHAS[vapply(RD_ALPHAS, function(a) any(abs(a - x) < 1e-12), TRUE)]           # frozen order
}
T2P <- T2_REFINE_POLICY                             # the refinement ladder; rung 0 is the primary GL(6,2)
T2_CACHE <- new.env(parent = emptyenv())            # refined rungs' tensors and truth, shared by every fit of the run
t2_diag <- list()                                   # per-fold refinement diagnostics (t2_refinement.csv)
REGIME <- arg_of("--regime", "many")
if (!REGIME %in% c("many", "fixed")) stop("--regime must be 'many' or 'fixed'", call. = FALSE)
E2C <- arg_of("--e2c", "none")                     # E2 construction (cell field e2c)
if (!E2C %in% c("none", "fpc", "strict", "wr"))
  stop("--e2c must be 'none', 'fpc', 'strict' or 'wr' (vanishing-fraction and unequal-allocation: runs/run_e2_census.R)", call. = FALSE)
if (E2C != "none" && REGIME != "many") stop("--e2c is a many-strata construction (experiment E2)", call. = FALSE)
E3 <- arg_of("--e3", "none")                       # E3 cell
if (!E3 %in% c("none", "panel", "vsens", "ctrl", "coarse")) stop("--e3 must be 'none', 'panel', 'vsens', 'ctrl' or 'coarse'", call. = FALSE)
if (E3 != "none" && REGIME != "fixed") stop("--e3 cells are fixed-strata", call. = FALSE)
ORACLE_ONLY <- arg_int("--oracle-only", 0L) == 1L
CALIB <- arg_int("--calib", 0L) == 1L
if (CALIB && !ORACLE_ONLY) stop("--calib 1 needs --oracle-only 1 (the calibration cell is oracle-only)", call. = FALSE)
U6 <- arg_int("--u6", 0L) == 1L
U10 <- arg_int("--u10", 0L) == 1L                    # the U10 cell
if (U10 && (REGIME != "fixed" || E3 != "none" || !ORACLE_ONLY || CALIB))
  stop("--u10 1 is the fixed-strata oracle-only U10 cell: --regime fixed --m 400 --oracle-only 1 --u10 1", call. = FALSE)
WR <- E2C == "wr"
if (WR && !ORACLE_ONLY) stop("the E2 WR code-check is oracle only: add --oracle-only 1", call. = FALSE)
if (E2C == "strict") {
  # E2 genuinely-strict trajectory, fitted cell: the first H
  # strata and the first 4 PSUs of the many-strata master [forced], f = 0.5; the fitted
  # R-a/R-b AIPW arm runs at m = 200 only. The census limb is runs/run_e2_census.R.
  MASTER_ID <- "MS"
  MASTER_H <- 800L; MASTER_J <- 20L
  H <- arg_int("--H", 100L)
  if (H != 100L) stop("the strict trajectory's fitted cell is m = 200 (H = 100); other rungs are census-only", call. = FALSE)
  M_H <- 2L; J <- 4L
} else if (E2C == "fpc") {
  # E2 FPC sweep: H = 50, m_h = 2, J in {100, 20, 8, 4} => f = m_h/J
  # in {0.02, 0.10, 0.25, 0.50}; prefixes (50, J) of the FPC-sweep master "FPC" (master
  # IDs and dimensions as in R/rng_streams.R). Only its 200 x 100 block is built: the master is 200 x 800,
  # and entity addressing makes every block's prefixes identical to the full master's.
  MASTER_ID <- "FPC"
  MASTER_H <- 200L; MASTER_J <- 100L; MASTER_J_FULL <- 800L
  H <- 50L; M_H <- 2L
  J <- arg_int("--J", 4L)
  if (!J %in% c(100L, 20L, 8L, 4L)) stop("--J must be 100, 20, 8 or 4 (the frozen FPC sweep of experiment E2)", call. = FALSE)
} else if (WR) {
  # E2 WR code-check: H = 50, m_h = 2, J = 4, f = 0.5, m = 100, prefix of
  # the FPC-sweep master; with-replacement first stage, selection-keyed stage 2.
  MASTER_ID <- "FPC"
  MASTER_H <- 200L; MASTER_J <- 100L; MASTER_J_FULL <- 800L            # the 200 x 100 block
  H <- 50L; M_H <- 2L; J <- 4L
} else if (REGIME == "many") {
  MASTER_ID <- "MS"
  MASTER_H <- 800L; MASTER_J <- 20L                 # many-strata master minimum,
  H <- arg_int("--H", 25L)                          # E1 many-strata rung (table)
  if (!H %in% c(25L, 50L, 100L, 200L, 400L, 800L))
    stop("--H must be one of 25, 50, 100, 200, 400, 800 (E1 many-strata rungs; m = 2H; H = 800 is the cheap-arm rung)", call. = FALSE)
  J <- 20L; M_H <- 2L
} else {
  # E1 fixed-strata (table): H = 5, m_h = m/5, J = 10 m_h; master "FS"
  # 5 x 6400 (RNG draft master IDs). DEV rungs: R-d's frozen grid starts at
  # m = 200 and every TMLE arm stops at m = 800, so m in {200, 400, 800}.
  MASTER_ID <- "FS"
  MASTER_H <- 5L; MASTER_J <- 6400L
  M_CELL <- arg_int("--m", 200L)
  if (!M_CELL %in% c(50L, 100L, 200L, 400L, 800L, 1600L, 3200L))
    stop("--m must be 50, 100, 200, 400, 800, 1600 or 3200 (the E1 fixed-strata ladder; R-d from 200; TMLE to 800)", call. = FALSE)
  if ((E3 %in% c("panel", "ctrl", "coarse") || U10) && M_CELL != 400L)
    stop("the E3 panel, its power controls, the leakage-neutral sensitivity and the U10 cell are m = 400", call. = FALSE)
  if (E3 == "vsens" && M_CELL != 200L) stop("the E3 V-sensitivity cell is m = 200", call. = FALSE)
  H <- 5L; M_H <- M_CELL %/% H; J <- 10L * M_H
}
N0 <- 20L                                           # (n0 = 20 [forced])
V <- arg_int("--V", 2L)                             # V = 2 fixed...
if (E3 == "vsens") {
  if (!V %in% c(2L, 5L)) stop("the V-sensitivity cell runs V = 2 or V = 5", call. = FALSE)
} else if (V != 2L) stop("V = 2 is fixed; only the E3 V-sensitivity cell runs V = 5", call. = FALSE)
M <- H * M_H
gcd2 <- function(a, b) if (b == 0L) a else gcd2(b, a %% b)
F_STR <- sprintf("%d/%d", M_H %/% gcd2(M_H, J), J %/% gcd2(M_H, J))   # f = m_h/J in lowest terms
# Cell-key construction label ('s E2-construction field): the U10 cell is its own cell at the
# panel's coordinates, so its sample stream is kept apart from the E1/E3 m = 400 cell.
CELL <- list(regime = REGIME, H = H, m_h = M_H, J = J, m = M, f = F_STR, e2c = if (U10) "u10" else E2C)
ARMS <- if (ORACLE_ONLY) { if (U10) c(ARM_ORACLE, ARM_ORACLE_U10) else ARM_ORACLE } else if (E2C == "fpc") {
  # oracle at all four fractions; fitted R-a & R-b PSU-CF AIPW at f = 0.50 only.
  if (J == 4L) c(ARM_ORACLE, ARM_RA_AIPW, ARM_RB_AIPW) else ARM_ORACLE
} else if (E2C == "strict") c(ARM_ORACLE, ARM_RA_AIPW, ARM_RB_AIPW) else
  if (E3 == "panel") c(ARM_ORACLE, ARM_RA_AIPW, ARM_RA_NOCF, ARM_RC_PSUCF, ARM_RC_NOCF, ARM_RC_INDCF, ARM_RC_INDCF_COARSE) else
  if (E3 == "ctrl") c(ARM_ORACLE, ARM_RA_AIPW, ARM_AE_AIPW) else
  if (E3 == "coarse") c(ARM_ORACLE, ARM_RC_INDCF, ARM_RC_INDCF_COARSE) else
  if (E3 == "vsens") c(ARM_ORACLE, ARM_RB_AIPW) else if (REGIME == "many") {
    if (H == 800L) c(ARM_ORACLE, ARM_RA_AIPW) else c(ARM_ORACLE, ARM_RA_AIPW, ARM_RA_TMLE, ARM_RB_AIPW, ARM_RB_TMLE)
  } else c(ARM_ORACLE,                                         # at m = 3200 the oracle is the internal psi_oracle
           if (M <= 1600L) ARM_RA_AIPW, if (M <= 800L) ARM_RA_TMLE,
           if (M >= 200L && !NO_RD) vapply(RD_RUN, rd_arm, ""))
HAS_RA <- ARM_RA_AIPW %in% ARMS; HAS_RA_TMLE <- ARM_RA_TMLE %in% ARMS
HAS_RB <- ARM_RB_AIPW %in% ARMS; HAS_RB_TMLE <- ARM_RB_TMLE %in% ARMS
HAS_RD <- any(startsWith(ARMS, "Rd-AIPW"))
HAS_RA_NOCF <- ARM_RA_NOCF %in% ARMS
HAS_RC_PSU <- ARM_RC_PSUCF %in% ARMS; HAS_RC_NOCF <- ARM_RC_NOCF %in% ARMS; HAS_RC_IND <- ARM_RC_INDCF %in% ARMS
HAS_RC <- HAS_RC_PSU || HAS_RC_NOCF || HAS_RC_IND
HAS_AE <- ARM_AE_AIPW %in% ARMS; HAS_RC_COARSE <- ARM_RC_INDCF_COARSE %in% ARMS
HAS_FITTED <- HAS_RA || HAS_RB || HAS_RD || HAS_RA_NOCF || HAS_RC
if (U6 && !HAS_RA) stop("--u6 1 needs the Ra-AIPW arm", call. = FALSE)
if (NO_RD && (REGIME != "fixed" || E3 != "none")) stop("--no-rd 1 is for the E1 fixed-strata cells", call. = FALSE)
PATH_LABEL <- paste0(if (WR) "E2 WR code-check (many-strata, f = 1/2; WR first stage, selection-keyed stage 2)" else
  if (U10) "G8 U10 cell (fixed-strata oracle; stage-2 weights omitted)" else
  if (E3 == "ctrl") "E3 power-control panel (fixed-strata; host R-a PSU-CF AIPW)" else
  if (E3 == "coarse") "E3 individual-CF leakage-neutral sensitivity (fixed-strata)" else
  if (E2C == "fpc") sprintf("E2 FPC sweep (many-strata), f = %s", F_STR) else
  if (E2C == "strict") "E2 genuinely-strict trajectory, fitted cell (many-strata, f = 1/2)" else
  if (E3 == "panel") "E3 R-c foil + R-a controls panel (fixed-strata)" else
  if (E3 == "vsens") sprintf("E3 V-sensitivity (fixed-strata), V = %d", V) else
  if (REGIME == "many") "E1 many-strata" else "E1 fixed-strata",
  if (ORACLE_ONLY && !CALIB && !U10 && !WR) " -- oracle (design-conditional anchor cell)" else "",
  if (CALIB && WR) " -- DEV calibration (R14: linearised known-N WR gate, U2 WR form)" else
  if (CALIB) " -- DEV calibration cell (U2, U8 limbs)" else "")

RM <- run_mode_setup(arg_of, repo_root, B_default = 2L, R_default = 10L, master_id = MASTER_ID)
MODE <- RM$mode; ROOT <- RM$root; POP_IDS <- RM$pops; REP_IDS <- RM$reps; B <- RM$B; R <- RM$R
RUN_LABEL <- RM$label
U9_BF <- arg_int("--u9-bf", if (MODE == "prod") 0L else 4L)     # U9 is a dry-run test: none in production
if (B < 1L || R < 1L) stop("need at least one population and one replicate", call. = FALSE)
if ((CALIB || U10) && R < 2L) stop("the calibration and U10 cells need R >= 2 (their gates are over redraws)", call. = FALSE)
# the r_hat for the forest arms: every replicate in dev and dry; in production each population's first
# replicate only (a display; R-b's rate gate and R-c's audit are dry-run objects).
RHAT_FOREST_ALL <- MODE != "prod"
out_root <- arg_of("--out-root", file.path(engine_dir, "runs", "out"))
RUN_TAG <- if (WR) sprintf("wr%s", if (CALIB) "-calib" else "") else if (U10) sprintf("u10-m%d", M) else
  if (E3 == "ctrl") sprintf("e3ctrl-m%d", M) else if (E3 == "coarse") sprintf("e3coarse-m%d", M) else
  if (E2C == "fpc") sprintf("fpc-J%d%s", J, if (CALIB) "-calib" else "") else
  if (E2C == "strict") sprintf("strict-m%d", M) else if (E3 == "panel") sprintf("e3panel-m%d", M) else
  if (E3 == "vsens") sprintf("e3vsens-V%d", V) else sprintf(if (REGIME == "many") "m%d" else "fs-m%d", M)
if (ORACLE_ONLY && !CALIB && !U10 && !WR) RUN_TAG <- paste0(RUN_TAG, "-oracle")
if (RM$task > 0L) RUN_TAG <- sprintf("%s-t%03d", RUN_TAG, RM$task)
out_dir <- file.path(out_root, sprintf("%s-e2e-%s-%s", MODE, RUN_TAG, stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Saved base quadrature rules, sha256-verified against grids/MANIFEST.json; required:
# a rule missing from the artifact is an error, never a recompute.
QR <- quad_rules_dev(engine_dir)
quad_rules_activate(QR, required = TRUE)
# The refinement ladder's finer base rules, a separately hashed artifact merged after the primary one
# (keys the frozen artifact lacks only; a shared key must agree bit for bit).
QRR <- t2_refine_rules(engine_dir)
t2_refine_rules_extend(QRR)

co <- dgp_coefficients(A_CPSU)
a_term <- co$a_cpsu != 0                            # the correct R-a model has A:Cpsu iff a != 0
learner <- function(train) fit_ra(train, a_cpsu_term = a_term)
cat(sprintf("[dev-e2e] %s  mode %s  cell %s e2c=%s e3=%s H=%d m_h=%d J=%d f=%s m=%d V=%d  arms %s  B=%d (pops %s) R=%d (reps %s)%s%s  a_cpsu=%.2f  out=%s\n",
            utc_start, MODE, REGIME, E2C, E3, H, M_H, J, F_STR, M, V, paste(ARMS, collapse = ","), B,
            paste(range(POP_IDS), collapse = ".."), R, paste(range(REP_IDS), collapse = ".."),
            if (U6) "  U6" else "", if (CALIB) "  CALIB" else "", A_CPSU, out_dir))

checks <- list(); results <- list(); census <- list(); census_strata <- list(); timing <- list()
add_checks <- function(df, pop_id, rep_id, arm, required) {
  if (is.null(df) || !nrow(df)) return(invisible())
  checks[[length(checks) + 1L]] <<- cbind(
    data.frame(pop_id = pop_id, rep_id = rep_id, arm = arm, stringsAsFactors = FALSE),
    df[, c("check", "value", "reference", "abs_diff", "tolerance", "pass")],
    required = required)
}
add_error <- function(what, e, pop_id, rep_id, arm) {
  cat("  ERROR in", what, ":", conditionMessage(e), "\n")
  add_checks(data.frame(check = paste(what, "completed without error --", conditionMessage(e)),
                        value = NA_real_, reference = NA_real_, abs_diff = NA_real_,
                        tolerance = NA_real_, pass = FALSE), pop_id, rep_id, arm, TRUE)
}

# Row-bind data frames whose columns differ (an arm that failed adds no row); missing -> NA.
rbind_fill <- function(dfs) {
  dfs <- Filter(Negate(is.null), dfs)
  if (!length(dfs)) return(data.frame())
  cols <- unique(unlist(lapply(dfs, names)))
  do.call(rbind, lapply(dfs, function(d) { for (cl in setdiff(cols, names(d))) d[[cl]] <- NA; d[cols] }))
}

# ---- P0 grid context: one per run (H-specific) ---------------------------
ctx <- timed("p0_context", p0_context(H, co))
Psi_grid <- sum(ctx$primary$w * (ctx$primary_truth$Q1 - ctx$primary_truth$Q0))   # Psi of P_{0,H},
# Grid identity: the hashes recorded for this H in grids/MANIFEST.json (required).
gh <- QR$manifest$expected_grid_hashes
gi <- if (is.data.frame(gh)) which(gh$H == H) else integer(0)
exp_hash <- if (length(gi) == 1L) c(gh$primary_sha256[gi], gh$second_sha256[gi]) else c(NA_character_, NA_character_)
add_checks(data.frame(check = sprintf("P0 grid hash (%s, H = %d) = grids/MANIFEST.json record (%s)",
                                      c("primary", "second"), H, ifelse(is.na(exp_hash), "NO RECORD", substr(exp_hash, 1, 12))),
                      value = NA_real_, reference = NA_real_, abs_diff = NA_real_, tolerance = NA_real_,
                      pass = c(identical(ctx$primary_hash, exp_hash[1]), identical(ctx$second_hash, exp_hash[2]))),
           NA, NA, "all", TRUE)
tryCatch({
  u9a <- timed("u9_grid", u9_psi(ctx))
  add_checks(transform(u9a, check = paste("U9 P0-grid limb:", check)), NA, NA, "all", TRUE)
}, error = function(e) add_error("U9 P0-grid limb", e, NA, NA, "all"))
cat(sprintf("[dev-e2e] P0 context: %d + %d grid points (hashes %s / %s), Psi(P0,H=%d) = %.10f, %.1f s; U9 grid limb %.1f s\n",
            nrow(ctx$primary), nrow(ctx$second), substr(ctx$primary_hash, 1, 12), substr(ctx$second_hash, 1, 12),
            H, Psi_grid, TIME$p0_context, TIME$u9_grid))

# U9 brute-force selection (u9_select's recorded convention): global
# replicate index g = (b - 1) * R + r, folds cycling. In the fixed-strata path the same
# (replicate, fold) pairs are also cross-checked for every R-d arm (A1 keeps U9's
# brute-force limb for R-d), so a run carries U9_BF x (1 + 4) brute-force checks.
u9_sel <- u9_select(B * R, V, U9_BF)

# G3 by arithmetic, fixed-strata path: the throttle per alpha and the master.
if (REGIME == "fixed") {
  if (HAS_RD) {
    g3 <- rd_g3_rows(M, M_H, V)
    add_checks(g3[g3$required, ], NA, NA, "all", TRUE)
    add_checks(g3[!g3$required, ], NA, NA, "all", FALSE)
  }
  add_checks(data.frame(check = sprintf("G3(b) master sufficiency: J = %d <= %d PSUs per stratum and H = %d <= %d strata of master %s",
                                        J, MASTER_J, H, MASTER_H, MASTER_ID),
                        value = J, reference = MASTER_J, abs_diff = NA_real_, tolerance = NA_real_,
                        pass = J <= MASTER_J && H <= MASTER_H), NA, NA, "all", TRUE)
}

# G3(b) by arithmetic, every other path: the cell fits inside its master (FPC is 200 x 800;
# the cell builds only its own (H, J) block, identical to the master's prefix by entity addressing).
if (REGIME == "many") {
  MJ <- if (E2C %in% c("fpc", "wr")) MASTER_J_FULL else MASTER_J
  add_checks(data.frame(check = sprintf("G3(b) master sufficiency: J = %d <= %d PSUs per stratum and H = %d <= %d strata of master %s",
                                        J, MJ, H, MASTER_H, MASTER_ID),
                        value = J, reference = MJ, abs_diff = NA_real_, tolerance = NA_real_,
                        pass = J <= MJ && H <= MASTER_H), NA, NA, "all", TRUE)
}

# ---- one replicate --------------------------------------------------------------
# r_hat: per fold fit, max(||Q_hat - Q0|| (max over a), ||g_hat - g0||) on the census, post-clip; equal
# weights over folds.
rhat_of <- function(fits, units) {
  r <- t(vapply(fits, function(f) nuisance_l2(f, units), numeric(5)))
  c(r_Q = mean(r[, "Q"]), r_g = mean(r[, "g"]), r_hat = mean(r[, "r"]))
}
# The replicate's fold structure as psu_cf_aipw() builds it: every arm shares it (R-d's shared checks use it, so
# they hold even where R-a FAILED or does not run, m = 3200).
fold_reference <- function(fold) {
  rows <- fold_rows(fold)
  list(fold = fold, fit_rows = lapply(seq_along(rows), function(v) sort(unlist(rows[-v], use.names = FALSE))),
       score_rows = rows)
}

run_rep <- function(b, r, bi, ri, pop, frame, cen) {
  t_rep <- now()
  g_idx <- (bi - 1L) * R + ri                                        # positions (b, r are the IDs)
  key <- list(cell = CELL, population_id = b, rep_id = r)
  s <- timed("sampling", if (WR) draw_sample_wr(frame, M_H, N0, key, ROOT) else
    draw_sample(frame, M_H, N0, key, ROOT))                       # level-2 stream,
  u1 <- if (WR) check_u1_wr(s, pop) else check_u1(s, pop)
  CL <- if (WR) s$j else s$k                                       # ultimate cluster: WR selection or PSU
  add_checks(transform(u1$checks, value = as.numeric(pass), reference = 1, abs_diff = NA_real_,
                       tolerance = NA_real_, check = paste("U1", check)), b, r, "all", TRUE)

  or <- timed("sampling", oracle_aipw(s, CL))
  u8o <- check_u8(s$phi0, s$w, s$h, CL)
  # The two-stage WOR SE has no meaning under a WR first stage (the WR cell declares no fpc).
  wor_o <- if (WR) list(se2_wor = NA_real_) else se2_two_stage_wor(s$phi0, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
  u8_rows <- function(u) data.frame(
    check = c("U8 se2 = survey::svymean, fpc not declared (rel 1e-12)",
              "U8 psi = survey::svymean (rel 1e-12)"),
    value = c(u$se2, u$psi), reference = c(u$se2_survey, u$psi_survey),
    abs_diff = c(abs(u$se2 - u$se2_survey), abs(u$psi - u$psi_survey)),
    tolerance = c(u$tol * abs(u$se2_survey), u$tol * abs(u$psi_survey)),
    pass = c(is.finite(u$rel_diff) && u$rel_diff <= u$tol, is.finite(u$psi_rel_diff) && u$psi_rel_diff <= u$tol))
  add_checks(u8_rows(u8o), b, r, ARM_ORACLE, TRUE)

  sel <- u9_sel$fold[u9_sel$rep == g_idx]
  sig_m <- sqrt(cen$sigma2_m); sig_wor <- sqrt(cen$sigma2_m_wor)
  m_s <- sum(attr(s, "design")$m_h)
  # One results row per arm. se2_two_stage_wor is the per-sample two-stage WOR estimate with
  # both stages' fpc declared (U8's two-stage limb, a dry-run object over redraws, not gated);
  # R_var_wor and se2_wor_over_se2 put it beside the UC/no-FPC se2.
  row <- function(arm, est, se2_wor, u8, extra) {
    Z <- (est$psi - cen$Psi_N) / est$se
    data.frame(pop_id = b, rep_id = r, arm = arm, psi_hat = est$psi, se_hat = est$se,
               Psi_N = cen$Psi_N, Psi = Psi_grid, sigma_m = sig_m, sigma_m_wor = sig_wor,
               Z = Z, cover95 = abs(Z) <= stats::qnorm(0.975),
               Z_wor = (est$psi - cen$Psi_N) / sig_wor, R_var = est$se2 / cen$sigma2_m,
               se2_two_stage_wor = se2_wor, R_var_wor = se2_wor / cen$sigma2_m_wor,
               se2_wor_over_se2 = se2_wor / est$se2, N_hat_over_N_minus_1 = sum(s$w) / cen$N - 1,
               n_units = nrow(s), n_psu = sum(attr(s, "design")$m_h),
               U1_pass = u1$pass, U8_rel_diff = u8$rel_diff, U8_pass = u8$pass, extra,
               arm_status = "OK", label = RUN_LABEL, stringsAsFactors = FALSE)
  }
  na_extra <- data.frame(T2 = NA_real_, T1 = NA_real_, T2_over_sigma_wor = NA_real_,
                         T1_over_sigma_wor = NA_real_, omega_1 = NA_real_, omega_2 = NA_real_,
                         U3_pass = NA, U4_pass = NA, U5_pass = NA, U7_lhs = NA_real_, U7_rhs = NA_real_,
                         U7_abs_diff = NA_real_, U7_tol = NA_real_, U7_pass = NA,
                         n_fits = NA_integer_, n_singular = NA_integer_, n_not_clean = NA_integer_,
                         n_warnings = NA_integer_, clip_hits = NA_integer_, fit_seconds = NA_real_)
  na_tmle <- data.frame(eps_1 = NA_real_, eps_0 = NA_real_, eps_norm = NA_real_, m14_eps_norm = NA_real_,
                        lambda_min = NA_real_, gram_premise_met = NA, conformance_max_abs = NA_real_,
                        score_max_abs = NA_real_, tmle_iter = NA_integer_, tmle_seconds = NA_real_,
                        paired_tmle_minus_aipw = NA_real_, sqrt_m_abs_paired = NA_real_,
                        se2_diff_over_sigma2_m = NA_real_, m14_abs_se2_diff_over_sigma2_m = NA_real_)
  rows <- list(row(ARM_ORACLE, or, wor_o$se2_wor, u8o, cbind(na_extra, na_tmle)))
  # A FAILED arm keeps a results row (psi NA) and an informational check row; it is counted, never dropped.
  fail_row <- function(arm, msg) {
    rows[[length(rows) + 1L]] <<- data.frame(pop_id = b, rep_id = r, arm = arm, Psi_N = cen$Psi_N, arm_status = "FAILED",
                                             failure = msg, label = RUN_LABEL, stringsAsFactors = FALSE)
    add_checks(data.frame(check = paste("R4 FAILED arm (counted, not dropped):", msg), value = NA_real_, reference = NA_real_,
                          abs_diff = NA_real_, tolerance = NA_real_, pass = TRUE), b, r, arm, FALSE)
  }
  rhat_cols <- function(rw, fits, forest = FALSE) {
    if (forest && !RHAT_FOREST_ALL && ri != 1L) return(rw)
    rh <- timed("rhat", rhat_of(fits, pop$units))
    rw$r_Q <- rh[["r_Q"]]; rw$r_g <- rh[["r_g"]]; rw$r_hat <- rh[["r_hat"]]
    rw
  }
  if (CALIB && WR) {                              # WR calibration, per draw: linearised known-N WR score
    lv <- lin_var_wr(s$phi0, cen$Psi_N, cen$N, s$w, s$h, s$j)
    rows[[1]]$lin_total <- lv$total; rows[[1]]$lin_se2_wr <- lv$se2
    rows[[1]]$R_var_wr <- or$se2 / cen$sigma2_m_wr
    ok_lv <- is.finite(lv$se2) && abs(lv$se2 - lv$se2_direct) <= 1e-12 * abs(lv$se2_direct) &&
      abs(lv$total - lv$total_direct) <= 1e-12 * max(abs(lv$total_direct), 1e-300) + 1e-18
    add_checks(data.frame(check = "WR linearised score: survey svytotal (ids = ~selection) = direct cluster-total formula (rel 1e-12)",
                          value = lv$se2, reference = lv$se2_direct, abs_diff = abs(lv$se2 - lv$se2_direct),
                          tolerance = 1e-12 * abs(lv$se2_direct), pass = ok_lv), b, r, ARM_ORACLE, TRUE)
  }
  if (CALIB && !WR) {                             # calibration, per draw (dry-run objects)
    f1 <- se2_single_stage_fpc(s$phi0, s$w, s$h, s$k, s$M_pop_h)            # Hajek, reported
    rows[[1]]$se2_single_stage_fpc <- f1$se2_fpc1
    lv <- lin_var_designs(s$phi0, cen$Psi_N, cen$N, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)   # gated
    rows[[1]]$lin_total <- lv$total
    rows[[1]]$lin_se2_uc <- lv$se2_uc; rows[[1]]$lin_se2_wor <- lv$se2_wor; rows[[1]]$lin_se2_fpc1 <- lv$se2_fpc1
  }
  add_checks(check_finite("oracle psi, se, Z", list(or$psi, or$se, rows[[1]]$Z, or$phi)), b, r, ARM_ORACLE, TRUE)
  line <- sprintf("oracle %.5f (se %.5f, se_wor %.5f)", or$psi, or$se, sqrt(wor_o$se2_wor))
  more <- character(0)
  if (U10) {                                      # U10 (G8): stage-2 weights omitted
    ue <- u10_estimate(s)
    u8u <- check_u8(s$phi0, ue$w1, s$h, s$k)
    add_checks(u8_rows(u8u), b, r, ARM_ORACLE_U10, TRUE)
    rows <- c(rows, list(row(ARM_ORACLE_U10, ue, NA_real_, u8u, cbind(na_extra, na_tmle))))
    add_checks(check_finite("oracle-U10 psi, se, Z", list(ue$psi, ue$se, rows[[length(rows)]]$Z)), b, r, ARM_ORACLE_U10, TRUE)
    line <- sprintf("%s  U10 %.5f (se %.5f, Z %+.2f)", line, ue$psi, ue$se, rows[[length(rows)]]$Z)
  }

  # Arms are switched by the path's flags (HAS_*); every block below is unchanged in what it
  # computes and checks when its arm runs. Oracle-only cells (E2 f < 0.5) draw no folds.
  if (HAS_FITTED) fold <- timed("sampling", dev_psu_folds(s, key, V, ROOT))      # level-3 stream,
  cf <- NULL; ra_failed <- NULL
  fref <- if (HAS_FITTED) fold_reference(fold) else NULL
  if (HAS_RA) {
    t_fit <- now()
    cf <- tryCatch(timed("ra_fits", psu_cf_aipw(s, fold, learner)),          # mixed-model fits + scoring
                   ra_fit_failed = function(e) { ra_failed <<- conditionMessage(e); NULL })
    fit_sec <- now() - t_fit
    if (!is.null(ra_failed)) {
      fail_row(ARM_RA_AIPW, ra_failed)
      if (HAS_RA_TMLE) fail_row(ARM_RA_TMLE, "paired with a FAILED Ra-AIPW (R4)")
      line <- sprintf("%s  Ra-AIPW FAILED (R4)", line)
    }
  }
  if (HAS_RA && !is.null(cf)) {
    u3 <- check_u3(s$w, cf$fold, cf$omega, cf$N_hat_v)
    u45 <- check_u4_u5(s, cf)
    u7 <- check_u7(s, cf$fold, cf$fits, cf$psi, or$psi)
    u8r <- check_u8(cf$phi, s$w, s$h, s$k)
    wor_r <- se2_two_stage_wor(cf$phi, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
    t2 <- timed("remainders", t2_refined(ctx, cf$fits, cf$omega, T2P, T2_CACHE))   # post-clip fits; refined T2
    t2c <- t2_refined_cols(t2, cf$psi, or$psi, sig_wor); T1 <- t2c$T1
    t2_diag[[length(t2_diag) + 1L]] <<- t2_refined_diag(t2, b, r, ARM_RA_AIPW)
    add_checks(u8_rows(u8r), b, r, ARM_RA_AIPW, TRUE)
    add_checks(u3, b, r, ARM_RA_AIPW, TRUE)
    add_checks(u45, b, r, ARM_RA_AIPW, TRUE)
    add_checks(u7, b, r, ARM_RA_AIPW, TRUE)
    # psi_CF as fold-weighted averages vs the pooled Hajek mean: equal in algebra (info).
    add_checks(data.frame(check = "psi_CF (fold-weighted) vs pooled Hajek mean of phi_hat (info, 64 eps)",
                          value = cf$psi, reference = cf$psi_pooled, abs_diff = abs(cf$psi - cf$psi_pooled),
                          tolerance = 64 * .Machine$double.eps * abs(cf$psi_pooled),
                          pass = abs(cf$psi - cf$psi_pooled) <= 64 * .Machine$double.eps * abs(cf$psi_pooled)),
               b, r, ARM_RA_AIPW, FALSE)
    if (length(sel)) add_checks(timed("u9_bruteforce", u9_cross_check_refined(ctx, cf$fits, t2, sel, T2P, T2_CACHE)), b, r, ARM_RA_AIPW, TRUE)
    d <- cf$diagnostics
    ra_extra <- data.frame(t2c, omega_1 = cf$omega[1], omega_2 = cf$omega[2],
                           U3_pass = all(u3$pass), U4_pass = u45$pass[1], U5_pass = all(u45$pass[2:3]),
                           U7_lhs = u7$value, U7_rhs = u7$reference, U7_abs_diff = u7$abs_diff,
                           U7_tol = u7$tolerance, U7_pass = u7$pass,
                           n_fits = nrow(d), n_singular = sum(d$singular), n_not_clean = sum(!d$clean),
                           n_warnings = sum(d$n_warnings), clip_hits = as.integer(cf$clip_hits),
                           fit_seconds = fit_sec)
    rows <- c(rows, list(rhat_cols(row(ARM_RA_AIPW, cf, wor_r$se2_wor, u8r, cbind(ra_extra, na_tmle)), cf$fits)))
    add_checks(t2_finite_check("Ra-AIPW psi, se, Z, T1, T2, phi_hat", t2, list(cf$psi, cf$se, rows[[length(rows)]]$Z, cf$phi), T1), b, r, ARM_RA_AIPW, TRUE)
    if (E3 == "ctrl") {
      # E3 power controls on the host: C1 shift Z (and the g-formula-keyed Z, reported) and the
      # three wrong-clustering se2 variants, each reproduced by survey's wrong design (required).
      ra_i <- length(rows)
      zc <- panel_shift_z(cf$psi, cf$se, cen$Psi_N, cen$Psi_N_gf, sig_wor)
      wc <- wrong_clustering_se2(cf$phi, s$w, s$h, s$k)
      sv <- wrong_clustering_survey(cf$phi, s$w, s$h, s$k)
      rel <- abs(wc - sv[names(wc)]) / abs(sv[names(wc)])
      add_checks(data.frame(check = c(sprintf("C2 wrong-clustering %s se2: direct = survey's wrong design (rel 1e-10)", names(wc)),
                                      "C2 wrong-clustering uc se2 = the host's UC se2 (rel 1e-12)"),
                            value = c(wc, wc[["uc"]]), reference = c(sv[names(wc)], cf$se2), abs_diff = NA_real_,
                            tolerance = c(rep(1e-10, 4), 1e-12),
                            pass = c(is.finite(rel) & rel <= 1e-10, abs(wc[["uc"]] - cf$se2) <= 1e-12 * cf$se2)), b, r, ARM_RA_AIPW, TRUE)
      rows[[ra_i]]$Z_shift <- zc[["Z_shift"]]; rows[[ra_i]]$Z_gf <- zc[["Z_gf"]]
      for (v in c("nomh", "unit", "nostrat")) rows[[ra_i]][[paste0("R_var_", v)]] <- wc[[v]] / cen$sigma2_m
      line <- sprintf("%s  Z_shift %+.2f  R_var nomh/unit/nostrat %.3f/%.3f/%.3f", line, zc[["Z_shift"]],
                      wc[["nomh"]] / cen$sigma2_m, wc[["unit"]] / cen$sigma2_m, wc[["nostrat"]] / cen$sigma2_m)
    }
    line <- sprintf("%s  AIPW %.5f (se %.5f)  T2 %+.2e  T1 %+.2e  fits %.1fs", line, cf$psi, cf$se, t2$T2, T1, fit_sec)

    if (HAS_RA_TMLE) {
      # Paired TMLE: targets Ra-AIPW's own clipped held-out nuisances (TC-T3).
      t_tm <- now()
      tm <- timed("targeting", psu_cf_tmle(s, cf, RA_Q_CLIP, RA_G_TRUNC))
      tmle_sec <- now() - t_tm
      cht <- check_tmle(s, cf, tm)
      u8t <- check_u8(tm$phi, s$w, s$h, s$k)
      wor_t <- se2_two_stage_wor(tm$phi, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
      # T2^TMLE on the targeted fold fits; the trigger is the AIPW fold's (the tilt is smooth: the kinks are the initial fit's)
      t2t <- timed("remainders", t2_refined(ctx, tm$targeted_fits, tm$omega, T2P, T2_CACHE, triggers = lapply(t2$folds, `[[`, "trigger")))
      t2tc <- t2_refined_cols(t2t, tm$psi, or$psi, sig_wor); T1t <- t2tc$T1
      t2_diag[[length(t2_diag) + 1L]] <<- t2_refined_diag(t2t, b, r, ARM_RA_TMLE)
      add_checks(u8_rows(u8t), b, r, ARM_RA_TMLE, TRUE)
      add_checks(cht[1:5, ], b, r, ARM_RA_TMLE, TRUE)                 # conformance + score: required
      add_checks(cht[6, ], b, r, ARM_RA_TMLE, FALSE)                  # one-step = plug-in (info)
      add_checks(cht[7, ], b, r, ARM_RA_TMLE, FALSE)                  # Gram premise (T2): recorded, not a build check
      eps_norm <- sqrt(sum(tm$eps^2))
      tm_tmle <- data.frame(eps_1 = tm$eps[1], eps_0 = tm$eps[2], eps_norm = eps_norm,
                            m14_eps_norm = m_s^0.25 * eps_norm, lambda_min = tm$lambda_min,
                            gram_premise_met = cht$pass[7], conformance_max_abs = cht$value[4],
                            score_max_abs = cht$value[5], tmle_iter = as.integer(tm$iter), tmle_seconds = tmle_sec,
                            paired_tmle_minus_aipw = tm$psi - cf$psi,
                            sqrt_m_abs_paired = sqrt(m_s) * abs(tm$psi - cf$psi),
                            se2_diff_over_sigma2_m = (tm$se2 - cf$se2) / cen$sigma2_m,
                            m14_abs_se2_diff_over_sigma2_m = m_s^0.25 * abs(tm$se2 - cf$se2) / cen$sigma2_m)
      tmle_extra <- na_extra
      for (cl in names(t2tc)) tmle_extra[[cl]] <- t2tc[[cl]]
      tmle_extra$omega_1 <- tm$omega[1]; tmle_extra$omega_2 <- tm$omega[2]
      tmle_extra$fit_seconds <- tmle_sec
      rows <- c(rows, list(row(ARM_RA_TMLE, tm, wor_t$se2_wor, u8t, cbind(tmle_extra, tm_tmle))))
      add_checks(t2_finite_check("Ra-TMLE psi, se, Z, eps, T1, T2, phi_star", t2t, list(tm$psi, tm$se, rows[[length(rows)]]$Z, tm$eps, tm$phi), T1t),
                 b, r, ARM_RA_TMLE, TRUE)
      line <- sprintf("%s  TMLE %.5f (se %.5f)  lmin %.3f", line, tm$psi, tm$se, tm$lambda_min)
    }

    if (U6) {
      # U6. Invariance limb (oracle and R-a only; R-a is a deterministic function of its
      # training set): a global relabel of fold indices leaves psi_CF and se2 bit-identical.
      # Sensitivity limb: a fresh within-stratum partition (fold structure "psu-s2", m_vh kept)
      # changes psi_CF. Both required. The oracle has no folds, so its invariance is vacuous.
      fold_rl <- relabel_folds(fold, rev(seq_len(V)))
      cf_rl <- timed("ra_fits", psu_cf_aipw(s, fold_rl, learner))
      fold_s2 <- timed("sampling", dev_psu_folds(s, key, V, ROOT, fold_structure_id = "psu-s2"))
      cf_s2 <- tryCatch(timed("ra_fits", psu_cf_aipw(s, fold_s2, learner)), ra_fit_failed = function(e) list(psi = NA_real_))
      new_part <- !identical(split(paste(s$h, s$k), as.vector(fold_s2)), split(paste(s$h, s$k), as.vector(fold)))
      u6 <- c(identical(cf_rl$psi, cf$psi), identical(cf_rl$se2, cf$se2), identical(cf_rl$phi, cf$phi),
              new_part, new_part && !identical(cf_s2$psi, cf$psi))
      add_checks(data.frame(check = c("U6 invariance limb (Ra-AIPW): global fold relabel leaves psi_CF bit-identical",
                                      "U6 invariance limb (Ra-AIPW): global fold relabel leaves se2 bit-identical",
                                      "U6 invariance limb (Ra-AIPW): global fold relabel leaves phi_hat bit-identical",
                                      "U6 sensitivity limb: the fresh partition (psu-s2) differs from the primary one",
                                      "U6 sensitivity limb (Ra-AIPW): re-assigning PSUs (psu-s2) changes psi_CF"),
                            value = c(cf_rl$psi - cf$psi, cf_rl$se2 - cf$se2, max(abs(cf_rl$phi - cf$phi)),
                                      as.numeric(new_part), cf_s2$psi - cf$psi),
                            reference = c(0, 0, 0, 1, NA_real_), abs_diff = NA_real_, tolerance = 0, pass = u6),
                 b, r, ARM_RA_AIPW, TRUE)
    }
  }

  # ---- engineered aligned-error control (R/controls.R C-AE-1..4): PSU-CF AIPW on the host's
  # sample and folds with the deterministic kappa_AE shift of eta0. T2 is a closed-form constant; the
  # generic T2 code must reproduce the solved map (required), and U9's brute-force limb applies.
  if (HAS_AE) tryCatch({
    t_ae <- now()
    cfa <- timed("remainders", psu_cf_aipw(s, fold, function(train) ae_fit(AE$kappa, co, RA_Q_CLIP, RA_G_TRUNC), arm = ARM_AE_AIPW))
    ae_sec <- now() - t_ae
    u3a <- check_u3(s$w, cfa$fold, cfa$omega, cfa$N_hat_v); u45a <- check_u4_u5(s, cfa)
    u7a <- check_u7(s, cfa$fold, cfa$fits, cfa$psi, or$psi); u8a <- check_u8(cfa$phi, s$w, s$h, s$k)
    wor_a <- se2_two_stage_wor(cfa$phi, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
    t2a <- timed("remainders", t2_closed_form(ctx, cfa$fits, cfa$omega))
    T1a <- t1_by_subtraction(cfa$psi, or$psi, t2a$T2)
    add_checks(u8_rows(u8a), b, r, ARM_AE_AIPW, TRUE); add_checks(u3a, b, r, ARM_AE_AIPW, TRUE)
    add_checks(u45a, b, r, ARM_AE_AIPW, TRUE); add_checks(u7a, b, r, ARM_AE_AIPW, TRUE)
    add_checks(data.frame(check = c("AE T2: the generic closed-form code reproduces the solved kappa_AE map (rel 1e-10; C-AE-3)",
                                    "AE clip bound: no scored unit's perturbed nuisance hits the fitted clip (C-AE-2)"),
                          value = c(t2a$T2, cfa$clip_hits), reference = c(AE$T2, 0), abs_diff = c(abs(t2a$T2 - AE$T2), NA),
                          tolerance = c(1e-10 * abs(AE$T2), 0),
                          pass = c(abs(t2a$T2 - AE$T2) <= 1e-10 * abs(AE$T2), cfa$clip_hits == 0L)), b, r, ARM_AE_AIPW, TRUE)
    if (length(sel)) {
      xca <- timed("u9_bruteforce", u9_cross_check(ctx, cfa$fits[sel]))
      add_checks(data.frame(check = paste0("U9 brute-force P0 f_v vs closed form (fold ", sel, ", 1e-6)"),
                            value = xca$closed_form, reference = xca$brute_force, abs_diff = xca$abs_diff,
                            tolerance = xca$tolerance, pass = xca$pass), b, r, ARM_AE_AIPW, TRUE)
    }
    ae_extra <- na_extra
    ae_extra$T2 <- t2a$T2; ae_extra$T1 <- T1a; ae_extra$T2_over_sigma_wor <- t2a$T2 / sig_wor
    ae_extra$T1_over_sigma_wor <- T1a / sig_wor; ae_extra$omega_1 <- cfa$omega[1]; ae_extra$omega_2 <- cfa$omega[2]
    ae_extra$U3_pass <- all(u3a$pass); ae_extra$U4_pass <- u45a$pass[1]; ae_extra$U5_pass <- all(u45a$pass[2:3])
    ae_extra$U7_lhs <- u7a$value; ae_extra$U7_rhs <- u7a$reference; ae_extra$U7_abs_diff <- u7a$abs_diff
    ae_extra$U7_tol <- u7a$tolerance; ae_extra$U7_pass <- u7a$pass; ae_extra$clip_hits <- as.integer(cfa$clip_hits)
    ae_extra$fit_seconds <- ae_sec
    rows <- c(rows, list(row(ARM_AE_AIPW, cfa, wor_a$se2_wor, u8a, cbind(ae_extra, na_tmle))))
    add_checks(check_finite("AE-AIPW psi, se, Z, T1, T2, phi_hat", list(cfa$psi, cfa$se, rows[[length(rows)]]$Z, T1a, t2a$T2, cfa$phi)),
               b, r, ARM_AE_AIPW, TRUE)
    more <- c(more, sprintf("AE-AIPW %.5f (se %.5f)  T2 %+.3e (= %.3f sigma_wor)  T1 %+.2e", cfa$psi, cfa$se, t2a$T2, t2a$T2 / sig_wor, T1a))
    rm(cfa)
  }, error = function(e) add_error(sprintf("AE arm (pop %d, rep %d)", b, r), e, b, r, ARM_AE_AIPW))

  # ---- R-a no-CF (E3 negative control): one full-sample R-a fit scores its own
  # training set; U3/U4/U5/U7 skipped (vacuous or false by construction), U8 applies.
  if (HAS_RA_NOCF) tryCatch(ra_fit_failed = function(e) fail_row(ARM_RA_NOCF, conditionMessage(e)), expr = {
    t_nc <- now()
    cn <- timed("ra_fits", nocf_aipw(s, learner, ARM_RA_NOCF))
    nc_sec <- now() - t_nc
    u8n <- check_u8(cn$phi, s$w, s$h, s$k)
    wor_n <- se2_two_stage_wor(cn$phi, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
    add_checks(u8_rows(u8n), b, r, ARM_RA_NOCF, TRUE)
    dn <- cn$diagnostics
    nc_extra <- na_extra
    nc_extra$n_fits <- nrow(dn); nc_extra$n_singular <- sum(dn$singular); nc_extra$n_not_clean <- sum(!dn$clean)
    nc_extra$n_warnings <- sum(dn$n_warnings); nc_extra$clip_hits <- as.integer(cn$clip_hits); nc_extra$fit_seconds <- nc_sec
    rows <- c(rows, list(rhat_cols(row(ARM_RA_NOCF, cn, wor_n$se2_wor, u8n, cbind(nc_extra, na_tmle)), cn$fits)))
    add_checks(check_finite("Ra-AIPW-nocf psi, se, Z, phi_hat", list(cn$psi, cn$se, rows[[length(rows)]]$Z, cn$phi)),
               b, r, ARM_RA_NOCF, TRUE)
    more <- c(more, sprintf("Ra-nocf %.5f (se %.5f)", cn$psi, cn$se))
  }, error = function(e) add_error(sprintf("Ra no-CF arm (pop %d, rep %d)", b, r), e, b, r, ARM_RA_NOCF))

  # ---- R-b (forest, DEV): AIPW (+ paired TMLE where the path has it) on the same sample
  # and folds. No forest T1/T2, so no U9 either. A learner error stops only this
  # arm for this replicate and is recorded as a failed required check (C-RB-7).
  if (HAS_RB) tryCatch({
    t_rb <- now()
    learner_b <- function(train, fold) fit_rb(train, rb_seeds(ROOT, CELL, b, r, fold), RA_Q_CLIP, RA_G_TRUNC)
    cfb <- timed("rb_forests", psu_cf_aipw(s, fold, learner_b, arm = ARM_RB_AIPW))   # forests + scoring
    rb_fit_sec <- now() - t_rb
    u3b <- check_u3(s$w, cfb$fold, cfb$omega, cfb$N_hat_v)
    u45b <- check_u4_u5(s, cfb)
    u7b <- check_u7(s, cfb$fold, cfb$fits, cfb$psi, or$psi)
    u8b <- check_u8(cfb$phi, s$w, s$h, s$k)
    wor_b <- se2_two_stage_wor(cfb$phi, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
    add_checks(u8_rows(u8b), b, r, ARM_RB_AIPW, TRUE)
    add_checks(u3b, b, r, ARM_RB_AIPW, TRUE)
    add_checks(u45b, b, r, ARM_RB_AIPW, TRUE)
    add_checks(u7b, b, r, ARM_RB_AIPW, TRUE)
    # The served predictions (explicit training seed, C-RB-8) are identical() to
    # the seed = NULL predictions, on every fold's scored units (required).
    pc <- do.call(rbind, lapply(seq_len(V), function(v)
      cbind(fold = v, rb_prediction_check(cfb$fits[[v]], s[cfb$score_rows[[v]], , drop = FALSE]))))
    add_checks(data.frame(check = sprintf("Rb predictions: explicit seed identical() to seed = NULL (OI-RB-3; fold %d, %s, %d units)",
                                          pc$fold, pc$component, pc$n),
                          value = pc$max_abs_diff, reference = 0, abs_diff = pc$max_abs_diff, tolerance = 0,
                          pass = pc$identical), b, r, ARM_RB_AIPW, TRUE)
    db <- cfb$diagnostics
    rb_common <- data.frame(omega_1 = cfb$omega[1], omega_2 = cfb$omega[2],
                            n_fits = nrow(db), n_not_clean = sum(db$n_warnings > 0), n_warnings = sum(db$n_warnings),
                            clip_hits = as.integer(cfb$clip_hits),
                            rb_oob_mse_Q1 = mean(db$oob_mse[db$component == "Q1"]),
                            rb_oob_mse_Q0 = mean(db$oob_mse[db$component == "Q0"]),
                            rb_oob_mse_g = mean(db$oob_mse[db$component == "g"]))
    rb_aipw_extra <- cbind(rb_common, data.frame(U3_pass = all(u3b$pass), U4_pass = u45b$pass[1],
                           U5_pass = all(u45b$pass[2:3]), U7_lhs = u7b$value, U7_rhs = u7b$reference,
                           U7_abs_diff = u7b$abs_diff, U7_tol = u7b$tolerance, U7_pass = u7b$pass,
                           fit_seconds = rb_fit_sec))
    rb_new <- list(rhat_cols(row(ARM_RB_AIPW, cfb, wor_b$se2_wor, u8b, rb_aipw_extra), cfb$fits, forest = TRUE))
    add_checks(check_finite("Rb-AIPW psi, se, Z, phi_hat", list(cfb$psi, cfb$se, rb_new[[1]]$Z, cfb$phi)),
               b, r, ARM_RB_AIPW, TRUE)
    rb_line <- sprintf("Rb-AIPW %.5f (se %.5f)  forests %.1fs", cfb$psi, cfb$se, rb_fit_sec)
    if (HAS_RB_TMLE) {
      t_tb <- now()
      tmb <- timed("targeting", psu_cf_tmle(s, cfb, RA_Q_CLIP, RA_G_TRUNC, arm = ARM_RB_TMLE))
      rb_tmle_sec <- now() - t_tb
      chtb <- check_tmle(s, cfb, tmb)
      u8tb <- check_u8(tmb$phi, s$w, s$h, s$k)
      wor_tb <- se2_two_stage_wor(tmb$phi, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
      add_checks(u8_rows(u8tb), b, r, ARM_RB_TMLE, TRUE)
      add_checks(chtb[1:5, ], b, r, ARM_RB_TMLE, TRUE)
      add_checks(chtb[6, ], b, r, ARM_RB_TMLE, FALSE)
      add_checks(chtb[7, ], b, r, ARM_RB_TMLE, FALSE)
      eps_b <- sqrt(sum(tmb$eps^2))
      rb_tmle_extra <- cbind(rb_common, data.frame(
        fit_seconds = rb_tmle_sec, eps_1 = tmb$eps[1], eps_0 = tmb$eps[2], eps_norm = eps_b,
        m14_eps_norm = m_s^0.25 * eps_b, lambda_min = tmb$lambda_min, gram_premise_met = chtb$pass[7],
        conformance_max_abs = chtb$value[4], score_max_abs = chtb$value[5], tmle_iter = as.integer(tmb$iter),
        tmle_seconds = rb_tmle_sec, paired_tmle_minus_aipw = tmb$psi - cfb$psi,
        sqrt_m_abs_paired = sqrt(m_s) * abs(tmb$psi - cfb$psi),
        se2_diff_over_sigma2_m = (tmb$se2 - cfb$se2) / cen$sigma2_m,
        m14_abs_se2_diff_over_sigma2_m = m_s^0.25 * abs(tmb$se2 - cfb$se2) / cen$sigma2_m))
      rb_new <- c(rb_new, list(row(ARM_RB_TMLE, tmb, wor_tb$se2_wor, u8tb, rb_tmle_extra)))
      add_checks(check_finite("Rb-TMLE psi, se, Z, eps, phi_star", list(tmb$psi, tmb$se, rb_new[[2]]$Z, tmb$eps, tmb$phi)),
                 b, r, ARM_RB_TMLE, TRUE)
      rb_line <- sprintf("Rb-AIPW %.5f (se %.5f)  Rb-TMLE %.5f  forests %.1fs", cfb$psi, cfb$se, tmb$psi, rb_fit_sec)
    }
    rows <- c(rows, rb_new)
    more <- c(more, rb_line)
  }, error = function(e) add_error(sprintf("Rb arm (pop %d, rep %d)", b, r), e, b, r, ARM_RB_AIPW))

  # ---- R-c (E3 foil, C-RC-1..3): PSU-CF, no-CF (fold 0) and individual-CF AIPW. No
  # T1/T2 or U9. A learner error stops only that arm for the replicate.
  rc_learner <- function(arm_key) function(train, fold)
    fit_rc(train, rb_seeds(ROOT, CELL, b, r, fold, arm = arm_key), RA_Q_CLIP, RA_G_TRUNC)
  rc_row <- function(arm_key, est, fit_sec, u8, u3 = NULL, u7 = NULL, extra_checks = NULL) {
    dc <- est$diagnostics
    ex <- data.frame(n_fits = nrow(dc), n_not_clean = sum(dc$n_warnings > 0), n_warnings = sum(dc$n_warnings),
                     clip_hits = as.integer(est$clip_hits), fit_seconds = fit_sec,
                     omega_1 = if (is.null(est$omega)) NA_real_ else est$omega[1],
                     omega_2 = if (is.null(est$omega)) NA_real_ else est$omega[2],
                     U3_pass = if (is.null(u3)) NA else all(u3$pass), U7_lhs = if (is.null(u7)) NA_real_ else u7$value,
                     U7_rhs = if (is.null(u7)) NA_real_ else u7$reference, U7_abs_diff = if (is.null(u7)) NA_real_ else u7$abs_diff,
                     U7_pass = if (is.null(u7)) NA else u7$pass)
    wor <- se2_two_stage_wor(est$phi, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
    rw <- rhat_cols(row(arm_key, est, wor$se2_wor, u8, ex), est$fits, forest = TRUE)
    add_checks(check_finite(paste(arm_key, "psi, se, Z, phi_hat"), list(est$psi, est$se, rw$Z, est$phi)), b, r, arm_key, TRUE)
    rw
  }
  if (HAS_RC_PSU) tryCatch({
    t_c <- now()
    cfc <- timed("rc_forests", psu_cf_aipw(s, fold, rc_learner(ARM_RC_PSUCF), arm = ARM_RC_PSUCF))
    c_sec <- now() - t_c
    u3c <- check_u3(s$w, cfc$fold, cfc$omega, cfc$N_hat_v); u45c <- check_u4_u5(s, cfc)
    u7c <- check_u7(s, cfc$fold, cfc$fits, cfc$psi, or$psi); u8c <- check_u8(cfc$phi, s$w, s$h, s$k)
    add_checks(u8_rows(u8c), b, r, ARM_RC_PSUCF, TRUE); add_checks(u3c, b, r, ARM_RC_PSUCF, TRUE)
    add_checks(u45c, b, r, ARM_RC_PSUCF, TRUE); add_checks(u7c, b, r, ARM_RC_PSUCF, TRUE)
    pcc <- do.call(rbind, lapply(seq_len(V), function(v)
      cbind(fold = v, rb_prediction_check(cfc$fits[[v]], s[cfc$score_rows[[v]], , drop = FALSE]))))
    add_checks(data.frame(check = sprintf("Rc predictions: explicit seed identical() to seed = NULL (C-RB-8; fold %d, %s, %d units)",
                                          pcc$fold, pcc$component, pcc$n),
                          value = pcc$max_abs_diff, reference = 0, abs_diff = pcc$max_abs_diff, tolerance = 0,
                          pass = pcc$identical), b, r, ARM_RC_PSUCF, TRUE)
    rows <- c(rows, list(rc_row(ARM_RC_PSUCF, cfc, c_sec, u8c, u3c, u7c)))
    more <- c(more, sprintf("Rc-psucf %.5f (se %.5f, %.1fs)", cfc$psi, cfc$se, c_sec))
    rm(cfc)
  }, error = function(e) add_error(sprintf("Rc arm psucf (pop %d, rep %d)", b, r), e, b, r, ARM_RC_PSUCF))
  if (HAS_RC_NOCF) tryCatch({
    t_c <- now()
    cnc <- timed("rc_forests", nocf_aipw(s, rc_learner(ARM_RC_NOCF), ARM_RC_NOCF))
    c_sec <- now() - t_c
    u8c <- check_u8(cnc$phi, s$w, s$h, s$k)
    add_checks(u8_rows(u8c), b, r, ARM_RC_NOCF, TRUE)
    rows <- c(rows, list(rc_row(ARM_RC_NOCF, cnc, c_sec, u8c)))
    more <- c(more, sprintf("Rc-nocf %.5f (se %.5f, %.1fs)", cnc$psi, cnc$se, c_sec))
    rm(cnc)
  }, error = function(e) add_error(sprintf("Rc arm nocf (pop %d, rep %d)", b, r), e, b, r, ARM_RC_NOCF))
  ufold <- if (HAS_RC_IND) timed("sampling", dev_unit_folds(s, key, V, ROOT)) else NULL   # level-3, unit-s1
  ind_ref <- NULL
  if (HAS_RC_IND) tryCatch({
    t_c <- now()
    cic <- timed("rc_forests", indiv_cf_aipw(s, ufold, rc_learner(ARM_RC_INDCF), ARM_RC_INDCF))
    c_sec <- now() - t_c
    u3c <- check_u3(s$w, cic$fold, cic$omega, cic$N_hat_v)
    u7c <- check_u7(s, cic$fold, cic$fits, cic$psi, or$psi); u8c <- check_u8(cic$phi, s$w, s$h, s$k)
    add_checks(u8_rows(u8c), b, r, ARM_RC_INDCF, TRUE); add_checks(u3c, b, r, ARM_RC_INDCF, TRUE)
    add_checks(check_indcf(s, cic), b, r, ARM_RC_INDCF, TRUE); add_checks(u7c, b, r, ARM_RC_INDCF, TRUE)
    if (HAS_RC_COARSE) {                          # its prediction identity, as for PSU-CF (required here)
      pci <- do.call(rbind, lapply(seq_len(V), function(v)
        cbind(fold = v, rb_prediction_check(cic$fits[[v]], s[cic$score_rows[[v]], , drop = FALSE]))))
      add_checks(data.frame(check = sprintf("Rc predictions: explicit seed identical() to seed = NULL (C-RB-8; indcf fold %d, %s, %d units)",
                                            pci$fold, pci$component, pci$n),
                            value = pci$max_abs_diff, reference = 0, abs_diff = pci$max_abs_diff, tolerance = 0,
                            pass = pci$identical), b, r, ARM_RC_INDCF, TRUE)
    }
    rows <- c(rows, list(rc_row(ARM_RC_INDCF, cic, c_sec, u8c, u3c, u7c)))
    more <- c(more, sprintf("Rc-indcf %.5f (se %.5f, %.1fs)", cic$psi, cic$se, c_sec))
    ind_ref <- list(fold = cic$fold, fit_rows = cic$fit_rows, score_rows = cic$score_rows, psi = cic$psi, se2 = cic$se2,
                    seeds = lapply(cic$fits, function(f) f$diagnostics$seed))
    rm(cic)
  }, error = function(e) add_error(sprintf("Rc arm indcf (pop %d, rep %d)", b, r), e, b, r, ARM_RC_INDCF))

  # ---- leakage-neutral sensitivity: the same individual-CF arm, C_psu coarsened to quartile bins
  # for the learner only (C-RC-4), on the same unit folds and level-4 seeds (C-RC-5). A sensitivity
  # analysis that ALSO REMOVES PREDICTIVE INFORMATION: C_psu is a true confounder of Q0 and g0.
  if (HAS_RC_COARSE) tryCatch({
    if (is.null(ind_ref)) stop("the paired Rc-AIPW-indcf arm did not complete")
    coarse_learner <- function(train, fold)
      fit_rc_coarse(train, rb_seeds(ROOT, CELL, b, r, fold, arm = ARM_RC_INDCF), RA_Q_CLIP, RA_G_TRUNC)
    t_c <- now()
    cco <- timed("rc_forests", indiv_cf_aipw(s, ufold, coarse_learner, ARM_RC_INDCF_COARSE))
    c_sec <- now() - t_c
    u3c <- check_u3(s$w, cco$fold, cco$omega, cco$N_hat_v)
    u7c <- check_u7(s, cco$fold, cco$fits, cco$psi, or$psi); u8c <- check_u8(cco$phi, s$w, s$h, s$k)
    add_checks(u8_rows(u8c), b, r, ARM_RC_INDCF_COARSE, TRUE); add_checks(u3c, b, r, ARM_RC_INDCF_COARSE, TRUE)
    add_checks(check_indcf(s, cco), b, r, ARM_RC_INDCF_COARSE, TRUE); add_checks(u7c, b, r, ARM_RC_INDCF_COARSE, TRUE)
    seeds_c <- lapply(cco$fits, function(f) f$diagnostics$seed)
    ncp <- vapply(cco$fits, function(f) f$n_cpsu_values_train, integer(1))
    shared <- c(identical(cco$fold, ind_ref$fold), identical(cco$fit_rows, ind_ref$fit_rows),
                identical(cco$score_rows, ind_ref$score_rows), identical(seeds_c, ind_ref$seeds), all(ncp <= 4L))
    add_checks(data.frame(check = c("R15 pairing: unit fold vector identical() to Rc-AIPW-indcf's",
                                    "R15 pairing: every fold's training rows identical() to Rc-AIPW-indcf's",
                                    "R15 pairing: every fold's scored rows identical() to Rc-AIPW-indcf's",
                                    "R15 pairing: the same three ranger seeds per fold as Rc-AIPW-indcf (C-RC-5)",
                                    "R15 coarsening: every fold's learner saw at most 4 C_psu values (quartile bins, C-RC-4)"),
                          value = c(as.numeric(shared[1:4]), max(ncp)), reference = c(1, 1, 1, 1, 4), abs_diff = NA_real_,
                          tolerance = NA_real_, pass = shared), b, r, ARM_RC_INDCF_COARSE, TRUE)
    pco <- do.call(rbind, lapply(seq_len(V), function(v)
      cbind(fold = v, rb_prediction_check(cco$fits[[v]]$inner, coarsen_cpsu(s[cco$score_rows[[v]], , drop = FALSE])))))
    add_checks(data.frame(check = sprintf("Rc predictions (coarsened C_psu): explicit seed identical() to seed = NULL (C-RB-8; fold %d, %s, %d units)",
                                          pco$fold, pco$component, pco$n),
                          value = pco$max_abs_diff, reference = 0, abs_diff = pco$max_abs_diff, tolerance = 0,
                          pass = pco$identical), b, r, ARM_RC_INDCF_COARSE, TRUE)
    rwc <- rc_row(ARM_RC_INDCF_COARSE, cco, c_sec, u8c, u3c, u7c)
    rwc$paired_minus_indcf_psi <- cco$psi - ind_ref$psi
    rwc$paired_minus_indcf_se2_over_sigma2_m <- (cco$se2 - ind_ref$se2) / cen$sigma2_m
    rows <- c(rows, list(rwc))
    more <- c(more, sprintf("Rc-indcf-coarse %.5f (se %.5f, %.1fs; minus indcf %+.2e)", cco$psi, cco$se, c_sec, cco$psi - ind_ref$psi))
    rm(cco)
  }, error = function(e) add_error(sprintf("Rc arm indcf-coarse (pop %d, rep %d)", b, r), e, b, r, ARM_RC_INDCF_COARSE))

  # ---- R-d (fixed-strata): the four throttled AIPW rungs on the same sample and folds.
  # alpha enters at level 4 only (C-RD-1), so every rung shares Ra-AIPW's sample and fold
  # structure; each draws its own q_m(alpha) whole-PSU subsample of each training fold. An
  # error (including q_m > pool) stops only that rung for this replicate and is recorded as a
  # failed required check; there is never a fallback to the full-training R-a fit (C-RD-4).
  rd_rows <- list(); rd_subs <- list()
  if (HAS_RD) {
    for (alpha in RD_RUN) {
      arm_d <- rd_arm(alpha)
      tryCatch(ra_fit_failed = function(e) fail_row(arm_d, conditionMessage(e)), expr = {
        t_rd <- now()
        learner_d <- function(train, fold) fit_rd(train, alpha, M, ROOT, CELL, b, r, fold, a_cpsu_term = a_term)
        cfd <- timed("rd_fits", psu_cf_aipw(s, fold, learner_d, arm = arm_d))   # R-a on the subsample + scoring
        rd_fit_sec <- now() - t_rd
        u3d <- check_u3(s$w, cfd$fold, cfd$omega, cfd$N_hat_v)
        u45d <- check_u4_u5(s, cfd)
        u7d <- check_u7(s, cfd$fold, cfd$fits, cfd$psi, or$psi)
        u8d <- check_u8(cfd$phi, s$w, s$h, s$k)
        wor_d <- se2_two_stage_wor(cfd$phi, s$w, s$h, s$k, s$i, s$M_pop_h, s$M_k)
        t2d <- timed("remainders", t2_refined(ctx, cfd$fits, cfd$omega, T2P, T2_CACHE))   # refined T2
        t2dc <- t2_refined_cols(t2d, cfd$psi, or$psi, sig_wor); T1d <- t2dc$T1
        t2_diag[[length(t2_diag) + 1L]] <<- t2_refined_diag(t2d, b, r, arm_d)
        add_checks(u8_rows(u8d), b, r, arm_d, TRUE)
        add_checks(u3d, b, r, arm_d, TRUE)
        add_checks(u45d, b, r, arm_d, TRUE)
        add_checks(u7d, b, r, arm_d, TRUE)
        add_checks(check_rd_throttle(s, cfd, alpha, M), b, r, arm_d, TRUE)
        shared <- c(identical(cfd$fold, fref$fold), identical(cfd$fit_rows, fref$fit_rows),
                    identical(cfd$score_rows, fref$score_rows))
        add_checks(data.frame(check = paste("Rd shared", arm_d, c(
                                ": fold vector identical() to the replicate's fold structure (shared with Ra-AIPW)",
                                ": every fold's training pool identical() to the replicate's (Ra-AIPW's)",
                                ": every fold's scored units identical() to the replicate's (Ra-AIPW's)")),
                              value = as.numeric(shared), reference = 1, abs_diff = NA_real_, tolerance = NA_real_,
                              pass = shared), b, r, arm_d, TRUE)
        if (length(sel)) add_checks(timed("u9_bruteforce", u9_cross_check_refined(ctx, cfd$fits, t2d, sel, T2P, T2_CACHE)), b, r, arm_d, TRUE)
        dd <- cfd$diagnostics
        th <- lapply(cfd$fits, `[[`, "throttle")
        rd_extra <- data.frame(t2dc, omega_1 = cfd$omega[1], omega_2 = cfd$omega[2],
                               U3_pass = all(u3d$pass), U4_pass = u45d$pass[1], U5_pass = all(u45d$pass[2:3]),
                               U7_lhs = u7d$value, U7_rhs = u7d$reference, U7_abs_diff = u7d$abs_diff,
                               U7_tol = u7d$tolerance, U7_pass = u7d$pass,
                               n_fits = nrow(dd), n_singular = sum(dd$singular), n_not_clean = sum(!dd$clean),
                               n_warnings = sum(dd$n_warnings), clip_hits = as.integer(cfd$clip_hits),
                               fit_seconds = rd_fit_sec, rd_alpha = alpha, rd_q = th[[1]]$q, rd_pool = th[[1]]$pool,
                               rd_share = th[[1]]$q / th[[1]]$pool,
                               rd_units_fit_mean = mean(vapply(th, `[[`, numeric(1), "n_units_fit")),
                               rd_units_pool_mean = mean(vapply(th, `[[`, numeric(1), "n_units_pool")))
        rd_rows[[arm_d]] <- rhat_cols(row(arm_d, cfd, wor_d$se2_wor, u8d, cbind(rd_extra, na_tmle)), cfd$fits)
        add_checks(t2_finite_check(paste(arm_d, "psi, se, Z, T1, T2, phi_hat"), t2d, list(cfd$psi, cfd$se, rd_rows[[arm_d]]$Z, cfd$phi), T1d),
                   b, r, arm_d, TRUE)
        rd_subs[[arm_d]] <- lapply(th, function(x) paste(x$psus$h, x$psus$k))
        more <- c(more, sprintf("%s %.5f (T2 %+.1e, q %d)", arm_d, cfd$psi, t2d$T2, th[[1]]$q))
      }, error = function(e) add_error(sprintf("Rd arm %s (pop %d, rep %d)", arm_d, b, r), e, b, r, arm_d))
    }
    # The four rungs draw distinct subsamples (distinct level-4 keys), fold by fold (info).
    if (length(rd_subs) == length(RD_RUN) && length(RD_RUN) > 1L) {
      distinct <- vapply(seq_len(V), function(v) length(unique(lapply(rd_subs, `[[`, v))) == length(RD_RUN), logical(1))
      add_checks(data.frame(check = sprintf("Rd subsamples: the %d alpha rungs run drew distinct PSU sets (fold %d; info)", length(RD_RUN), seq_len(V)),
                            value = as.numeric(distinct), reference = 1, abs_diff = NA_real_, tolerance = NA_real_,
                            pass = distinct), b, r, "all", FALSE)
    }
  }
  rows <- rbind_fill(c(rows, unname(rd_rows)))
  results[[length(results) + 1L]] <<- rows
  cat(sprintf("  pop %d rep %2d: %s  rep %.1fs\n", b, r, line, now() - t_rep))
  if (length(more)) cat(paste0("                  ", more, "\n"), sep = "")
  invisible()
}

# ---- E3 ctrl: the FROZEN kappa_AE record (C-AE-3). kappa_AE was derived once, by census algebra,
# from the complete production panel population set (FS populations 1..160 at the panel cell) by
# runs/run_prep.R --job panel; this run reads it, verifies its provenance (rule, root, master, population set,
# cell, specification identifier, P0 grid hash, content hash) and never recalibrates. Not buildable -> the AE arm is skipped and
# that is reported, never resized.
AE <- NULL; panel_pred <- list()
if (HAS_AE) {
  AE <- ae_frozen_load(file.path(engine_dir, AE_FROZEN_FILE), spec_sha = if (MODE == "dev") NULL else RM$spec$actual,
                       grid_hash = ctx$primary_hash)
  add_checks(data.frame(check = c(sprintf("AE frozen record verified (C-AE-3, R23): kappa_AE = %.10f from FS populations 1..160, content sha256 %s",
                                          AE$kappa, substr(AE$content_sha256, 1, 12)),
                                  sprintf("AE clip bound: eta0 inside the fitted clip on the covariate box (recorded: g in [%.4f, %.4f], Q in [%.4f, %.4f]; C-AE-2)",
                                          AE$clip$min_g, AE$clip$max_g, AE$clip$min_Q, AE$clip$max_Q),
                                  "AE map: T2(kappa) strictly increasing on [0, kappa_max] (recorded at calibration; C-AE-3)"),
                        value = c(AE$kappa, as.numeric(isTRUE(AE$clip$inside)), as.numeric(isTRUE(AE$monotone))),
                        reference = c(NA, 1, 1), abs_diff = NA_real_, tolerance = NA_real_,
                        pass = c(is.finite(AE$kappa) || !isTRUE(AE$buildable), isTRUE(AE$clip$inside), isTRUE(AE$monotone))),
             NA, NA, ARM_AE_AIPW, TRUE)
  if (!isTRUE(AE$buildable)) { ARMS <- setdiff(ARMS, ARM_AE_AIPW); HAS_AE <- FALSE }
  cat(sprintf("[dev-e2e] AE frozen record: %s\n", if (isTRUE(AE$buildable)) sprintf("kappa_AE = %.6f, T2 = %.4e, Q_0.9 = %.4f (populations 1..160)", AE$kappa, AE$T2, AE$stat) else
    paste("NOT BUILDABLE:", AE$reason)))
}

# ---- populations and replicates ---------------------------------------------------
for (bi in seq_along(POP_IDS)) {
  b <- POP_IDS[bi]
  t_pop <- now()
  pop <- NULL
  tryCatch({
    # The cell's (H, J) block of master MASTER_ID, built directly (level-1 stream). By entity addressing
    # it equals the master's prefix view; REQUIRED nesting check on a sub-block (A2's method).
    pop <- build_population(H, J, co, MASTER_ID, b, ROOT)
    Hs <- min(H, 10L); Js <- min(J, 40L)
    prefix_ok <- identical(population_prefix(pop, Hs, Js), build_population(Hs, Js, co, MASTER_ID, b, ROOT))
    add_checks(data.frame(check = sprintf("nesting: prefix (%d,%d) of the cell's (%d,%d) block identical() to a direct build (master %s; entity addressing)", Hs, Js, H, J, MASTER_ID),
                          value = as.numeric(prefix_ok), reference = 1, abs_diff = NA_real_,
                          tolerance = NA_real_, pass = prefix_ok), b, NA, "all", TRUE)
    cen <- census_targets(pop, m_h = M_H, n0 = N0)                             # census algebra,
    # U2's analytic SD of N_hat/N, for the record: Var(N_hat) = sum_h M_pop,h^2 (1 - f_h) S2_M,h / m_h.
    M_pop_h <- as.vector(table(pop$psu$h)); S2_M <- as.vector(tapply(pop$psu$M, pop$psu$h, stats::var))
    sd_Nhat_rel <- if (WR) sqrt(sum(M_pop_h^2 / M_H * (1 - 1 / M_pop_h) * S2_M)) / cen$N else   # U2, WR form
      sqrt(sum(M_pop_h^2 * (1 - M_H / M_pop_h) * S2_M / M_H)) / cen$N
    if (WR) { wrt <- census_wr_target(cen); cen$sigma2_m_wr <- wrt$sigma2_m_wr; cen$wr_gap <- wrt$wr_gap }
    ps <- cen$per_stratum
    fs <- census_fpc_summary(cen)                                              # E2 census algebra
    census_strata[[bi]] <- cbind(data.frame(pop_id = b, KB = kb_of_stratum(fs$per_stratum$h)), fs$per_stratum)
    census[[bi]] <- data.frame(
      pop_id = b, master_id = MASTER_ID, root = ROOT, H = H, J = J, N = cen$N, M_pop = cen$M_pop,
      m = cen$m, n0 = cen$n0, Psi_N = cen$Psi_N, Psi_N_gf = cen$Psi_N_gf, Psi = Psi_grid,
      Psi_N_minus_Psi = cen$Psi_N - Psi_grid, sigma2_m = cen$sigma2_m, sigma2_m_wor = cen$sigma2_m_wor,
      Delta_m = cen$Delta_m, sigma_m = sqrt(cen$sigma2_m), sigma_m_wor = sqrt(cen$sigma2_m_wor),
      Delta_over_sigma2_m = cen$Delta_m / cen$sigma2_m, mean_f_h = mean(ps$f_h), mean_beta_h = mean(ps$beta_h),
      sd_Nhat_over_N_analytic = sd_Nhat_rel, min_g0 = min(pop$units$g), max_g0 = max(pop$units$g),
      min_Q = min(pop$units$Q1, pop$units$Q0), max_Q = max(pop$units$Q1, pop$units$Q0),
      n_size_ties = pop$meta$n_size_ties,
      f = fs$f_bar, fpc = fs$fpc_bar, sigma_ratio = fs$sigma_ratio, beta_pooled = fs$beta_pooled,
      single_stage_fpc_deficit = fs$single_stage_fpc_deficit,
      deficit_over_sigma2_m_wor = fs$single_stage_fpc_deficit / cen$sigma2_m_wor,
      c3b_overlay_080 = fs$overlay[["0.80"]], c3b_overlay_090 = fs$overlay[["0.90"]], c3b_overlay_095 = fs$overlay[["0.95"]],
      min_M = min(pop$psu$M), max_M = max(pop$psu$M), min_pi_unit = N0 / max(pop$psu$M),
      label = RUN_LABEL, stringsAsFactors = FALSE)
    if (WR) { census[[bi]]$sigma2_m_wr <- cen$sigma2_m_wr; census[[bi]]$wr_gap <- cen$wr_gap }
    if (U10) {                                    # U10 census algebra (C-U10-1) and G6 limb (a), boundary reading
      u10c <- u10_census(pop, M_H, N0)
      pred_z <- u10c$bias / sqrt(u10c$sigma2_m_u10)
      g6u <- g6_boundary(pred_z, -0.10, 0.10)
      census[[bi]]$Psi_U10 <- u10c$Psi_U10; census[[bi]]$sigma2_m_u10 <- u10c$sigma2_m_u10
      census[[bi]]$u10_bias <- u10c$bias; census[[bi]]$u10_pred_mean_Z <- pred_z; census[[bi]]$u10_g6_distance <- g6u$distance
      add_checks(data.frame(check = sprintf("U10 G6 limb (a) (report, R5): predicted mean Z_U10 = (Psi_U10 - Psi_N)/sigma_m,U10 vs [-0.10, 0.10] -- %s", g6u$status),
                            value = pred_z, reference = 0.10, abs_diff = g6u$distance, tolerance = NA_real_, pass = TRUE),
                 b, NA, ARM_ORACLE_U10, FALSE)
    }
    if (E3 == "ctrl") {                           # the panel's census predictions and G6 limb (a) distances
      wcc <- wrong_clustering_census(pop, M_H, N0)
      pp <- panel_predictions(cen, wcc, ae_stat = if (!is.null(AE) && AE$buildable) AE$T2 / sqrt(cen$sigma2_m_wor) else NA_real_)
      g6p <- g6_boundary(pp$predicted, pp$lo, pp$hi)
      panel_pred[[bi]] <- cbind(data.frame(pop_id = b), pp, g6p[, c("distance", "inside", "status")])
      gf_off <- (cen$Psi_N - cen$Psi_N_gf) / sqrt(cen$sigma2_m)
      census[[bi]]$Z_gf_pred_offset <- gf_off
      for (v in c("nomh", "unit", "nostrat")) census[[bi]][[paste0("R_var_", v, "_pred")]] <- wcc$R_var[[v]]
      add_checks(data.frame(check = sprintf("panel G6 limb (a) (report, R5): %s, predicted %s vs [%s, %s] -- %s", pp$control,
                                            formatC(pp$predicted, format = "f", digits = 4), pp$lo, pp$hi, g6p$status),
                            value = pp$predicted, reference = NA_real_, abs_diff = g6p$distance, tolerance = NA_real_, pass = TRUE),
                 b, NA, ARM_RA_AIPW, FALSE)
    }
    add_checks(check_finite("census targets (Psi_N, sigma2_m, sigma2_m_wor, Delta_m)",
                            list(cen$Psi_N, cen$sigma2_m, cen$sigma2_m_wor, cen$Delta_m)), b, NA, "all", TRUE)
    # Census FPC algebra (required): Delta_m = sigma2_m - sigma2_m_wor = sum_h f_h B_h.
    add_checks(data.frame(check = "census FPC algebra: sigma2_m - sigma2_m_wor = Delta_m = sum_h f_h B_h (64 eps sigma2_m)",
                          value = cen$sigma2_m - cen$sigma2_m_wor, reference = cen$Delta_m, abs_diff = fs$identity_gap,
                          tolerance = 64 * .Machine$double.eps * cen$sigma2_m,
                          pass = fs$identity_gap <= 64 * .Machine$double.eps * cen$sigma2_m), b, NA, "all", TRUE)
    if (E2C %in% c("fpc", "wr"))                                                # G2, E2 only
      add_checks(data.frame(check = c("G2 feasibility: n0/M_k >= 0.40 for every PSU of the population (min n0/M_k in value)",
                                      "G2 feasibility: n0 <= min M_k (pi_{i|k} <= 1)"),
                            value = c(N0 / max(pop$psu$M), min(pop$psu$M)), reference = c(0.40, N0),
                            abs_diff = NA_real_, tolerance = NA_real_,
                            pass = c(N0 / max(pop$psu$M) >= 0.40, N0 <= min(pop$psu$M))), b, NA, "all", TRUE)
    frame <- sampling_frame(pop, H, J)
    TIME$populations <- TIME$populations + now() - t_pop
    cat(sprintf("[dev-e2e] population %d: N = %d, Psi_N = %.6f, sigma_m = %.5f, sigma_m_wor = %.5f (%.1f s)\n",
                b, cen$N, cen$Psi_N, sqrt(cen$sigma2_m), sqrt(cen$sigma2_m_wor), now() - t_pop))
    for (ri in seq_along(REP_IDS)) {
      r <- REP_IDS[ri]
      t_r <- now(); before <- unlist(mget(COMPONENTS, envir = TIME)); rng_r <- ambient_rng()
      tryCatch(run_rep(b, r, bi, ri, pop, frame, cen),
               error = function(e) add_error(sprintf("replicate (pop %d, rep %d)", b, r), e, b, r, "all"))
      rep_total <- now() - t_r
      d <- unlist(mget(COMPONENTS, envir = TIME)) - before
      timing[[length(timing) + 1L]] <- data.frame(
        pop_id = b, rep_id = r, sampling = d[["sampling"]], ra_fits = d[["ra_fits"]], rb_forests = d[["rb_forests"]],
        rc_forests = d[["rc_forests"]], rd_fits = d[["rd_fits"]], remainders = d[["remainders"]], u9_bruteforce = d[["u9_bruteforce"]], targeting = d[["targeting"]],
        rhat = d[["rhat"]],
        checks = rep_total - sum(d), replicate_total = rep_total)
      # One progress line per replicate, so a run stopped at its time limit still shows its cost.
      cat(sprintf("[rep] %s m=%d J=%d pop %d rep %d: %.1f s (Ra fits %.1f, %s %.1f, remainders %.1f, U9 bf %.1f, targeting %.2f, sampling %.1f, checks %.1f); peak RSS %s MB\n",
                  REGIME, M, J, b, r, rep_total, d[["ra_fits"]], if (HAS_RD) "Rd fits" else if (HAS_RC) "Rb/Rc forests" else "Rb forests",
                  d[["rb_forests"]] + d[["rc_forests"]] + d[["rd_fits"]], d[["remainders"]], d[["u9_bruteforce"]],
                  d[["targeting"]], d[["sampling"]], rep_total - sum(d),
                  format(round(peak_rss_kb_now() / 1024))))
      flush(stdout())
      rng_rep_ok <- identical(rng_r, ambient_rng())
      add_checks(data.frame(check = "ambient RNG state and kind unchanged by the replicate (OI-RB-3)",
                            value = as.numeric(rng_rep_ok), reference = 1, abs_diff = NA_real_, tolerance = NA_real_,
                            pass = rng_rep_ok), b, r, "all", TRUE)
    }
  }, error = function(e) add_error(sprintf("population %d", b), e, b, NA, "all"))
}

# ---- calibration cell: U2 and U8's three limbs
# over the redraws of each population, |z| <= 3.
# U8 is GATED on the census-centred linearised score with known N, whose design expectations are
# exactly sigma2_m, sigma2_m_wor and sigma2_m_wor - sum_h f_h V_h (exact-enumeration test). The
# actual Hajek estimators target the same quantities only up to O(1/m), so more redraws could
# fail correct code: they are REPORTED beside the gate (required = FALSE). C3a's both-levels
# quantities are reported without a grade. z = (MC estimate - target) / MCSE; the MCSE of an SD
# is SD / sqrt(2 (R - 1)).
zrow_add <- function(bb, arm, lab, est, target, mcse, n, required = TRUE) {
  z <- (est - target) / mcse
  add_checks(data.frame(check = sprintf("%s (|z| <= 3, R = %d; z = %+.2f)", lab, n, z), value = est, reference = target,
                        abs_diff = abs(est - target), tolerance = 3 * mcse, pass = is.finite(z) && abs(z) <= 3),
             bb, NA, arm, required)
}
# ---- the E2 WR code-check over redraws. GATED: the census-centred linearised score
# with known N, E[total] = 0 and E[UC/WR se2] = sigma2_m,WR (exact enumeration), and U2's WR form
# (E[N_hat/N - 1] = 0, SD = sqrt(sum_h M^2/m_h (1 - 1/M) S2_M,h)/N). REPORTED: the Hajek WR se2 against
# sigma2_m,WR (O(1/m) apart; ungated), and the WOR-defined sigma2_m, which correct WR code misses by
# sum_h B_h/M_pop,h.
if (CALIB && WR && length(results)) {
  rdf <- rbind_fill(results); cdf <- if (length(census)) do.call(rbind, census) else data.frame()
  for (bb in unique(rdf$pop_id)) {
    o <- rdf[rdf$arm == ARM_ORACLE & rdf$pop_id == bb, , drop = FALSE]; cr <- cdf[cdf$pop_id == bb, ]
    if (nrow(o) < 2L || !nrow(cr)) next
    n <- nrow(o); x <- o$N_hat_over_N_minus_1; se2 <- o$se_hat^2
    zrow_add(bb, ARM_ORACLE, "U2 (WR form) MC mean of N_hat/N - 1 = 0", mean(x), 0, stats::sd(x) / sqrt(n), n)
    zrow_add(bb, ARM_ORACLE, "U2 (WR form) MC SD of N_hat/N - 1 = analytic sqrt(Var_WR(N_hat))/N", stats::sd(x),
             cr$sd_Nhat_over_N_analytic, stats::sd(x) / sqrt(2 * (n - 1)), n)
    zrow_add(bb, ARM_ORACLE, "R14 WR code-check (linearised, known N): MC mean of sum_i w_i (phi_i - Psi_N)/N = 0",
             mean(o$lin_total), 0, stats::sd(o$lin_total) / sqrt(n), n)
    zrow_add(bb, ARM_ORACLE, "R14 WR code-check (linearised, known N): MC mean of UC/WR se2 = sigma2_m,WR = sum_h [(1 - 1/M_pop,h) B_h + V_h]",
             mean(o$lin_se2_wr), cr$sigma2_m_wr, stats::sd(o$lin_se2_wr) / sqrt(n), n)
    zrow_add(bb, ARM_ORACLE, "R14 WR code-check (Hajek, report): MC mean of UC/WR se2 vs sigma2_m,WR", mean(se2), cr$sigma2_m_wr,
             stats::sd(se2) / sqrt(n), n, required = FALSE)
    zl <- (mean(o$lin_se2_wr) - cr$sigma2_m) / (stats::sd(o$lin_se2_wr) / sqrt(n))
    add_checks(data.frame(check = sprintf("R14 mismatch (report): the linearised WR se2's mean vs the WOR-defined sigma2_m, which it misses by sum_h B_h/M_pop,h = %.4e (z = %+.1f)",
                                          cr$wr_gap, zl),
                          value = mean(o$lin_se2_wr), reference = cr$sigma2_m, abs_diff = abs(mean(o$lin_se2_wr) - cr$sigma2_m),
                          tolerance = NA_real_, pass = TRUE), bb, NA, ARM_ORACLE, FALSE)
  }
}
# ---- U10: its exact linear identities over redraws, GATED at |z| <= 3;
# the power reading (mean Z_U10 keyed to Psi_N against its census prediction) is REPORTED.
if (U10 && length(results)) {
  rdf <- rbind_fill(results); cdf <- if (length(census)) do.call(rbind, census) else data.frame()
  for (bb in unique(rdf$pop_id)) {
    o <- rdf[rdf$arm == ARM_ORACLE_U10 & rdf$pop_id == bb, , drop = FALSE]; cr <- cdf[cdf$pop_id == bb, ]
    if (nrow(o) < 2L || !nrow(cr)) next
    n <- nrow(o); se2 <- o$se_hat^2
    zrow_add(bb, ARM_ORACLE_U10, "U10 (exact, linear): MC mean of psi_U10 = Psi_U10 (the unweighted mean of PSU means)",
             mean(o$psi_hat), cr$Psi_U10, stats::sd(o$psi_hat) / sqrt(n), n)
    zrow_add(bb, ARM_ORACLE_U10, "U10 (exact): MC mean of se2_U10 = sigma2_m,U10", mean(se2), cr$sigma2_m_u10, stats::sd(se2) / sqrt(n), n)
    add_checks(data.frame(check = sprintf("U10 power (report): mean Z_U10 keyed to Psi_N over %d draws vs predicted %+.3f (first order)", n, cr$u10_pred_mean_Z),
                          value = mean(o$Z), reference = cr$u10_pred_mean_Z, abs_diff = abs(mean(o$Z) - cr$u10_pred_mean_Z),
                          tolerance = NA_real_, pass = TRUE), bb, NA, ARM_ORACLE_U10, FALSE)
  }
}
if (CALIB && !WR && length(results)) {
  rdf <- rbind_fill(results)
  cdf <- if (length(census)) do.call(rbind, census) else data.frame()
  for (bb in unique(rdf$pop_id)) {
    o <- rdf[rdf$arm == ARM_ORACLE & rdf$pop_id == bb, , drop = FALSE]; cr <- cdf[cdf$pop_id == bb, ]
    if (nrow(o) < 2L || !nrow(cr)) next
    n <- nrow(o); x <- o$N_hat_over_N_minus_1
    se2 <- o$se_hat^2; wor <- o$se2_two_stage_wor; fp1 <- o$se2_single_stage_fpc
    tgt_neg <- cr$sigma2_m_wor - cr$single_stage_fpc_deficit
    zrow <- function(lab, est, target, mcse, required = TRUE) {
      z <- (est - target) / mcse
      add_checks(data.frame(check = sprintf("%s (|z| <= 3, R = %d; z = %+.2f)", lab, n, z), value = est, reference = target,
                            abs_diff = abs(est - target), tolerance = 3 * mcse, pass = is.finite(z) && abs(z) <= 3),
                 bb, NA, ARM_ORACLE, required)
    }
    zrow("U2 MC mean of N_hat/N - 1 = 0", mean(x), 0, stats::sd(x) / sqrt(n))
    zrow("U2 MC SD of N_hat/N - 1 = analytic sqrt(Var(N_hat))/N", stats::sd(x), cr$sd_Nhat_over_N_analytic,
         stats::sd(x) / sqrt(2 * (n - 1)))
    # gated: the linearised score with known N (exact targets)
    zrow("U8 linearised (known N): MC mean of sum_i w_i (phi_i - Psi_N) / N = 0", mean(o$lin_total), 0,
         stats::sd(o$lin_total) / sqrt(n))
    zrow("U8 single-stage limb (linearised, known N): MC mean of UC se2 = sigma2_m", mean(o$lin_se2_uc), cr$sigma2_m,
         stats::sd(o$lin_se2_uc) / sqrt(n))
    zrow("U8 two-stage limb (linearised, known N): MC mean of two-stage WOR se2 = sigma2_m_wor", mean(o$lin_se2_wor),
         cr$sigma2_m_wor, stats::sd(o$lin_se2_wor) / sqrt(n))
    zrow("U8 negative limb (linearised, known N): MC mean of single-stage-fpc se2 = sigma2_m_wor - sum_h f_h V_h",
         mean(o$lin_se2_fpc1), tgt_neg, stats::sd(o$lin_se2_fpc1) / sqrt(n))
    # reported, not gated: the actual Hajek estimators (estimated psi and N_hat; O(1/m) apart)
    zrow("U8 single-stage limb (Hajek, report): MC mean of UC se2 vs sigma2_m", mean(se2), cr$sigma2_m,
         stats::sd(se2) / sqrt(n), required = FALSE)
    zrow("U8 two-stage limb (Hajek, report): MC mean of two-stage WOR se2 vs sigma2_m_wor", mean(wor), cr$sigma2_m_wor,
         stats::sd(wor) / sqrt(n), required = FALSE)
    zrow("U8 negative limb (Hajek, report): MC mean of single-stage-fpc se2 vs sigma2_m_wor - sum_h f_h V_h", mean(fp1),
         tgt_neg, stats::sd(fp1) / sqrt(n), required = FALSE)
    add_checks(data.frame(check = "U8 negative limb: census deficit sum_h f_h V_h > 0 (the blanket single-stage fpc is strictly below sigma2_m_wor)",
                          value = cr$single_stage_fpc_deficit, reference = 0, abs_diff = NA_real_, tolerance = NA_real_,
                          pass = cr$single_stage_fpc_deficit > 0), bb, NA, ARM_ORACLE, TRUE)
    # C3a's two levels, reported (no grade; C3a is graded on its frozen host after).
    v_psi <- stats::var(o$psi_hat)
    add_checks(data.frame(check = c(sprintf("C3a level 1 (report): Var(psi_oracle) over %d draws vs sigma2_m_wor (ratio in value)", n),
                                    sprintf("C3a level 2 (report): mean UC se2 over %d draws vs sigma2_m (ratio in value)", n),
                                    "C3a target (report): 1 - Delta_m/sigma2_m (census)"),
                          value = c(v_psi / cr$sigma2_m_wor, mean(se2) / cr$sigma2_m, 1 - cr$Delta_m / cr$sigma2_m),
                          reference = c(1, 1, NA_real_), abs_diff = NA_real_, tolerance = NA_real_, pass = TRUE),
               bb, NA, ARM_ORACLE, FALSE)
  }
}

# ---- ambient RNG untouched (REQUIRED): every stream above is scoped ----
rng_ok <- identical(rng_start, ambient_rng())
add_checks(data.frame(check = "ambient RNG state and kind unchanged by the whole run (OI-RB-3)",
                      value = as.numeric(rng_ok), reference = 1, abs_diff = NA_real_, tolerance = NA_real_,
                      pass = rng_ok), NA, NA, "all", TRUE)

# ---- outputs --------------------------------------------------------------------
checks_df <- do.call(rbind, checks)
results_df <- rbind_fill(results)
census_df <- if (length(census)) do.call(rbind, census) else data.frame()
census_strata_df <- if (length(census_strata)) do.call(rbind, census_strata) else data.frame()
status <- required_status(checks_df)

# Full binary64 precision in the CSVs (write.csv would round to 15 digits).
write_csv17 <- function(df, path) {
  out <- df
  for (cl in names(out)) if (is.double(out[[cl]]))
    out[[cl]] <- ifelse(is.na(out[[cl]]), NA_character_, sprintf("%.17g", out[[cl]]))
  utils::write.csv(out, path, row.names = FALSE, na = "NA")
}
t_out <- now()
write_csv17(results_df, file.path(out_dir, "results.csv"))
write_csv17(census_df, file.path(out_dir, "census.csv"))
write_csv17(census_strata_df, file.path(out_dir, "census_strata.csv"))
write_csv17(checks_df, file.path(out_dir, "checks.csv"))
panel_df <- if (length(panel_pred)) do.call(rbind, panel_pred) else data.frame()
if (nrow(panel_df)) write_csv17(panel_df, file.path(out_dir, "panel_predictions.csv"))
timing_df <- if (length(timing)) do.call(rbind, timing) else data.frame()
write_csv17(timing_df, file.path(out_dir, "timing.csv"))
# Every fold fit with a T2 -- its trigger, every rung's value and step, the status (written even when empty).
t2_diag_df <- if (length(t2_diag)) do.call(rbind, t2_diag) else
  data.frame(pop_id = integer(0), rep_id = integer(0), arm = character(0), fold = integer(0), triggered = logical(0), status = character(0))
write_csv17(t2_diag_df, file.path(out_dir, "t2_refinement.csv"))
t_outputs <- now() - t_out

wall <- now() - t_start
# Run time by component, no overlap: the named components, "checks" (everything else
# inside the replicates: U1/U3-U8, the WOR SE, conformance refits, prediction and
# RNG checks, row assembly), "outputs" (CSV writes) and "other" (start-up, rule loading,
# configuration, U9 selection) partition the wall time measured here; manifest and
# SUMMARY writing come after it.
comp_sec <- unlist(mget(COMPONENTS, envir = TIME))
checks_sec <- if (nrow(timing_df)) sum(timing_df$checks) else 0
run_timing <- c(comp_sec, checks = checks_sec, outputs = t_outputs)
run_timing <- c(run_timing, other = wall - sum(run_timing), wall = wall)
# Peak resident memory of this R process (Linux /proc; NA elsewhere) and R's own peak heap.
peak_rss_kb <- peak_rss_kb_now()
gc_max_mb <- sum(gc()[, 6])
utc_end <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
sha <- function(p) if (file.exists(p)) digest::digest(file = p, algo = "sha256") else NA_character_
git <- unname(Sys.which("git"))
git_out <- function(a) if (nzchar(git))
  tryCatch(suppressWarnings(system2(git, c("-C", shQuote(repo_root), a), stdout = TRUE, stderr = TRUE)),
           error = function(e) paste("ERROR:", conditionMessage(e))) else NA_character_
lock <- jsonlite::fromJSON(file.path(engine_dir, "renv.lock"), simplifyVector = FALSE)
lock_ver <- vapply(lock$Packages, `[[`, "", "Version")
loaded <- sort(loadedNamespaces())
pkg_ver <- lapply(stats::setNames(nm = loaded), function(p) {
  # The DESCRIPTION string, as renv.lock records it (packageVersion() would
  # normalise "1.1-37" to "1.1.37" and never match the lock).
  v <- utils::packageDescription(p, fields = "Version")
  list(version = v, lock = if (p %in% names(lock_ver)) unname(lock_ver[p]) else NA,
       match_lock = if (p %in% names(lock_ver)) identical(v, unname(lock_ver[p])) else NA)
})
engine_hashes <- c(
  stats::setNames(lapply(file.path(engine_dir, "R", ENGINE_FILES), sha), paste0("R/", ENGINE_FILES)),
  list("runs/run_e2e_dev.R" = sha(getOption("sim.engine.script")),
       "env/run.R" = sha(file.path(engine_dir, "env", "run.R")),
       "DESCRIPTION" = sha(file.path(engine_dir, "DESCRIPTION"))))
KP <- head(POP_IDS, 5L); KR <- head(REP_IDS, 5L)          # key examples (first 5 IDs on each axis); the rest follow the encoding
keys <- list(
  note = "examples: the first populations and replicates of this run; every other key follows the same encoding",
  population = vapply(KP, function(b)
    rng_key("population", list(root = ROOT, master_id = MASTER_ID, population_id = b)), ""),
  sample = as.vector(outer(KR, KP, Vectorize(function(r, b)
    rng_key("sample", list(root = ROOT, cell = CELL, population_id = b, rep_id = r))))),
  fold = as.vector(outer(KR, KP, Vectorize(function(r, b)
    rng_key("fold", list(root = ROOT, cell = CELL, population_id = b, rep_id = r,
                         fold_structure_id = FOLD_STRUCTURE_PRIMARY))))),
  fold_draw_convention = attr(dev_fold_draw(ROOT), "convention"),
  learner = c(list(Ra = "none: R-a is a deterministic function of its training set (TC-E4); no level-4 stream drawn"),
              stats::setNames(lapply(c(if (HAS_RB) RB_ARM_KEY, if (HAS_RC_PSU) ARM_RC_PSUCF, if (HAS_RC_IND) ARM_RC_INDCF,
                                       if (HAS_RD) vapply(RD_RUN, rd_arm, "")), function(ak)
                as.vector(outer(seq_len(V), seq_len(length(KR) * length(KP)), Vectorize(function(v, g) {
                  b <- KP[(g - 1L) %/% length(KR) + 1L]; r <- KR[(g - 1L) %% length(KR) + 1L]
                  rng_key("learner", list(root = ROOT, cell = CELL, population_id = b, rep_id = r,
                                          arm = ak, fold = v)) })))),
                c(if (HAS_RB) "Rb", if (HAS_RC_PSU) ARM_RC_PSUCF, if (HAS_RC_IND) ARM_RC_INDCF,
                  if (HAS_RD) vapply(RD_RUN, rd_arm, ""))),
              if (HAS_RC_NOCF) list(Rc_nocf_fold0 = as.vector(outer(KR, KP, Vectorize(function(r, b)
                rng_key("learner", list(root = ROOT, cell = CELL, population_id = b, rep_id = r, arm = ARM_RC_NOCF, fold = 0L))))))))
n_req <- sum(checks_df$required); n_req_fail <- sum(checks_df$required & !checks_df$pass)
manifest <- list(
  label = RUN_LABEL, status = RUN_LABEL, mode = MODE,
  spec_freeze = list(version = RM$spec$version, frozen_sha256 = RM$spec$expected, sha256 = RM$spec$actual, verified = RM$spec$ok),
  created_utc = utc_start, finished_utc = utc_end, wall_seconds = wall,
  phase_seconds = as.list(run_timing),
  timing_note = "components partition the in-R wall time (no overlap); per replicate in timing.csv",
  memory = list(peak_rss_kb = peak_rss_kb, peak_rss_source = "/proc/self/status VmHWM (Linux)",
                r_heap_max_used_mb = gc_max_mb),
  git = list(head = git_out(c("rev-parse", "HEAD"))[1], branch = git_out(c("rev-parse", "--abbrev-ref", "HEAD"))[1],
             engine_status_porcelain = git_out(c("status", "--porcelain", "--", "simulation/engine")),
             package_commit = if (file.exists(file.path(repo_root, "PACKAGE-COMMIT.txt")))
               readLines(file.path(repo_root, "PACKAGE-COMMIT.txt"), warn = FALSE)[1] else NA_character_,
             note = "engine files may be uncommitted; their sha256s below identify what ran"),
  spec = list(path = "letter/SIM-DESIGN-SPEC.md", sha256 = sha(file.path(repo_root, "letter", "SIM-DESIGN-SPEC.md"))),
  engine_files_sha256 = engine_hashes,
  renv_lock = list(sha256 = sha(file.path(engine_dir, "renv.lock")), R_locked = lock$R$Version,
                   n_packages = length(lock_ver)),
  runtime = list(R = R.version.string, platform = R.version$platform, os = utils::sessionInfo()$running,
                 in_wrapper = via_wrapper, libPaths = .libPaths(), threads = getOption("sim.engine.threads"),
                 RNGkind_at_start = rng_kind_start, La_version = La_version(),
                 double_eps = .Machine$double.eps, loaded_packages = pkg_ver),
  config = list(root = ROOT, master_id = MASTER_ID, master_dims = c(H = MASTER_H, J = MASTER_J),
                a_cpsu = A_CPSU, a_cpsu_status = "frozen by R1 (spec v9.8 §12.3)",
                cell = CELL, cell_encoded = rng_encode_cell(CELL), B = B, R = R, task = RM$task, block = RM$block,
                population_ids = POP_IDS, rep_ids = if (length(REP_IDS) <= 50L) REP_IDS else c(first = REP_IDS[1], last = REP_IDS[length(REP_IDS)], n = length(REP_IDS)),
                no_rd = NO_RD, rhat_forest_all_reps = RHAT_FOREST_ALL, n0 = N0, V = V,
                fold_structure_id = FOLD_STRUCTURE_PRIMARY,
                regime = REGIME, e2c = E2C, e3 = E3, path = PATH_LABEL, arms = ARMS, V_run = V,
                flags = list(oracle_only = ORACLE_ONLY, calib = CALIB, u6 = U6, u10 = U10, wr = WR),
                round2 = list(
                  rulings = "R1-R24 DECIDED 2026-10-11 (simulation/AUTHOR-DECISION-LOG-2026-10-04.md §8); spec v9.8",
                  wr = if (WR) list(sampler = "draw_sample_wr: WR first stage, stage 2 keyed by selection (stage2wr; RNG draft §7, PROPOSED)",
                                    target = "sigma2_m,WR = sum_h [(1 - 1/M_pop,h) B_h + V_h] (R14)",
                                    gate = "census-centred linearised score, known N (lin_var_wr); Hajek reported") else NULL,
                  u10 = if (U10) list(cell_label = "e2c = u10 (own cell at the panel coordinates, R10)",
                                      estimator = "oracle Hajek with stage-1 weights w1 = M_pop,h/m_h (C-U10-1)",
                                      gate = "E[psi_U10] = Psi_U10, E[se2_U10] = sigma2_m,U10 (exact, linear)") else NULL,
                  ctrl = if (E3 == "ctrl") list(host = ARM_RA_AIPW, kappa_shift = PANEL_KAPPA,
                                                ae = if (is.null(AE)) NULL else AE[c("kappa", "T2", "stat", "kappa_max", "buildable", "content_sha256", "spec_sha256", "p0_primary_hash")],
                                                ae_clip = if (is.null(AE)) NULL else AE$clip,
                                                ae_construction = "C-AE-1..4 (R/controls.R), confirmed by R23; kappa_AE frozen (C-AE-3)",
                                                g6 = "boundary reading (R5); 3 MCSE at the planned count by runs/run_prep.R") else NULL,
                  coarse = if (HAS_RC_COARSE) list(arm = ARM_RC_INDCF_COARSE, bins = RC_CPSU_BREAKS, codes = RC_CPSU_MIDS,
                                                   pairing = "unit folds and level-4 seeds of Rc-AIPW-indcf (C-RC-5)",
                                                   description = "a sensitivity analysis that also removes predictive information (R15)") else NULL),
                learner_rc = if (HAS_RC) c(RC_SETTINGS, list(conventions = "C-RC-1..3 (R/learner_rb.R), PROVISIONAL (ruling R18)",
                                                            unit_folds = "TC-E6 (R/estimators.R), PROPOSED (ruling R15)")) else NULL,
                learner_rd = if (HAS_RD) list(
                  name = "R-d (fit_rd, §1.7a): R-a on q_m(alpha) whole PSUs of the training fold",
                  alphas = RD_RUN, alphas_frozen = RD_ALPHAS, alphas_stopped_not_run = setdiff(RD_ALPHAS, RD_RUN), m0 = RD_M0, q0 = RD_Q0, m_grid = RD_M_GRID,
                  q_m = stats::setNames(as.list(vapply(RD_ALPHAS, function(a) rd_q(M, a), integer(1))), vapply(RD_ALPHAS, rd_arm, "")),
                  training_pool = M * (V - 1) / V, no_clipping = "q_m > pool is a hard error (§1.7a); no fallback",
                  conventions = "C-RD-1..4 (R/learner_rd.R), PROPOSED level-4 stream; R-a provisional (OI-RA-1..4)",
                  rates = "not measured or claimed (G5 needs the full ladder, §12.5)") else NULL,
                learner_rb = if (!HAS_RB) NULL else list(name = "R-b (fit_rb, §12.11)", num.trees = RB_NUM_TREES, mtry = RB_MTRY,
                                  min.node.size = RB_MIN_NODE, replace = TRUE, sample.fraction = 1,
                                  respect.unordered.factors = "order", features = RB_FEATURES,
                                  conventions = "C-RB-1..8 (R/learner_rb.R): Q by arm, regression forests and shared AIPW/TMLE fits (C-RB-1/2/6, author ruling 2026-10-07; DEC-2 grading open); KB factor, unweighted, PROPOSED level-4 seeds (provisional); explicit training seed at prediction (C-RB-8, OI-RB-3)",
                                  forest_remainders = "not computed (Amendment A1)"),
                quad_rules = list(file = QR$manifest$file, sha256 = QR$sha256, status = QR$manifest$status,
                                  required = TRUE, built_on = QR$manifest$built_on,
                                  refine = list(file = QRR$manifest$file, sha256 = QRR$sha256, status = QRR$manifest$status)),
                t2_refine = t2_policy_record(T2P, QRR),
                rd_alphas_run = RD_RUN,
                tmle = list(clever_covariate = "H = (A/g, (1-A)/(1-g)), 2 columns (TC-T1)",
                            targeting = "stats::glm.fit quasibinomial, offset logit Q_init, weights w/mean(w), epsilon 1e-14 (TC-T2)",
                            conformance = "independent damped Newton refit; Q* agrees to 1e-8 (§12.9)",
                            c_H = TMLE_C_H, pairing = "shares Ra-AIPW's sample, folds and fits (TC-T3)"),
                learner = list(name = "R-a (fit_ra)", a_cpsu_term = a_term, nAGQ = RA_NAGQ,
                               gh_nodes = RA_GH_NODES, q0 = RA_Q_CLIP, delta = RA_G_TRUNC),
                u9_bruteforce = list(k = U9_BF, selected = u9_sel), rng_version = RNG_DRAFT_VERSION),
  keys = keys,
  grids = list(primary_hash = ctx$primary_hash, second_hash = ctx$second_hash,
               primary_points = nrow(ctx$primary), second_points = nrow(ctx$second),
               primary_definition = attr(ctx$primary, "definition"),
               second_definition = attr(ctx$second, "definition"), Psi = Psi_grid),
  checks = list(n = nrow(checks_df), n_required = n_req, n_required_failed = n_req_fail,
                exit_status = status),
  outputs = c("results.csv", "census.csv", "census_strata.csv", "checks.csv", "timing.csv", "t2_refinement.csv", "manifest.json", "SUMMARY.md",
              if (nrow(panel_df)) "panel_predictions.csv"))
jsonlite::write_json(manifest, file.path(out_dir, "manifest.json"), auto_unbox = TRUE, pretty = TRUE,
                     digits = I(17), null = "null", na = "null")

# ---- SUMMARY.md -----------------------------------------------------------------
f6 <- function(x) formatC(x, format = "f", digits = 6)
fe <- function(x) formatC(x, format = "e", digits = 2)
md <- c(
  sprintf("# %s end-to-end run: %s", c(dev = "DEV", dry = "DRY-RUN", prod = "PRODUCTION")[[MODE]], stamp), "",
  sprintf("> **%s.**", RUN_LABEL),
  sprintf("> Spec %s (frozen SHA-256 %s; this file's %s; verified: %s). Adopted RNG contract (`%s`), root `%s`; a_cpsu = %.2f (R1).",
          RM$spec$version, substr(RM$spec$expected, 1, 12), substr(RM$spec$actual, 1, 12), RM$spec$ok, RNG_DRAFT_VERSION, ROOT, A_CPSU),
  sprintf("> %d population(s) (IDs %s) x %d replicate(s) (IDs %s)%s. Nothing below is a grade: grading is the frozen §6 plan's post-processing of production output.",
          B, paste(unique(range(POP_IDS)), collapse = ".."), R, paste(unique(range(REP_IDS)), collapse = ".."),
          if (RM$task > 0L) sprintf(", array task %d (block %d)", RM$task, RM$block) else ""), "",
  "## Configuration", "",
  sprintf("- Cell: %s, H = %d, m_h = %d, J = %d, f = %s, m = %d, V = %d, n0 = %d (`%s`).", PATH_LABEL,
          H, M_H, J, CELL$f, M, V, N0, rng_encode_cell(CELL)),
  sprintf("- Populations: the (%d, %d) block of master `%s` (%d x %d), built directly (nesting checked on a sub-block); population IDs %s; %d draws each.", H, J, MASTER_ID, MASTER_H, if (E2C %in% c("fpc", "wr")) MASTER_J_FULL else MASTER_J, paste(unique(range(POP_IDS)), collapse = ".."), R),
  if (E2C == "fpc") paste0("- Arms: ", paste0("`", ARMS, "`", collapse = ", "),
    " (§12.1 :1146: the oracle at every fraction; fitted R-a and R-b PSU-CF AIPW at f = 0.50 only, with R-a and R-b exactly as in the E1 path; no TMLE). The UC/no-FPC se2 and the two-stage WOR se2 (both stages' fpc declared) are saved for every arm.") else
  paste0("- Arms: `oracle` (eta0 known), `Ra-AIPW` (PSU-CF, V = 2, R-a per K_B class, nAGQ = 9, GH64 marginal, clip q0 = 0.005, truncation delta = 0.01), `Ra-TMLE` (one pooled 2-column logistic fluctuation of Ra-AIPW's clipped held-out fits, no post-clip; §12.9), ",
         if (REGIME == "many") "`Rb-AIPW` (PSU-CF, the §12.11 forest: ranger, 500 trees, mtry 3, min.node.size 25; Q by treatment group, regression forests; PROPOSED level-4 seeds; explicit prediction seed) and `Rb-TMLE` (targets Rb-AIPW's fits). R-b conventions: author ruling 2026-10-07; DEC-2 grading open. No forest T1/T2 (Amendment A1)."
         else if (!HAS_RD) "and no R-d rung in this run (R-d's grid is m = 200 ... 3200, §1.7a, or --no-rd 1). No R-b in fixed-strata (Amendment A1)."
         else sprintf("and the R-d AIPW rung(s) %s (§1.7a: R-a, as above, fitted on q_m(alpha) = floor(30 (m/200)^(2 alpha)) whole PSUs drawn from each training fold by the level-4 stream of the adopted RNG contract; here q_m = %s of a %d-PSU pool; no clipping, no fallback)%s. No R-b in fixed-strata (Amendment A1).",
                      paste0("`", vapply(RD_RUN, rd_arm, ""), "`", collapse = ", "), paste(vapply(RD_RUN, function(a) rd_q(M, a), integer(1)), collapse = "/"), as.integer(M * (V - 1) / V),
                      if (length(RD_RUN) < length(RD_ALPHAS)) sprintf("; not run: %s (STOPPED after G5, ruling R27)", paste0("`", vapply(setdiff(RD_ALPHAS, RD_RUN), rd_arm, ""), "`", collapse = ", ")) else "")),
  sprintf("- Quadrature: every base rule from the saved artifact `%s` (sha256 `%s`, required store, %s).", QR$manifest$file, substr(QR$sha256, 1, 12), QR$manifest$status),
  sprintf("- Wall time in R: %.1f s (by component below). Peak resident memory of the R process: %s; R heap peak %.0f MB.",
          wall, if (is.na(peak_rss_kb)) "not available on this OS" else sprintf("%.0f MB", peak_rss_kb / 1024), gc_max_mb),
  sprintf("- P0 grid: primary %d points (`%s`), second %d points; Psi(P0, H = %d) = %s.", nrow(ctx$primary),
          substr(ctx$primary_hash, 1, 12), nrow(ctx$second), H, formatC(Psi_grid, format = "f", digits = 10)), "",
  "## Checks", "")
fam <- sub("^(U[0-9]+|finite outputs|prefix view|nesting|R4 FAILED arm|AE frozen record|ambient RNG|psi_CF|population|replicate|TMLE premise|TMLE|Rb arm|Rb predictions|P0 grid hash|census FPC algebra|G2 feasibility|Rc predictions \\(coarsened C_psu\\)|Rc predictions|Rc arm|Ra no-CF arm|C3a level|C3a target|Rd arm|Rd throttle|Rd shared|Rd subsamples|G3\\(a-ii\\)|G3\\(a-i\\)|G3\\(a\\)|G3\\(b\\)|G3\\(c\\)|WR linearised score|R14 WR code-check \\(Hajek, report\\)|R14 WR code-check|R14 mismatch|AE T2|AE clip bound|AE map|AE control|AE arm|C2 wrong-clustering|panel G6 limb \\(a\\)|R15 pairing|R15 coarsening).*$", "\\1", checks_df$check)
tab <- do.call(rbind, lapply(split(seq_len(nrow(checks_df)), list(fam, checks_df$required), drop = TRUE), function(i)
  data.frame(family = fam[i[1]], required = checks_df$required[i[1]], rows = length(i), passed = sum(checks_df$pass[i]))))
tab <- tab[order(!tab$required, tab$family), ]
md <- c(md, "| check family | required | rows | passed |", "|---|---|---|---|",
        sprintf("| %s | %s | %d | %d |", tab$family, ifelse(tab$required, "yes", "no (info)"), tab$rows, tab$passed), "",
        sprintf("**Required checks: %d of %d passed. Exit status %d (%s).**", n_req - n_req_fail, n_req, status,
                if (status == 0L) "all required checks pass" else "at least one required check FAILED"))
if (n_req_fail) md <- c(md, "", "Failed required checks:", "",
                        sprintf("- pop %s, rep %s, %s: %s", checks_df$pop_id, checks_df$rep_id, checks_df$arm,
                                checks_df$check)[checks_df$required & !checks_df$pass])
u7r <- checks_df[grepl("^U7", checks_df$check), ]
if (nrow(u7r)) md <- c(md, "", sprintf("U7: the largest |LHS - RHS| was %s against tolerances from %s to %s (64 eps max(|psi_oracle|, 1), R17).",
                                      fe(max(u7r$abs_diff)), fe(min(u7r$tolerance)), fe(max(u7r$tolerance))))
u8r <- checks_df[grepl("^U8 se2", checks_df$check), ]
if (nrow(u8r)) md <- c(md, sprintf("U8: the largest relative se2 gap to survey::svymean was %s (tolerance 1e-12).",
                                   fe(max(u8r$abs_diff / abs(u8r$reference)))))
# ---- failures, runtime by component and memory ----
n_err <- function(pat) sum(grepl(pat, checks_df$check) & !checks_df$pass)
ra_d <- if (nrow(results_df)) results_df[results_df$arm == ARM_RA_AIPW, ] else results_df
rb_d <- if (nrow(results_df)) results_df[results_df$arm == ARM_RB_AIPW, ] else results_df
rd_d <- if (nrow(results_df)) results_df[startsWith(results_df$arm, "Rd-AIPW"), ] else results_df
n_failed_arm <- if (nrow(results_df) && "arm_status" %in% names(results_df)) table(results_df$arm[results_df$arm_status %in% "FAILED"]) else table(character(0))
md <- c(md, "", "## Failures", "",
  "| what | count |", "|---|---|",
  sprintf("| required checks failed | %d of %d |", n_req_fail, n_req),
  sprintf("| R4 FAILED arms (R-a fit failed after one Nelder-Mead refit; counted, never dropped) | %s |",
          if (length(n_failed_arm)) paste(sprintf("%s %d", names(n_failed_arm), as.integer(n_failed_arm)), collapse = "; ") else "none"),
  sprintf("| population errors | %d of %d populations |", n_err("^population [0-9]+ completed"), B),
  sprintf("| replicate errors | %d of %d replicates |", n_err("^replicate \\("), B * R),
  if (HAS_RB) sprintf("| R-b arm errors (learner failures) | %d of %d replicates |", n_err("^Rb arm"), B * R),
  if (HAS_RC) sprintf("| R-c arm errors (learner failures) | %d of %d arm-replicates |", n_err("^Rc arm"),
                      (HAS_RC_PSU + HAS_RC_NOCF + HAS_RC_IND + HAS_RC_COARSE) * B * R),
  if (HAS_RA_NOCF) sprintf("| R-a no-CF arm errors | %d of %d replicates |", n_err("^Ra no-CF arm"), B * R),
  if (U6) sprintf("| U6 limbs failed | %d of %d |", n_err("^U6"), sum(grepl("^U6", checks_df$check))),
  if (CALIB && !WR) { cal_i <- grepl("^U2 |^U8 (linearised|(two-stage|negative|single-stage) limb)", checks_df$check) & checks_df$required
    sprintf("| calibration gate (U2; U8 limbs, linearised known-N) failed | %d of %d |", sum(cal_i & !checks_df$pass), sum(cal_i)) },
  if (HAS_RD) sprintf("| R-d arm errors (learner failures or a throttle hard error) | %d of %d arm-replicates |", n_err("^Rd arm"), 4L * B * R),
  if (CALIB && WR) { wr_i <- grepl("^U2 \\(WR form\\)|^R14 WR code-check \\(linearised", checks_df$check) & checks_df$required
    sprintf("| R14 WR calibration gate (U2 WR form; linearised known-N WR se2) failed | %d of %d |", sum(wr_i & !checks_df$pass), sum(wr_i)) },
  if (U10) { u_i <- grepl("^U10 \\(exact", checks_df$check)
    sprintf("| U10 identities (exact, |z| <= 3) failed | %d of %d |", sum(u_i & !checks_df$pass), sum(u_i)) },
  if (E3 == "ctrl") sprintf("| AE arm errors | %d of %d replicates |", n_err("^AE arm"), B * R),
  if (E3 == "ctrl") sprintf("| C2 wrong-clustering cross-checks failed | %d of %d |", n_err("^C2 wrong-clustering"), sum(grepl("^C2 wrong-clustering", checks_df$check))),
  if (HAS_RC_COARSE) sprintf("| R15 pairing / coarsening checks failed | %d of %d |", n_err("^R15 "), sum(grepl("^R15 ", checks_df$check))),
  if (HAS_RA) sprintf("| R-a glmer fits singular / not clean / with warnings | %s |", if (nrow(ra_d))
    sprintf("%d / %d / %d of %d", sum(ra_d$n_singular, na.rm = TRUE), sum(ra_d$n_not_clean, na.rm = TRUE), sum(ra_d$n_warnings, na.rm = TRUE), sum(ra_d$n_fits, na.rm = TRUE)) else "none run"),
  if (HAS_RB) sprintf("| R-b forests with a ranger warning | %s |", if (nrow(rb_d))
    sprintf("%d of %d", sum(rb_d$n_not_clean, na.rm = TRUE), sum(rb_d$n_fits, na.rm = TRUE)) else "none run"),
  if (HAS_RD) sprintf("| R-d glmer fits singular / not clean / with warnings | %s |", if (nrow(rd_d))
    sprintf("%d / %d / %d of %d", sum(rd_d$n_singular, na.rm = TRUE), sum(rd_d$n_not_clean, na.rm = TRUE), sum(rd_d$n_warnings, na.rm = TRUE), sum(rd_d$n_fits, na.rm = TRUE)) else "none run"),
  if (HAS_RD) sprintf("| R-d throttle / shared-fold checks failed | %d of %d / %d of %d |",
    n_err("^Rd throttle"), sum(grepl("^Rd throttle", checks_df$check)), n_err("^Rd shared"), sum(grepl("^Rd shared", checks_df$check))),
  sprintf("| ambient RNG moved (replicates; whole run) | %d; %s |", n_err("^ambient RNG .* replicate"),
          if (rng_ok) "no" else "YES"),
  if (HAS_RB) sprintf("| R-b predictions not identical() to seed = NULL | %d of %d comparisons |", n_err("^Rb predictions"),
          sum(grepl("^Rb predictions", checks_df$check))))
lab <- c(p0_context = "P0 grid construction (grids + eta0 on them)",
         u9_grid = "U9 verification: P0-grid limb",
         u9_bruteforce = "U9 verification: brute-force P0 f_v",
         populations = sprintf("populations ((%d, %d) block of the %d-stratum master, sub-block nesting check, census%s)", H, J, MASTER_H,
                               if (E2C == "fpc") ", per-stratum FPC components" else ""),
         sampling = "sampling, folds and the oracle",
         ra_fits = "mixed-model fits (R-a glmer, with scoring)",
         rb_forests = "forests (R-b ranger, with scoring)",
         rc_forests = "forests (R-c full-depth ranger: PSU-CF, no-CF, individual-CF, with scoring)",
         rd_fits = "R-d fits (R-a on the throttled subsample, 4 rungs, with scoring)",
         remainders = if (HAS_RD) "remainder evaluation (T2 closed form, T1; Ra-AIPW, Ra-TMLE and the four R-d rungs)" else
           if (HAS_RA_TMLE) "remainder evaluation (T2 closed form, T1; Ra-AIPW and Ra-TMLE)" else "remainder evaluation (T2 closed form, T1; Ra-AIPW)",
         targeting = "targeting (TMLE fluctuation)",
         rhat = "r_hat: census L2 errors of every fold fit (R16)",
         checks = if (HAS_FITTED) "checks (U1, U3-U8, WOR SE, §12.9 conformance refit, RNG and prediction checks, rows)" else
           "checks (U1, U8, the two-stage WOR SE, RNG, rows)",
         outputs = "outputs (CSV writes)", other = "other (start-up, rule loading, configuration)")
ord <- c(if (HAS_RA || HAS_RA_NOCF) "ra_fits", if (HAS_RB) "rb_forests", if (HAS_RC) "rc_forests", if (HAS_RD) "rd_fits",
         if (HAS_RA || HAS_RD) "remainders",
         "u9_grid", if (HAS_RA || HAS_RD) "u9_bruteforce", if (HAS_RA_TMLE || HAS_RB_TMLE) "targeting", if (HAS_FITTED) "rhat",
         "p0_context", "populations", "sampling", "checks", "outputs", "other")
nrep <- max(1L, nrow(timing_df))
md <- c(md, "", "## Runtime by component", "",
  sprintf("Each second of the in-R wall time (%.1f s) is counted once. Per replicate = component / %d replicates (run-level components shown for scale). U9 verification = the two U9 rows (%.1f s).",
          wall, nrep, run_timing[["u9_grid"]] + run_timing[["u9_bruteforce"]]), "",
  "| component | seconds | share of wall | per replicate |", "|---|---|---|---|",
  sprintf("| %s | %.1f | %.1f%% | %.2f |", lab[ord], run_timing[ord], 100 * run_timing[ord] / wall, run_timing[ord] / nrep),
  sprintf("| **total (in R)** | **%.1f** | 100%% | %.2f |", wall, wall / nrep))
if (E2C == "fpc" && nrow(census_df)) {
  ar <- if (nrow(results_df)) do.call(rbind, lapply(split(results_df, list(results_df$arm, results_df$pop_id), drop = TRUE), function(d)
    data.frame(arm = d$arm[1], pop = d$pop_id[1], n = nrow(d), r_var = mean(d$R_var), r_var_wor = mean(d$R_var_wor),
               wor_over_uc = mean(d$se2_wor_over_se2), stringsAsFactors = FALSE))) else NULL
  md <- c(md, "", sprintf("## E2 FPC sweep: census variance components and FPC (J = %d, f = %s)", J, F_STR), "",
    "> **Implementation check only.** These are exact census quantities of each population (G1, :200-201) and per-sample estimates from 5 draws. No coverage grade and no theorem claim: C3a/C3b are graded on >= 100 populations (§8), G1's 3·MCSE adequacy gate needs the planned rep count, and the two-stage WOR estimator's design mean (U8's two-stage limb) needs redraws. The C3b overlay value is census algebra, reported, not graded. Per-stratum components (B_h, V_h, f_h, 1 - f_h, beta_h, ...) are in `census_strata.csv`.", "",
    "| pop | N | f | 1 - f | sigma_m | sigma_m_wor | Delta_m / sigma2_m | sigma_m / sigma_m_wor | pooled beta | mean beta_h | sum f_h V_h / sigma2_m_wor | overlay 0.80 / 0.90 / 0.95 | min n0/M_k |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    sprintf("| %d | %d | %.2f | %.2f | %s | %s | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f / %.4f / %.4f | %.3f |", census_df$pop_id, census_df$N,
            census_df$f, census_df$fpc, f6(census_df$sigma_m), f6(census_df$sigma_m_wor), census_df$Delta_over_sigma2_m,
            census_df$sigma_ratio, census_df$beta_pooled, census_df$mean_beta_h, census_df$deficit_over_sigma2_m_wor,
            census_df$c3b_overlay_080, census_df$c3b_overlay_090, census_df$c3b_overlay_095, census_df$min_pi_unit),
    if (!is.null(ar)) c("", "Per arm and population, 5-draw means (descriptive): R_var = se2 / sigma2_m (UC, no FPC); R_var_wor = se2_wor / sigma2_m_wor (two-stage WOR, both fpc declared); se2_wor / se2.", "",
      "| arm | pop | draws | mean R_var | mean R_var_wor | mean se2_wor / se2 |", "|---|---|---|---|---|---|",
      sprintf("| %s | %d | %d | %.3f | %.3f | %.4f |", ar$arm, ar$pop, ar$n, ar$r_var, ar$r_var_wor, ar$wor_over_uc)))
}
if (nrow(census_df)) {
  md <- c(md, "", "## Census targets (per population)", "",
          "| pop | N | Psi_N | Psi_N - Psi | sigma_m | sigma_m_wor | Delta_m / sigma_m^2 | min g0 | max g0 | min Q | max Q |",
          "|---|---|---|---|---|---|---|---|---|---|---|",
          sprintf("| %d | %d | %s | %s | %s | %s | %s | %.4f | %.4f | %.4f | %.4f |", census_df$pop_id, census_df$N,
                  f6(census_df$Psi_N), fe(census_df$Psi_N_minus_Psi), f6(census_df$sigma_m), f6(census_df$sigma_m_wor),
                  formatC(census_df$Delta_over_sigma2_m, format = "f", digits = 4),
                  census_df$min_g0, census_df$max_g0, census_df$min_Q, census_df$max_Q))
}
ok_df <- if (nrow(results_df) && "arm_status" %in% names(results_df)) results_df[results_df$arm_status %in% "OK", , drop = FALSE] else results_df
if (nrow(ok_df)) {
  results_df_all <- results_df; results_df <- ok_df          # per-arm summaries over the arms that ran (FAILED counted above)
  per <- do.call(rbind, lapply(split(results_df, list(results_df$arm, results_df$pop_id), drop = TRUE), function(d)
    data.frame(arm = d$arm[1], pop = d$pop_id[1], n = nrow(d), bias = mean(d$psi_hat - d$Psi_N),
               mean_se = mean(d$se_hat), sd_psi = stats::sd(d$psi_hat), sigma_m = d$sigma_m[1],
               sigma_wor = d$sigma_m_wor[1], cover = sum(d$cover95), mean_Z = mean(d$Z), sd_Z = stats::sd(d$Z),
               R_var = mean(d$R_var), stringsAsFactors = FALSE)))
  per <- per[order(per$arm != ARM_ORACLE, per$pop), ]
  md <- c(md, "", "## Arms, per population (descriptive only)", "",
          "SD(psi) is the SD of psi_hat over that population's draws. mean SE / SD compares the average estimated SE with that SD.", "",
          "| arm | pop | draws | mean(psi - Psi_N) | mean SE | SD(psi) | mean SE / SD | sigma_m | sigma_m_wor | mean Z | SD(Z) | cover95 | mean R_var |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|",
          sprintf("| %s | %d | %d | %s | %s | %s | %.3f | %s | %s | %+.3f | %.3f | %d/%d | %.3f |", per$arm, per$pop, per$n,
                  fe(per$bias), f6(per$mean_se), f6(per$sd_psi), per$mean_se / per$sd_psi, f6(per$sigma_m),
                  f6(per$sigma_wor), per$mean_Z, per$sd_Z, per$cover, per$n, per$R_var))
  pooled <- do.call(rbind, lapply(split(results_df, results_df$arm), function(d) data.frame(
    arm = d$arm[1], n = nrow(d), bias = mean(d$psi_hat - d$Psi_N),
    bias_mcse = stats::sd(d$psi_hat - d$Psi_N) / sqrt(nrow(d)), mean_se = mean(d$se_hat),
    sd_within = mean(tapply(d$psi_hat, d$pop_id, stats::sd)), cover = sum(d$cover95))))
  pooled <- pooled[order(pooled$arm != ARM_ORACLE), ]
  md <- c(md, "", "Pooled over populations (SD within population, averaged):", "",
          "| arm | draws | mean(psi - Psi_N) | naive MC SE | mean SE | mean within-pop SD(psi) | cover95 |",
          "|---|---|---|---|---|---|---|",
          sprintf("| %s | %d | %s | %s | %s | %s | %d/%d |", pooled$arm, pooled$n, fe(pooled$bias), fe(pooled$bias_mcse),
                  f6(pooled$mean_se), f6(pooled$sd_within), pooled$cover, pooled$n))
  ra <- results_df[results_df$arm == ARM_RA_AIPW, ]
  if (nrow(ra)) md <- c(md, "", "## R-a remainders and fits (descriptive only)", "",
    sprintf("- T2 (closed form; frozen primary grid, refined where a fold triggered -- round 3d): mean %s, range [%s, %s]; mean |T2| / sigma_m_wor = %.4f.",
            fe(mean(ra$T2, na.rm = TRUE)), fe(suppressWarnings(min(ra$T2, na.rm = TRUE))), fe(suppressWarnings(max(ra$T2, na.rm = TRUE))), mean(abs(ra$T2_over_sigma_wor), na.rm = TRUE)),
    sprintf("- T2 refinement (PROPOSED, round 3d): %d of %d replicates had a triggered fold; refined %d, unresolved %d (T1, T2 recorded NA; rungs in t2_refinement.csv).",
            sum(ra$T2_folds_triggered > 0, na.rm = TRUE), nrow(ra), sum(ra$T2_status %in% "refined"), sum(ra$T2_status %in% "unresolved")),
    sprintf("- T1 (by subtraction): mean %s, range [%s, %s]; mean |T1| / sigma_m_wor = %.4f.",
            fe(mean(ra$T1, na.rm = TRUE)), fe(suppressWarnings(min(ra$T1, na.rm = TRUE))), fe(suppressWarnings(max(ra$T1, na.rm = TRUE))), mean(abs(ra$T1_over_sigma_wor), na.rm = TRUE)),
    sprintf("- Paired psi_CF - psi_oracle = T1 + T2: mean %s.", fe(mean(ra$T1 + ra$T2, na.rm = TRUE))),
    sprintf("- glmer fits: %d (4 per replicate: Q and g per K_B class, each fold). Singular: %d. Not clean (warning or lme4 message): %d. Warnings: %d. Every fit was kept (C-RA-3).",
            sum(ra$n_fits), sum(ra$n_singular), sum(ra$n_not_clean), sum(ra$n_warnings)),
    sprintf("- Scored units with a clipped or truncated nuisance: %d of %d.", sum(ra$clip_hits), sum(ra$n_units)),
    sprintf("- Fit time: %.1f s in total, %.2f s per replicate on average.", sum(ra$fit_seconds), mean(ra$fit_seconds)))
}
for (tmle_arm in c(ARM_RA_TMLE, ARM_RB_TMLE)) {
tmr <- if (nrow(results_df)) results_df[results_df$arm == tmle_arm, ] else results_df
if (nrow(tmr)) md <- c(md, "", sprintf("## %s targeting and the paired TMLE - AIPW objects (descriptive only)", tmle_arm), "",
  sprintf("- eps (pooled, 2 components): eps_1 in [%s, %s], eps_0 in [%s, %s]; m^(1/4)||eps|| mean %.4f, max %.4f.",
          fe(min(tmr$eps_1)), fe(max(tmr$eps_1)), fe(min(tmr$eps_0)), fe(max(tmr$eps_0)),
          mean(tmr$m14_eps_norm), max(tmr$m14_eps_norm)),
  sprintf("- §12.9 conformance: largest |Q* - independent refit| %s (tolerance 1e-8); largest |score| %s (tolerance 1e-9); glm.fit iterations %d to %d.",
          fe(max(tmr$conformance_max_abs)), fe(max(tmr$score_max_abs)), min(tmr$tmle_iter), max(tmr$tmle_iter)),
  sprintf("- Gram premise (T2): lambda_min %.4f to %.4f; met (>= c_H = %.2f) in %d of %d replicates. A breach makes C6 UNGRADABLE (premise); it does not fail the build.",
          min(tmr$lambda_min), max(tmr$lambda_min), TMLE_C_H, sum(tmr$gram_premise_met), nrow(tmr)),
  sprintf("- Paired psi_TMLE - psi_AIPW: mean %s, range [%s, %s]; sqrt(m)|diff| max %.4f. (se2_TMLE - se2_AIPW)/sigma2_m mean %+.4f; m^(1/4)|.| max %.4f.",
          fe(mean(tmr$paired_tmle_minus_aipw)), fe(min(tmr$paired_tmle_minus_aipw)), fe(max(tmr$paired_tmle_minus_aipw)),
          max(tmr$sqrt_m_abs_paired), mean(tmr$se2_diff_over_sigma2_m), max(tmr$m14_abs_se2_diff_over_sigma2_m)),
  if (tmle_arm == ARM_RB_TMLE) "- T2^TMLE / T1^TMLE: not computed for a forest arm (Amendment A1)." else
    sprintf("- T2^TMLE (targeted fold fits; refined where the AIPW fold triggered): mean %s; T1^TMLE by subtraction: mean %s; unresolved %d.",
            fe(mean(tmr$T2, na.rm = TRUE)), fe(mean(tmr$T1, na.rm = TRUE)), sum(tmr$T2_status %in% "unresolved")),
  sprintf("- Targeting time: %.2f s in total, %.3f s per replicate on average.", sum(tmr$tmle_seconds), mean(tmr$tmle_seconds)))
}
rbr <- if (nrow(results_df)) results_df[results_df$arm == ARM_RB_AIPW, ] else results_df
rb_fail <- checks_df[grepl("^Rb arm", checks_df$check), , drop = FALSE]
if (HAS_RD) {
  rd_tab <- do.call(rbind, lapply(RD_RUN, function(a) {
    ak <- rd_arm(a); d <- results_df[results_df$arm == ak, , drop = FALSE]
    ck <- checks_df[checks_df$arm == ak, , drop = FALSE]
    u9 <- ck[grepl("^U9 brute-force", ck$check), , drop = FALSE]
    data.frame(arm = ak, alpha = a, q = rd_q(M, a), share = rd_q(M, a) / (M * (V - 1) / V), n = nrow(d),
               req = sum(ck$required), req_fail = sum(ck$required & !ck$pass),
               errors = sum(grepl("^Rd arm", checks_df$check) & checks_df$arm == ak & !checks_df$pass),
               singular = if (nrow(d)) sum(d$n_singular) else NA, fits = if (nrow(d)) sum(d$n_fits) else NA,
               T2 = if (nrow(d)) mean(d$T2, na.rm = TRUE) else NA, T2_lo = if (nrow(d)) suppressWarnings(min(d$T2, na.rm = TRUE)) else NA, T2_hi = if (nrow(d)) suppressWarnings(max(d$T2, na.rm = TRUE)) else NA,
               aT2 = if (nrow(d)) mean(abs(d$T2_over_sigma_wor), na.rm = TRUE) else NA, T1 = if (nrow(d)) mean(d$T1, na.rm = TRUE) else NA,
               aT1 = if (nrow(d)) mean(abs(d$T1_over_sigma_wor), na.rm = TRUE) else NA,
               u7 = if (nrow(d)) max(d$U7_abs_diff) else NA, u9 = if (nrow(u9)) max(u9$abs_diff) else NA, n_u9 = nrow(u9),
               stringsAsFactors = FALSE)
  }))
  md <- c(md, "", "## R-d throttled rungs (implementation verification; descriptive only)", "",
    "> No nuisance rate is measured or claimed here, and no remainder is graded: G5 (§12.5) needs the full ladder m in {200, 400, 800, 1600, 3200}, and C4/C5 are graded only on it. The remainder numbers below show the retained calculations run and agree with their identities.", "",
    sprintf("Training pool %d PSUs per fold (m(V-1)/V). Every R-d check family is in the table above; U7 is the P0-free fold identity (tolerance 64 eps |psi_oracle|) and U9 the brute-force P0 f_v cross-check of the closed-form T2 (tolerance 1e-6).", as.integer(M * (V - 1) / V)), "",
    "| arm | alpha | q_m | share of pool | draws | required failed | errors | glmer singular | mean T2 | T2 range | mean abs(T2)/sigma_wor | mean T1 | mean abs(T1)/sigma_wor | max U7 gap | U9 checks | max U9 gap |",
    "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    sprintf("| %s | %.2f | %d | %.3f | %d | %d of %d | %d | %s | %s | [%s, %s] | %.4f | %s | %.4f | %s | %d | %s |",
            rd_tab$arm, rd_tab$alpha, rd_tab$q, rd_tab$share, rd_tab$n, rd_tab$req_fail, rd_tab$req, rd_tab$errors,
            ifelse(is.na(rd_tab$fits), "-", sprintf("%d of %d", rd_tab$singular, rd_tab$fits)),
            fe(rd_tab$T2), fe(rd_tab$T2_lo), fe(rd_tab$T2_hi), rd_tab$aT2, fe(rd_tab$T1), rd_tab$aT1,
            fe(rd_tab$u7), rd_tab$n_u9, fe(rd_tab$u9)))
}
if (HAS_RB) md <- c(md, "", "## R-b forests (descriptive only)", "",
  sprintf("- Learner failures (an Rb arm that stopped with an error): %d of %d replicates.%s", nrow(rb_fail), B * R,
          if (nrow(rb_fail)) paste0(" First: ", rb_fail$check[1]) else ""),
  if (nrow(rbr)) c(
    sprintf("- Forests: %d (6 per replicate: Q1, Q0 and g in each of 2 folds), each 500 trees. Fits with a ranger warning: %d; warnings: %d. Every fit was kept (C-RB-7).",
            sum(rbr$n_fits), sum(rbr$n_not_clean), sum(rbr$n_warnings)),
    sprintf("- OOB MSE (mean over folds): Q1 %.4f to %.4f, Q0 %.4f to %.4f, g %.4f to %.4f.",
            min(rbr$rb_oob_mse_Q1), max(rbr$rb_oob_mse_Q1), min(rbr$rb_oob_mse_Q0), max(rbr$rb_oob_mse_Q0),
            min(rbr$rb_oob_mse_g), max(rbr$rb_oob_mse_g)),
    sprintf("- Scored units with a clipped or truncated nuisance: %d of %d.", sum(rbr$clip_hits), sum(rbr$n_units)),
    sprintf("- Forest fit time: %.1f s in total, %.2f s per replicate on average.", sum(rbr$fit_seconds), mean(rbr$fit_seconds))))
if (WR && nrow(census_df)) {
  wrr <- checks_df[grepl("^U2 \\(WR form\\)|^R14 ", checks_df$check), , drop = FALSE]
  md <- c(md, "", "## E2 WR code-check (ruling R14)", "",
    "> Implementation validation only (\"off-theorem\", §3/E2 :655). Target: the literal with-replacement variance sigma2_m,WR = sum_h [(1 - 1/M_pop,h) B_h + V_h]. GATED: the census-centred linearised score with known N, whose design expectation is exactly that target (exact enumeration, tests/testthat/test-sampling-wr.R). REPORTED: the Hajek WR estimator, and the WOR-defined sigma2_m that correct WR code misses by sum_h B_h / M_pop,h.", "",
    "| pop | sigma2_m (WOR) | sigma2_m,WR | gap sum_h B_h/M_pop,h | gap / sigma2_m |", "|---|---|---|---|---|",
    sprintf("| %d | %.6e | %.6e | %.6e | %.4f |", census_df$pop_id, census_df$sigma2_m, census_df$sigma2_m_wr, census_df$wr_gap,
            census_df$wr_gap / census_df$sigma2_m),
    if (nrow(wrr)) c("", "| pop | check | value | reference | required | pass |", "|---|---|---|---|---|---|",
      sprintf("| %s | %s | %.6g | %.6g | %s | %s |", wrr$pop_id, wrr$check, wrr$value, wrr$reference, wrr$required, wrr$pass)))
}
if (U10 && nrow(census_df)) {
  u1r <- checks_df[grepl("^U10 ", checks_df$check), , drop = FALSE]
  md <- c(md, "", "## U10: stage-2 weights omitted (G8 :348; ruling R10, its own cell)", "",
    "> U10 shows the stage-2 weights are load-bearing: without them the oracle estimates the unweighted mean of PSU means, Psi_U10, not Psi_N. Its exact linear identities are GATED (|z| <= 3, R17 PROPOSED); its Z keyed to Psi_N and the G6 limb (a) boundary distance (R5) are REPORTED -- the 3 MCSE test needs the planned count (R6/R7 OPEN).", "",
    "| pop | Psi_N | Psi_U10 | bias | sigma_m,U10 | predicted mean Z | G6 distance from [-0.10, 0.10] |", "|---|---|---|---|---|---|---|",
    sprintf("| %d | %.6f | %.6f | %+.3e | %.6f | %+.3f | %+.3f |", census_df$pop_id, census_df$Psi_N, census_df$Psi_U10, census_df$u10_bias,
            sqrt(census_df$sigma2_m_u10), census_df$u10_pred_mean_Z, census_df$u10_g6_distance),
    if (nrow(u1r)) c("", "| pop | check | value | reference | required | pass |", "|---|---|---|---|---|---|",
      sprintf("| %s | %s | %.6g | %.6g | %s | %s |", u1r$pop_id, u1r$check, u1r$value, u1r$reference, u1r$required, u1r$pass)))
}
if (E3 == "ctrl" && nrow(panel_df)) {
  ra <- results_df[results_df$arm == ARM_RA_AIPW, , drop = FALSE]
  ae <- results_df[results_df$arm == ARM_AE_AIPW, , drop = FALSE]
  obs <- c(if (nrow(ra)) mean(ra$Z_shift) else NA, if (nrow(ra)) mean(ra$R_var_nomh) else NA, if (nrow(ra)) mean(ra$R_var_unit) else NA,
           if (nrow(ra)) mean(ra$R_var_nostrat) else NA,
           if (nrow(ae)) unname(stats::quantile(abs(ae$T2_over_sigma_wor), 0.9, type = 7)) else NA)
  dev_sd <- c(if (nrow(ra)) stats::sd(ra$Z_shift) else NA, if (nrow(ra)) stats::sd(ra$R_var_nomh) else NA,
              if (nrow(ra)) stats::sd(ra$R_var_unit) else NA, if (nrow(ra)) stats::sd(ra$R_var_nostrat) else NA, NA)
  ctl <- unique(panel_df$control)
  pred <- vapply(ctl, function(cc) mean(panel_df$predicted[panel_df$control == cc]), numeric(1))
  dist <- vapply(ctl, function(cc) mean(panel_df$distance[panel_df$control == cc]), numeric(1))
  stat <- vapply(ctl, function(cc) panel_df$status[panel_df$control == cc][1], "")
  md <- c(md, "", "## E3 power-control panel (rulings R5, R10; descriptive)", "",
    sprintf("> Host: R-a PSU-CF AIPW. Shift constant 0.30 (R10; R8's graded form and symbol stay OPEN, so the g-formula-keyed Z is reported too). Aligned-error control: %s. G6 limb (a) is read from the compliant region's nearest boundary (R5): a control whose census-predicted value lies inside its region fails G6 and is a diagnostic, never graded. The 3 MCSE test needs the planned count and uncertainty procedure (R6/R7 OPEN), so nothing here is graded; the DEV SD is over these B x R = %d replicates.",
            if (is.null(AE)) "not run" else if (AE$buildable) sprintf("kappa_AE = %.5f (kappa_max %.3f), T2 = %.4e, predicted Q_0.9(|T2|/sigma_m,wor) = %.4f on this run's populations", AE$kappa, AE$clip$kappa_max, AE$T2, AE$stat) else paste("NOT BUILDABLE --", AE$reason), B * R), "",
    "| control | region | predicted (mean over pops) | G6 distance | G6 status (R5) | DEV observed | DEV SD |", "|---|---|---|---|---|---|---|",
    sprintf("| %s | [%s, %s] | %.4f | %+.4f | %s | %.4f | %s |", ctl,
            vapply(ctl, function(cc) format(panel_df$lo[panel_df$control == cc][1]), ""),
            vapply(ctl, function(cc) format(panel_df$hi[panel_df$control == cc][1]), ""),
            pred, dist, stat, obs, ifelse(is.na(dev_sd), "-", sprintf("%.4f", dev_sd))),
    if (nrow(ra)) sprintf("- g-formula-keyed Z (reported; R8 OPEN): mean %+.3f over %d replicates; census offsets (Psi_N - Psi_N^gf)/sigma_m = %s.",
                          mean(ra$Z_gf), nrow(ra), paste(sprintf("%+.3f", census_df$Z_gf_pred_offset), collapse = ", ")),
    if (nrow(ae)) sprintf("- AE-AIPW: T2 = %.4e in every replicate (a closed-form constant); |T2|/sigma_m,wor %.4f to %.4f; T1 by subtraction mean %+.3e; U7 max gap %s.",
                          ae$T2[1], min(abs(ae$T2_over_sigma_wor)), max(abs(ae$T2_over_sigma_wor)), mean(ae$T1), fe(max(ae$U7_abs_diff))))
}
if (HAS_RC_COARSE && nrow(results_df)) {
  ci <- results_df[results_df$arm == ARM_RC_INDCF, , drop = FALSE]; cc <- results_df[results_df$arm == ARM_RC_INDCF_COARSE, , drop = FALSE]
  md <- c(md, "", "## Individual-CF leakage-neutral sensitivity (ruling R15; descriptive)", "",
    "> **A sensitivity analysis that also removes predictive information.** The coarsened arm replaces the learner's C_psu by its quartile bin; C_psu is a true confounder of Q0 and g0 (§12.3), so the coarsening removes genuine predictive signal together with the PSU-identifying signal that individual-CF can leak. A difference between the two arms therefore cannot be attributed to leakage alone. Same sample, unit folds and level-4 seeds (C-RC-5); R-c is the PROPOSED full-depth foil (R18 OPEN).", "",
    "| arm | draws | mean(psi - Psi_N) | mean SE | mean R_var | cover95 |", "|---|---|---|---|---|---|",
    sprintf("| %s | %d | %s | %s | %.3f | %d/%d |", c(ARM_RC_INDCF, ARM_RC_INDCF_COARSE), c(nrow(ci), nrow(cc)),
            c(fe(mean(ci$psi_hat - ci$Psi_N)), fe(mean(cc$psi_hat - cc$Psi_N))), c(f6(mean(ci$se_hat)), f6(mean(cc$se_hat))),
            c(mean(ci$R_var), mean(cc$R_var)), c(sum(ci$cover95), sum(cc$cover95)), c(nrow(ci), nrow(cc))),
    if (nrow(cc)) sprintf("- Paired (coarse - indcf), per replicate: psi mean %+.3e (range %+.3e to %+.3e); (se2 difference)/sigma2_m mean %+.4f.",
                          mean(cc$paired_minus_indcf_psi), min(cc$paired_minus_indcf_psi), max(cc$paired_minus_indcf_psi),
                          mean(cc$paired_minus_indcf_se2_over_sigma2_m)))
}
if (HAS_RC || HAS_RA_NOCF || E3 == "vsens" || U6 || (CALIB && !WR)) {
  fam_rows <- function(pat, required_only = FALSE) {
    i <- grepl(pat, checks_df$check) & (!required_only | checks_df$required); c(sum(i), sum(checks_df$pass[i])) }
  ex <- rbind(
    if (U6) c("U6 invariance and sensitivity limbs (Ra-AIPW, required)", fam_rows("^U6")),
    if (HAS_RC_IND) c("individual-CF: U4 unit analogue and U5 false by design (required)", fam_rows("^U4 \\(unit|^U5 \\(unit|^U5 false")),
    if (HAS_RC_PSU) c("R-c prediction identity, explicit seed (required)", fam_rows("^Rc predictions")),
    if (CALIB && !WR) c("U2 and U8 limbs (linearised score, known N) over redraws, |z| <= 3 (required; tolerance PROPOSED, R17)",
                 fam_rows("^U2 |^U8 (linearised|(two-stage|negative|single-stage) limb)|^U8 negative limb: census", required_only = TRUE)),
    if (CALIB && !WR) c("U8 limbs on the actual Hajek estimators (reported, not gated: O(1/m) from the census targets)",
                 fam_rows("\\(Hajek, report\\)")))
  if (HAS_RC_COARSE) ex <- rbind(ex, c("R15 pairing and coarsening (required)", fam_rows("^R15 ")),
                                  c("R-c prediction identity on coarsened C_psu (required)", fam_rows("^Rc predictions \\(coarsened")))
  md <- c(md, "", "## Construction checks: tally", "",
          "> Implementation checks only (synthetic DEV populations). R-c and the calibration tolerances are PROPOSED conventions awaiting rulings R18 and R17; the individual-CF fold rule and its coarsening are DECIDED (R15, 2026-10-11). No grade is attached to any arm.", "",
          "| family | rows | passed |", "|---|---|---|",
          if (!is.null(ex)) sprintf("| %s | %s | %s |", ex[, 1], ex[, 2], ex[, 3]))
  if (CALIB && !WR) {
    cal <- checks_df[grepl("^U2 |^U8 (linearised|(two-stage|negative|single-stage) limb)|^C3a", checks_df$check), , drop = FALSE]
    md <- c(md, "", "Calibration rows (per population; `required` = gated):", "",
            "| pop | check | value | reference | required | pass |", "|---|---|---|---|---|---|",
            sprintf("| %s | %s | %.6g | %s | %s | %s |", cal$pop_id, cal$check, cal$value,
                    ifelse(is.na(cal$reference), "-", sprintf("%.6g", cal$reference)), cal$required, cal$pass))
  }
}
md <- c(md, "", "## Not covered", "",
        if (WR) "Every grade; the WR cell is an off-theorem implementation check (§3/E2 :655). The tolerance |z| <= 3 is PROPOSED (R17). See `docs/E2E-DEV.md`."
        else if (U10) "U10's G6 limb (a) 3 MCSE test and verdict (the planned count and the uncertainty procedure, R6/R7, are OPEN), every grade. See `docs/E2E-DEV.md`."
        else if (E3 == "ctrl") "Every control's verdict: G6 limb (a)'s 3 MCSE test at the planned count (R6/R7 OPEN), §6.1 dispositions, R8 (C1's graded form and symbol). kappa_AE is re-solved on the production panel populations before the freeze. See `docs/E2E-DEV.md`."
        else if (HAS_RC_COARSE) "Any reading of the coarsened arm as a leakage measurement (coarsening also removes predictive information, R15), R-c's settings (R18 OPEN), every grade. See `docs/E2E-DEV.md`."
        else if (HAS_RC || E3 == "vsens") "Any C2/R_var reading of R-c or the V comparison (diagnostic roles only, §3/E3), the leakage-neutral individual-CF sensitivity (ruling R15), the power-control panel (R5, R10), R-c rate audit, every grade. See `docs/E2E-DEV.md`."
        else if (E2C == "strict") "The strict trajectory's census limb (runs/run_e2_census.R), C3a/C3b grading, G1's 3·MCSE gate at planned counts, every grade. See `docs/E2E-DEV.md`."
        else if (E2C == "fpc") "C3a/C3b grading and every coverage statement (>= 100 populations, §8), G1's 3·MCSE adequacy gate at the planned rep count, the U8 two-stage limb's design mean and the U8 negative limb as gates (redraws), the WR code-check cell, the vanishing-fraction, genuinely-strict and unequal-allocation constructions, TMLE in E2, the E2 population-ID allocation (RNG draft §12 item 3, OPEN), production RNG and production populations. See `docs/E2E-DEV.md`."
        else if (REGIME == "many") "R-b's §12.11 rate gates (they need the m-ladder), forest T1/T2 (Amendment A1), R-c/R-d, every other cell, C6 grading (needs the m-sequence and the AIPW cell's C1/C2 grades), U2/U6/U10 gates, the U8 two-stage and negative limbs as gates, production RNG and production populations, and dry-run counts. See `docs/E2E-DEV.md`."
        else "G5 R-d rate validation and any achieved-rate claim (they need the full ladder m in {200, ..., 3200}), C4a/C4b/C5 grading, the R-d m = 400/1600/3200 rungs, the engineered aligned-error control, R-c, the fixed-strata oracle branch's design-CLT grading, U2/U6/U10 gates, the U8 two-stage and negative limbs as gates, production RNG and production populations, and dry-run counts. See `docs/E2E-DEV.md`.")
writeLines(md, file.path(out_dir, "SUMMARY.md"), useBytes = TRUE)

cat(sprintf("[dev-e2e] done in %.1f s: %d checks, %d required, %d required failed; exit %d\n[dev-e2e] wrote %s\n",
            wall, nrow(checks_df), n_req, n_req_fail, status, out_dir))
cat(sprintf("[dev-e2e] timing (s): %s\n", paste(sprintf("%s %.1f", names(run_timing), run_timing), collapse = ", ")))
cat(sprintf("[dev-e2e] peak RSS %s MB; R heap peak %.0f MB\n",
            if (is.na(peak_rss_kb)) "NA" else sprintf("%.0f", peak_rss_kb / 1024), gc_max_mb))
if (status != 0L) quit(save = "no", status = status)
