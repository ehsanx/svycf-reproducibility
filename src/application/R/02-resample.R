# ---------------------------------------------------------------------------
# 02-resample.R -- the SPA census-resampling study, behind the frozen settings in
# config/application/protocol.yaml.
#
# Exactly one mode flag is required:
#   --selftest                 synthetic data only: sampler, inclusion weights, design
#                              covariates, arm-failure handling, metrics, frame coding,
#                              protocol parsing. No real data.
#   --reference --survey <S>   the stabilised census pseudo-reference: each estimator
#                              family on the FULL frame (weights 1), R_ref splits,
#                              median; plus the census overlap diagnostic.
#   --run --survey <S>         R resamples of the imposed two-stage design; every arm on
#                              every resample (arm-level failures kept); per-replicate
#                              weight diagnostics.
#   --summarise --survey <S>   calibration per arm, with denominators.
# --reference and --run refuse to open any data file unless the protocol is FROZEN and
# complete, the installed svytmle matches the pinned version and build hash, and
# the protocol, the gate and both scripts are unchanged (git checkout: committed; copy without
# git: md5 as in config/application/code-manifest.tsv).
#
# Nuisance weighting. The imposed design selects facilities within strata and records
# with probability min(c, N_k)/N_k, so sampling is informative given the substantive
# covariates. PRIMARY arms use svytmle >= 0.1.2 with nuisance_weights = "survey" (the
# cross-fitted and CV nuisances are fit with the survey weights, including the
# training-fold PSU-selection correction) and the substantive C only. SENSITIVITY arms
# keep the package default (unweighted nuisances) and restore ignorability by adding
# the design stratum and log N_k to C; they target a contrast also adjusted for those.
#
# Claims: the variance ratio is primary; bias and coverage relative to the
# census PSEUDO-reference are descriptive. Never truth, coverage of Psi_N, or theorem
# validation. DATA POLICY: record-level frames stay in memory; outputs are aggregates.
# ---------------------------------------------------------------------------

local({   # reuse the count audit's utilities (YAML reader, preflight, readers, strata)
  here <- if (file.exists(file.path("src", "application", "R", "01-count-audit.R")))
    file.path("src", "application", "R") else "."
  sys.source(file.path(here, "01-count-audit.R"), envir = globalenv())
})
suppressPackageStartupMessages({ library(survey); library(svytmle); library(parallel) })
options(survey.lonely.psu = "adjust")

PROTOCOL_REL <- "config/application/protocol.yaml"
ARM_CODES  <- c(fa = "Fully-Aware", pa = "Partially-Aware", na = "Non-Aware", cv = "Fully-Aware-CV", cf = "Fully-Aware-CF")
SENS_CF   <- "Fully-Aware-CF [unw+design C]"          # sensitivity labels
SENS_AIPW <- "AIPW-CF [unw+design C]"
ARM_FAMILY <- c("Fully-Aware" = "tmle_sf", "Partially-Aware" = "tmle_sf", "Non-Aware" = "tmle_sf_unw",
                "Fully-Aware-CV" = "tmle_cv", "Fully-Aware-CF" = "tmle_cf", "AIPW-CF" = "aipw_cf",
                setNames(c("tmle_cf_unwd", "aipw_cf_unwd"), c(SENS_CF, SENS_AIPW)))
ALL_ARMS <- names(ARM_FAMILY)
DESIGN_COVARIATES <- c("design_stratum", "logNk")   # sensitivity arms only (see header)
RES_COLS <- c("method", "b", "se", "df", "lcl", "ucl")

# --- protocol ---------------------------------------------------------------------
validate_protocol <- function(p) {
  bad <- character(0); add <- function(...) bad <<- c(bad, paste0(...))
  if (!identical(p$status, "FROZEN")) add("status must be FROZEN (is ", format(p$status %||% "missing"), ")")
  for (nm in c("frozen_by", "frozen_on", "anchor", "spec_id", "survey")) if (is_placeholder(p[[nm]])) add(nm, " unset")
  d <- p$design %||% list()
  fv <- unlist(d$f)
  if (!is.numeric(fv) || !length(fv) || any(fv <= 0 | fv >= 1) || anyDuplicated(fv)) add("design$f must be distinct values in (0, 1)")
  if (length(fv) > 1 && (is.null(d$primary_f) || !isTRUE(d$primary_f %in% fv))) add("design$primary_f must be one of design$f")
  if (!is.numeric(d$c) || d$c < 1 || d$c != round(d$c)) add("design$c must be an integer >= 1")
  e <- p$estimation %||% list()
  if (!is.numeric(e$R) || e$R < 10) add("estimation$R must be >= 10")
  if (!is.numeric(e$R_ref) || e$R_ref < 1) add("estimation$R_ref must be >= 1")
  if (!is.numeric(e$cf_folds) || e$cf_folds < 2) add("estimation$cf_folds must be >= 2")
  if (!valid_set(e$learners, e$learners %||% "")) add("estimation$learners must be a non-empty list")
  if (!is.numeric(e$seed)) add("estimation$seed must be a number")
  for (nm in c("svytmle_version", "svytmle_packaged")) if (is_placeholder(e[[nm]])) add("estimation$", nm, " unset (pin)")
  if (length(bad)) fail("PROTOCOL NOT USABLE:\n  - ", paste(bad, collapse = "\n  - "))
  p$estimation$learners <- unlist(e$learners)
  p$design$f <- unlist(d$f)
  if (is.null(p$design$primary_f)) p$design$primary_f <- p$design$f[1]
  p
}

# The installed svytmle must be the pinned source build: same version AND the same
# "Packaged" stamp that R CMD build wrote into the vendored tarball
# (vendor/svytmle_0.1.2.tar.gz). Unlike an md5 of the compiled .rdb, this holds across
# platforms (Windows, Linux). The .rdb md5 is recorded, not enforced.
svytmle_build_md5 <- function() unname(tools::md5sum(system.file("R", "svytmle.rdb", package = "svytmle")))
check_pins <- function(prot) {
  e <- prot$estimation
  if (!identical(as.character(packageVersion("svytmle")), as.character(e$svytmle_version)))
    fail("installed svytmle ", as.character(packageVersion("svytmle")), " != pinned ", e$svytmle_version)
  pk <- utils::packageDescription("svytmle")$Packaged %||% "(none: not installed from a built tarball)"
  if (!identical(pk, as.character(e$svytmle_packaged)))
    fail("installed svytmle is not the pinned build (Packaged '", pk, "' != '", e$svytmle_packaged, "')")
  c(R = R.version.string, platform = R.version$platform, svytmle = as.character(packageVersion("svytmle")),
    svytmle_packaged = pk, svytmle_rdb_md5 = svytmle_build_md5(),
    SuperLearner = as.character(packageVersion("SuperLearner")), tmle = as.character(packageVersion("tmle")),
    survey = as.character(packageVersion("survey")), earth = as.character(packageVersion("earth")),
    glmnet = as.character(packageVersion("glmnet")), ranger = as.character(packageVersion("ranger")))
}

# --- analysis frame (anchor-specific coding + mandatory design covariates) ---------
# One row per analysable client record. Population rule = the frozen gate's: recode
# files, completed observation AND completed exit (c020a == 1 & c020b == 1), completed
# facility (v010a == 1); facilities with no such record are absent by construction
# (gate: drop). Records with A or Y missing are excluded and counted. N_k is then the
# facility's record count IN THIS FRAME -- exactly the N_k the sampler uses.
SPECS <- list(
  # Anchors Y and the shared covariate codings: Turcotte-Tremblay AM, Lee HY, Kruk ME.
  # BMJ Open 2023;13(12):e071037, doi:10.1136/bmjopen-2022-071037. The exposure A is the
  # urban versus rural facility. Applied pooled over ANC + FP + sick child.
  TT2023_pooled = list(
    covariates = c("age_grp", "educ", "service"),
    build = function(cli, s) {
      age <- as.numeric(cli$c511)
      data.frame(
        A = ifelse(cli$c003 == 1, 1L, ifelse(cli$c003 == 2, 0L, NA_integer_)),     # urban (1) vs rural (0)
        Y = ifelse(cli$c520 == 1, 1L, ifelse(cli$c520 %in% c(2, 3), 0L, NA_integer_)),  # very satisfied
        age_grp = factor(ifelse(is.na(age) | age >= 98, "not specified",
                         ifelse(age <= 19, "<=19", ifelse(age <= 35, "20-35", ">=36"))),
                         levels = c("<=19", "20-35", ">=36", "not specified")),
        educ = factor(ifelse(is.na(cli$c512) | cli$c512 == 8, "not specified",
                      c("none", "primary", "secondary", "tertiary+")[as.numeric(cli$c512) + 1]),
                      levels = c("none", "primary", "secondary", "tertiary+", "not specified")),
        service = factor(s, levels = SERVICES))
    },
    vars = c("c003", "c520", "c511", "c512"))
)
covariates_of <- function(spec, design = FALSE) c(spec$covariates, if (design) DESIGN_COVARIATES)

add_design_covariates <- function(fr) {
  N <- as.vector(table(fr$fac)[as.character(fr$fac)])
  fr$logNk <- log(N)
  fr$design_stratum <- factor(fr$stratum)
  fr
}

build_frame <- function(data_root, survey, spec_id, reader = read_dta_lower) {
  spec <- SPECS[[spec_id]]; if (is.null(spec)) fail("unknown spec_id ", spec_id)
  v <- VARS$recode
  fc <- reader(find_dta(data_root, survey, "recode", "FC"), cols = c(v$fac_id, v$fac_done, v$type))
  done <- !is.na(fc[[v$fac_done]]) & fc[[v$fac_done]] == DONE_CODE
  fac <- data.frame(fac = as.vector(fc[[v$fac_id]])[done], stratum = type_group_of(fc[[v$type]][done], survey))
  if (anyNA(fac$stratum)) fail("unmapped facility type among completed facilities")
  rows <- lapply(SERVICES, function(s) {
    d <- reader(find_dta(data_root, survey, "recode", s), cols = c(v$cli_id, v$cli_obs[[s]], v$cli_exit[[s]], spec$vars))
    keep <- !is.na(d[[v$cli_obs[[s]]]]) & d[[v$cli_obs[[s]]]] == DONE_CODE &
            !is.na(d[[v$cli_exit[[s]]]]) & d[[v$cli_exit[[s]]]] == DONE_CODE
    d <- as.data.frame(lapply(d[keep, , drop = FALSE], as.vector))
    cbind(fac = d[[v$cli_id]], spec$build(d, s))
  })
  fr <- merge(do.call(rbind, rows), fac, by = "fac")          # inner join = completed facilities only
  n_all <- nrow(fr); fr <- fr[!is.na(fr$A) & !is.na(fr$Y), ]
  fr <- add_design_covariates(fr)                              # N_k AFTER the A/Y exclusion
  attr(fr, "n_excluded_missing_AY") <- n_all - nrow(fr)
  fr
}

# --- the imposed two-stage design (pure; no I/O) ------------------------------------
# Stage 1: SRS-WOR of m_h = max(min_m, round(f * M_h)) facilities per stratum.
# Stage 2: SRS-WOR of min(c, N_k) records per sampled facility (N_k = frame count).
# Returns the sample with exact pi_fac, pi_rec and w = 1 / (pi_fac * pi_rec).
draw_two_stage <- function(frame, f, c, min_m = 2L) {
  fac <- unique(frame[, c("fac", "stratum")])
  pick <- unlist(lapply(split(fac$fac, fac$stratum), function(ks) {
    m <- max(min_m, round(f * length(ks))); if (m > length(ks)) m <- length(ks)
    ks[sample.int(length(ks), m)]
  }), use.names = FALSE)
  M_h <- table(fac$stratum); m_h <- table(fac$stratum[fac$fac %in% pick])
  s <- frame[frame$fac %in% pick, , drop = FALSE]
  s <- do.call(rbind, lapply(split(s, s$fac), function(d) {
    N <- nrow(d); n <- min(c, N)
    d <- d[sample.int(N, n), , drop = FALSE]
    d$N_k <- N; d$n_k <- n; d
  }))
  s$pi_fac <- as.numeric(m_h[as.character(s$stratum)] / M_h[as.character(s$stratum)])
  s$pi_rec <- s$n_k / s$N_k
  s$w <- 1 / (s$pi_fac * s$pi_rec)
  rownames(s) <- NULL
  s
}

# Per-replicate design diagnostics (aggregates only).
sample_diagnostics <- function(s, c) {
  fk <- s[!duplicated(s$fac), ]
  cv <- sd(s$w) / mean(s$w)
  data.frame(n_psu = nrow(fk), n_rec = nrow(s), w_min = min(s$w), w_median = median(s$w), w_max = max(s$w),
             w_cv = cv, kish_deff_w = 1 + cv^2, n_distinct_pi_rec = length(unique(round(s$pi_rec, 10))),
             share_fac_subsampled = mean(fk$N_k > c), share_A1_weighted = sum(s$w * s$A) / sum(s$w))
}

# --- fitting one sample: arm-level failures are kept, never drop the replicate --------
fit_arms <- function(s, spec, prot, tmle_fun = svytmle, aipw_fun = svyaipw) {
  e <- prot$estimation
  des <- svydesign(ids = ~fac, strata = ~stratum, weights = ~w, data = s, nest = TRUE)
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  failed <- function(method, err) data.frame(method = method, b = NA_real_, se = NA_real_, df = NA_real_,
                                             lcl = NA_real_, ucl = NA_real_, status = paste0("failed: ", class(err)[1]))
  run_t <- function(arms, nw, covs) q(tmle_fun(design = des, outcome = "Y", treatment = "A", covariates = covs,
                                      learners = e$learners, arms = arms, family = "binomial", variance = "linearization",
                                      cf_folds = e$cf_folds, g_bound = e$g_bound %||% 0.05, level = 0.95,
                                      nuisance_weights = nw))$results
  run_a <- function(nw, covs) q(aipw_fun(design = des, outcome = "Y", treatment = "A", covariates = covs,
                                 learners = e$learners, arms = "cf", family = "binomial", cf_folds = e$cf_folds,
                                 g_bound = e$g_bound %||% 0.05, level = 0.95, nuisance_weights = nw))$results
  C  <- covariates_of(spec); Cd <- covariates_of(spec, design = TRUE)
  # PRIMARY: survey-weighted nuisances, substantive C
  t <- tryCatch(cbind(run_t("all", "survey", C)[, RES_COLS], status = "ok"), error = function(err) NULL)
  if (is.null(t))                                  # salvage: each arm on its own
    t <- do.call(rbind, lapply(names(ARM_CODES), function(a)
      tryCatch(cbind(run_t(a, "survey", C)[, RES_COLS], status = "ok (single-arm refit)"), error = function(err) failed(ARM_CODES[[a]], err))))
  a <- tryCatch(cbind(run_a("survey", C)[, RES_COLS], status = "ok"), error = function(err) failed("AIPW-CF", err))
  # SENSITIVITY: package-default unweighted nuisances, C + design stratum + log N_k
  sc <- tryCatch({ r <- run_t("cf", "unweighted", Cd)[, RES_COLS]; r$method <- SENS_CF; cbind(r, status = "ok") },
                 error = function(err) failed(SENS_CF, err))
  sa <- tryCatch({ r <- run_a("unweighted", Cd)[, RES_COLS]; r$method <- SENS_AIPW; cbind(r, status = "ok") },
                 error = function(err) failed(SENS_AIPW, err))
  out <- rbind(t, a, sc, sa)
  missing <- setdiff(ALL_ARMS, out$method)          # an arm the engine silently omitted
  if (length(missing)) out <- rbind(out, do.call(rbind, lapply(missing, function(m) failed(m, structure(list(), class = "omitted")))))
  out[match(ALL_ARMS, out$method), ]
}

# The stabilised census pseudo-reference: each estimator family on the FULL frame
# (w = 1), one split per call; callers take the median over R_ref splits.
reference_split <- function(frame, spec, prot, seed) {
  e <- prot$estimation
  C <- covariates_of(spec); Cd <- covariates_of(spec, design = TRUE)
  frame$w <- 1
  des <- svydesign(ids = ~fac, strata = ~stratum, weights = ~w, data = frame, nest = TRUE)
  q <- function(expr) suppressWarnings(suppressMessages(expr))
  tm <- function(arms, nw, covs) q(svytmle(design = des, outcome = "Y", treatment = "A", covariates = covs,
                 learners = e$learners, arms = arms, family = "binomial", variance = "linearization",
                 cf_folds = e$cf_folds, g_bound = e$g_bound %||% 0.05, nuisance_weights = nw))$results
  ai <- function(nw, covs) q(svyaipw(design = des, outcome = "Y", treatment = "A", covariates = covs,
                 learners = e$learners, arms = "cf", family = "binomial", cf_folds = e$cf_folds,
                 g_bound = e$g_bound %||% 0.05, nuisance_weights = nw))$results
  set.seed(seed)
  t <- tm("all", "survey", C); a <- ai("survey", C)
  tu <- tm("cf", "unweighted", Cd); au <- ai("unweighted", Cd)
  fams <- c(tmle_sf = t$b[t$method == "Fully-Aware"], tmle_sf_unw = t$b[t$method == "Non-Aware"],
            tmle_cv = t$b[t$method == "Fully-Aware-CV"], tmle_cf = t$b[t$method == "Fully-Aware-CF"],
            aipw_cf = a$b[a$method == "AIPW-CF"],
            tmle_cf_unwd = tu$b[tu$method == "Fully-Aware-CF"], aipw_cf_unwd = au$b[au$method == "AIPW-CF"])
  data.frame(family = names(fams), value = unname(fams))
}

# Census overlap diagnostic: out-of-fold SuperLearner propensity P(A = 1 | C) on the
# frame (the Z matrix holds each learner's cross-validated predictions). Approximates
# what the estimators' propensity models see; reported, never used to respecify.
census_overlap <- function(frame, spec, prot) {
  e <- prot$estimation; b <- e$g_bound %||% 0.05
  one <- function(design, label) {
    X <- frame[, covariates_of(spec, design = design)]
    set.seed(e$seed + 200000L)
    sl <- suppressWarnings(suppressMessages(SuperLearner::SuperLearner(
      Y = frame$A, X = X, family = binomial(), SL.library = e$learners, cvControl = list(V = 5L))))
    g <- as.vector(sl$Z %*% sl$coef)
    data.frame(covariate_set = label,
               stat = c("p01", "p05", "p25", "p50", "p75", "p95", "p99", "share_below_bound", "share_above_1_minus_bound",
                        "share_truncated"),
               value = c(quantile(g, c(.01, .05, .25, .5, .75, .95, .99), names = FALSE), mean(g < b), mean(g > 1 - b),
                         mean(g < b | g > 1 - b)))
  }
  rbind(one(FALSE, "primary C"), one(TRUE, "C + design (sensitivity)"))
}

# --- calibration metrics (pure) ------------------------------------------------------
# reps: data.frame(rep, method, b, se, lcl, ucl, status); theta: by estimator family.
# Denominators: R_attempted (all replicates) and R_ok (arm succeeded).
calibration <- function(reps, theta) {
  do.call(rbind, lapply(split(reps, reps$method), function(d) {
    m <- d$method[1]                                 # captured BEFORE failed rows are dropped
    R_att <- length(unique(d$rep))
    th <- unname(theta[ARM_FAMILY[[m]]])            # NA if the family has no reference
    d <- d[is.finite(d$b) & is.finite(d$se), ]
    R <- nrow(d)
    if (R < 2) return(data.frame(method = m, R_attempted = R_att, R_ok = R,
                                 reference = th, bias = NA, bias_mcse = NA, mc_sd = NA, mean_se = NA, var_ratio = NA,
                                 se_ratio_sqrt = NA, mean_se_over_sd = NA, coverage = NA, coverage_mcse = NA, mean_ci_width = NA))
    v_hat <- mean(d$se^2); v_mc <- var(d$b); cover <- mean(d$lcl <= th & th <= d$ucl)
    data.frame(method = m, R_attempted = R_att, R_ok = R, reference = th,
               bias = mean(d$b) - th, bias_mcse = sd(d$b) / sqrt(R),
               mc_sd = sqrt(v_mc), mean_se = mean(d$se),
               var_ratio = v_hat / v_mc, se_ratio_sqrt = sqrt(v_hat / v_mc),
               mean_se_over_sd = mean(d$se) / sqrt(v_mc),
               coverage = cover, coverage_mcse = sqrt(cover * (1 - cover) / R),
               mean_ci_width = mean(d$ucl - d$lcl))
  }))
}

# --- self-test ---------------------------------------------------------------------------
selftest_resample <- function() {
  n_ok <- 0; check <- function(cond, what) { if (!isTRUE(cond)) fail("SELFTEST FAILED: ", what); n_ok <<- n_ok + 1 }
  set.seed(7)
  # Synthetic census: 3 strata, 60 facilities, N_k from 1 to 12.
  fac <- data.frame(fac = 1:60, stratum = rep(c("h1", "h2", "h3"), c(10, 30, 20)), N = sample(1:12, 60, TRUE))
  frame <- do.call(rbind, lapply(seq_len(nrow(fac)), function(i)
    data.frame(fac = fac$fac[i], stratum = fac$stratum[i], id = seq_len(fac$N[i]), A = i %% 2)))
  frame <- add_design_covariates(frame)
  s <- draw_two_stage(frame, f = 0.3, c = 3)
  m_h <- tapply(s$fac, s$stratum, function(x) length(unique(x)))
  check(identical(as.integer(m_h[c("h1", "h2", "h3")]), c(3L, 9L, 6L)), "stage-1 m_h = round(f * M_h)")
  check(all(s$n_k == pmin(3, s$N_k)) && all(tapply(s$id, s$fac, length) == tapply(s$n_k, s$fac, `[`, 1)), "stage-2 take min(c, N_k)")
  check(all(abs(s$w - 1 / ((m_h[s$stratum] / c(h1 = 10, h2 = 30, h3 = 20)[s$stratum]) * (s$n_k / s$N_k))) < 1e-12), "exact weights")
  check(all(abs(exp(s$logNk) - s$N_k) < 1e-9) && identical(as.character(s$design_stratum), s$stratum),
        "design covariates equal the sampler's N_k and stratum (sensitivity arms: ignorability given C + design)")
  est <- replicate(4000, sum(draw_two_stage(frame, f = 0.3, c = 3)$w))
  check(abs(mean(est) - nrow(frame)) < 4 * sd(est) / sqrt(4000), "HT estimate of N is unbiased (MC)")
  dg <- sample_diagnostics(s, 3)
  check(dg$n_psu == 18 && dg$n_rec == nrow(s) && dg$share_fac_subsampled == mean(s$N_k[!duplicated(s$fac)] > 3), "sample diagnostics")
  # Arm-level failure handling with mock fitters: "all" fails, "cv" fails alone, AIPW fails.
  fake <- function(arms) { m <- if (identical(arms, "all")) unname(ARM_CODES) else ARM_CODES[[arms]]
    list(results = data.frame(method = m, b = 0.1, se = 0.05, df = 10, lcl = 0, ucl = 0.2)) }
  seen <- list()
  mock_t <- function(..., covariates, arms, nuisance_weights) {
    seen[[length(seen) + 1]] <<- list(f = "tmle", arms = arms, nw = nuisance_weights, covs = covariates)
    if (identical(arms, "all") || identical(arms, "cv")) stop("boom"); fake(arms) }
  mock_a <- function(..., covariates, nuisance_weights) {
    seen[[length(seen) + 1]] <<- list(f = "aipw", nw = nuisance_weights, covs = covariates)
    if (identical(nuisance_weights, "survey")) stop("boom")
    list(results = data.frame(method = "AIPW-CF", b = 0.2, se = 0.05, df = 10, lcl = 0, ucl = 0.3)) }
  pr <- list(estimation = list(learners = "SL.glm", cf_folds = 2, g_bound = 0.05))
  sp <- list(covariates = c("age_grp"))
  fa <- fit_arms(s, sp, pr, tmle_fun = mock_t, aipw_fun = mock_a)
  check(identical(fa$method, ALL_ARMS) && sum(grepl("^ok", fa$status)) == 6 &&
        is.na(fa$b[fa$method == "Fully-Aware-CV"]) && is.na(fa$b[fa$method == "AIPW-CF"]) &&
        fa$b[fa$method == SENS_AIPW] == 0.2, "arm-level failures kept, others salvaged (8 arms)")
  prim <- Filter(function(z) identical(z$nw, "survey"), seen); sens <- Filter(function(z) identical(z$nw, "unweighted"), seen)
  check(length(prim) >= 2 && all(vapply(prim, function(z) identical(z$covs, "age_grp"), TRUE)) &&
        length(sens) == 2 && all(vapply(sens, function(z) identical(z$covs, c("age_grp", DESIGN_COVARIATES)), TRUE)),
        "primary = survey-weighted nuisances with substantive C; sensitivity = unweighted with C + design")
  # Calibration metrics and denominators.
  reps <- data.frame(rep = rep(1:4, 2), method = rep(c("AIPW-CF", "Fully-Aware-CF"), each = 4),
                     b = c(0.9, 1.1, 1.0, 1.0, 1, NA, 1.2, 0.8), se = 0.2)
  reps$lcl <- reps$b - 1.96 * reps$se; reps$ucl <- reps$b + 1.96 * reps$se
  cal <- calibration(reps, c(aipw_cf = 1, tmle_cf = 1))
  a1 <- cal[cal$method == "AIPW-CF", ]; c1 <- cal[cal$method == "Fully-Aware-CF", ]
  check(abs(a1$bias) < 1e-12 && abs(a1$var_ratio - 0.04 / var(c(0.9, 1.1, 1, 1))) < 1e-12 && a1$coverage == 1, "calibration metrics")
  check(c1$R_attempted == 4 && c1$R_ok == 3, "denominators: R_attempted vs R_ok")
  # an arm that failed in EVERY replicate is reported (R_ok = 0), not a crash
  dead <- rbind(reps, data.frame(rep = 1:4, method = SENS_AIPW, b = NA_real_, se = NA_real_, lcl = NA_real_, ucl = NA_real_))
  cz <- tryCatch(calibration(dead, c(aipw_cf = 1, tmle_cf = 1, aipw_cf_unwd = 0.5)), error = function(e) NULL)
  check(!is.null(cz) && identical(cz$method[cz$R_ok == 0], SENS_AIPW) && cz$R_attempted[cz$method == SENS_AIPW] == 4 &&
        cz$reference[cz$method == SENS_AIPW] == 0.5 && is.na(cz$var_ratio[cz$method == SENS_AIPW]),
        "entirely failed arm reported with R_ok = 0")
  # build_frame on synthetic recode files: population rule, codings, exclusions, N_k.
  root <- file.path(tempfile("rs-"), "raw"); sv <- "MW_SPA_2013-14"; v <- VARS$recode
  wd <- function(mod, df) { d <- file.path(root, sv, "recode", FILES[[sv]]$recode[[mod]]); dir.create(d, recursive = TRUE)
                            haven::write_dta(df, file.path(d, "X.dta")) }
  wd("FC", setNames(data.frame(c(11, 12, 13), c(1, 1, 2), c(2, 5, 8)), c(v$fac_id, v$fac_done, v$type)))
  cl <- data.frame(inv_id = c(11, 11, 12, 12, 13, 12), c020a = c(1, 1, 1, 0, 1, 1), c020b = c(1, 1, 1, 1, 1, 1),
                   c003 = c(1, 1, 2, 2, 2, 2), c520 = c(1, 2, 3, 1, 1, NA), c511 = c(18, 40, 98, 25, 30, 22), c512 = c(0, 3, 8, 1, 2, 1))
  for (s_ in SERVICES) wd(s_, cl)
  fr <- build_frame(root, sv, "TT2023_pooled")
  # per service: rows 1-3 kept (row 4 obs incomplete, row 5 facility 13 not completed, row 6 Y missing)
  check(nrow(fr) == 9 && attr(fr, "n_excluded_missing_AY") == 3, "population rule and A/Y exclusion")
  f1 <- fr[fr$service == "AN", ]; f1 <- f1[order(f1$fac, -f1$Y), ]
  check(identical(f1$A, c(1L, 1L, 0L)) && identical(f1$Y, c(1L, 0L, 0L)), "A = urban, Y = very satisfied")
  check(identical(as.character(f1$age_grp), c("<=19", ">=36", "not specified")) &&
        identical(as.character(f1$educ), c("none", "tertiary+", "not specified")), "age group and education coding")
  check(identical(sort(unique(fr$stratum)), c("health centre", "hospital")), "strata from v007 type groups")
  check(all(abs(exp(fr$logNk[fr$fac == 11]) - 6) < 1e-9) && all(abs(exp(fr$logNk[fr$fac == 12]) - 3) < 1e-9),
        "frame N_k counted after the A/Y exclusion")
  # Job plumbing: f grid in the protocol, ranges, the serial map path, chunk collection.
  okp <- list(status = "FROZEN", frozen_by = "T", frozen_on = "2026-01-01", anchor = "a", spec_id = "TT2023_pooled",
              survey = "MW_SPA_2013-14", design = list(f = list(0.1, 0.25, 0.5), primary_f = 0.25, c = 3),
              estimation = list(R = 20, R_ref = 5, cf_folds = 5, learners = list("SL.glm"), seed = 1,
                                svytmle_version = "0.1.2", svytmle_packaged = "x"))
  vp <- validate_protocol(okp)
  check(identical(vp$design$f, c(0.1, 0.25, 0.5)) && vp$design$primary_f == 0.25, "f grid parsed with primary_f")
  np <- okp; np$design$primary_f <- NULL
  check(grepl("primary_f", tryCatch({ validate_protocol(np); "" }, error = function(e) conditionMessage(e))),
        "multi-f protocol without primary_f rejected")
  check(identical(ftag(0.25), "f250") && identical(ftag(0.1), "f100") && identical(parse_range("3:7", 10), 3:7) &&
        identical(parse_range(NA, 4), 1:4) && grepl("range", tryCatch({ parse_range("5:2", 10); "" }, error = function(e) conditionMessage(e))),
        "ftag and parse_range")
  old <- Sys.getenv("SPA_WORKERS"); Sys.setenv(SPA_WORKERS = "1")
  check(identical(par_map(1:3, function(i, k) i * k, k = 2), list(2, 4, 6)), "par_map serial path")
  if (nzchar(old)) Sys.setenv(SPA_WORKERS = old) else Sys.unsetenv("SPA_WORKERS")
  cdir <- tempfile("chunks-"); dir.create(cdir)
  wc <- function(nm, ids) write.csv(data.frame(rep = ids, b = ids / 10), file.path(cdir, nm), row.names = FALSE)
  wc("a.csv", 1:3); wc("b.csv", 4:6)
  m <- collect_chunks(file.path(cdir, c("b.csv", "a.csv")), "rep", 1:6, "t")
  check(identical(m$rep, 1:6), "collect_chunks merges complete chunks in order")
  check(grepl("missing", tryCatch({ collect_chunks(file.path(cdir, "a.csv"), "rep", 1:6, "t"); "" }, error = function(e) conditionMessage(e))),
        "collect_chunks refuses missing replicates")
  wc("c.csv", 3:4)
  check(grepl("more than one chunk", tryCatch({ collect_chunks(file.path(cdir, c("a.csv", "b.csv", "c.csv")), "rep", 1:6, "t"); "" },
                                              error = function(e) conditionMessage(e))), "collect_chunks refuses duplicated replicates")
  expect <- tryCatch({ validate_protocol(list(status = "DRAFT")); "no error" }, error = function(e) conditionMessage(e))
  check(grepl("status must be FROZEN", expect) && grepl("svytmle_version", expect), "DRAFT or unpinned protocol rejected")
  message(sprintf("selftest (resample): %d checks passed (synthetic only)", n_ok))
}

# --- gated real-data modes -----------------------------------------------------------------
SCRIPT_RELS <- c("src/application/R/02-resample.R", "src/application/R/01-count-audit.R")
RAW_DIR     <- function() Sys.getenv("SPA_RAW_DIR", file.path(REPO, "data", "spa"))   # a private (chmod 700) dir
RESULTS_DIR <- function() Sys.getenv("SPA_RESULTS_DIR", file.path(REPO, "outputs", "application"))
CHUNK_DIR   <- function() file.path(RESULTS_DIR(), "chunks")                        # intermediate (gitignored)

gated_open <- function(survey) {
  ca <- commandArgs(FALSE); ca <- ca[seq_len(match("--args", ca, nomatch = length(ca) + 1L) - 1L)]
  self <- sub("^--file=", "", grep("^--file=", ca, value = TRUE))
  if (length(self) != 1 || !identical(normalizePath(self, winslash = "/", mustWork = FALSE),
                                      normalizePath(file.path(REPO, SCRIPT_RELS[1]), winslash = "/")))
    fail("run directly: Rscript src/application/R/02-resample.R <mode> --survey <S>")
  prov <- preflight_code(REPO, c(PROTOCOL_REL, "config/application/gate.yaml", SCRIPT_RELS))
  prot <- validate_protocol(read_gate_yaml(file.path(REPO, PROTOCOL_REL)))
  env <- check_pins(prot)
  if (!survey %in% c(prot$survey, unlist(prot$replication_surveys))) fail("survey not covered by the frozen protocol")
  list(prot = prot, prov = prov, env = env)
}

ftag <- function(f) sprintf("f%03d", round(1000 * f))                       # 0.25 -> "f250"
out_path <- function(survey, what, ext = "csv", dir = RESULTS_DIR()) file.path(dir, sprintf("resample_%s_%s.%s", survey, what, ext))
write_once <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) fail("refusing to overwrite ", basename(path), " (delete it deliberately and record why)")
  write.csv(df, path, row.names = FALSE)
}
log_run <- function(survey, mode, g, note) {
  runs <- file.path(RESULTS_DIR(), "RESAMPLE-RUNS.md")
  dir.create(dirname(runs), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(runs)) writeLines(c("# Resampling runs (append-only)", "",
    "| run (UTC) | host | survey | mode | HEAD | protocol blob | script blob | svytmle (Packaged) | note |",
    "|---|---|---|---|---|---|---|---|---|"), runs)
  cat(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s (%s) | %s |\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
              Sys.info()[["nodename"]], survey, mode, substr(g$prov[["head"]], 1, 12),
              substr(g$prov[[paste0("blob:", PROTOCOL_REL)]], 1, 12), substr(g$prov[[paste0("blob:", SCRIPT_RELS[1])]], 1, 12),
              g$env[["svytmle"]], g$env[["svytmle_packaged"]], note), file = runs, append = TRUE)
}

# Parallel map. Linux: forked mclapply over the Slurm allocation. Windows: a PSOCK
# pool. Results never depend on the worker count: every replicate / split sets its own seed.
n_workers <- function() {
  n <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", Sys.getenv("SPA_WORKERS", NA))))
  if (is.na(n)) n <- max(1L, parallel::detectCores() - 2L)
  n
}
par_map <- function(X, FUN, ...) {
  n <- min(n_workers(), length(X))
  if (n <= 1) return(lapply(X, FUN, ...))
  if (.Platform$OS.type == "unix") {
    out <- parallel::mclapply(X, FUN, ..., mc.cores = n, mc.preschedule = FALSE)
    bad <- vapply(out, function(o) inherits(o, "try-error") || is.null(o), logical(1))
    if (any(bad)) fail(sum(bad), " worker task(s) died (e.g. out of memory); nothing written for this chunk")
    return(out)
  }
  cl <- makeCluster(n); on.exit(stopCluster(cl))
  clusterCall(cl, function(repo) { setwd(repo); sys.source("src/application/R/02-resample.R", envir = globalenv()); TRUE }, REPO)
  parLapply(cl, X, FUN, ...)
}

parse_range <- function(x, max_n) {
  if (is.na(x) || !nzchar(x)) return(seq_len(max_n))
  ab <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
  if (length(ab) != 2 || anyNA(ab) || ab[1] < 1 || ab[2] < ab[1] || ab[2] > max_n) fail("range must be a:b within 1:", max_n)
  seq.int(ab[1], ab[2])
}

# Reference splits a:b (default all). Chunk files are collected by --collect.
mode_reference <- function(survey, splits_arg) {
  g <- gated_open(survey); prot <- g$prot; spec <- SPECS[[prot$spec_id]]
  sp <- parse_range(splits_arg, prot$estimation$R_ref)
  frame <- build_frame(RAW_DIR(), survey, prot$spec_id)
  out <- do.call(rbind, par_map(sp, function(r, frame, spec, prot)
    cbind(split = r, reference_split(frame, spec, prot, seed = prot$estimation$seed + 100000L + r)), frame, spec, prot))
  write_once(out, file.path(CHUNK_DIR(), sprintf("resample_%s_refsplits_%04d-%04d.csv", survey, min(sp), max(sp))))
  log_run(survey, "reference", g, sprintf("splits %d-%d; %d records; %d facilities; %d excluded (missing A/Y)", min(sp), max(sp),
          nrow(frame), length(unique(frame$fac)), attr(frame, "n_excluded_missing_AY")))
}

# Census overlap + environment (one cheap task per survey).
mode_overlap <- function(survey) {
  g <- gated_open(survey); prot <- g$prot; spec <- SPECS[[prot$spec_id]]
  frame <- build_frame(RAW_DIR(), survey, prot$spec_id)
  ov <- census_overlap(frame, spec, prot)
  write_once(ov, out_path(survey, "census_overlap"))
  write_once(data.frame(item = names(g$env), value = unname(g$env)), out_path(survey, "environment"))
  tr <- ov$value[ov$stat == "share_truncated"]
  log_run(survey, "overlap", g, sprintf("census truncation share %.3f (primary C), %.3f (C + design)", tr[1], tr[2]))
}

# Replicates a:b at one stage-1 fraction f. Replicate r uses seed + r at every f
# (common random numbers across the f grid); sampling AND folds are fresh per replicate.
mode_run <- function(survey, f_arg, reps_arg) {
  g <- gated_open(survey); prot <- g$prot; spec <- SPECS[[prot$spec_id]]
  f <- if (is.na(f_arg)) prot$design$primary_f else as.numeric(f_arg)
  if (!isTRUE(any(abs(prot$design$f - f) < 1e-12))) fail("f = ", f, " is not in the frozen design$f")
  rr <- parse_range(reps_arg, prot$estimation$R)
  frame <- build_frame(RAW_DIR(), survey, prot$spec_id)
  out <- par_map(rr, function(r, frame, spec, prot, f) {
    set.seed(prot$estimation$seed + r)
    s <- draw_two_stage(frame, f, prot$design$c)
    list(arms = cbind(rep = r, f = f, fit_arms(s, spec, prot)), diag = cbind(rep = r, f = f, sample_diagnostics(s, prot$design$c)))
  }, frame, spec, prot, f)
  reps <- do.call(rbind, lapply(out, `[[`, "arms")); diag <- do.call(rbind, lapply(out, `[[`, "diag"))
  stem <- sprintf("resample_%s_%s_reps_%05d-%05d", survey, ftag(f), min(rr), max(rr))
  write_once(reps, file.path(CHUNK_DIR(), paste0(stem, "_replicates.csv")))
  write_once(diag, file.path(CHUNK_DIR(), paste0(stem, "_diagnostics.csv")))
  nf <- tapply(!grepl("^ok", reps$status), reps$method, sum)
  log_run(survey, "run", g, sprintf("f=%g reps %d-%d; c=%d; arm failures: %s", f, min(rr), max(rr), prot$design$c,
                                    paste(sprintf("%s=%d", names(nf), nf), collapse = ", ")))
}

# Merge chunks after checking completeness: every replicate 1..R at every f, and every
# reference split 1..R_ref, present EXACTLY once. Refuses to write anything otherwise.
collect_chunks <- function(files, id_col, expected, what) {
  if (!length(files)) fail("no chunk files for ", what)
  d <- do.call(rbind, lapply(files, read.csv))
  per_file <- lapply(files, function(fl) unique(read.csv(fl)[[id_col]]))
  dup <- unlist(per_file)[duplicated(unlist(per_file))]
  miss <- setdiff(expected, unlist(per_file))
  if (length(dup)) fail(what, ": ", length(unique(dup)), " id(s) appear in more than one chunk (e.g. ", paste(head(unique(dup), 5), collapse = ", "), ")")
  if (length(miss)) fail(what, ": ", length(miss), " id(s) missing (e.g. ", paste(head(miss, 5), collapse = ", "), ")")
  d[order(d[[id_col]]), , drop = FALSE]
}
mode_collect <- function(survey) {
  g <- gated_open(survey); prot <- g$prot
  cd <- CHUNK_DIR()
  ref_files <- list.files(cd, sprintf("^resample_%s_refsplits_.*[.]csv$", survey), full.names = TRUE)
  splits <- collect_chunks(ref_files, "split", seq_len(prot$estimation$R_ref), "reference splits")
  agg <- aggregate(value ~ family, splits, median)
  write_once(splits, out_path(survey, "reference_splits"))
  write_once(data.frame(family = agg$family, pseudo_reference = agg$value), out_path(survey, "reference"))
  for (f in prot$design$f) {
    pat <- sprintf("^resample_%s_%s_reps_.*_replicates[.]csv$", survey, ftag(f))
    rf <- list.files(cd, pat, full.names = TRUE)
    reps <- collect_chunks(rf, "rep", seq_len(prot$estimation$R), paste0("replicates at f = ", f))
    diag <- collect_chunks(sub("_replicates[.]csv$", "_diagnostics.csv", rf), "rep", seq_len(prot$estimation$R),
                           paste0("diagnostics at f = ", f))
    write_once(reps, out_path(survey, paste0(ftag(f), "_replicates")))
    write_once(diag, out_path(survey, paste0(ftag(f), "_replicate_diagnostics")))
  }
  log_run(survey, "collect", g, sprintf("R_ref=%d splits; R=%d replicates at f = %s: complete", prot$estimation$R_ref,
                                        prot$estimation$R, paste(prot$design$f, collapse = ", ")))
}

mode_summarise <- function(survey) {
  g <- gated_open(survey); prot <- g$prot
  ref <- read.csv(out_path(survey, "reference"))
  all_cal <- list()
  for (f in prot$design$f) {
    reps <- read.csv(out_path(survey, paste0(ftag(f), "_replicates")))
    cal <- calibration(reps, setNames(ref$pseudo_reference, ref$family))
    set.seed(1)                                          # replicate-bootstrap MCSE for the variance ratio
    cal$var_ratio_mcse <- vapply(cal$method, function(m) {
      d <- reps[reps$method == m & is.finite(reps$b) & is.finite(reps$se), ]
      if (nrow(d) < 3) return(NA_real_)
      sd(replicate(500, { i <- sample.int(nrow(d), replace = TRUE); mean(d$se[i]^2) / var(d$b[i]) }))
    }, numeric(1))
    cal <- cbind(f = f, primary_f = isTRUE(abs(f - prot$design$primary_f) < 1e-12), cal[match(ALL_ARMS, cal$method), ])
    write_once(cal, out_path(survey, paste0(ftag(f), "_calibration")))
    diag <- read.csv(out_path(survey, paste0(ftag(f), "_replicate_diagnostics")))
    ds <- do.call(rbind, lapply(setdiff(names(diag), c("rep", "f")), function(v)
      data.frame(f = f, stat = v, min = min(diag[[v]]), median = median(diag[[v]]), max = max(diag[[v]]))))
    write_once(ds, out_path(survey, paste0(ftag(f), "_diagnostics_summary")))
    all_cal[[length(all_cal) + 1]] <- cal
  }
  allc <- do.call(rbind, all_cal)
  write_once(allc, out_path(survey, "calibration_all_f"))
  log_run(survey, "summarise", g, sprintf("f = %s", paste(prot$design$f, collapse = ", ")))
  print(allc[, c("f", "method", "R_ok", "var_ratio", "var_ratio_mcse", "mean_se_over_sd", "bias", "coverage", "coverage_mcse")], row.names = FALSE)
}

# --- entry point ---------------------------------------------------------------------------
if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  MODES <- c("--selftest", "--reference", "--overlap", "--run", "--collect", "--summarise")
  modes <- intersect(args, MODES)
  if (length(modes) != 1) fail("exactly one of ", paste(MODES, collapse = " | "), " is required")
  opt <- function(k) if (k %in% args) args[which(args == k) + 1] else NA_character_
  survey <- opt("--survey")
  if (modes != "--selftest" && (is.na(survey) || !survey %in% SURVEYS)) fail("--survey must be one of ", paste(SURVEYS, collapse = ", "))
  switch(modes, "--selftest" = selftest_resample(),
         "--reference" = mode_reference(survey, opt("--splits")), "--overlap" = mode_overlap(survey),
         "--run" = mode_run(survey, opt("--f"), opt("--reps")), "--collect" = mode_collect(survey),
         "--summarise" = mode_summarise(survey))
}
