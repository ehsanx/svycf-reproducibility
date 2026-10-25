# Pre-run gate machinery: G3 cell enumeration, G1's resize rule, G5's rate
# regression and R-b's nested rate gates, the G5 nuisance error r_m, and G6 limb (a)'s
# boundary distance and MCSE at the planned count.
#
# Every cell coordinate and count is fixed. G3 is pure
# arithmetic over the fixed cell list; G5 and R-b's gates ran on the full ladders of the dry run
# (its summary: results/simulation/execution-check); G1 and G6 run on the production populations' census before any fitted
# cell (runs/run_prep.R).

# FPC is 200 x 800: the vanishing cell's bridge-only terminal rung needs J = 800.
GATE_MASTERS <- list(MS = c(H = 800L, J = 20L), FPC = c(H = 200L, J = 800L), FS = c(H = 5L, J = 6400L))

# Every frozen cell: one row per (cell, rung, V). R-d rows carry the
# rung's q_m(alpha) for every alpha. E4's deterministic array has no master.
g3_cells <- function() {
  row <- function(cell, regime, H, m_h, J, master, V = 2L, arms = "", pending = "") {
    data.frame(cell = cell, regime = regime, H = as.integer(H), m_h = as.character(m_h), J = as.integer(J),
               master = master, V = as.integer(V), arms = arms, pending = pending, stringsAsFactors = FALSE)
  }
  r <- list()
  for (H in c(25L, 50L, 100L, 200L, 400L, 800L))
    r[[length(r) + 1L]] <- row("E1 many-strata", "many", H, 2L, 20L, "MS",
                               arms = if (H == 800L) "oracle, R-a AIPW (cheap rung)" else "oracle, R-a AIPW/TMLE, R-b AIPW/TMLE")
  for (mh in c(10L, 20L, 40L, 80L, 160L, 320L, 640L))
    r[[length(r) + 1L]] <- row("E1 fixed-strata", "fixed", 5L, mh, 10L * mh, "FS",
                               arms = paste0(if (mh <= 320L) "oracle, R-a AIPW" else "oracle (internal)",
                                             if (mh <= 160L) ", R-a TMLE" else "",
                                             if (5L * mh >= 200L) ", R-d x4" else ""))
  for (J in c(100L, 20L, 8L, 4L))
    r[[length(r) + 1L]] <- row("E2 FPC sweep", "many", 50L, 2L, J, "FPC",
                               arms = if (J == 4L) "oracle, R-a AIPW, R-b AIPW" else "oracle")
  for (HJ in list(c(25L, 10L), c(50L, 20L), c(100L, 40L), c(200L, 80L)))
    r[[length(r) + 1L]] <- row("E2 vanishing-fraction", "many", HJ[1], 2L, HJ[2], "FPC", arms = "census limb")
  r[[length(r) + 1L]] <- row("E2 vanishing-fraction (bridge-only terminal rung)", "many", 200L, 2L, 800L, "FPC",
                             arms = "census limb; the graded bridge (Amendment A2, R22)")
  for (H in c(25L, 50L, 100L, 200L, 400L))
    r[[length(r) + 1L]] <- row("E2 genuinely-strict", "many", H, 2L, 4L, "MS",
                               arms = if (H == 100L) "census limb; oracle, R-a AIPW, R-b AIPW" else "census limb")
  r[[length(r) + 1L]] <- row("E2 WR code-check", "many", 50L, 2L, 4L, "FPC", arms = "oracle (WR draw)")
  r[[length(r) + 1L]] <- row("E2 unequal-allocation pair", "many", 50L, "6/2 by class", 20L, "FPC", arms = "oracle")
  r[[length(r) + 1L]] <- row("E3 R-c foil + R-a controls", "fixed", 5L, 80L, 800L, "FS",
                             arms = "R-c PSU-CF / no-CF / individual-CF (+ coarsened C_psu, R15); R-a PSU-CF / no-CF")
  for (V in c(2L, 5L))
    r[[length(r) + 1L]] <- row("E3 V-sensitivity", "fixed", 5L, 40L, 400L, "FS", V = V, arms = "oracle, R-b AIPW")
  r[[length(r) + 1L]] <- row("E3 power-control panel", "fixed", 5L, 80L, 800L, "FS",
                             arms = "host R-a PSU-CF AIPW: C1 shift, wrong clustering, aligned-error T2 (R10); 160 x 5 (R7)")
  r[[length(r) + 1L]] <- row("G8 U10 cell", "fixed", 5L, 80L, 800L, "FS", arms = "oracle, stage-2 weights omitted (R10, R23); dry run 12 x 15")
  for (H in c(25L, 50L, 100L, 200L, 400L))
    r[[length(r) + 1L]] <- row("E4 deterministic array", "e4", H, 2L, 20L, "none", arms = "oracle (one unit per PSU, R11)")
  do.call(rbind, r)
}

# G3 on one enumeration: (a) R-d headroom for every fixed-strata rung on R-d's grid,
# (b) master sufficiency, (c) integrality. Returns one row per assertion; (a-ii) is reported
# (required = FALSE). Uses rd_q() and RD_ALPHAS from R/learner_rd.R.
g3_check <- function(cells = g3_cells()) {
  out <- list()
  add <- function(cell, check, value, reference, pass, required = TRUE)
    out[[length(out) + 1L]] <<- data.frame(cell = cell, check = check, value = value, reference = reference,
                                           pass = pass, required = required, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cells))) {
    c <- cells[i, ]
    lab <- sprintf("%s (H = %d, m_h = %s, J = %d, V = %d)", c$cell, c$H, c$m_h, c$J, c$V)
    mh <- as.integer(strsplit(c$m_h, "[^0-9]+")[[1]]); mh <- mh[!is.na(mh)]
    if (c$master != "none") {
      ms <- GATE_MASTERS[[c$master]]
      add(lab, sprintf("G3(b) master %s sufficient (J <= %d, H <= %d)", c$master, ms[["J"]], ms[["H"]]),
          c$J, ms[["J"]], c$J <= ms[["J"]] && c$H <= ms[["H"]])
    }
    add(lab, "G3(c) integrality: V divides every m_h and m_vh = m_h/V is whole", min(mh %% c$V), 0, all(mh %% c$V == 0L))
    add(lab, "G3(c) m_h <= J (a stratum cannot sample more PSUs than it has)", max(mh), c$J, all(mh <= c$J))
    if (c$regime == "fixed" && grepl("R-d", c$arms)) {
      m <- 5L * mh[1]; pool <- m * (c$V - 1L) / c$V
      for (a in RD_ALPHAS) {
        q <- rd_q(m, a)
        add(lab, sprintf("G3(a-i) %s at m = %d: q_m = %d <= pool %g", rd_arm(a), m, q, pool), q, pool, q <= pool)
        add(lab, sprintf("G3(a-ii) %s at m = %d: share %.3f <= 0.5 (design margin, reported)", rd_arm(a), m, q / pool),
            q / pool, 0.5, q / pool <= 0.5, required = FALSE)
      }
    }
  }
  do.call(rbind, out)
}

# G1's deterministic resize rule: the cell's planned population count B is doubled, to
# a cap of 4x, until the target contrast exceeds 3 MCSE; otherwise UNGRADABLE. mcse_at(B) gives
# the MCSE at B populations under the cell's MCSE convention. The census contrast is
# known before any fitted run, so the rule cannot use estimator output.
g1_resize <- function(contrast, mcse_at, B_planned, cap = 4L) {
  if (!is.finite(contrast) || B_planned < 1) stop("contrast must be finite and B_planned >= 1", call. = FALSE)
  for (k in c(1L, 2L, 4L)[c(1L, 2L, 4L) <= cap]) {
    B <- B_planned * k
    if (abs(contrast) > 3 * mcse_at(B)) return(list(B = B, multiplier = k, status = "GRADEABLE", mcse = mcse_at(B)))
  }
  list(B = B_planned * cap, multiplier = cap, status = "UNGRADABLE", mcse = mcse_at(B_planned * cap))
}

# CR2 cluster-robust covariance (Bell-McCaffrey) for OLS y ~ X, clusters g.
.cr2_vcov <- function(X, e, cluster) {
  XtX_inv <- solve(crossprod(X))
  meat <- matrix(0, ncol(X), ncol(X))
  for (cl in unique(cluster)) {
    i <- which(cluster == cl)
    Xg <- X[i, , drop = FALSE]
    Hgg <- Xg %*% XtX_inv %*% t(Xg)
    ev <- eigen(diag(length(i)) - Hgg, symmetric = TRUE)
    if (any(ev$values <= 1e-12)) stop("CR2: I - H_gg is singular for cluster ", cl, call. = FALSE)
    Ag <- ev$vectors %*% diag(1 / sqrt(ev$values), length(i)) %*% t(ev$vectors)
    u <- t(Xg) %*% Ag %*% e[i]
    meat <- meat + u %*% t(u)
  }
  XtX_inv %*% meat %*% XtX_inv
}

# G5: OLS of log r_hat on log m across the rungs, CR2 SE clustered on population,
# critical value from t_{G-1}; PASS iff the two-sided 95% Wald interval for the slope contains
# -alpha. Inputs are one row per (population, rung): r_hat already averaged over that
# population's draws at that rung.
g5_fit <- function(r_hat, m, population, alpha, level = 0.95) {
  if (any(!is.finite(r_hat)) || any(r_hat <= 0)) stop("r_hat must be finite and > 0", call. = FALSE)
  y <- log(r_hat); X <- cbind(1, log(m))
  fit <- stats::lm.fit(X, y)
  b <- unname(fit$coefficients)
  V <- .cr2_vcov(X, fit$residuals, population)
  G <- length(unique(population))
  se <- sqrt(V[2, 2]); crit <- stats::qt(1 - (1 - level) / 2, df = G - 1L)
  ci <- b[2] + c(-1, 1) * crit * se
  list(slope = b[2], intercept = b[1], se_cr2 = se, df = G - 1L, crit = crit, ci = ci,
       alpha = alpha, pass = ci[1] <= -alpha && -alpha <= ci[2])
}

# G5's nuisance error for one fold fit:
# ||Q_hat - Q0|| = max over a in {0,1} of the census L2 error of the post-clip fit, and
# ||g_hat - g0|| likewise for g; r = max of the two. units: the census units with the marginal
# truth columns Q1, Q0, g (R/dgp.R).
nuisance_l2 <- function(fit, units) {
  q1 <- sqrt(mean((fit$Q(1, units) - units$Q1)^2))
  q0 <- sqrt(mean((fit$Q(0, units) - units$Q0)^2))
  g <- sqrt(mean((fit$g(units) - units$g)^2))
  c(Q1 = q1, Q0 = q0, Q = max(q1, q0), g = g, r = max(q1, q0, g))
}

# G6 limb (a), boundary reading: the analytic effect of a stochastic power control is measured
# from the NEAREST BOUNDARY of the compliant region [lo, hi] against which its row is graded, not
# from zero. distance > 0 outside the region, <= 0 inside. A control predicted inside its region fails G6
# and is a diagnostic, never graded as a power control; outside, it is graded only if distance > 3 MCSE
# at the planned count (mcse; NA reports the distance only).
g6_boundary <- function(pred, lo, hi, mcse = NA_real_) {
  dist <- pmax(lo - pred, pred - hi)
  status <- ifelse(!is.finite(pred), "no prediction",
            ifelse(dist <= 0, "inside the compliant region: fails G6 limb (a); diagnostic, never graded (R5)",
            ifelse(is.na(mcse), "outside the region: distance reported (no MCSE given)",
            ifelse(dist > 3 * mcse, "outside by more than 3 MCSE: graded power control",
                   "outside by at most 3 MCSE: fails G6 limb (a); implementation check"))))
  data.frame(predicted = pred, lo = lo, hi = hi, distance = dist, mcse = mcse,
             inside = is.finite(dist) & dist <= 0,
             graded = is.finite(dist) & dist > 0 & is.finite(mcse) & dist > 3 * mcse,
             status = status, stringsAsFactors = FALSE)
}

# MCSE at the planned count (operational), for a mean-type statistic of a joint cell under the
# population-block bootstrap: [Var_b(p_b)/B + v/(B R)]^(1/2), p_b the census-predicted per-population
# value over the planned population set, v a conservative per-draw variance.
g6_mcse_joint <- function(pred_b, v_draw, R) {
  B <- length(pred_b)
  sqrt(stats::var(pred_b) / B + v_draw / (B * R))
}
# ... and for a design-conditional cell, per population: (v / R)^(1/2).
g6_mcse_design <- function(v_draw, R) sqrt(v_draw / R)

# R-b's two nested verdicts on its own ladder (one regime's rungs): OLS of log r_hat on log m,
# CR2 SE clustered on population, t_{G-1}; the ONE-SIDED 95% upper bound on the slope. L2 gate: bound < 0;
# quarter-rate gate: bound < -1/4.
rb_rate_gates <- function(r_hat, m, population, level = 0.95) {
  if (any(!is.finite(r_hat)) || any(r_hat <= 0)) stop("r_hat must be finite and > 0", call. = FALSE)
  y <- log(r_hat); X <- cbind(1, log(m))
  fit <- stats::lm.fit(X, y)
  b <- unname(fit$coefficients)
  V <- .cr2_vcov(X, fit$residuals, population)
  G <- length(unique(population))
  se <- sqrt(V[2, 2]); upper <- b[2] + stats::qt(level, df = G - 1L) * se
  list(slope = b[2], se_cr2 = se, df = G - 1L, upper = upper,
       L2 = if (upper < 0) "r_m -> 0 established (L2 gate)" else "L2 gate fails or is unresolved",
       quarter = if (upper < -0.25) "quarter-rate established" else "quarter-rate fails or is unresolved",
       pass_L2 = upper < 0, pass_quarter = upper < -0.25)
}
