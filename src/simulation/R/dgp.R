# Population generator for the DGP (a_cpsu = 1.25).
#
# Every draw comes from the entity-addressed streams of R/rng_streams.R, the random-number contract,
# under a SYNTHETIC root (tests) or the
# production root (dry and production runs). co comes from dgp_coefficients(a_cpsu)
# with a_cpsu passed explicitly.
#
# Requires, sourced beforehand: R/quadrature.R (gauss_hermite), R/nuisance_truth.R
# (lp_g, lp_Q, Q0_marginal, g0_marginal), R/remainders.R (aipw_phi), R/rng_streams.R.
#
# LAW:
#   K_B(h) = 1 for odd h, 2 for even h; no stratum-level channel, u^C = b_h = 0.
#   per PSU:  C_psu ~ U(-1,1); M_k = round(34 + 14 C_psu), support {20..48};
#             u^A_k ~ N(0, sigma_A[K_B]^2), u^Y_k ~ N(0, sigma_Y[K_B]^2).
#   per unit: C1..C4 ~ U(-1,1) i.i.d.;
#             A ~ Bern(expit(lp_g + u^A)), Y ~ Bern(expit(lp_Q(A) + u^Y)), with lp_Q
#             carrying co$a_cpsu * A * C_psu (S11; a = 1.25 frozen).
#   eta0 = (Q1, Q0, g) is the MARGINAL truth, Gauss-Hermite over u^A, u^Y at 64 nodes
#   phi0 = aipw_phi(Y, A, Q1, Q0, g), the uncentred EIF.
#   M_k depends on C_psu only, never on u^A, u^Y.
#
# TECHNICAL CONVENTIONS:
#   TC-1 streams: one population-level stream per PSU, key (root, master_id, pop_id,
#        psu(h,k)) (draft). A PSU's values never depend on H, J or any other PSU, so
#        the (H,J) population equals the prefix view of any larger (H',J') build with the
#        same ids (prefix rule; population stream cell-free).
#   TC-2 within-PSU draw order (draft): U_c, U_A, U_Y, then 6*M_k uniforms read
#        column-major into an M_k x 6 matrix [C1, C2, C3, C4, V^A, V^Y].
#   TC-3 transforms: C = -1 + 2U; u = sigma * qnorm(U) (inversion); A = 1{V^A < p_A},
#        Y = 1{V^Y < p_Y} (uniform thresholds, draft), so the coefficient-free identity
#        tier (h,k,i,M,K_B,C,U_A,U_Y,V^A,V^Y) cannot depend on any coefficient, and a
#        coefficient change can only flip A or Y, never move a draw.
#   TC-4 size rounding (draft): x = 34 + 14*C_psu in binary64 (product first), then
#        round HALF TO EVEN (size_from_driver). Ties have probability zero under the
#        continuous law, so the frozen discrete law (E M = 34.000) is unaffected;
#        the number of exact ties is recorded in meta$n_size_ties (draft: expected 0).

DEV_ROOT_DEFAULT <- "SYNTHETIC-DEV-E2E"
DEV_LABEL <- "population (spec v9.8 DGP, adopted RNG contract); root in meta$root"
POP_IDENTITY_COLS <- c("h", "k", "i", "M", "KB", "Cpsu", "UA", "UY", "C1", "C2", "C3", "C4", "VA", "VY")   # contract identity tier
POP_DERIVED_COLS <- c("h", "k", "i", "A", "Y", "Q1", "Q0", "g", "phi0")

.assert_dev_root <- function(root) assert_engine_root(root, "population builder")

.is_whole <- function(x, lo = 1) is.numeric(x) && length(x) == 1L && is.finite(x) && x == floor(x) && x >= lo

# K_B by stratum parity.
kb_of_stratum <- function(h) ifelse(as.integer(h) %% 2L == 1L, 1L, 2L)

# One PSU's coefficient-free uniforms (TC-2). Must run inside rng_scoped() with the
# PSU's state set. The draw count 3 + 6*M_k depends on this PSU's own first draw only.
.dev_psu_draws <- function() {
  u0 <- stats::runif(3L)
  M <- size_from_driver(-1 + 2 * u0[1], 34, 14)     # TC-4,
  list(u0 = u0, U = matrix(stats::runif(6L * M), nrow = M))
}

# keep_identity = TRUE also keeps the coefficient-free uniforms (psu$UA, psu$UY; units$VA, units$VY)
# for the identity-tier store hash (contract); the default leaves the frame as before.
build_population <- function(H, J, co, master_id, pop_id, root = DEV_ROOT_DEFAULT, keep_identity = FALSE) {
  .assert_dev_root(root)
  if (!.is_whole(H) || !.is_whole(J)) stop("H and J must be whole numbers >= 1", call. = FALSE)
  if (missing(co) || !is.list(co) || is.null(co$a_cpsu))
    stop("co must come from dgp_coefficients(a_cpsu) with a_cpsu passed explicitly (S11)", call. = FALSE)
  H <- as.integer(H); J <- as.integer(J)
  ctx <- list(root = root, master_id = master_id, population_id = pop_id)

  # TC-1: entity-addressed states, one per PSU, in canonical (h, k) order.
  hk <- data.frame(h = rep(seq_len(H), each = J), k = rep(seq_len(J), times = H))
  states <- rng_states(rng_entity_keys("population", ctx, "psu", hk))
  dr <- rng_scoped(lapply(seq_len(nrow(hk)), function(r) {
    rng_set_state(states[r, ])
    .dev_psu_draws()
  }))

  u0 <- matrix(unlist(lapply(dr, `[[`, "u0"), use.names = FALSE), ncol = 3L, byrow = TRUE)
  KB <- kb_of_stratum(hk$h)
  Cpsu <- -1 + 2 * u0[, 1]
  M <- vapply(dr, function(d) nrow(d$U), integer(1))
  psu <- data.frame(h = as.integer(hk$h), k = as.integer(hk$k), KB = KB, Cpsu = Cpsu, M = M,
                    uA = co$sigma_A[KB] * stats::qnorm(u0[, 2]),
                    uY = co$sigma_Y[KB] * stats::qnorm(u0[, 3]))

  U <- do.call(rbind, lapply(dr, `[[`, "U"))
  idx <- rep(seq_len(nrow(psu)), psu$M)
  units <- data.frame(h = psu$h[idx], k = psu$k[idx], i = sequence(psu$M), KB = psu$KB[idx],
                      C1 = -1 + 2 * U[, 1], C2 = -1 + 2 * U[, 2], C3 = -1 + 2 * U[, 3],
                      C4 = -1 + 2 * U[, 4], Cpsu = psu$Cpsu[idx])

  # Conditional laws include the PSU effects; thresholds are the stored uniforms (TC-3).
  pA <- stats::plogis(lp_g(units, units$KB, co) + psu$uA[idx])
  units$A <- as.integer(U[, 5] < pA)
  pY <- stats::plogis(lp_Q(units$A, units, units$KB, co) + psu$uY[idx])
  units$Y <- as.integer(U[, 6] < pY)

  # Marginal eta0 at GH64; Q0 here is Q0(0, c), as in truth_at().
  gh <- gauss_hermite(64L)
  units$Q1 <- Q0_marginal(1, units, units$KB, co, gh)
  units$Q0 <- Q0_marginal(0, units, units$KB, co, gh)
  units$g <- g0_marginal(units, units$KB, co, gh)
  units$phi0 <- aipw_phi(units$Y, units$A, units$Q1, units$Q0, units$g)

  if (keep_identity) {
    psu$UA <- u0[, 2]; psu$UY <- u0[, 3]
    units$VA <- U[, 5]; units$VY <- U[, 6]
  }
  x <- 34 + 14 * psu$Cpsu
  list(psu = psu, units = units,
       meta = list(H = H, J = J, a_cpsu = co$a_cpsu, master_id = master_id, pop_id = pop_id,
                   root = root, rng = "adopted RNG contract (docs/RNG-CONVENTIONS-DRAFT.md, R2)",
                   rng_version = RNG_DRAFT_VERSION, gh_nodes = 64L, N = nrow(units),
                   size_rounding = "half-to-even on binary64 (contract §8, adopted)",
                   n_size_ties = sum(x - floor(x) == 0.5),
                   label = DEV_LABEL))
}

# Prefix view: the first H strata and the first J PSUs of each. By TC-1 it
# equals build_population(H, J, ...) with the same ids; tested, not assumed.
population_prefix <- function(pop, H, J) {
  if (!.is_whole(H) || !.is_whole(J) || H > pop$meta$H || J > pop$meta$J)
    stop("prefix (H, J) must lie within the population's (H, J)", call. = FALSE)
  keep_p <- pop$psu$h <= H & pop$psu$k <= J
  keep_u <- pop$units$h <= H & pop$units$k <= J
  psu <- pop$psu[keep_p, , drop = FALSE]; rownames(psu) <- NULL
  units <- pop$units[keep_u, , drop = FALSE]; rownames(units) <- NULL
  meta <- pop$meta
  meta$H <- as.integer(H); meta$J <- as.integer(J); meta$N <- nrow(units)
  x <- 34 + 14 * psu$Cpsu
  meta$n_size_ties <- sum(x - floor(x) == 0.5)
  list(psu = psu, units = units, meta = meta)
}

# Identity- and derived-tier store hashes of a population (contract item 4), streamed
# through a temporary file so a 5.4-million-unit frame never holds its bytes twice in memory. The
# bytes are exactly rng_store_bytes()'s (tested). The identity tier needs keep_identity = TRUE.
population_store_hashes <- function(pop, tag, tmpdir = tempdir()) {
  psu <- pop$psu; u <- pop$units
  if (is.null(psu$UA) || is.null(u$VA)) stop("population_store_hashes needs build_population(keep_identity = TRUE)", call. = FALSE)
  idx <- rep(seq_len(nrow(psu)), psu$M)
  get_col <- function(cl) switch(cl, M = psu$M[idx], KB = as.integer(psu$KB[idx]), Cpsu = psu$Cpsu[idx],
                                 UA = psu$UA[idx], UY = psu$UY[idx], A = as.integer(u$A), Y = as.integer(u$Y),
                                 h = as.integer(u$h), k = as.integer(u$k), i = as.integer(u$i), u[[cl]])
  one <- function(cols, tier) {
    f <- tempfile(paste0("store-", tier, "-"), tmpdir = tmpdir); on.exit(unlink(f))
    con <- file(f, "wb")
    types <- ifelse(cols %in% c("h", "k", "i", "M", "KB", "A", "Y"), "i4", "f8")
    writeBin(charToRaw(paste0(tag, ";tier=", tier, ";cols=", paste(cols, collapse = ","), ";types=",
                              paste(types, collapse = ","), ";nrow=", nrow(u), "\n")), con)
    for (j in seq_along(cols)) {
      x <- get_col(cols[j])
      if (anyNA(x)) { close(con); stop("store column '", cols[j], "' contains NA/NaN", call. = FALSE) }
      if (types[j] == "i4") writeBin(as.integer(x), con, size = 4L, endian = "little") else writeBin(as.double(x) + 0, con, size = 8L, endian = "little")
    }
    close(con)
    digest::digest(file = f, algo = "sha256")
  }
  if (any(order(u$h, u$k, u$i, method = "radix") != seq_len(nrow(u)))) stop("units must be in canonical (h, k, i) order", call. = FALSE)
  c(identity = one(POP_IDENTITY_COLS, "identity"), derived = one(POP_DERIVED_COLS, "derived"))
}
