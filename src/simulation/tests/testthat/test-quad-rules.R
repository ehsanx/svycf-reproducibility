# Saved quadrature-rule artifacts (R/quad_rules_io.R, grids/): exact round trip,
# checksum and integrity verification, the required-store contract, and the DEV
# artifact's cross-platform identity -- the grids built from it (H = 5, 25, 50, 100, 400) must
# hash to the values recorded in grids/MANIFEST.json on every platform.

for (f in c("quadrature.R", "quad_rules_io.R", "nuisance_truth.R", "p0_grid.R"))
  source(test_path("..", "..", "R", f), local = TRUE)
ENGINE <- normalizePath(test_path("..", ".."), winslash = "/")
withr::defer(quad_rules_deactivate())

SPEC <- data.frame(rule = c("gauss-legendre", "clenshaw-curtis", "gauss-hermite"), n = c(5L, 7L, 12L),
                   stringsAsFactors = FALSE)

test_that("save -> load is an exact round trip, verified by sha256", {
  quad_rules_deactivate()
  r <- quad_rules_build(SPEC)
  p <- withr::local_tempfile(fileext = ".tsv")
  sha <- quad_rules_save(r, p)
  expect_match(sha, "^[0-9a-f]{64}$")
  expect_identical(sha, digest::digest(file = p, algo = "sha256"))
  ld <- quad_rules_load(p, sha)
  expect_identical(names(ld$rules), names(r))
  for (k in names(r)) expect_identical(ld$rules[[k]], r[[k]])
  expect_false(any(grepl("\r", readLines(p))))                  # LF only
})

test_that("a wrong checksum, a tampered value or a bad header is a hard error", {
  quad_rules_deactivate()
  p <- withr::local_tempfile(fileext = ".tsv")
  sha <- quad_rules_save(quad_rules_build(SPEC), p)
  expect_error(quad_rules_load(p, strrep("0", 64)), "has sha256")
  expect_error(quad_rules_load(p), "expected sha256")
  l <- readLines(p)
  l[4] <- sub("\t[^\t]+$", "\t0.5", l[4])                         # decimal no longer matches hex
  q <- withr::local_tempfile(fileext = ".tsv"); con <- file(q, "wb"); writeLines(l, con, sep = "\n"); close(con)
  expect_error(quad_rules_load(q, digest::digest(file = q, algo = "sha256")), "disagree")
  l2 <- readLines(p); l2[1] <- "# something else"
  q2 <- withr::local_tempfile(fileext = ".tsv"); con <- file(q2, "wb"); writeLines(l2, con, sep = "\n"); close(con)
  expect_error(quad_rules_load(q2, digest::digest(file = q2, algo = "sha256")), "bad header")
})

test_that("an active required store serves stored base rules and refuses missing ones", {
  p <- withr::local_tempfile(fileext = ".tsv")
  quad_rules_deactivate()
  r <- quad_rules_build(SPEC)
  ld <- quad_rules_load(p, quad_rules_save(r, p))
  quad_rules_activate(ld, required = TRUE)
  expect_identical(quad_rules_active(), ld$sha256)
  b <- ld$rules[["gauss-legendre:5"]]
  expect_identical(gauss_legendre(5, -0.3, 0.7), list(x = -0.3 + 0.5 * (b$x + 1), w = 0.5 * b$w))
  expect_identical(gauss_hermite(12), ld$rules[["gauss-hermite:12"]])
  expect_error(gauss_legendre(4), "not in the active saved-rule artifact")
  expect_error(clenshaw_curtis(3, 0, 1), "not in the active saved-rule artifact")
  quad_rules_deactivate()
  expect_true(is.na(quad_rules_active()))
  expect_length(gauss_legendre(4)$x, 4L)                            # computed again once inactive
})

test_that("the DEV artifact verifies and reproduces the recorded H = 25 grid hashes on this platform", {
  quad_rules_deactivate()
  dev <- quad_rules_dev(ENGINE)
  man <- dev$manifest
  expect_identical(dev$sha256, man$sha256)
  expect_match(man$status, "^FROZEN with spec v9.8")
  quad_rules_activate(dev, required = TRUE)
  on.exit(quad_rules_deactivate(), add = TRUE)
  e <- man$expected_grid_hash
  gp <- p0_grid(e$H, rule = e$primary$definition$rule, n_c = e$primary$definition$n_c, n_psu = e$primary$definition$n_psu)
  gs <- p0_grid(e$H, rule = e$second$definition$rule, n_c = e$second$definition$n_c, n_psu = e$second$definition$n_psu)
  expect_identical(p0_grid_hash(gp), e$primary$sha256)
  expect_identical(p0_grid_hash(gs), e$second$sha256)
})

test_that("the DEV artifact reproduces the recorded grid hashes at every DEV H (5, 25, 50, 100, 400)", {
  quad_rules_deactivate()
  dev <- quad_rules_dev(ENGINE)
  quad_rules_activate(dev, required = TRUE)
  on.exit(quad_rules_deactivate(), add = TRUE)
  e <- dev$manifest$expected_grid_hash
  hs <- dev$manifest$expected_grid_hashes
  expect_identical(as.integer(hs$H), c(5L, 25L, 50L, 100L, 200L, 400L, 800L))   # H = 800 is the m = 1600 rung
  expect_identical(hs$primary_sha256[hs$H == e$H], e$primary$sha256)
  expect_identical(hs$second_sha256[hs$H == e$H], e$second$sha256)
  for (i in seq_len(nrow(hs))) {
    gp <- p0_grid(hs$H[i], rule = e$primary$definition$rule, n_c = e$primary$definition$n_c, n_psu = e$primary$definition$n_psu)
    gs <- p0_grid(hs$H[i], rule = e$second$definition$rule, n_c = e$second$definition$n_c, n_psu = e$second$definition$n_psu)
    expect_identical(p0_grid_hash(gp), hs$primary_sha256[i])
    expect_identical(p0_grid_hash(gs), hs$second_sha256[i])
  }
})
