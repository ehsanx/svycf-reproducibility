# Saved quadrature-rule artifacts: build, save, load with checksum verification, and
# activate for R/quadrature.R's store ("save the
# actual nodes and weights with a checksum, then load that artifact").
#
# The saved artifact grids/quad-rules-dev-v1.tsv is the one the recorded runs use;
# nothing here changes it.
#
# Format (one TSV, LF line endings, UTF-8/ASCII):
#   # svytmle-eif-clt quadrature rules v1        <- line 1, fixed
#   rule  n  i  x_hex  w_hex  x_dec  w_dec        <- header
#   one row per node of each BASE rule (Gauss-Legendre and Clenshaw-Curtis on [-1,1]
#   before scaling; Gauss-Hermite as used). x_hex/w_hex are C99 hex floats ("%a"),
#   exact. x_dec/w_dec are "%.17g", also exact, kept for reading; load requires the
#   two to agree bit for bit.
# The artifact's identity is the sha256 of the file's bytes, recorded in grids/MANIFEST.json.
#
# Requires R/quadrature.R (.gh_base, .gl_base, .cc_base, .QUAD_STORE) and digest.

QUAD_RULES_HEADER <- "# svytmle-eif-clt quadrature rules v1"
QUAD_RULES_COLS <- c("rule", "n", "i", "x_hex", "w_hex", "x_dec", "w_dec")

# Base rules for a spec data.frame(rule, n), computed on this platform.
quad_rules_build <- function(spec) {
  out <- list()
  for (j in seq_len(nrow(spec))) {
    rule <- spec$rule[j]; n <- as.integer(spec$n[j])
    b <- switch(rule, "gauss-hermite" = .gh_base(n), "gauss-legendre" = .gl_base(n),
                "clenshaw-curtis" = .cc_base(n), stop("unknown rule: ", rule))
    out[[paste0(rule, ":", n)]] <- b
  }
  out
}

quad_rules_save <- function(rules, path) {
  rows <- unlist(lapply(names(rules), function(key) {
    rn <- strsplit(key, ":", fixed = TRUE)[[1]]
    b <- rules[[key]]
    sprintf("%s\t%s\t%d\t%a\t%a\t%.17g\t%.17g", rn[1], rn[2], seq_along(b$x), b$x, b$w, b$x, b$w)
  }), use.names = FALSE)
  con <- file(path, open = "wb")                    # "wb": no CRLF translation on Windows
  writeLines(c(QUAD_RULES_HEADER, paste(QUAD_RULES_COLS, collapse = "\t"), rows), con, sep = "\n")
  close(con)                                        # flush before hashing the bytes
  invisible(digest::digest(file = path, algo = "sha256"))
}

# Load and verify. sha256 = the expected file checksum (from grids/MANIFEST.json);
# a mismatch, a malformed row or a hex/decimal disagreement is a hard error.
quad_rules_load <- function(path, sha256 = NULL) {
  if (!is.character(sha256) || length(sha256) != 1L || !grepl("^[0-9a-f]{64}$", sha256))
    stop("quad_rules_load: an expected sha256 (64 hex digits) is required", call. = FALSE)
  got <- digest::digest(file = path, algo = "sha256")
  if (!identical(got, sha256))
    stop("quadrature artifact ", basename(path), " has sha256 ", got, ", expected ", sha256, call. = FALSE)
  l <- readLines(path, warn = FALSE, encoding = "UTF-8")
  if (length(l) < 3L || l[1] != QUAD_RULES_HEADER || l[2] != paste(QUAD_RULES_COLS, collapse = "\t"))
    stop("quadrature artifact ", basename(path), ": bad header", call. = FALSE)
  f <- strsplit(l[-(1:2)], "\t", fixed = TRUE)
  if (any(lengths(f) != 7L)) stop("quadrature artifact: every row needs 7 fields", call. = FALSE)
  m <- do.call(rbind, f)
  x <- as.numeric(m[, 4]); w <- as.numeric(m[, 5])
  if (anyNA(x) || anyNA(w) || !identical(x, as.numeric(m[, 6])) || !identical(w, as.numeric(m[, 7])))
    stop("quadrature artifact: hex and decimal values disagree", call. = FALSE)
  key <- paste0(m[, 1], ":", m[, 2])
  rules <- lapply(split(seq_len(nrow(m)), factor(key, levels = unique(key))), function(i) {
    if (!identical(as.integer(m[i, 3]), seq_along(i)) || length(i) != as.integer(m[i[1], 2]))
      stop("quadrature artifact: rule ", key[i[1]], " rows are incomplete or out of order", call. = FALSE)
    list(x = x[i], w = w[i])
  })
  list(rules = rules, sha256 = got, path = normalizePath(path, winslash = "/"))
}

quad_rules_activate <- function(loaded, required = TRUE) {
  .QUAD_STORE$rules <- loaded$rules
  .QUAD_STORE$sha256 <- loaded$sha256
  .QUAD_STORE$required <- isTRUE(required)
  invisible(loaded$sha256)
}

quad_rules_deactivate <- function() {
  rm(list = ls(.QUAD_STORE, all.names = TRUE), envir = .QUAD_STORE)
  invisible(NULL)
}

quad_rules_active <- function() if (is.null(.QUAD_STORE$rules)) NA_character_ else .QUAD_STORE$sha256

# The DEV artifact recorded in grids/MANIFEST.json (path relative to the engine dir).
quad_rules_dev <- function(engine_dir) {
  man <- jsonlite::fromJSON(file.path(engine_dir, "grids", "MANIFEST.json"), simplifyVector = TRUE)
  c(quad_rules_load(file.path(engine_dir, "grids", man$file), man$sha256), list(manifest = man))
}
