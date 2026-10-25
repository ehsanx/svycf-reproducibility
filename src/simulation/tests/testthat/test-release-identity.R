# The reproduction identity check of the summaries (release_identity_check in R/grade_run.R). Output made under this
# release's identity passes only when its runs used this release's code and that code has the executable code the
# execution check passed (same bytes, or comments edited only), apart from the adapted files; output made at a commit
# keeps the original check. Constructed fixtures only: no population, sample or fit.

for (f in c("rng_streams.R", "run_mode.R", "learner_rd.R", "estimators.R", "tmle.R", "learner_rb.R", "controls.R",
            "grading.R", "grade_plan.R", "grade_io.R", "grade_run.R"))
  source(test_path("..", "..", "R", f), local = TRUE)

ID <- "release:test-release+000000000000"
fixture <- function() {
  d <- tempfile("release-"); dir.create(file.path(d, "src", "simulation", "R"), recursive = TRUE)
  dir.create(file.path(d, "config"))
  writeLines(c("# a comment", "a <- 1  # inline"), file.path(d, "src", "simulation", "R", "a.R"))
  writeLines("b <- 2", file.path(d, "src", "simulation", "R", "run_mode.R"))
  p <- function(x) file.path(d, "src", "simulation", x)
  h <- function(x) digest::digest(file = p(x), algo = "sha256")
  rec <- list(release_id = sub("^release:", "", ID),
              execution_check = list(commit = "c0ffee", adapted = list("R/run_mode.R"),
                                     engine_files = list("R/a.R" = h("R/a.R"), "R/run_mode.R" = strrep("0", 64)),
                                     engine_code = list("R/a.R" = code_sha256(p("R/a.R")))))
  jsonlite::write_json(rec, file.path(d, "config", "execution-identity.json"), auto_unbox = TRUE)
  runs <- function() rep(list(list("R/a.R" = h("R/a.R"), "R/run_mode.R" = h("R/run_mode.R"))), 3)
  list(dir = d, path = p, runs = runs)
}

test_that("code_sha256 ignores comments and blank lines but not code or strings", {
  f <- tempfile(fileext = ".R")
  writeLines(c("x <- \"# not a comment\"  # a comment", "", "y <- 2"), f); h1 <- code_sha256(f)
  writeLines(c("# new header", "x <- \"# not a comment\"", "y <- 2  # another"), f); expect_identical(code_sha256(f), h1)
  writeLines(c("x <- \"# changed string\"", "y <- 2"), f); expect_false(identical(code_sha256(f), h1))
  writeLines(c("x <- \"# not a comment\"", "y <- 3"), f); expect_false(identical(code_sha256(f), h1))
})

test_that("output made at a commit keeps the original check", {
  fx <- fixture()
  expect_false(release_identity_check("c0ffee", fx$runs(), fx$dir, "c0ffee")$applies)
  expect_false(release_identity_check(character(0), fx$runs(), fx$dir, "c0ffee")$applies)
})

test_that("output made under this release passes when every run used this release's code", {
  fx <- fixture()
  r <- release_identity_check(ID, fx$runs(), fx$dir, "c0ffee")
  expect_true(r$applies); expect_true(r$ok)
  expect_match(r$detail, "3 run folders; 1 engine files with the code of c0ffee (0 with comments edited), adapted: R/run_mode.R",
               fixed = TRUE)
})

test_that("a comment-only edit keeps the check; a change to code or strings fails it", {
  fx <- fixture()
  writeLines(c("# a rewritten comment", "", "a <- 1"), fx$path("R/a.R"))
  r <- release_identity_check(ID, fx$runs(), fx$dir, "c0ffee")              # runs made with the edited file
  expect_true(r$ok); expect_match(r$detail, "(1 with comments edited)", fixed = TRUE)
  writeLines(c("# a rewritten comment", "a <- 2"), fx$path("R/a.R"))
  r <- release_identity_check(ID, fx$runs(), fx$dir, "c0ffee")
  expect_false(r$ok); expect_match(r$detail, "differ from the execution-check commit: R/a.R", fixed = TRUE)
})

test_that("it fails on another identity, another execution-check commit, runs of other code or missing records", {
  fx <- fixture()
  expect_false(release_identity_check("release:other", fx$runs(), fx$dir, "c0ffee")$ok)
  expect_false(release_identity_check(c(ID, "c0ffee"), fx$runs(), fx$dir, "c0ffee")$ok)     # mixed identities
  expect_match(release_identity_check(ID, fx$runs(), fx$dir, "beef00")$detail, "execution-check commit c0ffee")
  runs <- fx$runs(); runs[[2]][["R/a.R"]] <- strrep("1", 64)
  expect_match(release_identity_check(ID, runs, fx$dir, "c0ffee")$detail, "other code for: R/a.R", fixed = TRUE)
  runs <- fx$runs(); runs[[3]][["R/z.R"]] <- strrep("2", 64)
  expect_match(release_identity_check(ID, runs, fx$dir, "c0ffee")$detail, "not in the release record: R/z.R", fixed = TRUE)
  expect_match(release_identity_check(ID, c(fx$runs(), list(NULL)), fx$dir, "c0ffee")$detail, "1 run folders recorded no")
  expect_false(release_identity_check(ID, fx$runs(), NULL, "c0ffee")$ok)                      # no release record
})
