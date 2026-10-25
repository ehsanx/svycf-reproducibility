# make_tasks.R -- turn the FROZEN protocol into the SLURM array's task table (tasks.tsv).
# One line per work unit:  <mode> <survey> <f or -> <range or ->
#   overlap   <survey> - -          census overlap + environment (cheap)
#   reference <survey> - a:b        census pseudo-reference splits a..b
#   run       <survey> <f> a:b      resamples a..b at stage-1 fraction f
# Chunk sizes come from the environment (defaults sized for 32-core tasks).
# Nothing here opens data; it only reads the protocol.
# Run from the repository root:  Rscript scripts/application/make_tasks.R > outputs/application/tasks.tsv
suppressPackageStartupMessages(library(yaml))
txt <- readLines(file.path("config", "application", "protocol.yaml"), warn = FALSE)
fm  <- which(trimws(txt) == "---")
p   <- yaml.load(paste(txt[(fm[1] + 1):(fm[2] - 1)], collapse = "\n"))
if (!identical(p$status, "FROZEN")) stop("the protocol is not FROZEN; refusing to build tasks")
rep_chunk <- as.integer(Sys.getenv("SPA_REP_CHUNK", "64"))     # 2 replicates per core on a 32-core task
ref_chunk <- as.integer(Sys.getenv("SPA_REF_CHUNK", "20"))   # one task per survey (REV4-PANEL-RECOMMENDATION.md:51)
surveys <- c(p$survey, unlist(p$replication_surveys))
fs <- unlist(p$design$f); R <- as.integer(p$estimation$R); Rr <- as.integer(p$estimation$R_ref)
ranges <- function(n, k) { a <- seq(1L, n, by = k); sprintf("%d:%d", a, pmin(a + k - 1L, n)) }
out <- character(0)
for (s in surveys) {
  out <- c(out, sprintf("overlap\t%s\t-\t-", s))
  out <- c(out, sprintf("reference\t%s\t-\t%s", s, ranges(Rr, ref_chunk)))
  for (f in fs) out <- c(out, sprintf("run\t%s\t%s\t%s", s, format(f), ranges(R, rep_chunk)))
}
writeLines(out)
message(length(out), " tasks: ", length(surveys), " surveys x (1 overlap + ", length(ranges(Rr, ref_chunk)),
        " reference chunks + ", length(fs), " f x ", length(ranges(R, rep_chunk)), " run chunks)")
