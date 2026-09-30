# Run from this directory: Rscript --vanilla run.R [summary|check|test|fit|score]
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "summary"
if (length(args) > 1L || !mode %in% c("summary", "check", "test", "fit", "score"))
  stop("Usage: Rscript --vanilla run.R [summary|check|test|fit|score]")
if (!file.exists("concrete.csv")) stop("Run from the bslmAODS directory.")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
Sys.setenv(BSLM_INPUT = if (mode %in% c("fit", "score")) "new" else "reference")
steps <- switch(mode,
  summary = c("summarize.R", "verify.R"),
  check = c("prepare.R", "summarize.R", "verify.R"),
  test = c("prepare.R", "summarize.R", "verify.R", "tests.R", "test-bart.R"),
  fit = c("prepare.R", "fit-models.R", "sensitivity.R", "fit-bart.R",
          "evaluate.R", "summarize.R", "tests.R"),
  score = c("evaluate.R", "summarize.R", "tests.R"))
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
for (script in steps) {
  cat("\n===", script, "===\n")
  status <- system2(rscript, c("--vanilla", shQuote(script)))
  if (status != 0L) stop(script, " failed; no further steps were run.")
}
cat("\nCompleted:", mode, "\nAll generated analysis files are in this directory.\n")
