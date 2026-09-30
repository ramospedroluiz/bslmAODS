source("config.R")
compare_csv <- function(name, keys) {
  expected <- read.csv(paste0("reference-", name))
  actual <- read.csv(path_result(name))
  expected <- expected[do.call(order, expected[keys]), , drop = FALSE]
  actual <- actual[do.call(order, actual[keys]), , drop = FALSE]
  rownames(expected) <- rownames(actual) <- NULL
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-10)))
  cat("PASS: reproduced", name, "\n")
}
compare_csv("summary.csv", "model")
compare_csv("metric-intervals.csv", c("model", "metric"))
compare_csv("paired-differences.csv", c("first", "second", "metric"))
