source("config.R")
need_packages(c("readxl", "digest"))
stopifnot(digest::digest(file = "Concrete_Data.xls", algo = "sha256") ==
  "710076c66b9ca3f8050e7942f3dcbdbe04013534daeb0077ffd3079a52d8e0c4")
d <- as.data.frame(readxl::read_xls("Concrete_Data.xls"))
names(d) <- c("cement", "slag", "flyash", "water", "superplasticizer",
              "coarse", "fine", "age", "strength")
stopifnot(nrow(d) == 1030L, ncol(d) == 9L, !anyNA(d), all(d$age > 0))
key <- apply(d[1:7], 1, function(x) paste(sprintf("%.6f", x), collapse = "|"))
d$recipe_id <- match(key, unique(key))
d$row_id <- seq_len(nrow(d))
sizes <- tabulate(d$recipe_id)
set.seed(cfg$seed)
order_groups <- order(-sizes, runif(length(sizes)))
assignment <- integer(length(sizes))
loads <- integer(cfg$folds)
for (g in order_groups) {
  candidates <- which(loads == min(loads))
  f <- candidates[sample.int(length(candidates), 1L)]
  assignment[g] <- f
  loads[f] <- loads[f] + sizes[g]
}
d$fold <- assignment[d$recipe_id]
stopifnot(isTRUE(all.equal(d, read_data(), check.attributes = FALSE, tolerance = 1e-12)),
          length(sizes) == 428L, all(loads == 206L))
cat("PASS: the original UCI data reproduce all 1,030 rows and the five grouped folds.\n")
