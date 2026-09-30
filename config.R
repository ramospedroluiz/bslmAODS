options(stringsAsFactors = FALSE, timeout = 600, width = 110)
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
library_path <- Sys.getenv("BSLM_LIBRARY", tools::R_user_dir("bslmAODS", "data"))
if (dir.exists(library_path)) .libPaths(c(library_path, .libPaths()))
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

cfg <- list(seed = 20260909L, folds = 5L, chains = 4L, warmup = 2000L,
            draws = 2000L, bart_trees = 1000L, bart_k = 2,
            gam_k = 10L, nn_neurons = 8L, nn_epochs = 3000L,
            gp_starts = 3L, gp_maxit = 150L,
            prediction_draws = 8000L, bootstrap = 2000L)
primary_models <- c("OLS", "Ridge", "Lasso", "Horseshoe", "BMA", "GAM", "BART-1000", "GP", "NN")
required <- c("bayesreg", "dbarts", "BAS", "mgcv", "DiceKriging", "brnn",
              "posterior", "readxl", "jsonlite", "ggplot2", "scoringRules", "digest", "foreach")
need_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "),
                            ". Run Rscript --vanilla setup.R first.")
}
path_result <- function(name) file.path(root, paste0("output-", name))
path_input <- function(name) {
  prefix <- if (Sys.getenv("BSLM_INPUT", "reference") == "new") "output-" else "reference-"
  file.path(root, paste0(prefix, name))
}
fit_path <- function(model, fold) file.path(root, sprintf("fit-%s-fold%d.rds", model, fold))
read_data <- function() read.csv(file.path(root, "concrete.csv"))

prepare_fold <- function(d, fold) {
  train <- which(d$fold != fold)
  test <- which(d$fold == fold)
  x <- as.matrix(d[, c("cement", "slag", "flyash", "water", "superplasticizer",
                      "coarse", "fine", "age")])
  x[, "age"] <- log(x[, "age"])
  colnames(x)[8] <- "log_age"
  center <- colMeans(x[train, , drop = FALSE])
  scale <- apply(x[train, , drop = FALSE], 2, sd)
  x <- sweep(sweep(x, 2, center), 2, scale, "/")
  yc <- mean(d$strength[train])
  ys <- sd(d$strength[train])
  list(train = train, test = test, xtrain = x[train, , drop = FALSE],
       xtest = x[test, , drop = FALSE], ytrain = (d$strength[train] - yc)/ys,
       ytest = d$strength[test], yc = yc, ys = ys,
       preprocessing = list(center = center, scale = scale, response_center = yc, response_scale = ys))
}

# Only reuse new checkpoints when their data, code and numerical environment match.
cache_stamp <- function(files) {
  need_packages(c("digest", required))
  list(data = digest::digest(file = "concrete.csv", algo = "sha256"),
       source = setNames(vapply(files, function(f) digest::digest(file = f, algo = "sha256"),
                               character(1)), files),
       R = R.version.string, platform = R.version$platform,
       packages = setNames(vapply(required, function(p) as.character(packageVersion(p)),
                                  character(1)), required))
}
check_cache <- function(z, stamp) {
  if (!identical(z$cache_stamp, stamp))
    stop("Checkpoint does not match the current data, scripts or environment. ",
         "Move previous fit-*.rds and chain-*.rds files aside before refitting.")
}
