source("models.R")
d <- read_data()
data_file <- file.path(root,"concrete.csv")
stopifnot(isTRUE(all.equal(d,read.csv(data_file),check.attributes=FALSE,tolerance=1e-12)))
data_hash <- digest::digest(file=data_file,algo="sha256")
models <- c("Ridge","Lasso","Horseshoe","BMA","GAM","NN","OLS","GP")
stamp <- cache_stamp(c("config.R", "models.R", "fit-models.R"))
for (fold in 1:5) {
  f <- prepare_fold(d,fold)
  for (model in models) {
    file <- fit_path(model, fold)
    if (file.exists(file)) {
      z <- readRDS(file)
      check_cache(z, stamp)
      stopifnot(identical(z$row_id,d$row_id[f$test]),identical(z$preprocessing,f$preprocessing))
      cat("CACHED",model,fold,"\n")
      next
    }
    cat("START",model,fold,format(Sys.time()),"\n")
    z <- if (model %in% c("Ridge","Lasso","Horseshoe")) fit_sampled(f,model,fold) else fit_other(f,model,fold)
    z$model <- model
    z$fold <- fold
    z$row_id <- d$row_id[f$test]
    z$cache_stamp <- stamp
    z$data_sha256 <- data_hash
    saveRDS(z,file,compress=FALSE)
  }
}
cat("NON-TREE PRIMARY FITS COMPLETE; RUN fit-bart.R FOR THE PRIMARY BART FITS.\n")
