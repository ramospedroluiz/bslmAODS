source("config.R")
need_packages(required)
stamp <- cache_stamp(c("config.R", "sensitivity.R"))
d <- read_data()
g_factors <- c(.1,10)
for (fold in seq_len(cfg$folds)) {
  f <- prepare_fold(d,fold)
  for (factor in g_factors) {
    model <- paste0("BMA-g",factor)
    target <- fit_path(model, fold)
    if (file.exists(target)) {
      check_cache(readRDS(target), stamp)
      cat("CACHED", model, fold, "\n")
      next
    }
    set.seed(cfg$seed+fold*10000L+as.integer(factor*100))
    g <- nrow(f$xtrain)*factor
    tm <- proc.time()[[3]]
    obj <- BAS::bas.lm(y~.,data.frame(y=f$ytrain,f$xtrain),prior="g-prior",alpha=g,
                       modelprior=BAS::uniform(),n.models=256,method="BAS")
    fit_seconds <- proc.time()[[3]]-tm
    tm <- proc.time()[[3]]
    pr <- predict(obj,as.data.frame(f$xtest),se.fit=TRUE)
    stopifnot(length(obj$which)==256L,abs(sum(pr$postprobs)-1)<1e-12)
    z <- list(model=model,fold=fold,row_id=d$row_id[f$test],type="t_mixture",
               mu=pr$Ypred,sigma=pr$se.pred,df=pr$df,weights=pr$postprobs,
               fit_seconds=fit_seconds,construction_seconds=proc.time()[[3]]-tm,
               settings=list(g=g,model_prior="uniform"),
               diagnostics=list(models=length(obj$which),posterior_sum=sum(pr$postprobs)),
               preprocessing=f$preprocessing,fit=obj,
               cache_stamp=stamp)
    saveRDS(z,target,compress=FALSE)
    cat("SAVED",model,fold,"g",g,"\n")
  }
}
