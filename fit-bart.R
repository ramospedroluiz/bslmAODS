source("models.R")

pack_bart <- function(sampler,raw) {
  z <- dbarts:::packageBartResults(sampler,raw,NULL,NULL,TRUE,TRUE)
  list(mu=z$yhat.test,sigma=as.numeric(z$sigma))
}

new_bart_chain <- function(f,seed,warmup=8000L,draws=16000L) {
  set.seed(seed)
  time <- proc.time()[[3]]
  sampler <- dbarts::bart(f$xtrain,f$ytrain,f$xtest,ntree=1000L,k=2,
    usequants=FALSE,numcut=100L,nchain=1L,nthread=1L,
    nskip=warmup,ndpost=draws,keepevery=1L,keeptrainfits=FALSE,
    keeptrees=FALSE,verbose=FALSE,sampleronly=TRUE)
  raw <- sampler$run(warmup,draws,updateState=TRUE)
  seconds <- proc.time()[[3]]-time
  z <- pack_bart(sampler,raw)
  invisible(sampler$state)
  c(z,list(sampler=sampler,rng=.Random.seed,seed=seed,fit_seconds=seconds,warmup=warmup))
}

extend_bart_chain <- function(z,draws) {
  assign(".Random.seed",z$rng,envir=.GlobalEnv)
  time <- proc.time()[[3]]
  raw <- z$sampler$run(0L,draws,updateState=TRUE)
  z$fit_seconds <- z$fit_seconds+proc.time()[[3]]-time
  extra <- pack_bart(z$sampler,raw)
  z$mu <- rbind(z$mu,extra$mu)
  z$sigma <- c(z$sigma,extra$sigma)
  z$rng <- .Random.seed
  invisible(z$sampler$state)
  z
}

pending_bart_targets <- function(chains,histories=list(),budgets=c(16000L,32000L,64000L)) {
  have <- vapply(chains,function(z) if(is.null(z)) 0L else nrow(z$mu),integer(1))
  stopifnot(all(have %in% c(0L,budgets)))
  targets <- budgets[budgets>=max(have)]
  if (length(histories)) {
    last <- tail(histories,1)[[1]]
    if (all(have==last$draws)) targets <- targets[targets>last$draws]
  }
  targets
}

if (sys.nframe()==0L) {
  d <- read_data()
  stamp <- cache_stamp(c("config.R", "models.R", "fit-bart.R"))
  hashes <- stamp$source
  data_hash <- stamp$data
  for (fold in 1:5) {
    target <- fit_path("BART-1000", fold)
    histories <- list()
    if (file.exists(target)) {
      saved <- readRDS(target)
      check_cache(saved, stamp)
      if (tail(saved$diagnostics,1)[[1]]$pass) {cat("CACHED CONVERGED FOLD",fold,"\n"); next}
      histories <- saved$diagnostics
      rm(saved)
    }
    f <- prepare_fold(d,fold)
    chains <- vector("list",4)
    for (ch in 1:4) {
      checkpoint <- file.path(root,sprintf("chain-bart-fold%d-chain%d.rds",fold,ch))
      if (file.exists(checkpoint)) {
        chains[[ch]] <- readRDS(checkpoint)
        check_cache(chains[[ch]], stamp)
      }
    }
    if (length(histories)) stopifnot(!any(vapply(chains,is.null,logical(1))))
    for (draw_target in pending_bart_targets(chains,histories)) {
      for (ch in 1:4) {
        checkpoint <- file.path(root,sprintf("chain-bart-fold%d-chain%d.rds",fold,ch))
        before <- if(is.null(chains[[ch]])) 0L else nrow(chains[[ch]]$mu)
        sources <- if(before==0L) list() else chains[[ch]]$source_history
        if (before>0L && is.null(sources)) {
          sources <- list(list(first=1L,last=before,
            sha256=unname(chains[[ch]]$script_hashes["fit-bart.R"])))
        }
        cat("BART1000",fold,"CHAIN",ch,"TARGET",draw_target,format(Sys.time()),"\n")
        if (is.null(chains[[ch]])) chains[[ch]] <- new_bart_chain(f,202609100L+fold*100L+ch)
        have <- nrow(chains[[ch]]$mu)
        if (have<draw_target) chains[[ch]] <- extend_bart_chain(chains[[ch]],draw_target-have)
        if (nrow(chains[[ch]]$mu)>before) sources[[length(sources)+1L]] <-
          list(first=before+1L,last=nrow(chains[[ch]]$mu),sha256=unname(hashes["fit-bart.R"]))
        chains[[ch]]$source_history <- sources
        chains[[ch]]$script_hashes <- hashes
        chains[[ch]]$data_sha256 <- data_hash
        chains[[ch]]$cache_stamp <- stamp
        saveRDS(chains[[ch]],checkpoint,compress=FALSE)
      }
      stopifnot(length(unique(vapply(chains,function(z)nrow(z$mu),integer(1))))==1L,
        all(vapply(chains,function(z)length(z$sigma)==nrow(z$mu),logical(1))))
      nd <- nrow(chains[[1]]$mu)
      monitored <- lapply(chains,function(z) {
        monitor <- cbind(sigma=z$sigma,z$mu)
        colnames(monitor) <- c("sigma",paste0("f_",seq_len(ncol(z$mu))))
        list(monitor=monitor)
      })
      diagnostic <- diagnose_chains(monitored)
      histories[[length(histories)+1L]] <- c(list(draws=nd,warmup=8000L,
        seeds=vapply(chains,`[[`,numeric(1),"seed")),diagnostic)
      cat("DIAGNOSTIC",fold,"DRAWS",nd,"Rhat",diagnostic$max_rhat,
          "bulk",diagnostic$min_bulk,"tail",diagnostic$min_tail,"pass",diagnostic$pass,"\n")
      result <- list(model="BART-1000",fold=fold,row_id=d$row_id[f$test],
        type="normal_mixture",mu=do.call(rbind,lapply(chains,`[[`,"mu")),
        sigma=unlist(lapply(chains,`[[`,"sigma")),
        fit_seconds=sum(vapply(chains,`[[`,numeric(1),"fit_seconds")),
        diagnostics=histories,preprocessing=f$preprocessing,script_hashes=hashes,cache_stamp=stamp,
        chain_source_history=lapply(chains,`[[`,"source_history"),
        data_sha256=data_hash,settings=list(package="dbarts",version=as.character(packageVersion("dbarts")),
          ntree=1000,k=2,depth_base=.95,depth_power=2,numcut=100,usequants=FALSE,
          sigdf=3,sigquant=.9,continuation="same chains and RNG states",threads=1))
      saveRDS(result,target,compress=FALSE)
      rm(result,monitored)
      gc()
      if (diagnostic$pass) break
    }
    rm(chains)
    gc()
  }
  cat("BART1000 FITTING COMPLETE; NO RESPONSE SCORES CALCULATED\n")
}
