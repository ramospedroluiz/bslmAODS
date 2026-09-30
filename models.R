source("config.R")
need_packages(required)

diagnose_chains <- function(chains) {
  a <- array(NA_real_, c(nrow(chains[[1]]$monitor), length(chains), ncol(chains[[1]]$monitor)),
             dimnames = list(NULL, NULL, colnames(chains[[1]]$monitor)))
  for (j in seq_along(chains)) a[, j, ] <- chains[[j]]$monitor
  z <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(a),
                                              "rhat", "ess_bulk", "ess_tail"))
  list(detail = z, max_rhat = max(z$rhat), min_bulk = min(z$ess_bulk),
       min_tail = min(z$ess_tail), pass = all(is.finite(z$rhat)) &&
         max(z$rhat) <= 1.01 && min(z$ess_bulk) >= 400 && min(z$ess_tail) >= 400)
}

fit_sampled <- function(f, model, fold, k = cfg$bart_k) {
  stopifnot(model %in% c("Ridge", "Lasso", "Horseshoe"))
  histories <- list()
  elapsed <- 0
  for (attempt in 1:3) {
    multiplier <- c(1L, 4L, 8L)[attempt]
    nd <- cfg$draws * multiplier
    nb <- cfg$warmup * multiplier
    chains <- vector("list", cfg$chains)
    seeds <- cfg$seed + fold*10000L + match(model, c("Ridge", "Lasso", "Horseshoe", "BART"))*1000L +
      as.integer(k*100) + attempt*10L + seq_len(cfg$chains)
    for (ch in seq_len(cfg$chains)) {
      set.seed(seeds[ch])
      cat("SAMPLING", model, "fold", fold, "k", k, "attempt", attempt, "chain", ch, "draws", nd, "\n")
      tm <- proc.time()[[3]]
      obj <- bayesreg::bayesreg(y ~ ., data.frame(y = f$ytrain, f$xtrain),
                   model = "normal", prior = tolower(model), n.samples = nd,
                   burnin = nb, thin = 1, n.cores = 1)
      mu <- sweep(t(obj$beta) %*% t(f$xtest), 1, as.numeric(obj$beta0), "+")
      sigma <- sqrt(obj$sigma2)
      monitor <- cbind(sigma = sigma, t(obj$beta), intercept = as.numeric(obj$beta0),
                       log_tau2 = log(obj$tau2), mu)
      info <- list(prior = obj$prior, model = obj$model)
      elapsed <- elapsed + proc.time()[[3]] - tm
      colnames(monitor) <- make.unique(c(colnames(monitor)[seq_len(ncol(monitor)-ncol(mu))],
                                         paste0("f_", seq_len(ncol(mu)))))
      chains[[ch]] <- list(mu = mu, sigma = sigma, monitor = monitor)
    }
    dg <- diagnose_chains(chains)
    histories[[attempt]] <- c(list(attempt = attempt, draws = nd, warmup = nb, seeds = seeds), dg)
    cat("DIAGNOSTIC", model, fold, "Rhat", dg$max_rhat, "bulk", dg$min_bulk,
        "tail", dg$min_tail, "pass", dg$pass, "\n")
    if (dg$pass) break
  }
  list(type = "normal_mixture", mu = do.call(rbind, lapply(chains, `[[`, "mu")),
       sigma = unlist(lapply(chains, `[[`, "sigma")),
       fit_seconds = elapsed, diagnostics = histories, settings = info,
       preprocessing = f$preprocessing)
}

nn_jacobian <- function(x, fit) {
  brnn:::jacobian(as.double(x), nrow(x), ncol(x), fit$npar, fit$theta, fit$neurons, cores = 1)
}

fit_other <- function(f, model, fold) {
  set.seed(cfg$seed + fold*10000L + match(model, c("BMA", "GAM", "GP", "NN", "OLS"))*100L)
  dat <- data.frame(y = f$ytrain, f$xtrain)
  newdat <- as.data.frame(f$xtest)
  tm <- proc.time()[[3]]
  if (model == "BMA") {
    obj <- BAS::bas.lm(y ~ ., dat, prior = "g-prior", alpha = nrow(dat),
                        modelprior = BAS::uniform(), n.models = 256, method = "BAS")
  } else if (model == "GAM") {
    form <- as.formula(paste("y ~", paste(sprintf("s(%s, bs='cr', k=%d)",
                                             colnames(f$xtrain), cfg$gam_k), collapse = "+")))
    obj <- mgcv::gam(form, data = dat, method = "REML", select = TRUE)
  } else if (model == "GP") {
    foreach::registerDoSEQ()
    obj <- DiceKriging::km(~1, design = as.data.frame(f$xtrain), response = f$ytrain,
                           covtype = "matern5_2", nugget.estim = TRUE, nugget = .05,
                           lower = rep(.1, 8), upper = rep(20, 8), parinit = rep(2, 8),
                           multistart = cfg$gp_starts, control = list(trace = FALSE, maxit = cfg$gp_maxit))
  } else if (model == "NN") {
    obj <- brnn::brnn(f$xtrain, f$ytrain, neurons = cfg$nn_neurons,
                     normalize = FALSE, epochs = cfg$nn_epochs, change = 1e-7, cores = 1)
  } else if (model == "OLS") obj <- lm(y ~ ., dat)
  elapsed <- proc.time()[[3]] - tm
  tm <- proc.time()[[3]]
  if (model == "BMA") {
    pr <- predict(obj, newdat, se.fit = TRUE)
    stopifnot(length(obj$which) == 256L, abs(sum(pr$postprobs)-1) < 1e-12)
    ans <- list(type = "t_mixture", mu = pr$Ypred, sigma = pr$se.pred,
                df = pr$df, weights = pr$postprobs,
                diagnostics = list(models = length(obj$which), posterior_sum = sum(pr$postprobs)),
                settings = list(g = nrow(dat), model_prior = "uniform", df = unique(pr$df)))
  } else if (model == "GAM") {
    xp <- predict(obj, newdat, type = "lpmatrix")
    stopifnot(!is.null(obj$Vc))
    v <- vcov(obj, unconditional = TRUE)
    ans <- list(type = "normal", mu = as.numeric(xp %*% coef(obj)),
                sigma = sqrt(obj$sig2 + rowSums((xp %*% v)*xp)),
                sigma_conditional = sqrt(obj$sig2 + rowSums((xp %*% obj$Vp)*xp)),
                noise_variance = obj$sig2,
                diagnostics = list(converged = obj$converged, outer = obj$outer.info, edf = sum(obj$edf)),
                settings = list(k = cfg$gam_k, sp = obj$sp, covariance = "Vc"))
  } else if (model == "GP") {
    pr <- predict(obj, newdat, type = "UK", checkNames = TRUE)
    stopifnot(all(pr$sd^2 >= obj@covariance@nugget - 1e-8))
    ans <- list(type = "normal", mu = as.numeric(pr$mean), sigma = as.numeric(pr$sd),
                noise_variance = obj@covariance@nugget,
                diagnostics = list(logLik = obj@logLik, optimizer = obj@optim.method,
                  convergence_code_available = FALSE,
                  range = obj@covariance@range.val, nugget = obj@covariance@nugget),
                settings = list(covariance = "matern5_2", multistart = cfg$gp_starts,
                                lower = .1, upper = 20, maxit = cfg$gp_maxit))
  } else if (model == "NN") {
    j <- nn_jacobian(f$xtrain, obj)
    jp <- nn_jacobian(f$xtest, obj)
    precision <- 2*obj$beta*crossprod(j) + diag(2*obj$alpha, obj$npar)
    v <- chol2inv(chol(precision))
    ans <- list(type = "normal", mu = as.numeric(predict(obj, f$xtest)),
                sigma = sqrt(1/(2*obj$beta) + rowSums((jp %*% v)*jp)),
                noise_variance = 1/(2*obj$beta),
                diagnostics = list(reason = obj$reason, epoch = obj$epoch, alpha = obj$alpha,
                  beta = obj$beta, gamma = obj$gamma, precision_condition = kappa(precision)),
                settings = list(neurons = cfg$nn_neurons, epochs = cfg$nn_epochs,
                                approximation = "evidence/Gauss-Newton/linearized"))
  } else {
    pr <- predict(obj, newdat, se.fit = TRUE)
    ans <- list(type = "t", mu = as.numeric(pr$fit),
                sigma = sqrt(pr$se.fit^2 + pr$residual.scale^2), df = pr$df,
                diagnostics = list(rank = obj$rank), settings = list(formula = "y~all eight inputs"))
  }
  ans$fit_seconds <- elapsed
  ans$construction_seconds <- proc.time()[[3]] - tm
  ans$preprocessing <- f$preprocessing
  ans$fit <- obj
  ans
}
