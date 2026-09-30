source("config.R")
need_packages(c("scoringRules", "digest"))

log_sum_exp <- function(x) {
  z <- max(x)
  z + log(sum(exp(x-z)))
}

forecast_rows <- function(z, d) {
  f <- prepare_fold(d, z$fold)
  stopifnot(identical(z$row_id, d$row_id[f$test]))
  n <- length(f$test)
  levels <- c(.50, .80, .90, .95)
  probs <- sort(unique(c(.5, (1-levels)/2, (1+levels)/2)))
  qnames <- paste0("q", sprintf("%03d", round(probs*1000)))
  out <- data.frame(model = z$model, fold = z$fold, row_id = z$row_id,
                    recipe_id = d$recipe_id[f$test], y = f$ytest)
  set.seed(cfg$seed + z$fold*10000L + sum(utf8ToInt(z$model)))
  out$mean <- out$variance <- out$nlpd <- out$pit <- out$crps <- out$negative_mass <- NA_real_
  out$crps_half_difference <- 0
  for (nm in qnames) out[[nm]] <- NA_real_
  for (i in seq_len(n)) {
    y <- f$ytest[i]
    if (z$type %in% c("normal_mixture", "t_mixture")) {
      mu <- z$mu[, i]*f$ys + f$yc
      sig <- if (is.matrix(z$sigma)) z$sigma[, i]*f$ys else z$sigma*f$ys
      w <- if (is.null(z$weights)) rep(1/length(mu), length(mu)) else z$weights
      stopifnot(length(mu) == length(sig), length(w) == length(mu), all(sig > 0))
      idx <- sample.int(length(mu), cfg$prediction_draws, replace = TRUE, prob = w)
      if (z$type == "normal_mixture") {
        lp <- dnorm(y, mu, sig, log = TRUE)
        pit <- pnorm(y, mu, sig)
        neg <- pnorm(0, mu, sig)
        draws <- rnorm(cfg$prediction_draws, mu[idx], sig[idx])
        within <- sig^2
      } else {
        lp <- dt((y-mu)/sig, z$df, log = TRUE) - log(sig)
        pit <- pt((y-mu)/sig, z$df)
        neg <- pt(-mu/sig, z$df)
        draws <- mu[idx] + sig[idx]*rt(cfg$prediction_draws, z$df[idx])
        within <- sig^2*z$df/(z$df-2)
      }
      out$mean[i] <- sum(w*mu)
      out$variance[i] <- sum(w*(within + (mu-out$mean[i])^2))
      out$nlpd[i] <- -log_sum_exp(log(w) + lp)
      out$pit[i] <- sum(w*pit)
      out$negative_mass[i] <- sum(w*neg)
      out$crps[i] <- scoringRules::crps_sample(y, draws, method = "edf", show_messages = FALSE)
      h <- cfg$prediction_draws/2
      out$crps_half_difference[i] <- scoringRules::crps_sample(y, draws[seq_len(h)], show_messages = FALSE) -
        scoringRules::crps_sample(y, draws[h+seq_len(h)], show_messages = FALSE)
      out[i, qnames] <- quantile(draws, probs, names = FALSE)
    } else {
      mu <- z$mu[i]*f$ys + f$yc
      sig <- z$sigma[i]*f$ys
      stopifnot(is.finite(mu), is.finite(sig), sig > 0)
      out$mean[i] <- mu
      if (z$type == "normal") {
        out$variance[i] <- sig^2
        out$nlpd[i] <- -dnorm(y, mu, sig, log = TRUE)
        out$pit[i] <- pnorm(y, mu, sig)
        out$negative_mass[i] <- pnorm(0, mu, sig)
        out$crps[i] <- scoringRules::crps_norm(y, mu, sig)
        out[i, qnames] <- qnorm(probs, mu, sig)
      } else {
        out$variance[i] <- sig^2*z$df/(z$df-2)
        out$nlpd[i] <- -dt((y-mu)/sig, z$df, log = TRUE) + log(sig)
        out$pit[i] <- pt((y-mu)/sig, z$df)
        out$negative_mass[i] <- pt(-mu/sig, z$df)
        out$crps[i] <- scoringRules::crps_t(y, df = z$df, location = mu, scale = sig)
        out[i, qnames] <- mu + sig*qt(probs, z$df)
      }
    }
  }
  out$squared_error <- (out$mean-out$y)^2
  out$absolute_error <- abs(out$q500-out$y)
  for (lev in levels) {
    lo <- paste0("q", sprintf("%03d", round(1000*(1-lev)/2)))
    hi <- paste0("q", sprintf("%03d", round(1000*(1+lev)/2)))
    out[[paste0("cover", lev*100)]] <- as.integer(out$y >= out[[lo]] & out$y <= out[[hi]])
    out[[paste0("width", lev*100)]] <- out[[hi]] - out[[lo]]
  }
  stopifnot(all(is.finite(as.matrix(out[, vapply(out, is.numeric, logical(1))]))),
            all(out$crps >= 0), all(out$pit >= 0 & out$pit <= 1), all(out$width90 > 0))
  out
}

if (sys.nframe() == 0L) {
  d <- read_data()
  files <- list.files(root, pattern = "^fit-.*[.]rds$", full.names = TRUE)
  primary <- c("OLS","Ridge","Lasso","Horseshoe","BMA","GAM","BART-1000","GP","NN")
  expected <- unlist(lapply(c(primary,"BMA-g0.1","BMA-g10"),function(m)
    sprintf("fit-%s-fold%d.rds",m,1:5)))
  stopifnot(all(expected %in% basename(files)))
  files <- files[basename(files) %in% expected]
  for (fold in 1:5) {
    for (model in c("Ridge", "Lasso", "Horseshoe", "BART-1000")) {
      sampled <- readRDS(fit_path(model, fold))
      stopifnot(isTRUE(tail(sampled$diagnostics,1)[[1]]$pass))
    }
  }
  output <- timing <- list()
  diagnostic <- list()
  for (file in files) {
    z <- readRDS(file)
    source_files <- if (z$model == "BART-1000") c("config.R", "models.R", "fit-bart.R") else
      if (grepl("^BMA-g", z$model)) c("config.R", "sensitivity.R") else
        c("config.R", "models.R", "fit-models.R")
    check_cache(z, cache_stamp(source_files))
    cat("SCORING", z$model, z$fold, "\n")
    tm <- proc.time()[[3]]
    pred <- forecast_rows(z, d)
    seconds <- proc.time()[[3]]-tm
    output[[length(output)+1L]] <- pred
    timing[[length(timing)+1L]] <- data.frame(model = z$model, fold = z$fold,
      fit_seconds = z$fit_seconds, predictive_seconds = seconds +
        if (is.null(z$construction_seconds)) 0 else z$construction_seconds,
      result_size_bytes = as.numeric(object.size(z)))
    diagnostic[[paste(z$model, z$fold, sep = "-")]] <- z$diagnostics
    if (z$model == "GAM") {
      z$model <- "GAM-conditional"
      z$sigma <- z$sigma_conditional
      output[[length(output)+1L]] <- forecast_rows(z, d)
    }
  }
  predictions <- do.call(rbind, output)
  predictions$role <- ifelse(predictions$model=="BART-1000","primary",
    ifelse(grepl("^BART",predictions$model),"diagnostic_only",
    ifelse(predictions$model %in% c("GAM-conditional","BMA-g0.1","BMA-g10"),"sensitivity","primary"))
  )
  write.csv(predictions, path_result("predictions.csv"), row.names = FALSE)
  write.csv(do.call(rbind, timing), path_result("timings.csv"), row.names = FALSE)
  saveRDS(diagnostic, path_result("diagnostics.rds"))
  writeLines(capture.output(str(diagnostic, max.level = 4)), path_result("diagnostics.txt"))
  cat("PREDICTIVE EVALUATION COMPLETE\n")
}
