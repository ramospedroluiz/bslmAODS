source("models.R")
source("evaluate.R")
d <- read_data()
stopifnot(nrow(d) == 1030L, length(unique(d$recipe_id)) == 428L,
          nrow(unique(d[1:7])) == 428L,
          all(table(d$fold) == 206L), !anyNA(d))
stopifnot(all(vapply(split(d[1:7], d$recipe_id), function(x) nrow(unique(x)) == 1L, logical(1))))
for (fold in 1:5) {
  f <- prepare_fold(d, fold)
  stopifnot(length(intersect(d$recipe_id[f$train], d$recipe_id[f$test])) == 0L,
            max(abs(colMeans(f$xtrain))) < 1e-12,
            max(abs(apply(f$xtrain, 2, sd)-1)) < 1e-12,
            abs(mean(f$ytrain)) < 1e-12, abs(sd(f$ytrain)-1) < 1e-12)
  changed <- d
  changed$strength[f$test] <- changed$strength[f$test] + 1000
  alt <- prepare_fold(changed, fold)
  stopifnot(identical(f$xtrain, alt$xtrain), identical(f$ytrain, alt$ytrain),
            identical(f$preprocessing, alt$preprocessing))
}

# A degenerate mixture must agree with the corresponding analytic forecast.
f <- prepare_fold(d, 1L)
z <- list(model = "unit", fold = 1L, row_id = d$row_id[f$test],
           type = "normal", mu = rep(0,206), sigma = rep(1,206))
a <- forecast_rows(z, d)
z$type <- "normal_mixture"
z$mu <- matrix(0,1,206)
z$sigma <- 1
b <- forecast_rows(z, d)
stopifnot(max(abs(a$nlpd-b$nlpd)) < 1e-12, max(abs(a$pit-b$pit)) < 1e-12,
          max(abs(a$variance-b$variance)) < 1e-10,
          abs(mean(a$crps-b$crps)) < .04)
z$type <- "t"
z$mu <- rep(0,206)
z$sigma <- rep(1,206)
z$df <- 20
a <- forecast_rows(z, d)
z$type <- "t_mixture"
z$mu <- matrix(0,1,206)
z$sigma <- matrix(1,1,206)
z$weights <- 1
b <- forecast_rows(z, d)
stopifnot(max(abs(a$nlpd-b$nlpd)) < 1e-12, max(abs(a$pit-b$pit)) < 1e-12,
          max(abs(a$variance-b$variance)) < 1e-10,
          abs(mean(a$crps-b$crps)) < .04)

nnfile <- fit_path("NN", 1L)
if (file.exists(nnfile)) {
  nn <- readRDS(nnfile)$fit
  j <- nn_jacobian(f$xtrain[1:7, , drop = FALSE], nn)
  theta <- unlist(nn$theta)
  for (pos in c(1,2,3,10,11,40,80)) {
    plus <- minus <- nn
    tp <- tm <- theta
    tp[pos] <- tp[pos]+1e-5
    tm[pos] <- tm[pos]-1e-5
    plus$theta <- split(tp, rep(seq_len(nn$neurons), each = ncol(f$xtrain)+2))
    minus$theta <- split(tm, rep(seq_len(nn$neurons), each = ncol(f$xtrain)+2))
    fd <- (predict(plus,f$xtrain[1:7, , drop = FALSE])-predict(minus,f$xtrain[1:7, , drop = FALSE]))/2e-5
    stopifnot(max(abs(fd+j[,pos])) < 1e-6)
  }
}

if (file.exists(path_input("predictions.csv"))) {
  p <- read.csv(path_input("predictions.csv"))
  primary <- c("OLS","Ridge","Lasso","Horseshoe","BMA","GAM","BART-1000","GP","NN")
  complete_families <- c(primary,intersect("BART",unique(p$model)))
  for (m in complete_families) {
    q <- p[p$model == m, ]
    stopifnot(nrow(q) == 1030L, !anyDuplicated(q$row_id),
              all(q$q025 < q$q050), all(q$q050 < q$q100), all(q$q100 < q$q250),
              all(q$q250 < q$q500), all(q$q500 < q$q750), all(q$q750 < q$q900),
              all(q$q900 < q$q950), all(q$q950 < q$q975))
  }
  stopifnot(all(p$role[p$model == "BART"] == "diagnostic_only"),
            all(p$role[p$model %in% primary] == "primary"))
  if (Sys.getenv("BSLM_INPUT", "reference") == "new") {
    diag <- readRDS(path_result("diagnostics.rds"))
    sample_names <- grep("^(Ridge|Lasso|Horseshoe|BART-1000)-", names(diag), value = TRUE)
    stopifnot(length(sample_names)==20L,
      all(vapply(diag[sample_names], function(x) tail(x,1)[[1]]$pass, logical(1))))
  } else {
    diag <- read.csv("reference-mcmc-diagnostics.csv")
    diag <- diag[grepl("^(Ridge|Lasso|Horseshoe|BART-1000)-", diag$fit), ]
    stopifnot(nrow(diag) == 20L, all(diag$pass), all(diag$max_rhat <= 1.01),
              all(diag$min_bulk >= 400), all(diag$min_tail >= 400))
  }
}
cat("PASS: grouped folds, training-only preprocessing, forecast units, mixture density/CDF/variance,\n")
cat("Monte Carlo CRPS, NN finite-difference Jacobian (when fitted), and row-level interval integrity.\n")
