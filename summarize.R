source("config.R")
need_packages("ggplot2")
p <- read.csv(path_input("predictions.csv"))
timing <- read.csv(path_input("timings.csv"))
models <- c("OLS", "Ridge", "Lasso", "Horseshoe", "BMA", "GAM", "BART-1000", "GP", "NN")
metric_cols <- c(RMSE = "squared_error", MAE = "absolute_error", NLPD = "nlpd", CRPS = "crps",
                 Coverage90 = "cover90", Width90 = "width90")
summarize_one <- function(x) {
  v <- vapply(metric_cols, function(n) mean(x[[n]]), numeric(1))
  v["RMSE"] <- sqrt(v["RMSE"])
  v
}
summary <- do.call(rbind, lapply(split(p, p$model), function(x)
  data.frame(model = x$model[1], n = nrow(x), t(summarize_one(x)))))
summary <- merge(summary, aggregate(cbind(fit_seconds, predictive_seconds) ~ model, timing, sum),
                 all.x = TRUE, sort = FALSE)
summary$role <- ifelse(summary$model %in% models,"primary",
  ifelse(grepl("^BART",summary$model),"diagnostic_only","sensitivity"))
write.csv(summary, path_result("summary-all-forecasts.csv"), row.names = FALSE)
summary <- summary[summary$model %in% models,]
write.csv(summary, path_result("summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(split(p, interaction(p$model, p$fold, drop = TRUE)), function(x)
  data.frame(model = x$model[1], fold = x$fold[1], t(summarize_one(x))))),
  path_result("fold-metrics.csv"), row.names = FALSE)

d <- read_data()
groups <- sort(unique(d$recipe_id))
set.seed(cfg$seed + 777L)
counts <- rmultinom(cfg$bootstrap, length(groups), rep(1/length(groups), length(groups)))
boot <- vector("list", length(models))
names(boot) <- models
ci <- list()
for (model in models) {
  x <- p[p$model == model, ]
  stopifnot(nrow(x) == nrow(d), !anyDuplicated(x$row_id))
  sums <- rowsum(cbind(n = 1, as.matrix(x[, metric_cols])), x$recipe_id, reorder = TRUE)
  stopifnot(identical(as.integer(rownames(sums)), groups))
  b <- t(counts) %*% sums
  vals <- b[, -1, drop = FALSE]/b[, 1]
  colnames(vals) <- names(metric_cols)
  vals[, "RMSE"] <- sqrt(vals[, "RMSE"])
  boot[[model]] <- vals
  ci[[model]] <- data.frame(model = model, metric = names(metric_cols),
    estimate = summarize_one(x), lower = apply(vals, 2, quantile, .025),
    upper = apply(vals, 2, quantile, .975))
}
write.csv(do.call(rbind, ci), path_result("metric-intervals.csv"), row.names = FALSE)
pairs <- combn(models, 2, simplify = FALSE)
paired <- lapply(pairs, function(pair) {
  delta <- boot[[pair[1]]] - boot[[pair[2]]]
  data.frame(first = pair[1], second = pair[2], metric = names(metric_cols),
    difference = summarize_one(p[p$model == pair[1], ])-summarize_one(p[p$model == pair[2], ]),
    lower = apply(delta, 2, quantile, .025), upper = apply(delta, 2, quantile, .975))
})
write.csv(do.call(rbind, paired), path_result("paired-differences.csv"), row.names = FALSE)

coverage <- do.call(rbind, lapply(models, function(model) do.call(rbind, lapply(c(50,80,90,95), function(l)
  data.frame(model = model, nominal = l/100, observed = mean(p[p$model==model,paste0("cover",l)]))))))
write.csv(coverage, path_result("coverage.csv"), row.names = FALSE)
mc <- do.call(rbind, lapply(split(p, p$model), function(x) data.frame(model = x$model[1],
  crps_half_mean_difference = mean(x$crps_half_difference),
  mean_negative_mass = mean(x$negative_mass), maximum_negative_mass = max(x$negative_mass))))
write.csv(mc, path_result("predictive-checks.csv"), row.names = FALSE)
print(summary, row.names = FALSE)

labels <- c(OLS="OLS", Ridge="Ridge", Lasso="Lasso", Horseshoe="Horseshoe", BMA="BMA",
            GAM="GAM (REML)", "BART-1000"="BART", GP="GP (EB)", NN="NN (EB/GN)")
colors <- c(OLS="#707070",Ridge="#827600",Lasso="#91500E",Horseshoe="#AA4499",BMA="#3B7897",
            GAM="#007D5A","BART-1000"="#005B91",GP="#AD4500",NN="#222222")
ltypes <- c(OLS="longdash",Ridge="dashed",Lasso="dotdash",Horseshoe="dotted",BMA="twodash",
             GAM="F282","BART-1000"="solid",GP="88",NN="3313")
pit <- do.call(rbind, lapply(models, function(m) {
  v <- sort(p$pit[p$model == m])
  data.frame(model=m,x=c(0,v,1),y=c(0,seq_along(v)/length(v),1), panel="(b) Probability integral transforms")
}))
cal <- data.frame(model=coverage$model,x=coverage$nominal,y=coverage$observed,
                  panel="(a) Central predictive coverage")
plotdata <- rbind(cal, pit)
plotdata$model <- factor(plotdata$model, levels=models)
cal$model <- factor(cal$model,levels=models)
library(ggplot2)
g <- ggplot(plotdata, aes(x,y,color=model,linetype=model)) +
  geom_abline(slope=1,intercept=0,color="grey75",linewidth=.35) +
  geom_line(data=cal,linewidth=.55) + geom_step(data=pit,linewidth=.55) +
  geom_point(data=cal,size=1.3) + facet_wrap(~panel,nrow=1,scales="free") +
  scale_color_manual(values=colors,labels=labels,name=NULL) +
  scale_linetype_manual(values=ltypes,labels=labels,name=NULL) +
  scale_x_continuous() + scale_y_continuous() +
  labs(x="Nominal probability / PIT value",y="Observed proportion") +
  theme_bw(base_size=10,base_family="sans") +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(linewidth=.2,color="grey93"),
        strip.background=element_blank(),strip.text=element_text(face="bold",size=10),
        legend.position="bottom",legend.key.width=grid::unit(1.3,"cm"),legend.text=element_text(size=9)) +
  guides(color=guide_legend(nrow=3,byrow=TRUE),linetype=guide_legend(nrow=3,byrow=TRUE))
ggsave(path_result("calibration.pdf"),g,width=7.0,height=4.5,device=cairo_pdf)
ggsave(path_result("calibration.png"),g,width=7.0,height=4.5,dpi=220)
cat("SUMMARIES COMPLETE\n")
