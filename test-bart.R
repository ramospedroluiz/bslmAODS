source("fit-bart.R")
mock <- function(n) list(mu=matrix(0,n,1))
stopifnot(identical(pending_bart_targets(list(mock(30L),mock(15L)),budgets=c(15L,30L,60L)),c(30L,60L)),
  identical(pending_bart_targets(list(mock(30L),mock(30L)),list(list(draws=30L)),
    budgets=c(15L,30L,60L)),60L),
  identical(pending_bart_targets(list(NULL,mock(15L)),budgets=c(15L,30L,60L)),c(15L,30L,60L)))
f <- prepare_fold(read_data(),1L)
whole <- new_bart_chain(f,14789L,warmup=10L,draws=30L)
set.seed(14789L)
native <- dbarts::bart(f$xtrain,f$ytrain,f$xtest,ntree=1000L,k=2,
  usequants=FALSE,numcut=100L,nchain=1L,nthread=1L,nskip=10L,ndpost=30L,
  keepevery=1L,keeptrainfits=FALSE,keeptrees=FALSE,verbose=FALSE)
stopifnot(max(abs(whole$mu-native$yhat.test))<1e-12,
          max(abs(whole$sigma-tail(native$sigma,30L)))<1e-12)
part <- new_bart_chain(f,14789L,warmup=10L,draws=15L)
part <- extend_bart_chain(part,15L)
stopifnot(identical(whole$mu,part$mu),identical(whole$sigma,part$sigma))
part <- new_bart_chain(f,14789L,warmup=10L,draws=15L)
file <- tempfile(fileext=".rds")
saveRDS(part,file)
part <- readRDS(file)
part <- extend_bart_chain(part,15L)
stopifnot(max(abs(whole$mu-part$mu))<1e-12,max(abs(whole$sigma-part$sigma))<1e-12,
          identical(whole$rng,part$rng))
cat("PASS: native bart output agrees within 1e-12; in-memory continuation is identical;\n")
cat("serialized continuation agrees within 1e-12 with identical RNG state.\n")
unlink(file)
cat("PASS: partial resumes select a common chain length.\n")
