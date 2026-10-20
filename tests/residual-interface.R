library(frpp)
source(file.path("helpers","kernel-api.R"))
fails <- function(expr,text) {
 error <- tryCatch(force(expr),error=identity)
 stopifnot(inherits(error,"error"),grepl(text,conditionMessage(error),fixed=TRUE))
}
primary <- c("dresidual","presidual","qresidual","rresidual")
aliases <- c("density","cdf","quantile","rng")
functions <- list(dexp,pexp,qexp,rexp)
make_law <- function(selected=rep(TRUE,4),names=primary) {
 args <- functions[selected];names(args) <- names[selected]
 do.call(residual_law,c(list(param=c(rate=1)),args))
}
full <- make_law()
stopifnot(identical(full,make_law(names=aliases)),
 identical(full,residual_law(c(rate=1),dexp,pexp,qexp,rexp)),
 identical(residual_law(c(rate=1),dexp,pexp,qexp,rexp,"Exponential"),
   residual_law(param=c(rate=1),dresidual=dexp,presidual=pexp,
     qresidual=qexp,rresidual=rexp,name="Exponential")))
set.seed(763);original_seed <- .Random.seed
for(mask in 0:15) {
 selected <- as.logical(intToBits(mask)[1:4])
 mixed <- aliases;mixed[selected] <- primary[selected]
 stopifnot(identical(full,make_law(names=mixed)))
}
for(i in 1:4) {
 args <- setNames(list(functions[[i]],functions[[i]]),c(primary[i],aliases[i]))
 if(i %in% 2:3) args$rresidual <- rexp
 fails(do.call(residual_law,args),paste0("'",primary[i],"' or '",aliases[i],"'"))
 args[[aliases[i]]] <- NULL
 stopifnot(inherits(do.call(residual_law,args),"frpp_residual_law"))
 args <- setNames(functions,primary);args[aliases[i]] <- list(NULL)
 stopifnot(identical(full,do.call(residual_law,c(list(param=c(rate=1)),args))))
 args[primary[i]] <- list(NULL);args[aliases[i]] <- list(functions[[i]])
 stopifnot(identical(full,do.call(residual_law,c(list(param=c(rate=1)),args))))
}
fails(residual_law(rresidual=1),"functions or NULL")
fails(residual_law(rng=1),"functions or NULL")
stopifnot(identical(original_seed,.Random.seed))
# Neither constructing a law nor validating a model runs callbacks/defaults.
callback <- function(n,param=stop("default must remain lazy")) stop("callback must not run")
law <- residual_law(rresidual=callback)
stopifnot(identical(law$functions$rng,callback),
 identical(environment(law$functions$rng),environment(callback)),
 identical(formals(law$functions$rng),formals(callback)))
kernel <- methods::new("hspec",mu=1,alpha=0,beta=1)
invisible(fspec(kernel,law))
stopifnot(identical(original_seed,.Random.seed))
rm(".Random.seed",envir=.GlobalEnv)
invisible(fspec(kernel,residual_law(rresidual=callback)))
stopifnot(!exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
assign(".Random.seed",original_seed,envir=.GlobalEnv)

fits <- list()
providers <- c("stats",if(requireNamespace("maxLik",quietly=TRUE))"maxLik")
fit_routes <- c("bridge",if(frpp_has_parameter_api)"independent")
for(provider in providers)for(route in fit_routes)
 fits[[paste(provider,route)]] <- ffit(fspec(kernel,full),c(0,.2,.3,.5),
  optimizer=provider,fit_backend=route,lambda_component0=0,fixed=c(FALSE,TRUE,TRUE))

# Exhaust all 16 subsets. Each operation demands only its own callbacks.
for(mask in 0:15) {
 selected <- as.logical(intToBits(mask)[1:4])
 if(!any(selected[c(1,4)])) {
  fails(make_law(selected),"dresidual")
  fails(make_law(selected,names=aliases),"dresidual")
  next
 }
 law <- make_law(selected)
 stopifnot(identical(law,make_law(selected,names=aliases)))
 for(d in 1:2) {
  k <- methods::new("hspec",mu=rep(1,d),alpha=matrix(0,d,d),beta=matrix(1,d,d))
  model <- fspec(k,law);reference <- fspec(k,full)
  for(backend in frpp_backends) {
   set.seed(72);before <- .Random.seed
   if(selected[4]) {
    a <- fsim(model,size=9,lambda_component0=0,backend=backend)
    rng <- .Random.seed;following <- runif(2)
    set.seed(72);b <- fsim(reference,size=9,lambda_component0=0,backend=backend)
    fields <- setdiff(names(a),"model")
    stopifnot(identical(a[fields],b[fields]),identical(rng,.Random.seed),identical(following,runif(2)))
   } else {
    fails(fsim(model,size=9,lambda_component0=0,backend=backend),"rresidual")
    stopifnot(identical(before,.Random.seed))
   }
   args <- list(inter_arrival=c(0,.2,.3,.5),type=c(0,rep(seq_len(d),length.out=3)),
     lambda_component0=0,backend=backend)
   needed <- if(!selected[1])"dresidual" else if(d>1 && !selected[2])"presidual" else NULL
   if(is.null(needed)) {
    stopifnot(identical(do.call(logLik,c(list(object=model),args)),
      do.call(logLik,c(list(object=reference),args))))
   } else fails(do.call(logLik,c(list(object=model),args)),needed)
   needed <- if(!selected[1])"dresidual" else if(!selected[2])"presidual" else if(!selected[3])"qresidual" else NULL
   if(is.null(needed)) {
    stopifnot(identical(do.call(logLik,c(list(object=model,infer=TRUE),args)),
      do.call(logLik,c(list(object=reference,infer=TRUE),args))))
   } else fails(do.call(logLik,c(list(object=model,infer=TRUE),args)),needed)
  }
  if(!selected[1] || (d>1 && !selected[2])) {
   set.seed(81);before <- .Random.seed
   fails(ffit(model,c(0,.2,.3,.5),type=c(0,rep(seq_len(d),length.out=3)),
     optimizer="stats",valid=function(param)stop("validator must not run")),
     if(!selected[1])"dresidual" else "presidual")
   stopifnot(identical(before,.Random.seed))
  }
 }
 if(selected[1])for(provider in providers)for(route in fit_routes) {
  set.seed(82);before <- .Random.seed
  fit <- ffit(fspec(kernel,law),c(0,.2,.3,.5),optimizer=provider,fit_backend=route,
    lambda_component0=0,fixed=c(FALSE,TRUE,TRUE))
  expected <- fits[[paste(provider,route)]]
  stopifnot(identical(coef(fit),coef(expected)),identical(logLik(fit),logLik(expected)),
    identical(vcov(fit),vcov(expected)),identical(fdiagnostics(fit),fdiagnostics(expected)),
    identical(finitial_state(fit),finitial_state(expected)),identical(before,.Random.seed))
 }
}
cat("Residual names/aliases, positional calls, mixed aliases, conflicts, all 16 subsets, partial fits and RNG: OK\n")
