# Optional independent replay, after tools/validate-legacy-fit.R.
# Rscript tools/replay-fit-validation.R <output-directory>
# Requires the emhawkes development shared-kernel API.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply the fit-validation output directory.")
out <- args[1L]
library(frpp)
if (!all(c("as_exp_kernel", "exp_kernel_replay") %in%
    getNamespaceExports("emhawkes"))) stop("Independent replay requires the emhawkes development shared-kernel API.")
x <- readRDS(file.path(out,"inputs.rds"))
results <- readRDS(file.path(out,"all-fits.rds"))
set.seed(451); rng <- .Random.seed
rows <- lapply(names(results),function(name) {
 fit <- results[[name]]
 model <- if (inherits(fit,"ffit")) fitted_fspec(fit) else emhawkes::fitted_hspec(fit)
 initial <- if (inherits(fit,"ffit")) finitial_state(fit) else emhawkes::fit_initial_state(fit)
 if (!inherits(model,"frpp_spec")) model <- as_fspec(model)
 bridge <- do.call(stats::logLik,c(list(object=model,inter_arrival=x$path$inter_arrival,
   type=x$path$type,infer=TRUE),initial))
 state <- do.call(stats::logLik,c(list(object=model,inter_arrival=x$path$inter_arrival,
   type=x$path$type,infer=TRUE,backend="state"),initial))
 if (name == names(results)[1L]) {
   law <- as_legacy_hspec(model)
   driver_total <- residual_total <- event_total <- 0
   for (n in seq_len(length(x$path$type) - 1L) + 1L) {
     terms <- numeric(ncol(state$integrated_driver))
     for (k in seq_along(terms)) {
       terms[k] <- if (k == x$path$type[n])
         law@dresidual(state$integrated_driver[n, k], log = TRUE) else
         law@presidual(state$integrated_driver[n, k], lower.tail = FALSE, log.p = TRUE)
       residual_total <- residual_total + terms[k]
     }
     value <- log(state$lambda[n, x$path$type[n]])
     driver_total <- driver_total + value
     event_total <- event_total + value + sum(terms)
   }
   orders <- data.frame(bridge = bridge$loglikelihood,
     bridge_sum_order = driver_total + residual_total, state = state$loglikelihood,
     state_sum_order = event_total)
   stopifnot(identical(unname(driver_total + residual_total), unname(bridge$loglikelihood)),
     identical(unname(event_total), unname(state$loglikelihood)))
   write.csv(orders, file.path(out, "summation-order.csv"), row.names = FALSE)
 }
 stopifnot(isTRUE(all.equal(bridge,state[names(bridge)],tolerance=1e-10)),
   isTRUE(all.equal(as.numeric(stats::logLik(fit)),bridge$loglikelihood,tolerance=1e-10)),
   identical(rng,.Random.seed))
 data.frame(fit=name,optimized=as.numeric(stats::logLik(fit)),bridge=bridge$loglikelihood,
   state=state$loglikelihood,state_difference=state$loglikelihood-bridge$loglikelihood,
   maximum_state_output_difference=max(vapply(names(bridge),function(nm)
      if(length(bridge[[nm]]))max(abs(bridge[[nm]]-state[[nm]]))else 0,numeric(1))))
})
rows <- do.call(rbind,rows)
write.csv(rows,file.path(out,"saved-and-state-replay.csv"),row.names=FALSE)
print(rows,digits=15)
cat("PASS: saved fits, restored FRPP laws, independent shared-kernel inference and unchanged RNG.\n")