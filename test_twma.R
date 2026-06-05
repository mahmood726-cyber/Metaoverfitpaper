#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## Smoke test for twma_analysis.R
##
## Sources the TWMA function definitions (without running the full, slow
## analysis) and asserts a few basic correctness properties with stopifnot().
## Exits 0 on success, non-zero (via stop) on failure.
##
## Run:  Rscript test_twma.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
    stopifnot(requireNamespace("metafor", quietly = TRUE))
    stopifnot(requireNamespace("metadat", quietly = TRUE))
    library(metafor)
    library(metadat)
})

## --- Source ONLY the function definitions from twma_analysis.R -------------
## We parse the script and evaluate just the top-level `<- function(...)`
## assignments, so sourcing is fast and side-effect free (no datasets built,
## no bootstrap run).
stopifnot(file.exists("twma_analysis.R"))
exprs <- parse("twma_analysis.R")
fn_env <- new.env()
for (e in exprs) {
    if (is.call(e) &&
        as.character(e[[1]])[1] %in% c("<-", "=") &&
        length(e) == 3 &&
        is.call(e[[3]]) &&
        identical(as.character(e[[3]][[1]])[1], "function")) {
        eval(e, envir = fn_env)
    }
}

cat("Loaded functions:", paste(ls(fn_env), collapse = ", "), "\n")
stopifnot("calculate_t_squared_robust" %in% ls(fn_env))
stopifnot("t2_optimism_corrected"      %in% ls(fn_env))

calculate_t_squared_robust <- get("calculate_t_squared_robust", fn_env)
t2_optimism_corrected      <- get("t2_optimism_corrected",      fn_env)

## --- Property 1: a known metadat dataset loads -----------------------------
data(dat.bcg)
stopifnot(nrow(dat.bcg) == 13L)
es <- escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg,
             data = dat.bcg)
y <- es$yi; v <- es$vi
X <- cbind(intercept = 1,
           latitude  = as.numeric(scale(dat.bcg$ablat)),
           year      = as.numeric(scale(dat.bcg$year)))
stopifnot(length(y) == nrow(X))

## --- Property 2: a fitted rma object has finite tau^2 >= 0 ------------------
fit <- rma(yi = y, vi = v, method = "REML")
stopifnot(is.finite(fit$tau2))
stopifnot(fit$tau2 >= 0)

## --- Property 3: T^2 from the core function is in [0, 1] and finite --------
res <- suppressWarnings(calculate_t_squared_robust(y, v, X))
stopifnot(is.finite(res$t_squared))
stopifnot(res$t_squared >= 0, res$t_squared <= 1)
stopifnot(is.finite(res$tau2_null), res$tau2_null >= 0)

## --- Property 4: optimism correction returns non-negative weights/quantities
oc <- suppressWarnings(t2_optimism_corrected(y, v, X))
stopifnot(is.finite(oc$tau2_cv), oc$tau2_cv >= 0)        # transportability var >= 0
stopifnot(oc$t2_corrected >= 0, oc$t2_corrected <= 1)
stopifnot(oc$t2_original  >= 0)
## corrected T^2 should not exceed apparent T^2 (optimism is non-negative here)
stopifnot(oc$optimism >= -1e-8)

cat("\nALL SMOKE TESTS PASSED\n")
cat(sprintf("  BCG: tau^2(null)=%.4f  T^2=%.3f  T~^2=%.3f  optimism=%.3f\n",
            res$tau2_null, res$t_squared, oc$t2_corrected, oc$optimism))
