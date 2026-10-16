suppressPackageStartupMessages(library(mrgsolve))
arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file_argument), ".."), winslash = "/", mustWork = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
source(file.path(APP_ROOT, "R", "ddi_engine.R"))
source(file.path(APP_ROOT, "R", "pd_engine.R"))
fails <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))
base <- pd_config(list(type = "emax", parameters = list(E0 = 10, EMAX = 50, EC50 = 2)))
stopifnot(abs(pd_simulate(base)$effect[[1]] - (10 + 50 * 10 / 12)) < 1e-8)
for (type in PD_TYPES) for (delay in c(FALSE, TRUE)) {
  config <- pd_config(list(type = type, delay = delay, parameters = list(E0 = 10, EMAX = 0.8)))
  result <- pd_simulate(config)
  stopifnot(nrow(result) == 501, all(is.finite(as.matrix(result))))
  config$c0 <- 0
  stopifnot(max(abs(pd_simulate(config)$effect - 10)) < 1e-6)
  model <- mrgsolve::mcode(paste0("pd_", type, "_", delay), pd_model_code(config), soloc = tempdir(), quiet = TRUE)
  compiled <- as.data.frame(mrgsolve::mrgsim(model, end = 24, delta = 0.25))
  stopifnot(max(abs(compiled$PRED - 10)) < 1e-6)
  # A constant concentration isolates equation/export agreement from grid interpolation.
  config$c0 <- 5; config$kel <- 0
  model <- mrgsolve::param(model, CP = 5)
  compiled <- as.data.frame(mrgsolve::mrgsim(model, end = 24, delta = 0.25))
  expected <- pd_simulate(config, compiled$time)
  stopifnot(max(abs(compiled$PRED - expected$effect)) < 1e-5)
}
for (type in c("linear", "emax", "hill", "inhibit_in")) {
  truth <- pd_config(list(type = type, parameters = list(E0 = 10, SLOPE = 2, EMAX = 0.8, EC50 = 2, KOUT = 0.15)))
  data <- pd_simulate(truth, seq(0, 24, length.out = 30))[, c("time", "effect")]
  initial <- truth
  estimate <- if (type == "linear") c("E0", "SLOPE") else c("E0", "EMAX", "EC50")
  initial$parameters[estimate] <- initial$parameters[estimate] * 1.2
  fit <- pd_fit(initial, data, "time", estimate)
  stopifnot(fit$rmse < 1e-5, max(abs(fit$config$parameters[estimate] - truth$parameters[estimate])) < 0.01)
}
data <- data.frame(concentration = seq(0, 20, length.out = 30))
data$effect <- pd_predict(base, data, "concentration")
initial <- base; initial$parameters["EC50"] <- 3
stopifnot(pd_fit(initial, data, "concentration", "EC50")$rmse < 1e-7)
fails(pd_predict(modifyList(base, list(delay = TRUE)), data, "concentration"))
fails(pd_config(list(type = "inhibit_in", parameters = list(EMAX = 2))))
fails(pd_config(list(exposure = "table", profile = data.frame(time = c(0, 0), concentration = c(1, 2)))))
fails(pd_config(list(exposure = "table", profile = data.frame(time = c(1, 24), concentration = c(1, 2)))))
fails(pd_read_table("time,effect\n2,not-a-number", c("time", "effect")))
fails(pd_fit(base, data[1:3, ], "concentration"))
pulse <- pd_config(list(type = "linear", delay = TRUE, horizon = 1000, exposure = "table",
  parameters = list(E0 = 10, SLOPE = 2, KE0 = 0.5),
  profile = data.frame(time = c(0, 10, 10.001, 10.002, 1000), concentration = c(0, 0, 5, 0, 0))))
expected <- sum(vapply(c(10, 10.001), function(from) stats::integrate(function(time) 0.5 * pd_exposure(time, pulse) * exp(-0.5 * (11 - time)), from, from + 0.001)$value, numeric(1)))
stopifnot(abs(pd_simulate(pulse, c(0, 11, 1000))$driver[[2]] - expected) < 1e-7)
export <- new.env(parent = globalenv())
pdf(file = NULL)
eval(parse(text = pd_export_script(base)), envir = export)
dev.off()
stopifnot(isTRUE(all.equal(export$result, pd_simulate(base))))
cat("PD: 14 model/delay exports compiled, equation agreement, parameter recovery, validation and R export passed.\n")
