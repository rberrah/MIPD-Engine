suppressPackageStartupMessages(library(mrgsolve))
arguments <- commandArgs(trailingOnly = FALSE)
path <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(path), ".."), winslash = "/")
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
for (name in c("ddi_engine", "pd_engine", "pd_pk", "onco_engine", "workbench_bridge")) source(file.path(APP_ROOT, "R", paste0(name, ".R")))
fails <- function(expression) stopifnot(inherits(try(force(expression), silent = TRUE), "try-error"))
config <- onco_config()
comparison <- onco_compare(config, fraction = 0.5, delay = 7)
before <- function(curve) curve[curve$time < config$decision, c("time", "TUMOR", "ANC", "concentration")]
stopifnot(isTRUE(all.equal(before(comparison$curves$maintain), before(comparison$curves$change), tolerance = 1e-6, check.attributes = FALSE)))
stopifnot(comparison$metrics$final_tumor_ratio[2] > comparison$metrics$final_tumor_ratio[1], comparison$metrics$future_dose_mg[2] < comparison$metrics$future_dose_mg[1])
same <- onco_compare(config, fraction = 1)
stopifnot(identical(same$curves$maintain, same$curves$change))
no_drug <- config; no_drug$dose <- 0; no_drug$history <- config$history[FALSE, ]
control <- onco_simulate(no_drug)
stopifnot(max(abs(control$ANC - config$parameters$ANC0)) < 1e-8,
  max(abs(control$TUMOR / (config$parameters$T0 * exp(config$parameters$KG * control$time)) - 1)) < 1e-6)
for (growth in c("exponential", "logistic", "gompertz")) for (toxicity in c(FALSE, TRUE)) {
  candidate <- config; candidate$growth <- growth; candidate$toxicity <- toxicity; candidate$horizon <- 28; candidate$decision <- 21; candidate$history <- config$history[1, ]
  doses <- onco_schedule(candidate)
  curve <- onco_simulate(candidate)
  model <- mrgsolve::mcode(paste0("onco_", growth, toxicity), onco_model_code(candidate), soloc = tempdir(), quiet = TRUE)
  events <- data.frame(ID = 1, time = doses$time, amt = doses$amount, rate = doses$amount / (doses$infusion / 24), cmt = 1, evid = 1)
  cpp <- as.data.frame(mrgsolve::mrgsim_d(model, data = events, tgrid = curve$time, end = candidate$horizon, recsort = 3))
  cpp <- cpp[!duplicated(cpp$time, fromLast = TRUE), ]
  for (endpoint in c("TUMOR", if (toxicity) "ANC")) {
    expected <- approx(cpp$time, cpp[[endpoint]], xout = curve$time, rule = 2)$y
    stopifnot(max(abs(expected - curve[[endpoint]]) / pmax(1, curve[[endpoint]])) < 2e-4)
  }
}
times <- seq(0, 42, by = 3)
truth <- onco_simulate(config, doses = config$history, times = times)
truth <- truth[match(times, truth$time), ]
data <- rbind(data.frame(time = times, endpoint = "tumor", value = truth$TUMOR), data.frame(time = times, endpoint = "anc", value = truth$ANC))
start <- config; start$parameters$KILL <- 0.2; start$parameters$SLOPE <- 0.12
fit <- onco_fit(start, data)
stopifnot(abs(fit$config$parameters$KILL - config$parameters$KILL) < 1e-3, abs(fit$config$parameters$SLOPE - config$parameters$SLOPE) < 1e-3)
fails(onco_fit(config, transform(data, time = time + 43)))
fails(onco_fit(config, subset(data, endpoint == "tumor"), c("KILL", "SLOPE")))
fails(onco_read_observations("time,endpoint,value\n0,tumor,0"))
fails(onco_config(list(history = data.frame(time = 42, amount = 100, infusion = 1))))
fails(onco_config(list(parameters = list(unsafe_code = 1))))
payload <- list(version = 1, id = "test-123", view = "onco", config = config)
stopifnot(identical(workbench_payload(payload)$config, config))
fails(workbench_payload(utils::modifyList(payload, list(version = 2))))
export <- new.env(parent = globalenv()); pdf(file = NULL)
eval(parse(text = onco_export_script(config, fraction = 0.5, delay = 7)), export)
dev.off()
stopifnot(isTRUE(all.equal(comparison$metrics, export$result$metrics)))
cat("Oncology: growth, equilibrium, history preservation, scenario separation, six C++ variants, parameter recovery, validation and standalone export passed.\n")
