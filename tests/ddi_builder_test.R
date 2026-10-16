suppressPackageStartupMessages(library(mrgsolve))
arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file_argument), ".."), winslash = "/", mustWork = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
source(file.path(APP_ROOT, "R", "model_library.R"))
source(file.path(APP_ROOT, "R", "engine.R"))
source(file.path(APP_ROOT, "R", "ddi_engine.R"))
source(file.path(APP_ROOT, "R", "ddi_builder.R"))
code <- "$PARAM CL=5, V=30\n$CMT CENT\n$ODE dxdt_CENT = -CL/V*CENT;\n$TABLE double CP=CENT/V;\n$CAPTURE CP"
cache <- new.env()
affected <- ddi_custom_context(code, 1, tempdir(), cache, TRUE)
driver <- ddi_custom_context(code, 2, tempdir(), cache, TRUE)
config <- list(affected = list(dose = 100, interval = 12, infusion = 0), driver = list(dose = 100, interval = 12, infusion = 0),
  start_day = 1, stop_day = 3, followup_days = 1, target = "CL", type = "factor", factor = 0.5, strength = 0.8, c50 = 1)
for (type in c("factor", "inhibition", "induction", "reversible", "hill_inhibition")) {
  config$type <- type
  cpp <- ddi_model_code(affected, config)
  model <- mrgsolve::mcode(paste0("ddi_", type), cpp, soloc = tempdir(), quiet = TRUE)
  control <- as.data.frame(mrgsolve::mrgsim(mrgsolve::ev(model, amt = 100), end = 24, delta = 0.1))
  original <- as.data.frame(mrgsolve::mrgsim(mrgsolve::ev(affected$model, amt = 100), end = 24, delta = 0.1))
  stopifnot(max(abs(control$CP - original$CP)) < 1e-8)
  model <- mrgsolve::param(model, DDI_CP = 2, DDI_ACTIVE = 1)
  modified <- as.data.frame(mrgsolve::mrgsim(mrgsolve::ev(model, amt = 100), end = 24, delta = 0.1))
  multiplier <- ddi_modifier(type, 2, 1, 0, 2, 0.5, 0.8, 1)
  expected <- 100/30 * exp(-5/30 * multiplier * modified$time)
  stopifnot(max(abs(expected - modified$CP)[modified$time > 0]) < 1e-6)
}
for (type in c("tdi", "turnover_induction")) {
  config$type <- type
  model <- mrgsolve::mcode(paste0("ddi_", type), ddi_model_code(affected, config), soloc = tempdir(), quiet = TRUE)
  time <- seq(0, 24, by = 0.1)
  activity <- ddi_modifier(type, rep(2, length(time)), time, 0, 24, 1, config$strength, config$c50)
  data <- data.frame(ID = 1, time = time, evid = 0, cmt = 1, amt = 0, DDI_ACTIVITY = activity)
  data$evid[1] <- 1; data$amt[1] <- 100
  simulation <- as.data.frame(mrgsolve::mrgsim_d(model, data = data, end = 24, delta = 0.1, nocb = FALSE, recsort = 3))
  simulation <- simulation[!duplicated(simulation$time, fromLast = TRUE), ]
  expected <- 100/30 * exp(-5/30 * c(0, cumsum(head(activity, -1) * 0.1)))
  stopifnot(max(abs(expected - simulation$CP)) < 2e-5)
}
multi <- config
multi$type <- "factor"
multi$factor <- 0.5
multi$targets <- list(list(parameter = "CL", fraction = 0.5), list(parameter = "V", fraction = 0.25))
multi_model <- mrgsolve::mcode("ddi_multi", ddi_model_code(affected, multi), soloc = tempdir(), quiet = TRUE)
stopifnot(all(c("DDI_BASE_CL", "DDI_BASE_V", "DDI_FRACTION_1", "DDI_FRACTION_2") %in% model_param_names(multi_model)))
multi_model <- mrgsolve::param(multi_model, DDI_ACTIVE = 1)
multi_simulation <- as.data.frame(mrgsolve::mrgsim(mrgsolve::ev(multi_model, amt = 100), end = 24, delta = 0.1))
expected <- 100/(30 * 0.875) * exp(-(5 * 0.75)/(30 * 0.875) * multi_simulation$time)
stopifnot(max(abs(expected - multi_simulation$CP)[multi_simulation$time > 0]) < 1e-6)
affected_library <- ddi_library_context("tacrolimus_woillard_ddi")
config$target <- "TVCL_TAC"
compiled <- mrgsolve::mcode("ddi_library", ddi_model_code(affected_library, config), soloc = tempdir(), quiet = TRUE)
stopifnot("DDI_BASE_TVCL_TAC" %in% model_param_names(compiled))
invalid <- config; invalid$affected$infusion <- 1
stopifnot(inherits(try(ddi_simulate(invalid, affected_library, driver), silent = TRUE), "try-error"))
config$target <- "CL"
config$type <- "tdi"
result <- ddi_simulate(config, affected, driver)
export <- new.env(parent = globalenv())
pdf(file = NULL)
eval(parse(text = ddi_export_script(config, affected, driver)), envir = export)
dev.off()
stopifnot(isTRUE(all.equal(result$affected_profile, export$result$affected_profile)))
cat("DDI builder: pasted models, C++ exports, analytical agreement, library model and standalone R export passed.\n")
