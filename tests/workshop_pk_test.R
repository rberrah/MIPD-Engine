args <- commandArgs(FALSE)
path <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(path), ".."), winslash = "/")
setwd(APP_ROOT)
source("app.R", local = TRUE)
cases <- jsonlite::fromJSON(paste(system2("node", "../scripts/test_workshop_pk.mjs --json", stdout = TRUE), collapse = ""), simplifyVector = FALSE)
original_exposure <- pd_exposure
pd_exposure <- function(time, config) {
  if (is.null(config$test_doses)) return(original_exposure(time, config))
  onco_exposure(time / 24, config$test_doses, list(V = config$v, CL = config$cl * 24))
}
for (case in cases) {
  config <- case$config
  config$parameters <- unlist(config$parameters)
  rows <- do.call(rbind, lapply(case$rows, as.data.frame))
  # Independent R exposure and LSODA integration, not a replay of the JS profile.
  doses <- data.frame(time = seq(0, 24 - 1e-8, by = config$regimen$interval), amount = config$regimen$dose, infusion = config$regimen$infusion)
  time <- sort(unique(c(0, 24, doses$time, doses$time + doses$infusion)))
  config$test_doses <- transform(doses, time = time / 24)
  config$exposure <- "pk"
  config$profile <- data.frame(time = time, concentration = onco_exposure(time / 24, transform(doses, time = time / 24), list(V = config$v, CL = config$cl * 24)))
  result <- pd_simulate(config, unique(rows$time))
  stopifnot(max(abs(result$effect - rows$effect[match(result$time, rows$time)])) < .002)
}
pd_exposure <- original_exposure
workdir <- tempfile("workshop-pk-"); dir.create(workdir); cache <- new.env()
context <- pd_builtin_iv_context(30, 6, workdir, cache)
stopifnot(mrgsolve::param(context$model)$TV_v_L1_CENT == 30, mrgsolve::param(context$model)$TV_cl_L1_CENT == 6)
profile <- pd_pk_profile(context, data.frame(time = 0, amount = 300, infusion = 0), 12)
stopifnot(max(abs(profile$concentration - 10 * exp(-.2 * profile$time))) < 1e-5)
fit <- infection_population_fit(context, workdir, cache)
stopifnot(all(diag(as.matrix(mrgsolve::omat(fit$model))) == .09))
invisible(pd_builtin_iv_context(40, 8, workdir, cache))
stopifnot(length(ls(cache)) == 1)

config <- onco_config(list(toxicity = FALSE))
grid <- list(min = 50, max = 100, step = 50, intervals = c(14, 21))
result <- onco_compare_grid(config, grid)
stopifnot(nrow(result$metrics) == 5, length(result$curves) == 6, !anyNA(result$metrics$targets_met))
same <- result$regimens$scenario[result$regimens$dose == 100 & result$regimens$interval == 21]
stopifnot(identical(result$curves[[same[1]]], result$curves[[same[2]]]))
pair <- onco_grid_selection(result, "candidate_1")
stopifnot(identical(names(pair$curves), c("maintain", "change", "untreated")), nrow(pair$metrics) == 2)
before <- result$curves$maintain$time < config$decision
stopifnot(max(abs(result$curves$maintain$TUMOR[before] - approx(pair$curves$change$time, pair$curves$change$TUMOR, result$curves$maintain$time[before])$y)) < .001)
stopifnot(inherits(try(onco_compare_grid(config, list(min = 1, max = 100, step = 1, intervals = 21)), silent = TRUE), "try-error"))
script <- onco_export_script(config, grid = grid)
env <- new.env(parent = globalenv()); eval(parse(text = script), envir = env)
stopifnot(identical(result$metrics, env$result$metrics))
cat("15 browser PK/PD cases agree with independent R/LSODA; IV PTA context and oncology grid/export passed.\n")
