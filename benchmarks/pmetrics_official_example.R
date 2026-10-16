options(warn = 1)
suppressPackageStartupMessages({
  library(Pmetrics)
  library(jsonlite)
})

full_run <- "--full" %in% commandArgs(trailingOnly = TRUE)
cycles <- if (full_run) 1000L else 75L
initial_points <- if (full_run) 2028L else 128L
tolerance_percent <- if (full_run) 2 else 8

data_path <- system.file("Examples", "src", "bimodal_ke.csv", package = "Pmetrics")
if (!nzchar(data_path) || !file.exists(data_path)) {
  stop("The official bimodal_ke.csv example is unavailable in this Pmetrics installation.")
}

data <- PM_data$new(data_path, quiet = TRUE)
example_rows <- utils::read.csv(data_path, na.strings = ".", check.names = FALSE)
id_column <- grep("ID$", names(example_rows), ignore.case = TRUE, value = TRUE)[1]
if (is.na(id_column)) stop("The official example has no identifiable subject ID column.")
common <- list(
  pri = list(ke = ab(0.001, 3.0), v = ab(25.0, 250.0)),
  out = function() { y[1] <- x[1] / v },
  err = list(additive(1, c(0.0, 0.5, 0.0, 0.0)))
)

analytical_model <- do.call(PM_model$new, c(common, list(
  eqn = function() { one_comp_iv }
)))
ode_model <- do.call(PM_model$new, c(common, list(
  eqn = function() { dx[1] <- -ke * x[1] + rateiv[1] + b[1] },
  solver = "TSIT45"
)))

fit_model <- function(model, label) {
  path <- tempfile(paste0("pmetrics-official-", label, "-"))
  dir.create(path, recursive = TRUE)
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)
  elapsed <- system.time(capture.output({
    fit <- model$fit(
      data = data, path = path, run = 1L, algorithm = "NPAG",
      cycles = cycles, prior = "sobol", points = initial_points,
      seed = 22L, idelta = 0.01, report = "none", quiet = TRUE,
      overwrite = TRUE
    )
  }))[["elapsed"]]

  points <- as.data.frame(fit$final$popPoints)
  required <- c("ke", "v", "prob")
  if (!all(required %in% names(points))) stop(label, " output is missing ke, v or prob.")
  if (any(!is.finite(as.matrix(points[, required])))) stop(label, " output contains non-finite values.")
  if (any(points$prob <= 0) || abs(sum(points$prob) - 1) > 1e-8) stop(label, " weights are invalid.")
  if (any(points$ke < 0.001 | points$ke > 3.0) || any(points$v < 25 | points$v > 250)) {
    stop(label, " support points fall outside the official parameter bounds.")
  }

  mean_table <- as.data.frame(fit$final$popMean)
  means <- setNames(as.numeric(mean_table[1, c("ke", "v")]), c("ke", "v"))
  list(
    label = label,
    status = as.character(fit$cycle$data$status %||% "unknown"),
    elapsed_seconds = unname(elapsed),
    support_points = nrow(points),
    weights_sum = sum(points$prob),
    mean_ke_per_h = unname(means[["ke"]]),
    mean_v_L = unname(means[["v"]]),
    minimum_ke_per_h = min(points$ke),
    maximum_ke_per_h = max(points$ke),
    minimum_v_L = min(points$v),
    maximum_v_L = max(points$v)
  )
}

`%||%` <- function(value, fallback) if (is.null(value) || !length(value)) fallback else value

analytical <- fit_model(analytical_model, "official-analytical-shortcut")
ode <- fit_model(ode_model, "official-ode-reference")
relative_difference <- c(
  ke = 100 * abs(analytical$mean_ke_per_h - ode$mean_ke_per_h) / max(abs(ode$mean_ke_per_h), .Machine$double.eps),
  v = 100 * abs(analytical$mean_v_L - ode$mean_v_L) / max(abs(ode$mean_v_L), .Machine$double.eps)
)
passed <- all(relative_difference <= tolerance_percent)

description <- utils::packageDescription("Pmetrics")
result <- list(
  status = if (passed) "passed" else "failed",
  generated_on = "2026-10-11",
  mode = if (full_run) "official-full" else "ci-smoke",
  provenance = list(
    package = "Pmetrics",
    package_version = as.character(utils::packageVersion("Pmetrics")),
    package_remote_sha = description$RemoteSha %||% NA_character_,
    source_file = "Examples/src/bimodal_ke.csv",
    source_md5 = unname(tools::md5sum(data_path)),
    source_url = "https://github.com/LAPKB/Pmetrics/tree/v3.2.6/inst/Examples",
    data_are_loaded_from_installed_package = TRUE
  ),
  design = list(
    subjects = length(unique(example_rows[[id_column]])),
    algorithm = "NPAG",
    seed = 22L,
    cycles = cycles,
    initial_points = initial_points,
    parameter_bounds = list(ke_per_h = c(0.001, 3.0), v_L = c(25, 250)),
    error_model = "additive(1, c(0.0, 0.5, 0.0, 0.0))"
  ),
  analytical = analytical,
  ode_reference = ode,
  comparison = list(
    relative_difference_percent = as.list(relative_difference),
    tolerance_percent = tolerance_percent,
    passed = passed,
    interpretation = "Numerical equivalence check between the two official Pmetrics formulations; not clinical validation."
  )
)

output <- Sys.getenv("PMETRICS_BENCHMARK_OUTPUT", unset = "")
if (!nzchar(output)) {
  output <- if (full_run) {
    file.path(getwd(), "tdm-engine", "pmetrics", "benchmarks", "pmetrics-official-example.json")
  } else {
    tempfile("pmetrics-official-smoke-", fileext = ".json")
  }
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write_json(result, output, pretty = TRUE, auto_unbox = TRUE, na = "null")
if (!full_run && !nzchar(Sys.getenv("PMETRICS_BENCHMARK_OUTPUT", unset = ""))) {
  unlink(output, force = TRUE)
}

cat(sprintf(
  "Official Pmetrics example: %s; max relative mean difference %.3f%% (tolerance %.1f%%).\n",
  result$status, max(relative_difference), tolerance_percent
))
if (!passed) stop("The analytical and ODE official example fits differ beyond the documented tolerance.")
