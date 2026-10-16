options(warn = 1)
suppressPackageStartupMessages({
  library(Pmetrics)
  library(jsonlite)
})

output <- file.path(getwd(), "tdm-engine", "pmetrics", "benchmarks", "pmetrics-npod-smoke.json")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

model <- PM_model$new(
  pri = list(cl = ab(0.5, 12)),
  eqn = function() { dx[1] <- -(cl / 30) * x[1] + rateiv[1] + b[1] },
  out = function() { y[1] <- x[1] / 30 },
  err = list(proportional(1, c(0.02, 0.1, 0, 0), fixed = TRUE))
)

set.seed(20261011)
ids <- seq_len(24)
clearance <- ifelse(ids <= 12, 3, 8) * exp(stats::rnorm(length(ids), 0, 0.12))
times <- c(1, 4, 8)
rows <- do.call(rbind, lapply(ids, function(id) {
  pred <- 100 / 30 * exp(-(clearance[[id]] / 30) * times)
  rbind(
    data.frame(id = id, time = 0, dose = 100, input = 1L, out = NA_real_),
    data.frame(id = id, time = times, dose = NA_real_, input = NA_integer_,
      out = pmax(0.001, stats::rnorm(length(times), pred, 0.02 + 0.1 * pred)))
  )
}))
data <- PM_data$new(data = rows, quiet = TRUE)

run_algorithm <- function(algorithm, seed) {
  path <- tempfile(paste0("pmetrics-", tolower(algorithm), "-"))
  dir.create(path, recursive = TRUE)
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)
  started <- proc.time()[[3]]
  log_output <- capture.output({
    fit <- model$fit(
      data = data, path = path, run = 1L, algorithm = algorithm,
      cycles = 35, prior = "sobol", points = 30, seed = seed,
      report = "none", quiet = TRUE, overwrite = TRUE
    )
  })
  objective_lines <- grep("Objective function =", log_output, value = TRUE, fixed = TRUE)
  final_objective <- if (length(objective_lines)) {
    suppressWarnings(as.numeric(sub(".*Objective function = ([0-9.+-]+).*", "\\1", tail(objective_lines, 1L))))
  } else NA_real_
  points <- as.data.frame(fit$final$popPoints)
  final_values <- if (is.environment(fit$final)) as.list.environment(fit$final) else as.list(fit$final)
  final_scalars <- final_values[vapply(final_values, function(value) is.numeric(value) && length(value) == 1L && is.finite(value), logical(1))]
  list(
    algorithm = algorithm,
    seed = seed,
    elapsed_seconds = unname(proc.time()[[3]] - started),
    final_objective = final_objective,
    support_points = nrow(points),
    mean_clearance = sum(points$cl * points$prob) / sum(points$prob),
    weights_sum = sum(points$prob),
    finite_positive_weights = all(is.finite(points$prob) & points$prob > 0),
    final_numeric_fields = final_scalars
  )
}

result <- tryCatch({
  runs <- list(
    run_algorithm("NPAG", 23L), run_algorithm("NPAG", 71L),
    run_algorithm("NPOD", 23L), run_algorithm("NPOD", 71L)
  )
  list(
    status = "available",
    pmetrics_version = as.character(utils::packageVersion("Pmetrics")),
    generated_on = "2026-10-11",
    note = "Experimental smoke comparison only; algorithms use different searches and this is not a superiority benchmark.",
    runs = runs
  )
}, error = function(error) {
  list(
    status = "unavailable",
    pmetrics_version = as.character(utils::packageVersion("Pmetrics")),
    generated_on = "2026-10-11",
    error = conditionMessage(error),
    note = "Non-blocking experimental check; no public NPOD control should be exposed from this result."
  )
})

write_json(result, output, pretty = TRUE, auto_unbox = TRUE)
cat("Pmetrics NPOD experimental status:", result$status, "\n")
