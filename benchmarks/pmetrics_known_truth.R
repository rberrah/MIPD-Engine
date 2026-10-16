options(warn = 1)
suppressPackageStartupMessages({
  library(Pmetrics)
  library(jsonlite)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
output_dir <- file.path(root, "tdm-engine", "pmetrics", "benchmarks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

config <- list(
  benchmark = "known-truth parametric versus NPAG",
  generated_on = "2026-10-11",
  population_seed = 20261011L,
  train_subjects = 90L,
  test_subjects = 180L,
  dose_mg = 100,
  volume_L = 30,
  train_times_h = c(1, 4, 8),
  heldout_conditioning_time_h = 1,
  heldout_prediction_time_h = 8,
  residual_model = "Gaussian SD = 0.02 mg/L + 0.10 * prediction",
  npag_runs = list(list(seed = 23L, points = 40L), list(seed = 71L, points = 80L)),
  npag_cycles = 50L,
  note = "Synthetic methodological benchmark; not clinical validation or evidence of universal superiority."
)
write_json(config, file.path(output_dir, "pmetrics-known-truth-config.json"), pretty = TRUE, auto_unbox = TRUE)

dose <- config$dose_mg
volume <- config$volume_L
times <- config$train_times_h
residual_sd <- function(prediction) 0.02 + 0.10 * prediction
prediction <- function(cl, time) dose / volume * exp(-(cl / volume) * time)

sample_clearance <- function(scenario, n) {
  if (scenario == "lognormal") return(rlnorm(n, log(5) - 0.5 * 0.28^2, 0.28))
  if (scenario == "asymmetric") {
    high <- runif(n) < 0.15
    return(ifelse(high, rlnorm(n, log(9) - 0.5 * 0.14^2, 0.14), rlnorm(n, log(4) - 0.5 * 0.2^2, 0.2)))
  }
  high <- runif(n) < 0.5
  ifelse(high, rlnorm(n, log(8) - 0.5 * 0.14^2, 0.14), rlnorm(n, log(3) - 0.5 * 0.14^2, 0.14))
}

simulate_population <- function(scenario, n, set_name, id_offset) {
  cl <- sample_clearance(scenario, n)
  do.call(rbind, lapply(seq_len(n), function(index) {
    pred <- prediction(cl[[index]], times)
    data.frame(
      scenario = scenario, set = set_name, id = id_offset + index, cl_true = cl[[index]],
      time = times, prediction_true = pred,
      concentration = pmax(0.001, stats::rnorm(length(times), pred, residual_sd(pred))),
      stringsAsFactors = FALSE
    )
  }))
}

set.seed(config$population_seed)
scenarios <- c("lognormal", "asymmetric", "bimodal")
synthetic <- do.call(rbind, lapply(seq_along(scenarios), function(index) {
  scenario <- scenarios[[index]]
  rbind(
    simulate_population(scenario, config$train_subjects, "train", index * 10000L),
    simulate_population(scenario, config$test_subjects, "test", index * 10000L + 1000L)
  )
}))
utils::write.csv(synthetic, file.path(output_dir, "pmetrics-known-truth-data.csv"), row.names = FALSE)

log_sum_exp <- function(values) {
  maximum <- max(values)
  maximum + log(sum(exp(values - maximum)))
}

subject_log_likelihood <- function(rows, cl_grid) {
  vapply(cl_grid, function(cl) {
    pred <- prediction(cl, rows$time)
    sum(stats::dnorm(rows$concentration, pred, residual_sd(pred), log = TRUE))
  }, numeric(1))
}

fit_parametric <- function(rows) {
  ids <- split(rows, rows$id)
  z_grid <- seq(log(0.3), log(15), length.out = 181)
  cl_grid <- exp(z_grid)
  likelihood <- vapply(ids, subject_log_likelihood, numeric(length(cl_grid)), cl_grid = cl_grid)
  objective <- function(parameters) {
    prior <- stats::dnorm(z_grid, parameters[[1]], exp(parameters[[2]]), log = TRUE)
    -sum(apply(likelihood, 2, function(column) log_sum_exp(prior + column) - log_sum_exp(prior)))
  }
  optimum <- stats::optim(
    c(log(5), log(0.35)), objective, method = "L-BFGS-B",
    lower = c(log(0.5), log(0.05)), upper = c(log(12), log(2))
  )
  weights <- stats::dnorm(z_grid, optimum$par[[1]], exp(optimum$par[[2]]))
  list(
    cl = cl_grid, prob = weights / sum(weights),
    mu_log_cl = optimum$par[[1]], sd_log_cl = exp(optimum$par[[2]]),
    objective = optimum$value, convergence = optimum$convergence
  )
}

as_pmetrics_data <- function(rows) {
  ids <- unique(rows$id)
  data <- do.call(rbind, lapply(ids, function(id) {
    observations <- rows[rows$id == id, , drop = FALSE]
    rbind(
      data.frame(id = id, time = 0, dose = dose, input = 1L, out = NA_real_),
      data.frame(id = id, time = observations$time, dose = NA_real_, input = NA_integer_, out = observations$concentration)
    )
  }))
  PM_data$new(data = data, quiet = TRUE)
}

npag_model <- PM_model$new(
  pri = list(cl = ab(0.3, 15)),
  eqn = function() { dx[1] <- -(cl / 30) * x[1] + rateiv[1] + b[1] },
  out = function() { y[1] <- x[1] / 30 },
  err = list(proportional(1, c(0.02, 0.10, 0, 0), fixed = TRUE))
)

fit_npag <- function(rows, scenario, run_config) {
  path <- tempfile(paste0("npag-benchmark-", scenario, "-"))
  dir.create(path, recursive = TRUE)
  on.exit(unlink(path, recursive = TRUE, force = TRUE), add = TRUE)
  fit <- suppressMessages(npag_model$fit(
    data = as_pmetrics_data(rows), path = path, run = 1L,
    cycles = config$npag_cycles, prior = "sobol", points = run_config$points,
    seed = run_config$seed, algorithm = "NPAG", report = "none", quiet = TRUE, overwrite = TRUE
  ))
  points <- as.data.frame(fit$final$popPoints)
  points$prob <- points$prob / sum(points$prob)
  list(cl = points$cl, prob = points$prob, objective = fit$final$logLik %||% NA_real_)
}

`%||%` <- function(value, fallback) if (is.null(value) || !length(value)) fallback else value

mixture_quantile <- function(probability, cl, weights, time) {
  means <- prediction(cl, time)
  sds <- residual_sd(means)
  cdf <- function(value) sum(weights * stats::pnorm(value, means, sds)) - probability
  stats::uniroot(cdf, c(0, max(means + 8 * sds)), extendInt = "yes")$root
}

evaluate_distribution <- function(distribution, test_rows) {
  keep <- is.finite(distribution$cl) & distribution$cl > 0 & is.finite(distribution$prob) & distribution$prob > 0
  distribution$cl <- distribution$cl[keep]
  distribution$prob <- distribution$prob[keep] / sum(distribution$prob[keep])
  if (!length(distribution$cl) || any(!is.finite(distribution$prob))) stop("Invalid fitted population distribution.")
  subjects <- split(test_rows, test_rows$id)
  evaluated <- lapply(subjects, function(rows) {
    first <- rows[which.min(abs(rows$time - config$heldout_conditioning_time_h)), ]
    last <- rows[which.min(abs(rows$time - config$heldout_prediction_time_h)), ]
    first_mean <- prediction(distribution$cl, first$time)
    log_weights <- log(distribution$prob) + stats::dnorm(first$concentration, first_mean, residual_sd(first_mean), log = TRUE)
    posterior <- exp(log_weights - log_sum_exp(log_weights))
    future_mean <- prediction(distribution$cl, last$time)
    estimate <- sum(posterior * future_mean)
    predictive_log_likelihood <- log_sum_exp(log(posterior) + stats::dnorm(last$concentration, future_mean, residual_sd(future_mean), log = TRUE))
    lower <- mixture_quantile(0.05, distribution$cl, posterior, last$time)
    upper <- mixture_quantile(0.95, distribution$cl, posterior, last$time)
    data.frame(
      id = last$id, true = last$prediction_true, observed = last$concentration,
      estimate = estimate, predictive_log_likelihood = predictive_log_likelihood,
      covered = last$concentration >= lower && last$concentration <= upper
    )
  })
  values <- do.call(rbind, evaluated)
  error <- values$estimate - values$true
  data.frame(
    mean_predictive_log_likelihood = mean(values$predictive_log_likelihood),
    bias_mg_L = mean(error),
    relative_bias_percent = mean(error / values$true) * 100,
    rmse_mg_L = sqrt(mean(error^2)),
    relative_rmse_percent = sqrt(mean((error / values$true)^2)) * 100,
    mae_mg_L = mean(abs(error)),
    coverage_90_percent = mean(values$covered) * 100
  )
}

mode_recovery <- function(scenario, distribution, method) {
  if (method == "parametric_lognormal") {
    return(if (scenario == "lognormal") "one imposed mode, matching the generating family" else "one imposed mode; separate modes cannot be represented")
  }
  if (scenario == "lognormal") return("support-point count is not a mode count; no extra biology inferred")
  low_mass <- sum(distribution$prob[distribution$cl < 5.5])
  high_mass <- sum(distribution$prob[distribution$cl >= 5.5])
  threshold <- if (scenario == "asymmetric") 0.05 else 0.20
  if (low_mass >= threshold && high_mass >= threshold) "both known clearance regions represented" else "one known clearance region underrepresented"
}

results <- list()
for (scenario in scenarios) {
  train <- synthetic[synthetic$scenario == scenario & synthetic$set == "train", , drop = FALSE]
  test <- synthetic[synthetic$scenario == scenario & synthetic$set == "test", , drop = FALSE]
  parametric <- fit_parametric(train)
  metrics <- evaluate_distribution(parametric, test)
  results[[length(results) + 1L]] <- cbind(
    data.frame(scenario = scenario, method = "parametric_lognormal", run = "MLE", seed = NA_integer_, initial_points = NA_integer_,
      support_points = NA_integer_, distribution_mean_cl = sum(parametric$cl * parametric$prob),
      mode_recovery = mode_recovery(scenario, parametric, "parametric_lognormal"), stringsAsFactors = FALSE),
    metrics
  )
  for (index in seq_along(config$npag_runs)) {
    run_config <- config$npag_runs[[index]]
    npag <- fit_npag(train, scenario, run_config)
    metrics <- evaluate_distribution(npag, test)
    results[[length(results) + 1L]] <- cbind(
      data.frame(scenario = scenario, method = "NPAG", run = paste0("NPAG-", index), seed = run_config$seed,
        initial_points = run_config$points, support_points = length(npag$cl),
        distribution_mean_cl = sum(npag$cl * npag$prob), mode_recovery = mode_recovery(scenario, npag, "NPAG"), stringsAsFactors = FALSE),
      metrics
    )
  }
}

results <- do.call(rbind, results)
utils::write.csv(results, file.path(output_dir, "pmetrics-known-truth-results.csv"), row.names = FALSE)
write_json(results, file.path(output_dir, "pmetrics-known-truth-results.json"), dataframe = "rows", pretty = TRUE, auto_unbox = TRUE, na = "null")

stability <- do.call(rbind, lapply(scenarios, function(scenario) {
  rows <- results[results$scenario == scenario & results$method == "NPAG", , drop = FALSE]
  data.frame(
    scenario = scenario,
    mean_cl_range = diff(range(rows$distribution_mean_cl)),
    rmse_range = diff(range(rows$rmse_mg_L)),
    predictive_log_likelihood_range = diff(range(rows$mean_predictive_log_likelihood))
  )
}))
utils::write.csv(stability, file.path(output_dir, "pmetrics-known-truth-stability.csv"), row.names = FALSE)

print(results[, c("scenario", "method", "run", "support_points", "relative_bias_percent", "relative_rmse_percent", "coverage_90_percent")])
cat("Known-truth Pmetrics benchmark complete. Results are synthetic and conditional on this design.\n")
