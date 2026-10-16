trap_auc <- function(time, concentration) {
  keep <- is.finite(time) & is.finite(concentration)
  time <- time[keep]
  concentration <- concentration[keep]
  if (length(time) < 2) return(NA_real_)
  order_index <- order(time)
  time <- time[order_index]
  concentration <- concentration[order_index]
  sum(diff(time) * (head(concentration, -1) + tail(concentration, -1)) / 2)
}

carry_covariates <- function(event_times, covariate_history, covariates) {
  history_names <- if (is.null(covariate_history)) character() else setdiff(names(covariate_history), "time")
  names_to_apply <- unique(c(names(covariates), history_names))
  if (!length(names_to_apply)) return(data.frame(row.names = seq_along(event_times)))

  output <- list()
  for (name in names_to_apply) {
    fallback <- suppressWarnings(as.numeric(covariates[[name]] %||% NA_real_))
    if (length(fallback) != 1L || !is.finite(fallback)) fallback <- NA_real_
    history <- data.frame()
    if (!is.null(covariate_history) && all(c("time", name) %in% names(covariate_history))) {
      history <- data.frame(
        time = suppressWarnings(as.numeric(covariate_history$time)),
        value = suppressWarnings(as.numeric(covariate_history[[name]]))
      )
      history <- history[is.finite(history$time) & is.finite(history$value), , drop = FALSE]
      history <- history[order(history$time), , drop = FALSE]
    }

    if (!nrow(history)) {
      output[[name]] <- rep(fallback, length(event_times))
      next
    }
    positions <- findInterval(event_times, history$time)
    positions[positions == 0L] <- 1L
    output[[name]] <- history$value[positions]
  }
  as.data.frame(output, check.names = FALSE)
}

parse_dosing_intervals <- function(selected = numeric(), custom = "") {
  selected_tokens <- as.character(unlist(selected, use.names = FALSE))
  custom_text <- paste(as.character(custom %||% ""), collapse = " ")
  custom_tokens <- if (nzchar(trimws(custom_text))) {
    strsplit(trimws(custom_text), "[,;[:space:]]+", perl = TRUE)[[1]]
  } else {
    character()
  }
  tokens <- c(selected_tokens[nzchar(trimws(selected_tokens))], custom_tokens)
  if (!length(tokens)) stop("Select or enter at least one dosing interval.")
  values <- suppressWarnings(as.numeric(tokens))
  if (any(!is.finite(values) | values <= 0)) {
    stop("Dosing intervals must be positive numbers separated by spaces, commas, or semicolons.")
  }
  sort(unique(values))
}

normalize_steady_state_doses <- function(doses) {
  if (!nrow(doses)) return(doses)
  dose_ss <- if ("ss" %in% names(doses)) suppressWarnings(as.integer(doses$ss)) else rep(0L, nrow(doses))
  doses$time[dose_ss == 1L] <- 0
  doses
}

expand_steady_state_doses <- function(doses, end_time) {
  doses <- normalize_steady_state_doses(doses)
  if (!nrow(doses) || !"ss" %in% names(doses) || !is.finite(end_time)) return(doses)

  rows <- lapply(seq_len(nrow(doses)), function(index) {
    dose <- doses[index, , drop = FALSE]
    steady_state <- suppressWarnings(as.integer(dose$ss[[1]])) == 1L
    interval <- suppressWarnings(as.numeric(dose$interval[[1]]))
    if (!steady_state || !is.finite(interval) || interval <= 0) return(dose)

    repetitions <- max(0L, floor((end_time - dose$time[[1]] + 1e-8) / interval))
    expanded <- dose[rep(1L, repetitions + 1L), , drop = FALSE]
    expanded$time <- dose$time[[1]] + seq.int(0L, repetitions) * interval
    expanded$count <- 1L
    expanded$ss <- c(1L, rep(0L, repetitions))
    expanded
  })
  output <- do.call(rbind, rows)
  output <- output[order(output$time), , drop = FALSE]
  rownames(output) <- NULL
  output
}

LEGO_STEADY_STATE_WARMUP_DOSES <- 50L

is_split_lego_model <- function(model) {
  "LEGO_INPUT" %in% model@cmtL
}

dose_event_rows <- function(doses, adm_cmt, end_time, split_lego = FALSE) {
  doses <- expand_steady_state_doses(doses, end_time)
  if (isTRUE(split_lego) && nrow(doses)) {
    rows <- lapply(seq_len(nrow(doses)), function(index) {
      dose <- doses[index, , drop = FALSE]
      if (as.integer(dose$ss[[1]] %||% 0L) != 1L) return(dose)
      interval <- as.numeric(dose$interval[[1]])
      warmup <- dose[rep(1L, LEGO_STEADY_STATE_WARMUP_DOSES + 1L), , drop = FALSE]
      warmup$time <- dose$time[[1]] + seq.int(-LEGO_STEADY_STATE_WARMUP_DOSES, 0L) * interval
      warmup$count <- 1L
      warmup$ss <- 0L
      warmup
    })
    doses <- do.call(rbind, rows)
    doses <- doses[order(doses$time), , drop = FALSE]
  }
  dose_ss <- if ("ss" %in% names(doses)) as.integer(doses$ss) else rep(0L, nrow(doses))
  output <- data.frame(
    ID = 1,
    time = doses$time,
    evid = 1,
    cmt = adm_cmt,
    amt = doses$amount,
    rate = ifelse(doses$infusion > 0, doses$amount / doses$infusion, 0),
    ii = doses$interval,
    addl = ifelse(dose_ss == 1L, 0, pmax(0, doses$count - 1)),
    ss = dose_ss,
    stringsAsFactors = FALSE
  )
  administration_names <- setdiff(names(doses), c("time", "amount", "interval", "count", "infusion", "ss", "status", "time_uncertainty"))
  for (name in administration_names) output[[name]] <- doses[[name]]
  output
}

build_map_data <- function(doses, observations, adm_cmt, obs_cmt, covariates, covariate_history = NULL, split_lego = FALSE) {
  if (!nrow(doses)) stop("At least one administered dose is required.")
  doses <- normalize_steady_state_doses(doses)
  dose_ss <- if ("ss" %in% names(doses)) as.integer(doses$ss) else rep(0L, nrow(doses))
  if (any(!dose_ss %in% c(0L, 1L)) || any(dose_ss == 1L & doses$interval <= 0)) {
    stop("Steady-state doses require ss = 1 and a positive interval.")
  }

  dose_horizon <- if (nrow(observations)) max(observations$time, na.rm = TRUE) else max(doses$time, na.rm = TRUE)
  dose_rows <- dose_event_rows(doses, adm_cmt, dose_horizon, split_lego = split_lego)
  event_columns <- c("ID", "time", "evid", "cmt", "amt", "rate", "ii", "addl", "ss")
  administration_names <- setdiff(names(dose_rows), event_columns)
  administration_history <- if (length(administration_names)) dose_rows[, c("time", administration_names), drop = FALSE] else NULL
  dose_rows <- dose_rows[, event_columns, drop = FALSE]
  dose_rows$DV <- NA_real_
  dose_rows$mdv <- 1

  observation_rows <- data.frame(
    ID = rep(1, nrow(observations)),
    time = observations$time,
    evid = rep(0, nrow(observations)),
    cmt = rep(obs_cmt, nrow(observations)),
    amt = rep(0, nrow(observations)),
    rate = rep(0, nrow(observations)),
    ii = rep(0, nrow(observations)),
    addl = rep(0, nrow(observations)),
    ss = rep(0L, nrow(observations)),
    DV = observations$concentration,
    mdv = rep(0, nrow(observations)),
    stringsAsFactors = FALSE
  )

  data <- rbind(dose_rows, observation_rows)
  data <- data[order(data$time, -data$evid), , drop = FALSE]
  time_covariates <- carry_covariates(data$time, covariate_history, covariates)
  for (name in names(time_covariates)) data[[name]] <- time_covariates[[name]]
  if (length(administration_names)) {
    administration_covariates <- carry_covariates(data$time, administration_history, as.list(administration_history[1, administration_names, drop = FALSE]))
    for (name in administration_names) data[[name]] <- administration_covariates[[name]]
  }
  rownames(data) <- NULL
  data
}

validate_model_route_set <- function(specifications) {
  routes <- unique(vapply(specifications, function(item) as.character(item$route %||% ""), character(1)))
  routes <- routes[nzchar(routes)]
  if (length(routes) > 1L) stop("Model averaging is only allowed for models using the same administration route.")
  modes <- unique(vapply(specifications, function(item) as.character(item$mode %||% ""), character(1)))
  modes <- modes[nzchar(modes)]
  if (length(modes) > 1L) stop("Model averaging is only allowed for models using the same administration mode.")
  invisible(if (length(routes)) routes[[1]] else "")
}

resolve_administration_cmt <- function(model, contract, specification) {
  name <- as.character(specification$adm_cmt_name %||% "")
  if (!nzchar(name)) return(contract$adm_cmt)
  position <- match(name, model@cmtL)
  if (!is.finite(position)) stop("Administration compartment ", name, " is absent from model ", specification$label, ".")
  position
}

apply_residual_error_setting <- function(model, mode = "model", fixed_cv = 1) {
  if (identical(mode, "model")) return(model)
  if (!identical(mode, "fixed_cv")) stop("Unknown residual error setting: ", mode, ".")

  fixed_cv <- suppressWarnings(as.numeric(fixed_cv))
  if (length(fixed_cv) != 1L || !is.finite(fixed_cv) || fixed_cv <= 0 || fixed_cv > 100) {
    stop("The fixed proportional residual error must be between 0 and 100%.")
  }

  sigma <- as.matrix(mrgsolve::smat(model))
  outcome_count <- max(1L, ceiling(nrow(sigma) / 2L))
  cv <- fixed_cv / 100
  uses_log_error <- isTRUE(getFromNamespace("log_transformation", "mapbayr")(model))
  sigma_pair <- if (uses_log_error) c(0, log1p(cv^2)) else c(cv^2, 0)
  diagonal <- rep(sigma_pair, outcome_count)
  mrgsolve::smat(model, diag(diagonal, nrow = length(diagonal)))
}

fit_one_model <- function(
  specification,
  doses,
  observations,
  covariates,
  allow_custom,
  covariate_history = NULL,
  custom_soloc = NULL,
  custom_cache = NULL,
  residual_error_mode = "model",
  fixed_residual_cv = 1
) {
  if (identical(tolower(as.character(specification$engine %||% "mapbayr")), "pmetrics")) {
    return(fit_one_pmetrics_model(
      specification = specification,
      doses = doses,
      observations = observations,
      covariates = covariates,
      covariate_history = covariate_history,
      custom_soloc = custom_soloc,
      custom_cache = custom_cache,
      residual_error_mode = residual_error_mode,
      fixed_residual_cv = fixed_residual_cv
    ))
  }
  model <- compile_model(
    model_id = specification$id,
    custom_code = specification$code,
    allow_custom = allow_custom,
    custom_soloc = custom_soloc,
    custom_cache = custom_cache
  )
  model <- apply_residual_error_setting(model, residual_error_mode, fixed_residual_cv)
  accepted_covariates <- model_param_names(model)
  model_covariates <- covariates[names(covariates) %in% accepted_covariates]
  current_covariates <- model_covariates
  if (!is.null(covariate_history)) {
    history_names <- intersect(setdiff(names(covariate_history), "time"), accepted_covariates)
    covariate_history <- covariate_history[, c("time", history_names), drop = FALSE]
    history_order <- order(suppressWarnings(as.numeric(covariate_history$time)))
    for (name in history_names) {
      values <- suppressWarnings(as.numeric(covariate_history[[name]][history_order]))
      values <- values[is.finite(values)]
      if (length(values)) current_covariates[[name]] <- tail(values, 1)
    }
  }
  model <- safe_param(model, model_covariates)
  contract <- validate_model_contract(model)
  if (!contract$ok) stop(paste(contract$errors, collapse = " "))
  contract$adm_cmt <- resolve_administration_cmt(model, contract, specification)
  contract$route <- specification$route %||% "Unspecified"
  contract$mode <- specification$mode %||% ""
  split_lego <- is_split_lego_model(model)

  data <- build_map_data(
    doses = doses,
    observations = observations,
    adm_cmt = contract$adm_cmt,
    obs_cmt = contract$obs_cmt,
    covariates = model_covariates,
    covariate_history = covariate_history,
    split_lego = split_lego
  )

  estimate <- NULL
  if (nrow(observations)) {
    estimate <- mapbayr::mapbayest(
      x = model,
      data = data,
      hessian = stats::optimHess,
      verbose = FALSE,
      progress = FALSE
    )
  }

  list(
    id = specification$id %||% "custom",
    label = specification$label,
    model = model,
    estimate = estimate,
    data = data,
    contract = contract,
    split_lego = split_lego,
    current_covariates = current_covariates,
    source_doses = doses,
    source_observations = observations,
    source_covariates = model_covariates,
    source_covariate_history = covariate_history,
    residual_error = list(
      mode = residual_error_mode,
      fixed_cv = if (identical(residual_error_mode, "fixed_cv")) fixed_residual_cv else NA_real_,
      scale = if (identical(residual_error_mode, "fixed_cv") &&
        isTRUE(getFromNamespace("log_transformation", "mapbayr")(model))) "log" else "standard"
    )
  )
}

SCENARIO_MAINTAIN <- "Poursuite de la derni\u00e8re posologie"
SCENARIO_RECOMMENDED <- "Application de la recommandation"

fit_model_set <- function(
  specifications,
  doses,
  observations,
  covariates,
  allow_custom,
  covariate_history = NULL,
  custom_soloc = NULL,
  custom_cache = NULL,
  residual_error_mode = "model",
  fixed_residual_cv = 1
) {
  validate_model_route_set(specifications)
  doses <- normalize_steady_state_doses(doses)
  fits <- lapply(specifications, function(specification) {
    tryCatch(
      fit_one_model(
        specification,
        doses,
        observations,
        covariates,
        allow_custom,
        covariate_history,
        custom_soloc,
        custom_cache,
        residual_error_mode,
        fixed_residual_cv
      ),
      error = function(error) structure(list(
        id = specification$id %||% "custom",
        label = specification$label,
        message = conditionMessage(error)
      ), class = "tdm_fit_error")
    )
  })
  names(fits) <- vapply(specifications, function(item) item$id %||% "custom", character(1))
  fits
}

successful_fits <- function(fits) {
  Filter(function(item) !inherits(item, "tdm_fit_error"), fits)
}

fit_numeric_diagnostics <- function(fit, boundary_fraction = 0.01, condition_limit = 1e8) {
  observations <- fit$source_observations %||% data.frame()
  output <- list(
    id = fit$id %||% "unknown",
    label = fit$label %||% fit$id %||% "Model",
    engine = fit$engine %||% "mapbayr",
    observations = nrow(observations),
    status = "fitted",
    converged = NA,
    sparse = FALSE,
    eta_count = as.integer(fit$contract$n_eta %||% 0L),
    at_boundary = FALSE,
    uncertainty_available = FALSE,
    uncertainty_stable = NA,
    covariance_condition = NA_real_,
    effective_support = NA_real_
  )

  if (inherits(fit, "tdm_fit_error")) {
    output$status <- "error"
    return(output)
  }

  if (identical(output$engine, "pmetrics")) {
    output$status <- if (isTRUE(fit$pmetrics$mapped)) "fitted" else "population"
    points <- fit$posterior_points %||% data.frame()
    probabilities <- suppressWarnings(as.numeric(points$prob %||% numeric()))
    probabilities <- probabilities[is.finite(probabilities) & probabilities > 0]
    if (length(probabilities)) {
      probabilities <- probabilities / sum(probabilities)
      output$effective_support <- 1 / sum(probabilities^2)
      output$uncertainty_available <- length(probabilities) > 1L
    }
    return(output)
  }

  estimate <- fit$estimate
  if (is.null(estimate)) {
    output$status <- "population"
    return(output)
  }

  selected <- as.integer(estimate$arg.optim$select_eta %||% seq_len(output$eta_count))
  output$eta_count <- length(selected)
  output$sparse <- output$observations < output$eta_count

  convergence <- suppressWarnings(as.numeric(estimate$opt.value$convergence %||% NA_real_))
  if (length(convergence) && any(is.finite(convergence))) {
    output$converged <- all(convergence[is.finite(convergence)] == 0)
  }

  eta <- suppressWarnings(as.numeric(estimate$final_eta[[1]] %||% numeric()))
  lower <- suppressWarnings(as.numeric(estimate$arg.optim$lower %||% numeric()))
  upper <- suppressWarnings(as.numeric(estimate$arg.optim$upper %||% numeric()))
  if (length(selected) && length(eta) >= max(selected) && length(lower) == length(selected) && length(upper) == length(selected)) {
    width <- upper - lower
    margin <- pmin(eta[selected] - lower, upper - eta[selected]) / width
    output$at_boundary <- any(is.finite(margin) & margin <= boundary_fraction)
  }

  covariance <- tryCatch(as.matrix(estimate$covariance[[1]]), error = function(error) matrix(NA_real_))
  output$uncertainty_available <- nrow(covariance) > 0L && nrow(covariance) == ncol(covariance) && all(is.finite(covariance))
  if (isTRUE(output$uncertainty_available)) {
    eigenvalues <- tryCatch(eigen((covariance + t(covariance)) / 2, symmetric = TRUE, only.values = TRUE)$values, error = function(error) numeric())
    if (length(eigenvalues) && all(is.finite(eigenvalues)) && min(eigenvalues) > 0) {
      output$covariance_condition <- max(eigenvalues) / min(eigenvalues)
      output$uncertainty_stable <- output$covariance_condition <= condition_limit
    } else {
      output$uncertainty_stable <- FALSE
    }
  }
  output
}

compute_model_weights <- function(fits, scheme = "AIC") {
  valid <- successful_fits(fits)
  if (!length(valid)) return(numeric())
  estimates <- lapply(valid, `[[`, "estimate")
  if (length(valid) == 1 || any(vapply(estimates, is.null, logical(1)))) {
    return(stats::setNames(rep(1 / length(valid), length(valid)), names(valid)))
  }

  weights <- tryCatch(
    mapbayr::compute_weights(estlist = estimates, scheme = scheme),
    error = function(error) NULL
  )
  if (is.null(weights)) {
    return(stats::setNames(rep(1 / length(valid), length(valid)), names(valid)))
  }
  values <- as.numeric(weights[1, ])
  stats::setNames(values / sum(values), colnames(weights))
}

individual_model <- function(fit) {
  if (identical(fit$engine %||% "mapbayr", "pmetrics")) {
    model <- safe_param(fit$model, as.list(fit$posterior_parameters %||% numeric()))
  } else {
    model <- if (is.null(fit$estimate)) fit$model else mapbayr::use_posterior(fit$estimate)
  }
  if (length(fit$ml_eta_override %||% numeric())) {
    model <- safe_param(model, as.list(fit$ml_eta_override))
  }
  model <- safe_param(model, fit$current_covariates %||% list())
  mrgsolve::zero_re(model)
}

simulate_regimen <- function(fit, dose, interval, infusion, horizon = 24, delta = 0.1) {
  model <- individual_model(fit)
  extra_doses <- max(0, floor((horizon - 1e-8) / interval))
  warmup_doses <- if (isTRUE(fit$split_lego)) LEGO_STEADY_STATE_WARMUP_DOSES else 0L
  warmup_time <- warmup_doses * interval
  event <- mrgsolve::ev(
    amt = dose,
    ii = interval,
    cmt = fit$contract$adm_cmt,
    ss = if (isTRUE(fit$split_lego)) 0 else 1,
    addl = warmup_doses + extra_doses,
    rate = if (infusion > 0) dose / infusion else 0
  )
  simulation <- model |>
    mrgsolve::ev(event) |>
    mrgsolve::mrgsim(start = warmup_time, end = warmup_time + horizon, delta = delta, recsort = 3) |>
    as.data.frame()
  simulation$time <- simulation$time - warmup_time
  simulation <- simulation[simulation$time >= -1e-8 & simulation$time <= horizon + 1e-8, , drop = FALSE]
  column <- pick_concentration_column(simulation)
  data.frame(
    time = simulation$time,
    concentration = pmax(0, simulation[[column]]),
    model = fit$id,
    stringsAsFactors = FALSE
  )
}

average_profiles <- function(profiles, weights) {
  if (!length(profiles)) stop("No simulation profile is available.")
  grid <- sort(unique(unlist(lapply(profiles, `[[`, "time"))))
  averaged <- rep(0, length(grid))
  used_weight <- 0
  for (name in names(profiles)) {
    weight <- weights[[name]] %||% 0
    if (!is.finite(weight) || weight <= 0) next
    profile <- profiles[[name]]
    averaged <- averaged + stats::approx(
      profile$time,
      profile$concentration,
      xout = grid,
      rule = 2,
      ties = "ordered"
    )$y * weight
    used_weight <- used_weight + weight
  }
  if (used_weight <= 0) stop("No valid model weight is available.")
  data.frame(time = grid, concentration = averaged / used_weight)
}

time_above_threshold_pct <- function(time, concentration, threshold) {
  threshold <- suppressWarnings(as.numeric(threshold))
  if (length(threshold) != 1L || !is.finite(threshold)) return(NA_real_)
  keep <- is.finite(time) & is.finite(concentration)
  time <- as.numeric(time[keep])
  concentration <- as.numeric(concentration[keep])
  if (length(time) < 2L) return(NA_real_)
  order_index <- order(time)
  time <- time[order_index]
  concentration <- concentration[order_index]
  duration <- diff(time)
  if (sum(duration) <= 0) return(NA_real_)
  left <- head(concentration, -1L)
  right <- tail(concentration, -1L)
  above <- ifelse(
    left > threshold & right > threshold,
    duration,
    ifelse(
      left <= threshold & right <= threshold,
      0,
      duration * ifelse(left > threshold, left - threshold, right - threshold) / abs(right - left)
    )
  )
  100 * sum(above) / sum(duration)
}

normalized_range_distance <- function(value, low, high) {
  span <- max(high - low, .Machine$double.eps)
  ifelse(value < low, (low - value) / span, ifelse(value > high, (value - high) / span, 0))
}

score_target_component <- function(value, operator = "range", low, high = low) {
  operator <- as.character(operator %||% "range")[[1]]
  scale <- max(abs(c(low, high)), 1, na.rm = TRUE)
  if (identical(operator, "above")) {
    distance <- ifelse(value >= low, 0, (low - value) / scale)
  } else if (identical(operator, "below")) {
    distance <- ifelse(value <= high, 0, (value - high) / scale)
  } else if (identical(operator, "range")) {
    distance <- normalized_range_distance(value, low, high)
  } else {
    stop("Unknown target operator: ", operator, ".")
  }
  list(distance = distance, in_target = distance == 0)
}

score_target_metrics <- function(data, metric, target_low, target_high, target_options = list()) {
  if (identical(metric, "CminCmax")) {
    cmin_operator <- as.character(target_options$cmin_operator %||% "range")[[1]]
    cmax_operator <- as.character(target_options$cmax_operator %||% "range")[[1]]
    cmax_low <- suppressWarnings(as.numeric(target_options$cmax_low %||% NA_real_))
    cmax_high <- suppressWarnings(as.numeric(target_options$cmax_high %||% NA_real_))
    if (!cmin_operator %in% c("range", "above", "below") ||
        !cmax_operator %in% c("range", "above", "below") ||
        !is.finite(target_low) || !is.finite(target_high) ||
        !is.finite(cmax_low) || !is.finite(cmax_high) ||
        (identical(cmin_operator, "range") && target_high <= target_low) ||
        (identical(cmax_operator, "range") && cmax_high <= cmax_low)) {
      stop("Combined Cmin/Cmax target constraints are invalid.")
    }
    cmin_score <- score_target_component(data$cmin, cmin_operator, target_low, target_high)
    cmax_score <- score_target_component(data$cmax, cmax_operator, cmax_low, cmax_high)
    distance <- pmax(cmin_score$distance, cmax_score$distance)
    return(list(
      value = distance,
      in_target = cmin_score$in_target & cmax_score$in_target,
      distance = distance,
      center_distance = distance
    ))
  }

  selected <- switch(
    metric,
    AUC24 = data$auc24,
    Cmin = data$cmin,
    Cmax = data$cmax,
    TimeAbove = data$time_above_pct,
    stop("Unknown target metric: ", metric, ".")
  )
  center <- (target_low + target_high) / 2
  list(
    value = selected,
    in_target = selected >= target_low & selected <= target_high,
    distance = normalized_range_distance(selected, target_low, target_high),
    center_distance = abs(selected - center) / max(abs(center), .Machine$double.eps)
  )
}

profile_metrics <- function(profile, interval = 24, concentration_threshold = NA_real_) {
  auc_profile <- profile[profile$time <= 24 + 1e-8, , drop = FALSE]
  interval_profile <- profile[profile$time <= interval + 1e-8, , drop = FALSE]
  data.frame(
    auc24 = trap_auc(auc_profile$time, auc_profile$concentration),
    cmin = min(interval_profile$concentration, na.rm = TRUE),
    cmax = max(interval_profile$concentration, na.rm = TRUE),
    time_above_pct = time_above_threshold_pct(
      interval_profile$time,
      interval_profile$concentration,
      concentration_threshold
    )
  )
}

simulate_averaged_regimen <- function(
  fits,
  weights,
  dose,
  interval,
  infusion,
  horizon = NULL,
  delta = 0.1,
  concentration_threshold = NA_real_
) {
  valid <- successful_fits(fits)
  horizon <- max(24, interval, suppressWarnings(as.numeric(horizon %||% 0)), na.rm = TRUE)
  profiles <- lapply(
    valid,
    simulate_regimen,
    dose = dose,
    interval = interval,
    infusion = infusion,
    horizon = horizon,
    delta = delta
  )
  names(profiles) <- names(valid)
  averaged <- average_profiles(profiles, weights)
  list(
    profile = averaged,
    metrics = profile_metrics(averaged, interval, concentration_threshold),
    per_model = profiles
  )
}

simulate_known_history <- function(fit, end_time, delta = 0.05) {
  model <- individual_model(fit)
  event_columns <- c("ID", "time", "evid", "cmt", "amt", "rate", "ii", "addl", "ss")
  history <- dose_event_rows(fit$source_doses, fit$contract$adm_cmt, max(0, end_time), split_lego = isTRUE(fit$split_lego))
  time_covariates <- carry_covariates(
    history$time,
    fit$source_covariate_history %||% NULL,
    fit$source_covariates %||% list()
  )
  covariate_names <- intersect(names(time_covariates), model_param_names(model))
  columns <- c(event_columns, covariate_names)
  for (name in names(time_covariates)) history[[name]] <- time_covariates[[name]]
  for (name in setdiff(columns, names(history))) history[[name]] <- 0
  history <- history[, columns, drop = FALSE]
  simulation <- mrgsolve::mrgsim_d(
    model,
    data = history,
    start = min(0, history$time, na.rm = TRUE),
    end = max(0, end_time),
    delta = delta,
    recsort = 3
  ) |>
    as.data.frame()
  simulation <- simulation[simulation$time >= -1e-8, , drop = FALSE]
  column <- pick_concentration_column(simulation)
  data.frame(
    time = simulation$time,
    concentration = pmax(0, simulation[[column]]),
    model = fit$id,
    stringsAsFactors = FALSE
  )
}

historical_exposure <- function(fits, weights, doses, observations = data.frame(), delta = 0.05) {
  dose_ss <- if ("ss" %in% names(doses)) as.integer(doses$ss) else rep(0L, nrow(doses))
  last_administration <- doses$time + ifelse(dose_ss == 1L, 0, doses$interval * pmax(0, doses$count - 1))
  observation_end <- if (nrow(observations)) max(observations$time, na.rm = TRUE) else -Inf
  decision_time <- max(max(last_administration, na.rm = TRUE), observation_end)
  context <- regimen_context(doses, observations)
  current_c0_time <- context$future_start
  window_start <- max(0, decision_time - 24)
  valid <- successful_fits(fits)
  profiles <- lapply(valid, simulate_known_history, end_time = current_c0_time, delta = delta)
  names(profiles) <- names(valid)
  average <- average_profiles(profiles, weights)
  rows <- average$time >= window_start - 1e-8 & average$time <= decision_time + 1e-8
  current <- if (nrow(average) < 2L) {
    average$concentration[[nrow(average)]]
  } else {
    stats::approx(
      average$time,
      average$concentration,
      xout = decision_time,
      rule = 2,
      ties = "ordered"
    )$y[[1]]
  }
  current_c0 <- if (nrow(average) < 2L) {
    average$concentration[[nrow(average)]]
  } else {
    stats::approx(
      average$time,
      average$concentration,
      xout = current_c0_time,
      rule = 2,
      ties = "ordered"
    )$y[[1]]
  }
  auc <- if (sum(rows) < 2L) 0 else trap_auc(average$time[rows], average$concentration[rows])
  list(
    auc24 = auc,
    concentration = current,
    c0 = current_c0,
    c0_time = current_c0_time,
    decision_time = decision_time,
    window_start = window_start,
    coverage_hours = decision_time - window_start,
    profile = average
  )
}

current_regimen_exposure <- function(fits, weights, doses, observations = data.frame()) {
  if (!nrow(doses)) stop("At least one administered dose is required.")
  dose_ss <- if ("ss" %in% names(doses)) as.integer(doses$ss) else rep(0L, nrow(doses))
  last_administration <- doses$time + ifelse(dose_ss == 1L, 0, doses$interval * pmax(0, doses$count - 1))
  regimen <- doses[which.max(last_administration), , drop = FALSE]
  steady_state <- "ss" %in% names(regimen) && isTRUE(as.integer(regimen$ss[[1]]) == 1L)
  recorded_interval <- as.numeric(regimen$interval[[1]])
  single_dose <- !steady_state && (!is.finite(recorded_interval) || recorded_interval <= 0)
  interval <- if (single_dose) 24 else recorded_interval
  dose <- as.numeric(regimen$amount[[1]])
  infusion <- as.numeric(regimen$infusion[[1]])
  valid <- successful_fits(fits)
  if (!length(valid)) stop("No fitted model is available for current exposure.")

  auc_profiles <- lapply(
    valid,
    simulate_regimen,
    dose = dose,
    interval = interval,
    infusion = infusion,
    horizon = 24,
    delta = 0.05
  )
  names(auc_profiles) <- names(valid)
  average_auc_profile <- average_profiles(auc_profiles, weights)
  auc_rows <- average_auc_profile$time >= 0 & average_auc_profile$time <= 24
  auc24 <- trap_auc(
    average_auc_profile$time[auc_rows],
    average_auc_profile$concentration[auc_rows]
  )

  trough_profiles <- lapply(
    valid,
    simulate_regimen,
    dose = dose,
    interval = interval,
    infusion = infusion,
    horizon = interval,
    delta = min(0.05, interval / 100)
  )
  names(trough_profiles) <- names(valid)
  average_trough_profile <- average_profiles(trough_profiles, weights)
  c0 <- stats::approx(
    average_trough_profile$time,
    average_trough_profile$concentration,
    xout = interval,
    rule = 2,
    ties = "ordered"
  )$y[[1]]

  history <- historical_exposure(fits, weights, doses, observations)

  list(
    auc24 = history$auc24,
    c0 = history$c0,
    historical_auc24 = history$auc24,
    historical_concentration = history$concentration,
    historical_c0 = history$c0,
    historical_c0_time = history$c0_time,
    historical_window_start = history$window_start,
    historical_decision_time = history$decision_time,
    historical_coverage_hours = history$coverage_hours,
    steady_state_auc24 = auc24,
    steady_state_c0 = c0,
    dose = dose,
    interval = interval,
    recorded_interval = recorded_interval,
    infusion = infusion,
    steady_state = steady_state,
    single_dose = single_dose
  )
}

recommend_regimens <- function(
  fits,
  weights,
  dose_min,
  dose_max,
  dose_step,
  intervals,
  infusion,
  metric,
  target_low,
  target_high,
  target_options = list(),
  delta = 0.1
) {
  doses <- seq(dose_min, dose_max, by = dose_step)
  scenarios <- expand.grid(
    dose = doses,
    interval = intervals,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  scenarios$infusion <- infusion
  if (nrow(scenarios) > 500) stop("The dose grid is too large; use a wider dose step.")

  metrics <- lapply(seq_len(nrow(scenarios)), function(index) {
    result <- simulate_averaged_regimen(
      fits,
      weights,
      dose = scenarios$dose[[index]],
      interval = scenarios$interval[[index]],
      infusion = infusion,
      delta = delta,
      concentration_threshold = target_options$concentration_threshold %||% NA_real_
    )$metrics
    result
  })
  metrics <- do.call(rbind, metrics)
  output <- cbind(scenarios, metrics)
  scores <- score_target_metrics(output, metric, target_low, target_high, target_options)
  output$value <- scores$value
  output$in_target <- scores$in_target
  output$distance <- scores$distance
  output$center_distance <- scores$center_distance
  output <- output[order(output$distance, output$center_distance, output$dose), , drop = FALSE]
  rownames(output) <- NULL
  output
}

regimen_context <- function(doses, observations = data.frame()) {
  doses <- normalize_steady_state_doses(doses)
  dose_ss <- if ("ss" %in% names(doses)) as.integer(doses$ss) else rep(0L, nrow(doses))
  last_administration <- doses$time + ifelse(dose_ss == 1L, 0, doses$interval * pmax(0, doses$count - 1))
  regimen <- doses[which.max(last_administration), , drop = FALSE]
  observation_end <- if (nrow(observations)) max(observations$time, na.rm = TRUE) else -Inf
  decision_time <- max(max(last_administration, na.rm = TRUE), observation_end)
  interval <- as.numeric(regimen$interval[[1]])
  if (!is.finite(interval) || interval <= 0) interval <- 24
  origin <- as.numeric(regimen$time[[1]])
  next_index <- floor((decision_time - origin) / interval) + 1
  future_start <- origin + max(1, next_index) * interval
  if (future_start <= decision_time + 1e-8) future_start <- future_start + interval
  list(regimen = regimen, decision_time = decision_time, future_start = future_start, interval = interval)
}

simulate_projected_history <- function(fit, dose, interval, infusion, future_start, dose_count, end_time, delta = 0.1) {
  model <- individual_model(fit)
  event_columns <- c("ID", "time", "evid", "cmt", "amt", "rate", "ii", "addl", "ss")
  covariate_names <- intersect(names(fit$current_covariates %||% list()), model_param_names(model))
  columns <- c(event_columns, covariate_names)

  history <- fit$data[fit$data$evid == 1, , drop = FALSE]
  for (name in setdiff(columns, names(history))) history[[name]] <- 0
  history <- history[, columns, drop = FALSE]

  future <- data.frame(
    ID = 1,
    time = future_start,
    evid = 1,
    cmt = fit$contract$adm_cmt,
    amt = dose,
    rate = if (infusion > 0) dose / infusion else 0,
    ii = interval,
    addl = max(0L, as.integer(dose_count) - 1L),
    ss = 0L,
    stringsAsFactors = FALSE
  )
  for (name in covariate_names) future[[name]] <- as.numeric(fit$current_covariates[[name]])
  future <- future[, columns, drop = FALSE]
  events <- rbind(history, future)
  events <- events[order(events$time), , drop = FALSE]

  simulation <- mrgsolve::mrgsim_d(
    model,
    data = events,
    start = min(0, events$time, na.rm = TRUE),
    end = end_time,
    delta = delta,
    recsort = 3
  ) |>
    as.data.frame()
  simulation <- simulation[simulation$time >= -1e-8, , drop = FALSE]
  column <- pick_concentration_column(simulation)
  data.frame(
    time = simulation$time,
    concentration = pmax(0, simulation[[column]]),
    model = fit$id,
    stringsAsFactors = FALSE
  )
}

compare_future_regimens <- function(
  fits,
  weights,
  doses,
  observations,
  recommended,
  additional_doses = 6,
  delta = 0.1
) {
  context <- regimen_context(doses, observations)
  current <- context$regimen
  current_interval <- context$interval
  current_infusion <- as.numeric(current$infusion[[1]])
  scenarios <- list()
  scenarios[[SCENARIO_MAINTAIN]] <- list(
    dose = as.numeric(current$amount[[1]]),
    interval = current_interval,
    infusion = current_infusion
  )
  scenarios[[SCENARIO_RECOMMENDED]] <- list(
    dose = as.numeric(recommended$dose[[1]]),
    interval = as.numeric(recommended$interval[[1]]),
    infusion = as.numeric(recommended$infusion[[1]] %||% current_infusion)
  )
  maximum_interval <- max(vapply(scenarios, `[[`, numeric(1), "interval"))
  end_time <- context$future_start + max(1L, as.integer(additional_doses)) * maximum_interval
  valid <- successful_fits(fits)

  output <- lapply(names(scenarios), function(label) {
    scenario <- scenarios[[label]]
    profiles <- lapply(valid, simulate_projected_history,
      dose = scenario$dose,
      interval = scenario$interval,
      infusion = scenario$infusion,
      future_start = context$future_start,
      dose_count = additional_doses,
      end_time = end_time,
      delta = delta
    )
    names(profiles) <- names(valid)
    average <- average_profiles(profiles, weights)
    average$scenario <- label
    average
  })

  list(
    profiles = do.call(rbind, output),
    decision_time = context$decision_time,
    future_start = context$future_start,
    end_time = end_time,
    additional_doses = as.integer(additional_doses),
    scenarios = scenarios
  )
}

draw_multivariate_normal <- function(count, center, covariance) {
  center <- as.numeric(center)
  covariance <- as.matrix(covariance)
  if (!length(center)) return(matrix(numeric(), nrow = count, ncol = 0))
  covariance <- (covariance + t(covariance)) / 2
  decomposition <- eigen(covariance, symmetric = TRUE)
  values <- pmax(decomposition$values, 0)
  transform <- decomposition$vectors %*% diag(sqrt(values), nrow = length(values))
  draws <- matrix(stats::rnorm(count * length(center)), nrow = count) %*% t(transform)
  sweep(draws, 2, center, "+")
}

timing_eta_offsets <- function(fit, count, seed = 419) {
  if (is.null(fit$estimate) || count <= 0) return(matrix(numeric(), nrow = 0, ncol = 0))
  doses <- fit$source_doses
  observations <- fit$source_observations
  dose_uncertainty <- suppressWarnings(as.numeric(doses$time_uncertainty %||% rep(0, nrow(doses))))
  observation_uncertainty <- suppressWarnings(as.numeric(observations$time_uncertainty %||% rep(0, nrow(observations))))
  dose_uncertainty[!is.finite(dose_uncertainty)] <- 0
  observation_uncertainty[!is.finite(observation_uncertainty)] <- 0
  if (!any(c(dose_uncertainty, observation_uncertainty) > 0)) {
    return(matrix(numeric(), nrow = 0, ncol = 0))
  }

  center_source <- fit$estimate$final_eta[[1]]
  center <- as.numeric(center_source)
  names(center) <- names(center_source)
  set.seed(seed)
  rows <- lapply(seq_len(min(20L, as.integer(count))), function(index) {
    jittered_doses <- doses
    jittered_observations <- observations
    jittered_doses$time <- pmax(0, doses$time + stats::runif(nrow(doses), -dose_uncertainty, dose_uncertainty))
    if (nrow(observations)) {
      shifts <- stats::runif(nrow(observations), -observation_uncertainty, observation_uncertainty)
      jittered_observations$time <- pmax(0, observations$time + shifts)
    }
    data <- build_map_data(
      doses = jittered_doses,
      observations = jittered_observations,
      adm_cmt = fit$contract$adm_cmt,
      obs_cmt = fit$contract$obs_cmt,
      covariates = fit$source_covariates,
      covariate_history = fit$source_covariate_history,
      split_lego = isTRUE(fit$split_lego)
    )
    eta <- tryCatch(
      mapbayr::mapbayest(
        x = fit$model,
        data = data,
        hessian = FALSE,
        output = "eta",
        verbose = FALSE,
        progress = FALSE
      ),
      error = function(error) NULL
    )
    if (is.null(eta)) return(NULL)
    as.numeric(eta[1, ]) - center
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) return(matrix(numeric(), nrow = 0, ncol = 0))
  output <- do.call(rbind, rows)
  colnames(output) <- names(center)
  output
}

posterior_eta_draws <- function(fit, count, include_posterior = TRUE, include_timing = FALSE, timing_refits = 20L, seed = 381) {
  if (identical(fit$engine %||% "mapbayr", "pmetrics")) {
    points <- fit$posterior_points
    parameter_names <- names(fit$posterior_parameters %||% numeric())
    if (!is.data.frame(points) || !nrow(points) || !length(parameter_names)) {
      return(list(draws = NULL, available = FALSE, mode = "pmetrics_unavailable", timing_available = FALSE))
    }
    set.seed(seed)
    selected <- if (isTRUE(include_posterior)) {
      sample(seq_len(nrow(points)), count, replace = TRUE, prob = points$prob)
    } else {
      rep(NA_integer_, count)
    }
    draws <- if (isTRUE(include_posterior)) {
      as.matrix(points[selected, parameter_names, drop = FALSE])
    } else {
      matrix(rep(as.numeric(fit$posterior_parameters), each = count), nrow = count,
        dimnames = list(NULL, parameter_names))
    }
    colnames(draws) <- parameter_names
    return(list(
      draws = draws,
      available = isTRUE(fit$pmetrics$mapped) && nrow(points) > 1L,
      mode = if (include_posterior) "pmetrics_posterior_support" else "pmetrics_posterior_mean",
      timing_available = FALSE
    ))
  }
  if (is.null(fit$estimate)) {
    return(list(draws = NULL, available = FALSE, mode = if (include_posterior) "population_iiv" else "population_fixed", timing_available = FALSE))
  }
  center <- fit$ml_eta_override %||% fit$estimate$final_eta[[1]]
  covariance <- fit$estimate$covariance[[1]]
  covariance_valid <- is.matrix(covariance) && all(dim(covariance) == length(center)) && all(is.finite(covariance))
  set.seed(seed)
  draws <- if (isTRUE(include_posterior) && covariance_valid) {
    draw_multivariate_normal(count, center, covariance)
  } else {
    matrix(rep(as.numeric(center), each = count), nrow = count)
  }
  colnames(draws) <- names(center)

  timing_offsets <- if (isTRUE(include_timing)) timing_eta_offsets(fit, min(count, timing_refits), seed = seed + 97L) else matrix(numeric(), 0, 0)
  if (nrow(timing_offsets)) {
    selected <- sample(seq_len(nrow(timing_offsets)), count, replace = TRUE)
    draws <- draws + timing_offsets[selected, , drop = FALSE]
  }
  list(
    draws = draws,
    available = covariance_valid,
    mode = if (include_posterior && covariance_valid) "posterior_map" else "map_fixed",
    timing_available = nrow(timing_offsets) > 0
  )
}

simulate_model_distribution <- function(
  fit,
  dose,
  interval,
  infusion,
  replicates,
  delta,
  include_posterior,
  include_residual,
  include_timing = FALSE,
  timing_refits = 20L,
  concentration_threshold = NA_real_,
  seed = 381,
  return_profiles = FALSE,
  hours_per_model_time = 1,
  concentration_column = NULL,
  concentration_scale = 1
) {
  if (!is.finite(hours_per_model_time) || hours_per_model_time <= 0 ||
      !is.finite(concentration_scale) || concentration_scale <= 0) stop("Invalid profile unit conversion.")
  if (return_profiles && include_residual) stop("Exposure profiles require residual error disabled.")
  eta_draws <- posterior_eta_draws(
    fit,
    replicates,
    include_posterior = include_posterior,
    include_timing = include_timing,
    timing_refits = timing_refits,
    seed = seed
  )
  model <- fit$model
  model <- safe_param(model, fit$current_covariates %||% list())
  pmetrics_fit <- identical(fit$engine %||% "mapbayr", "pmetrics")
  if (pmetrics_fit && !is.null(eta_draws$draws)) {
    idata <- data.frame(ID = seq_len(replicates), eta_draws$draws, check.names = FALSE)
    model <- mrgsolve::idata_set(model, idata)
  } else if (!is.null(fit$estimate)) {
    model <- mrgsolve::zero_re(model, omega)
    idata <- data.frame(ID = seq_len(replicates), eta_draws$draws, check.names = FALSE)
    model <- mrgsolve::idata_set(model, idata)
  } else if (!include_posterior) {
    model <- mrgsolve::zero_re(model, omega)
  }
  fixed_residual_cv <- suppressWarnings(as.numeric(fit$residual_error$fixed_cv %||% NA_real_))
  use_fixed_residual <- include_residual && identical(fit$residual_error$mode %||% "model", "fixed_cv") && is.finite(fixed_residual_cv)
  if (!pmetrics_fit && (!include_residual || use_fixed_residual)) model <- mrgsolve::zero_re(model, sigma)
  horizon <- max(24, interval)
  warmup_doses <- if (isTRUE(fit$split_lego)) LEGO_STEADY_STATE_WARMUP_DOSES else 0L
  warmup_time <- warmup_doses * interval
  event <- mrgsolve::ev(
    amt = dose,
    ii = interval / hours_per_model_time,
    cmt = fit$contract$adm_cmt,
    ss = if (isTRUE(fit$split_lego)) 0 else 1,
    addl = warmup_doses + max(0, floor((horizon - 1e-8) / interval)),
    rate = if (infusion > 0) dose / infusion * hours_per_model_time else 0
  )
  simulation <- model |>
    mrgsolve::ev(event) |>
    mrgsolve::mrgsim(start = warmup_time / hours_per_model_time, end = (warmup_time + horizon) / hours_per_model_time,
      delta = delta / hours_per_model_time, nid = replicates, recsort = 3) |>
    as.data.frame()
  simulation$time <- simulation$time * hours_per_model_time - warmup_time
  simulation <- simulation[simulation$time >= -1e-8 & simulation$time <= horizon + 1e-8, , drop = FALSE]
  column <- concentration_column %||% pick_concentration_column(simulation)
  if (!column %in% names(simulation)) stop("Concentration output is unavailable.")
  if (return_profiles) {
    values <- simulation[[column]] * concentration_scale
    if (any(!is.finite(values)) || any(values < -1e-6)) stop("Invalid simulated concentrations.")
    return(list(profiles = data.frame(ID = simulation$ID, time = simulation$time, concentration = pmax(0, values)),
      posterior_available = eta_draws$available, uncertainty_mode = eta_draws$mode))
  }
  use_pmetrics_residual <- pmetrics_fit && isTRUE(include_residual)
  if (use_fixed_residual || use_pmetrics_residual) set.seed(seed + 41L)
  rows <- lapply(split(simulation, simulation$ID), function(profile) {
    values <- pmax(0, profile[[column]])
    if (use_pmetrics_residual) {
      values <- pmetrics_residual_values(values, fit$residual_error$pmetrics)
    } else if (use_fixed_residual) {
      values <- pmax(0, values * (1 + stats::rnorm(length(values), 0, fixed_residual_cv / 100)))
    }
    profile_metrics(
      data.frame(time = profile$time, concentration = values),
      interval = interval,
      concentration_threshold = concentration_threshold
    )
  })
  output <- do.call(rbind, rows)
  output$model <- fit$id
  output$uncertainty_mode <- eta_draws$mode
  rownames(output) <- NULL
  list(
    data = output,
    posterior_available = eta_draws$available,
    timing_available = eta_draws$timing_available,
    uncertainty_mode = eta_draws$mode
  )
}

simulate_regimen_distribution <- function(
  fits,
  weights,
  dose,
  interval,
  infusion,
  metric = "AUC24",
  replicates = 250,
  interval_level = 90,
  delta = 0.1,
  include_posterior = TRUE,
  include_residual = FALSE,
  include_timing = FALSE,
  timing_refits = 20L,
  target_low = NA_real_,
  target_high = NA_real_,
  target_options = list(),
  seed = 381
) {
  valid <- successful_fits(fits)
  replicates <- max(20L, min(1000L, as.integer(replicates)))
  probabilities <- weights[names(valid)]
  probabilities[!is.finite(probabilities) | probabilities < 0] <- 0
  if (sum(probabilities) <= 0) probabilities[] <- 1
  probabilities <- probabilities / sum(probabilities)
  set.seed(seed)
  selected <- sample(names(valid), size = replicates, replace = TRUE, prob = probabilities)
  allocations <- table(factor(selected, levels = names(valid)))

  simulations <- lapply(seq_along(valid), function(index) {
    name <- names(valid)[[index]]
    count <- as.integer(allocations[[name]])
    if (count <= 0) return(NULL)
    simulate_model_distribution(
      valid[[name]], dose, interval, infusion, count, delta,
      include_posterior = include_posterior,
      include_residual = include_residual,
      include_timing = include_timing,
      timing_refits = timing_refits,
      concentration_threshold = target_options$concentration_threshold %||% NA_real_,
      seed = seed + index * 101L
    )
  })
  simulations <- Filter(Negate(is.null), simulations)
  distribution <- do.call(rbind, lapply(simulations, `[[`, "data"))
  scores <- score_target_metrics(distribution, metric, target_low, target_high, target_options)
  distribution$value <- scores$value
  distribution$in_target <- scores$in_target
  alpha <- (100 - interval_level) / 200
  quantiles <- stats::quantile(distribution$value, probs = c(alpha, 0.5, 1 - alpha), na.rm = TRUE, names = FALSE)
  component_quantiles <- lapply(c("cmin", "cmax"), function(name) {
    stats::setNames(
      stats::quantile(distribution[[name]], probs = c(alpha, 0.5, 1 - alpha), na.rm = TRUE, names = FALSE),
      c("lower", "median", "upper")
    )
  })
  names(component_quantiles) <- c("cmin", "cmax")
  list(
    data = distribution,
    lower = quantiles[[1]],
    median = quantiles[[2]],
    upper = quantiles[[3]],
    interval_level = interval_level,
    metric = metric,
    component_quantiles = component_quantiles,
    include_posterior = include_posterior,
    include_residual = include_residual,
    residual_error_mode = if (include_residual && length(valid)) valid[[1]]$residual_error$mode %||% "model" else "none",
    fixed_residual_cv = if (include_residual && length(valid)) valid[[1]]$residual_error$fixed_cv %||% NA_real_ else NA_real_,
    include_timing = include_timing,
    posterior_available = all(vapply(simulations, `[[`, logical(1), "posterior_available")),
    timing_available = any(vapply(simulations, `[[`, logical(1), "timing_available")),
    uncertainty_modes = unique(vapply(simulations, `[[`, character(1), "uncertainty_mode"))
  )
}

rank_regimens_by_pta <- function(
  recommendations,
  fits,
  weights,
  metric,
  target_low,
  target_high,
  target_options = list(),
  replicates = 100,
  delta = 0.1,
  include_posterior = TRUE,
  include_residual = FALSE,
  top_n = 12L,
  seed = 773
) {
  output <- recommendations
  output$p_under <- NA_real_
  output$p_target <- NA_real_
  output$p_over <- NA_real_
  output$pta_evaluated <- FALSE
  candidates <- seq_len(min(nrow(output), max(1L, as.integer(top_n))))
  for (index in candidates) {
    distribution <- simulate_regimen_distribution(
      fits = fits,
      weights = weights,
      dose = output$dose[[index]],
      interval = output$interval[[index]],
      infusion = output$infusion[[index]],
      metric = metric,
      replicates = replicates,
      interval_level = 90,
      delta = delta,
      include_posterior = include_posterior,
      include_residual = include_residual,
      include_timing = FALSE,
      target_low = target_low,
      target_high = target_high,
      target_options = target_options,
      seed = seed
    )
    values <- distribution$data$value
    output$p_target[[index]] <- mean(distribution$data$in_target, na.rm = TRUE)
    if (!identical(metric, "CminCmax")) {
      output$p_under[[index]] <- mean(values < target_low, na.rm = TRUE)
      output$p_over[[index]] <- mean(values > target_high, na.rm = TRUE)
    }
    output$pta_evaluated[[index]] <- TRUE
  }
  evaluated <- output[output$pta_evaluated, , drop = FALSE]
  unevaluated <- output[!output$pta_evaluated, , drop = FALSE]
  evaluated <- evaluated[order(-evaluated$p_target, evaluated$p_over, evaluated$p_under, evaluated$distance, evaluated$center_distance), , drop = FALSE]
  output <- rbind(evaluated, unevaluated)
  rownames(output) <- NULL
  output
}

model_averaging_sensitivity <- function(
  fits,
  selected_weights,
  dose_min,
  dose_max,
  dose_step,
  intervals,
  infusion,
  metric,
  target_low,
  target_high,
  target_options = list(),
  delta = 0.1
) {
  valid <- successful_fits(fits)
  if (!length(valid)) return(data.frame())
  normalize <- function(values) {
    values <- values[names(valid)]
    values[!is.finite(values) | values < 0] <- 0
    if (sum(values) <= 0) values[] <- 1
    values / sum(values)
  }
  sets <- list()
  sets[["Pond\u00e9ration s\u00e9lectionn\u00e9e"]] <- normalize(selected_weights)
  if (length(valid) > 1L) {
    sets[["Poids égaux"]] <- stats::setNames(rep(1 / length(valid), length(valid)), names(valid))
    sets[["AIC"]] <- normalize(compute_model_weights(fits, "AIC"))
    sets[["Log-vraisemblance"]] <- normalize(compute_model_weights(fits, "LL"))
    for (excluded in names(valid)) {
      remaining <- setdiff(names(valid), excluded)
      if (!length(remaining)) next
      values <- selected_weights
      values[excluded] <- 0
      sets[[paste0("Sans ", valid[[excluded]]$label)]] <- normalize(values)
    }
  }
  keys <- vapply(sets, function(values) paste(round(values, 8), collapse = ":"), character(1))
  sets <- sets[!duplicated(keys)]
  rows <- lapply(names(sets), function(label) {
    weights <- sets[[label]]
    recommendation <- recommend_regimens(
      fits, weights, dose_min, dose_max, dose_step, intervals, infusion,
      metric, target_low, target_high, target_options, delta
    )[1, , drop = FALSE]
    data.frame(
      scenario = label,
      dose = recommendation$dose,
      interval = recommendation$interval,
      auc24 = recommendation$auc24,
      cmin = recommendation$cmin,
      cmax = recommendation$cmax,
      time_above_pct = recommendation$time_above_pct,
      target_value = recommendation$value,
      in_target = recommendation$in_target,
      weights = paste(paste(names(weights), round(weights, 3), sep = "="), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

fit_profiles <- function(fits, weights, end_time, delta = 0.25) {
  valid <- successful_fits(fits)
  profiles <- list()
  for (name in names(valid)) {
    fit <- valid[[name]]
    has_individual_fit <- !is.null(fit$estimate) ||
      (identical(fit$engine %||% "mapbayr", "pmetrics") && isTRUE(fit$pmetrics$mapped))
    if (!has_individual_fit) next
    profile <- simulate_known_history(fit, end_time = end_time, delta = delta)[, c("time", "concentration"), drop = FALSE]
    profile$model <- name
    profiles[[name]] <- profile
  }
  if (!length(profiles)) return(list(per_model = data.frame(), average = data.frame()))
  list(
    per_model = do.call(rbind, profiles),
    average = average_profiles(profiles, weights)
  )
}

model_summary <- function(fits, weights) {
  rows <- lapply(names(fits), function(name) {
    fit <- fits[[name]]
    if (inherits(fit, "tdm_fit_error")) {
      return(data.frame(
        model = fit$label,
        status = "Error",
        weight = 0,
        clearance = NA_real_,
        eta = "",
        ml = "Non appliqué",
        detail = fit$message,
        stringsAsFactors = FALSE
      ))
    }
    pmetrics_fit <- identical(fit$engine %||% "mapbayr", "pmetrics")
    clearance <- tryCatch({
      if (pmetrics_fit) {
        values <- fit$posterior_parameters
        if ("cl" %in% names(values)) values[["cl"]] else if (all(c("ke", "v") %in% names(values))) values[["ke"]] * values[["v"]] else NA_real_
      } else if (is.null(fit$estimate)) NA_real_ else as.numeric(mapbayr::get_param(fit$estimate, "CL"))[[1]]
    }, error = function(error) NA_real_)
    eta <- tryCatch({
      if (pmetrics_fit) {
        paste(paste(names(fit$posterior_parameters), round(fit$posterior_parameters, 4), sep = "="), collapse = ", ")
      } else if (is.null(fit$estimate)) "Population" else paste(round(mapbayr::get_eta(fit$estimate), 4), collapse = ", ")
    }, error = function(error) "")
    ml <- fit$ml_correction %||% list(applied = FALSE, reason = "not_evaluated")
    ml_label <- if (isTRUE(ml$applied)) {
      if (identical(ml$type %||% "", "auc24_direct")) {
        paste0(ml$artifact_id, " · AUC24 ", round(ml$auc24, 1), " ", ml$unit %||% "mg.h/L")
      } else {
        paste0(ml$artifact_id, " · ", ml$eta, " ", if (ml$applied_delta >= 0) "+" else "", round(ml$applied_delta, 4))
      }
    } else {
      "Repli MAP"
    }
    data.frame(
      model = fit$label,
      status = if (pmetrics_fit) {
        if (isTRUE(fit$pmetrics$mapped)) "MAP Pmetrics" else "Population Pmetrics"
      } else if (is.null(fit$estimate)) "Population" else "MAP",
      weight = weights[[name]] %||% 0,
      clearance = clearance,
      eta = eta,
      ml = ml_label,
      detail = paste0(fit$contract$n_eta, " ETA; ", fit$contract$n_sigma, " residual terms"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
