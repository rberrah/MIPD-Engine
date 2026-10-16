infection_config <- function(config = list()) {
  defaults <- list(source = "pk", metric = "time", basis = "free", fu = 1, mic = 1,
    multiple = 1, target = 100, replicates = 250, delta = 0.05,
    dose = 1000, interval = 12, infusion = 1, dose2 = 1000, interval2 = 8, infusion2 = 1, grid = NULL, pta_target = 90,
    exposure = "pk", v = 20, cl = 4)
  x <- utils::modifyList(defaults, config)
  x$mic_kind <- NULL
  if (!x$exposure %in% c("iv1", "pk")) stop("Unknown PK exposure.")
  x$v <- ddi_numeric(x$v, "V", 1e-6, 1e6); x$cl <- ddi_numeric(x$cl, "CL", 1e-6, 1e6)
  for (name in c("source", "metric", "basis")) {
    choices <- switch(name, source = c("pk", "tdm"), metric = c("time", "auc", "peak"),
      basis = c("free", "total"))
    if (length(x[[name]]) != 1 || !x[[name]] %in% choices) stop("Invalid infectiology option: ", name)
  }
  limits <- list(fu = c(1e-6, 1), mic = c(1e-6, 1e5), multiple = c(0.01, 1000), target = c(0.001, if (x$metric == "time") 100 else 1e6),
    replicates = c(50, 1000), delta = c(0.01, 0.25), dose = c(1e-6, 1e7), interval = c(0.25, 8760), infusion = c(0, 168),
    dose2 = c(1e-6, 1e7), interval2 = c(0.25, 8760), infusion2 = c(0, 168), pta_target = c(.001, 100))
  for (name in names(limits)) x[[name]] <- ddi_numeric(x[[name]], name, limits[[name]][1], limits[[name]][2])
  if (x$replicates != floor(x$replicates)) stop("Replicates must be an integer.")
  if (x$infusion > x$interval || x$infusion2 > x$interval2) stop("Infusion duration must not exceed the dosing interval.")
  if (x$replicates * (max(24, x$interval, x$interval2) / x$delta + 1) > 3000000)
    stop("Simulation exceeds 3 million samples per regimen. Reduce draws or increase the time step.")
  infection_regimens(x)
  x
}

infection_regimens <- function(config) {
  current <- data.frame(regimen = "current", dose = config$dose, interval = config$interval, infusion = config$infusion)
  if (is.null(config$grid)) return(data.frame(regimen = c("1", "2"), dose = c(config$dose, config$dose2),
    interval = c(config$interval, config$interval2), infusion = c(config$infusion, config$infusion2)))
  g <- config$grid
  minimum <- ddi_numeric(g$min, "Dose min", 1e-6, 1e7)
  maximum <- ddi_numeric(g$max, "Dose max", minimum, 1e7)
  step <- ddi_numeric(g$step, "Dose step", 1e-6, 1e7)
  intervals <- sort(unique(as.numeric(unlist(g$intervals))))
  infusion <- ddi_numeric(g$infusion, "Candidate infusion", 0, 168)
  if (!length(intervals) || any(!is.finite(intervals) | intervals < .25 | intervals > 8760)) stop("Select valid dosing intervals between 0.25 and 8760 hours.")
  count <- (floor((maximum - minimum)/step) + 1) * length(intervals)
  if (count > 24) stop("Maximum 24 candidate regimens. Narrow the dose grid.")
  if (infusion > min(intervals)) stop("Infusion duration must not exceed any candidate dosing interval.")
  candidates <- expand.grid(dose = seq(minimum, maximum, by = step), interval = intervals)
  candidates$infusion <- if (isTRUE(g$continuous)) candidates$interval else infusion
  candidates$regimen <- paste0("candidate-", seq_len(nrow(candidates)))
  regimens <- rbind(current, candidates[names(current)])
  points <- config$replicates * (pmax(24, regimens$interval) / config$delta + 1)
  if (any(points > 3000000)) stop("Maximum 3 million samples per regimen. Reduce draws or increase the time step.")
  if (sum(points) > 12000000)
    stop("Grid exceeds 12 million samples. Reduce candidate regimens or Monte Carlo draws.")
  regimens
}

infection_rank <- function(result) {
  grid <- result$regimens
  at_mic <- result$curves[result$curves$mic == result$config$mic, ]
  grid$daily_dose <- grid$dose * 24 / grid$interval
  grid$pta <- at_mic$pta[match(grid$regimen, at_mic$regimen)]
  grid$target_met <- grid$pta >= result$config$pta_target
  grid[order(grid$regimen != "current", !grid$target_met, grid$daily_dose, -grid$pta), ]
}

infection_index <- function(profile, mic, config, window) {
  profile <- profile[c("time", "concentration")]
  if (!nrow(profile) || any(!is.finite(as.matrix(profile[c("time", "concentration")])))) stop("Invalid concentration profile.")
  if (any(!is.finite(mic) | mic <= 0)) stop("MIC values must be exact positive numbers, not censored limits.")
  profile <- profile[order(profile$time), ]
  if (min(profile$time) > window[1] + 1e-7 || max(profile$time) < window[2] - 1e-7 || diff(window) <= 0) stop("Incomplete exposure window.")
  # Keep both sides of instantaneous dose events; their zero-duration segments add no area.
  for (edge in window) if (!any(abs(profile$time - edge) < 1e-8)) {
    profile <- rbind(profile, data.frame(time = edge,
      concentration = stats::approx(profile$time, profile$concentration, edge, ties = "ordered")$y))
  }
  profile <- profile[order(profile$time), ]
  profile <- profile[profile$time >= window[1] - 1e-8 & profile$time <= window[2] + 1e-8, ]
  concentration <- profile$concentration * if (config$basis == "free") config$fu else 1
  switch(config$metric,
    time = vapply(mic, function(m) time_above_threshold_pct(profile$time, concentration, config$multiple * m), numeric(1)),
    auc = {
      if (abs(diff(window) - 24) > 1e-6) return(rep(NA_real_, length(mic)))
      trap_auc(profile$time, concentration) / mic
    },
    peak = max(concentration) / mic)
}

infection_population_fit <- function(context, soloc, cache) {
  # The shared PD selector intentionally zeros random effects. Restore the original OMEGA for PTA.
  model <- if (context$id == "custom-pk") compile_model(custom_code = context$code, allow_custom = ALLOW_CUSTOM_MODELS,
    custom_soloc = soloc, custom_cache = cache) else compile_model(model_id = context$id)
  model <- safe_param(model, as.list(mrgsolve::param(context$model)))
  list(id = context$id, label = context$label, model = model, estimate = NULL,
    contract = list(adm_cmt = context$adm_cmt), split_lego = is_split_lego_model(model),
    current_covariates = list(), concentration_column = context$concentration,
    administration_modes = context$administration_modes %||% character(),
    concentration_scale = context$concentration_scale, hours_per_model_time = if (context$time_unit == "day") 24 else 1)
}

infection_simulate <- function(fits, weights, config, route, seed = 831L) {
  config <- infection_config(config)
  valid <- successful_fits(fits)
  if (!length(valid) || length(valid) != length(fits)) stop("All selected PK models must be available.")
  if (!route %in% c("IV", "Oral")) stop("Unknown administration route.")
  regimens <- infection_regimens(config)
  for (fit in valid) for (r in seq_len(nrow(regimens))) {
    modes <- fit$administration_modes %||% fit$contract$mode %||% character()
    if (route == "Oral" && regimens$infusion[r] > 0) stop("Oral administration requires infusion = 0.")
    if (route == "IV" && length(modes)) {
      continuous <- "IV_CONTINUOUS" %in% modes && abs(regimens$infusion[r] - regimens$interval[r]) < 1e-8
      intermittent <- "IV_INTERMITTENT" %in% modes && regimens$infusion[r] < regimens$interval[r] - 1e-8
      if (!continuous && !intermittent) stop("Regimen ", r, " is incompatible with the model's IV administration mode.")
    }
  }
  posterior <- vapply(valid, function(f) !is.null(f$estimate), logical(1))
  if (any(posterior) && !all(posterior)) stop("Do not mix population and posterior distributions.")
  for (fit in valid[posterior]) {
    covariance <- fit$estimate$covariance[[1]]
    if (!is.matrix(covariance) || any(!is.finite(covariance)) || !nrow(covariance) ||
        min(eigen((covariance + t(covariance))/2, symmetric = TRUE, only.values = TRUE)$values) <= 0)
      stop("A positive-definite posterior covariance is required for individual target probabilities.")
    if (!is.null(fit$ml_eta_override)) stop("Experimental ML-adjusted ETA cannot replace the MAP posterior in this module.")
  }
  if (!any(posterior) && any(vapply(valid, function(f) sum(abs(as.matrix(mrgsolve::omat(f$model)))) <= 0, logical(1))))
    stop("PTA requires a PK model with nonzero interindividual variability (OMEGA). A fixed profile is not a population distribution.")
  probabilities <- weights[names(valid)]
  if (length(probabilities) != length(valid) || any(!is.finite(probabilities) | probabilities < 0) || sum(probabilities) <= 0)
    stop("Invalid model averaging weights.")
  probabilities <- probabilities / sum(probabilities)
  set.seed(seed)
  allocations <- table(factor(sample(names(valid), config$replicates, replace = TRUE, prob = probabilities), levels = names(valid)))
  mic <- sort(unique(c(2^seq(-9, 9), config$mic / 2, config$mic, config$mic * 2)))
  curves <- list(); exposures <- list()
  for (r in seq_len(nrow(regimens))) {
    reg <- regimens[r, ]; indices <- list(); trace <- list()
    for (i in seq_along(valid)) {
      count <- as.integer(allocations[[i]]); if (!count) next
      fit <- valid[[i]]
      set.seed(seed + i * 101L)
      result <- simulate_model_distribution(fit, reg$dose, reg$interval, reg$infusion, count, config$delta,
        include_posterior = TRUE, include_residual = FALSE, seed = seed + i * 101L, return_profiles = TRUE,
        hours_per_model_time = fit$hours_per_model_time %||% 1,
        concentration_column = fit$concentration_column, concentration_scale = fit$concentration_scale %||% 1)
      if (any(posterior) && !isTRUE(result$posterior_available)) stop("Posterior covariance unavailable: a MAP point estimate cannot be displayed as a PTA.")
      profiles <- split(result$profiles[c("time", "concentration")], result$profiles$ID)
      window <- c(0, if (config$metric == "auc") 24 else reg$interval)
      indices[[i]] <- do.call(rbind, lapply(profiles, infection_index, mic = mic, config = config, window = window))
      trace[[i]] <- result$profiles
      trace[[i]]$subject <- paste(i, trace[[i]]$ID, sep = "-")
    }
    values <- do.call(rbind, indices)
    if (any(!is.finite(values))) stop("Target index is not calculable for all simulated subjects.")
    n <- nrow(values); p <- colMeans(values >= config$target - 1e-9)
    # Wilson bounds describe Monte Carlo precision, not uncertainty in MIC, fu or model structure.
    z <- stats::qnorm(.975); centre <- (p + z^2 / (2*n)) / (1 + z^2/n)
    width <- z * sqrt(p*(1-p)/n + z^2/(4*n^2)) / (1 + z^2/n)
    curves[[r]] <- data.frame(regimen = reg$regimen, mic = mic, pta = 100*p,
      lower = 100*pmax(0, centre-width), upper = 100*pmin(1, centre+width), n = n)
    all <- do.call(rbind, trace)
    all <- all[!duplicated(all[c("subject", "time")], fromLast = TRUE), ]
    summaries <- lapply(split(all$concentration, all$time), stats::quantile, probs = c(.05, .5, .95), names = FALSE)
    exposures[[r]] <- data.frame(regimen = reg$regimen, time = as.numeric(names(summaries)), do.call(rbind, summaries))
    names(exposures[[r]])[3:5] <- c("lower", "median", "upper")
    exposures[[r]] <- exposures[[r]][order(exposures[[r]]$time), ]
  }
  list(curves = do.call(rbind, curves), exposure = do.call(rbind, exposures), regimens = regimens,
    config = config, posterior = all(posterior), models = vapply(valid, `[[`, character(1), "label"), weights = probabilities,
    parameters = lapply(valid, function(f) list(parameters = as.list(mrgsolve::param(f$model)), covariates = f$current_covariates,
      hours_per_model_time = f$hours_per_model_time %||% 1, concentration_scale = f$concentration_scale %||% 1,
      concentration_column = f$concentration_column %||% "automatic")),
    route = route, seed = seed, created = format(Sys.time(), tz = "UTC", usetz = TRUE))
}

infection_tdm_summary <- function(result, config, scale = 1) {
  config <- infection_config(config)
  scale <- ddi_numeric(scale, "Concentration conversion to mg/L", 1e-9, 1e9)
  exposure <- historical_exposure(result$fits, result$weights, result$doses, result$observations, delta = config$delta)
  profile <- exposure$profile
  profile$concentration <- profile$concentration * scale
  window <- c(exposure$window_start, exposure$decision_time)
  value <- if (diff(window) > 0) infection_index(profile, config$mic, config, window) else NA_real_
  list(value = value, window = window, profile = profile, config = config,
    concentration = exposure$concentration * scale, achieved = if (is.finite(value)) value >= config$target else NA)
}
