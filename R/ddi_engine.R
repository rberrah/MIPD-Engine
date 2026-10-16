ddi_numeric <- function(value, label, minimum = -Inf, maximum = Inf) {
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value) || value < minimum || value > maximum) {
    stop(label, " is outside the allowed range.")
  }
  value
}

ddi_library_context <- function(model_id, route = NULL) {
  record <- model_record(model_id)
  routes <- model_routes(record)
  route <- as.character(route %||% routes[[1]])
  if (!route %in% routes) stop("The selected route is not supported by this model.")

  model <- compile_model(model_id = model_id)
  model <- mrgsolve::zero_re(model)
  adm_cmt <- match(model_administration_cmt(record, route), model@cmtL)
  if (!is.finite(adm_cmt)) stop("The administration compartment is missing from the compiled model.")

  list(
    id = model_id,
    label = record$label[[1]],
    drug = record$drug[[1]],
    model = model,
    code = read_library_code(model_id),
    route = route,
    administration_modes = model_administration_modes(record, route),
    adm_cmt = adm_cmt,
    split_lego = is_split_lego_model(model),
    covariate_names = parse_covariates(read_library_code(model_id))$name,
    source = "population",
    citation = record$citation[[1]] %||% "",
    doi = record$doi[[1]] %||% ""
  )
}

ddi_fit_context <- function(fit) {
  library_model <- (fit$id %||% "") %in% MODEL_CATALOG$id
  record <- if (library_model) model_record(fit$id) else NULL
  list(
    id = fit$id %||% "custom",
    label = fit$label %||% if (library_model) record$label[[1]] else "Lego session model",
    drug = if (library_model) record$drug[[1]] else "Lego",
    model = individual_model(fit),
    code = if (library_model) read_library_code(fit$id) else paste(fit$model@code, collapse = "\n"),
    route = fit$contract$route %||% "Unspecified",
    administration_modes = fit$contract$mode %||% character(),
    adm_cmt = fit$contract$adm_cmt,
    split_lego = isTRUE(fit$split_lego),
    covariate_names = names(fit$current_covariates %||% list()),
    source = "tdm",
    citation = if (library_model) record$citation[[1]] %||% "" else "Session Lego model",
    doi = if (library_model) record$doi[[1]] %||% "" else ""
  )
}

ddi_filter_parameters <- function(values, covariate_names = character()) {
  excluded <- unique(c(
    covariate_names,
    grep("^(ETA|BETA|EPS|SIGMA|OMEGA|RUV|ERR|RES)", names(values), value = TRUE, ignore.case = TRUE),
    grep("(^|_)(PROP|ADD|CV)(_|$)", names(values), value = TRUE, ignore.case = TRUE)
  ))
  keep <- is.finite(values) & values > 0 & !names(values) %in% excluded
  data.frame(name = names(values)[keep], value = as.numeric(values[keep]), stringsAsFactors = FALSE)
}

ddi_parameter_table <- function(context) {
  values <- unlist(as.list(mrgsolve::param(context$model)), use.names = TRUE)
  ddi_filter_parameters(values, context$covariate_names %||% character())
}

ddi_code_parameter_table <- function(code) {
  lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  values <- numeric()
  in_parameters <- FALSE
  covariate_block <- FALSE
  for (line in lines) {
    if (grepl("^\\s*(\\$PARAM|\\[PARAM\\])", line, ignore.case = TRUE, perl = TRUE)) {
      in_parameters <- TRUE
      covariate_block <- grepl("@covariates", tolower(line), fixed = TRUE)
      next
    }
    if (grepl("^\\s*(\\$|\\[)[A-Za-z]", line, perl = TRUE)) {
      in_parameters <- FALSE
      covariate_block <- FALSE
      next
    }
    if (!in_parameters || covariate_block) next
    clean <- trimws(sub("//.*$", "", line))
    match <- regexec("^([A-Za-z_][A-Za-z0-9_]*)\\s*:\\s*([-+]?(?:[0-9]*\\.?[0-9]+)(?:[eE][-+]?[0-9]+)?)", clean, perl = TRUE)
    fields <- regmatches(clean, match)[[1]]
    if (length(fields) == 3L) values[[fields[[2]]]] <- as.numeric(fields[[3]])
  }
  ddi_filter_parameters(values, parse_covariates(code)$name)
}

ddi_context_with_route <- function(context, route) {
  if (!(context$id %||% "") %in% MODEL_CATALOG$id) return(context)
  record <- model_record(context$id)
  route <- as.character(route %||% context$route)
  if (!route %in% model_routes(record)) stop("The selected route is not supported by this model.")
  context$route <- route
  context$adm_cmt <- match(model_administration_cmt(record, route), context$model@cmtL)
  if (!is.finite(context$adm_cmt)) stop("The administration compartment is missing from the compiled model.")
  context
}

ddi_validate_regimen <- function(regimen, label) {
  regimen$dose <- ddi_numeric(regimen$dose, paste(label, "dose"), 0.001, 1e7)
  regimen$interval <- ddi_numeric(regimen$interval, paste(label, "interval"), 0.25, 720)
  regimen$infusion <- ddi_numeric(regimen$infusion %||% 0, paste(label, "infusion"), 0, regimen$interval)
  regimen
}

ddi_dose_data <- function(context, regimen, start, stop, steady_state = FALSE) {
  times <- if (stop < start) numeric() else seq(start, stop + 1e-8, by = regimen$interval)
  if (!length(times)) return(data.frame())

  if (isTRUE(steady_state) && isTRUE(context$split_lego)) {
    times <- c(seq(-LEGO_STEADY_STATE_WARMUP_DOSES * regimen$interval, -regimen$interval, by = regimen$interval), times)
  }
  rows <- data.frame(
    ID = 1,
    time = times,
    evid = 1,
    cmt = context$adm_cmt,
    amt = regimen$dose,
    rate = if (regimen$infusion > 0) regimen$dose / regimen$infusion else 0,
    ii = 0,
    addl = 0,
    ss = 0,
    stringsAsFactors = FALSE
  )
  if (isTRUE(steady_state) && !isTRUE(context$split_lego)) {
    index <- which.min(abs(rows$time))
    rows$ii[index] <- regimen$interval
    rows$ss[index] <- 1
  }
  rows
}

ddi_simulate_profile <- function(context, regimen, start, stop, horizon, delta, steady_state = FALSE,
                                 parameter_name = NULL, parameter_times = NULL, parameter_values = NULL) {
  doses <- ddi_dose_data(context, regimen, start, stop, steady_state)
  simulation_start <- min(0, doses$time, na.rm = TRUE)
  grid <- seq(simulation_start, horizon, by = delta)
  observations <- data.frame(
    ID = 1,
    time = grid,
    evid = 0,
    cmt = 0,
    amt = 0,
    rate = 0,
    ii = 0,
    addl = 0,
    ss = 0,
    stringsAsFactors = FALSE
  )
  data <- rbind(doses, observations)
  data <- data[order(data$time, -data$evid), , drop = FALSE]

  if (length(parameter_name)) {
    if (any(!parameter_name %in% model_param_names(context$model))) stop("Unknown interaction target parameter.")
    values <- if (is.list(parameter_values)) parameter_values else list(parameter_values)
    if (length(values) != length(parameter_name)) stop("One parameter trajectory is required per interaction target.")
    for (i in seq_along(parameter_name)) {
      data[[parameter_name[[i]]]] <- stats::approx(
        parameter_times,
        values[[i]],
        xout = data$time,
        rule = 2,
        ties = "ordered"
      )$y
    }
  }

  simulation <- mrgsolve::mrgsim_d(
    context$model,
    data = data,
    start = simulation_start,
    end = horizon,
    delta = delta,
    recsort = 3,
    nocb = FALSE
  ) |>
    as.data.frame()
  simulation <- simulation[simulation$time >= -1e-8, , drop = FALSE]
  simulation <- simulation[!duplicated(simulation$time, fromLast = TRUE), , drop = FALSE]
  concentration <- pick_concentration_column(simulation)
  data.frame(
    time = simulation$time,
    day = simulation$time / 24,
    concentration = pmax(0, simulation[[concentration]]),
    stringsAsFactors = FALSE
  )
}

ddi_validate_targets <- function(config, available = NULL) {
  targets <- config$targets
  if (is.null(targets) || !length(targets)) {
    targets <- list(list(parameter = config$target %||% "", fraction = config$target_fraction %||% 1))
  } else if (is.data.frame(targets)) {
    targets <- lapply(seq_len(nrow(targets)), function(i) as.list(targets[i, , drop = FALSE]))
  } else if (!is.null(targets$parameter) || !is.null(targets$target) || !is.null(targets$name)) {
    targets <- list(targets)
  }
  if (!is.list(targets) || !length(targets) || length(targets) > 6L) stop("Select between one and six interaction targets.")
  targets <- lapply(seq_along(targets), function(i) {
    target <- targets[[i]]
    parameter <- as.character(target$parameter %||% target$target %||% target$name %||% "")
    if (length(parameter) != 1L || !grepl("^[A-Za-z_][A-Za-z0-9_]*$", parameter)) stop("Invalid interaction target.")
    list(parameter = parameter, fraction = ddi_numeric(target$fraction %||% 1, paste("Affected fraction", i), 0, 1))
  })
  parameters <- vapply(targets, `[[`, character(1), "parameter")
  if (anyDuplicated(parameters)) stop("Each interaction target must be unique.")
  if (!is.null(available) && any(!parameters %in% available)) stop("Select positive structural parameters from model 1.")
  config$targets <- targets
  config$target <- targets[[1]]$parameter
  config$target_fraction <- targets[[1]]$fraction
  config
}

ddi_validate_mechanism <- function(config) {
  if (!config$type %in% c("factor", "inhibition", "induction", "reversible", "hill_inhibition", "tdi", "turnover_induction")) stop("Unknown interaction type.")
  config$factor <- ddi_numeric(config$factor %||% 1, "Interaction factor", 0.01, 20)
  config$strength <- ddi_numeric(config$strength %||% 1, "Imax / Emax", 0, if (config$type %in% c("inhibition", "hill_inhibition")) 1 else 20)
  config$c50 <- ddi_numeric(config$c50 %||% 1, "IC50 / EC50 / Ki / KI", 1e-9, 1e9)
  config$hill <- ddi_numeric(config$hill %||% 1, "Hill", 0.1, 10)
  config$kdeg <- ddi_numeric(config$kdeg %||% 0.02, "kdeg (1/h)", 1e-6, 10)
  config$kinact <- ddi_numeric(config$kinact %||% 0.1, "kinact (1/h)", 0, 10)
  ddi_validate_targets(config)
}

ddi_modifier <- function(type, concentration, time, start, stop, factor, strength, c50,
                         hill = 1, kdeg = 0.02, kinact = 0.1) {
  concentration <- pmax(0, concentration)
  if (identical(type, "factor")) {
    return(ifelse(time >= start & time < stop, factor, 1))
  }
  occupancy <- concentration / (c50 + concentration)
  if (type == "reversible") return(pmax(0.01, 1 / (1 + concentration / c50)))
  if (type == "hill_inhibition") return(pmax(0.01, 1 - strength * plogis(hill * (log(concentration) - log(c50)))))
  if (type %in% c("tdi", "turnover_induction")) {
    if (length(time) != length(concentration) || any(!is.finite(time)) || any(diff(time) <= 0)) stop("Enzyme turnover requires an ordered concentration time grid.")
    activity <- rep(1, length(time))
    # Exact turnover solution over each left-held concentration step, matching PK coupling.
    for (i in seq_len(length(time) - 1L)) {
      loss <- kdeg + if (type == "tdi") kinact * occupancy[i] else 0
      synthesis <- kdeg * if (type == "turnover_induction") 1 + strength * occupancy[i] else 1
      equilibrium <- synthesis / loss
      activity[i + 1L] <- equilibrium + (activity[i] - equilibrium) * exp(-loss * (time[i + 1L] - time[i]))
    }
    return(pmax(0.01, activity))
  }
  effect <- strength * occupancy
  switch(type, inhibition = pmax(0.01, 1 - effect), induction = 1 + effect, stop("Unknown interaction type."))
}

ddi_profile_window <- function(profile, from, to) {
  inside <- profile$time > from & profile$time < to
  time <- c(from, profile$time[inside], to)
  concentration <- stats::approx(profile$time, profile$concentration, xout = time, rule = 2, ties = "ordered")$y
  list(
    auc = trap_auc(time, concentration),
    cmin = min(concentration, na.rm = TRUE),
    cmax = max(concentration, na.rm = TRUE)
  )
}

ddi_simulate <- function(config, affected_context, driver_context, delta = 0.1) {
  config$affected <- ddi_validate_regimen(config$affected, "Model 1")
  config$driver <- ddi_validate_regimen(config$driver, "Model 2")
  if ((identical(affected_context$route, "Oral") && config$affected$infusion > 0) ||
      (identical(driver_context$route, "Oral") && config$driver$infusion > 0)) stop("Oral dosing requires infusion duration = 0.")
  config$start_day <- ddi_numeric(config$start_day, "Interaction start", 0.25, 180)
  config$stop_day <- ddi_numeric(config$stop_day, "Interaction stop", config$start_day + 0.25, 365)
  config$followup_days <- ddi_numeric(config$followup_days, "Follow-up", 0, 180)
  config <- ddi_validate_mechanism(config)

  parameters <- ddi_parameter_table(affected_context)
  config <- ddi_validate_targets(config, parameters$name)
  target_names <- vapply(config$targets, `[[`, character(1), "parameter")
  target_fractions <- vapply(config$targets, `[[`, numeric(1), "fraction")
  target_indexes <- match(target_names, parameters$name)
  baseline_parameters <- stats::setNames(parameters$value[target_indexes], target_names)
  start <- config$start_day * 24
  stop <- config$stop_day * 24
  horizon <- (config$stop_day + config$followup_days) * 24

  driver <- ddi_simulate_profile(
    driver_context,
    config$driver,
    start = start,
    stop = stop - 1e-8,
    horizon = horizon,
    delta = delta
  )
  modifier <- ddi_modifier(
    config$type,
    driver$concentration,
    driver$time,
    start,
    stop,
    config$factor,
    config$strength,
    config$c50, config$hill, config$kdeg, config$kinact
  )
  target_modifiers <- stats::setNames(lapply(target_fractions, function(fraction) (1 - fraction) + fraction * modifier), target_names)
  target_values <- stats::setNames(lapply(target_names, function(target) pmax(1e-12, baseline_parameters[[target]] * target_modifiers[[target]])), target_names)

  baseline <- ddi_simulate_profile(
    affected_context,
    config$affected,
    start = 0,
    stop = horizon,
    horizon = horizon,
    delta = delta,
    steady_state = TRUE
  )
  interaction <- ddi_simulate_profile(
    affected_context,
    config$affected,
    start = 0,
    stop = horizon,
    horizon = horizon,
    delta = delta,
    steady_state = TRUE,
    parameter_name = target_names,
    parameter_times = driver$time,
    parameter_values = target_values
  )

  window_to <- min(stop, horizon)
  window_from <- max(0, window_to - 24)
  baseline_metrics <- ddi_profile_window(baseline, window_from, window_to)
  interaction_metrics <- ddi_profile_window(interaction, window_from, window_to)
  ratio <- function(name) interaction_metrics[[name]] / baseline_metrics[[name]]

  affected_profile <- rbind(
    transform(baseline, scenario = "baseline"),
    transform(interaction, scenario = "interaction")
  )
  effect <- do.call(rbind, lapply(seq_along(target_names), function(i) data.frame(
    time = driver$time,
    day = driver$day,
    driver_concentration = driver$concentration,
    modifier = target_modifiers[[i]],
    parameter_value = target_values[[i]],
    target = target_names[[i]],
    fraction = target_fractions[[i]],
    stringsAsFactors = FALSE
  )))
  schedule <- data.frame(
    model = c(affected_context$label, driver_context$label),
    source = c(affected_context$source, driver_context$source),
    route = c(affected_context$route, driver_context$route),
    start_day = c(0, config$start_day),
    stop_day = c(config$stop_day + config$followup_days, config$stop_day),
    dose = c(config$affected$dose, config$driver$dose),
    interval = c(config$affected$interval, config$driver$interval),
    infusion = c(config$affected$infusion, config$driver$infusion),
    stringsAsFactors = FALSE
  )

  list(
    config = config,
    affected = affected_context,
    driver = driver_context,
    affected_profile = affected_profile,
    driver_profile = driver,
    effect = effect,
    schedule = schedule,
    metrics = list(
      baseline = baseline_metrics,
      interaction = interaction_metrics,
      auc_ratio = ratio("auc"),
      cmin_ratio = ratio("cmin"),
      cmax_ratio = ratio("cmax"),
      minimum_modifier = min(effect$modifier, na.rm = TRUE),
      maximum_modifier = max(effect$modifier, na.rm = TRUE),
      baseline_parameter = baseline_parameters[[1]],
      baseline_parameters = baseline_parameters
    )
  )
}
