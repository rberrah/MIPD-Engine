# Generic research structures, not drug-specific clinical models or dosing protocols.
ONCO_DEFAULTS <- list(V = 20, CL = 10, T0 = 60, KG = 0.01, CAP = 300, KILL = 0.25,
  EC50 = 1, RES = 0, ANC0 = 4, MTT = 5, GAMMA = 0.17, SLOPE = 0.15)

onco_config <- function(config = list()) {
  config$growth <- config$growth %||% "exponential"
  if (!config$growth %in% c("exponential", "logistic", "gompertz")) stop("Unknown tumor growth law.")
  config$toxicity <- config$toxicity %||% TRUE
  config$free_pk <- isTRUE(config$free_pk)
  if (!is.logical(config$toxicity) || length(config$toxicity) != 1 || is.na(config$toxicity)) stop("Invalid toxicity setting.")
  if (!is.null(config$parameters) && (!is.list(config$parameters) || any(!names(config$parameters) %in% names(ONCO_DEFAULTS)))) stop("Unknown oncology parameters.")
  config$parameters <- utils::modifyList(ONCO_DEFAULTS, config$parameters %||% list())
  for (name in names(ONCO_DEFAULTS)) {
    minimum <- if (name %in% c("KG", "KILL", "RES", "GAMMA", "SLOPE")) 0 else 1e-6
    config$parameters[[name]] <- ddi_numeric(config$parameters[[name]], name, minimum, 1e5)
  }
  config$horizon <- ddi_numeric(config$horizon %||% 84, "Horizon (days)", 1, 730)
  if (!config$free_pk && config$parameters$CL / config$parameters$V > 100) stop("This exploratory engine supports CL/V up to 100 per day.")
  config$decision <- ddi_numeric(config$decision %||% 42, "Decision day", 0, config$horizon - 0.1)
  config$dose <- ddi_numeric(config$dose %||% 100, "Future dose (mg)", 0, 1e5)
  config$interval <- ddi_numeric(config$interval %||% 21, "Cycle interval (days)", 0.25, 180)
  config$infusion <- ddi_numeric(config$infusion %||% 1, "Infusion (h)", 0, min(168, config$interval * 24))
  config$anc_floor <- ddi_numeric(config$anc_floor %||% 1, "Exploratory ANC floor", 0, 100)
  config$tumor_goal <- ddi_numeric(config$tumor_goal %||% 0.8, "Final tumor / baseline ceiling", 0.001, 100)
  config$sigma_tumor <- ddi_numeric(config$sigma_tumor %||% 0.15, "Tumor log-error SD", 0.001, 10)
  config$sigma_anc <- ddi_numeric(config$sigma_anc %||% 0.2, "ANC log-error SD", 0.001, 10)
  history <- config$history %||% data.frame(time = c(0, 21), amount = c(100, 100), infusion = c(1, 1))
  if (!is.data.frame(history) || !all(c("time", "amount", "infusion") %in% names(history))) stop("Dose history requires time,amount,infusion.")
  if (nrow(history) > 500 || !all(vapply(history[c("time", "amount", "infusion")], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(history[c("time", "amount", "infusion")])))) stop("Invalid dose history.")
  if (any(history$time < 0 | history$time >= config$decision | history$amount <= 0 | history$amount > 1e5 | history$infusion < 0 | history$infusion > 168)) {
    stop("Historical doses must be positive and strictly before the decision day; infusion is in hours (0 = IV bolus).")
  }
  config$history <- history[order(history$time), c("time", "amount", "infusion"), drop = FALSE]
  config
}

onco_schedule <- function(config, fraction = 1, delay = 0, interval = config$interval) {
  fraction <- ddi_numeric(fraction, "Dose fraction", 0, 2)
  delay <- ddi_numeric(delay, "Delay (days)", 0, 180)
  interval <- ddi_numeric(interval, "Cycle interval (days)", max(0.25, config$infusion / 24), 180)
  first <- config$decision + delay
  times <- if (first < config$horizon && fraction * config$dose > 0) seq(first, config$horizon - 1e-8, by = interval) else numeric()
  future <- data.frame(time = times, amount = rep(config$dose * fraction, length(times)), infusion = rep(config$infusion, length(times)))
  rbind(config$history, future)
}

onco_exposure <- function(time, doses, p) {
  kel <- p$CL / p$V
  result <- rep(0, length(time))
  for (i in seq_len(nrow(doses))) {
    elapsed <- time - doses$time[i]
    active <- elapsed >= 0
    duration <- doses$infusion[i] / 24
    if (duration == 0) result[active] <- result[active] + doses$amount[i] / p$V * exp(-kel * elapsed[active])
    else {
      u <- elapsed[active]
      result[active] <- result[active] + doses$amount[i] / (duration * p$CL) *
        (-expm1(-kel * pmin(u, duration))) * exp(-kel * pmax(0, u - duration))
    }
  }
  result
}

onco_rhs <- function(time, state, parameters) {
  p <- parameters$p
  tumor <- max(1e-12, state[["TUMOR"]])
  c <- if (is.null(parameters$exposure)) onco_exposure(time, parameters$doses, p) else
    stats::approx(parameters$exposure$time, parameters$exposure$concentration, xout = time, rule = 2)$y
  growth <- switch(parameters$growth, exponential = p$KG, logistic = p$KG * (1 - tumor / p$CAP), gompertz = p$KG * log(p$CAP / tumor))
  kill <- p$KILL * c / (p$EC50 + c) * exp(-p$RES * time)
  derivative <- c(TUMOR = (growth - kill) * tumor)
  if (parameters$toxicity) {
    cells <- pmax(state[c("PROL", "TR1", "TR2", "TR3", "ANC")], 1e-10)
    ktr <- 4 / p$MTT
    effect <- min(1, p$SLOPE * c)
    feedback <- (p$ANC0 / cells[["ANC"]])^p$GAMMA
    derivative <- c(derivative, PROL = ktr * cells[["PROL"]] * ((1 - effect) * feedback - 1),
      TR1 = ktr * (cells[["PROL"]] - cells[["TR1"]]), TR2 = ktr * (cells[["TR1"]] - cells[["TR2"]]),
      TR3 = ktr * (cells[["TR2"]] - cells[["TR3"]]), ANC = ktr * (cells[["TR3"]] - cells[["ANC"]]))
  }
  list(derivative)
}

onco_simulate <- function(config, doses = NULL, times = NULL, pk_context = NULL, exposure = NULL) {
  config <- onco_config(config)
  p <- config$parameters
  if (is.null(doses)) doses <- onco_schedule(config)
  if (config$free_pk && is.null(exposure)) {
    if (is.null(pk_context)) stop("Load the free PK model before oncology simulation.")
    exposure <- pd_pk_profile(pk_context, doses, config$horizon, hours_per_time = 24)
  }
  if (is.null(times)) times <- seq(0, config$horizon, length.out = ceiling(config$horizon * 20) + 1)
  if (any(!is.finite(times)) || any(times < 0 | times > config$horizon)) stop("Simulation times outside horizon.")
  knots <- unique(round(c(doses$time, doses$time + doses$infusion / 24, if (!is.null(exposure)) exposure$time), 10))
  knots <- knots[knots >= 0 & knots <= max(times)]
  # Merge floating-point aliases of the same hour/day knot before restarting LSODA.
  grid <- sort(unique(round(c(0, times, knots), 10)))
  initial <- c(TUMOR = p$T0)
  if (config$toxicity) initial <- c(initial, stats::setNames(rep(p$ANC0, 5), c("PROL", "TR1", "TR2", "TR3", "ANC")))
  solution <- as.data.frame(deSolve::ode(initial, grid, onco_rhs,
    list(p = p, doses = doses, growth = config$growth, toxicity = config$toxicity, exposure = exposure), method = "lsoda",
    events = if (length(knots)) list(func = function(t, y, p) y, time = knots) else NULL,
    hmax = if (config$free_pk) 0.1 else min(0.1, p$V / p$CL), maxsteps = 50000, rtol = 1e-7, atol = 1e-9))
  if (nrow(solution) != length(grid) || any(!is.finite(as.matrix(solution))) || any(as.matrix(solution[, -1, drop = FALSE]) < -1e-6) ||
      any(as.matrix(solution[, -1, drop = FALSE]) > 1e12)) stop("Unstable oncology simulation; check parameters.")
  solution$concentration <- if (is.null(exposure)) onco_exposure(solution$time, doses, p) else
    stats::approx(exposure$time, exposure$concentration, xout = solution$time, rule = 2)$y
  solution$tumor_ratio <- solution$TUMOR / p$T0
  solution
}

onco_read_observations <- function(text) {
  if (!is.character(text) || length(text) != 1 || nchar(text, type = "bytes") > 200000) stop("Observation CSV limit: 200 kB.")
  data <- utils::read.csv(text = text, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("time", "endpoint", "value") %in% names(data)) || !nrow(data) || nrow(data) > 1000) stop("Expected CSV: time,endpoint,value (endpoint: tumor or anc).")
  if (!is.numeric(data$time) || !is.numeric(data$value) || any(!is.finite(data$time)) || any(!is.finite(data$value)) ||
      any(data$time < 0 | data$value <= 0) || any(!data$endpoint %in% c("tumor", "anc"))) stop("Observations require nonnegative days and positive values; endpoints: tumor, anc.")
  data[order(data$time), c("time", "endpoint", "value"), drop = FALSE]
}

onco_fit <- function(config, data, estimate = c("KILL", "SLOPE"), pk_context = NULL) {
  config <- onco_config(config)
  if (!is.data.frame(data) || !all(c("time", "endpoint", "value") %in% names(data)) || any(!is.finite(data$time)) || any(!is.finite(data$value)) ||
      any(data$time < 0 | data$value <= 0) || any(!data$endpoint %in% c("tumor", "anc"))) stop("Complete the oncology observations with nonnegative times and positive values before fitting.")
  allowed <- c("T0", "KG", "KILL", "EC50", "RES", if (config$growth != "exponential") "CAP",
    if (config$toxicity) c("ANC0", "MTT", "GAMMA", "SLOPE"))
  if (!length(estimate) || anyDuplicated(estimate) || any(!estimate %in% allowed) || length(estimate) > 5 || nrow(data) < length(estimate) + 2) stop("Select 1-5 active PD parameters with at least two more observations than parameters.")
  if (any(data$time > config$decision)) stop("Observations after the decision day cannot be used for this adaptation.")
  if (!config$toxicity && any(data$endpoint == "anc")) stop("Enable myelosuppression for ANC observations.")
  toxicity_parameters <- c("ANC0", "MTT", "GAMMA", "SLOPE")
  if ((any(estimate %in% toxicity_parameters) && !any(data$endpoint == "anc")) ||
      (any(!estimate %in% toxicity_parameters) && !any(data$endpoint == "tumor"))) stop("Each fitted endpoint needs corresponding observations.")
  exposure <- if (config$free_pk) pd_pk_profile(pk_context, config$history, config$horizon, hours_per_time = 24) else NULL
  predict <- function(values) {
    candidate <- config
    candidate$parameters[estimate] <- as.list(values)
    simulation <- onco_simulate(candidate, doses = candidate$history, times = sort(unique(c(0, data$time, config$decision))), exposure = exposure)
    result <- stats::approx(simulation$time, simulation$TUMOR, data$time, rule = 2)$y
    if (config$toxicity) result[data$endpoint == "anc"] <- stats::approx(simulation$time, simulation$ANC, data$time[data$endpoint == "anc"], rule = 2)$y
    result
  }
  sigma <- ifelse(data$endpoint == "anc", config$sigma_anc, config$sigma_tumor)
  residual <- function(values) (log(pmax(1e-12, predict(values))) - log(data$value)) / sigma
  initial <- unlist(config$parameters[estimate])
  lower <- rep(1e-8, length(initial)); upper <- rep(1e5, length(initial))
  fit <- minpack.lm::nls.lm(initial, lower = lower, upper = upper, fn = residual, control = minpack.lm::nls.lm.control(maxiter = 120))
  if (!fit$info %in% 1:4) stop("Oncology fit did not converge.")
  config$parameters[estimate] <- as.list(fit$par)
  data$prediction <- predict(fit$par)
  data$residual <- data$value - data$prediction
  singular <- tryCatch(!is.finite(kappa(fit$hessian)) || kappa(fit$hessian) > 1e10, error = function(e) TRUE)
  boundary <- any(fit$par < lower * 10 | fit$par > upper * 0.99)
  list(config = config, data = data, parameters = data.frame(parameter = estimate, initial = initial, fitted = fit$par),
    objective = sum(residual(fit$par)^2), ill_conditioned = singular, boundary = boundary)
}

onco_metrics <- function(config, curves, schedules) {
  do.call(rbind, lapply(names(curves), function(name) {
    data <- curves[[name]]
    future <- data[data$time >= config$decision, ]
    ratio <- tail(data$tumor_ratio, 1)
    nadir <- if (config$toxicity) min(future$ANC) else NA_real_
    below <- if (config$toxicity) sum(diff(future$time) * (head(future$ANC, -1) < config$anc_floor)) else NA_real_
    data.frame(scenario = name, final_tumor_ratio = ratio, future_anc_nadir = nadir, days_below_floor = below,
      tumor_goal_met = ratio <= config$tumor_goal, anc_goal_met = if (config$toxicity) nadir >= config$anc_floor else NA,
      future_dose_mg = sum(schedules[[name]]$amount[schedules[[name]]$time >= config$decision]))
  }))
}

onco_compare <- function(config, fraction = 0.75, delay = 0, interval = config$interval, pk_context = NULL) {
  config <- onco_config(config)
  schedules <- list(maintain = onco_schedule(config), change = onco_schedule(config, fraction, delay, interval))
  curves <- lapply(schedules, function(doses) onco_simulate(config, doses, pk_context = pk_context))
  metrics <- onco_metrics(config, curves, schedules)
  curves$untreated <- onco_simulate(config, config$history[FALSE, ], exposure = data.frame(time = c(0, config$horizon), concentration = 0))
  pk_metadata <- if (is.null(pk_context)) NULL else c(pk_context[intersect(names(pk_context), c("label", "route", "adm_cmt", "time_unit", "concentration", "concentration_scale"))],
    list(parameters = as.list(mrgsolve::param(pk_context$model))))
  list(config = config, curves = curves, schedules = schedules, metrics = metrics, pk = pk_metadata,
    change = list(fraction = fraction, delay = delay, interval = interval))
}

onco_compare_grid <- function(config, grid, delay = 0, pk_context = NULL) {
  config <- onco_config(config)
  minimum <- ddi_numeric(grid$min, "Dose min", 0, 1e5)
  maximum <- ddi_numeric(grid$max, "Dose max", minimum, 1e5)
  step <- ddi_numeric(grid$step, "Dose step", 1e-6, 1e5)
  delay <- ddi_numeric(delay, "Delay (days)", 0, min(180, config$horizon - config$decision - 1e-8))
  intervals <- sort(unique(suppressWarnings(as.numeric(unlist(grid$intervals)))))
  if (!length(intervals) || any(!is.finite(intervals) | intervals < max(.25, config$infusion / 24) | intervals > 180)) stop("Choose valid cycle intervals (days), not shorter than the infusion.")
  count <- floor((maximum - minimum) / step + 1e-9) + 1
  if (!is.finite(count) || count * length(intervals) > 24) stop("Oncology grid: maximum 24 candidates.")
  candidates <- expand.grid(dose = minimum + seq.int(0, count - 1) * step, interval = intervals)
  candidates$scenario <- paste0("candidate_", seq_len(nrow(candidates)))
  regimens <- rbind(data.frame(dose = config$dose, interval = config$interval, scenario = "maintain"), candidates)
  schedules <- setNames(lapply(seq_len(nrow(regimens)), function(i) {
    candidate <- config; candidate$dose <- regimens$dose[i]
    onco_schedule(candidate, delay = if (i == 1) 0 else delay, interval = regimens$interval[i])
  }), regimens$scenario)
  curves <- lapply(schedules, function(doses) onco_simulate(config, doses, pk_context = pk_context))
  metrics <- cbind(regimens, onco_metrics(config, curves, schedules)[, -1, drop = FALSE])
  metrics$targets_met <- metrics$tumor_goal_met & (!config$toxicity | metrics$anc_goal_met)
  curves$untreated <- onco_simulate(config, config$history[FALSE, ], exposure = data.frame(time = c(0, config$horizon), concentration = 0))
  pk_metadata <- if (is.null(pk_context)) NULL else c(pk_context[intersect(names(pk_context), c("label", "route", "adm_cmt", "time_unit", "concentration", "concentration_scale"))], list(parameters = as.list(mrgsolve::param(pk_context$model))))
  list(config = config, grid = grid, delay = delay, regimens = regimens, curves = curves, schedules = schedules, metrics = metrics, pk = pk_metadata)
}

onco_grid_selection <- function(result, selected) {
  if (!selected %in% result$regimens$scenario || selected == "maintain") stop("Choose a candidate regimen.")
  pair <- result
  pair$curves <- setNames(result$curves[c("maintain", selected, "untreated")], c("maintain", "change", "untreated"))
  pair$schedules <- setNames(result$schedules[c("maintain", selected)], c("maintain", "change"))
  pair$metrics <- result$metrics[match(c("maintain", selected), result$metrics$scenario), setdiff(names(result$metrics), c("dose", "interval", "targets_met"))]
  pair$metrics$scenario <- c("maintain", "change")
  pair$change <- result$regimens[result$regimens$scenario == selected, ]
  pair
}

onco_model_code <- function(config) {
  config <- onco_config(config)
  external <- config$free_pk
  growth <- switch(config$growth, exponential = "KG", logistic = "KG*(1.0-TUMOR/CAP)", gompertz = "KG*log(CAP/fmax(1e-12,TUMOR))")
  paste(c("// Generic oncology research model; not a validated drug-specific dosing tool.",
    if (external) "// TIME IN DAYS. CP mg/L, tumor SLD mm, ANC 10^9/L. PK units are defined in the separate PK code." else "// TIME IN DAYS. Doses mg, concentrations mg/L, CL L/day, V L, tumor SLD mm, ANC 10^9/L.",
    if (external) "// PD ONLY: supply external CP in PD concentration units; the R export couples the selected PK." else "// Infusion RATE = mg/day (convert infusion hours to days). Initial exposure is zero.",
    "$PARAM", paste(paste0(names(config$parameters)[if (external) !names(config$parameters) %in% c("V", "CL") else TRUE], "=", unlist(config$parameters)[if (external) !names(config$parameters) %in% c("V", "CL") else TRUE]), collapse = ", "),
    if (external) c("$PARAM CP=0", "$CMT TUMOR") else "$CMT CENT TUMOR", if (config$toxicity) "$CMT PROL TR1 TR2 TR3 ANC",
    "$MAIN", "TUMOR_0=T0;", if (config$toxicity) "PROL_0=ANC0; TR1_0=ANC0; TR2_0=ANC0; TR3_0=ANC0; ANC_0=ANC0;",
    "$ODE", if (external) "double C=fmax(0.0,CP);" else c("double C=fmax(0.0,CENT/V);", "dxdt_CENT=-CL/V*CENT;"),
    paste0("dxdt_TUMOR=(", growth, "-KILL*C/(EC50+C)*exp(-RES*SOLVERTIME))*fmax(1e-12,TUMOR);"),
    if (config$toxicity) c("double ktr=4.0/MTT;", "double drug=fmin(1.0,SLOPE*C);", "double feedback=pow(ANC0/fmax(1e-10,ANC),GAMMA);",
      "dxdt_PROL=ktr*fmax(1e-10,PROL)*((1.0-drug)*feedback-1.0);",
      "dxdt_TR1=ktr*(fmax(1e-10,PROL)-fmax(1e-10,TR1));", "dxdt_TR2=ktr*(fmax(1e-10,TR1)-fmax(1e-10,TR2));",
      "dxdt_TR3=ktr*(fmax(1e-10,TR2)-fmax(1e-10,TR3));", "dxdt_ANC=ktr*(fmax(1e-10,TR3)-fmax(1e-10,ANC));"),
    if (!external) "$TABLE double CP=CENT/V;", "$CAPTURE CP"), collapse = "\n")
}

onco_export_script <- function(config, fraction = 0.75, delay = 0, interval = config$interval, pk_context = NULL, grid = NULL) {
  dump <- function(x) paste(capture.output(dput(x)), collapse = "\n")
  helpers <- c("ddi_numeric", "pd_pk_profile", "onco_config", "onco_schedule", "onco_exposure", "onco_rhs", "onco_simulate", "onco_metrics", "onco_compare", "onco_compare_grid")
  paste(c("# Exploratory oncology simulation. Requires deSolve; time in days, dose in mg.",
    "# Selected PK + tumor growth/inhibition + modified Friberg structure.",
    "# Not a validated anticancer-drug model. Contains the explicitly exported dose history.",
    "`%||%` <- function(x,y) if (is.null(x) || !length(x)) y else x", paste0("ONCO_DEFAULTS <- ", dump(ONCO_DEFAULTS)),
    vapply(helpers, function(name) paste0(name, " <- ", paste(deparse(get(name, mode = "function")), collapse = "\n")), character(1)),
    pd_pk_export_setup(pk_context), paste0("config <- ", dump(config)),
    if (is.null(grid)) paste0("result <- onco_compare(config, fraction=", fraction, ", delay=", delay, ", interval=", interval, ", pk_context=pk_context)") else
      paste0("result <- onco_compare_grid(config, grid=", dump(grid), ", delay=", delay, ", pk_context=pk_context)"),
    "print(result$metrics)", "with(result$curves$maintain, plot(time,TUMOR,type='l',lty=2,xlab='Day',ylab='Tumor size'))",
    if (is.null(grid)) "with(result$curves$change, lines(time,TUMOR,col='#196f76'))" else
      "for (name in setdiff(names(result$curves),c('maintain','untreated'))) with(result$curves[[name]], lines(time,TUMOR,col='#196f76'))"), collapse = "\n\n")
}
