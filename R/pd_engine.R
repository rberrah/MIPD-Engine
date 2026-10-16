PD_TYPES <- c("linear", "emax", "hill", "inhibit_in", "stimulate_in", "inhibit_out", "stimulate_out")
PD_DEFAULTS <- c(E0 = 100, SLOPE = 2, EMAX = 0.8, EC50 = 2, HILL = 1.5, KE0 = 0.5, KOUT = 0.15)

pd_parameters <- function(type, delay = FALSE) {
  if (!type %in% PD_TYPES) stop("Unknown PD model.")
  c("E0", if (type == "linear") "SLOPE" else c("EMAX", "EC50"),
    if (type == "hill") "HILL", if (delay) "KE0", if (grepl("_(in|out)$", type)) "KOUT")
}

pd_config <- function(config) {
  config$type <- as.character(config$type %||% "emax")
  config$delay <- isTRUE(config$delay)
  names <- pd_parameters(config$type, config$delay)
  p <- PD_DEFAULTS
  for (name in names) p[[name]] <- ddi_numeric(config$parameters[[name]] %||% p[[name]], name,
    if (name %in% c("E0", "EC50", "HILL", "KE0", "KOUT")) 1e-8 else if (name == "EMAX" && grepl("^inhibit", config$type)) 0 else -1e6,
    if (name == "EMAX" && grepl("^inhibit", config$type)) 1 else 1e6)
  if (grepl("^stimulate", config$type) && p[["EMAX"]] < 0) stop("Stimulation requires EMAX >= 0.")
  config$parameters <- p
  config$horizon <- ddi_numeric(config$horizon %||% 24, "Horizon", 0.01, 2400)
  config$c0 <- ddi_numeric(config$c0 %||% 10, "Initial concentration", 0, 1e6)
  config$kel <- ddi_numeric(config$kel %||% 0.1, "Exposure elimination rate", 0, 100)
  config$exposure <- as.character(config$exposure %||% "exponential")
  if (!config$exposure %in% c("exponential", "table", "pk")) stop("Unknown exposure source.")
  if (config$exposure %in% c("table", "pk")) {
    data <- config$profile
    if (!is.data.frame(data) || !all(c("time", "concentration") %in% names(data)) || nrow(data) < 2 || nrow(data) > if (config$exposure == "pk") 200000 else 2000) stop("Invalid exposure profile row count.")
    if (any(!is.finite(data$time) | !is.finite(data$concentration) | data$time < 0 | data$concentration < 0) || anyDuplicated(data$time)) stop("Invalid or duplicate exposure times.")
    config$profile <- data[order(data$time), c("time", "concentration")]
    if (min(data$time) != 0 || max(data$time) < config$horizon) stop("The exposure profile must cover time 0 through the simulation horizon.")
  }
  config
}

pd_exposure <- function(time, config) {
  if (config$exposure %in% c("table", "pk")) stats::approx(config$profile$time, config$profile$concentration, xout = time, rule = 2)$y
  else config$c0 * exp(-config$kel * time)
}

pd_signal <- function(concentration, config) {
  p <- config$parameters
  c <- pmax(0, concentration)
  if (config$type == "linear") return(p[["SLOPE"]] * c)
  h <- if (config$type == "hill") p[["HILL"]] else 1
  # Logistic form avoids overflow in C^h at extreme concentrations.
  p[["EMAX"]] * stats::plogis(h * (log(c) - log(p[["EC50"]])))
}

pd_simulate <- function(config, times = NULL) {
  config <- pd_config(config)
  p <- config$parameters
  if (is.null(times)) times <- seq(0, config$horizon, length.out = 501)
  if (!length(times) || any(!is.finite(times) | times < 0 | times > config$horizon)) stop("Invalid PD output times.")
  grid <- sort(unique(round(c(0, times), 10)))
  indirect <- grepl("_(in|out)$", config$type)
  if (config$delay || indirect) {
    knots <- if (config$exposure %in% c("table", "pk")) unique(round(config$profile$time[config$profile$time > 0 & config$profile$time < max(grid)], 10)) else numeric()
    grid <- sort(unique(c(grid, knots)))
    # Restart at exposure knots, without changing states, so a sparse observation grid cannot skip a short pulse.
    events <- if (length(knots)) list(func = function(time, state, parameters) state, time = knots) else NULL
    rhs <- function(time, state, parameters) {
      cp <- pd_exposure(time, config)
      signal <- pd_signal(if (config$delay) state[[1]] else cp, config)
      production <- switch(config$type, inhibit_in = 1 - signal, stimulate_in = 1 + signal, 1)
      loss <- switch(config$type, inhibit_out = 1 - signal, stimulate_out = 1 + signal, 1)
      list(c(p[["KE0"]] * (cp - state[[1]]), p[["KOUT"]] * (p[["E0"]] * production - loss * state[[2]])))
    }
    if (length(grid) == 1) states <- data.frame(time = 0, CE = 0, RESPONSE = p[["E0"]])
    else states <- as.data.frame(deSolve::ode(y = c(CE = 0, RESPONSE = p[["E0"]]), times = grid,
      func = rhs, parms = NULL, method = "lsoda", events = events, rtol = 1e-8, atol = 1e-9))
    ce <- if (config$delay) states$CE else pd_exposure(grid, config)
    effect <- if (indirect) states$RESPONSE else p[["E0"]] + pd_signal(ce, config)
  } else {
    ce <- pd_exposure(grid, config)
    effect <- p[["E0"]] + pd_signal(ce, config)
  }
  output <- data.frame(time = grid, concentration = pd_exposure(grid, config), driver = ce, effect = effect, baseline = p[["E0"]])
  if (any(!is.finite(as.matrix(output)))) stop("Non-finite PD prediction.")
  output <- output[match(round(times, 10), output$time), , drop = FALSE]
  output$time <- times
  output
}

pd_read_table <- function(text, required) {
  if (!is.character(text) || length(text) != 1 || nchar(text, type = "bytes") > 200000) stop("Table is empty or too large.")
  data <- utils::read.csv(text = text, check.names = FALSE, na.strings = c("", "NA"))
  if (!all(required %in% names(data)) || !nrow(data) || nrow(data) > 2000) stop("Expected CSV columns: ", paste(required, collapse = ", "), "; 1-2000 rows.")
  data <- data[, required, drop = FALSE]
  for (name in required) data[[name]] <- suppressWarnings(as.numeric(data[[name]]))
  if (any(!is.finite(as.matrix(data)))) stop("Every supplied value must be numeric and finite.")
  data
}

pd_predict <- function(config, data, mode) {
  config <- pd_config(config)
  if (!mode %in% c("time", "concentration")) stop("Unknown observation mode.")
  if (mode == "concentration") {
    if (config$delay || grepl("_(in|out)$", config$type)) stop("Delayed/indirect effects require timed data and an exposure history.")
    if (any(data$concentration < 0)) stop("Concentrations cannot be negative.")
    return(config$parameters[["E0"]] + pd_signal(data$concentration, config))
  }
  pd_simulate(config, times = data$time)$effect
}

pd_fit <- function(config, data, mode = "time", estimate = c("E0", "EMAX", "EC50")) {
  config <- pd_config(config)
  allowed <- pd_parameters(config$type, config$delay)
  if (!length(estimate) || anyDuplicated(estimate) || !all(estimate %in% allowed)) stop("Select active parameters to estimate.")
  required <- c(if (mode == "concentration") "concentration" else "time", "effect")
  if (!is.data.frame(data) || !all(required %in% names(data)) || nrow(data) <= length(estimate) + 1 || any(!is.finite(as.matrix(data[, required])))) stop("More observations than estimated parameters + 1 are required.")
  residual <- function(par) {
    candidate <- config
    candidate$parameters[estimate] <- par
    pd_predict(candidate, data, mode) - data$effect
  }
  lower <- ifelse(estimate %in% c("E0", "EC50", "HILL", "KE0", "KOUT"), 1e-8, -1e6)
  upper <- rep(1e6, length(estimate))
  if (grepl("^(inhibit|stimulate)", config$type)) lower[estimate == "EMAX"] <- 0
  if (grepl("^inhibit", config$type)) upper[estimate == "EMAX"] <- 1
  initial <- config$parameters[estimate]
  fit <- minpack.lm::nls.lm(par = initial, fn = residual, lower = lower, upper = upper,
    control = minpack.lm::nls.lm.control(maxiter = 150, maxfev = 3000))
  if (!fit$info %in% 1:4) stop("PD optimization did not converge: ", fit$message)
  fitted <- config
  fitted$parameters[estimate] <- fit$par
  prediction <- pd_predict(fitted, data, mode)
  sse <- sum((prediction - data$effect)^2)
  sensitivity <- vapply(seq_along(fit$par), function(i) {
    step <- max(1e-7, abs(fit$par[[i]]) * 1e-5)
    trial <- fit$par
    trial[[i]] <- min(upper[[i]], trial[[i]] + step)
    if (trial[[i]] == fit$par[[i]]) trial[[i]] <- max(lower[[i]], trial[[i]] - step)
    (residual(trial) - residual(fit$par)) / (trial[[i]] - fit$par[[i]])
  }, numeric(nrow(data)))
  ill_conditioned <- !is.finite(kappa(sensitivity)) || kappa(sensitivity) > 1e8
  list(config = fitted, initial = config, data = transform(data, prediction = prediction, residual = data$effect - prediction),
    mode = mode, parameters = data.frame(parameter = estimate, initial = as.numeric(initial), estimate = as.numeric(fit$par)),
    rmse = sqrt(mean((prediction - data$effect)^2)), sse = sse, ill_conditioned = ill_conditioned,
    boundary = any(abs(fit$par - lower) < 1e-7 | abs(fit$par - upper) < 1e-7))
}

pd_model_code <- function(config) {
  config <- pd_config(config)
  p <- config$parameters
  driver <- if (config$delay) "CE" else "CP"
  signal <- if (config$type == "linear") paste0("SLOPE*", driver) else paste0("EMAX/(1.0+pow(EC50/fmax(", driver, ",1e-300),", if (config$type == "hill") "HILL" else "1.0", "))")
  production <- switch(config$type, inhibit_in = "(1.0-SIGNAL)", stimulate_in = "(1.0+SIGNAL)", "1.0")
  loss <- switch(config$type, inhibit_out = "(1.0-SIGNAL)", stimulate_out = "(1.0+SIGNAL)", "1.0")
  paste(c("$PROB PD builder: research only; CP is an external concentration input.",
    "$PARAM", paste(paste(names(p), format(p, digits = 15, trim = TRUE), sep = " = "), collapse = ", "), "CP = 0",
    "$CMT CE RESPONSE", "$MAIN", "CE_0 = 0; RESPONSE_0 = E0;", "$ODE",
    "dxdt_CE = KE0*(CP-CE);", paste0("double SIGNAL = ", signal, ";"),
    paste0("dxdt_RESPONSE = KOUT*(E0*", production, "-", loss, "*RESPONSE);"),
    "$TABLE", paste0("double PRED = ", if (grepl("_(in|out)$", config$type)) "RESPONSE" else paste0("E0+", signal), ";"),
    "$CAPTURE PRED CP"), collapse = "\n")
}

pd_export_script <- function(config) {
  config <- pd_config(config)
  helpers <- c("ddi_numeric", "pd_parameters", "pd_config", "pd_exposure", "pd_signal", "pd_simulate")
  paste(c("# PD simulation; deSolve required. Includes the selected exposure profile, without identifiers.",
    "`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x",
    paste0("PD_TYPES <- ", paste(capture.output(dput(PD_TYPES)), collapse = "\n")),
    paste0("PD_DEFAULTS <- ", paste(capture.output(dput(PD_DEFAULTS)), collapse = "\n")),
    vapply(helpers, function(name) paste0(name, " <- ", paste(deparse(get(name, mode = "function")), collapse = "\n")), character(1)),
    paste0("config <- ", paste(capture.output(dput(config)), collapse = "\n")),
    "result <- pd_simulate(config)", "plot(result$time, result$effect, type='l', xlab='Time', ylab='Effect')"), collapse = "\n\n")
}
