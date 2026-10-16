PMETRICS_ARTIFACT_ROOT <- file.path(APP_ROOT, "pmetrics", "artifacts")
PMETRICS_ARTIFACT_VERSION <- 1L
PMETRICS_STEADY_STATE_TOLERANCE <- 1e-8

pmetrics_scalar <- function(value, field, type = c("character", "numeric", "logical")) {
  type <- match.arg(type)
  if (is.null(value) || length(value) != 1L) stop("Pmetrics artifact field `", field, "` must be scalar.")
  output <- switch(type,
    character = as.character(value),
    numeric = suppressWarnings(as.numeric(value)),
    logical = as.logical(value)
  )
  if (type == "character" && !nzchar(output)) stop("Pmetrics artifact field `", field, "` cannot be empty.")
  if (type == "numeric" && !is.finite(output)) stop("Pmetrics artifact field `", field, "` must be finite.")
  if (type == "logical" && is.na(output)) stop("Pmetrics artifact field `", field, "` must be boolean.")
  output
}

pmetrics_artifact_files <- function() {
  if (!dir.exists(PMETRICS_ARTIFACT_ROOT)) return(character())
  sort(list.files(PMETRICS_ARTIFACT_ROOT, pattern = "\\.json$", full.names = TRUE))
}

pmetrics_artifact_path <- function(id) {
  id <- pmetrics_scalar(id, "id")
  if (!grepl("^[a-z0-9][a-z0-9_-]{0,79}$", id)) stop("Invalid Pmetrics artifact id.")
  path <- normalizePath(file.path(PMETRICS_ARTIFACT_ROOT, paste0(id, ".json")), winslash = "/", mustWork = TRUE)
  root <- paste0(normalizePath(PMETRICS_ARTIFACT_ROOT, winslash = "/", mustWork = TRUE), "/")
  if (!startsWith(path, root)) stop("Pmetrics artifact path escapes the trusted artifact directory.")
  path
}

pmetrics_parameter_table <- function(artifact) {
  rows <- artifact$structure$parameters
  if (!is.list(rows) || !length(rows)) stop("Pmetrics artifact requires structural parameters.")
  output <- do.call(rbind, lapply(seq_along(rows), function(index) {
    row <- rows[[index]]
    data.frame(
      name = pmetrics_scalar(row$name, paste0("structure.parameters[", index, "].name")),
      lower = pmetrics_scalar(row$lower, paste0("structure.parameters[", index, "].lower"), "numeric"),
      upper = pmetrics_scalar(row$upper, paste0("structure.parameters[", index, "].upper"), "numeric"),
      stringsAsFactors = FALSE
    )
  }))
  if (anyDuplicated(output$name)) stop("Pmetrics structural parameter names must be unique.")
  if (any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", output$name))) stop("Invalid Pmetrics structural parameter name.")
  if (any(output$upper <= output$lower)) stop("Each Pmetrics parameter upper bound must exceed its lower bound.")
  output
}

pmetrics_support_table <- function(artifact, parameters = pmetrics_parameter_table(artifact)) {
  rows <- artifact$supportPoints
  if (!is.list(rows) || !length(rows)) stop("Pmetrics artifact requires weighted support points.")
  required <- c(parameters$name, "prob")
  output <- do.call(rbind, lapply(seq_along(rows), function(index) {
    row <- rows[[index]]
    missing <- setdiff(required, names(row))
    extra <- setdiff(names(row), required)
    if (length(missing)) stop("Support point ", index, " is missing: ", paste(missing, collapse = ", "), ".")
    if (length(extra)) stop("Support point ", index, " has unknown fields: ", paste(extra, collapse = ", "), ".")
    values <- vapply(required, function(name) pmetrics_scalar(row[[name]], paste0("supportPoints[", index, "].", name), "numeric"), numeric(1))
    as.data.frame(as.list(values), check.names = FALSE)
  }))
  output[] <- lapply(output, as.numeric)
  if (any(output$prob <= 0) || sum(output$prob) <= 0) stop("Pmetrics support-point probabilities must be positive and sum to a positive value.")
  for (index in seq_len(nrow(parameters))) {
    name <- parameters$name[[index]]
    if (any(output[[name]] < parameters$lower[[index]] | output[[name]] > parameters$upper[[index]])) {
      stop("Support-point values for `", name, "` must remain inside its search bounds.")
    }
  }
  output$prob <- output$prob / sum(output$prob)
  output
}

pmetrics_error_spec <- function(artifact, residual_error_mode = "model", fixed_residual_cv = 1) {
  if (identical(residual_error_mode, "fixed_cv")) {
    cv <- pmetrics_scalar(fixed_residual_cv, "fixed_residual_cv", "numeric") / 100
    if (cv <= 0 || cv > 1) stop("The fixed proportional residual error must be between 0 and 100%.")
    return(list(type = "proportional", initial = 1, coefficients = c(0, cv, 0, 0), fixed = TRUE))
  }
  if (!identical(residual_error_mode, "model")) stop("Unknown residual error setting: ", residual_error_mode, ".")
  error <- artifact$error
  type <- pmetrics_scalar(error$type, "error.type")
  if (!type %in% c("proportional", "additive")) stop("Pmetrics error type must be `proportional` or `additive`.")
  coefficients <- suppressWarnings(as.numeric(unlist(error$coefficients, use.names = FALSE)))
  if (length(coefficients) != 4L || any(!is.finite(coefficients))) {
    stop("Pmetrics error coefficients must contain four finite values.")
  }
  initial <- pmetrics_scalar(error$initial %||% 1, "error.initial", "numeric")
  if (initial <= 0) stop("Pmetrics error initial value must be positive.")
  list(
    type = type,
    initial = initial,
    coefficients = coefficients,
    fixed = pmetrics_scalar(error$fixed %||% TRUE, "error.fixed", "logical")
  )
}

validate_pmetrics_artifact <- function(artifact) {
  if (!is.list(artifact)) stop("Pmetrics artifact must be a JSON object.")
  if (pmetrics_scalar(artifact$schemaVersion, "schemaVersion", "numeric") != PMETRICS_ARTIFACT_VERSION) {
    stop("Unsupported Pmetrics artifact schema version.")
  }
  id <- pmetrics_scalar(artifact$id, "id")
  if (!grepl("^[a-z0-9][a-z0-9_-]{0,79}$", id)) stop("Invalid Pmetrics artifact id.")
  if (!identical(tolower(pmetrics_scalar(artifact$engine, "engine")), "pmetrics")) stop("Artifact engine must be Pmetrics.")
  template <- pmetrics_scalar(artifact$structure$template, "structure.template")
  supported <- list(
    one_compartment_iv = c("ke", "v"),
    one_compartment_oral = c("ka", "ke", "v"),
    two_compartment_iv = c("ke", "k12", "k21", "v")
  )
  if (!template %in% names(supported)) stop("Unsupported Pmetrics structural template: ", template, ".")
  parameters <- pmetrics_parameter_table(artifact)
  if (!setequal(parameters$name, supported[[template]])) {
    stop("Template `", template, "` requires parameters: ", paste(supported[[template]], collapse = ", "), ".")
  }
  pmetrics_error_spec(artifact)
  support <- pmetrics_support_table(artifact, parameters)
  if (any(as.matrix(support[, parameters$name, drop = FALSE]) <= 0)) stop("Pmetrics structural support-point parameters must be positive.")
  if (!isTRUE(pmetrics_error_spec(artifact)$fixed)) stop("Pmetrics patient artifacts require a fixed, previously estimated error model.")
  route <- pmetrics_scalar(artifact$administration$route, "administration.route")
  mode <- pmetrics_scalar(artifact$administration$mode, "administration.mode")
  valid_administration <- switch(template,
    one_compartment_oral = identical(route, "Oral") && identical(mode, "ORAL"),
    one_compartment_iv = identical(route, "IV") && mode %in% c("IV_INTERMITTENT", "IV_CONTINUOUS"),
    two_compartment_iv = identical(route, "IV") && mode %in% c("IV_INTERMITTENT", "IV_CONTINUOUS"),
    FALSE
  )
  if (!valid_administration) stop("Pmetrics artifact route and administration mode do not match its structural template.")
  artifact$schemaVersion <- PMETRICS_ARTIFACT_VERSION
  artifact$parameters <- parameters
  artifact$support <- support
  artifact$route <- route
  artifact$mode <- mode
  artifact
}

read_pmetrics_artifact <- function(id) {
  artifact <- jsonlite::fromJSON(pmetrics_artifact_path(id), simplifyVector = FALSE)
  artifact <- validate_pmetrics_artifact(artifact)
  if (!identical(artifact$id, id)) stop("Pmetrics artifact filename and id do not match.")
  artifact
}

pmetrics_artifact_catalog <- function() {
  files <- pmetrics_artifact_files()
  if (!length(files)) return(data.frame(id = character(), label = character(), labelEn = character()))
  rows <- lapply(files, function(path) {
    id <- tools::file_path_sans_ext(basename(path))
    artifact <- read_pmetrics_artifact(id)
    data.frame(id = id, label = pmetrics_scalar(artifact$label, "label"),
      labelEn = pmetrics_scalar(artifact$labelEn %||% artifact$label, "labelEn"), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

build_pmetrics_model <- function(artifact, error = pmetrics_error_spec(artifact)) {
  if (!requireNamespace("Pmetrics", quietly = TRUE)) stop("The Pmetrics package is not installed on this R server.")
  bounds <- stats::setNames(lapply(seq_len(nrow(artifact$parameters)), function(index) {
    Pmetrics::ab(artifact$parameters$lower[[index]], artifact$parameters$upper[[index]])
  }), artifact$parameters$name)
  error_model <- switch(error$type,
    proportional = Pmetrics::proportional(error$initial, error$coefficients, fixed = error$fixed),
    additive = Pmetrics::additive(error$initial, error$coefficients, fixed = error$fixed)
  )
  switch(artifact$structure$template,
    one_compartment_iv = Pmetrics::PM_model$new(
      pri = bounds,
      eqn = function() {
        dx[1] <- -ke * x[1] + rateiv[1] + b[1]
      },
      out = function() {
        y[1] <- x[1] / v
      },
      err = list(error_model)
    ),
    one_compartment_oral = Pmetrics::PM_model$new(
      pri = bounds,
      eqn = function() {
        dx[1] <- -ka * x[1] + b[1]
        dx[2] <- ka * x[1] - ke * x[2]
      },
      out = function() {
        y[1] <- x[2] / v
      },
      err = list(error_model)
    ),
    two_compartment_iv = Pmetrics::PM_model$new(
      pri = bounds,
      eqn = function() {
        dx[1] <- -(ke + k12) * x[1] + k21 * x[2] + rateiv[1] + b[1]
        dx[2] <- k12 * x[1] - k21 * x[2]
      },
      out = function() {
        y[1] <- x[1] / v
      },
      err = list(error_model)
    ),
    stop("Unsupported Pmetrics structural template.")
  )
}

pmetrics_bridge_code <- function(artifact) {
  switch(artifact$structure$template,
    one_compartment_iv = paste(
      "$PARAM ke = 0.15, v = 30",
      "$CMT @annotated",
      "CENT : central compartment [ADM, OBS]",
      "$ODE",
      "dxdt_CENT = -ke * CENT;",
      "$TABLE",
      "double CONC = CENT / v;",
      "$CAPTURE CONC",
      sep = "\n"
    ),
    one_compartment_oral = paste(
      "$PARAM ka = 1, ke = 0.15, v = 30",
      "$CMT @annotated",
      "DEPOT : oral depot [ADM]",
      "CENT : central compartment [OBS]",
      "$ODE",
      "dxdt_DEPOT = -ka * DEPOT;",
      "dxdt_CENT = ka * DEPOT - ke * CENT;",
      "$TABLE",
      "double CONC = CENT / v;",
      "$CAPTURE CONC",
      sep = "\n"
    ),
    two_compartment_iv = paste(
      "$PARAM ke = 0.15, k12 = 0.3, k21 = 0.2, v = 30",
      "$CMT @annotated",
      "CENT : central compartment [ADM, OBS]",
      "PERIPH : peripheral compartment",
      "$ODE",
      "dxdt_CENT = -(ke + k12) * CENT + k21 * PERIPH;",
      "dxdt_PERIPH = k12 * CENT - k21 * PERIPH;",
      "$TABLE",
      "double CONC = CENT / v;",
      "$CAPTURE CONC",
      sep = "\n"
    ),
    stop("Unsupported Pmetrics structural template.")
  )
}

pmetrics_compartment_contract <- function(artifact) {
  switch(artifact$structure$template,
    one_compartment_oral = list(adm_cmt = 1L, obs_cmt = 2L),
    one_compartment_iv = list(adm_cmt = 1L, obs_cmt = 1L),
    two_compartment_iv = list(adm_cmt = 1L, obs_cmt = 1L),
    stop("Unsupported Pmetrics structural template.")
  )
}

pmetrics_terminal_rate <- function(artifact) {
  support <- artifact$support[artifact$support$prob > 0, , drop = FALSE]
  rate <- switch(artifact$structure$template,
    one_compartment_iv = support$ke,
    one_compartment_oral = pmin(support$ka, support$ke),
    two_compartment_iv = {
      total <- support$ke + support$k12 + support$k21
      (total - sqrt(pmax(0, total^2 - 4 * support$ke * support$k21))) / 2
    },
    stop("Unsupported Pmetrics structural template.")
  )
  minimum <- min(rate)
  if (!is.finite(minimum) || minimum <= 0) stop("Pmetrics support points imply an invalid terminal rate.")
  minimum
}

compile_pmetrics_bridge <- function(artifact, soloc, cache = NULL) {
  if (is.null(soloc) || !dir.exists(soloc)) stop("A session-scoped compilation directory is required for Pmetrics models.")
  code <- pmetrics_bridge_code(artifact)
  key <- safe_model_id(paste0("pmetrics_", artifact$id, "_", short_hash(code)))
  if (!is.null(cache) && !is.null(cache[[key]])) return(cache[[key]])
  model <- mrgsolve::mcode(key, code, soloc = soloc, quiet = TRUE)
  if (!is.null(cache)) cache[[key]] <- model
  model
}

expand_pmetrics_doses <- function(doses, end_time = max(doses$time), min_ke = 0.08) {
  if (!nrow(doses)) stop("At least one administered dose is required.")
  doses <- expand_steady_state_doses(doses, end_time)
  output <- lapply(seq_len(nrow(doses)), function(index) {
    dose <- doses[index, , drop = FALSE]
    count <- max(1L, as.integer(dose$count[[1]] %||% 1L))
    interval <- as.numeric(dose$interval[[1]] %||% 0)
    steady_state <- as.integer(dose$ss[[1]] %||% 0L) == 1L
    offsets <- if (steady_state) {
      if (!is.finite(interval) || interval <= 0) stop("Steady state requires a positive dosing interval.")
      warmup <- ceiling(-log(PMETRICS_STEADY_STATE_TOLERANCE) / (min_ke * interval))
      if (!is.finite(warmup) || warmup > 10000) stop("Pmetrics steady-state warmup exceeds 10000 doses.")
      seq.int(-warmup, 0L) * interval
    } else if (count > 1L) {
      seq.int(0L, count - 1L) * interval
    } else {
      0
    }
    expanded <- dose[rep(1L, length(offsets)), , drop = FALSE]
    expanded$time <- as.numeric(dose$time[[1]]) + offsets
    expanded
  })
  doses <- do.call(rbind, output)
  doses[order(doses$time), , drop = FALSE]
}

build_pmetrics_data <- function(doses, observations, min_ke = 0.08) {
  doses <- expand_pmetrics_doses(doses, max(c(doses$time, observations$time)), min_ke)
  shift <- max(0, -min(c(doses$time, observations$time %||% numeric()), na.rm = TRUE))
  dose_rows <- data.frame(
    id = 1L, evid = 1L, time = as.numeric(doses$time) + shift,
    dur = pmax(0, as.numeric(doses$infusion)), dose = as.numeric(doses$amount),
    addl = 0L, ii = 0, input = 1L, out = NA_real_, outeq = NA_integer_
  )
  observation_rows <- data.frame(
    id = rep(1L, nrow(observations)), evid = rep(0L, nrow(observations)),
    time = as.numeric(observations$time) + shift, dur = NA_real_, dose = NA_real_,
    addl = NA_integer_, ii = NA_real_, input = NA_integer_,
    out = as.numeric(observations$concentration), outeq = rep(1L, nrow(observations))
  )
  data <- rbind(dose_rows, observation_rows)
  data <- data[order(data$time, -data$evid), , drop = FALSE]
  rownames(data) <- NULL
  list(data = data, time_shift = shift)
}

pmetrics_posterior_summary <- function(points, parameter_names) {
  if (!is.data.frame(points) || !nrow(points)) stop("Pmetrics returned no posterior support points.")
  missing <- setdiff(c(parameter_names, "prob"), names(points))
  if (length(missing)) stop("Pmetrics posterior is missing: ", paste(missing, collapse = ", "), ".")
  points <- as.data.frame(points, check.names = FALSE)
  points$prob <- suppressWarnings(as.numeric(points$prob))
  if (any(!is.finite(points$prob) | points$prob < 0) || sum(points$prob) <= 0) stop("Invalid Pmetrics posterior probabilities.")
  points$prob <- points$prob / sum(points$prob)
  values <- as.matrix(points[, parameter_names, drop = FALSE])
  storage.mode(values) <- "double"
  if (any(!is.finite(values))) stop("Invalid Pmetrics posterior parameter values.")
  mean <- colSums(values * points$prob)
  mode <- values[which.max(points$prob), ]
  list(points = points, mean = mean, mode = mode)
}

run_pmetrics_map <- function(artifact, doses, observations, residual_error_mode = "model", fixed_residual_cv = 1) {
  if (!nrow(observations)) {
    summary <- pmetrics_posterior_summary(artifact$support, artifact$parameters$name)
    return(list(mapped = FALSE, posterior = summary, warnings = character(), data = NULL, time_shift = 0))
  }
  error <- pmetrics_error_spec(artifact, residual_error_mode, fixed_residual_cv)
  model <- build_pmetrics_model(artifact, error)
  input <- build_pmetrics_data(doses, observations, pmetrics_terminal_rate(artifact))
  data <- Pmetrics::PM_data$new(data = input$data, quiet = TRUE)
  run_dir <- tempfile(paste0("pmetrics-map-", artifact$id, "-"))
  dir.create(run_dir, recursive = TRUE)
  on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
  prior_path <- file.path(run_dir, "prior.csv")
  utils::write.csv(artifact$support, prior_path, row.names = FALSE)
  warnings <- character()
  result <- withCallingHandlers(
    suppressMessages(model$map(
      data = data, path = run_dir, run = 1L, prior = prior_path,
      idelta = 0.05, report = "none", quiet = TRUE, overwrite = TRUE
    )),
    warning = function(warning) {
      warnings <<- c(warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  posterior <- pmetrics_posterior_summary(result$final$postPoints, artifact$parameters$name)
  list(mapped = TRUE, posterior = posterior, warnings = unique(warnings), data = input$data, time_shift = input$time_shift)
}

fit_one_pmetrics_model <- function(
  specification, doses, observations, covariates, covariate_history = NULL,
  custom_soloc = NULL, custom_cache = NULL,
  residual_error_mode = "model", fixed_residual_cv = 1
) {
  artifact <- read_pmetrics_artifact(specification$artifact_id %||% specification$id)
  if (length(covariates) || (!is.null(covariate_history) && ncol(covariate_history) > 1L)) {
    stop("This Pmetrics artifact does not define covariates.")
  }
  if (identical(artifact$route, "Oral") && any(as.numeric(doses$infusion) > 0)) {
    stop("The one-compartment oral Pmetrics template does not accept infusion durations.")
  }
  model <- compile_pmetrics_bridge(artifact, custom_soloc, custom_cache)
  mapped <- run_pmetrics_map(artifact, doses, observations, residual_error_mode, fixed_residual_cv)
  posterior_parameters <- mapped$posterior$mean
  model <- safe_param(model, as.list(posterior_parameters))
  compartments <- pmetrics_compartment_contract(artifact)
  contract <- list(ok = TRUE, errors = character(), warnings = mapped$warnings,
    adm_cmt = compartments$adm_cmt, obs_cmt = compartments$obs_cmt, route = artifact$route, mode = artifact$mode,
    n_eta = 0L, n_sigma = 0L, engine = "pmetrics")
  map_data <- build_map_data(doses, observations, contract$adm_cmt, contract$obs_cmt, list())
  list(
    id = artifact$id,
    label = specification$label %||% artifact$label,
    engine = "pmetrics",
    model = model,
    estimate = NULL,
    data = map_data,
    contract = contract,
    split_lego = FALSE,
    current_covariates = list(),
    source_doses = doses,
    source_observations = observations,
    source_covariates = list(),
    source_covariate_history = NULL,
    posterior_parameters = posterior_parameters,
    posterior_mode = mapped$posterior$mode,
    posterior_points = mapped$posterior$points,
    pmetrics = list(mapped = mapped$mapped, artifact_id = artifact$id,
      time_shift = mapped$time_shift, warnings = mapped$warnings),
    residual_error = list(
      mode = if (identical(residual_error_mode, "fixed_cv")) "fixed_cv" else "pmetrics_model",
      fixed_cv = if (identical(residual_error_mode, "fixed_cv")) fixed_residual_cv else NA_real_,
      pmetrics = pmetrics_error_spec(artifact, residual_error_mode, fixed_residual_cv)
    )
  )
}

pmetrics_residual_sd <- function(prediction, error) {
  coefficients <- as.numeric(error$coefficients)
  assay <- coefficients[[1]] + coefficients[[2]] * prediction +
    coefficients[[3]] * prediction^2 + coefficients[[4]] * prediction^3
  # Match pharmsol's AssayErrorModel, including the additive variance component.
  sd <- if (identical(error$type, "additive")) sqrt(assay^2 + error$initial^2) else error$initial * assay
  if (any(!is.finite(sd) | sd < 0)) stop("Invalid Pmetrics error SD over the simulated concentration range.")
  sd
}

pmetrics_residual_values <- function(prediction, error) {
  pmax(0, prediction + stats::rnorm(length(prediction), 0, pmetrics_residual_sd(prediction, error)))
}
