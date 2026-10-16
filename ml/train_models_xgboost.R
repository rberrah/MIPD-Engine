suppressPackageStartupMessages({
  library(mrgsolve)
  library(xgboost)
  library(jsonlite)
})

`%||%` <- function(value, fallback) {
  if (is.null(value) || !length(value)) fallback else value
}

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(Sys.getenv("PMX_TRAIN_ROOT", file.path(dirname(file_argument), "..")), winslash = "/", mustWork = TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "engine.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "ml_engine.R"), local = TRUE)

args <- if (exists(".training_args")) .training_args else commandArgs(trailingOnly = TRUE)
option_value <- function(prefix, default = NULL) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (length(hit)) sub(paste0("^", prefix, "="), "", hit[[1]]) else default
}

smoke <- "--smoke" %in% args
publish_artifacts <- "--publish" %in% args || "--publish-research" %in% args
n_patients <- as.integer(option_value("--n", if (smoke) "30" else "1000"))
base_argument <- option_value("--base", "all")
drug_argument <- option_value("--drug", "all")
mode_argument <- toupper(option_value("--mode", "all"))
report_path <- option_value("--report", "")
seed <- as.integer(option_value("--seed", "20260906"))
full_library_publish <- publish_artifacts && identical(base_argument, "all") &&
  identical(tolower(drug_argument), "all") && identical(mode_argument, "ALL")

if (!is.finite(n_patients) || n_patients < 24L) stop("--n must be at least 24 simulated patients per model and mode.")
if (publish_artifacts && (smoke || n_patients < 1000L)) {
  stop("Artifact publication requires a non-smoke run with at least 1000 simulated patients per model and mode.")
}

config_path <- file.path(APP_ROOT, "ml", "training-regimens.json")
training_config <- jsonlite::fromJSON(config_path, simplifyVector = FALSE)
if (!identical(as.integer(training_config$version %||% 0L), 1L)) stop("Unsupported ML training regimen version.")

population_config_path <- file.path(APP_ROOT, "ml", "training-populations.json")
population_config <- jsonlite::fromJSON(population_config_path, simplifyVector = FALSE)
if (!identical(as.integer(population_config$version %||% 0L), 1L)) stop("Unsupported ML training population version.")

MODE_ROUTE <- c(ORAL = "Oral", IV_INTERMITTENT = "IV", IV_CONTINUOUS = "IV")
MODE_ID <- c(ORAL = "oral", IV_INTERMITTENT = "intermittent", IV_CONTINUOUS = "continuous")
VALID_MODES <- names(MODE_ROUTE)

catalog_scopes <- function() {
  rows <- lapply(seq_len(nrow(MODEL_CATALOG)), function(index) {
    record <- MODEL_CATALOG[index, , drop = FALSE]
    modes <- model_administration_modes(record)
    do.call(rbind, lapply(modes, function(mode) data.frame(
      model_id = record$id[[1]],
      drug_key = record$drugKey[[1]],
      drug = record$drug[[1]],
      mode = mode,
      route = unname(MODE_ROUTE[[mode]]),
      stringsAsFactors = FALSE
    )))
  })
  do.call(rbind, rows)
}

scopes <- catalog_scopes()
if (!identical(base_argument, "all")) {
  selected <- strsplit(base_argument, ",", fixed = TRUE)[[1]]
  scopes <- scopes[scopes$model_id %in% selected, , drop = FALSE]
  missing <- setdiff(selected, scopes$model_id)
  if (length(missing)) stop("Unknown or filtered model(s): ", paste(missing, collapse = ", "))
}
if (!identical(tolower(drug_argument), "all")) {
  selected_drugs <- tolower(strsplit(drug_argument, ",", fixed = TRUE)[[1]])
  scopes <- scopes[tolower(scopes$drug_key) %in% selected_drugs | tolower(scopes$drug) %in% selected_drugs, , drop = FALSE]
}
if (!identical(mode_argument, "ALL")) {
  selected_modes <- strsplit(mode_argument, ",", fixed = TRUE)[[1]]
  if (any(!selected_modes %in% VALID_MODES)) stop("--mode accepts: all, ", paste(VALID_MODES, collapse = ", "))
  scopes <- scopes[scopes$mode %in% selected_modes, , drop = FALSE]
}
if (!nrow(scopes)) stop("No model matches the requested training scope.")

continuous_scopes <- scopes$mode == "IV_CONTINUOUS"
if (any(continuous_scopes) && !"--paired-benchmark" %in% args) {
  message("Skipping ", sum(continuous_scopes), " continuous-infusion scope(s): C0 plus one hour after infusion end is not defined for an ongoing infusion.")
  scopes <- scopes[!continuous_scopes, , drop = FALSE]
}
if (!nrow(scopes)) stop("The requested scopes only contain continuous infusions, for which this sampling benchmark is not applicable.")

regimen_config <- function(scope) {
  profile <- training_config$profiles[[scope$drug_key[[1]]]][[scope$mode[[1]]]]
  if (is.null(profile)) {
    stop("Missing regimen configuration for ", scope$drug_key[[1]], " / ", scope$mode[[1]], ".")
  }
  profile
}

invisible(lapply(seq_len(nrow(scopes)), function(index) regimen_config(scopes[index, , drop = FALSE])))

OUTER_V <- if (smoke) 3L else 5L
OUTER_REPEATS <- if (smoke) 1L else 2L
MAX_RELATIVE_RMSE_PCT <- 15
MAX_ABSOLUTE_BIAS_PCT <- 10
MIN_WITHIN_20_PCT <- 80
DALEX_BACKGROUND_SIZE <- 200L
MIN_AUC_RATIO <- 0.05
MAX_AUC_RATIO <- 20
SAMPLING_PROTOCOL <- list(
  minimumObservations = 2L,
  steadyStateRequired = TRUE,
  strategy = "c0_plus_one_hour_post_infusion",
  description = "Pre-dose C0 plus one concentration one hour after infusion end; infusion end is time zero for oral and bolus administrations.",
  target = "AUC24",
  residualError = "Model SIGMA retained for sparse concentrations; OMEGA supplied through idata."
)
C0_SAMPLING_PROTOCOL <- list(
  minimumObservations = 1L,
  steadyStateRequired = TRUE,
  strategy = "c0_only",
  description = "One pre-dose steady-state concentration.",
  target = "AUC24",
  residualError = "Model SIGMA retained for the sparse concentration; OMEGA supplied through idata."
)
XGB_PARAMETERS <- list(
  max_depth = 4L,
  eta = 0.03,
  min_child_weight = 5,
  subsample = 0.85,
  colsample_bytree = 0.9,
  nrounds = if (smoke) 100L else 350L
)

set.seed(seed)
.compiled_models <- new.env(parent = emptyenv())
.covariate_defaults <- new.env(parent = emptyenv())

compiled_model <- function(model_id) {
  if (!exists(model_id, envir = .compiled_models, inherits = FALSE)) {
    assign(model_id, compile_model(model_id = model_id), envir = .compiled_models)
  }
  get(model_id, envir = .compiled_models, inherits = FALSE)
}

model_covariate_defaults <- function(model_id) {
  if (!exists(model_id, envir = .covariate_defaults, inherits = FALSE)) {
    parsed <- parse_covariates(read_library_code(model_id))
    values <- if (nrow(parsed)) stats::setNames(as.numeric(parsed$value), parsed$name) else numeric()
    assign(model_id, values, envir = .covariate_defaults)
  }
  get(model_id, envir = .covariate_defaults, inherits = FALSE)
}

population_profile <- function(model_id) {
  population_config$models[[model_id]] %||% list(covariates = list())
}

population_sampling_record <- function(model_id) {
  defaults <- model_covariate_defaults(model_id)
  profile <- population_profile(model_id)
  specs <- profile$covariates %||% list()
  varied <- intersect(names(defaults), names(specs)[vapply(specs, function(spec) !identical(spec$type %||% "", "fixed"), logical(1))])
  list(
    source = profile$source %||% "Model reference values",
    varied = varied,
    fixed = setdiff(names(defaults), varied)
  )
}

validate_population_config <- function() {
  unknown_models <- setdiff(names(population_config$models), MODEL_CATALOG$id)
  if (length(unknown_models)) stop("Unknown model(s) in training-populations.json: ", paste(unknown_models, collapse = ", "))
  allowed_types <- c("fixed", "uniform", "normal", "lognormal", "bernoulli", "cockcroft_creatinine", "schwartz_height", "bmi_height")
  for (model_id in names(population_config$models)) {
    specs <- population_profile(model_id)$covariates %||% list()
    unknown_covariates <- setdiff(names(specs), names(model_covariate_defaults(model_id)))
    if (length(unknown_covariates)) stop("Unknown covariate(s) for ", model_id, ": ", paste(unknown_covariates, collapse = ", "))
    invalid_types <- setdiff(vapply(specs, function(spec) spec$type %||% "", character(1)), allowed_types)
    if (length(invalid_types)) stop("Unsupported population distribution(s) for ", model_id, ": ", paste(invalid_types, collapse = ", "))
  }
  invisible(TRUE)
}

validate_population_config()

sample_truncated_normal <- function(n, mean, sd, min = -Inf, max = Inf) {
  values <- stats::rnorm(n, mean, sd)
  outside <- values < min | values > max
  while (any(outside)) {
    values[outside] <- stats::rnorm(sum(outside), mean, sd)
    outside <- values < min | values > max
  }
  values
}

sample_population_spec <- function(spec, n) {
  type <- spec$type %||% stop("Population covariate specification has no type.")
  if (identical(type, "fixed")) return(rep(as.numeric(spec$value), n))
  if (identical(type, "uniform")) return(stats::runif(n, as.numeric(spec$min), as.numeric(spec$max)))
  if (identical(type, "bernoulli")) return(stats::rbinom(n, 1, as.numeric(spec$probability)))
  if (identical(type, "normal")) {
    return(sample_truncated_normal(
      n, as.numeric(spec$mean), as.numeric(spec$sd),
      as.numeric(spec$min %||% -Inf), as.numeric(spec$max %||% Inf)
    ))
  }
  if (identical(type, "lognormal")) {
    mean <- as.numeric(spec$mean)
    sd <- as.numeric(spec$sd)
    sigma2 <- log1p((sd / mean)^2)
    return(stats::rlnorm(n, log(mean) - sigma2 / 2, sqrt(sigma2)))
  }
  stop("Unsupported population covariate distribution: ", type)
}

derive_population_covariate <- function(name, spec, sampled, n) {
  type <- spec$type
  if (identical(type, "cockcroft_creatinine")) {
    required <- c("WT", "AGE", "SEX")
    if (identical(spec$weightMode %||% "total_body_weight", "lean_body_weight")) required <- c(required, "HT")
    if (!all(required %in% names(sampled))) stop("Cockcroft creatinine requires WT, AGE, and SEX before CREAT.")
    weight <- sampled$WT
    if (identical(spec$weightMode %||% "total_body_weight", "lean_body_weight")) {
      weight <- ifelse(
        sampled$SEX < 0.5,
        0.407 * sampled$WT + 0.267 * sampled$HT - 19.2,
        0.252 * sampled$WT + 0.473 * sampled$HT - 48.3
      )
    }
    factor <- ifelse(
      sampled$SEX < 0.5,
      as.numeric(spec$maleFactor %||% 1.25),
      as.numeric(spec$femaleFactor %||% 1.04)
    )
    numerator <- factor * weight * (140 - sampled$AGE)
    crcl_min <- rep(as.numeric(spec$crclMin), n)
    crcl_max <- rep(as.numeric(spec$crclMax), n)
    if (!is.null(spec$valueMax)) crcl_min <- pmax(crcl_min, numerator / as.numeric(spec$valueMax))
    if (!is.null(spec$valueMin)) crcl_max <- pmin(crcl_max, numerator / as.numeric(spec$valueMin))
    if (any(crcl_min > crcl_max)) stop("Cockcroft creatinine constraints are incompatible with the sampled population.")
    target_crcl <- crcl_min + stats::runif(n) * (crcl_max - crcl_min)
    return(numerator / target_crcl)
  }
  if (identical(type, "schwartz_height")) {
    if (!"CREAT" %in% names(sampled)) stop("Schwartz height requires CREAT before HT.")
    target_egfr <- stats::runif(n, as.numeric(spec$egfrMin), as.numeric(spec$egfrMax))
    return(target_egfr * (sampled$CREAT / 88.4) / 0.413)
  }
  if (identical(type, "bmi_height")) {
    if (!"WT" %in% names(sampled)) stop("BMI-derived height requires WT before HT.")
    bmi <- sample_truncated_normal(n, as.numeric(spec$bmiMean), as.numeric(spec$bmiSd), as.numeric(spec$bmiMin %||% 10))
    return(100 * sqrt(sampled$WT / bmi))
  }
  stop("Unsupported derived population covariate ", name, ": ", type)
}

sample_covariates <- function(base_id, generator_id, n) {
  base <- model_covariate_defaults(base_id)
  generator <- model_covariate_defaults(generator_id)
  names_all <- union(names(base), names(generator))
  if (!length(names_all)) return(data.frame(row.names = seq_len(n)))
  generator_specs <- population_profile(generator_id)$covariates %||% list()
  base_specs <- population_profile(base_id)$covariates %||% list()
  values <- list()
  derived <- list()
  for (name in names_all) {
    from_generator <- name %in% names(generator)
    center <- if (from_generator) generator[[name]] else base[[name]]
    specs <- if (from_generator) generator_specs else base_specs
    spec <- specs[[name]]
    if (is.null(spec)) {
      values[[name]] <- rep(center, n)
    } else if ((spec$type %||% "") %in% c("cockcroft_creatinine", "schwartz_height", "bmi_height")) {
      derived[[name]] <- spec
    } else {
      values[[name]] <- sample_population_spec(spec, n)
    }
  }
  for (name in names(derived)) values[[name]] <- derive_population_covariate(name, derived[[name]], values, n)
  as.data.frame(values[names_all], check.names = FALSE)
}

sample_eta_matrix <- function(model, n) {
  omega <- as.matrix(mrgsolve::omat(model))
  if (!nrow(omega)) return(data.frame(row.names = seq_len(n)))
  decomposition <- eigen((omega + t(omega)) / 2, symmetric = TRUE)
  root <- decomposition$vectors %*% diag(sqrt(pmax(decomposition$values, 0)), nrow = nrow(omega))
  values <- matrix(stats::rnorm(n * nrow(omega)), nrow = n) %*% t(root)
  limits <- 2 * sqrt(pmax(diag(omega), 0))
  for (index in seq_len(ncol(values))) values[, index] <- pmax(-limits[[index]], pmin(limits[[index]], values[, index]))
  colnames(values) <- paste0("ETA", seq_len(ncol(values)))
  as.data.frame(values, check.names = FALSE)
}

sample_regimens <- function(scope, n) {
  config <- regimen_config(scope)
  draw <- function(values, size) values[sample.int(length(values), size, replace = TRUE)]
  amount <- draw(as.numeric(unlist(config$amounts)), n)
  interval <- draw(as.numeric(unlist(config$intervals)), n)
  mode <- scope$mode[[1]]
  infusion <- if (identical(mode, "ORAL")) {
    rep(0, n)
  } else if (identical(mode, "IV_CONTINUOUS")) {
    interval
  } else {
    candidates <- as.numeric(unlist(config$infusions))
    vapply(interval, function(ii) draw(candidates[candidates < ii], 1), numeric(1))
  }
  data.frame(amount = amount, interval = interval, infusion = infusion)
}

# Paired research benchmarks can override these hooks without changing the
# published model files or the historical artifact-training pipeline.
align_regimen_covariates <- function(base_scope, generator_scope, covariates, regimens) covariates
observation_model <- function(model) model

sample_times <- function(regimens, mode) {
  do.call(rbind, lapply(seq_len(nrow(regimens)), function(index) {
    interval <- regimens$interval[[index]]
    infusion_end <- if (mode %in% c("ORAL", "IV_CONTINUOUS") || regimens$infusion[[index]] <= 0) 0 else regimens$infusion[[index]]
    post_time <- infusion_end + 1
    if (post_time >= interval) stop(
      "The post-infusion sample falls outside the dosing interval: mode=", mode,
      ", interval=", interval, ", infusion=", regimens$infusion[[index]], "."
    )
    data.frame(
      ID = index,
      sample = c("C0", "POST"),
      time = c(if (mode == "IV_CONTINUOUS") 0 else interval - 1e-6, post_time),
      feature_time = c(if (mode == "IV_CONTINUOUS") 0 else interval, post_time),
      stringsAsFactors = FALSE
    )
  }))
}

trap_auc_vector <- function(time, concentration) {
  sum(diff(time) * (head(concentration, -1L) + tail(concentration, -1L)) / 2)
}

simulate_batch <- function(base_scope, generator_scope, n) {
  base_id <- base_scope$model_id[[1]]
  generator_id <- generator_scope$model_id[[1]]
  generator_record <- model_record(generator_id)
  model <- compiled_model(generator_id)
  covariates <- sample_covariates(base_id, generator_id, n)
  etas <- sample_eta_matrix(model, n)
  regimens <- sample_regimens(base_scope, n)
  covariates <- align_regimen_covariates(base_scope, generator_scope, covariates, regimens)
  samples <- sample_times(regimens, base_scope$mode[[1]])
  continuous <- identical(base_scope$mode[[1]], "IV_CONTINUOUS")

  parameter_names <- model_param_names(model)
  individual <- cbind(ID = seq_len(n), covariates, etas)
  individual <- individual[, c("ID", intersect(setdiff(names(individual), "ID"), parameter_names)), drop = FALSE]
  compartment <- match(model_administration_cmt(generator_record, base_scope$route[[1]]), model@cmtL)
  events <- data.frame(
    ID = seq_len(n), time = 0, evid = 1, amt = if (continuous) 0 else regimens$amount,
    ii = if (continuous) 0 else regimens$interval,
    addl = if (continuous) 0 else pmax(0, floor((24 - 1e-8) / regimens$interval)),
    ss = 1, cmt = compartment,
    rate = if (continuous) regimens$amount / regimens$interval else ifelse(regimens$infusion > 0, regimens$amount / regimens$infusion, 0)
  )
  ss_n <- if (continuous) 5000 else 500
  ss_tolerance <- if (continuous) 1e-6 else 1e-8

  output <- model |>
    mrgsolve::zero_re() |>
    mrgsolve::idata_set(individual) |>
    mrgsolve::data_set(events) |>
    mrgsolve::mrgsim(start = 0, end = 24, delta = 0.1, recsort = 3,
                     ss_n = ss_n, ss_rtol = ss_tolerance, ss_atol = ss_tolerance) |>
    as.data.frame()
  if (!"DV" %in% names(output)) stop("Model ", generator_id, " does not capture DV.")

  base_model <- compiled_model(base_id)
  base_record <- model_record(base_id)
  base_parameters <- model_param_names(base_model)
  base_individual <- cbind(ID = seq_len(n), covariates)
  base_individual <- base_individual[, c("ID", intersect(setdiff(names(base_individual), "ID"), base_parameters)), drop = FALSE]
  base_events <- events
  base_events$cmt <- match(model_administration_cmt(base_record, base_scope$route[[1]]), base_model@cmtL)
  population_horizon <- max(24, samples$time)
  population_output <- base_model |>
    mrgsolve::zero_re() |>
    mrgsolve::idata_set(base_individual) |>
    mrgsolve::data_set(base_events) |>
    mrgsolve::mrgsim(start = 0, end = population_horizon, delta = 0.1, recsort = 3,
                     ss_n = ss_n, ss_rtol = ss_tolerance, ss_atol = ss_tolerance) |>
    as.data.frame()
  if (!"DV" %in% names(population_output)) stop("Model ", base_id, " does not capture DV.")

  omega <- as.matrix(mrgsolve::omat(model))
  stochastic_model <- observation_model(model)
  observation_records <- data.frame(
    ID = samples$ID,
    time = samples$time,
    evid = 0,
    amt = 0,
    ii = 0,
    addl = 0,
    ss = 0,
    cmt = tagged_compartment(stochastic_model, "OBS"),
    rate = 0
  )
  stochastic_data <- rbind(events, observation_records)
  stochastic_data <- stochastic_data[order(stochastic_data$ID, stochastic_data$time, -stochastic_data$evid), , drop = FALSE]
  stochastic_output <- stochastic_model |>
    mrgsolve::zero_re(omega) |>
    mrgsolve::idata_set(individual) |>
    mrgsolve::data_set(stochastic_data) |>
    mrgsolve::mrgsim(obsonly = TRUE, recsort = 3,
                     ss_n = ss_n, ss_rtol = ss_tolerance, ss_atol = ss_tolerance) |>
    as.data.frame()
  if (!"DV" %in% names(stochastic_output)) stop("Model ", generator_id, " does not capture stochastic DV.")

  rows <- lapply(seq_len(n), function(index) {
    profile <- output[output$ID == index, c("time", "DV"), drop = FALSE]
    profile <- profile[!duplicated(profile$time, fromLast = TRUE), , drop = FALSE]
    concentration <- pmax(0, profile$DV)
    true_auc24 <- trap_auc_vector(profile$time, concentration)
    population_profile <- population_output[population_output$ID == index, c("time", "DV"), drop = FALSE]
    population_profile <- population_profile[!duplicated(population_profile$time, fromLast = TRUE), , drop = FALSE]
    population_auc_profile <- population_profile[population_profile$time <= 24 + 1e-8, , drop = FALSE]
    population_auc24 <- trap_auc_vector(population_auc_profile$time, pmax(0, population_auc_profile$DV))
    sample_rows <- samples[samples$ID == index, , drop = FALSE]
    times <- sample_rows$time
    feature_times <- sample_rows$feature_time
    sparse_rows <- stochastic_output[stochastic_output$ID == index, c("time", "DV"), drop = FALSE]
    sparse_rows <- sparse_rows[order(sparse_rows$time), , drop = FALSE]
    sparse <- pmax(0.001, vapply(times, function(sample_time) {
      sparse_rows$DV[[which.min(abs(sparse_rows$time - sample_time))]]
    }, numeric(1)))
    population_sparse <- stats::approx(
      population_profile$time,
      pmax(0, population_profile$DV),
      xout = times,
      rule = 2,
      ties = "ordered"
    )$y
    base_covariates <- covariates[index, intersect(names(model_covariate_defaults(base_id)), names(covariates)), drop = FALSE]
    data.frame(
      patient_id = index,
      TRUE_AUC24 = true_auc24,
      base_covariates,
      POP_AUC24 = population_auc24,
      DOSE = regimens$amount[[index]],
      INTERVAL = regimens$interval[[index]],
      INFUSION = regimens$infusion[[index]],
      PREV_CONC = sparse[[1]],
      PREV_TIME = feature_times[[1]],
      LAST_CONC = sparse[[2]],
      LAST_TIME = feature_times[[2]],
      PREV_POP_CONC = population_sparse[[1]],
      LAST_POP_CONC = population_sparse[[2]],
      PREV_CONC_RATIO = max(1e-4, min(1e4, sparse[[1]] / max(population_sparse[[1]], 1e-8))),
      LAST_CONC_RATIO = max(1e-4, min(1e4, sparse[[2]] / max(population_sparse[[2]], 1e-8))),
      CONC_DIFF = sparse[[2]] - sparse[[1]],
      TIME_DIFF = feature_times[[2]],
      N_OBS = 2,
      HAS_PREV = 1,
      check.names = FALSE
    )
  })
  data <- do.call(rbind, rows)
  numeric_columns <- vapply(data, is.numeric, logical(1))
  auc_ratio <- data$TRUE_AUC24 / data$POP_AUC24
  keep <- stats::complete.cases(data) & apply(data[, numeric_columns, drop = FALSE], 1, function(row) all(is.finite(row)) && all(abs(row) < 1e12)) &
    is.finite(auc_ratio) & auc_ratio >= MIN_AUC_RATIO & auc_ratio <= MAX_AUC_RATIO &
    data$TRUE_AUC24 > 0 & data$PREV_CONC > 0 & data$LAST_CONC > 0
  data[keep, , drop = FALSE]
}

safe_simulate_batch <- function(base_scope, generator_scope, n) {
  tryCatch(
    simulate_batch(base_scope, generator_scope, n),
    error = function(error) {
      if (n <= 1L) {
        message("  profile skipped: ", conditionMessage(error))
        return(NULL)
      }
      left <- safe_simulate_batch(base_scope, generator_scope, floor(n / 2))
      right <- safe_simulate_batch(base_scope, generator_scope, ceiling(n / 2))
      rbind(left, right)
    }
  )
}

make_cohort <- function(base_scope, generator_scope, n, label) {
  cat("  ", label, ": ", n, " profiles from ", generator_scope$model_id[[1]], "\n", sep = "")
  rows <- list()
  attempts <- 0L
  while (sum(vapply(rows, nrow, integer(1))) < n && attempts < ceiling(n / 100) * 3L) {
    attempts <- attempts + 1L
    current <- sum(vapply(rows, nrow, integer(1)))
    requested <- min(100L, ceiling((n - current) * 1.1))
    batch <- safe_simulate_batch(base_scope, generator_scope, requested)
    if (!is.null(batch) && nrow(batch)) rows[[length(rows) + 1L]] <- batch
  }
  cohort <- do.call(rbind, rows)
  if (is.null(cohort) || nrow(cohort) < n) stop("Too few valid profiles for ", label, ".")
  cohort <- cohort[seq_len(n), , drop = FALSE]
  cohort$patient_id <- seq_len(n)
  cohort
}

make_alternate_cohort <- function(base_scope, peers, n) {
  if (!nrow(peers)) return(NULL)
  counts <- rep(floor(n / nrow(peers)), nrow(peers))
  counts[seq_len(n %% nrow(peers))] <- counts[seq_len(n %% nrow(peers))] + 1L
  rows <- lapply(seq_len(nrow(peers)), function(index) {
    make_cohort(base_scope, peers[index, , drop = FALSE], counts[[index]], "alternate PopPK")
  })
  data <- do.call(rbind, rows)
  data$patient_id <- seq_len(nrow(data))
  data
}

NON_PREDICTORS <- c("patient_id", "TRUE_AUC24")
predictor_names <- function(data) {
  candidates <- setdiff(names(data), NON_PREDICTORS)
  candidates[vapply(data[candidates], function(values) stats::sd(values) > sqrt(.Machine$double.eps), logical(1))]
}

c0_only_dataset <- function(data) {
  data$LAST_CONC <- data$PREV_CONC
  data$LAST_TIME <- data$PREV_TIME
  data$LAST_POP_CONC <- data$PREV_POP_CONC
  data$LAST_CONC_RATIO <- data$PREV_CONC_RATIO
  data$PREV_CONC <- 0
  data$PREV_TIME <- 0
  data$PREV_POP_CONC <- 0
  data$PREV_CONC_RATIO <- 0
  data$CONC_DIFF <- 0
  data$TIME_DIFF <- 0
  data$N_OBS <- 1
  data$HAS_PREV <- 0
  data
}

fold_ids <- function(n, v, seed_value) {
  set.seed(seed_value)
  sample(rep(seq_len(min(v, n)), length.out = n))
}

fit_booster <- function(data, predictors) {
  matrix <- as.matrix(data[, predictors, drop = FALSE])
  dtrain <- xgboost::xgb.DMatrix(matrix, label = log(data$TRUE_AUC24 / data$POP_AUC24))
  xgboost::xgb.train(
    params = list(
      objective = "reg:squarederror",
      eval_metric = "rmse",
      base_score = mean(log(data$TRUE_AUC24 / data$POP_AUC24)),
      max_depth = XGB_PARAMETERS$max_depth,
      eta = XGB_PARAMETERS$eta,
      min_child_weight = XGB_PARAMETERS$min_child_weight,
      subsample = XGB_PARAMETERS$subsample,
      colsample_bytree = XGB_PARAMETERS$colsample_bytree,
      nthread = 1
    ),
    data = dtrain,
    nrounds = XGB_PARAMETERS$nrounds,
    verbose = 0
  )
}

auc_metrics <- function(truth, prediction) {
  keep <- is.finite(truth) & truth > 0 & is.finite(prediction)
  truth <- truth[keep]
  prediction <- prediction[keep]
  errors <- prediction - truth
  c(
    relative_rmse_pct = 100 * sqrt(mean(errors^2)) / mean(truth),
    relative_bias_pct = 100 * mean(errors) / mean(truth),
    within_20_pct = 100 * mean(abs(errors / truth) <= 0.2)
  )
}

repeated_cross_validate <- function(data, predictors, seed_value) {
  truth <- prediction <- numeric()
  for (repeat_index in seq_len(OUTER_REPEATS)) {
    ids <- fold_ids(nrow(data), OUTER_V, seed_value + repeat_index)
    for (fold in sort(unique(ids))) {
      booster <- fit_booster(data[ids != fold, , drop = FALSE], predictors)
      truth <- c(truth, data$TRUE_AUC24[ids == fold])
      prediction <- c(prediction, exp(stats::predict(booster, as.matrix(data[ids == fold, predictors, drop = FALSE]))) * data$POP_AUC24[ids == fold])
    }
  }
  auc_metrics(truth, prediction)
}

validation_record <- function(metrics) {
  list(
    passed = metrics[["relative_rmse_pct"]] <= MAX_RELATIVE_RMSE_PCT &&
      abs(metrics[["relative_bias_pct"]]) <= MAX_ABSOLUTE_BIAS_PCT &&
      metrics[["within_20_pct"]] >= MIN_WITHIN_20_PCT,
    relativeRmsePct = unname(metrics[["relative_rmse_pct"]]),
    relativeBiasPct = unname(metrics[["relative_bias_pct"]]),
    within20Pct = unname(metrics[["within_20_pct"]])
  )
}

rmse_gain_pct <- function(without_ml, with_ml) {
  baseline <- without_ml[["relative_rmse_pct"]]
  if (!is.finite(baseline) || baseline <= 0) return(NA_real_)
  100 * (baseline - with_ml[["relative_rmse_pct"]]) / baseline
}

feature_schema <- function(base_id, predictors) {
  covariates <- names(model_covariate_defaults(base_id))
  regimen <- c("DOSE", "INTERVAL", "INFUSION")
  observation <- c("PREV_CONC", "PREV_TIME", "LAST_CONC", "LAST_TIME", "CONC_DIFF", "TIME_DIFF", "N_OBS", "HAS_PREV")
  population_observation <- c("PREV_POP_CONC", "LAST_POP_CONC", "PREV_CONC_RATIO", "LAST_CONC_RATIO")
  lapply(predictors, function(name) {
    source <- if (name %in% covariates) "covariate" else if (name %in% regimen) "regimen" else if (name == "POP_AUC24") "population_auc24" else if (name %in% population_observation) "population_observation" else if (name %in% observation) "observation" else stop("Unknown feature ", name)
    list(name = name, source = source, key = name)
  })
}

evaluate_scope <- function(scope, index) {
  base_id <- scope$model_id[[1]]
  population <- population_sampling_record(base_id)
  mode <- scope$mode[[1]]
  peers <- catalog_scopes()
  peers <- peers[peers$drug_key == scope$drug_key[[1]] & peers$mode == mode & peers$model_id != base_id, , drop = FALSE]
  cat("\n=== ", base_id, " / ", mode, " ===\n", sep = "")
  development <- make_cohort(scope, scope, n_patients, "development")
  ratio_quantiles <- stats::quantile(development$TRUE_AUC24 / development$POP_AUC24, c(0.01, 0.5, 0.99), na.rm = TRUE)
  cat(sprintf("  AUC24 ratio p01/median/p99: %.3f / %.3f / %.3f\n", ratio_quantiles[[1]], ratio_quantiles[[2]], ratio_quantiles[[3]]))

  set.seed(seed + index)
  holdout_indices <- sample(seq_len(nrow(development)), max(6L, floor(0.25 * nrow(development))))
  holdout <- development[holdout_indices, , drop = FALSE]
  training <- development[-holdout_indices, , drop = FALSE]
  if (anyDuplicated(development$patient_id) ||
      length(intersect(training$patient_id, holdout$patient_id))) {
    stop("Patient leakage detected between ML training and untouched holdout partitions.")
  }
  predictors <- predictor_names(training)
  repeated <- repeated_cross_validate(training, predictors, seed + index * 1000L)
  repeated_without_ml <- auc_metrics(training$TRUE_AUC24, training$POP_AUC24)
  booster <- fit_booster(training, predictors)
  holdout_prediction <- exp(stats::predict(booster, as.matrix(holdout[, predictors, drop = FALSE]))) * holdout$POP_AUC24
  holdout_metrics <- auc_metrics(holdout$TRUE_AUC24, holdout_prediction)
  holdout_without_ml <- auc_metrics(holdout$TRUE_AUC24, holdout$POP_AUC24)
  alternate <- make_alternate_cohort(scope, peers, if (smoke) n_patients else min(n_patients, 500L))
  alternate_metrics <- if (is.null(alternate)) NULL else auc_metrics(
    alternate$TRUE_AUC24,
    exp(stats::predict(booster, as.matrix(alternate[, predictors, drop = FALSE]))) * alternate$POP_AUC24
  )
  alternate_without_ml <- if (is.null(alternate)) NULL else auc_metrics(alternate$TRUE_AUC24, alternate$POP_AUC24)

  c0_development <- c0_only_dataset(development)
  c0_holdout <- c0_development[holdout_indices, , drop = FALSE]
  c0_training <- c0_development[-holdout_indices, , drop = FALSE]
  c0_predictors <- predictor_names(c0_training)
  c0_repeated <- repeated_cross_validate(c0_training, c0_predictors, seed + index * 1000L)
  c0_booster <- fit_booster(c0_training, c0_predictors)
  c0_holdout_prediction <- exp(stats::predict(c0_booster, as.matrix(c0_holdout[, c0_predictors, drop = FALSE]))) * c0_holdout$POP_AUC24
  c0_holdout_metrics <- auc_metrics(c0_holdout$TRUE_AUC24, c0_holdout_prediction)
  c0_alternate <- if (is.null(alternate)) NULL else c0_only_dataset(alternate)
  c0_alternate_metrics <- if (is.null(c0_alternate)) NULL else auc_metrics(
    c0_alternate$TRUE_AUC24,
    exp(stats::predict(c0_booster, as.matrix(c0_alternate[, c0_predictors, drop = FALSE]))) * c0_alternate$POP_AUC24
  )
  improves_holdout <- holdout_metrics[["relative_rmse_pct"]] < holdout_without_ml[["relative_rmse_pct"]]
  gates <- isTRUE(validation_record(repeated)$passed) && isTRUE(validation_record(holdout_metrics)$passed) && improves_holdout
  print_metrics <- function(label, metrics) cat(sprintf(
    "  %-18s rRMSE %6.2f%% | bias %+6.2f%% | within 20%% %5.1f%%\n",
    label, metrics[["relative_rmse_pct"]], metrics[["relative_bias_pct"]], metrics[["within_20_pct"]]
  ))
  print_metrics("CV without ML", repeated_without_ml)
  print_metrics("CV with ML", repeated)
  print_metrics("test without ML", holdout_without_ml)
  print_metrics("test with ML", holdout_metrics)
  print_metrics("test C0 only", c0_holdout_metrics)
  cat(sprintf("  holdout rRMSE gain: %+6.2f%%\n", rmse_gain_pct(holdout_without_ml, holdout_metrics)))
  if (!is.null(alternate_metrics)) {
    print_metrics("other PopPK no ML", alternate_without_ml)
    print_metrics("other PopPK + ML", alternate_metrics)
  }
  cat("  research gates: ", if (gates) "PASS" else "FAIL", "\n", sep = "")

  set.seed(seed + index * 10000L)
  background_rows <- sample(seq_len(nrow(training)), min(DALEX_BACKGROUND_SIZE, nrow(training)))
  list(
    row = data.frame(
      base_model = base_id,
      drug = scope$drug[[1]],
      route = scope$route[[1]],
      administration_mode = mode,
      population_source = population$source,
      varied_covariates = paste(population$varied, collapse = "+"),
      reference_covariates = paste(population$fixed, collapse = "+"),
      peers = paste(peers$model_id, collapse = "+"),
      n_training = nrow(training),
      n_holdout = nrow(holdout),
      n_alternate = if (is.null(alternate)) 0L else nrow(alternate),
      cv_without_ml_relative_rmse_pct = repeated_without_ml[["relative_rmse_pct"]],
      cv_without_ml_relative_bias_pct = repeated_without_ml[["relative_bias_pct"]],
      cv_with_ml_relative_rmse_pct = repeated[["relative_rmse_pct"]],
      cv_with_ml_relative_bias_pct = repeated[["relative_bias_pct"]],
      cv_relative_rmse_pct = repeated[["relative_rmse_pct"]],
      cv_relative_bias_pct = repeated[["relative_bias_pct"]],
      cv_within_20_pct = repeated[["within_20_pct"]],
      holdout_without_ml_relative_rmse_pct = holdout_without_ml[["relative_rmse_pct"]],
      holdout_without_ml_relative_bias_pct = holdout_without_ml[["relative_bias_pct"]],
      holdout_with_ml_relative_rmse_pct = holdout_metrics[["relative_rmse_pct"]],
      holdout_with_ml_relative_bias_pct = holdout_metrics[["relative_bias_pct"]],
      holdout_relative_rmse_gain_pct = rmse_gain_pct(holdout_without_ml, holdout_metrics),
      holdout_relative_rmse_pct = holdout_metrics[["relative_rmse_pct"]],
      holdout_relative_bias_pct = holdout_metrics[["relative_bias_pct"]],
      holdout_within_20_pct = holdout_metrics[["within_20_pct"]],
      c0_cv_relative_rmse_pct = c0_repeated[["relative_rmse_pct"]],
      c0_cv_relative_bias_pct = c0_repeated[["relative_bias_pct"]],
      c0_holdout_relative_rmse_pct = c0_holdout_metrics[["relative_rmse_pct"]],
      c0_holdout_relative_bias_pct = c0_holdout_metrics[["relative_bias_pct"]],
      c0_holdout_within_20_pct = c0_holdout_metrics[["within_20_pct"]],
      c0_alternate_relative_rmse_pct = if (is.null(c0_alternate_metrics)) NA_real_ else c0_alternate_metrics[["relative_rmse_pct"]],
      c0_alternate_relative_bias_pct = if (is.null(c0_alternate_metrics)) NA_real_ else c0_alternate_metrics[["relative_bias_pct"]],
      alternate_without_ml_relative_rmse_pct = if (is.null(alternate_without_ml)) NA_real_ else alternate_without_ml[["relative_rmse_pct"]],
      alternate_without_ml_relative_bias_pct = if (is.null(alternate_without_ml)) NA_real_ else alternate_without_ml[["relative_bias_pct"]],
      alternate_with_ml_relative_rmse_pct = if (is.null(alternate_metrics)) NA_real_ else alternate_metrics[["relative_rmse_pct"]],
      alternate_with_ml_relative_bias_pct = if (is.null(alternate_metrics)) NA_real_ else alternate_metrics[["relative_bias_pct"]],
      alternate_relative_rmse_pct = if (is.null(alternate_metrics)) NA_real_ else alternate_metrics[["relative_rmse_pct"]],
      alternate_relative_bias_pct = if (is.null(alternate_metrics)) NA_real_ else alternate_metrics[["relative_bias_pct"]],
      alternate_within_20_pct = if (is.null(alternate_metrics)) NA_real_ else alternate_metrics[["within_20_pct"]],
      research_gates_pass = gates,
      artifact_saved = FALSE,
      stringsAsFactors = FALSE
    ),
    booster = booster,
    background = training[background_rows, predictors, drop = FALSE],
    predictors = predictors,
    schema = feature_schema(base_id, predictors),
    feature_bounds = stats::setNames(lapply(predictors, function(name) list(min = min(training[[name]]), max = max(training[[name]]))), predictors),
    auc_bounds = list(min = min(training$TRUE_AUC24), max = max(training$TRUE_AUC24)),
    concentration_max = max(training$PREV_CONC, training$LAST_CONC),
    repeated = repeated,
    repeated_without_ml = repeated_without_ml,
    holdout = holdout_metrics,
    holdout_without_ml = holdout_without_ml,
    alternate = alternate_metrics,
    alternate_without_ml = alternate_without_ml,
    c0 = list(
      protocol = C0_SAMPLING_PROTOCOL,
      repeated = c0_repeated,
      holdout = c0_holdout_metrics,
      alternate = c0_alternate_metrics
    ),
    scope = scope,
    peers = peers,
    population = population
  )
}

publish_candidate <- function(evaluation) {
  scope <- evaluation$scope
  base_id <- scope$model_id[[1]]
  mode_id <- unname(MODE_ID[[scope$mode[[1]]]])
  artifact_id <- paste0(base_id, "-", mode_id, "-c0-post-auc24-xgb-v4")
  artifact_directory <- file.path(APP_ROOT, "ml", "artifacts")
  dir.create(artifact_directory, recursive = TRUE, showWarnings = FALSE)
  artifact_file <- file.path(artifact_directory, paste0(artifact_id, ".rds"))
  background_file <- file.path(artifact_directory, paste0(artifact_id, "-dalex-background.rds"))
  saveRDS(evaluation$booster, artifact_file)
  saveRDS(evaluation$background, background_file)

  config <- regimen_config(scope)
  manifest <- read_ml_manifest()
  alternate_validation <- if (is.null(evaluation$alternate)) {
    list(status = "not_applicable", passed = NULL)
  } else {
    validation_record(evaluation$alternate)
  }
  artifact <- list(
    id = artifact_id,
    releaseLevel = if (isTRUE(evaluation$row$research_gates_pass[[1]])) "research" else "experimental",
    drug = scope$drug[[1]],
    route = scope$route[[1]],
    administrationMode = mode_id,
    samplingProtocol = SAMPLING_PROTOCOL,
    baseModelId = base_id,
    baseModelSha256 = model_sha256(base_id),
    artifactPath = paste0("artifacts/", basename(artifact_file)),
    artifactSha256 = digest::digest(artifact_file, algo = "sha256", file = TRUE, serialize = FALSE),
    explanation = list(
      type = "dalex_break_down",
      backgroundPath = paste0("artifacts/", basename(background_file)),
      backgroundSha256 = digest::digest(background_file, algo = "sha256", file = TRUE, serialize = FALSE),
      sampleSize = nrow(evaluation$background),
      synthetic = TRUE
    ),
    featureSchema = evaluation$schema,
    prediction = list(type = "auc24_direct", metric = "AUC24", unit = config$aucUnit, horizonHours = 24, transform = "exp_times_feature", scaleFeature = "POP_AUC24"),
    methodology = list(
      label = "Simulation-trained hybrid AUC prediction from C0 plus one-hour post-infusion sampling",
      doi = "10.1016/j.phrs.2021.105578",
      developmentGenerators = base_id,
      developmentModelSha256 = stats::setNames(list(model_sha256(base_id)), base_id),
      unseenValidationGenerators = evaluation$peers$model_id,
      regimenConfiguration = "ml/training-regimens.json",
      populationConfiguration = "ml/training-populations.json",
      sourcePopulation = evaluation$population$source,
      variedCovariates = evaluation$population$varied,
      referenceCovariates = evaluation$population$fixed,
      target = "log(individual AUC24 / population AUC24)",
      acceptedAucRatio = list(min = MIN_AUC_RATIO, max = MAX_AUC_RATIO),
      realPatientValidation = "Pending"
    ),
    training = list(
      seed = seed,
      nVirtualDevelopment = evaluation$row$n_training[[1]] + evaluation$row$n_holdout[[1]],
      nTraining = evaluation$row$n_training[[1]],
      nUntouchedHoldout = evaluation$row$n_holdout[[1]],
      nUnseenPopPk = evaluation$row$n_alternate[[1]],
      split = list(
        unit = "simulated_patient",
        method = "seeded_random_75_25_before_model_fitting",
        holdoutUntouched = TRUE
      ),
      hyperparameters = XGB_PARAMETERS,
      software = list(
        R = R.version.string,
        mrgsolve = as.character(utils::packageVersion("mrgsolve")),
        xgboost = as.character(utils::packageVersion("xgboost"))
      )
    ),
    trainingDomain = list(
      auc24 = evaluation$auc_bounds,
      concentrationMax = evaluation$concentration_max,
      features = evaluation$feature_bounds
    ),
    validation = list(
      thresholds = list(relativeRmsePct = MAX_RELATIVE_RMSE_PCT, absoluteBiasPct = MAX_ABSOLUTE_BIAS_PCT, within20Pct = MIN_WITHIN_20_PCT, requiresImprovementOverPopulationBaseline = TRUE),
      repeatedCv = validation_record(evaluation$repeated),
      untouchedHoldout = validation_record(evaluation$holdout),
      alternatePopPk = alternate_validation,
      populationBaseline = list(
        repeatedCv = validation_record(evaluation$repeated_without_ml),
        untouchedHoldout = validation_record(evaluation$holdout_without_ml),
        alternatePopPk = if (is.null(evaluation$alternate_without_ml)) list(status = "not_applicable", passed = NULL) else validation_record(evaluation$alternate_without_ml)
      ),
      comparison = list(
        repeatedCvRelativeRmseGainPct = rmse_gain_pct(evaluation$repeated_without_ml, evaluation$repeated),
        untouchedHoldoutRelativeRmseGainPct = rmse_gain_pct(evaluation$holdout_without_ml, evaluation$holdout),
        alternatePopPkRelativeRmseGainPct = if (is.null(evaluation$alternate)) NULL else rmse_gain_pct(evaluation$alternate_without_ml, evaluation$alternate)
      ),
      samplingStrategies = list(
        c0Only = list(
          protocol = C0_SAMPLING_PROTOCOL,
          repeatedCv = validation_record(evaluation$c0$repeated),
          untouchedHoldout = validation_record(evaluation$c0$holdout),
          alternatePopPk = if (is.null(evaluation$c0$alternate)) list(status = "not_applicable", passed = NULL) else validation_record(evaluation$c0$alternate)
        ),
        c0PlusOneHourPostInfusion = list(
          protocol = SAMPLING_PROTOCOL,
          repeatedCv = validation_record(evaluation$repeated),
          untouchedHoldout = validation_record(evaluation$holdout),
          alternatePopPk = alternate_validation
        )
      ),
      realPatient = list(status = "pending", gainPct = NULL)
    )
  )
  same_scope <- vapply(manifest$artifacts, function(item) {
    identical(item$baseModelId %||% "", base_id) &&
      identical(item$route %||% "", scope$route[[1]]) &&
      identical(ml_administration_mode(item$administrationMode), scope$mode[[1]]) &&
      identical(ml_prediction_type(item), "auc24_direct")
  }, logical(1))
  old_paths <- unlist(lapply(manifest$artifacts[same_scope], function(item) c(
    item$artifactPath %||% "",
    (item$explanation %||% list())$backgroundPath %||% ""
  )), use.names = FALSE)
  manifest$artifacts <- c(manifest$artifacts[!same_scope], list(artifact))
  manifest$version <- max(3L, as.integer(manifest$version %||% 0L))
  manifest$benchmarkDate <- format(Sys.Date(), "%Y-%m-%d")
  manifest$samplingDesign <- list(
    primary = SAMPLING_PROTOCOL,
    comparator = C0_SAMPLING_PROTOCOL,
    notApplicable = "Continuous infusion: no infusion end or pre-dose C0 within an ongoing infusion interval."
  )
  jsonlite::write_json(manifest, ML_MANIFEST_PATH, auto_unbox = TRUE, pretty = TRUE, null = "null")
  for (relative_path in old_paths[nzchar(old_paths)]) {
    old_path <- file.path(APP_ROOT, "ml", relative_path)
    if (file.exists(old_path) && normalizePath(old_path, winslash = "/") != normalizePath(artifact_file, winslash = "/") && normalizePath(old_path, winslash = "/") != normalizePath(background_file, winslash = "/")) unlink(old_path)
  }
  TRUE
}

prune_non_applicable_continuous_artifacts <- function() {
  manifest <- read_ml_manifest()
  remove <- vapply(manifest$artifacts, function(item) {
    identical(ml_administration_mode(item$administrationMode), "IV_CONTINUOUS") &&
      identical(ml_prediction_type(item), "auc24_direct")
  }, logical(1))
  if (!any(remove)) return(invisible(0L))
  old_paths <- unlist(lapply(manifest$artifacts[remove], function(item) c(
    item$artifactPath %||% "",
    (item$explanation %||% list())$backgroundPath %||% ""
  )), use.names = FALSE)
  manifest$artifacts <- manifest$artifacts[!remove]
  jsonlite::write_json(manifest, ML_MANIFEST_PATH, auto_unbox = TRUE, pretty = TRUE, null = "null")
  for (relative_path in old_paths[nzchar(old_paths)]) {
    old_path <- file.path(APP_ROOT, "ml", relative_path)
    if (file.exists(old_path)) unlink(old_path)
  }
  invisible(sum(remove))
}

if (!identical(Sys.getenv("PMX_TRAIN_FUNCTIONS_ONLY"), "1")) {
if ("--paired-benchmark" %in% args) {
  source(file.path(APP_ROOT, "ml", "paired_auc_benchmark.R"), local = TRUE)
  run_paired_benchmark()
} else {
cat("Generic simulation-trained direct AUC24 evaluation\n")
cat("Scopes: ", nrow(scopes), " | patients per scope: ", n_patients, "\n", sep = "")
cat("No patient data are read or written. Only synthetic profiles are used.\n")

evaluations <- vector("list", nrow(scopes))
for (index in seq_len(nrow(scopes))) {
  evaluations[[index]] <- evaluate_scope(scopes[index, , drop = FALSE], index)
  if (publish_artifacts) {
    saved <- publish_candidate(evaluations[[index]])
    evaluations[[index]]$row$artifact_saved <- saved
    cat("  artifact: ", if (saved) "published" else "not published", "\n", sep = "")
  }
}
if (full_library_publish) {
  removed <- prune_non_applicable_continuous_artifacts()
  if (removed) cat("Removed ", removed, " non-applicable continuous-infusion artifact(s).\n", sep = "")
}

results <- do.call(rbind, lapply(evaluations, `[[`, "row"))
if (nzchar(report_path)) {
  dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(results, report_path, row.names = FALSE)
  cat("Validation report written to ", normalizePath(report_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
}
cat(
  "Completed: ", nrow(results), " scope(s) | research: ", sum(results$research_gates_pass),
  " | experimental: ", sum(!results$research_gates_pass), "\n",
  sep = ""
)
if (!publish_artifacts) cat("\nNo artifact saved. Add --publish after a full evaluation run.\n")
}
}
