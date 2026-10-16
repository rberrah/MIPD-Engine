ML_MANIFEST_VERSION <- 3L
if (!exists("APP_ROOT", inherits = TRUE)) {
  APP_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
ML_ROOT <- file.path(APP_ROOT, "ml")
ML_MANIFEST_PATH <- file.path(ML_ROOT, "registry.json")
.ml_explainer_cache <- new.env(parent = emptyenv())

model_sha256 <- function(model_id) {
  if (!requireNamespace("digest", quietly = TRUE)) return(NA_character_)
  record <- model_record(model_id)
  path <- file.path(MODEL_ROOT, record$file[[1]])
  size <- file.info(path)$size
  content <- readBin(path, what = "raw", n = size)
  canonical <- content[content != as.raw(0x0d)]
  digest::digest(canonical, algo = "sha256", serialize = FALSE)
}

ml_contract_value <- function(object, path) {
  value <- object
  for (key in strsplit(path, ".", fixed = TRUE)[[1]]) {
    if (!is.list(value) || is.null(value[[key]])) return(NULL)
    value <- value[[key]]
  }
  value
}

ml_contract_present <- function(value) {
  if (is.null(value) || !length(value)) return(FALSE)
  if (is.character(value) && any(!nzchar(value))) return(FALSE)
  TRUE
}

validate_ml_manifest <- function(manifest) {
  if (!identical(as.integer(manifest$version %||% 0), ML_MANIFEST_VERSION)) {
    stop("Unsupported ML manifest version.")
  }
  artifacts <- manifest$artifacts %||% list()
  required <- c(
    "id", "releaseLevel", "baseModelId", "baseModelSha256", "artifactPath", "artifactSha256",
    "route", "administrationMode", "samplingProtocol.strategy", "samplingProtocol.minimumObservations",
    "samplingProtocol.target", "prediction.type", "prediction.metric", "prediction.unit",
    "prediction.horizonHours", "featureSchema", "methodology.sourcePopulation",
    "methodology.regimenConfiguration", "methodology.populationConfiguration", "training.seed",
    "training.nVirtualDevelopment", "training.nTraining", "training.nUntouchedHoldout",
    "training.nUnseenPopPk", "training.split.unit", "training.split.method",
    "training.split.holdoutUntouched", "training.software", "trainingDomain.features",
    "validation.untouchedHoldout.relativeRmsePct", "validation.untouchedHoldout.relativeBiasPct",
    "validation.untouchedHoldout.within20Pct"
  )
  ids <- vapply(artifacts, function(artifact) as.character(artifact$id %||% ""), character(1))
  if (any(!nzchar(ids)) || anyDuplicated(ids)) stop("ML artifact identifiers must be present and unique.")
  for (artifact in artifacts) {
    label <- as.character(artifact$id %||% "<unknown>")
    missing <- required[!vapply(required, function(path) ml_contract_present(ml_contract_value(artifact, path)), logical(1))]
    if (length(missing)) stop("Incomplete ML artifact ", label, ": ", paste(missing, collapse = ", "))
    if (!artifact$route %in% c("IV", "Oral")) stop("Unsupported ML route in ", label, ".")
    if (!artifact$administrationMode %in% c("intermittent", "oral", "continuous")) {
      stop("Unsupported ML administration mode in ", label, ".")
    }
    if (!identical(artifact$samplingProtocol$target, artifact$prediction$metric)) {
      stop("Sampling and prediction targets differ in ", label, ".")
    }
    horizon <- suppressWarnings(as.numeric(artifact$prediction$horizonHours))
    if (length(horizon) != 1L || !is.finite(horizon) || horizon <= 0) stop("Invalid ML horizon in ", label, ".")
    counts <- suppressWarnings(as.integer(c(
      artifact$training$nVirtualDevelopment,
      artifact$training$nTraining,
      artifact$training$nUntouchedHoldout
    )))
    if (anyNA(counts) || counts[[2]] + counts[[3]] != counts[[1]]) stop("Inconsistent ML split in ", label, ".")
    if (!identical(artifact$training$split$unit, "simulated_patient") ||
        !isTRUE(artifact$training$split$holdoutUntouched)) {
      stop("Invalid ML split contract in ", label, ".")
    }
    features <- vapply(artifact$featureSchema, function(feature) as.character(feature$name %||% ""), character(1))
    if (any(!nzchar(features)) || anyDuplicated(features)) stop("Invalid ML feature schema in ", label, ".")
    domains <- names(artifact$trainingDomain$features)
    if (!setequal(features, domains)) stop("ML features and training domains differ in ", label, ".")
    for (name in features) {
      bounds <- artifact$trainingDomain$features[[name]]
      range <- suppressWarnings(as.numeric(c(bounds$min, bounds$max)))
      if (length(range) != 2L || any(!is.finite(range)) || range[[1]] > range[[2]]) {
        stop("Invalid ML training domain for ", name, " in ", label, ".")
      }
    }
    metrics <- artifact$validation$untouchedHoldout
    values <- suppressWarnings(as.numeric(c(metrics$relativeRmsePct, metrics$relativeBiasPct, metrics$within20Pct)))
    if (length(values) != 3L || any(!is.finite(values)) || values[[3]] < 0 || values[[3]] > 100) {
      stop("Invalid ML holdout metrics in ", label, ".")
    }
  }
  manifest$artifacts <- artifacts
  manifest
}

read_ml_manifest <- function() {
  if (!file.exists(ML_MANIFEST_PATH)) return(list(version = ML_MANIFEST_VERSION, artifacts = list()))
  manifest <- jsonlite::fromJSON(ML_MANIFEST_PATH, simplifyVector = FALSE)
  validate_ml_manifest(manifest)
}

ml_prediction_type <- function(artifact) {
  as.character((artifact$prediction %||% list())$type %||% (artifact$correction %||% list())$type %||% "")
}

ml_prediction_value <- function(value, artifact, feature_data = NULL) {
  transform <- as.character((artifact$prediction %||% list())$transform %||% "identity")
  if (identical(transform, "identity")) return(as.numeric(value))
  if (identical(transform, "exp")) return(exp(as.numeric(value)))
  if (identical(transform, "exp_times_feature")) {
    feature <- as.character((artifact$prediction %||% list())$scaleFeature %||% "")
    if (is.null(feature_data) || !nzchar(feature) || !feature %in% names(feature_data)) {
      stop("ML prediction scale feature is unavailable.")
    }
    return(exp(as.numeric(value)) * as.numeric(feature_data[[feature]]))
  }
  stop("Unsupported ML prediction transform: ", transform)
}

ml_administration_mode <- function(value) {
  value <- toupper(as.character(value %||% ""))
  if (value %in% c("INTERMITTENT", "IV_INTERMITTENT")) return("IV_INTERMITTENT")
  if (value %in% c("CONTINUOUS", "IV_CONTINUOUS")) return("IV_CONTINUOUS")
  if (value == "ORAL") return("ORAL")
  ""
}

ml_artifact_eligibility <- function(artifact, model_id, drug, route, administration_mode = "") {
  validation <- artifact$validation %||% list()
  expected_mode <- ml_administration_mode(administration_mode)
  mode_ok <- !nzchar(expected_mode) || identical(ml_administration_mode(artifact$administrationMode), expected_mode)
  identity_ok <- identical(artifact$baseModelId %||% "", model_id) &&
    identical(artifact$drug %||% "", drug) &&
    identical(artifact$route %||% "", route) && mode_ok
  expected_hash <- model_sha256(model_id)
  hash_ok <- is.character(expected_hash) && length(expected_hash) == 1L && !is.na(expected_hash) &&
    identical(tolower(artifact$baseModelSha256 %||% ""), tolower(expected_hash))
  supported_type <- ml_prediction_type(artifact) %in% c("auc24_direct", "eta_additive")
  validation_passed <- function(name, legacy_name) {
    block <- validation[[name]] %||% list()
    if (is.list(block) && length(block)) return(isTRUE(block$passed))
    legacy_value <- suppressWarnings(as.numeric(validation[[legacy_name]] %||% NA_real_))
    is.finite(legacy_value) && legacy_value > 0
  }
  validation_evaluated <- function(name) {
    block <- validation[[name]] %||% list()
    is.list(block) && is.finite(suppressWarnings(as.numeric(block$relativeRmsePct %||% NA_real_))) &&
      is.finite(suppressWarnings(as.numeric(block$relativeBiasPct %||% NA_real_))) &&
      is.finite(suppressWarnings(as.numeric(block$within20Pct %||% NA_real_)))
  }
  nested_ok <- validation_passed("repeatedCv", "repeatedCvGainPct") ||
    validation_passed("repeatedNestedCv", "repeatedNestedCvGainPct")
  holdout_ok <- validation_passed("untouchedHoldout", "untouchedHoldoutGainPct")
  alternate_ok <- validation_passed("alternatePopPk", "alternatePopPkGainPct")
  repeated_evaluated <- validation_evaluated("repeatedCv") || validation_evaluated("repeatedNestedCv")
  experimental_ok <- identity_ok && hash_ok && supported_type && repeated_evaluated && validation_evaluated("untouchedHoldout")
  improvement_required <- isTRUE((validation$thresholds %||% list())$requiresImprovementOverPopulationBaseline)
  holdout_gain <- suppressWarnings(as.numeric((validation$comparison %||% list())$untouchedHoldoutRelativeRmseGainPct %||% NA_real_))
  improvement_ok <- !improvement_required || (is.finite(holdout_gain) && holdout_gain > 0)
  research_ok <- experimental_ok && nested_ok && holdout_ok && improvement_ok
  real_patient <- validation$realPatient %||% list()
  clinical_ok <- research_ok && identical(real_patient$status %||% "", "validated") &&
    is.finite(suppressWarnings(as.numeric(real_patient$gainPct %||% NA_real_))) &&
    as.numeric(real_patient$gainPct) > 0
  reasons <- c(
    if (!improvement_ok) "absence d'amelioration de la RMSE par rapport a l'AUC populationnelle sans ML",
    if (!identity_ok) "molécule, voie, mode d'administration ou modèle de base incompatible",
    if (!hash_ok) "empreinte SHA-256 du modèle incompatible ou indisponible",
    if (!supported_type) "type de prédiction ML non pris en charge",
    if (!nested_ok) "validation croisée répétée/imbriquée non favorable",
    if (!holdout_ok) "jeu de test interne non touché non favorable",
    if (!clinical_ok) "validation sur patients réels absente ou non favorable"
  )
  list(experimental = experimental_ok, research = research_ok, transportability = alternate_ok, clinical = clinical_ok, reasons = unique(reasons))
}

compatible_ml_artifacts <- function(model_id, drug, route, administration_mode = "") {
  if (model_id %in% MODEL_CATALOG$id && !model_analysis_eligible(model_record(model_id))) return(list())
  artifacts <- read_ml_manifest()$artifacts
  if (!length(artifacts)) return(list())
  Filter(function(artifact) {
    eligibility <- ml_artifact_eligibility(artifact, model_id, drug, route, administration_mode)
    isTRUE(eligibility$experimental)
  }, artifacts)
}

ml_latest_regimen <- function(fit) {
  doses <- fit$source_doses %||% data.frame()
  if (!nrow(doses)) stop("ML AUC24 requires an administered dose.")
  dose_ss <- if ("ss" %in% names(doses)) as.integer(doses$ss) else rep(0L, nrow(doses))
  interval <- suppressWarnings(as.numeric(doses$interval))
  count <- suppressWarnings(as.numeric(doses$count %||% rep(1, nrow(doses))))
  last_time <- doses$time + ifelse(dose_ss == 1L, 0, pmax(0, interval) * pmax(0, count - 1))
  regimen <- doses[which.max(last_time), , drop = FALSE]
  values <- c(
    DOSE = suppressWarnings(as.numeric(regimen$amount[[1]])),
    INTERVAL = suppressWarnings(as.numeric(regimen$interval[[1]])),
    INFUSION = suppressWarnings(as.numeric(regimen$infusion[[1]]))
  )
  if (any(!is.finite(values)) || values[["DOSE"]] <= 0 || values[["INTERVAL"]] <= 0 || values[["INFUSION"]] < 0) {
    stop("The current regimen is incompatible with ML AUC24.")
  }
  values
}

ml_latest_regimen_is_steady_state <- function(fit) {
  doses <- fit$source_doses %||% data.frame()
  if (!nrow(doses)) return(FALSE)
  dose_ss <- if ("ss" %in% names(doses)) as.integer(doses$ss) else rep(0L, nrow(doses))
  interval <- suppressWarnings(as.numeric(doses$interval))
  count <- suppressWarnings(as.numeric(doses$count %||% rep(1, nrow(doses))))
  last_time <- doses$time + ifelse(dose_ss == 1L, 0, pmax(0, interval) * pmax(0, count - 1))
  dose_ss[[which.max(last_time)]] == 1L
}

ml_dose_occurrences <- function(doses, end_time) {
  occurrences <- numeric()
  for (index in seq_len(nrow(doses))) {
    start <- suppressWarnings(as.numeric(doses$time[[index]]))
    interval <- suppressWarnings(as.numeric(doses$interval[[index]]))
    count <- suppressWarnings(as.integer((doses$count %||% rep(1L, nrow(doses)))[[index]]))
    steady_state <- suppressWarnings(as.integer((doses$ss %||% rep(0L, nrow(doses)))[[index]])) == 1L
    if (!is.finite(start)) next
    if (is.finite(interval) && interval > 0) {
      repetitions <- if (steady_state) {
        max(0L, floor((end_time - start) / interval))
      } else {
        max(0L, count - 1L)
      }
      occurrences <- c(occurrences, start + seq.int(0L, repetitions) * interval)
    } else {
      occurrences <- c(occurrences, start)
    }
  }
  sort(unique(occurrences[is.finite(occurrences) & occurrences <= end_time + 1e-8]))
}

ml_c0_post_observation_values <- function(fit, observations, strategy) {
  regimen <- ml_latest_regimen(fit)
  interval <- regimen[["INTERVAL"]]
  infusion_end <- if (regimen[["INFUSION"]] > 0) regimen[["INFUSION"]] else 0
  occurrences <- ml_dose_occurrences(fit$source_doses, max(observations$time))

  if (identical(strategy, "c0_only")) {
    last <- observations[nrow(observations), , drop = FALSE]
    exact_dose <- any(abs(occurrences - last$time[[1]]) <= 1e-8)
    previous_doses <- occurrences[occurrences < last$time[[1]] - 1e-8]
    relative_time <- if (exact_dose) interval else if (length(previous_doses)) last$time[[1]] - max(previous_doses) else last$time[[1]] - min(occurrences)
    return(c(
      PREV_CONC = 0, PREV_TIME = 0,
      LAST_CONC = last$concentration[[1]], LAST_TIME = relative_time,
      CONC_DIFF = 0, TIME_DIFF = 0, N_OBS = 1, HAS_PREV = 0
    ))
  }

  for (anchor in rev(occurrences)) {
    c0_candidates <- observations[
      observations$time <= anchor + 1e-8 & observations$time > anchor - interval - 1e-8,
      , drop = FALSE
    ]
    post_candidates <- observations[
      observations$time > anchor + 1e-8 & observations$time < anchor + interval - 1e-8,
      , drop = FALSE
    ]
    if (!nrow(c0_candidates) || !nrow(post_candidates)) next
    c0 <- c0_candidates[nrow(c0_candidates), , drop = FALSE]
    expected_post <- anchor + infusion_end + 1
    post <- post_candidates[["time"]]
    post <- post_candidates[which.min(abs(post - expected_post)), , drop = FALSE]
    c0_time <- interval - max(0, anchor - c0$time[[1]])
    post_time <- post$time[[1]] - anchor
    return(c(
      PREV_CONC = c0$concentration[[1]], PREV_TIME = c0_time,
      LAST_CONC = post$concentration[[1]], LAST_TIME = post_time,
      CONC_DIFF = post$concentration[[1]] - c0$concentration[[1]],
      TIME_DIFF = post$time[[1]] - c0$time[[1]],
      N_OBS = 2, HAS_PREV = 1
    ))
  }
  stop("le prédicteur AUC24 ML exige un C0 pré-dose suivi d'une concentration post-dose dans l'intervalle suivant")
}

ml_observation_values <- function(fit, sampling_protocol = list()) {
  observations <- fit$source_observations %||% data.frame()
  keep <- if (nrow(observations)) {
    is.finite(suppressWarnings(as.numeric(observations$time))) &
      is.finite(suppressWarnings(as.numeric(observations$concentration)))
  } else {
    logical()
  }
  observations <- observations[keep, c("time", "concentration"), drop = FALSE]
  if (!nrow(observations)) stop("ML AUC24 requires at least one measured concentration.")
  observations <- observations[order(observations$time), , drop = FALSE]
  strategy <- sampling_protocol$strategy %||% ""
  if (strategy %in% c("c0_only", "c0_plus_one_hour_post_infusion")) {
    return(ml_c0_post_observation_values(fit, observations, strategy))
  }
  last_observation_time <- observations$time[[nrow(observations)]]
  occurrences <- ml_dose_occurrences(fit$source_doses, last_observation_time)
  if (!length(occurrences)) stop("No dose precedes the concentration used by ML AUC24.")
  # A concentration sampled exactly at the next scheduled dose is a pre-dose
  # concentration and belongs to the dosing interval that just ended.
  interval_index <- findInterval(observations$time - 1e-8, occurrences)
  observations <- observations[interval_index > 0L, , drop = FALSE]
  interval_index <- interval_index[interval_index > 0L]
  if (!nrow(observations)) stop("No concentration follows an available dose for ML AUC24.")
  interval_counts <- table(interval_index)
  eligible_intervals <- as.integer(names(interval_counts)[interval_counts >= 2L])
  selected_interval <- if (length(eligible_intervals)) max(eligible_intervals) else max(interval_index)
  observations <- observations[interval_index == selected_interval, , drop = FALSE]
  anchor <- occurrences[[selected_interval]]
  observations$relative_time <- observations$time - anchor
  last <- observations[nrow(observations), , drop = FALSE]
  has_previous <- nrow(observations) >= 2L
  previous <- if (has_previous) observations[nrow(observations) - 1L, , drop = FALSE] else NULL
  c(
    PREV_CONC = if (has_previous) previous$concentration[[1]] else 0,
    PREV_TIME = if (has_previous) previous$relative_time[[1]] else 0,
    LAST_CONC = last$concentration[[1]],
    LAST_TIME = last$relative_time[[1]],
    CONC_DIFF = if (has_previous) last$concentration[[1]] - previous$concentration[[1]] else 0,
    TIME_DIFF = if (has_previous) last$relative_time[[1]] - previous$relative_time[[1]] else 0,
    N_OBS = nrow(observations),
    HAS_PREV = as.integer(has_previous)
  )
}

ml_population_profile <- function(fit, horizon = 24) {
  regimen <- ml_latest_regimen(fit)
  population_fit <- fit
  population_fit$estimate <- NULL
  population_fit$ml_eta_override <- NULL
  simulate_regimen(
    population_fit,
    dose = regimen[["DOSE"]],
    interval = regimen[["INTERVAL"]],
    infusion = regimen[["INFUSION"]],
    horizon = max(24, horizon),
    delta = 0.1
  )
}

ml_population_observation_values <- function(profile, observations) {
  predicted <- stats::approx(
    profile$time,
    profile$concentration,
    xout = c(observations[["PREV_TIME"]], observations[["LAST_TIME"]]),
    rule = 2,
    ties = "ordered"
  )$y
  ratios <- pmax(1e-4, pmin(1e4, c(
    observations[["PREV_CONC"]] / max(predicted[[1]], 1e-8),
    observations[["LAST_CONC"]] / max(predicted[[2]], 1e-8)
  )))
  c(
    PREV_POP_CONC = predicted[[1]],
    LAST_POP_CONC = predicted[[2]],
    PREV_CONC_RATIO = ratios[[1]],
    LAST_CONC_RATIO = ratios[[2]]
  )
}

ml_feature_value <- function(feature, fit, context = list()) {
  source <- feature$source %||% ""
  key <- feature$key %||% feature$name %||% ""
  if (!nzchar(key)) stop("ML feature key is missing.")
  value <- switch(
    source,
    covariate = fit$current_covariates[[key]] %||% NA_real_,
    eta = if (is.null(fit$estimate)) NA_real_ else fit$estimate$final_eta[[1]][[key]] %||% NA_real_,
    map_parameter = if (is.null(fit$estimate)) NA_real_ else tryCatch(
      as.numeric(mapbayr::get_param(fit$estimate, key))[[1]],
      error = function(error) NA_real_
    ),
    population_auc24 = context$population_auc24 %||% NA_real_,
    population_observation = (context$population_observations %||% numeric())[[key]] %||% NA_real_,
    regimen = ml_latest_regimen(fit)[[key]] %||% NA_real_,
    observation = (context$observations %||% ml_observation_values(fit))[[key]] %||% NA_real_,
    stop("Unsupported ML feature source: ", source)
  )
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value)) stop("ML feature is unavailable: ", key)
  value
}

ml_feature_row <- function(fit, artifact) {
  features <- artifact$featureSchema %||% list()
  if (!length(features)) stop("ML feature schema is empty.")
  feature_names <- vapply(features, function(feature) as.character(feature$name %||% ""), character(1))
  if (any(!nzchar(feature_names)) || anyDuplicated(feature_names)) stop("ML feature names are invalid.")
  sources <- vapply(features, function(feature) as.character(feature$source %||% ""), character(1))
  context <- list()
  if (any(sources %in% c("observation", "population_observation"))) {
    context$observations <- ml_observation_values(fit, artifact$samplingProtocol %||% list())
  }
  if (any(sources %in% c("population_auc24", "population_observation"))) {
    observations <- if ("population_observation" %in% sources) context$observations else NULL
    horizon <- if (is.null(observations)) 24 else max(24, observations[["PREV_TIME"]], observations[["LAST_TIME"]])
    profile <- ml_population_profile(fit, horizon)
    auc_profile <- profile[profile$time <= 24 + 1e-8, , drop = FALSE]
    context$population_auc24 <- trap_auc(auc_profile$time, auc_profile$concentration)
    if ("population_observation" %in% sources) {
      context$population_observations <- ml_population_observation_values(profile, observations)
    }
  }
  values <- vapply(features, ml_feature_value, numeric(1), fit = fit, context = context)
  names(values) <- feature_names
  list(
    values = values,
    data = as.data.frame(as.list(values), check.names = FALSE),
    matrix = matrix(values, nrow = 1L, dimnames = list(NULL, feature_names))
  )
}

resolve_ml_artifact_path <- function(relative_path) {
  relative_path <- as.character(relative_path %||% "")
  if (!nzchar(relative_path) || !grepl("\\.rds$", relative_path, ignore.case = TRUE)) {
    stop("ML artifact must be an RDS file.")
  }
  root <- normalizePath(ML_ROOT, winslash = "/", mustWork = TRUE)
  path <- normalizePath(file.path(ML_ROOT, relative_path), winslash = "/", mustWork = TRUE)
  if (!startsWith(path, paste0(root, "/"))) stop("ML artifact path escapes its controlled directory.")
  path
}

verified_ml_rds_path <- function(relative_path, expected_hash, label = "ML artifact") {
  path <- resolve_ml_artifact_path(relative_path)
  expected_hash <- tolower(as.character(expected_hash %||% ""))
  if (!nzchar(expected_hash) || !requireNamespace("digest", quietly = TRUE)) {
    stop(label, " fingerprint is unavailable.")
  }
  actual_hash <- digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
  if (!identical(expected_hash, tolower(actual_hash))) {
    stop(label, " fingerprint does not match the registry.")
  }
  path
}

build_ml_explanation <- function(artifact, booster, feature_row) {
  configuration <- artifact$explanation %||% list()
  if (!identical(configuration$type %||% "", "dalex_break_down")) {
    stop("DALEX explanation is not configured for this predictor.")
  }
  if (!requireNamespace("DALEX", quietly = TRUE)) stop("Package DALEX is unavailable.")
  background_path <- verified_ml_rds_path(
    configuration$backgroundPath,
    configuration$backgroundSha256,
    "DALEX background"
  )
  background <- readRDS(background_path)
  feature_names <- names(feature_row$values)
  if (!is.data.frame(background) || !all(feature_names %in% names(background))) {
    stop("DALEX background does not match the ML feature schema.")
  }
  background <- background[, feature_names, drop = FALSE]
  background[] <- lapply(background, function(value) suppressWarnings(as.numeric(value)))
  if (!nrow(background) || any(!vapply(background, function(value) all(is.finite(value)), logical(1)))) {
    stop("DALEX background contains invalid values.")
  }

  cache_key <- paste(
    artifact$id %||% "unnamed",
    artifact$artifactSha256 %||% "",
    configuration$backgroundSha256 %||% "",
    sep = ":"
  )
  explainer <- .ml_explainer_cache[[cache_key]]
  if (is.null(explainer)) {
    explainer <- DALEX::explain(
      model = booster,
      data = background,
      y = NULL,
      predict_function = function(model, newdata) {
        ml_prediction_value(stats::predict(model, as.matrix(newdata)), artifact, newdata)
      },
      label = "AUC24 ML",
      verbose = FALSE
    )
    .ml_explainer_cache[[cache_key]] <- explainer
  }
  parts <- as.data.frame(DALEX::predict_parts(
    explainer = explainer,
    new_observation = feature_row$data,
    type = "break_down"
  ))
  parts$variable_name <- as.character(parts$variable_name)
  baseline_row <- which(parts$variable_name == "intercept")
  contribution_rows <- which(nzchar(parts$variable_name) & parts$variable_name != "intercept")
  if (!length(baseline_row) || !length(contribution_rows)) stop("DALEX returned an incomplete explanation.")
  contributions <- data.frame(
    variable = parts$variable_name[contribution_rows],
    value = as.numeric(feature_row$values[parts$variable_name[contribution_rows]]),
    contribution = as.numeric(parts$contribution[contribution_rows]),
    stringsAsFactors = FALSE
  )
  prediction <- tail(parts$cumulative[is.finite(parts$cumulative)], 1L)
  list(
    available = TRUE,
    method = "DALEX break_down",
    baseline = as.numeric(parts$contribution[baseline_row[[1]]]),
    prediction = as.numeric(prediction),
    contributions = contributions,
    background_size = nrow(background),
    background_sha256 = configuration$backgroundSha256,
    synthetic = isTRUE(configuration$synthetic)
  )
}

ml_feature_domain_warnings <- function(values, artifact) {
  domain <- artifact$trainingDomain %||% list()
  bounds <- domain$features %||% list()
  warnings <- lapply(intersect(names(bounds), names(values)), function(name) {
    lower <- suppressWarnings(as.numeric(bounds[[name]]$min %||% NA_real_))
    upper <- suppressWarnings(as.numeric(bounds[[name]]$max %||% NA_real_))
    below <- is.finite(lower) && values[[name]] < lower
    above <- is.finite(upper) && values[[name]] > upper
    if (!below && !above) return(NULL)
    list(
      feature = name,
      value = as.numeric(values[[name]]),
      min = lower,
      max = upper,
      direction = if (below) "below" else "above",
      artifact_id = artifact$id %||% "unnamed"
    )
  })
  Filter(Negate(is.null), warnings)
}

apply_ml_artifact <- function(fit, artifact) {
  feature_row <- ml_feature_row(fit, artifact)
  values <- feature_row$values
  domain <- artifact$trainingDomain %||% list()
  domain_warnings <- ml_feature_domain_warnings(values, artifact)

  if (!requireNamespace("xgboost", quietly = TRUE)) stop("Package xgboost unavailable; MAP fallback used.")
  artifact_path <- verified_ml_rds_path(artifact$artifactPath, artifact$artifactSha256)
  booster <- readRDS(artifact_path)
  if (!inherits(booster, "xgb.Booster")) stop("ML artifact is not an xgboost booster.")
  prediction_value <- ml_prediction_value(stats::predict(booster, feature_row$matrix), artifact, feature_row$data)[[1]]
  prediction_type <- ml_prediction_type(artifact)

  if (identical(prediction_type, "auc24_direct")) {
    regimen <- ml_latest_regimen(fit)
    observation_values <- ml_observation_values(fit, artifact$samplingProtocol %||% list())
    minimum_observations <- suppressWarnings(as.integer(
      (artifact$samplingProtocol %||% list())$minimumObservations %||% 1L
    ))
    if (observation_values[["N_OBS"]] < minimum_observations) {
      strategy <- (artifact$samplingProtocol %||% list())$strategy %||% ""
      if (identical(strategy, "c0_plus_one_hour_post_infusion")) {
        stop("le prédicteur AUC24 ML exige un C0 pré-dose suivi d'une concentration post-dose dans l'intervalle suivant")
      }
      stop("le prédicteur AUC24 ML exige au moins ", minimum_observations, " concentrations dans un même intervalle posologique")
    }
    if (isTRUE((artifact$samplingProtocol %||% list())$steadyStateRequired) &&
        !ml_latest_regimen_is_steady_state(fit)) {
      stop("le prédicteur AUC24 ML exige une administration déclarée à l'état stationnaire (ss = 1)")
    }
    mode <- artifact$administrationMode %||% "intermittent"
    if (identical(mode, "intermittent") &&
        regimen[["INFUSION"]] >= regimen[["INTERVAL"]]) {
      stop("This ML AUC24 artifact only covers bolus or intermittent IV administrations.")
    }
    if (identical(mode, "continuous") && regimen[["INFUSION"]] < regimen[["INTERVAL"]]) {
      stop("This ML AUC24 artifact only covers continuous IV infusions.")
    }
    if (!is.finite(prediction_value) || prediction_value <= 0) stop("ML AUC24 prediction is invalid.")
    auc_domain <- domain$auc24 %||% list()
    auc_min <- suppressWarnings(as.numeric(auc_domain$min %||% NA_real_))
    auc_max <- suppressWarnings(as.numeric(auc_domain$max %||% NA_real_))
    if (is.finite(auc_min) && prediction_value < auc_min || is.finite(auc_max) && prediction_value > auc_max) {
      domain_warnings <- c(domain_warnings, list(list(
        feature = "ML_AUC24",
        value = as.numeric(prediction_value),
        min = auc_min,
        max = auc_max,
        direction = if (is.finite(auc_min) && prediction_value < auc_min) "below" else "above",
        artifact_id = artifact$id %||% "unnamed"
      )))
    }
    fit$ml_auc24 <- prediction_value
    fit$ml_features <- as.list(values)
    fit$ml_domain_warnings <- domain_warnings
    fit$ml_explanation <- tryCatch(
      build_ml_explanation(artifact, booster, feature_row),
      error = function(error) list(available = FALSE, reason = conditionMessage(error))
    )
    fit$ml_correction <- list(
      applied = TRUE,
      type = prediction_type,
      artifact_id = artifact$id %||% "unnamed",
      auc24 = prediction_value,
      unit = (artifact$prediction %||% list())$unit %||% "mg.h/L",
      extrapolated = length(domain_warnings) > 0L,
      research_validation = isTRUE(ml_artifact_eligibility(
        artifact,
        fit$id,
        model_record(fit$id)$drug[[1]],
        fit$contract$route,
        fit$contract$mode %||% ""
      )$research),
      clinical_validation = isTRUE(ml_artifact_eligibility(
        artifact,
        fit$id,
        model_record(fit$id)$drug[[1]],
        fit$contract$route,
        fit$contract$mode %||% ""
      )$clinical)
    )
    return(fit)
  }

  if (!identical(prediction_type, "eta_additive")) stop("Unsupported ML prediction type.")
  if (is.null(fit$estimate)) stop("ETA correction requires an individual MAP estimate.")
  correction <- artifact$correction %||% list()
  eta_name <- as.character(correction$eta %||% "")
  center <- fit$estimate$final_eta[[1]]
  if (!eta_name %in% names(center)) stop("ML correction references an unknown ETA: ", eta_name)
  maximum <- suppressWarnings(as.numeric(correction$maxAbsDelta %||% 0.5))
  if (!is.finite(maximum) || maximum <= 0 || maximum > 2) stop("ML correction bound is invalid.")
  bounded <- max(-maximum, min(maximum, prediction_value))
  center[[eta_name]] <- center[[eta_name]] + bounded
  fit$ml_eta_override <- center
  fit$ml_domain_warnings <- domain_warnings
  fit$ml_correction <- list(
    applied = TRUE,
    type = prediction_type,
    artifact_id = artifact$id %||% "unnamed",
    eta = eta_name,
    raw_delta = prediction_value,
    applied_delta = bounded,
    bounded = !isTRUE(all.equal(prediction_value, bounded)),
    extrapolated = length(domain_warnings) > 0L,
    clinical_validation = isTRUE(ml_artifact_eligibility(
      artifact,
      fit$id,
      model_record(fit$id)$drug[[1]],
      fit$contract$route,
      fit$contract$mode %||% ""
    )$clinical)
  )
  fit
}

apply_hybrid_ml_to_fits <- function(fits, route, enabled = FALSE) {
  output <- lapply(names(fits), function(id) {
    fit <- fits[[id]]
    if (inherits(fit, "tdm_fit_error") || !id %in% MODEL_CATALOG$id) return(fit)
    record <- model_record(id)
    artifacts <- compatible_ml_artifacts(id, record$drug[[1]], route, fit$contract$mode %||% "")
    if (!isTRUE(enabled)) {
      fit$ml_correction <- list(
        applied = FALSE,
        reason = if (length(artifacts)) "experimental_ml_disabled" else "no_compatible_artifact"
      )
      return(fit)
    }
    if (!length(artifacts)) {
      fit$ml_correction <- list(applied = FALSE, reason = "no_compatible_artifact")
      return(fit)
    }
    tryCatch(
      apply_ml_artifact(fit, artifacts[[1]]),
      error = function(error) {
        fit$ml_correction <- list(applied = FALSE, reason = conditionMessage(error))
        fit
      }
    )
  })
  stats::setNames(output, names(fits))
}

ml_application_summary <- function(fits, weights = NULL) {
  valid <- successful_fits(fits)
  if (!length(valid)) return(list(available = 0L, auc24 = NA_real_, clinical = FALSE, message = "Aucune estimation ML disponible."))
  applied <- vapply(valid, function(fit) isTRUE((fit$ml_correction %||% list())$applied), logical(1))
  direct <- vapply(valid, function(fit) identical((fit$ml_correction %||% list())$type %||% "", "auc24_direct"), logical(1))
  if (all(applied & direct)) {
    values <- vapply(valid, function(fit) as.numeric(fit$ml_auc24), numeric(1))
    units <- unique(vapply(valid, function(fit) as.character(fit$ml_correction$unit %||% "mg.h/L"), character(1)))
    if (length(units) != 1L) {
      return(list(
        available = length(valid),
        auc24 = NA_real_,
        clinical = FALSE,
        message = "Unités AUC incompatibles entre les artefacts ML sélectionnés."
      ))
    }
    if (is.null(weights) || !length(weights)) weights <- stats::setNames(rep(1 / length(valid), length(valid)), names(valid))
    selected_weights <- as.numeric(weights[names(valid)])
    selected_weights[!is.finite(selected_weights) | selected_weights < 0] <- 0
    if (sum(selected_weights) <= 0) selected_weights <- rep(1, length(valid))
    selected_weights <- selected_weights / sum(selected_weights)
    auc24 <- sum(values * selected_weights)
    clinical <- all(vapply(valid, function(fit) isTRUE(fit$ml_correction$clinical_validation), logical(1)))
    research <- all(vapply(valid, function(fit) isTRUE(fit$ml_correction$research_validation), logical(1)))
    labels <- vapply(valid, function(fit) fit$ml_correction$artifact_id, character(1))
    explanation <- if (length(valid) == 1L) {
      valid[[1]]$ml_explanation %||% list(available = FALSE, reason = "Explication DALEX indisponible.")
    } else {
      list(
        available = FALSE,
        reason = "L'explication locale n'est affichée que pour un prédicteur ML individuel, sans model averaging."
      )
    }
    domain_warnings <- unlist(
      lapply(valid, function(fit) fit$ml_domain_warnings %||% list()),
      recursive = FALSE,
      use.names = FALSE
    )
    return(list(
      available = length(valid),
      auc24 = auc24,
      unit = units[[1]],
      clinical = clinical,
      research = research,
      artifacts = labels,
      explanation = explanation,
      domain_warnings = domain_warnings,
      message = paste0(
        "AUC24 ML expérimentale : ", round(auc24, 1),
        " ", units[[1]], ". Les projections de doses restent calculées par MAP",
        if (length(valid) > 1L) " et model averaging." else ".",
        if (length(domain_warnings)) paste0(
          " Avertissement : extrapolation hors du domaine d'entraînement pour ",
          length(domain_warnings),
          " variable(s)."
        ) else "",
        if (!research) " Un ou plusieurs artefacts restent sous les seuils de performance internes et sont affichés uniquement comme expérimentaux." else ""
      )
    ))
  }
  if (any(applied)) {
    return(list(
      available = sum(applied),
      auc24 = NA_real_,
      clinical = FALSE,
      message = "Artefacts ML incomplets pour les modèles sélectionnés : estimation agrégée non affichée."
    ))
  }
  reasons <- unique(vapply(valid, function(fit) as.character((fit$ml_correction %||% list())$reason %||% ""), character(1)))
  reasons <- reasons[nzchar(reasons)]
  detailed_reasons <- setdiff(reasons, c("experimental_ml_disabled", "no_compatible_artifact"))
  if (length(detailed_reasons)) {
    return(list(
      available = 0L,
      auc24 = NA_real_,
      clinical = FALSE,
      message = paste0("AUC24 ML non appliquée : ", paste(detailed_reasons, collapse = " | "), ". Estimation MAP conservée.")
    ))
  }
  list(
    available = 0L,
    auc24 = NA_real_,
    clinical = FALSE,
    message = if ("experimental_ml_disabled" %in% reasons) {
      "AUC24 ML expérimentale non activée : estimation MAP conservée."
    } else {
      "Aucun prédicteur direct d'AUC24 évalué et compatible : estimation MAP conservée."
    }
  )
}

ml_status_summary <- function(model_ids = character(), route = "", administration_mode = "") {
  if (!length(model_ids)) {
    return(list(available = 0L, message = "ML : indisponible pour un modèle de session."))
  }
  counts <- vapply(model_ids, function(id) {
    record <- model_record(id)
    length(compatible_ml_artifacts(id, record$drug[[1]], route, administration_mode))
  }, integer(1))
  list(
    available = sum(counts),
    message = if (sum(counts)) {
      paste(sum(counts), "prédicteur(s) direct(s) d'AUC24 expérimental(aux) évalué(s) disponible(s).")
    } else {
      "Aucun prédicteur direct d'AUC24 évalué et compatible : estimation MAP conservée."
    }
  )
}
