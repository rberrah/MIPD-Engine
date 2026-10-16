suppressPackageStartupMessages({
  library(mrgsolve)
  library(mapbayr)
  library(dplyr)
  library(jsonlite)
})

`%||%` <- function(value, fallback) {
  if (is.null(value) || !length(value)) fallback else value
}

smoke_phase <- "startup"
if (identical(Sys.getenv("GITHUB_ACTIONS"), "true")) {
  options(error = function() {
    message <- paste0("Phase ", smoke_phase, ": ", geterrmessage())
    message <- gsub("%", "%25", message, fixed = TRUE)
    message <- gsub("\r", "%0D", message, fixed = TRUE)
    message <- gsub("\n", "%0A", message, fixed = TRUE)
    cat("::error title=TDM engine smoke test::", message, "\n", sep = "")
    quit(save = "no", status = 1L)
  })
}

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file_argument), ".."), winslash = "/", mustWork = TRUE)

source(file.path(APP_ROOT, "R", "model_library.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "i18n.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "clinical_presets.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "engine.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "ml_engine.R"), local = TRUE)

stopifnot(
  isTRUE(all.equal(time_above_threshold_pct(c(0, 1, 2), c(0, 2, 0), 1), 50)),
  identical(time_above_threshold_pct(c(0, 1), c(1, 1), 1), 0)
)
combined_scores <- score_target_metrics(
  data.frame(cmin = c(10, 4, 10), cmax = c(30, 30, 50)),
  "CminCmax",
  target_low = 8,
  target_high = 12,
  target_options = list(cmax_low = 25, cmax_high = 35)
)
stopifnot(
  identical(combined_scores$in_target, c(TRUE, FALSE, FALSE)),
  identical(combined_scores$value, combined_scores$distance),
  combined_scores$distance[[1]] == 0,
  all(combined_scores$distance[-1] > 0)
)
combined_operator_cases <- list(
  above_above = list(cmin = c(11, 9, 11), cmax = c(21, 21, 19), cmin_operator = "above", cmax_operator = "above"),
  above_below = list(cmin = c(11, 9, 11), cmax = c(19, 19, 21), cmin_operator = "above", cmax_operator = "below"),
  below_above = list(cmin = c(9, 11, 9), cmax = c(21, 21, 19), cmin_operator = "below", cmax_operator = "above"),
  below_below = list(cmin = c(9, 11, 9), cmax = c(19, 19, 21), cmin_operator = "below", cmax_operator = "below")
)
for (case in combined_operator_cases) {
  score <- score_target_metrics(
    data.frame(cmin = case$cmin, cmax = case$cmax),
    "CminCmax",
    target_low = 10,
    target_high = 10,
    target_options = list(
      cmin_operator = case$cmin_operator,
      cmax_operator = case$cmax_operator,
      cmax_low = 20,
      cmax_high = 20
    )
  )
  stopifnot(identical(score$in_target, c(TRUE, FALSE, FALSE)), score$distance[[1]] == 0, all(score$distance[-1] > 0))
}

stopifnot(
  identical(app_language_from_query("?model=vanco_roberts&lang=en"), "en"),
  identical(app_language_from_query("?lang=de"), "fr"),
  identical(app_t("en", "Analyse"), "Analysis"),
  identical(names(app_t("en", stats::setNames(c("Analyse", "Réglages"), c("analysis", "settings")))), c("analysis", "settings")),
  inherits(app_t("en", htmltools::HTML("<strong>test</strong>")), "html"),
  any(grepl("Vancomycin", names(catalog_choices_i18n(lang = "en")), fixed = TRUE))
)

smoke_phase <- "model compilation and contracts"
doses <- data.frame(time = 0, amount = 1000, interval = 12, count = 4, infusion = 1, ss = 0)
observations <- data.frame(time = 47.5, concentration = 18)
covariates <- list(WT = 70, AGE = 65, CREAT = 90, CREAT2 = 90, CRCL = 90, SEX = 0, HT = 175, DIAL = 0)
covariate_history <- data.frame(
  time = c(0, 47.5),
  WT = c(70, 74),
  AGE = c(65, 65),
  CREAT = c(90, 110),
  CREAT2 = c(90, 110),
  CRCL = c(90, 110),
  SEX = c(0, 0),
  HT = c(175, 175),
  DIAL = c(0, 0)
)
model_specification <- function(id, mode = NULL) {
  record <- model_record(id)
  if (is.null(mode)) mode <- model_administration_modes(record, "IV")[[1]]
  list(
    id = id,
    label = record$label[[1]],
    code = NULL,
    route = "IV",
    mode = model_administration_mode(record, "IV", mode),
    adm_cmt_name = model_administration_cmt(record, "IV")
  )
}
model_ids <- c("vanco_goti", "vanco_pkjust")
specifications <- lapply(model_ids, model_specification)
roberts_specification <- list(model_specification("vanco_roberts"))

mellon_record <- model_record("amox_mellon")
stopifnot(identical(model_routes(mellon_record), c("IV", "Oral")))
mellon_model <- compile_model(model_id = "amox_mellon")
mellon_contract <- validate_model_contract(mellon_model)
fixed_residual_model <- apply_residual_error_setting(mellon_model, "fixed_cv", 1)
fixed_log_model <- apply_residual_error_setting(compile_model(model_id = "amox_carlier"), "fixed_cv", 1.5)
stopifnot(
  isTRUE(all.equal(diag(as.matrix(mrgsolve::smat(fixed_residual_model))), c(0.0001, 0), tolerance = 1e-12)),
  isTRUE(all.equal(
    diag(as.matrix(mrgsolve::smat(fixed_log_model))),
    c(0, log1p(0.015^2)),
    tolerance = 1e-12
  )),
  identical(clinical_target_preset(model_record("vanco_roberts"), "IV_CONTINUOUS")$target_metric, "Cmin"),
  identical(clinical_target_preset(model_record("vanco_goti"), "IV_INTERMITTENT")$target_low, 400),
  is.null(clinical_target_preset(model_record("amox_carlier"), "IV_INTERMITTENT"))
)
stopifnot(resolve_administration_cmt(
  mellon_model,
  mellon_contract,
  list(label = "Mellon oral", adm_cmt_name = model_administration_cmt(mellon_record, "Oral"))
) == match("DEPOT", mellon_model@cmtL))

route_error <- tryCatch({
  validate_model_route_set(list(list(route = "IV"), list(route = "Oral")))
  NULL
}, error = identity)
stopifnot(inherits(route_error, "error"))
mode_error <- tryCatch({
  validate_model_route_set(list(
    list(route = "IV", mode = "IV_INTERMITTENT"),
    list(route = "IV", mode = "IV_CONTINUOUS")
  ))
  NULL
}, error = identity)
stopifnot(inherits(mode_error, "error"))

smoke_phase <- "MAP fits and residual error"
fits <- fit_model_set(
  specifications,
  doses,
  observations,
  covariates,
  allow_custom = FALSE,
  covariate_history = covariate_history
)
valid <- successful_fits(fits)
if (length(valid) != length(model_ids)) {
  messages <- vapply(fits, function(item) item$message %||% "OK", character(1))
  stop("Some vancomycin models failed: ", paste(messages, collapse = " | "))
}
stopifnot(identical(as.numeric(valid[[1]]$current_covariates$WT), 74))

fit_diagnostic <- fit_numeric_diagnostics(valid[[1]])
stopifnot(
  identical(fit_diagnostic$status, "fitted"),
  isTRUE(fit_diagnostic$converged),
  isTRUE(fit_diagnostic$sparse),
  isTRUE(fit_diagnostic$uncertainty_available),
  isFALSE(fit_diagnostic$at_boundary)
)
boundary_fit <- valid[[1]]
selected_eta <- boundary_fit$estimate$arg.optim$select_eta
boundary_fit$estimate$final_eta[[1]][selected_eta[[1]]] <- boundary_fit$estimate$arg.optim$lower[[1]]
stopifnot(isTRUE(fit_numeric_diagnostics(boundary_fit)$at_boundary))
no_uncertainty_fit <- valid[[1]]
no_uncertainty_fit$estimate$covariance[[1]] <- matrix(NA_real_)
stopifnot(isFALSE(fit_numeric_diagnostics(no_uncertainty_fit)$uncertainty_available))

fixed_residual_fits <- fit_model_set(
  specifications[1],
  doses,
  observations,
  covariates,
  allow_custom = FALSE,
  covariate_history = covariate_history,
  residual_error_mode = "fixed_cv",
  fixed_residual_cv = 1.5
)
fixed_residual_valid <- successful_fits(fixed_residual_fits)
if (length(fixed_residual_valid) != 1L) stop("The fixed residual-error MAP fit failed: ", fixed_residual_fits[[1]]$message %||% "unknown error")
stopifnot(
  identical(fixed_residual_valid[[1]]$residual_error$mode, "fixed_cv"),
  isTRUE(all.equal(diag(as.matrix(mrgsolve::smat(fixed_residual_valid[[1]]$model))), c(0.000225, 0), tolerance = 1e-12))
)

smoke_phase <- "steady-state dosing"
steady_doses <- data.frame(time = 36, amount = 1000, interval = 12, count = 1, infusion = 12, ss = 1)
steady_observations <- data.frame(time = 11.5, concentration = 18)
steady_fits <- fit_model_set(
  roberts_specification,
  steady_doses,
  steady_observations,
  covariates,
  allow_custom = FALSE,
  covariate_history = covariate_history
)
steady_valid <- successful_fits(steady_fits)
if (length(steady_valid) != 1L) stop("The steady-state MAP fit failed: ", steady_fits[[1]]$message %||% "unknown error")
steady_dose_row <- steady_valid[[1]]$data[steady_valid[[1]]$data$evid == 1, , drop = FALSE]
stopifnot(nrow(steady_dose_row) == 1L, steady_dose_row$time[[1]] == 0, steady_dose_row$ss[[1]] == 1L, steady_dose_row$addl[[1]] == 0)
steady_events <- dose_event_rows(steady_doses, steady_valid[[1]]$contract$adm_cmt, end_time = 36)
stopifnot(
  identical(as.numeric(steady_events$time), c(0, 12, 24, 36)),
  identical(as.integer(steady_events$ss), c(1L, 0L, 0L, 0L))
)
steady_profiles <- fit_profiles(steady_fits, c(vanco_roberts = 1), end_time = 36, delta = 0.25)
steady_cycle_concentrations <- stats::approx(
  steady_profiles$average$time,
  steady_profiles$average$concentration,
  xout = c(11.5, 35.5),
  ties = "ordered"
)$y
stopifnot(abs(diff(steady_cycle_concentrations)) < 0.1)
steady_exposure <- current_regimen_exposure(steady_fits, c(vanco_roberts = 1), steady_doses)
stopifnot(isTRUE(steady_exposure$steady_state), !isTRUE(steady_exposure$single_dose))
steady_context <- regimen_context(steady_doses, data.frame(time = 20, concentration = 18))
stopifnot(steady_context$decision_time == 20, steady_context$future_start == 24)
steady_comparison <- compare_future_regimens(
  steady_fits,
  c(vanco_roberts = 1),
  steady_doses,
  data.frame(time = 20, concentration = 18),
  recommended = data.frame(dose = 750, interval = 12, infusion = 12),
  additional_doses = 3,
  delta = 0.2
)
stopifnot(
  steady_comparison$future_start == 24,
  steady_comparison$additional_doses == 3L,
  max(steady_comparison$profiles$time) >= 60
)

smoke_phase <- "model averaging and exposure"
weights <- compute_model_weights(fits, scheme = "AIC")
stopifnot(length(weights) == length(model_ids), abs(sum(weights) - 1) < 1e-8)

long_interval <- simulate_averaged_regimen(
  fits,
  weights,
  dose = 1000,
  interval = 48,
  infusion = 1,
  delta = 0.2,
  concentration_threshold = 10
)
stopifnot(
  max(long_interval$profile$time) >= 48,
  all(is.finite(unlist(long_interval$metrics))),
  long_interval$metrics$time_above_pct >= 0,
  long_interval$metrics$time_above_pct <= 100
)

carried <- carry_covariates(
  event_times = c(0, 12, 47.5, 60),
  covariate_history = covariate_history[, c("time", "WT")],
  covariates = list(WT = 70)
)
stopifnot(identical(as.numeric(carried$WT), c(70, 70, 74, 74)))

profiles <- fit_profiles(fits, weights, end_time = 96)
stopifnot(nrow(profiles$per_model) > 0, nrow(profiles$average) > 0)
stopifnot(all(is.finite(profiles$average$concentration)))

summary_table <- model_summary(fits, weights)
stopifnot(nrow(summary_table) == length(model_ids))

current_exposure <- current_regimen_exposure(fits, weights, doses, observations)
stopifnot(is.finite(current_exposure$auc24), current_exposure$auc24 > 0)
stopifnot(is.finite(current_exposure$c0), current_exposure$c0 >= 0)
stopifnot(is.finite(current_exposure$steady_state_auc24), current_exposure$steady_state_auc24 > 0)
stopifnot(is.finite(current_exposure$steady_state_c0), current_exposure$steady_state_c0 >= 0)
stopifnot(identical(current_exposure$interval, 12))

smoke_phase <- "dose recommendations and PTA"
recommendations <- recommend_regimens(
  fits = fits,
  weights = weights,
  dose_min = 500,
  dose_max = 1500,
  dose_step = 500,
  intervals = c(8, 12, 24),
  infusion = 1,
  metric = "AUC24",
  target_low = 400,
  target_high = 600
)

recommendations <- rank_regimens_by_pta(
  recommendations,
  fits,
  weights,
  metric = "AUC24",
  target_low = 400,
  target_high = 600,
  replicates = 40,
  delta = 0.2,
  top_n = 3
)

stopifnot(nrow(recommendations) == 9)
stopifnot(all(is.finite(recommendations$auc24)))
stopifnot(all(recommendations$auc24 > 0))

time_recommendations <- recommend_regimens(
  fits = fits,
  weights = weights,
  dose_min = 500,
  dose_max = 1000,
  dose_step = 500,
  intervals = c(12, 24),
  infusion = 1,
  metric = "TimeAbove",
  target_low = 30,
  target_high = 80,
  target_options = list(concentration_threshold = 10)
)
time_recommendations <- rank_regimens_by_pta(
  time_recommendations,
  fits,
  weights,
  metric = "TimeAbove",
  target_low = 30,
  target_high = 80,
  target_options = list(concentration_threshold = 10),
  replicates = 20,
  delta = 0.2,
  top_n = 2
)
stopifnot(
  all(is.finite(time_recommendations$time_above_pct)),
  all(time_recommendations$time_above_pct >= 0 & time_recommendations$time_above_pct <= 100),
  all(is.finite(time_recommendations$p_target[time_recommendations$pta_evaluated]))
)

best <- recommendations[1, ]
comparison <- compare_future_regimens(
  fits,
  weights,
  doses,
  observations,
  recommended = best,
  additional_doses = 4,
  delta = 0.2
)
stopifnot(nrow(comparison$profiles) > 0)
stopifnot(identical(sort(unique(comparison$profiles$scenario)), sort(c(
  SCENARIO_MAINTAIN,
  SCENARIO_RECOMMENDED
))))

smoke_phase <- "uncertainty distribution"
distribution <- simulate_regimen_distribution(
  fits,
  weights,
  dose = best$dose,
  interval = best$interval,
  infusion = best$infusion,
  metric = "AUC24",
  replicates = 60,
  interval_level = 90,
  delta = 0.2,
  include_posterior = TRUE,
  include_residual = FALSE
)
stopifnot(nrow(distribution$data) == 60L)
stopifnot(all(is.finite(c(distribution$lower, distribution$median, distribution$upper))))
stopifnot(distribution$lower <= distribution$median, distribution$median <= distribution$upper)
stopifnot(distribution$lower < distribution$upper)

smoke_phase <- "model-averaging sensitivity"
sensitivity <- model_averaging_sensitivity(
  fits, weights,
  dose_min = 500, dose_max = 1500, dose_step = 500,
  intervals = c(8, 12, 24), infusion = 1,
  metric = "AUC24", target_low = 400, target_high = 600,
  delta = 0.2
)
stopifnot(nrow(sensitivity) >= 2L, all(is.finite(sensitivity$target_value)))

if (identical(Sys.getenv("TDM_SMOKE_CORE_ONLY"), "true")) {
  cat("TDM engine core smoke test OK (ML artifact checks skipped).\n")
  quit(save = "no", status = 0L)
}

ml_manifest <- read_ml_manifest()
invalid_ml_manifest <- ml_manifest
invalid_ml_manifest$artifacts[[1]]$prediction$unit <- NULL
stopifnot(inherits(try(validate_ml_manifest(invalid_ml_manifest), silent = TRUE), "try-error"))
artifact_for <- function(model_id, mode) {
  matches <- Filter(function(artifact) {
    identical(artifact$baseModelId, model_id) &&
      identical(ml_administration_mode(artifact$administrationMode), ml_administration_mode(mode))
  }, ml_manifest$artifacts)
  stopifnot(length(matches) == 1L)
  matches[[1]]
}
analysis_model_rows <- MODEL_CATALOG[vapply(seq_len(nrow(MODEL_CATALOG)), function(index) {
  model_analysis_eligible(MODEL_CATALOG[index, , drop = FALSE])
}, logical(1)), , drop = FALSE]
expected_ml_scopes <- unlist(lapply(seq_len(nrow(analysis_model_rows)), function(index) {
  record <- analysis_model_rows[index, , drop = FALSE]
  modes <- vapply(model_administration_modes(record), ml_administration_mode, character(1))
  modes <- modes[modes != "IV_CONTINUOUS"]
  if (!length(modes)) character() else paste(record$id[[1]], modes, sep = "::")
}), use.names = FALSE)
analysis_artifacts <- Filter(function(artifact) artifact$baseModelId %in% analysis_model_rows$id, ml_manifest$artifacts)
actual_ml_scopes <- vapply(analysis_artifacts, function(artifact) {
  paste(artifact$baseModelId, ml_administration_mode(artifact$administrationMode), sep = "::")
}, character(1))
published_artifact <- artifact_for("vanco_pkjust", "IV_INTERMITTENT")
stopifnot(
  identical(basename(ML_MANIFEST_PATH), "registry.json"),
  file.exists(ML_MANIFEST_PATH),
  identical(as.integer(ml_manifest$version), 3L),
  length(analysis_artifacts) == length(expected_ml_scopes),
  setequal(actual_ml_scopes, expected_ml_scopes),
  !anyDuplicated(actual_ml_scopes),
  published_artifact$trainingDomain$features$PREV_TIME$min <= 8,
  published_artifact$trainingDomain$features$PREV_TIME$max >= 8,
  abs(ml_prediction_value(
    log(2),
    list(prediction = list(transform = "exp_times_feature", scaleFeature = "POP_AUC24")),
    data.frame(POP_AUC24 = 50)
  ) - 100) < 1e-10
)
invisible(lapply(analysis_artifacts, function(artifact) {
  eligibility <- ml_artifact_eligibility(
    artifact,
    artifact$baseModelId,
    artifact$drug,
    artifact$route,
    artifact$administrationMode
  )
  if (!isTRUE(eligibility$experimental)) {
    stop(sprintf(
      "ML artifact %s is incompatible with its analysis model: %s",
      artifact$baseModelId,
      paste(eligibility$reasons, collapse = "; ")
    ))
  }
  stopifnot(
    identical(artifact$releaseLevel, if (isTRUE(eligibility$research)) "research" else "experimental"),
    is.finite(as.numeric(artifact$validation$populationBaseline$untouchedHoldout$relativeRmsePct)),
    is.finite(as.numeric(artifact$validation$populationBaseline$untouchedHoldout$relativeBiasPct)),
    is.finite(as.numeric(artifact$validation$comparison$untouchedHoldoutRelativeRmseGainPct))
  )
  verified_ml_rds_path(artifact$artifactPath, artifact$artifactSha256)
  verified_ml_rds_path(
    artifact$explanation$backgroundPath,
    artifact$explanation$backgroundSha256,
    "DALEX background"
  )
}))
published_eligibility <- ml_artifact_eligibility(
  published_artifact,
  "vanco_pkjust",
  "Vancomycine",
  "IV",
  "IV_INTERMITTENT"
)
stopifnot(
  isTRUE(published_eligibility$experimental),
  !isTRUE(published_eligibility$clinical),
  identical(
    published_artifact$releaseLevel,
    if (isTRUE(published_eligibility$research)) "research" else "experimental"
  )
)
experimental_contract <- published_artifact
experimental_contract$validation$repeatedCv$passed <- FALSE
experimental_contract$validation$untouchedHoldout$passed <- FALSE
experimental_eligibility <- ml_artifact_eligibility(
  experimental_contract,
  "vanco_pkjust",
  "Vancomycine",
  "IV",
  "IV_INTERMITTENT"
)
stopifnot(isTRUE(experimental_eligibility$experimental), !isTRUE(experimental_eligibility$research))
published_explanation <- published_artifact$explanation
stopifnot(
  identical(published_explanation$type, "dalex_break_down"),
  isTRUE(published_explanation$synthetic),
  published_explanation$sampleSize >= 100L
)
invisible(verified_ml_rds_path(
  published_explanation$backgroundPath,
  published_explanation$backgroundSha256,
  "DALEX background"
))

revilla_observations <- data.frame(
  time = c(36, 38),
  concentration = c(18, 30)
)
revilla_doses <- doses
revilla_doses$time <- 36
revilla_doses$count <- 1L
revilla_doses$ss <- 1L
revilla_fits <- fit_model_set(
  specifications[match("vanco_pkjust", model_ids)],
  revilla_doses,
  revilla_observations,
  covariates,
  allow_custom = FALSE,
  covariate_history = covariate_history
)
revilla_fits <- apply_hybrid_ml_to_fits(revilla_fits, "IV", enabled = TRUE)
stopifnot(
  isTRUE(revilla_fits[[1]]$ml_correction$applied),
  identical(revilla_fits[[1]]$ml_correction$type, "auc24_direct"),
  revilla_fits[[1]]$ml_auc24 >= 100,
  revilla_fits[[1]]$ml_auc24 <= 1200,
  isTRUE(revilla_fits[[1]]$ml_explanation$available),
  abs(
    revilla_fits[[1]]$ml_explanation$baseline +
      sum(revilla_fits[[1]]$ml_explanation$contributions$contribution) -
      revilla_fits[[1]]$ml_auc24
  ) < 1e-3
)
boundary_fits <- fit_model_set(
  specifications[match("vanco_pkjust", model_ids)],
  data.frame(time = 0, amount = 1000, interval = 12, count = 1, infusion = 1, ss = 1),
  data.frame(time = c(0, 2), concentration = c(15, 30)),
  list(WT = 73, AGE = 61.1, CREAT = 123.8, CRCL = 74.7),
  allow_custom = FALSE,
  covariate_history = data.frame(
    time = c(0, 2), WT = c(74.8, 73), AGE = 61.1,
    CREAT = 123.8, CRCL = c(90.7, 74.7)
  )
)
boundary_fits <- apply_hybrid_ml_to_fits(boundary_fits, "IV", enabled = TRUE)
boundary_ml_summary <- ml_application_summary(boundary_fits, c(vanco_pkjust = 1))
stopifnot(
  boundary_ml_summary$available == 1L,
  is.finite(boundary_ml_summary$auc24),
  !length(boundary_ml_summary$domain_warnings),
  identical(unname(boundary_ml_summary$artifacts), published_artifact$id)
)
revilla_ml_summary <- ml_application_summary(revilla_fits, c(vanco_pkjust = 1))
stopifnot(
  revilla_ml_summary$available == 1L,
  is.finite(revilla_ml_summary$auc24),
  isTRUE(revilla_ml_summary$explanation$available)
)
stopifnot(
  length(compatible_ml_artifacts("vanco_pkjust", "Vancomycine", "IV", "IV_INTERMITTENT")) == 1L,
  length(compatible_ml_artifacts("vanco_pkjust", "Vancomycine", "IV", "IV_CONTINUOUS")) == 0L
)

averaging_fits <- fit_model_set(
  specifications,
  revilla_doses,
  revilla_observations,
  covariates,
  allow_custom = FALSE,
  covariate_history = covariate_history
)
averaging_weights <- compute_model_weights(averaging_fits, "AIC")
averaging_fits <- apply_hybrid_ml_to_fits(averaging_fits, "IV", enabled = TRUE)
averaging_ml_summary <- ml_application_summary(averaging_fits, averaging_weights)
stopifnot(
  averaging_ml_summary$available == 2L,
  is.finite(averaging_ml_summary$auc24),
  length(averaging_ml_summary$artifacts) == 2L,
  !isTRUE(averaging_ml_summary$explanation$available),
  !grepl("Artefacts ML incomplets", averaging_ml_summary$message, fixed = TRUE)
)

incomplete_averaging_fits <- averaging_fits
incomplete_averaging_fits[[2]]$ml_correction <- list(
  applied = FALSE,
  reason = "no_compatible_artifact"
)
incomplete_averaging_fits[[2]]$ml_auc24 <- NULL
incomplete_averaging_summary <- ml_application_summary(incomplete_averaging_fits, averaging_weights)
stopifnot(
  incomplete_averaging_summary$available == 1L,
  !is.finite(incomplete_averaging_summary$auc24),
  grepl("Artefacts ML incomplets", incomplete_averaging_summary$message, fixed = TRUE)
)

default_fits <- fit_model_set(
  roberts_specification,
  doses,
  observations,
  list(WT = 74.8, CRCL = 90.7),
  allow_custom = FALSE,
  covariate_history = data.frame(time = 47.5, WT = 74.8, CRCL = 90.7)
)
default_weights <- compute_model_weights(default_fits, "AIC")
default_recommendations <- recommend_regimens(
  default_fits, default_weights,
  dose_min = 250, dose_max = 2000, dose_step = 250,
  intervals = c(8, 12, 24), infusion = 1,
  metric = "AUC24", target_low = 400, target_high = 600,
  delta = 0.1
)
default_recommendations <- rank_regimens_by_pta(
  default_recommendations, default_fits, default_weights,
  metric = "AUC24", target_low = 400, target_high = 600,
  replicates = 150, delta = 0.1
)
default_best <- default_recommendations[1, , drop = FALSE]
default_sensitivity <- model_averaging_sensitivity(
  default_fits, default_weights,
  dose_min = 250, dose_max = 2000, dose_step = 250,
  intervals = c(8, 12, 24), infusion = 1,
  metric = "AUC24", target_low = 400, target_high = 600,
  delta = 0.1
)
default_distribution <- simulate_regimen_distribution(
  default_fits, default_weights,
  dose = default_best$dose, interval = default_best$interval, infusion = 1,
  metric = "AUC24", replicates = 250, interval_level = 90, delta = 0.1,
  include_posterior = TRUE, include_residual = FALSE, include_timing = FALSE
)
stopifnot(nrow(default_sensitivity) == 1L, nrow(default_distribution$data) == 250L)
stopifnot(isTRUE(default_distribution$posterior_available))

default_fits <- apply_hybrid_ml_to_fits(default_fits, "IV")
stopifnot(!isTRUE(default_fits[[1]]$ml_correction$applied))
stopifnot(identical(default_fits[[1]]$ml_correction$reason, "no_compatible_artifact"))

runtime_regimen <- ml_latest_regimen(default_fits[[1]])
runtime_observations <- ml_observation_values(default_fits[[1]])
stopifnot(
  identical(unname(runtime_regimen[["DOSE"]]), 1000),
  identical(unname(runtime_regimen[["INTERVAL"]]), 12),
  identical(unname(runtime_regimen[["INFUSION"]]), 1),
  identical(unname(runtime_observations[["LAST_CONC"]]), 18),
  abs(runtime_observations[["LAST_TIME"]] - 11.5) < 1e-8,
  identical(unname(runtime_observations[["HAS_PREV"]]), 0)
)

long_interval_fit <- default_fits[[1]]
long_interval_fit$source_doses <- data.frame(
  time = 0, amount = 1000, interval = 48, count = 1, infusion = 1, ss = 1
)
long_interval_fit$source_observations <- data.frame(
  time = c(8, 36), concentration = c(18, 6)
)
long_interval_features <- ml_feature_row(long_interval_fit, list(featureSchema = list(
  list(name = "LAST_POP_CONC", source = "population_observation", key = "LAST_POP_CONC")
)))
long_interval_profile <- ml_population_profile(long_interval_fit, 36)
expected_long_concentration <- stats::approx(
  long_interval_profile$time,
  long_interval_profile$concentration,
  xout = 36,
  ties = "ordered"
)$y
stopifnot(
  max(long_interval_profile$time) >= 36,
  abs(long_interval_features$values[["LAST_POP_CONC"]] - expected_long_concentration) < 1e-8
)

same_interval_fit <- default_fits[[1]]
same_interval_fit$source_doses <- data.frame(time = 0, amount = 1000, interval = 12, count = 3, infusion = 1, ss = 0)
same_interval_fit$source_observations <- data.frame(time = c(3, 8, 25), concentration = c(31, 18, 12))
same_interval_values <- ml_observation_values(same_interval_fit)
stopifnot(
  same_interval_values[["N_OBS"]] == 2,
  same_interval_values[["PREV_TIME"]] == 3,
  same_interval_values[["LAST_TIME"]] == 8,
  same_interval_values[["LAST_CONC"]] == 18
)
split_interval_fit <- same_interval_fit
split_interval_fit$source_observations <- data.frame(time = c(3, 15), concentration = c(31, 18))
stopifnot(ml_observation_values(split_interval_fit)[["N_OBS"]] == 1)

steady_boundary_fit <- same_interval_fit
steady_boundary_fit$source_doses <- data.frame(time = 0, amount = 1000, interval = 12, count = 1, infusion = 1, ss = 1)
steady_boundary_fit$source_observations <- data.frame(time = c(8, 12), concentration = c(18, 15))
steady_boundary_values <- ml_observation_values(steady_boundary_fit)
stopifnot(
  steady_boundary_values[["N_OBS"]] == 2,
  steady_boundary_values[["PREV_TIME"]] == 8,
  steady_boundary_values[["LAST_TIME"]] == 12,
  steady_boundary_values[["LAST_CONC"]] == 15
)

c0_post_fit <- steady_boundary_fit
c0_post_fit$source_observations <- data.frame(time = c(0, 2), concentration = c(15, 31))
c0_post_values <- ml_observation_values(c0_post_fit, list(strategy = "c0_plus_one_hour_post_infusion"))
c0_only_fit <- steady_boundary_fit
c0_only_fit$source_observations <- data.frame(time = 0, concentration = 15)
c0_only_values <- ml_observation_values(c0_only_fit, list(strategy = "c0_only"))
stopifnot(
  c0_post_values[["N_OBS"]] == 2,
  c0_post_values[["PREV_CONC"]] == 15,
  c0_post_values[["PREV_TIME"]] == 12,
  c0_post_values[["LAST_CONC"]] == 31,
  c0_post_values[["LAST_TIME"]] == 2,
  c0_post_values[["TIME_DIFF"]] == 2,
  c0_only_values[["N_OBS"]] == 1,
  c0_only_values[["LAST_CONC"]] == 15,
  c0_only_values[["LAST_TIME"]] == 12
)

hash <- model_sha256("vanco_roberts")
if (!is.na(hash)) {
  eligible_artifact <- list(
    id = "test-artifact",
    baseModelId = "vanco_roberts",
    baseModelSha256 = hash,
    drug = "Vancomycine",
    route = "IV",
    administrationMode = "intermittent",
    prediction = list(type = "auc24_direct", metric = "AUC24", unit = "mg.h/L"),
    validation = list(
      repeatedCv = list(passed = TRUE, relativeRmsePct = 10, relativeBiasPct = 1, within20Pct = 90),
      untouchedHoldout = list(passed = TRUE, relativeRmsePct = 11, relativeBiasPct = 2, within20Pct = 88),
      alternatePopPk = list(passed = TRUE, relativeRmsePct = 12, relativeBiasPct = 2, within20Pct = 85),
      realPatient = list(status = "pending")
    )
  )
  eligibility <- ml_artifact_eligibility(eligible_artifact, "vanco_roberts", "Vancomycine", "IV", "IV_INTERMITTENT")
  stopifnot(isTRUE(eligibility$research), !isTRUE(eligibility$clinical))
  stopifnot(!isTRUE(ml_artifact_eligibility(eligible_artifact, "vanco_roberts", "Vancomycine", "IV", "IV_CONTINUOUS")$research))
  eligible_artifact$validation$untouchedHoldout$passed <- FALSE
  stopifnot(!isTRUE(ml_artifact_eligibility(eligible_artifact, "vanco_roberts", "Vancomycine", "IV")$research))

  test_ml_root <- tempfile("tdm-ml-")
  dir.create(test_ml_root)
  training_matrix <- matrix(
    c(500, 10, 1000, 18, 1500, 24, 2000, 30),
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("DOSE", "LAST_CONC"))
  )
  training_target <- c(300, 420, 550, 700)
  booster <- xgboost::xgb.train(
    params = list(objective = "reg:squarederror", max_depth = 2, eta = 0.3, min_child_weight = 1, nthread = 1),
    data = xgboost::xgb.DMatrix(training_matrix, label = training_target),
    nrounds = 8,
    verbose = 0
  )
  saveRDS(booster, file.path(test_ml_root, "direct-auc24.rds"))
  old_ml_root <- ML_ROOT
  ML_ROOT <- test_ml_root
  eligible_artifact$artifactPath <- "direct-auc24.rds"
  eligible_artifact$artifactSha256 <- digest::digest(
    file.path(test_ml_root, "direct-auc24.rds"),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  )
  eligible_artifact$featureSchema <- list(
    list(name = "DOSE", source = "regimen", key = "DOSE"),
    list(name = "LAST_CONC", source = "observation", key = "LAST_CONC")
  )
  eligible_artifact$trainingDomain <- list(
    auc24 = list(min = 1, max = 2000),
    features = list(
      DOSE = list(min = 500, max = 2000),
      LAST_CONC = list(min = 0.5, max = 100)
    )
  )
  eligible_artifact$validation$untouchedHoldout$passed <- TRUE
  ml_fit <- apply_ml_artifact(default_fits[[1]], eligible_artifact)
  extrapolation_artifact <- eligible_artifact
  extrapolation_artifact$trainingDomain$features$DOSE$max <- 100
  extrapolated_ml_fit <- apply_ml_artifact(default_fits[[1]], extrapolation_artifact)
  ML_ROOT <- old_ml_root
  unlink(test_ml_root, recursive = TRUE, force = TRUE)
  stopifnot(
    isTRUE(ml_fit$ml_correction$applied),
    identical(ml_fit$ml_correction$type, "auc24_direct"),
    is.finite(ml_fit$ml_auc24),
    ml_fit$ml_auc24 > 0,
    !length(ml_fit$ml_eta_override %||% numeric())
  )
  stopifnot(
    isTRUE(extrapolated_ml_fit$ml_correction$applied),
    isTRUE(extrapolated_ml_fit$ml_correction$extrapolated),
    length(extrapolated_ml_fit$ml_domain_warnings) == 1L,
    identical(extrapolated_ml_fit$ml_domain_warnings[[1]]$feature, "DOSE")
  )
  ml_summary <- ml_application_summary(list(vanco_roberts = ml_fit), c(vanco_roberts = 1))
  stopifnot(ml_summary$available == 1L, is.finite(ml_summary$auc24), !isTRUE(ml_summary$clinical))
  extrapolated_summary <- ml_application_summary(list(vanco_roberts = extrapolated_ml_fit), c(vanco_roberts = 1))
  stopifnot(
    extrapolated_summary$available == 1L,
    is.finite(extrapolated_summary$auc24),
    length(extrapolated_summary$domain_warnings) == 1L,
    grepl("extrapolation", extrapolated_summary$message, fixed = TRUE)
  )
}

cat("TDM engine smoke test OK\n")
cat("Weights:", paste(names(weights), round(weights, 4), collapse = " | "), "\n")
cat("Current exposure: AUC0-24", round(current_exposure$auc24, 1), "; C0", round(current_exposure$c0, 2), "\n")
cat("Best tested scenario:", best$dose, "mg every", best$interval, "h; AUC24", round(best$auc24, 1), "\n")
