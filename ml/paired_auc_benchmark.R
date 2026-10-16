# Sourced by train_models_xgboost.R --paired-benchmark. No clinical data are used.
PAIRED_RESIDUAL_SD_SCALE <- 0.10
PAIRED_AUC_TAIL_PROBABILITY <- 0.05

# Preserve each article's residual-error structure while reducing every
# residual standard deviation ten-fold for this sparse-sampling experiment.
observation_model <- function(model) {
  sigma <- as.matrix(mrgsolve::smat(model))
  if (!nrow(sigma)) return(model)
  mrgsolve::smat(model, sigma * PAIRED_RESIDUAL_SD_SCALE^2)
}

# In the Woillard model ST=1 is immediate-release Prograf and ST=0 is
# once-daily prolonged-release Advagraf. Keep formulation and interval paired.
align_regimen_covariates <- function(base_scope, generator_scope, covariates, regimens) {
  if (identical(base_scope$model_id[[1]], "tacrolimus_woillard_ddi") && "ST" %in% names(covariates)) {
    covariates$ST <- ifelse(regimens$interval == 12, 1, 0)
  }
  covariates
}

central_auc_filter <- function(data, tail_probability = PAIRED_AUC_TAIL_PROBABILITY) {
  bounds <- unname(stats::quantile(data$TRUE_AUC24, c(tail_probability, 1 - tail_probability), type = 7))
  keep <- data$TRUE_AUC24 >= bounds[[1]] & data$TRUE_AUC24 <= bounds[[2]]
  list(data = data[keep, , drop = FALSE], bounds = bounds, excluded = sum(!keep))
}

paired_metrics <- function(truth, prediction) {
  stopifnot(length(truth) == length(prediction), all(is.finite(truth) & truth > 0))
  valid <- is.finite(prediction) & prediction > 0
  relative <- (prediction[valid] - truth[valid]) / truth[valid]
  list(n = length(truth), nEstimated = sum(valid), nFailed = sum(!valid),
       relativeBiasPct = if (any(valid)) 100 * mean(relative) else NULL,
       relativeRmsePct = if (any(valid)) 100 * sqrt(mean(relative^2)) else NULL,
       within20Pct = 100 * sum(abs(relative) <= 0.2 + 1e-12) / length(truth))
}

paired_events <- function(data, model, scope) {
  continuous <- scope$mode[[1]] == "IV_CONTINUOUS"
  data.frame(ID = seq_len(nrow(data)), time = 0, evid = 1,
    amt = if (continuous) 0 else data$DOSE,
    ii = if (continuous) 0 else data$INTERVAL,
    addl = if (continuous) 0 else pmax(0, floor((24 - 1e-8) / data$INTERVAL)),
    ss = 1, cmt = match(model_administration_cmt(model_record(scope$model_id[[1]]), scope$route[[1]]), model@cmtL),
    rate = if (continuous) data$DOSE / data$INTERVAL else ifelse(data$INFUSION > 0, data$DOSE / data$INFUSION, 0),
    DV = NA_real_, mdv = 1)
}

paired_map_model <- function(scope) {
  model <- observation_model(compiled_model(scope$model_id[[1]]))
  reserved <- c("IPRED", "PRED", grep("^ETA[0-9]+$", model@capL, value = TRUE))
  model <- mrgsolve::update(model, outvars = c(model@cmtL, setdiff(model@capL, reserved)))
  sigma <- diag(as.matrix(mrgsolve::smat(model)))
  if (isTRUE(getFromNamespace("log_transformation", "mapbayr")(model))) {
    if (length(sigma) == 2L && sigma[[1]] > 0 && sigma[[2]] == 0) {
      # Exact lognormal likelihood: mapbayr expects log-scale variance in slot 2.
      # This fitting object is never used to generate noisy observations.
      model <- mrgsolve::smat(model, diag(c(0, sigma[[1]])))
    } else if (length(sigma) == 2L && all(sigma > 0)) {
      stop("Exponential plus additive residual error is not supported exactly by mapbayr; no Gaussian approximation applied.")
    }
  }
  mapbayr::check_mapbayr_model(model)
  model
}

paired_map <- function(data, scope, two_points) {
  model <- paired_map_model(scope)
  continuous <- scope$mode[[1]] == "IV_CONTINUOUS"
  model <- mrgsolve::update(model, ss_rtol = if (continuous) 1e-6 else 1e-8,
                            ss_atol = if (continuous) 1e-6 else 1e-8)
  events <- paired_events(data, model, scope)
  observations <- events
  observations$time <- if (continuous) 0 else data$INTERVAL - 1e-6
  observations$DV <- data$PREV_CONC
  observations$evid <- observations$amt <- observations$ii <- observations$addl <-
    observations$ss <- observations$rate <- observations$mdv <- 0
  observations$cmt <- tagged_compartment(model, "OBS")
  if (two_points) {
    second <- observations
    second$time <- data$LAST_TIME
    second$DV <- data$LAST_CONC
    observations <- rbind(observations, second)
  }
  fit_data <- rbind(events, observations)
  fit_data <- fit_data[order(fit_data$ID, fit_data$time, -fit_data$evid), ]
  cov_names <- intersect(names(model_covariate_defaults(scope$model_id[[1]])), names(data))
  for (name in cov_names) fit_data[[name]] <- data[[name]][fit_data$ID]
  # Only observed concentrations and known covariates enter MAP, never generating ETAs.
  estimate <- mapbayr::mapbayest(model, data = fit_data, hessian = NULL,
                                select_eta = which(diag(as.matrix(mrgsolve::omat(model))) > 0),
                                verbose = FALSE, progress = FALSE)
  etas <- as.data.frame(mapbayr::get_eta(estimate))
  if ("ID" %in% names(etas)) etas <- etas[match(seq_len(nrow(data)), etas$ID), setdiff(names(etas), "ID"), drop = FALSE]
  stopifnot(nrow(etas) == nrow(data))
  names(etas) <- paste0("ETA", seq_len(ncol(etas)))
  individual <- cbind(ID = seq_len(nrow(data)), data[, cov_names, drop = FALSE], etas)
  out <- model |> mrgsolve::zero_re() |> mrgsolve::idata_set(individual) |>
    mrgsolve::data_set(events) |> mrgsolve::mrgsim(end = 24, delta = 0.1, recsort = 3,
      ss_n = if (continuous) 5000 else 500) |> as.data.frame()
  auc <- vapply(seq_len(nrow(data)), function(id) {
    profile <- out[out$ID == id, c("time", "DV")]
    profile <- profile[!duplicated(profile$time, fromLast = TRUE), ]
    trap_auc_vector(profile$time, pmax(0, profile$DV))
  }, numeric(1))
  if ("convergence" %in% names(estimate$opt.value)) {
    failed <- as.integer(estimate$opt.value$ID[estimate$opt.value$convergence != 0])
    auc[failed] <- NA_real_
  }
  auc
}

safe_paired_map <- function(data, scope, two_points) {
  tryCatch(paired_map(data, scope, two_points), error = function(error) {
    if (nrow(data) == 1L) {
      message(scope$model_id, " MAP failed: ", conditionMessage(error))
      return(NA_real_)
    }
    middle <- floor(nrow(data) / 2)
    c(safe_paired_map(data[seq_len(middle), , drop = FALSE], scope, two_points),
      safe_paired_map(data[seq.int(middle + 1L, nrow(data)), , drop = FALSE], scope, two_points))
  })
}

paired_features <- function(data, two_points) {
  exclude <- c("patient_id", "TRUE_AUC24", "POP_AUC24", "PREV_POP_CONC", "LAST_POP_CONC",
               "PREV_CONC_RATIO", "LAST_CONC_RATIO", "N_OBS", "HAS_PREV")
  if (!two_points) exclude <- c(exclude, "LAST_CONC", "LAST_TIME", "CONC_DIFF", "TIME_DIFF")
  names(data)[!names(data) %in% exclude]
}

paired_xgboost <- function(training, testing, two_points, seed_value) {
  features <- paired_features(training, two_points)
  dtrain <- xgboost::xgb.DMatrix(as.matrix(training[, features]), label = training$TRUE_AUC24)
  set.seed(seed_value)
  fold <- fold_ids(nrow(training), if (smoke) 3 else 10, seed_value)
  folds <- lapply(sort(unique(fold)), function(id) which(fold == id))
  grid <- if (smoke) {
    data.frame(max_depth = 1L, eta = 0.0261, min_child_weight = 1)
  } else {
    expand.grid(max_depth = c(1L, 2L, 4L), eta = c(0.0261, 0.05), min_child_weight = 40)
  }
  mtry <- min(7L, length(features))
  colsample <- mtry / length(features)
  best <- NULL
  for (i in seq_len(nrow(grid))) {
    parameters <- list(objective = "reg:squarederror", eval_metric = "rmse", nthread = 1,
      max_depth = grid$max_depth[[i]], eta = grid$eta[[i]], min_child_weight = grid$min_child_weight[[i]],
      subsample = 1, colsample_bynode = colsample, base_score = mean(training$TRUE_AUC24), seed = seed_value)
    cv <- xgboost::xgb.cv(params = parameters, data = dtrain, folds = folds,
                         nrounds = if (smoke) 50 else 1000, early_stopping_rounds = 30, verbose = 0)
    iteration <- which.min(cv$evaluation_log$test_rmse_mean)
    score <- cv$evaluation_log$test_rmse_mean[[iteration]]
    if (is.null(best) || score < best$rmse) best <- list(params = parameters, rounds = iteration, rmse = score)
  }
  booster <- xgboost::xgb.train(params = best$params, data = dtrain, nrounds = best$rounds, verbose = 0)
  list(prediction = as.numeric(predict(booster, as.matrix(testing[, features]))),
       selection = best, features = features)
}

paired_scope <- function(index) {
  scope <- scopes[index, , drop = FALSE]
  key <- paste(scope$model_id, scope$mode, sep = "-")
  checkpoint <- file.path(paired_directory, paste0(key, ".rds"))
  if (file.exists(checkpoint)) return(readRDS(checkpoint))
  local_seed <- seed + sum(utf8ToInt(key) * seq_along(utf8ToInt(key)))
  set.seed(local_seed)
  invisible(compiled_model(scope$model_id[[1]]))
  unavailable <- tryCatch({ paired_map_model(scope); NULL }, error = function(error) conditionMessage(error))
  # The benchmark uses the complete normal ETA prior, without clipping its tails.
  sample_eta_matrix <- function(model, n) {
    omega <- as.matrix(mrgsolve::omat(model))
    if (!nrow(omega)) return(data.frame(row.names = seq_len(n)))
    eig <- eigen((omega + t(omega)) / 2, symmetric = TRUE)
    root <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = nrow(omega))
    values <- matrix(rnorm(n * nrow(omega)), nrow = n) %*% t(root)
    colnames(values) <- paste0("ETA", seq_len(ncol(values)))
    as.data.frame(values)
  }
  # Simulation helpers resolve these names in their shared source environment.
  assign("sample_eta_matrix", sample_eta_matrix, envir = environment(simulate_batch))
  assign("MIN_AUC_RATIO", 0, envir = environment(simulate_batch))
  assign("MAX_AUC_RATIO", Inf, envir = environment(simulate_batch))
  scope_n <- if (!smoke && !identical(scope$model_id[[1]], "tacrolimus_woillard_ddi")) min(n_patients, 3000L) else n_patients
  simulated <- make_cohort(scope, scope, scope_n, "paired MAP/ML")
  filtered <- central_auc_filter(simulated)
  data <- filtered$data
  set.seed(local_seed + 1L)
  train_ids <- sample.int(nrow(data), floor(0.75 * nrow(data)))
  training <- data[train_ids, ]
  holdout_pool <- data[-train_ids, ]
  evaluation_n <- if (identical(scope$model_id[[1]], "tacrolimus_woillard_ddi")) 1098L else 300L
  set.seed(local_seed + 2L)
  testing <- holdout_pool[sample.int(nrow(holdout_pool), min(evaluation_n, nrow(holdout_pool))), , drop = FALSE]
  results <- list()
  for (two_points in c(FALSE, TRUE)) {
    design <- if (two_points) "c0PlusOneHourPostInfusion" else "c0Only"
    cat(key, " ", design, " MAP and XGBoost\n")
    chunks <- split(seq_len(nrow(testing)), ceiling(seq_len(nrow(testing)) / 50))
    map <- if (!is.null(unavailable)) rep(NA_real_, nrow(testing)) else
      unlist(lapply(chunks, function(ids) safe_paired_map(testing[ids, , drop = FALSE], scope, two_points)), use.names = FALSE)
    ml <- paired_xgboost(training, testing, two_points, local_seed + 3L)
    paired <- is.finite(map) & map > 0 & is.finite(ml$prediction) & ml$prediction > 0
    results[[design]] <- list(map = paired_metrics(testing$TRUE_AUC24, map),
      ml = paired_metrics(testing$TRUE_AUC24, ml$prediction),
      paired = list(n = sum(paired), map = paired_metrics(testing$TRUE_AUC24[paired], map[paired]),
                    ml = paired_metrics(testing$TRUE_AUC24[paired], ml$prediction[paired])),
      tuning = ml$selection, features = ml$features)
  }
  result <- list(modelId = scope$model_id[[1]], drug = scope$drug[[1]],
    administrationMode = unname(MODE_ID[[scope$mode[[1]]]]),
    baseModelSha256 = model_sha256(scope$model_id[[1]]),
    nSimulated = nrow(simulated), nIncluded = nrow(data), nExcluded = filtered$excluded,
    aucInclusionBounds = list(lower = filtered$bounds[[1]], upper = filtered$bounds[[2]]),
    nTraining = nrow(training), nHoldoutPool = nrow(holdout_pool), nHoldout = nrow(testing), seed = local_seed,
    mapUnavailableReason = unavailable,
    strategies = results)
  saveRDS(result, checkpoint)
  result
}

run_paired_benchmark <- function() {
  # Checkpoints are isolated by inputs, including helper and model sources.
  fingerprint <- paste(tools::md5sum(c(file.path(APP_ROOT, "ml", c("paired_auc_benchmark.R", "train_models_xgboost.R", "training-populations.json", "training-regimens.json")),
                                      list.files(file.path(APP_ROOT, "models"), "\\.cpp$", full.names = TRUE))), collapse = "")
  versions <- setNames(lapply(c("mrgsolve", "mapbayr", "xgboost"), function(package) as.character(utils::packageVersion(package))), c("mrgsolve", "mapbayr", "xgboost"))
  fingerprint <- paste(fingerprint, R.version.string, paste(unlist(versions), collapse = ":"))
  signature_file <- tempfile()
  writeLines(fingerprint, signature_file)
  signature <- unname(tools::md5sum(signature_file))
  unlink(signature_file)
  paired_directory <<- file.path(APP_ROOT, "ml", ".benchmark-cache", paste0(signature, "-", seed, "-", n_patients, if (smoke) "-smoke" else "-full"))
  dir.create(paired_directory, recursive = TRUE, showWarnings = FALSE)
  workers <- max(1L, as.integer(option_value("--workers", "1")))
  if (workers > 1L) {
    cluster <- parallel::makePSOCKcluster(min(workers, nrow(scopes)), outfile = file.path(paired_directory, "workers.log"))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterExport(cluster, c("APP_ROOT", "args", "paired_directory"), envir = environment())
    parallel::clusterEvalQ(cluster, {
      Sys.setenv(PMX_TRAIN_ROOT = APP_ROOT, PMX_TRAIN_FUNCTIONS_ONLY = "1")
      .training_args <- args
      source(file.path(APP_ROOT, "ml", "train_models_xgboost.R"))
      source(file.path(APP_ROOT, "ml", "paired_auc_benchmark.R"))
      NULL
    })
    cat("Worker log: ", file.path(paired_directory, "workers.log"), "\n", sep = "")
    results <- parallel::parLapplyLB(cluster, seq_len(nrow(scopes)), paired_scope, chunk.size = 1L)
  } else results <- lapply(seq_len(nrow(scopes)), paired_scope)
  document <- list(version = 1L, benchmarkDate = format(Sys.Date(), "%Y-%m-%d"), smoke = smoke,
    software = c(list(R = R.version.string), versions),
    methodology = list(target = "AUC24", comparator = "MAP-BE fitted to the same observed concentrations",
      trainingFraction = 0.75, tuningFolds = if (smoke) 3L else 10L,
      cohortSize = "9000 simulated profiles for Woillard when --n=9000 and 3000 per other model/mode; central 90% retained; paired evaluation uses 1098 Woillard or 300 other untouched holdout patients",
      learningTarget = "Direct AUC24, squared-error loss; no population AUC predictor or ratio target",
      c0 = "Pre-dose trough at steady state; sampled at tau-epsilon in the equivalent preceding cycle",
      c1 = "Tinf+1 h for intermittent infusion; 1 h for oral, bolus and continuous infusion",
      continuous = "Ongoing constant-rate infusion at steady state, sampled at t=0 and t=1 h",
      truth = "Trapezoidal integration of the noise-free individual curve, 0-24 h, step 0.1 h",
      inclusion = "Central 90% of simulated true AUC24 retained per model and administration mode (5th to 95th percentiles), before the 75/25 split",
      residual = "Published SIGMA structure retained, with each residual standard deviation multiplied by 0.10 for both sparse-observation simulation and MAP-BE",
      woillardRegimen = "Woillard ST is aligned with regimen: ST=1 for Prograf q12h and ST=0 for Advagraf q24h",
      tuning = "Ten-fold CV on training patients; two-sample grid max_depth={1,2,4}, eta={0.0261,0.05}, min_child_weight=40, mtry=min(7,p), up to 1000 rounds with 30-round early stopping",
      metric = "e=(estimate-truth)/truth; bias=100*mean(e); rRMSE=100*sqrt(mean(e^2)); within20=100*count(abs(e)<=0.2)/N",
      failures = "Bias/RMSE use finite positive estimates; failures count as outside +/-20%. Paired metrics restrict both methods to the same successful patients.",
      limitations = "Not a replication of Woillard 2021: C0 or C0+C1 instead of C0+C3, AUC24 instead of AUC0-12, generic central-AUC filtering instead of tacrolimus-specific rules, ten-fold residual-SD reduction instead of the article's exact residual values, broader regimens and different covariate sampling. The complete article search grid is not published; this grid includes its reported two-sample optimum. Benchmark predictors are not the currently deployed Shiny predictors."),
    results = results)
  destination <- if (nzchar(report_path)) report_path else file.path(APP_ROOT, "ml", "validation", "paired-auc-benchmark.json")
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(document, destination, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 9)
  cat("Paired benchmark complete: ", length(results), " scopes; ", destination, "\n", sep = "")
}
