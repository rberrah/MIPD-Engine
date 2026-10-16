args <- commandArgs(FALSE)
file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
root <- normalizePath(file.path(dirname(file), ".."), winslash = "/")
Sys.setenv(PMX_TRAIN_ROOT = root, PMX_TRAIN_FUNCTIONS_ONLY = "1")
.training_args <- c("--paired-benchmark", "--smoke")
source(file.path(root, "ml", "train_models_xgboost.R"))
source(file.path(root, "ml", "paired_auc_benchmark.R"))

metrics <- paired_metrics(c(100, 200, 100), c(120, 160, NA))
stopifnot(abs(metrics$relativeBiasPct) < 1e-10,
          abs(metrics$relativeRmsePct - 20) < 1e-10,
          abs(metrics$within20Pct - 200 / 3) < 1e-10, metrics$nFailed == 1)
stopifnot(is.null(paired_metrics(100, NA_real_)$relativeRmsePct))

filtered <- central_auc_filter(data.frame(TRUE_AUC24 = 1:100))
stopifnot(nrow(filtered$data) == 90L, filtered$excluded == 10L,
          min(filtered$data$TRUE_AUC24) == 6, max(filtered$data$TRUE_AUC24) == 95)

woillard_scope <- data.frame(model_id = "tacrolimus_woillard_ddi")
aligned <- align_regimen_covariates(woillard_scope, woillard_scope,
  data.frame(ST = c(0, 1)), data.frame(interval = c(12, 24)))
stopifnot(identical(aligned$ST, c(1, 0)))

regimens <- data.frame(amount = 100, interval = 24, infusion = 2)
stopifnot(sample_times(regimens, "IV_INTERMITTENT")$time[[2]] == 3)
regimens$infusion <- 0
stopifnot(sample_times(regimens, "ORAL")$time[[2]] == 1,
          sample_times(regimens, "IV_INTERMITTENT")$time[[2]] == 1,
          sample_times(regimens, "IV_INTERMITTENT")$time[[1]] < 24)
regimens$infusion <- 24
stopifnot(identical(sample_times(regimens, "IV_CONTINUOUS")$time, c(0, 1)))

columns <- c("patient_id", "TRUE_AUC24", "POP_AUC24", "PREV_CONC", "LAST_CONC", "LAST_TIME",
             "PREV_POP_CONC", "LAST_POP_CONC", "PREV_CONC_RATIO", "LAST_CONC_RATIO", "CONC_DIFF", "TIME_DIFF", "WT", "DOSE")
fake <- as.data.frame(setNames(rep(list(1:5), length(columns)), columns))
stopifnot(identical(paired_features(fake, FALSE), c("PREV_CONC", "WT", "DOSE")),
          !any(c("TRUE_AUC24", "patient_id", "POP_AUC24") %in% paired_features(fake, TRUE)))
cat("Paired benchmark metrics, sampling and no-leakage checks passed.\n")
