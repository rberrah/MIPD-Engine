suppressPackageStartupMessages(library(mrgsolve))

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file_argument), ".."), winslash = "/", mustWork = TRUE)

`%||%` <- function(value, fallback) if (is.null(value) || !length(value)) fallback else value
source(file.path(APP_ROOT, "R", "model_library.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "engine.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "ddi_engine.R"), local = TRUE)

affected <- ddi_library_context("tacrolimus_woillard_ddi", "Oral")
driver <- ddi_library_context("voriconazole_vandenborn_ddi", "Oral")
code_parameters <- ddi_code_parameter_table(read_library_code("tacrolimus_woillard_ddi"))
result <- ddi_simulate(
  list(
    affected = list(dose = 3, interval = 12, infusion = 0),
    driver = list(dose = 200, interval = 12, infusion = 0),
    start_day = 3,
    stop_day = 10,
    followup_days = 3,
    target = "TVCL_TAC",
    type = "factor",
    factor = 0.5,
    strength = 1,
    c50 = 1
  ),
  affected,
  driver,
  delta = 0.25
)
multi_config <- result$config
multi_config$targets <- list(
  list(parameter = "TVCL_TAC", fraction = 0.5),
  list(parameter = "TVV1_TAC", fraction = 0.25)
)
multi <- ddi_simulate(multi_config, affected, driver, delta = 0.25)
active_effect <- subset(multi$effect, day >= multi$config$start_day & day < multi$config$stop_day)

personalized <- ddi_fit_context(list(
  id = affected$id,
  label = affected$label,
  model = affected$model,
  estimate = NULL,
  contract = list(route = affected$route, adm_cmt = affected$adm_cmt),
  split_lego = FALSE,
  current_covariates = list(HT = 40)
))
lego_context <- ddi_fit_context(list(
  id = "custom",
  label = "Imported Lego model",
  model = affected$model,
  estimate = NULL,
  contract = list(route = "Oral", adm_cmt = affected$adm_cmt),
  split_lego = FALSE,
  current_covariates = list(HT = 40)
))
inhibition <- ddi_modifier("inhibition", c(0, 1, 10), c(0, 1, 2), 0, 3, 1, 0.8, 1)
induction <- ddi_modifier("induction", c(0, 1, 10), c(0, 1, 2), 0, 3, 1, 2, 1)

stopifnot(
  nrow(result$affected_profile) > 100,
  nrow(result$driver_profile) > 50,
  nrow(result$schedule) == 2,
  "TVCL_TAC" %in% code_parameters$name,
  !"HT" %in% code_parameters$name,
  is.finite(result$metrics$baseline$auc),
  result$metrics$auc_ratio > 1,
  result$metrics$cmin_ratio > 1,
  abs(result$metrics$minimum_modifier - 0.5) < 1e-8,
  abs(result$metrics$maximum_modifier - 1) < 1e-8,
  identical(unique(multi$effect$target), c("TVCL_TAC", "TVV1_TAC")),
  abs(min(subset(active_effect, target == "TVCL_TAC")$modifier) - 0.75) < 1e-8,
  abs(min(subset(active_effect, target == "TVV1_TAC")$modifier) - 0.875) < 1e-8,
  length(multi$metrics$baseline_parameters) == 2,
  inherits(try(ddi_validate_targets(list(targets = list(list(parameter = "CL", fraction = 1.1)))), silent = TRUE), "try-error"),
  inherits(try(ddi_validate_targets(list(targets = list(list(parameter = "CL", fraction = 1), list(parameter = "CL", fraction = 0.5)))), silent = TRUE), "try-error"),
  personalized$source == "tdm",
  lego_context$drug == "Lego",
  lego_context$source == "tdm",
  !"HT" %in% ddi_parameter_table(personalized)$name,
  inhibition[[1]] == 1,
  tail(inhibition, 1) < 0.3,
  induction[[1]] == 1,
  tail(induction, 1) > 2
)

cat("Generic DDI engine test passed.\n")
