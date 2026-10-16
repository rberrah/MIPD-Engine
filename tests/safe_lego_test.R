suppressPackageStartupMessages({
  library(mrgsolve)
  library(mapbayr)
  library(dplyr)
  library(jsonlite)
})

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file_argument), ".."), winslash = "/", mustWork = TRUE)

source(file.path(APP_ROOT, "R", "model_library.R"), local = TRUE)

expect_error <- function(expression, pattern) {
  message <- tryCatch({
    force(expression)
    NULL
  }, error = function(error) conditionMessage(error))
  if (is.null(message) || !grepl(pattern, message, fixed = TRUE)) {
    stop("Expected error containing `", pattern, "`, got: ", message %||% "no error")
  }
}

`%||%` <- function(value, fallback) if (is.null(value) || !length(value)) fallback else value
source(file.path(APP_ROOT, "R", "engine.R"), local = TRUE)

stopifnot(identical(parse_dosing_intervals(c("8", "12"), "168, 336; 12"), c(8, 12, 168, 336)))
expect_error(parse_dosing_intervals(character(), "0, 168"), "positive numbers")

oral_one_compartment <- list(
  version = 1,
  nodes = list(
    list(id = 1, kind = "depot", name = "depot", dose = 100),
    list(id = 2, kind = "central", name = "centr", dose = 0, vol = 30)
  ),
  edges = list(
    list(from = 1, to = 2, k = 1),
    list(from = 2, to = "OUT", k = 0.2)
  ),
  covariates = list(
    list(name = "WT", type = "continuous", target = "v_centr", reference = 70, comparison = 90, beta = 0.75),
    list(name = "SEX", type = "categorical", target = "v_centr", reference = 0, comparison = 1, beta = 0.2)
  )
)

safe_code <- lego_model_code(oral_one_compartment)
stopifnot(startsWith(safe_code, LEGO_SPEC_PREFIX))
stopifnot(grepl("$PARAM @covariates", safe_code, fixed = TRUE))
stopifnot(grepl("pow(WT/70", safe_code, fixed = TRUE))
stopifnot(grepl("exp(BETA_SEX_2 * (SEX == 1))", safe_code, fixed = TRUE))
stopifnot(!grepl("LEGO_INPUT", safe_code, fixed = TRUE), !grepl("$PLUGIN evtools", safe_code, fixed = TRUE))

population_specification <- oral_one_compartment
population_specification$population <- list(
  iivVariances = list(k_centr_e = 0.16, v_centr = 0.25),
  iivDistributions = list(k_centr_e = "lognormal", v_centr = "normal"),
  iivCovariances = list(`k_centr_e::v_centr` = 0.04),
  residualError = list(type = "additive", additive = 0.3, proportional = 0.2)
)
population_code <- lego_model_code(population_specification)
stopifnot(
  grepl("$OMEGA @block\n0.16\n0.04 0.25", population_code, fixed = TRUE),
  grepl("double v_L2_centr = TV_v_L2_centr", population_code, fixed = TRUE),
  grepl("+ ETA2 + ETA(2);", population_code, fixed = TRUE),
  grepl("PROP : 0 : proportional residual variance", population_code, fixed = TRUE),
  grepl("ADD : 0.09", population_code, fixed = TRUE),
  grepl("double DV = IPRED + EPS(2);", population_code, fixed = TRUE)
)
population_dir <- tempfile("safe-lego-population-")
dir.create(population_dir, recursive = TRUE)
on.exit(unlink(population_dir, recursive = TRUE, force = TRUE), add = TRUE)
population_model <- compile_model(
  custom_code = population_code,
  allow_custom = FALSE,
  custom_soloc = population_dir,
  custom_cache = new.env(parent = emptyenv())
)
population_contract <- validate_model_contract(population_model)
stopifnot(isTRUE(population_contract$ok), population_contract$n_eta == 2L, population_contract$n_sigma == 2L)
mapbayr::check_mapbayr_model(population_model)

proportional_specification <- population_specification
proportional_specification$population$residualError$type <- "proportional"
proportional_code <- lego_model_code(proportional_specification)
stopifnot(
  grepl("PROP : 0.04", proportional_code, fixed = TRUE),
  grepl("ADD : 0 : additive residual variance", proportional_code, fixed = TRUE),
  grepl("double DV = IPRED * (1 + EPS(1));", proportional_code, fixed = TRUE)
)
population_roundtrip <- lego_spec_from_code(population_code)
stopifnot(
  identical(names(population_roundtrip$population$iivVariances), c("k_centr_e", "v_centr")),
  identical(unname(population_roundtrip$population$iivDistributions), c("lognormal", "normal"))
)
invalid_population <- population_specification
invalid_population$population$iivCovariances[[1]] <- 1
expect_error(lego_model_code(invalid_population), "positive definite")

logit_population <- population_specification
logit_population$nodes[[2]]$vol <- 0.5
logit_population$population$iivDistributions$v_centr <- "logit"
logit_code <- lego_model_code(logit_population)
stopifnot(grepl("1.0/(1.0 + exp(-(log(TV_v_L2_centr/(1.0-TV_v_L2_centr))", logit_code, fixed = TRUE))
invalid_logit <- population_specification
invalid_logit$population$iivDistributions$v_centr <- "logit"
expect_error(lego_model_code(invalid_logit), "strictly between 0 and 1")

legacy_specification <- oral_one_compartment
legacy_specification$covariates <- list(
  list(name = "WT", target = "v_centr", reference = 70, beta = 0.75)
)
legacy_code <- lego_model_code(legacy_specification)
stopifnot(grepl("pow(WT/70", legacy_code, fixed = TRUE))

many_covariates <- oral_one_compartment
many_covariates$covariates <- lapply(seq_len(12), function(index) list(
  name = paste0("COV", index), type = "continuous", target = "v_centr",
  reference = 1, comparison = 1.25, beta = 0.1
))
many_covariate_code <- lego_model_code(many_covariates)
stopifnot(nrow(parse_covariates(many_covariate_code)) == 12L)

q_specification <- list(
  version = 3,
  nodes = list(
    list(id = 1, kind = "central", name = "central", dose = 100, vol = 20),
    list(id = 2, kind = "periph", name = "periph", dose = 0, vol = 40)
  ),
  edges = list(
    list(from = 1, to = 2, transferParameterization = "clearance", q = 3.2),
    list(from = 2, to = 1, transferParameterization = "clearance", q = 3.2),
    list(from = 1, to = "OUT", eliminationParameterization = "clearance", cl = 4)
  ),
  covariates = list(
    list(name = "WT", type = "continuous", target = "q_central_periph", reference = 70, comparison = 90, beta = 0.5),
    list(name = "WT", type = "continuous", target = "v_central", reference = 70, comparison = 90, beta = 1)
  )
)
q_code <- lego_model_code(q_specification)
stopifnot(
  grepl("TV_q_L1_central_L2_periph : 3.2", q_code, fixed = TRUE),
  grepl("q_L1_central_L2_periph*L1_central/v_L1_central", q_code, fixed = TRUE),
  grepl("q_L1_central_L2_periph*L2_periph/v_L2_periph", q_code, fixed = TRUE),
  grepl("BETA_WT_1", q_code, fixed = TRUE),
  grepl("BETA_WT_2", q_code, fixed = TRUE)
)

renal_specification <- q_specification
renal_specification$covariates <- list(
  list(name = "WT", type = "continuous", form = "power", target = "cl_central", reference = 74.4, comparison = 90, beta = 0.727),
  list(name = "CRCL", type = "continuous", form = "linear-additive", target = "cl_central", reference = 0, comparison = 100, maximum = 150, beta = 0.024)
)
renal_code <- lego_model_code(renal_specification)
stopifnot(
  grepl("TV_cl_L1_central * pow(WT/74.4, BETA_WT_1) + BETA_CRCL_2 * (fmin(CRCL, 150) - 0)", renal_code, fixed = TRUE),
  identical(lego_spec_from_code(renal_code)$covariates[[2]]$form, "linear-additive"),
  identical(lego_spec_from_code(renal_code)$covariates[[2]]$maximum, 150)
)

session_dir <- tempfile("safe-lego-test-")
dir.create(session_dir, recursive = TRUE)
on.exit(unlink(session_dir, recursive = TRUE, force = TRUE), add = TRUE)

model <- compile_model(
  custom_code = safe_code,
  allow_custom = FALSE,
  custom_soloc = session_dir,
  custom_cache = new.env(parent = emptyenv())
)
contract <- validate_model_contract(model)
if (!isTRUE(contract$ok)) stop(paste(contract$errors, collapse = " | "))
stopifnot(all(c("WT", "SEX") %in% model_param_names(model)))
covariate_definition <- parse_covariates(safe_code)
stopifnot(nrow(covariate_definition) == 2L, all(c("WT", "SEX") %in% covariate_definition$name))

q_model <- compile_model(
  custom_code = q_code,
  allow_custom = FALSE,
  custom_soloc = session_dir,
  custom_cache = new.env(parent = emptyenv())
)
q_contract <- validate_model_contract(q_model)
if (!isTRUE(q_contract$ok)) stop(paste(q_contract$errors, collapse = " | "))
stopifnot("WT" %in% model_param_names(q_model))

renal_model <- compile_model(
  custom_code = renal_code,
  allow_custom = FALSE,
  custom_soloc = session_dir,
  custom_cache = new.env(parent = emptyenv())
)
renal_contract <- validate_model_contract(renal_model)
if (!isTRUE(renal_contract$ok)) stop(paste(renal_contract$errors, collapse = " | "))
stopifnot(all(c("WT", "CRCL") %in% model_param_names(renal_model)))

advanced_absorption <- list(
  version = 2,
  nodes = list(
    list(id = 1, kind = "depot", name = "rapid", dose = 100, inputType = "zero_order", inputDuration = 1.5, tlag = 0, doseFraction = 30),
    list(id = 2, kind = "depot", name = "slow", dose = 100, inputType = "bolus", inputDuration = 1, tlag = 1.5, doseFraction = 70),
    list(id = 3, kind = "central", name = "centr", dose = 0, vol = 30, inputType = "bolus", inputDuration = 1, tlag = 0, doseFraction = 100)
  ),
  edges = list(
    list(from = 1, to = 3, kinetics = "michaelis_menten", k = 1, vmax = 40, km = 20, eliminationParameterization = "rate", cl = 5),
    list(from = 2, to = 3, kinetics = "first_order", k = 0.35, vmax = 10, km = 10, eliminationParameterization = "rate", cl = 5),
    list(from = 3, to = "OUT", kinetics = "first_order", k = 0.17, vmax = 10, km = 10, eliminationParameterization = "clearance", cl = 5)
  ),
  covariates = list()
)
advanced_code <- lego_model_code(advanced_absorption)
stopifnot(
  grepl("$PLUGIN evtools", advanced_code, fixed = TRUE),
  grepl("split dose input [ADM]", advanced_code, fixed = TRUE),
  grepl("evt::infuse(AMT*f_L1_rapid", advanced_code, fixed = TRUE),
  grepl("evt::retime(route_2, TIME + tlag_L2_slow)", advanced_code, fixed = TRUE),
  grepl("vmax_L1_rapid_L3_centr", advanced_code, fixed = TRUE),
  grepl("cl_L3_centr*L3_centr/v_L3_centr", advanced_code, fixed = TRUE)
)
advanced_model <- compile_model(
  custom_code = advanced_code,
  allow_custom = FALSE,
  custom_soloc = session_dir,
  custom_cache = new.env(parent = emptyenv())
)
advanced_contract <- validate_model_contract(advanced_model)
if (!isTRUE(advanced_contract$ok)) stop(paste(advanced_contract$errors, collapse = " | "))
advanced_map_data <- build_map_data(
  doses = data.frame(time = 0, amount = 100, interval = 12, count = 1, infusion = 0, ss = 1),
  observations = data.frame(time = 6, concentration = 2),
  adm_cmt = advanced_contract$adm_cmt,
  obs_cmt = advanced_contract$obs_cmt,
  covariates = list(),
  split_lego = TRUE
)
stopifnot(min(advanced_map_data$time) == -600, all(advanced_map_data$ss == 0))
advanced_estimate <- mapbayr::mapbayest(
  advanced_model,
  data = advanced_map_data,
  hessian = FALSE,
  verbose = FALSE,
  progress = FALSE
)
stopifnot(inherits(advanced_estimate, "mapbayests"))
advanced_output <- as.data.frame(mrgsolve::mrgsim(
  mrgsolve::zero_re(advanced_model),
  events = mrgsolve::ev(amt = 100, cmt = advanced_contract$adm_cmt),
  end = 12,
  delta = 0.1
))
stopifnot(any(advanced_output$DV > 0), all(is.finite(advanced_output$DV)))
advanced_ss_output <- as.data.frame(mrgsolve::mrgsim(
  mrgsolve::zero_re(advanced_model),
  events = mrgsolve::ev(amt = 100, cmt = advanced_contract$adm_cmt, ii = 12, addl = 50),
  start = 600,
  end = 612,
  delta = 0.1,
  recsort = 3
))
advanced_ss_output$time <- advanced_ss_output$time - 600
concentration_column <- grep("^CONC_.*centr$", names(advanced_ss_output), value = TRUE)[[1]]
stopifnot(max(advanced_ss_output[[concentration_column]][advanced_ss_output$time == 0]) > 0, all(is.finite(advanced_ss_output[[concentration_column]])))
advanced_long_output <- as.data.frame(mrgsolve::mrgsim(
  mrgsolve::zero_re(advanced_model),
  events = mrgsolve::ev(amt = 100, cmt = advanced_contract$adm_cmt, ii = 12, addl = 99),
  end = 1200,
  delta = 0.1,
  recsort = 3
))
steady_profile <- advanced_ss_output[advanced_ss_output$time > 0, , drop = FALSE]
long_profile <- advanced_long_output[advanced_long_output$time > 1188, , drop = FALSE]
long_profile$time <- long_profile$time - 1188
for (check_time in c(3, 6, 9, 12)) {
  steady_value <- stats::approx(steady_profile$time, steady_profile[[concentration_column]], xout = check_time, ties = "ordered")$y
  long_value <- stats::approx(long_profile$time, long_profile[[concentration_column]], xout = check_time, ties = "ordered")$y
  stopifnot(abs(steady_value / long_value - 1) < 0.1)
}
advanced_fit <- list(
  id = "advanced_lego",
  model = advanced_model,
  estimate = NULL,
  contract = advanced_contract,
  current_covariates = list(),
  split_lego = TRUE
)
engine_profile <- simulate_regimen(advanced_fit, dose = 100, interval = 12, infusion = 0, horizon = 12, delta = 0.1)
stopifnot(abs(min(engine_profile$time)) < 1e-8, abs(max(engine_profile$time) - 12) < 1e-8, max(engine_profile$concentration) > 0)

koka <- list(
  version = 3,
  nodes = list(
    list(id = 1, kind = "central", name = "central", dose = 100, vol = 391, inputType = "zero_order", inputDuration = 319, inputDurationTlagOf = 2, tlag = 0, doseFraction = 16.8),
    list(id = 2, kind = "depot", name = "slow", dose = 100, inputType = "bolus", inputDuration = 1, tlag = 319, doseFraction = 83.2, fractionComplementOf = 1)
  ),
  edges = list(
    list(from = 2, to = 1, kinetics = "first_order", k = 0.000488),
    list(from = 1, to = "OUT", kinetics = "first_order", k = 0.005, eliminationParameterization = "clearance", cl = 4.95)
  ),
  covariates = list(
    list(name = "IVOL", type = "continuous", scope = "administration", target = "f_central", reference = 1.75, comparison = 1, beta = 0.2),
    list(name = "IVOL", type = "continuous", scope = "administration", target = "cl_central", reference = 1.75, comparison = 1, beta = 0.1)
  )
)
koka_code <- lego_model_code(koka)
stopifnot(
  grepl("evt::infuse(AMT*f_L1_central, 2, AMT*f_L1_central/tlag_L2_slow)", koka_code, fixed = TRUE),
  grepl("evt::bolus(AMT*(1-f_L1_central), 3)", koka_code, fixed = TRUE),
  grepl("[administration]", koka_code, fixed = TRUE),
  !grepl("TV_f_L2_slow", koka_code, fixed = TRUE),
  !grepl("TV_tk0_L1_central", koka_code, fixed = TRUE)
)
koka_model <- compile_model(custom_code = koka_code, allow_custom = FALSE, custom_soloc = session_dir, custom_cache = new.env(parent = emptyenv()))
koka_contract <- validate_model_contract(koka_model)
stopifnot(isTRUE(koka_contract$ok), identical(parse_covariates(koka_code)$scope, "administration"))
koka_data <- build_map_data(
  doses = data.frame(time = c(0, 24), amount = 100, interval = 24, count = 1, infusion = 0, ss = 0, IVOL = c(1.75, 2.5)),
  observations = data.frame(time = c(12, 36), concentration = c(0.01, 0.02)),
  adm_cmt = koka_contract$adm_cmt,
  obs_cmt = koka_contract$obs_cmt,
  covariates = list(IVOL = 1.75),
  split_lego = TRUE
)
stopifnot(koka_data$IVOL[koka_data$time == 12] == 1.75, koka_data$IVOL[koka_data$time == 36] == 2.5)

pp6m <- list(
  version = 3,
  nodes = list(
    list(id = 1, kind = "depot", name = "slow", dose = 1000, doseFraction = 79.1),
    list(id = 2, kind = "depot", name = "rapid", dose = 1000, doseFraction = 20.9, fractionComplementOf = 1),
    list(id = 3, kind = "central", name = "central", dose = 0, vol = 1960)
  ),
  edges = list(
    list(from = 1, to = 3, kinetics = "hill", k = 0.1, vmax = 0.0904, km = 120, gamma = 1.44),
    list(from = 2, to = 3, kinetics = "hill", k = 0.1, vmax = 0.149, km = 23.8, gamma = 1),
    list(from = 3, to = "OUT", kinetics = "first_order", k = 0.002, eliminationParameterization = "clearance", cl = 3.9)
  ),
  covariates = list()
)
pp6m_code <- lego_model_code(pp6m)
stopifnot(grepl("pow(L1_slow, gamma_L1_slow_L3_central)", pp6m_code, fixed = TRUE))
pp6m_model <- compile_model(custom_code = pp6m_code, allow_custom = FALSE, custom_soloc = session_dir, custom_cache = new.env(parent = emptyenv()))
stopifnot(isTRUE(validate_model_contract(pp6m_model)$ok))

invalid_fractions <- advanced_absorption
invalid_fractions$nodes[[2]]$doseFraction <- 60
expect_error(lego_model_code(invalid_fractions), "must total 100%")

expect_error(
  compile_model(
    custom_code = "$GLOBAL\n#include <cstdlib>\n$MAIN\nsystem(\"whoami\");",
    allow_custom = FALSE,
    custom_soloc = session_dir
  ),
  "specification marker is missing"
)

injected <- oral_one_compartment
injected$nodes[[2]]$name <- "centr;system_call"
expect_error(lego_model_code(injected), "Invalid Lego text field")

unknown_source <- oral_one_compartment
unknown_source$nodes[[3]] <- list(
  id = 3, kind = "effect", name = "Ce", dose = 0, ke0 = 0.4, source = 999
)
expect_error(lego_model_code(unknown_source), "does not exist")

invalid_covariate <- oral_one_compartment
invalid_covariate$covariates[[1]]$name <- "DV"
expect_error(lego_model_code(invalid_covariate), "Reserved Lego covariate name")

unknown_target <- oral_one_compartment
unknown_target$covariates[[1]]$target <- "CL_DOES_NOT_EXIST"
expect_error(lego_model_code(unknown_target), "Unknown Lego covariate target")

invalid_type <- oral_one_compartment
invalid_type$covariates[[1]]$type <- "free_cpp"
expect_error(lego_model_code(invalid_type), "Invalid Lego text field")

same_categories <- oral_one_compartment
same_categories$covariates[[2]]$comparison <- 0
expect_error(lego_model_code(same_categories), "reference and comparison must differ")

cat("Safe Lego compilation OK; arbitrary C++ and invalid specifications rejected.\n")
