suppressPackageStartupMessages({
  library(mrgsolve)
  library(jsonlite)
})

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file_argument), ".."), winslash = "/", mustWork = TRUE)

source(file.path(APP_ROOT, "R", "model_library.R"), local = TRUE)

failures <- list()

required_metadata <- c(
  "model", "citation", "doi", "population", "modelType", "sourceStatus",
  "administrationModes", "administrationCategories", "administrationCategoriesEn",
  "implementationStatus", "implementationStatusLabel", "implementationStatusLabelEn",
  "analysisEligible", "analysisEligibilityLabel", "analysisEligibilityLabelEn",
  "verificationStatus", "validationStatus", "contextOfUse"
)
missing_metadata <- setdiff(required_metadata, names(MODEL_CATALOG))
if (length(missing_metadata)) stop("Missing catalog metadata: ", paste(missing_metadata, collapse = ", "))
if (any(c("source", "provenance") %in% names(MODEL_CATALOG))) {
  stop("Application provenance and filename-derived source fields are forbidden in the model catalog.")
}
if (any(!nzchar(MODEL_CATALOG$citation)) || any(!nzchar(MODEL_CATALOG$doi))) {
  stop("Every library model must cite an article and a DOI.")
}
if (any(grepl("\\bDDI\\b", MODEL_CATALOG$model, ignore.case = TRUE))) {
  stop("Technical DDI suffixes are forbidden in displayed model names.")
}
woillard <- MODEL_CATALOG[MODEL_CATALOG$id == "tacrolimus_woillard_ddi", , drop = FALSE]
if (nrow(woillard) != 1L || !identical(woillard$model[[1]], "Woillard")) {
  stop("The tacrolimus Woillard model must be displayed as `Woillard`.")
}

analysis_ids <- unname(catalog_choices(analysis_only = TRUE))
expected_analysis_ids <- MODEL_CATALOG$id[!is.na(MODEL_CATALOG$analysisEligible) & MODEL_CATALOG$analysisEligible]
if (!setequal(analysis_ids, expected_analysis_ids)) {
  stop("Analysis choices do not match the models explicitly eligible for MAP estimation.")
}
if (any(c("amik_burdet", "cefazoline_schmitz", "amox_mellon", "linez_buerger", "voriconazole_vandenborn_ddi") %in% analysis_ids)) {
  stop("A documented article adaptation is incorrectly exposed as a faithful MAP prior.")
}
if (!"tacrolimus_woillard_ddi" %in% analysis_ids) {
  stop("The audited baseline Woillard prior must remain available for MAP estimation.")
}

for (model_id in MODEL_CATALOG$id) {
  result <- tryCatch({
    record <- model_record(model_id)
    routes <- model_routes(record)
    if (!length(routes) || any(!routes %in% c("IV", "Oral")) || anyDuplicated(routes)) {
      stop("Invalid administration routes: ", paste(routes, collapse = ", "))
    }
    modes <- model_administration_modes(record)
    if (!length(modes) || any(!modes %in% c("ORAL", "IV_INTERMITTENT", "IV_CONTINUOUS")) || anyDuplicated(modes)) {
      stop("Invalid administration modes: ", paste(modes, collapse = ", "))
    }
    model <- compile_model(model_id = model_id)
    contract <- validate_model_contract(model)
    if (!isTRUE(contract$ok)) stop(paste(contract$errors, collapse = " | "))
    for (route in routes) {
      if (!length(model_administration_modes(record, route))) stop("No administration mode declared for route ", route, ".")
      compartment <- model_administration_cmt(record, route)
      if (!compartment %in% model@cmtL) {
        stop(route, " administration compartment `", compartment, "` is absent from the compiled model.")
      }
    }
    contract
  }, error = function(error) error)

  if (inherits(result, "error")) {
    failures[[model_id]] <- conditionMessage(result)
    cat("FAIL ", model_id, ": ", conditionMessage(result), "\n", sep = "")
  } else {
    cat("OK   ", model_id, "\n", sep = "")
  }
}

vancomycin_routes <- lapply(c("vanco_goti", "vanco_pkjust", "vanco_roberts"), function(id) {
  model_routes(model_record(id))
})
if (!all(vapply(vancomycin_routes, identical, logical(1), "IV"))) {
  stop("All vancomycin models must be IV-only.")
}
if (!identical(model_routes(model_record("tacrolimus_woillard_ddi")), "Oral")) {
  stop("The Woillard tacrolimus model must be oral-only.")
}
if (!identical(model_administration_mode(model_record("vanco_roberts"), "IV"), "IV_CONTINUOUS") ||
    !identical(model_administration_mode(model_record("amox_rambaud"), "IV"), "IV_CONTINUOUS")) {
  stop("Roberts and Rambaud must be categorized as continuous IV models.")
}
if (!setequal(model_administration_modes(model_record("vanco_pkjust"), "IV"), c("IV_INTERMITTENT", "IV_CONTINUOUS")) ||
    !identical(model_administration_mode(model_record("vanco_goti"), "IV"), "IV_INTERMITTENT")) {
  stop("Revilla must support both IV modes and Goti must remain intermittent IV.")
}

expected_covariates <- list(
  amox_carlier = "CRCL",
  cefazoline_schmitz = c("WT", "AGE", "CLCR"),
  cefepime_jonckheere = c("CRCL", "IHD"),
  dapto_garreau = c("CLCR", "AGE", "WT", "SEX", "RIF"),
  levo_canoui = c("CLCR", "AGE"),
  dapto_dvorchik_IV_adults = c("WT", "SEX", "DIAL", "INFECT"),
  cefepime_an = c("AGE", "WT", "HT", "CREAT", "SEX", "CRRT"),
  vanco_pkjust = c("AGE", "CRCL", "CREAT", "WT"),
  rifampicine_marsot = "FUS"
)
for (model_id in names(expected_covariates)) {
  parsed <- parse_covariates(read_library_code(model_id))$name
  missing <- setdiff(expected_covariates[[model_id]], parsed)
  if (length(missing)) stop(model_id, " is missing parsed covariates: ", paste(missing, collapse = ", "))
}

reference_parameters <- list(
  amox_carlier = c(TVCL = 10.0, REF_CRCL = 102.0),
  cefazoline_lanoiselee = c(TVCL = 2.86, THETA_CLCR = 0.79, TVVC = 5.20, TVVP = 4.56, TVQ = 10.9, REF_CLCR = 80.0),
  cefepime_jonckheere = c(TVV1 = 18.3, TVV2 = 11.1, TVQ = 6.63, THETA1 = 2.88, THETA2 = 0.368, CLDIAL = 5.74, CLOTHER = 0.87),
  cefepime_seo = c(TVCL = 6.60, TVV1 = 13.3, TVQ = 16.5, TVV2 = 13.0, CL_EXP = 0.656, REF_CLCR = 77.21),
  dapto_garreau = c(TVCL = 0.365, TVV1 = 3.59, TVQ = 0.752, TVV2 = 4.71, RIF_V1 = -0.121),
  levo_canoui = c(TVCL = 5.57, TVV = 96.3, KA = 1.6, CLCR_CL = 0.684, AGE_CL = -0.312),
  rifampicine_jing = c(TVCL = 4.02, TVVD = 57.8, TVKA = 1.61),
  rifampicine_marsot = c(TVKa = 1.15, TVCL_NOFUS = 13.7, TVV_NOFUS = 61.1, TVCL_FUS = 5.1, TVV_FUS = 23.8),
  tacrolimus_woillard_ddi = c(TVCL_TAC = 21.2, TVV1_TAC = 486, TVQ_TAC = 79, TVV2_TAC = 271, TVKTR_TAC = 3.34),
  cefepime_an = c(TVCLR = 2, TVCLCRRT = 1.64, TVCLNR = 0.526, TVVC = 13.4, TVVP = 7.52, TVQ = 12),
  dapto_dvorchik_IV_adults = c(TVCL = 0.807, TVCL_DIAL = 0.269, FEMALE_CL = 0.801, INFECT_V2 = 1.93),
  vanco_pkjust = c(TVCLCR = 0.67, TVV1 = 0.82, AGE_CL = -0.24, CREAT_V = 2.49),
  vanco_goti = c(TVCL = 4.5, TVV1 = 58.4, TVV2 = 38.4, TVQ = 6.5, DIAL_CL = 0.7, DIAL_V1 = 0.5),
  vanco_roberts = c(TVCL = 4.58, TVV = 1.53)
)
for (model_id in names(reference_parameters)) {
  model <- compile_model(model_id = model_id)
  actual <- unlist(param(model))[names(reference_parameters[[model_id]])]
  if (!isTRUE(all.equal(actual, reference_parameters[[model_id]], tolerance = 1e-12, check.attributes = FALSE))) {
    stop(model_id, " reference parameters differ from the audited publication values.")
  }
}

source_checks <- list(
  cefazoline_schmitz = c("CLCR     : 90.0", "ADD  : 0.714025", "pow(clcr_safe / 90.0, POW_RCC)"),
  cefepime_an = c("ADD  : 41.4736", "TVCLCRRT : 1.64"),
  dapto_dvorchik_IV_adults = c("4.28 // additive variance", "TVCL_DIAL : 0.269"),
  rifampicine_marsot = c("TVKa       : 1.15", "2.256 //", "dxdt_DEPOT =  Ktr * TR"),
  tacrolimus_woillard_ddi = c("double CL_IND =", "double CL_TAC_DDI = CL_IND", "+ CL_OTH"),
  vanco_pkjust = c("PROP : 0", "ADD  : 17.8929", "CRCL  : 74.7")
)
for (model_id in names(source_checks)) {
  code <- read_library_code(model_id)
  if (any(!vapply(source_checks[[model_id]], grepl, logical(1), x = code, fixed = TRUE))) {
    stop(model_id, " no longer contains the audited article constants.")
  }
}

if (length(failures)) {
  stop(length(failures), " model(s) failed library validation.")
}

cat("TDM library validation OK: ", nrow(MODEL_CATALOG), "/", nrow(MODEL_CATALOG), " models.\n", sep = "")
