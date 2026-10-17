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

source_review <- jsonlite::fromJSON(
  file.path(APP_ROOT, "models", "source-review.json"),
  simplifyVector = FALSE
)
reviewed_ids <- vapply(source_review$reviewed, `[[`, character(1), "id")
if (anyDuplicated(reviewed_ids) || !setequal(reviewed_ids, MODEL_CATALOG$id)) {
  stop("Source-review coverage must contain every active catalog model exactly once.")
}
allowed_review_results <- c(
  "ARTICLE_ALIGNED",
  "DOCUMENTED_ADAPTATION_CONFIRMED",
  "SOURCE_AMBIGUITY_DOCUMENTED"
)
review_results <- vapply(source_review$reviewed, `[[`, character(1), "result")
if (any(!review_results %in% allowed_review_results)) {
  stop("Unknown source-review result.")
}
review_status <- setNames(review_results, reviewed_ids)
expected_implementation <- ifelse(review_status[MODEL_CATALOG$id] == "ARTICLE_ALIGNED", "ARTICLE", "ADAPTED")
# Implementation lineage (historical ARTICLE code) and source fidelity are
# separate: a newly identified source ambiguity is a conservative downgrade
# and must not require rewriting frozen model metadata just to warn users.
conservative_source_downgrade <- MODEL_CATALOG$implementationStatus == "ARTICLE" &
  review_status[MODEL_CATALOG$id] == "SOURCE_AMBIGUITY_DOCUMENTED"
if (any(MODEL_CATALOG$implementationStatus != expected_implementation &
        !conservative_source_downgrade)) {
  stop("Catalog implementation status must agree with the source-review conclusion.")
}
if (any(!nzchar(vapply(source_review$reviewed, `[[`, character(1), "note")))) {
  stop("Every reviewed model must have a source-review note.")
}
allowed_evidence <- names(source_review$evidenceDefinitions)
review_evidence <- vapply(source_review$reviewed, `[[`, character(1), "evidence")
if (!length(allowed_evidence) || any(!review_evidence %in% allowed_evidence)) {
  stop("Every source review must use a defined evidence level.")
}
if (any(!vapply(source_review$reviewed, function(review) {
  is.logical(review$fullTextReviewed) && length(review$fullTextReviewed) == 1L
}, logical(1)))) {
  stop("Every source review must state whether the primary full text was reviewed.")
}

pdf_inventory <- jsonlite::fromJSON(
  file.path(APP_ROOT, "models", "pdf-inventory.json"),
  simplifyVector = FALSE
)
pdf_category_counts <- unlist(pdf_inventory$counts[setdiff(names(pdf_inventory$counts), "documents")], use.names = FALSE)
if (length(pdf_inventory$documents) != 47L ||
    pdf_inventory$counts$documents != 47L ||
    sum(pdf_category_counts) != 47L) {
  stop("The private-source inventory must account for exactly 47 PDFs.")
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
  amiodarone_pollak_ddi = c(TVCL = 9.541667, TVV1 = 882, TVQ = 24.958333, TVV2 = 12700, TVKA = 0.6),
  cefazoline_lanoiselee = c(TVCL = 2.86, THETA_CLCR = 0.79, TVVC = 5.20, TVVP = 4.56, TVQ = 10.9, REF_CLCR = 80.0),
  cefepime_jonckheere = c(TVV1 = 18.3, TVV2 = 11.1, TVQ = 6.63, THETA1 = 2.88, THETA2 = 0.368, CLDIAL = 5.74, CLOTHER = 0.87),
  cefepime_seo = c(TVCL = 6.60, TVV1 = 13.3, TVQ = 16.5, TVV2 = 13.0, CL_EXP = 0.656, REF_CLCR = 77.21),
  dapto_garreau = c(TVCL = 0.365, TVV1 = 3.59, TVQ = 0.752, TVV2 = 4.71, RIF_V1 = -0.121),
  levo_canoui = c(TVCL = 5.57, TVV = 96.3, KA = 1.6, CLCR_CL = 0.684, AGE_CL = -0.312),
  rifampicine_jing = c(TVCL = 4.02, TVVD = 57.8, TVKA = 1.61),
  rifampicine_loos_ddi = c(TVCL = 9.03, TVV = 50, TVKA = 1.5, TVF = 0.68),
  rifampicine_marsot = c(TVKa = 1.15, TVCL_NOFUS = 13.7, TVV_NOFUS = 61.1, TVCL_FUS = 5.1, TVV_FUS = 23.8),
  ritonavir_kappelhoff_ddi = c(TVCL = 10.5, TVV = 96.6, TVKA = 0.871, TVLAG = 0.778, LPV_CL = 2.72),
  posaconazole_dolton_ddi = c(TVCL = 30.2, TVV = 1100, TVKA = 1.26, TVLAG = 1.79, PHE_CL = 7.21, RIF_CL = 7.21, FOS_CL = 1.342),
  ciclosporine_press_ddi = c(TVCL_CSA = 15, TVV1_CSA = 56, TVQ_CSA = 14, TVV2_CSA = 125, TVKA_CSA = 2, WTREF = 76, PREDNI_KA = 0.45, PREDNI_F = 0.78),
  isavuconazole_desai_ddi = c(TVCL = 2.36, TVV1 = 49.1, TVV2 = 417, TVQ = 26.6, KAMAX = 1.08, RA = 0.72, GAM1 = 4.88),
  tacrolimus_kamp_envarsus_ddi = c(TVCL = 19.6, TVV1 = 123, TVQ = 74.9, TVV2 = 500, TVKTR = 0.752, LAG = 2.29, WTREF = 70),
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
  amiodarone_pollak_ddi = c("229 L/day converted to L/h", "dxdt_PERI = (Q / TVV1) * CENT - (Q / V2) * PERI"),
  cefazoline_schmitz = c("CLCR     : 90.0", "ADD  : 0.714025", "pow(clcr_safe / 90.0, POW_RCC)"),
  cefazoline_lanoiselee = c("0.151392 0.324900", "IIV_Q   : 0.435600", "PROP : 0.014400"),
  cefepime_an = c("ADD  : 41.4736", "TVCLCRRT : 1.64"),
  cefepime_seo = c("ETA_Q  : 0.258064", "ETA_V2 : 0.160801"),
  ciclosporine_press_ddi = c("PREDNI_KA : 0.45", "PREDNI_F  : 0.78", "double f3A5   = 0.0"),
  dapto_dvorchik_IV_adults = c("4.28 // additive variance", "TVCL_DIAL : 0.269"),
  cefazoline_cheng = c("Q_WT_EXP    : 0.78", "Predicted Unbound Concentration"),
  rifampicine_marsot = c("5.089536 // additive variance", "0        // no proportional residual error", "dxdt_DEPOT =  Ktr * TR"),
  ritonavir_kappelhoff_ddi = c("ALAG_GUT = TVLAG", "LPV > 0.5 ? LPV_CL : 1.0"),
  posaconazole_dolton_ddi = c("pow(PHE_CL, PHE)", "pow(0.549, PPI)", "ALAG_GUT = TVLAG"),
  isavuconazole_desai_ddi = c("ASIAN > 0.5 ? 1.51 : TVCL", "PATIENT > 0.5 ? TVV2 : 260.0", "RA_IND"),
  tacrolimus_kamp_envarsus_ddi = c("IIV_V1 : 0.23 0.609", "pow(WT / WTREF, ALLOCL)", "LAG * MIX_LAGGED_ABS"),
  tacrolimus_woillard_ddi = c("double CL_IND =", "double CL_TAC_DDI = CL_IND", "+ CL_OTH"),
  vanco_goti = c("ETA_CL : 0.147043196", "ETA_V1 : 0.510339105", "ETA_V2 : 0.282197811", "age_safe > 60.0", "CRCL_RAW > 150.0", "published CV 39.8%"),
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
