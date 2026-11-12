# Drug-specific reference dossiers used by the Shiny analysis interface.
# One committed JSON file per library drugKey; no network dependency, no patient data.
if (!exists("APP_ROOT", inherits = TRUE)) {
  APP_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
TDM_DRUG_PROFILE_DIR <- file.path(APP_ROOT, "integration", "tdm-drugs")

tdm_drug_profiles_load <- function(directory = TDM_DRUG_PROFILE_DIR, catalog = MODEL_CATALOG) {
  paths <- list.files(directory, pattern = "\\.json$", full.names = TRUE)
  required_keys <- sort(unique(as.character(catalog$drugKey)))
  if (length(paths) != length(required_keys)) {
    stop("TDM drug profiles: expected exactly one JSON file for each library drug")
  }
  result <- list()
  for (path in paths) {
    d <- jsonlite::fromJSON(path, simplifyVector = FALSE)
    key <- as.character(d$drugKey %||% "")
    if (!identical(d$schema, "mipd.tdm.drug-profile.v1") ||
        !grepl("^[a-z][a-z0-9_-]*$", key) ||
        !identical(tools::file_path_sans_ext(basename(path)), key) ||
        !key %in% required_keys || !is.null(result[[key]]) ||
        !identical(d$intendedUse, "TEACHING_RESEARCH_ONLY") ||
        !identical(d$clinicalInterpretationAllowed, FALSE) ||
        !identical(d$automaticDoseAdjustmentAllowed, FALSE) ||
        !identical(d$patientDataIncluded, FALSE)) {
      stop("Invalid or unsafe TDM drug profile: ", basename(path))
    }
    model_rows <- catalog[as.character(catalog$drugKey) == key, , drop = FALSE]
    ids <- as.character(unlist(d$modelIds %||% list(), use.names = FALSE))
    if (anyDuplicated(ids) || !setequal(ids, as.character(model_rows$id)) ||
        !identical(as.character(d$concentrationUnit), unique(as.character(model_rows$concentrationUnit))) ||
        !identical(d$doseUnit, "mg") ||
        !grepl("^https://rberrah\\.github\\.io/pharmacometrie/cibles-tdm/\\?drug=", d$guideHref %||% "") ||
        !grepl("^[a-f0-9]{64}$", d$sourceSnapshotSha256 %||% "")) {
      stop("TDM drug profile identity or units mismatch: ", key)
    }
    if (!is.list(d$guideTargetPresets) || !is.list(d$simulationPresets) ||
        !is.list(d$sourceReferences) || !is.list(d$conditionalPdObjectives)) {
      stop("TDM drug profile has invalid reference lists: ", key)
    }
    for (entry in d$guideTargetPresets) {
      if (!is.list(entry$candidates)) stop("Invalid TDM target candidates: ", key)
      for (c in entry$candidates) {
        if (!identical(c$modelApplicabilityVerified, FALSE) ||
            !identical(c$automaticDoseAdjustmentAllowed, FALSE) ||
            !grepl("^https://", c$url %||% "") ||
            !is.numeric(c$target_low) || !is.numeric(c$target_high) ||
            c$target_low >= c$target_high) {
          stop("Unqualified TDM target in: ", key)
        }
      }
      if (!is.null(entry$doseGrid) &&
          (!identical(entry$doseGrid$signoff, "PENDING") ||
           !grepl("^https://", entry$doseGrid$source$url %||% ""))) {
        stop("Unreviewed TDM dose search grid in: ", key)
      }
    }
    result[[key]] <- d
  }
  if (!setequal(names(result), required_keys)) stop("TDM drug profile catalogue is incomplete")
  result
}

TDM_DRUG_PROFILES <- tdm_drug_profiles_load()

tdm_drug_profile <- function(record) {
  key <- as.character(record$drugKey[[1]] %||% "")
  TDM_DRUG_PROFILES[[key]]
}

tdm_drug_profile_guide_entry <- function(record, mode) {
  profile <- tdm_drug_profile(record)
  if (is.null(profile)) return(NULL)
  profile$guideTargetPresets[[as.character(mode %||% "")]]
}

tdm_drug_profile_search_grid <- function(record, mode) {
  entry <- tdm_drug_profile_guide_entry(record, mode)
  if (is.null(entry)) return(NULL)
  entry$doseGrid
}

tdm_drug_profile_metric_unit <- function(record, metric) {
  profile <- tdm_drug_profile(record)
  unit <- if (is.null(profile)) as.character(record$concentrationUnit[[1]] %||% "") else profile$concentrationUnit
  if (identical(metric, "TimeAbove")) return("%")
  if (identical(metric, "AUC24")) {
    if (identical(unit, "ng/mL")) return("ng·h/mL")
    if (identical(unit, "mg/L")) return("mg·h/L")
    return(paste0(unit, "·h"))
  }
  unit
}

tdm_drug_profile_formulation <- function(record, st = NA_real_) {
  profile <- tdm_drug_profile(record)
  id <- as.character(record$id[[1]] %||% "")
  data <- profile$formulationCovariates[[id]] %||% NULL
  if (is.null(data)) return(NULL)
  if (!is.null(data$values)) {
    if (!is.finite(st) || !st %in% c(0, 1)) return(NULL)
    return(data$values[[as.character(as.integer(st))]])
  }
  data
}

# A label-derived grid can be copied into a simulation ONLY on an explicit
# user action, after confirming its product/formulation. No clinical sign-off.
tdm_drug_profile_grid_eligible <- function(record, mode, confirmed_formulation, st = NA_real_) {
  grid <- tdm_drug_profile_search_grid(record, mode)
  if (is.null(grid) || !identical(grid$signoff, "PENDING")) return(FALSE)
  if (!is.character(confirmed_formulation) || length(confirmed_formulation) != 1L ||
      !nzchar(confirmed_formulation)) return(FALSE)
  allowed <- as.character(unlist(grid$requiresFormulationIds %||% list(), use.names = FALSE))
  if (!confirmed_formulation %in% allowed) return(FALSE)
  id <- as.character(record$id[[1]])
  if (identical(id, "tacrolimus_kamp_envarsus_ddi")) return(FALSE)
  if (identical(id, "tacrolimus_woillard_ddi")) {
    return(identical(confirmed_formulation, "tacrolimus-immediate-release") &&
             is.finite(st) && identical(as.numeric(st), 1))
  }
  # An oncology/transplant mixed source cannot receive a transplant-only grid.
  if (identical(id, "everolimus_terheine_ddi")) return(FALSE)
  declared <- tdm_drug_profile(record)$formulationByModel[[id]]$formulationId %||% ""
  if (nzchar(declared) && !declared %in% c("not-declared-by-source") &&
      !declared %in% allowed) return(FALSE)
  TRUE
}
