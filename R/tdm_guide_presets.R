# Target presets derived from the PMx TDM Guide, read from a versioned artifact.
# MIPD keeps its own curated presets in clinical_presets.R; those stay authoritative.
# This layer only ever OFFERS a published target. It writes nothing on its own:
# a reference is applied when the reader selects it, or, when a molecule and route
# have exactly one qualified reference, by the existing curated preset of that key.
if (!exists("APP_ROOT", inherits = TRUE)) {
  APP_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

TDM_GUIDE_PRESET_FILE <- file.path(APP_ROOT, "integration", "tdm-guide-clinical-presets.v1.json")
TDM_GUIDE_PRESET_METRICS <- c("AUC24", "Cmin", "Cmax")
# A dose grid is a posology statement: only a signed-off grid may reach the inputs.
TDM_GUIDE_GRID_SIGNOFFS <- c("PENDING", "APPROVED")

# Fail closed: anything unexpected disables the guide layer instead of guessing.
tdm_guide_presets_load <- function(path = TDM_GUIDE_PRESET_FILE) {
  if (!nzchar(path) || !file.exists(path)) return(NULL)
  object <- try(jsonlite::fromJSON(path, simplifyVector = FALSE), silent = TRUE)
  if (inherits(object, "try-error")) return(NULL)
  ok <- identical(object$schema, "mipd.tdm.clinical-presets.v1") &&
    identical(object$clinicalInterpretationAllowed, FALSE) &&
    identical(object$automaticDoseAdjustmentAllowed, FALSE) &&
    identical(object$patientDataIncluded, FALSE) &&
    is.character(object$sourceSnapshotSha256) &&
    grepl("^[a-f0-9]{64}$", object$sourceSnapshotSha256) &&
    is.list(object$presets)
  if (!ok) return(NULL)
  for (entry in object$presets) {
    if (!is.list(entry) || !is.list(entry$candidates) || !length(entry$candidates) ||
        !is.logical(entry$autoApply) || length(entry$autoApply) != 1L) return(NULL)
    # A single qualified reference is the only case the artifact may mark auto-applicable.
    if (!identical(entry$autoApply, length(entry$candidates) == 1L)) return(NULL)
    for (candidate in entry$candidates) {
      if (!identical(candidate$modelApplicabilityVerified, FALSE) ||
          !identical(candidate$automaticDoseAdjustmentAllowed, FALSE) ||
          !(as.character(candidate$target_metric %||% "") %in% TDM_GUIDE_PRESET_METRICS) ||
          !is.numeric(candidate$target_low) || !is.numeric(candidate$target_high) ||
          !is.finite(candidate$target_low) || !is.finite(candidate$target_high) ||
          candidate$target_high <= candidate$target_low ||
          !grepl("^https://", as.character(candidate$url %||% ""))) return(NULL)
    }
    grid <- entry$doseGrid
    if (!is.null(grid)) {
      if (!(as.character(grid$signoff %||% "") %in% TDM_GUIDE_GRID_SIGNOFFS) ||
          !is.numeric(grid$dose_min) || !is.numeric(grid$dose_max) ||
          !is.numeric(grid$dose_step) || grid$dose_max <= grid$dose_min ||
          grid$dose_step <= 0 || !length(grid$intervals) ||
          !grepl("^https://", as.character(grid$source$url %||% ""))) return(NULL)
      # The applicable flag must agree with the signoff it claims to summarise.
      if (!identical(entry$doseGridApplicable, identical(grid$signoff, "APPROVED"))) return(NULL)
    } else if (!identical(entry$doseGridApplicable, FALSE)) return(NULL)
  }
  object
}

TDM_GUIDE_PRESETS <- tdm_guide_presets_load()

tdm_guide_preset_key <- function(record, mode) {
  paste(as.character(record$drugKey[[1]] %||% ""), as.character(mode %||% ""), sep = "::")
}

#' Published references available for this model and route, or NULL.
tdm_guide_preset_entry <- function(record, mode) {
  if (is.null(TDM_GUIDE_PRESETS)) return(NULL)
  TDM_GUIDE_PRESETS$presets[[tdm_guide_preset_key(record, mode)]]
}

#' TRUE only when a signed-off grid exists for this model and route.
#' While the grid is PENDING it is displayed but never written to the inputs.
tdm_guide_preset_grid_applicable <- function(record, mode) {
  entry <- tdm_guide_preset_entry(record, mode)
  isTRUE(entry$doseGridApplicable)
}

tdm_guide_preset_candidate <- function(record, mode, id) {
  entry <- tdm_guide_preset_entry(record, mode)
  if (is.null(entry) || !nzchar(id %||% "")) return(NULL)
  match <- Filter(function(candidate) identical(as.character(candidate$id), as.character(id)),
    entry$candidates)
  if (!length(match)) NULL else match[[1]]
}

#' Published references that contradict MIPD's own curated preset for the same key.
#' Reported, never silently resolved: the curated value keeps precedence.
tdm_guide_preset_divergence <- function(record, mode) {
  entry <- tdm_guide_preset_entry(record, mode)
  if (is.null(entry)) return(list())
  curated <- clinical_target_preset(record, mode)
  if (is.null(curated)) return(list())
  Filter(function(candidate) {
    identical(as.character(candidate$target_metric), as.character(curated$target_metric)) &&
      !(isTRUE(all.equal(candidate$target_low, curated$target_low)) &&
        isTRUE(all.equal(candidate$target_high, curated$target_high)))
  }, entry$candidates)
}

tdm_guide_preset_choices <- function(record, mode, english = FALSE) {
  entry <- tdm_guide_preset_entry(record, mode)
  if (is.null(entry)) return(character())
  ids <- vapply(entry$candidates, function(candidate) as.character(candidate$id), character(1))
  labels <- vapply(entry$candidates, function(candidate) {
    scope <- if (english) candidate$labelEn %||% "" else candidate$labelFr %||% ""
    paste0(
      if (nzchar(scope)) paste0(scope, " · ") else "",
      candidate$target_metric, " ",
      format(candidate$target_low), "–", format(candidate$target_high), " ",
      candidate$target_unit %||% ""
    )
  }, character(1))
  stats::setNames(ids, labels)
}
