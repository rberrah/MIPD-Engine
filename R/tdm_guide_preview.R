# Experimental STP documentation only. This module is NOT active unless local opt-in
# ENABLE_TDM_GUIDE_PREVIEW=true and an explicit snapshot file are provided.
# It does not read observations, derive doses, or modify the model target fields.
tdm_guide_preview_load <- function(path, catalog) {
  if (!nzchar(path) || !file.exists(path)) stop("Missing explicit local STP preview JSON")
  obj <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  fail <- function() stop("Unsafe or incompatible STP preview file")
  if (!identical(obj$schema, "mipd.tdm.model-recommendations.preview.v1") ||
      !identical(obj$enabledInRuntime, FALSE) ||
      !identical(obj$patientDataIncluded, FALSE) ||
      !identical(obj$clinicalInterpretationAllowed, FALSE) ||
      !identical(obj$automaticDoseAdjustmentAllowed, FALSE) ||
      is.null(obj$sourceSnapshotSha256) ||
      !grepl("^[a-f0-9]{64}$", obj$sourceSnapshotSha256) ||
      !is.list(obj$models)) fail()
  ids <- vapply(obj$models, function(x) as.character(x$modelId %||% ""), character(1))
  expected <- as.character(catalog$id)
  if (length(ids) != length(expected) || anyDuplicated(ids) ||
      !setequal(ids, expected)) fail()
  for (model in obj$models) {
    r <- model$recommendation
    if (is.null(r) || !identical(r$interpretationAllowed, FALSE) ||
        !identical(r$automaticDoseAdjustmentAllowed, FALSE) ||
        !identical(r$patientContextMatched, FALSE) ||
        !identical(r$modelFormulationMatched, FALSE) ||
        !identical(r$analyteAndMatrixMatched, FALSE) ||
        !identical(r$assayMethodMatched, FALSE) ||
        !identical(r$evidenceSourceVersion, obj$sourceSnapshotSha256)) fail()
    refs <- r$numericalReferenceCandidates %||% list()
    pd <- r$pharmacodynamicCandidates %||% list()
    if (!is.list(refs) || !is.list(pd) ||
        any(vapply(refs, function(x) !identical(x$modelApplicabilityVerified, FALSE) ||
          !identical(x$automaticDoseAdjustmentAllowed, FALSE) ||
          !grepl("^https://", as.character(x$source$url %||% "")), logical(1))) ||
        any(vapply(pd, function(x) !identical(x$modelApplicabilityVerified, FALSE) ||
          !identical(x$automaticDoseAdjustmentAllowed, FALSE), logical(1)))) fail()
  }
  stats::setNames(obj$models, ids)
}

tdm_guide_preview_available <- function() {
  identical(tolower(trimws(Sys.getenv("ENABLE_TDM_GUIDE_PREVIEW", "false"))), "true")
}
# Only explicit local opt-in can load a file. The public deployment remains off.
TDM_GUIDE_PREVIEW_BY_MODEL <- if (tdm_guide_preview_available()) {
  tdm_guide_preview_load(Sys.getenv("TDM_GUIDE_PREVIEW_JSON", ""), MODEL_CATALOG)
} else NULL

tdm_guide_preview_model_card <- function(model_id, english = FALSE) {
  if (is.null(TDM_GUIDE_PREVIEW_BY_MODEL)) return(NULL)
  rec <- TDM_GUIDE_PREVIEW_BY_MODEL[[model_id]]
  if (is.null(rec)) return(NULL)
  r <- rec$recommendation
  translate <- function(fr, en) if (english) en else fr
  index <- if (english) r$index$en else r$index$fr
  sampling <- if (english) r$sampling$en else r$sampling$fr
  candidates <- r$numericalReferenceCandidates %||% list()
  pd <- r$pharmacodynamicCandidates %||% list()
  monitoring <- r$monitoringContext
  make_ref <- function(item) {
    # Numeric values remain bibliographic candidates, never patient-specific bounds.
    vmin <- item$min %||% NULL
    vmax <- item$max %||% NULL
    relation <- item$relation %||% "range"
    formatted_value <- if (!is.numeric(vmin) || length(vmin)!=1L || !is.finite(vmin)) "" else {
      lower <- as.character(vmin)
      if (identical(relation,"range") && is.numeric(vmax) && length(vmax)==1L && is.finite(vmax))
        paste0(lower,"–",as.character(vmax))
      else if (identical(relation,"above")) paste0("≥",lower)
      else if (identical(relation,"at_most")) paste0("≤",lower)
      else if (relation %in% c("target","approx")) paste0("≈",lower)
      else lower
    }
    if (nzchar(formatted_value)) formatted_value <- paste(formatted_value,item$unit %||% "")
    label <- paste0(item$metric %||% "", " · ",
      if (english) item$indicationEn %||% "" else item$indicationFr %||% "")
    context <- if (english) item$contextEn %||% "" else item$contextFr %||% ""
    shiny::tags$li(
      shiny::strong(label),
      if (nzchar(formatted_value)) shiny::tags$span(paste0(" · ", formatted_value)),
      if (nzchar(context)) shiny::tags$span(paste0(" — ", context)),
      shiny::tags$span(translate(
        " · plage bibliographique candidate (non applicable automatiquement)",
        " · bibliographic candidate (not automatically applicable)")),
      if (grepl("^https://", item$source$url %||% ""))
        shiny::tags$a(href = item$source$url, target = "_blank",
          rel = "noopener noreferrer", translate(" Source", " Source"))
    )
  }
  make_pd <- function(item) {
    label <- if (english) item$criterionEn %||% item$metric else item$criterionFr %||% item$metric
    shiny::tags$li(
      shiny::span(label),
      shiny::tags$span(translate(
        " · objectif contextuel (CMI, fraction libre, perfusion à vérifier)",
        " · contextual objective (MIC, unbound fraction and infusion must be checked)")),
      if (grepl("^https://", item$source$url %||% ""))
        shiny::tags$a(href = item$source$url, target = "_blank",
          rel = "noopener noreferrer", translate(" Source", " Source"))
    )
  }
  shiny::tags$details(
    class = "tdm-guide-preview-readonly",
    `data-model-id` = model_id,
    shiny::tags$summary(translate("Recommandation STP · prototype bibliographique", "TDM guidance · bibliographic preview")),
    shiny::tags$p(shiny::tags$strong(translate("Indice PK/PD : ", "PK/PD index: ")), index %||% translate("Non établi", "Not established")),
    shiny::tags$p(shiny::tags$strong(translate("Prélèvement pratique : ", "Practical sampling: ")), sampling %||% translate("À documenter", "Not documented")),
    if (!is.null(monitoring)) shiny::tags$p(if (english) monitoring$noteEn else monitoring$noteFr),
    if (length(candidates)) shiny::tags$div(
      shiny::tags$strong(translate("Références publiées possibles : ", "Published candidate references: ")),
      shiny::tags$ul(lapply(candidates, make_ref))
    ),
    if (length(pd)) shiny::tags$div(
      shiny::tags$strong(translate("Repères PK/PD conditionnels : ", "Conditional PK/PD objectives: ")),
      shiny::tags$ul(lapply(pd, make_pd))
    ),
    if (!length(candidates) && !length(pd)) shiny::tags$p(
      translate("Aucune cible chiffrée qualifiée ici. Ne pas inventer de plage.",
        "No qualified numerical target here. Do not infer a range.")
    ),
    shiny::tags$p(translate(
      "Références générales sans correspondance patient, modèle/formulation ni méthode analytique vérifiées. Recherche/enseignement uniquement. Aucune proposition posologique.",
      "General references; no verified patient, model/formulation or assay match. Research/education only. No dosing proposal."
    ))
  )
}
