# Experimental STP documentation only. This module is NOT active unless local opt-in
# ENABLE_TDM_GUIDE_PREVIEW=true and an explicit snapshot file are provided.
# It does not read observations, derive doses, or modify the model target fields.

# A preview file may only ever ask for human review. Any file claiming that a
# reference is applicable to a model is rejected.
TDM_GUIDE_ALLOWED_VERDICTS <- c("requires_human_context_review", "structural_context_mismatch")

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
          # A numeric reference must carry an explicit, review-only applicability verdict.
          !is.list(x$applicability) ||
          !(as.character(x$applicability$verdict %||% "") %in% TDM_GUIDE_ALLOWED_VERDICTS) ||
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

# Readable labels for the identity dimensions the guide currently leaves unstructured.
TDM_GUIDE_DIMENSION_LABELS <- list(
  analyteId = c("analyte", "analyte"),
  analyteBasis = c("base analyte parent/métabolite", "parent/metabolite basis"),
  matrixId = c("matrice", "matrix"),
  formulationId = c("formulation", "formulation"),
  assayMethod = c("méthode de dosage", "assay method"),
  samplingTime = c("horaire de prélèvement", "sampling time"),
  transplanted_organ = c("organe greffé", "grafted organ"),
  administration_mode = c("mode d'administration", "administration mode"),
  transplanted_organ_not_declared_by_model = c(
    "organe greffé non déclaré par le modèle", "grafted organ not declared by the model")
)

# The guide indication codes and the catalogue population tags are both closed
# vocabularies, so the conflicting contexts can be named in the reader's language.
TDM_GUIDE_ORGAN_LABELS <- list(
  kidney = c("greffe rénale", "kidney transplant"),
  liver = c("greffe hépatique", "liver transplant"),
  heart = c("greffe cardiaque", "heart transplant"),
  lung = c("greffe pulmonaire", "lung transplant")
)

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
  population <- if (english) rec$populationEn %||% rec$population else rec$population
  mismatch_n <- as.integer(r$structuralMismatchCount %||% 0L)
  mode_n <- as.integer(r$administrationModeWarningCount %||% 0L)

  label_dimension <- function(key) {
    entry <- TDM_GUIDE_DIMENSION_LABELS[[key]]
    if (is.null(entry)) key else if (english) entry[[2]] else entry[[1]]
  }
  # A context may combine several declared organs ("kidney+heart").
  label_context <- function(value) {
    parts <- strsplit(as.character(value), "+", fixed = TRUE)[[1]]
    if (!length(parts)) return(as.character(value))
    paste(vapply(parts, function(part) {
      entry <- TDM_GUIDE_ORGAN_LABELS[[part]]
      if (is.null(entry)) part else if (english) entry[[2]] else entry[[1]]
    }, character(1)), collapse = " + ")
  }
  # Dimensions missing for every reference are stated once instead of on each line.
  undecidable_of <- function(item)
    vapply(item$applicability$undecidable %||% list(), as.character, character(1))
  common_undecidable <- if (length(candidates))
    Reduce(intersect, lapply(candidates, undecidable_of)) else character(0)
  format_value <- function(item) {
    vmin <- item$min %||% NULL
    vmax <- item$max %||% NULL
    relation <- item$relation %||% "range"
    out <- if (!is.numeric(vmin) || length(vmin) != 1L || !is.finite(vmin)) "" else {
      lower <- as.character(vmin)
      if (identical(relation, "range") && is.numeric(vmax) && length(vmax) == 1L && is.finite(vmax))
        paste0(lower, "–", as.character(vmax))
      else if (identical(relation, "above")) paste0("≥", lower)
      else if (identical(relation, "at_most")) paste0("≤", lower)
      else if (relation %in% c("target", "approx")) paste0("≈", lower)
      else lower
    }
    if (nzchar(out)) paste(out, item$unit %||% "") else ""
  }
  make_ref <- function(item) {
    # Numeric values remain bibliographic candidates, never patient-specific bounds.
    formatted_value <- format_value(item)
    label <- paste0(item$metric %||% "", " · ",
      if (english) item$indicationEn %||% "" else item$indicationFr %||% "")
    context <- if (english) item$contextEn %||% "" else item$contextFr %||% ""
    mismatches <- item$applicability$mismatches %||% list()
    undecidable <- setdiff(undecidable_of(item), common_undecidable)
    shiny::tags$li(
      class = if (length(mismatches)) "tdm-guide-ref tdm-guide-ref-mismatch" else "tdm-guide-ref",
      shiny::strong(label),
      if (nzchar(formatted_value)) shiny::tags$span(paste0(" · ", formatted_value)),
      if (nzchar(context)) shiny::tags$span(paste0(" — ", context)),
      shiny::tags$span(translate(
        " · plage bibliographique candidate (non applicable automatiquement)",
        " · bibliographic candidate (not automatically applicable)")),
      # The declared context of the reference contradicts the declared context of the model.
      if (length(mismatches)) shiny::tags$div(
        class = "tdm-guide-mismatch",
        lapply(mismatches, function(m) shiny::tags$div(shiny::strong(translate(
          paste0("⚠ Contexte déclaré incompatible avec ce modèle : ",
            label_dimension(as.character(m$dimension %||% "")), " — référence « ",
            label_context(m$referenceContext %||% ""), " » vs modèle « ",
            label_context(m$modelContext %||% ""), " »."),
          paste0("⚠ Declared context conflicts with this model: ",
            label_dimension(as.character(m$dimension %||% "")), " — reference “",
            label_context(m$referenceContext %||% ""), "” vs model “",
            label_context(m$modelContext %||% ""), "”.")))))
      ),
      if (length(undecidable)) shiny::tags$div(
        class = "tdm-guide-undecidable",
        paste0(translate("Non vérifiable faute de métadonnée publiée : ",
          "Cannot be checked, metadata not published: "),
          paste(vapply(undecidable, label_dimension, character(1)), collapse = ", "), ".")
      ),
      if (grepl("^https://", item$source$url %||% ""))
        shiny::tags$a(href = item$source$url, target = "_blank",
          rel = "noopener noreferrer", translate(" Source", " Source"))
    )
  }
  make_pd <- function(item) {
    label <- if (english) item$criterionEn %||% item$metric else item$criterionFr %||% item$metric
    mode_mismatch <- identical(as.character(item$status %||% ""), "different_administration_mode")
    shiny::tags$li(
      class = if (mode_mismatch) "tdm-guide-ref tdm-guide-ref-mismatch" else "tdm-guide-ref",
      shiny::span(label),
      shiny::tags$span(translate(
        " · objectif contextuel (CMI, fraction libre, perfusion à vérifier)",
        " · contextual objective (MIC, unbound fraction and infusion must be checked)")),
      # Previously computed but never shown: a continuous-infusion objective must not
      # be read silently on an intermittent model, and vice versa.
      if (mode_mismatch) shiny::tags$div(
        class = "tdm-guide-mismatch",
        shiny::strong(translate(
          "⚠ Objectif établi pour un autre mode d’administration que celui du modèle.",
          "⚠ Objective established for an administration mode other than this model's."))
      ),
      if (grepl("^https://", item$source$url %||% ""))
        shiny::tags$a(href = item$source$url, target = "_blank",
          rel = "noopener noreferrer", translate(" Source", " Source"))
    )
  }

  warning_banner <- if (mismatch_n > 0L || mode_n > 0L) shiny::tags$p(
    class = "tdm-guide-warning-banner",
    shiny::strong(translate("⚠ À vérifier avant toute lecture : ",
      "⚠ Check before reading any value: ")),
    paste0(paste(c(
      if (mismatch_n > 0L) translate(
        paste0(mismatch_n, " référence(s) au contexte déclaré incompatible avec ce modèle"),
        paste0(mismatch_n, " reference(s) whose declared context conflicts with this model")),
      if (mode_n > 0L) translate(
        paste0(mode_n, " objectif(s) PK/PD d’un autre mode d’administration"),
        paste0(mode_n, " PK/PD objective(s) from another administration mode"))
    ), collapse = " ; "), ".")
  ) else NULL

  # Pedagogical order: the index and the sampling answer the question immediately;
  # numbers, sources and limits stay one click away.
  shiny::tags$div(
    class = "tdm-guide-preview-readonly",
    `data-model-id` = model_id,
    shiny::tags$p(class = "tdm-guide-preview-title",
      shiny::strong(translate("Recommandation STP · prototype bibliographique",
        "TDM guidance · bibliographic preview"))),
    shiny::tags$p(shiny::tags$strong(translate("Indice PK/PD : ", "PK/PD index: ")),
      index %||% translate("Non établi", "Not established")),
    shiny::tags$p(shiny::tags$strong(translate("Prélèvement pratique : ", "Practical sampling: ")),
      sampling %||% translate("À documenter", "Not documented")),
    if (!is.null(population)) shiny::tags$p(class = "tdm-guide-model-population",
      shiny::tags$strong(translate("Population du modèle : ", "Model population: ")), population),
    warning_banner,
    shiny::tags$details(
      shiny::tags$summary(translate("Références, contexte et limites d’applicabilité",
        "References, context and applicability limits")),
      if (!is.null(monitoring)) shiny::tags$p(if (english) monitoring$noteEn else monitoring$noteFr),
      if (length(candidates)) shiny::tags$div(
        shiny::tags$strong(translate("Références publiées possibles : ",
          "Published candidate references: ")),
        # Stated once: what the guide does not publish for any of these references.
        if (length(common_undecidable)) shiny::tags$p(
          class = "tdm-guide-undecidable tdm-guide-undecidable-common",
          paste0(translate(
            "Pour toutes les références ci-dessous, le guide ne publie pas encore : ",
            "For every reference below, the guide does not yet publish: "),
            paste(vapply(common_undecidable, label_dimension, character(1)), collapse = ", "),
            translate(". La comparaison de ces dimensions avec le modèle reste donc à faire par un pharmacologue.",
              ". Comparing these dimensions with the model therefore remains a pharmacologist's task."))
        ),
        shiny::tags$ul(lapply(candidates, make_ref))
      ),
      if (length(pd)) shiny::tags$div(
        shiny::tags$strong(translate("Repères PK/PD conditionnels : ",
          "Conditional PK/PD objectives: ")),
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
  )
}
