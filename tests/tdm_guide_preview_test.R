# Standalone R/Shiny presentation contract: no patient, dose or target mutation.
root <- normalizePath(file.path(getwd()), mustWork=TRUE)
APP_ROOT <- root
`%||%` <- function(x,y) if(is.null(x)||!length(x)) y else x
catalog <- jsonlite::fromJSON(file.path(root,"models","catalog.json"))$models
MODEL_CATALOG <- catalog
Sys.unsetenv("ENABLE_TDM_GUIDE_PREVIEW")
source(file.path(root,"R","tdm_guide_preview.R"),local=FALSE)
stopifnot(is.null(TDM_GUIDE_PREVIEW_BY_MODEL))
stopifnot(is.null(tdm_guide_preview_model_card(catalog$id[[1]], english=FALSE)))
testmodels <- lapply(seq_len(nrow(catalog)), function(i) {
  list(modelId=catalog$id[[i]],recommendation=list(
    interpretationAllowed=FALSE,automaticDoseAdjustmentAllowed=FALSE,
    patientContextMatched=FALSE,modelFormulationMatched=FALSE,
    analyteAndMatrixMatched=FALSE,assayMethodMatched=FALSE,
    evidenceSourceVersion=paste(rep("a",64),collapse=""),
    index=list(fr="AUC",en="AUC"),sampling=list(fr="C0",en="Trough"),
    numericalReferenceCandidates=list(),pharmacodynamicCandidates=list()))
})
fixture <- list(schema="mipd.tdm.model-recommendations.preview.v1",
  enabledInRuntime=FALSE,patientDataIncluded=FALSE,
  clinicalInterpretationAllowed=FALSE,automaticDoseAdjustmentAllowed=FALSE,
  sourceSnapshotSha256=paste(rep("a",64),collapse=""),models=testmodels)
tmp <- tempfile(fileext=".json")
on.exit(unlink(tmp), add=TRUE)
jsonlite::write_json(fixture,tmp,auto_unbox=TRUE,null="null",pretty=TRUE)
parsed <- tdm_guide_preview_load(tmp,catalog)
stopifnot(length(parsed)==nrow(catalog))
TDM_GUIDE_PREVIEW_BY_MODEL <- parsed
html <- as.character(htmltools::renderTags(tdm_guide_preview_model_card(catalog$id[[1]],english=FALSE))$html)
stopifnot(grepl("Recommandation STP",html,fixed=TRUE))
stopifnot(grepl("Aucune cible chiffrée",html,fixed=TRUE))
stopifnot(!grepl("automaticDoseAdjustmentAllowed",html,fixed=TRUE))
stopifnot(grepl("TDM guidance",as.character(htmltools::renderTags(
  tdm_guide_preview_model_card(catalog$id[[1]],english=TRUE))$html),fixed=TRUE))
with_range <- fixture
with_range$models[[1]]$recommendation$numericalReferenceCandidates <- list(list(
  id="test-synthetic",metric="Cmax",unit="mg/L",relation="range",min=10,max=15,
  indicationFr="Essai",indicationEn="Study",contextFr="Contexte spécifique",
  contextEn="Specific context",
  source=list(url="https://example.org/source"),
  applicability=list(verdict="requires_human_context_review",
    mismatches=list(),undecidable=list("formulationId","assayMethod")),
  modelApplicabilityVerified=FALSE,automaticDoseAdjustmentAllowed=FALSE
))
jsonlite::write_json(with_range,tmp,auto_unbox=TRUE,null="null")
TDM_GUIDE_PREVIEW_BY_MODEL <- tdm_guide_preview_load(tmp,catalog)
html_range <- as.character(htmltools::renderTags(
  tdm_guide_preview_model_card(catalog$id[[1]],english=FALSE))$html)
stopifnot(grepl("10–15 mg/L",html_range,fixed=TRUE))
stopifnot(grepl("Contexte spécifique",html_range,fixed=TRUE))
stopifnot(grepl("non applicable automatiquement",html_range,fixed=TRUE))
unsafe <- fixture
unsafe$models[[1]]$recommendation$interpretationAllowed <- TRUE
jsonlite::write_json(unsafe,tmp,auto_unbox=TRUE,null="null")
stopifnot(inherits(try(tdm_guide_preview_load(tmp,catalog),silent=TRUE),"try-error"))
unsafe <- fixture
unsafe$models <- unsafe$models[-1]
jsonlite::write_json(unsafe,tmp,auto_unbox=TRUE,null="null")
stopifnot(inherits(try(tdm_guide_preview_load(tmp,catalog),silent=TRUE),"try-error"))
cat("Experimental Shiny STP panels: 46/46 models, disabled by default, safe read-only display.\n")

# --- Applicability warnings must reach the panel, and none may claim applicability ---
# The undecidable dimensions declared by the adapter are shown, not dropped.
stopifnot(grepl("formulation", html_range, fixed=TRUE))
stopifnot(grepl("éthode de dosage", html_range, fixed=TRUE))
# Key information stays visible: index and sampling precede the collapsed section.
stopifnot(regexpr("Indice PK/PD", html_range, fixed=TRUE) <
  regexpr("<details", html_range, fixed=TRUE))
stopifnot(regexpr("èvement pratique", html_range, fixed=TRUE) <
  regexpr("<details", html_range, fixed=TRUE))

# A structural context mismatch is rendered as a visible warning plus a counted banner.
mismatched <- fixture
mismatched$models[[1]]$recommendation$structuralMismatchCount <- 1L
mismatched$models[[1]]$recommendation$numericalReferenceCandidates <- list(list(
  id="test-organ-mismatch", metric="C0", unit="ng/mL", relation="range", min=3, max=5,
  indicationFr="Greffe hépatique", indicationEn="Liver transplant",
  contextFr="Contexte hépatique", contextEn="Liver context",
  source=list(url="https://example.org/source"),
  applicability=list(verdict="structural_context_mismatch",
    mismatches=list(list(dimension="transplanted_organ", code="organ_indication_mismatch",
      referenceContext="liver", modelContext="kidney")),
    undecidable=list()),
  modelApplicabilityVerified=FALSE, automaticDoseAdjustmentAllowed=FALSE))
jsonlite::write_json(mismatched, tmp, auto_unbox=TRUE, null="null")
TDM_GUIDE_PREVIEW_BY_MODEL <- tdm_guide_preview_load(tmp, catalog)
html_mismatch <- as.character(htmltools::renderTags(
  tdm_guide_preview_model_card(catalog$id[[1]], english=FALSE))$html)
stopifnot(grepl("Contexte déclaré incompatible", html_mismatch, fixed=TRUE))
stopifnot(grepl("greffe hépatique", html_mismatch, fixed=TRUE))
stopifnot(grepl("greffe rénale", html_mismatch, fixed=TRUE))
stopifnot(grepl("tdm-guide-warning-banner", html_mismatch, fixed=TRUE))
stopifnot(grepl("tdm-guide-ref-mismatch", html_mismatch, fixed=TRUE))
# The banner is visible without opening the collapsed section.
stopifnot(regexpr("tdm-guide-warning-banner", html_mismatch, fixed=TRUE) <
  regexpr("<details", html_mismatch, fixed=TRUE))

# A PK/PD objective from another administration mode is no longer silently equivalent.
pd_mode <- fixture
pd_mode$models[[1]]$recommendation$administrationModeWarningCount <- 1L
pd_mode$models[[1]]$recommendation$pharmacodynamicCandidates <- list(list(
  id="pd-continuous", kind="beta_continuous", metric="fCss/MIC",
  criterionFr="4 CMI en perfusion continue", criterionEn="4 MIC under continuous infusion",
  status="different_administration_mode", source=list(url="https://example.org/pd"),
  modelApplicabilityVerified=FALSE, automaticDoseAdjustmentAllowed=FALSE))
jsonlite::write_json(pd_mode, tmp, auto_unbox=TRUE, null="null")
TDM_GUIDE_PREVIEW_BY_MODEL <- tdm_guide_preview_load(tmp, catalog)
html_pd <- as.character(htmltools::renderTags(
  tdm_guide_preview_model_card(catalog$id[[1]], english=FALSE))$html)
stopifnot(grepl("autre mode d’administration", html_pd, fixed=TRUE))
stopifnot(grepl("tdm-guide-ref-mismatch", html_pd, fixed=TRUE))

# Fail closed: a numeric reference without a verdict, or claiming applicability, is rejected.
no_verdict <- with_range
no_verdict$models[[1]]$recommendation$numericalReferenceCandidates[[1]]$applicability <- NULL
jsonlite::write_json(no_verdict, tmp, auto_unbox=TRUE, null="null")
stopifnot(inherits(try(tdm_guide_preview_load(tmp, catalog), silent=TRUE), "try-error"))
claims <- with_range
claims$models[[1]]$recommendation$numericalReferenceCandidates[[1]]$applicability$verdict <- "applicable"
jsonlite::write_json(claims, tmp, auto_unbox=TRUE, null="null")
stopifnot(inherits(try(tdm_guide_preview_load(tmp, catalog), silent=TRUE), "try-error"))
cat("Shiny STP applicability warnings: organ and administration-mode conflicts visible, applicability claims rejected.
")
