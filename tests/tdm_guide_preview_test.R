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
unsafe <- fixture
unsafe$models[[1]]$recommendation$interpretationAllowed <- TRUE
jsonlite::write_json(unsafe,tmp,auto_unbox=TRUE,null="null")
stopifnot(inherits(try(tdm_guide_preview_load(tmp,catalog),silent=TRUE),"try-error"))
unsafe <- fixture
unsafe$models <- unsafe$models[-1]
jsonlite::write_json(unsafe,tmp,auto_unbox=TRUE,null="null")
stopifnot(inherits(try(tdm_guide_preview_load(tmp,catalog),silent=TRUE),"try-error"))
cat("Experimental Shiny STP panels: 46/46 models, disabled by default, safe read-only display.\n")
