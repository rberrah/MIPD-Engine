# Shiny-only JSON runtime bridge. Keep the original imported files immutable.
# Source this after R/clinical_presets.R and R/tdm_drug_profiles.R.
if (!exists("TDM_DRUG_PROFILES", inherits = TRUE)) {
  stop("R/tdm_drug_profiles.R must be loaded before runtime drug presets")
}
CLINICAL_TARGET_PRESETS <- list()
for (drug_key in names(TDM_DRUG_PROFILES)) {
  presets <- TDM_DRUG_PROFILES[[drug_key]]$simulationPresets
  for (mode in names(presets)) {
    preset <- presets[[mode]]
    preset$intervals <- as.numeric(unlist(preset$intervals, use.names = FALSE))
    CLINICAL_TARGET_PRESETS[[paste(drug_key, mode, sep = "::")]] <- preset
  }
}
clinical_target_preset <- function(record, mode) {
  drug_key <- as.character(record$drugKey[[1]] %||% "")
  CLINICAL_TARGET_PRESETS[[paste(drug_key, mode, sep = "::")]]
}
