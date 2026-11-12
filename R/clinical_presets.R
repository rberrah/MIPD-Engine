# Historic MIPD curated simulation presets now come from the 23 per-drug JSON dossiers.
# This is the existing educational behavior, not automatic clinical prescribing.
if (!exists("TDM_DRUG_PROFILES", inherits = TRUE)) {
  source(file.path(APP_ROOT, "R", "tdm_drug_profiles.R"), local = FALSE)
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
