CLINICAL_TARGET_PRESETS <- list(
  "vanco::IV_INTERMITTENT" = list(
    target_metric = "AUC24",
    target_low = 400,
    target_high = 600,
    target_unit = "mg.h/L",
    dose_min = 250,
    dose_max = 3000,
    dose_step = 250,
    intervals = c(6, 8, 12, 24),
    infusion = 1,
    scope_fr = "Infection invasive grave à SARM, CMI supposée à 1 mg/L.",
    scope_en = "Serious invasive MRSA infection, assuming an MIC of 1 mg/L.",
    grid_note_fr = "La cible est issue de la recommandation; la grille de doses reste une plage exploratoire à adapter au poids et à la fonction rénale.",
    grid_note_en = "The target comes from the guideline; the dose grid remains an exploratory range to adapt to weight and renal function.",
    citation = "Rybak MJ et al. Am J Health Syst Pharm. 2020;77:835-864.",
    doi = "10.1093/ajhp/zxaa036"
  ),
  "vanco::IV_CONTINUOUS" = list(
    target_metric = "Cmin",
    target_low = 20,
    target_high = 25,
    target_unit = "mg/L",
    dose_min = 500,
    dose_max = 5000,
    dose_step = 250,
    intervals = 24,
    infusion = 24,
    scope_fr = "Perfusion continue chez le patient de réanimation; la concentration à l'état stationnaire est représentée par Cmin.",
    scope_en = "Continuous infusion in critically ill patients; the steady-state concentration is represented by Cmin.",
    grid_note_fr = "La cible est issue de la recommandation; la grille quotidienne reste exploratoire et doit être adaptée au poids et à la fonction rénale.",
    grid_note_en = "The target comes from the guideline; the daily dose grid remains exploratory and must be adapted to weight and renal function.",
    citation = "Rybak MJ et al. Am J Health Syst Pharm. 2020;77:835-864.",
    doi = "10.1093/ajhp/zxaa036"
  ),
  "linez::IV_INTERMITTENT" = list(
    target_metric = "Cmin",
    target_low = 2,
    target_high = 8,
    target_unit = "mg/L",
    dose_min = 300,
    dose_max = 1200,
    dose_step = 100,
    intervals = c(8, 12, 24),
    infusion = 1,
    scope_fr = "Cible de concentration résiduelle du consensus TDM linezolide.",
    scope_en = "Linezolid trough target from the TDM expert consensus.",
    grid_note_fr = "La cible est issue du consensus; la grille de doses est exploratoire et ne remplace pas l'adaptation à l'indication et à la tolérance.",
    grid_note_en = "The target comes from the consensus; the dose grid is exploratory and does not replace adjustment for indication and tolerability.",
    citation = "Lin B et al. Front Public Health. 2022;10:967311.",
    doi = "10.3389/fpubh.2022.967311"
  ),
  "linez::ORAL" = list(
    target_metric = "Cmin",
    target_low = 2,
    target_high = 8,
    target_unit = "mg/L",
    dose_min = 300,
    dose_max = 1200,
    dose_step = 100,
    intervals = c(8, 12, 24),
    infusion = 0,
    scope_fr = "Cible de concentration résiduelle du consensus TDM linezolide.",
    scope_en = "Linezolid trough target from the TDM expert consensus.",
    grid_note_fr = "La cible est issue du consensus; la grille de doses est exploratoire et ne remplace pas l'adaptation à l'indication et à la tolérance.",
    grid_note_en = "The target comes from the consensus; the dose grid is exploratory and does not replace adjustment for indication and tolerability.",
    citation = "Lin B et al. Front Public Health. 2022;10:967311.",
    doi = "10.3389/fpubh.2022.967311"
  ),
  "voriconazole::ORAL" = list(
    target_metric = "Cmin",
    target_low = 1,
    target_high = 5.5,
    target_unit = "mg/L",
    dose_min = 50,
    dose_max = 600,
    dose_step = 50,
    intervals = c(8, 12, 24),
    infusion = 0,
    scope_fr = "Traitement d'une infection fongique invasive; cible résiduelle d'efficacité et de toxicité.",
    scope_en = "Treatment of invasive fungal infection; trough efficacy and toxicity target.",
    grid_note_fr = "La cible est issue de la recommandation; la grille de doses est exploratoire et doit tenir compte de l'indication, des interactions et de la fonction hépatique.",
    grid_note_en = "The target comes from the guideline; the dose grid is exploratory and must account for indication, interactions, and liver function.",
    citation = "Ashbee HR et al. J Antimicrob Chemother. 2014;69:1162-1176.",
    doi = "10.1093/jac/dkt508"
  )
)

clinical_target_preset <- function(record, mode) {
  drug_key <- as.character(record$drugKey[[1]] %||% "")
  CLINICAL_TARGET_PRESETS[[paste(drug_key, mode, sep = "::")]]
}
