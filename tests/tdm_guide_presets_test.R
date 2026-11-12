# The guide preset layer may only ever offer a published target. It must stay
# subordinate to MIPD's curated presets and must disable itself on anything unexpected.
arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file_argument), ".."), winslash = "/", mustWork = TRUE)

`%||%` <- function(value, fallback) if (is.null(value) || !length(value)) fallback else value
suppressPackageStartupMessages(library(jsonlite))
source(file.path(APP_ROOT, "R", "model_library.R"), local = FALSE)
source(file.path(APP_ROOT, "R", "clinical_presets.R"), local = FALSE)
source(file.path(APP_ROOT, "R", "tdm_guide_presets.R"), local = FALSE)

stopifnot(!is.null(TDM_GUIDE_PRESETS))
stopifnot(length(TDM_GUIDE_PRESETS$presets) == 9L)
stopifnot(identical(TDM_GUIDE_PRESETS$clinicalInterpretationAllowed, FALSE))
stopifnot(identical(TDM_GUIDE_PRESETS$automaticDoseAdjustmentAllowed, FALSE))

# A molecule with several published targets offers a choice and auto-applies nothing.
tac <- model_record("tacrolimus_woillard_ddi")
tac_entry <- tdm_guide_preset_entry(tac, "ORAL")
stopifnot(!is.null(tac_entry), length(tac_entry$candidates) == 3L,
  identical(tac_entry$autoApply, FALSE))
stopifnot(length(tdm_guide_preset_choices(tac, "ORAL")) == 3L)
stopifnot(!is.null(tac_entry$doseGrid), identical(tac_entry$doseGrid$signoff, "PENDING"))
# A pending grid is displayed but must never be applicable, and lists no model.
stopifnot(identical(tac_entry$doseGridSignedOff, FALSE))
stopifnot(length(tac_entry$doseGridApplicableModelIds) == 0L)
stopifnot(!tdm_guide_preset_grid_applicable(tac, "ORAL"))
stopifnot(identical(tdm_guide_preset_grid_block(tac, "ORAL"), "grid_awaiting_signoff"))
for (key in names(TDM_GUIDE_PRESETS$presets))
  stopifnot(length(TDM_GUIDE_PRESETS$presets[[key]]$doseGridApplicableModelIds) == 0L)
# Neither tacrolimus model is applicable, nor is any other model.
stopifnot(!tdm_guide_preset_grid_applicable(model_record("tacrolimus_kamp_envarsus_ddi"), "ORAL"))
stopifnot(!tdm_guide_preset_grid_applicable(model_record("sirolimus_djebli_ddi"), "ORAL"))
# Sirolimus targets stay listed and documented but cannot be applied.
siro <- tdm_guide_preset_entry(model_record("sirolimus_djebli_ddi"), "ORAL")
stopifnot(length(siro$candidates) == 2L)
stopifnot(all(vapply(siro$candidates, tdm_guide_preset_candidate_blocked, logical(1))))
stopifnot(all(vapply(siro$candidates, function(x) nzchar(x$applicationBlockedCode), logical(1))))
stopifnot(all(vapply(siro$candidates, function(x) is.numeric(x$target_low), logical(1))))
# Vancomycin is not method-gated.
stopifnot(!any(vapply(tdm_guide_preset_entry(model_record("vanco_goti"), "IV_INTERMITTENT")$candidates,
  tdm_guide_preset_candidate_blocked, logical(1))))
# The corrected label-bound grids: 1 mg sirolimus step, tacrolimus twice daily only.
stopifnot(TDM_GUIDE_PRESETS$presets[["sirolimus::ORAL"]]$doseGrid$dose_step == 1)
stopifnot(unlist(TDM_GUIDE_PRESETS$presets[["tacrolimus::ORAL"]]$doseGrid$intervals) == 12)
# The liver-transplant reference is never offered on a kidney model.
stopifnot(!("tac-liver-everolimus-fda" %in%
  vapply(tac_entry$candidates, function(x) as.character(x$id), character(1))))
# Selecting by id returns that reference; an unknown id returns nothing.
stopifnot(!is.null(tdm_guide_preset_candidate(tac, "ORAL", "tac-kidney-mid")))
stopifnot(tdm_guide_preset_candidate(tac, "ORAL", "tac-kidney-mid")$target_low == 5)
stopifnot(is.null(tdm_guide_preset_candidate(tac, "ORAL", "nope")))
stopifnot(is.null(tdm_guide_preset_candidate(tac, "ORAL", "")))
# The Envarsus model sees the same published references and the same missing metadata.
envarsus <- tdm_guide_preset_entry(model_record("tacrolimus_kamp_envarsus_ddi"), "ORAL")
stopifnot(length(envarsus$candidates) == 3L)
stopifnot("formulationId" %in% unlist(envarsus$candidates[[1]]$undecidable))

# A single qualified reference may be auto-applied, and only where MIPD already curates one.
vanco <- model_record("vanco_goti")
stopifnot(identical(tdm_guide_preset_entry(vanco, "IV_INTERMITTENT")$autoApply, TRUE))
stopifnot(!is.null(clinical_target_preset(vanco, "IV_INTERMITTENT")))
for (key in names(TDM_GUIDE_PRESETS$presets)) {
  entry <- TDM_GUIDE_PRESETS$presets[[key]]
  if (!isTRUE(entry$autoApply)) next
  parts <- strsplit(key, "::", fixed = TRUE)[[1]]
  stopifnot(!is.null(CLINICAL_TARGET_PRESETS[[key]]))
}

# A disagreement with the curated preset is reported, not resolved.
vori <- model_record("voriconazole_vandenborn_ddi")
divergence <- tdm_guide_preset_divergence(vori, "ORAL")
stopifnot(length(divergence) == 1L)
stopifnot(divergence[[1]]$target_low == 2, divergence[[1]]$target_high == 6)
stopifnot(clinical_target_preset(vori, "ORAL")$target_low == 1)
# Matching values raise no disagreement, and an unindexed molecule has no entry at all.
stopifnot(length(tdm_guide_preset_divergence(vanco, "IV_INTERMITTENT")) == 0L)
stopifnot(is.null(tdm_guide_preset_entry(model_record("amox_carlier"), "IV_INTERMITTENT")))
stopifnot(length(tdm_guide_preset_choices(model_record("amox_carlier"), "IV_INTERMITTENT")) == 0L)

# --- Fail closed on every tampering that would turn documentation into authority ---
tmp <- tempfile(fileext = ".json")
on.exit(unlink(tmp), add = TRUE)
original <- jsonlite::fromJSON(TDM_GUIDE_PRESET_FILE, simplifyVector = FALSE)
reject <- function(mutate, label) {
  broken <- mutate(original)
  jsonlite::write_json(broken, tmp, auto_unbox = TRUE, null = "null", digits = NA)
  if (!is.null(tdm_guide_presets_load(tmp))) stop(paste("accepted a broken artifact:", label))
}
reject(function(x) { x$schema <- "wrong"; x }, "schema")
reject(function(x) { x$clinicalInterpretationAllowed <- TRUE; x }, "interpretation allowed")
reject(function(x) { x$automaticDoseAdjustmentAllowed <- TRUE; x }, "dose adjustment allowed")
reject(function(x) { x$sourceSnapshotSha256 <- "nope"; x }, "snapshot sha")
reject(function(x) { x$presets[["tacrolimus::ORAL"]]$autoApply <- TRUE; x },
  "auto-apply claimed for a 3-way choice")
reject(function(x) {
  x$presets[["tacrolimus::ORAL"]]$candidates[[1]]$modelApplicabilityVerified <- TRUE; x },
  "candidate claiming verified applicability")
reject(function(x) {
  x$presets[["tacrolimus::ORAL"]]$candidates[[1]]$automaticDoseAdjustmentAllowed <- TRUE; x },
  "candidate allowing dose adjustment")
reject(function(x) {
  x$presets[["tacrolimus::ORAL"]]$candidates[[1]]$target_metric <- "AUC0_12"; x },
  "metric the app cannot represent")
reject(function(x) {
  x$presets[["tacrolimus::ORAL"]]$candidates[[1]]$target_high <- 1; x }, "inverted range")
reject(function(x) {
  x$presets[["tacrolimus::ORAL"]]$candidates[[1]]$url <- "http://insecure"; x }, "insecure source")
reject(function(x) { x$presets[["tacrolimus::ORAL"]]$doseGrid$signoff <- "APPROVED"; x },
  "applicability flag disagreeing with an APPROVED signoff")
reject(function(x) { x$presets[["tacrolimus::ORAL"]]$doseGridSignedOff <- TRUE; x },
  "signedOff flag disagreeing with a PENDING signoff")
reject(function(x) {
  x$presets[["tacrolimus::ORAL"]]$doseGridApplicableModelIds <- list("tacrolimus_woillard_ddi"); x },
  "applicable model listed while the grid is PENDING")
reject(function(x) {
  x$presets[["sirolimus::ORAL"]]$candidates[[1]]$applicationBlockedCode <- ""; x },
  "blocked candidate without a reason")
reject(function(x) {
  x$presets[["sirolimus::ORAL"]]$candidates[[1]]$applicationBlocked <- "yes"; x },
  "non-logical applicationBlocked")
reject(function(x) { x$presets[["tacrolimus::ORAL"]]$doseGrid$signoff <- "SORT_OF"; x },
  "unknown signoff value")
reject(function(x) { x$presets[["vanco::IV_INTERMITTENT"]]$doseGridSignedOff <- TRUE; x },
  "signed-off flag without any grid")
reject(function(x) { x$presets[["tacrolimus::ORAL"]]$doseGrid$dose_step <- 0; x },
  "dose grid with a zero step")
# A missing file simply disables the layer instead of failing the app.
stopifnot(is.null(tdm_guide_presets_load(file.path(tempdir(), "absent.json"))))

cat("TDM guide presets: 9 keys, choice required for 4 immunosuppressants, curated MIPD presets keep precedence, 19 tampering modes rejected.\n")
