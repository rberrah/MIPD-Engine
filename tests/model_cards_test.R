suppressPackageStartupMessages({
  library(shiny)
  library(jsonlite)
  library(htmltools)
})
APP_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
'%||%' <- function(x, fallback) if (is.null(x) || !length(x)) fallback else x
source(file.path(APP_ROOT, "R", "model_library.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "model_cards.R"), local = TRUE)

stopifnot(nrow(MODEL_CATALOG) == 46L)
stopifnot(length(unique(.model_source_review$id)) == 46L)
stopifnot(setequal(.model_source_review$id, MODEL_CATALOG$id))
stopifnot(all(vapply(MODEL_CATALOG$id, function(id) {
  status <- model_card_status(id)
  if (status != "aligned") nzchar(model_card_detail(id)) else TRUE
}, logical(1))))
all_status <- vapply(MODEL_CATALOG$id, model_card_status, character(1))
stopifnot(sum(all_status == "aligned") == 10L)
stopifnot(sum(all_status == "adapted") == 25L)
stopifnot(model_card_status("vanco_goti") == "not verified")
stopifnot(sum(all_status == "not verified") == 11L)

for (i in seq_len(nrow(MODEL_CATALOG))) {
  row <- MODEL_CATALOG[i, , drop = FALSE]
  for (english in c(FALSE, TRUE)) {
    html <- as.character(htmltools::renderTags(model_catalog_card(row, english))$html)
    stopifnot(grepl(paste0("https://doi.org/", row$doi), html, fixed = TRUE))
    stopifnot(grepl(row$model, html, fixed = TRUE))
    stopifnot(grepl(if (english) row$drugEn else row$drug, html, fixed = TRUE))
    stopifnot(grepl(if (english) "Source:" else "Source :", html, fixed = TRUE))
    status <- model_card_status(row$id[[1]])
    if (status == "aligned") {
      stopifnot(!grepl("<details", html, fixed = TRUE))
    } else {
      stopifnot(grepl("<details", html, fixed = TRUE))
      stopifnot(grepl(if (english) "Details" else "Détails", html, fixed = TRUE))
    }
    # Software-maturity metadata must not be exposed in the concise card.
    stopifnot(!grepl("SOFTWARE_TESTED", html, fixed = TRUE))
    stopifnot(!grepl("TEACHING_RESEARCH", html, fixed = TRUE))
  }
}
stopifnot(grepl("Mellon", as.character(htmltools::renderTags(
  model_catalog_cards(MODEL_CATALOG, query = "obèse")
)$html), fixed = TRUE))
stopifnot(!grepl("Mellon", as.character(htmltools::renderTags(
  model_catalog_cards(MODEL_CATALOG, query = "Djebli")
)$html), fixed = TRUE))
stopifnot(grepl("Aucun modèle", as.character(htmltools::renderTags(
  model_catalog_cards(MODEL_CATALOG, query = "un-findable-model-xyz")
)$html), fixed = TRUE))
cat("Model cards, bibliography links and independent statuses: PASS (46, FR/EN).\n")
