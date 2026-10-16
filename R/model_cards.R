# Compact bilingual model cards. Scientific fields remain available internally.
.model_source_review <- jsonlite::fromJSON(
  file.path(APP_ROOT, "models", "source-review.json")
)$reviewed
.model_fr_notes <- jsonlite::fromJSON(
  file.path(APP_ROOT, "models", "model-details.fr.json")
)$details

model_card_review <- function(id) {
  index <- match(id, .model_source_review$id)
  if (is.na(index)) stop("Missing source review for model: ", id)
  .model_source_review[index, , drop = FALSE]
}

model_card_status <- function(id) {
  review <- model_card_review(id)
  # This status expresses scientific SOURCE FIDELITY, not DOI or CI status.
  if (!isTRUE(review$fullTextReviewed[[1]]) ||
      identical(review$result[[1]], "SOURCE_AMBIGUITY_DOCUMENTED")) return("not verified")
  if (identical(review$result[[1]], "ARTICLE_ALIGNED")) return("aligned")
  if (identical(review$result[[1]], "DOCUMENTED_ADAPTATION_CONFIRMED")) return("adapted")
  "not verified"
}

model_card_detail <- function(id, english = FALSE) {
  review <- model_card_review(id)
  if (english) return(as.character(review$note[[1]] %||% ""))
  as.character(.model_fr_notes[[id]] %||% review$note[[1]] %||% "")
}

model_catalog_card <- function(record, english = FALSE) {
  id <- as.character(record$id[[1]])
  drug <- as.character(if (english) record$drugEn[[1]] else record$drug[[1]])
  author <- as.character(record$model[[1]])
  population <- as.character(if (english) record$populationEn[[1]] else record$population[[1]])
  citation <- as.character(record$citation[[1]])
  doi <- as.character(record$doi[[1]])
  status <- model_card_status(id)
  detail <- model_card_detail(id, english)
  label <- function(fr, en) if (english) en else fr

  shiny::tags$article(
    class = "pk-model-card",
    shiny::div(
      class = "pk-model-card-heading",
      shiny::div(class = "pk-model-card-identity",
        shiny::strong(class = "pk-model-card-drug", drug),
        shiny::span(class = "pk-model-card-author", author)
      ),
      shiny::span(class = paste("pk-model-status", gsub(" ", "-", status, fixed = TRUE)),
                  status)
    ),
    shiny::p(class = "pk-model-card-population", population),
    shiny::p(class = "pk-model-card-citation",
      shiny::span(class = "pk-model-source-label", label("Source : ", "Source: ")),
      shiny::span(citation),
      " · ",
      shiny::a(paste0("DOI ", doi), href = paste0("https://doi.org/", doi),
               target = "_blank", rel = "noopener noreferrer")
    ),
    if (status != "aligned") shiny::tags$details(
      class = "pk-model-card-details",
      shiny::tags$summary(label("Détails", "Details")),
      shiny::p(detail)
    )
  )
}

model_catalog_cards <- function(catalog, query = "", english = FALSE) {
  query <- trimws(as.character(query %||% "")[[1]])
  if (nzchar(query)) {
    keep <- vapply(seq_len(nrow(catalog)), function(i) {
      item <- catalog[i, , drop = FALSE]
      fields <- c(item$drug[[1]], item$drugEn[[1]], item$model[[1]],
                  item$population[[1]], item$populationEn[[1]], item$citation[[1]],
                  item$doi[[1]], item$id[[1]])
      grepl(tolower(query), tolower(paste(fields, collapse = " ")), fixed = TRUE)
    }, logical(1))
    catalog <- catalog[keep, , drop = FALSE]
  }
  if (!nrow(catalog)) return(shiny::p(
    class = "pk-model-empty",
    if (english) "No model matches your search." else "Aucun modèle ne correspond à cette recherche."
  ))
  shiny::div(
    class = "pk-model-grid",
    lapply(seq_len(nrow(catalog)), function(i) model_catalog_card(catalog[i, , drop = FALSE], english))
  )
}
