arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file_argument), ".."), winslash = "/", mustWork = TRUE)
setwd(APP_ROOT)

environment <- new.env(parent = globalenv())
sys.source(file.path(APP_ROOT, "app.R"), envir = environment)

fr <- htmltools::renderTags(environment$app_ui(list(QUERY_STRING = "lang=fr")))$html
en <- htmltools::renderTags(environment$app_ui(list(QUERY_STRING = "lang=en")))$html

stopifnot(
  grepl("pmetrics_artifact_id", fr, fixed = TRUE),
  grepl("Pmetrics", fr, fixed = TRUE),
  grepl("Pmetrics", en, fixed = TRUE)
)

stopifnot(
  grepl("Atelier DDI", fr, fixed = TRUE),
  grepl("DDI builder", en, fixed = TRUE),
  grepl("Molécule affectée", fr, fixed = TRUE),
  grepl("Affected drug", en, fixed = TRUE),
  grepl("Reprendre TDM / Lego", fr, fixed = TRUE),
  grepl("Use TDM / Lego fit", en, fixed = TRUE),
  grepl("Facteur constant", fr, fixed = TRUE),
  grepl("Emax inhibition", en, fixed = TRUE),
  grepl("Time-dependent inhibition", en, fixed = TRUE),
  grepl("Induction with turnover", en, fixed = TRUE),
  grepl("Import workshop (.json)", en, fixed = TRUE)
)
for (type in c("factor", "inhibition", "induction", "reversible", "hill_inhibition", "tdi", "turnover_induction")) {
  stopifnot(grepl(paste0('value="', type, '"'), fr, fixed = TRUE), grepl(paste0('value="', type, '"'), en, fixed = TRUE))
}

cat("Generic DDI Shiny UI test passed.\n")
