# A/B dose comparison: deterministic model-based regression, no patient files.
suppressPackageStartupMessages({
  library(mrgsolve)
  library(mapbayr)
  library(ggplot2)
  library(jsonlite)
})
`%||%` <- function(value, fallback) if (is.null(value) || !length(value)) fallback else value
APP_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local = FALSE)
source(file.path(APP_ROOT, "R", "engine.R"), local = FALSE)
source(file.path(APP_ROOT, "R", "regimen_ab.R"), local = FALSE)

recommendations <- data.frame(
  dose = c(500, 600, 600, 750),
  interval = c(12, 12, 24, 336),
  infusion = c(1, 1, 1, 1)
)
minimal <- list(recommendations = recommendations)
candidates <- regimen_ab_candidates(minimal)
stopifnot(nrow(candidates) == 4L, identical(candidates$rank, seq_len(4L)),
          nrow(regimen_ab_candidates(list(recommendations = recommendations[1, , drop = FALSE]))) == 1L,
          nrow(regimen_ab_candidates(list(recommendations = rbind(recommendations[1, ], recommendations[1, ])))) == 1L,
          nrow(regimen_ab_candidates(list(recommendations = rbind(recommendations, data.frame(dose = 500, interval = 337, infusion = 1))))) == 4L,
          nrow(regimen_ab_candidates(NULL)) == 0L,
          is.null(regimen_ab_compare(minimal, 1L, 1L)),
          is.null(regimen_ab_compare(minimal, 1L, 99L)))

record <- model_record("vanco_roberts")
spec <- list(
  id = "vanco_roberts", label = "Roberts 2011", code = NULL,
  route = "IV", mode = model_administration_mode(record, "IV", "IV_CONTINUOUS"),
  adm_cmt_name = model_administration_cmt(record, "IV")
)
doses <- data.frame(time = 0, amount = 1000, interval = 12, count = 4L,
                    infusion = 12, ss = 0L, time_uncertainty = 0)
obs <- data.frame(time = 47.5, concentration = 20, time_uncertainty = 0)
fit <- fit_one_model(spec, doses, obs, covariates = list(CRCL = 90, WT = 74),
                     allow_custom = FALSE)
stopifnot(!is.null(fit$estimate))
result <- list(
  fits = list(vanco_roberts = fit),
  weights = c(vanco_roberts = 1),
  recommendations = recommendations,
  target_options = list(concentration_threshold = 5)
)
old <- result$recommendations
comparison <- regimen_ab_compare(result, 1L, 2L)
stopifnot(!is.null(comparison), nrow(comparison$metrics) == 2L,
          all(c("A", "B") %in% comparison$metrics$scenario),
          all(comparison$metrics$auc0_24 > 0),
          all(is.finite(comparison$metrics$cmin_interval)),
          all(is.finite(comparison$metrics$cmax_interval)),
          all(is.finite(comparison$profiles$concentration)),
          all(comparison$profiles$concentration >= 0),
          length(unique(comparison$profiles$scenario)) == 2L,
          isTRUE(all.equal(
            comparison$metrics$auc0_24[[2]] / comparison$metrics$auc0_24[[1]],
            600 / 500, tolerance = 0.025
          )),
          identical(result$recommendations, old))
app <- paste(readLines(file.path(APP_ROOT, "app.R"), encoding = "UTF-8", warn = FALSE), collapse = "\n")
stopifnot(grepl('source(file.path(APP_ROOT, "R", "regimen_ab.R")', app, fixed = TRUE),
          grepl('output$regimen_ab_plot <- renderPlot', app, fixed = TRUE),
          grepl('output$download_regimen_ab_profiles <- downloadHandler', app, fixed = TRUE),
          grepl('output$download_regimen_ab_metrics <- downloadHandler', app, fixed = TRUE))
cat("MIPD A/B: candidates, safe boundaries, independent overlays, linearity and UI hooks PASS\n")
