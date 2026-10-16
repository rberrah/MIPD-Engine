arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file_argument), ".."), winslash = "/", mustWork = TRUE)
setwd(APP_ROOT)

environment <- new.env(parent = globalenv())
sys.source(file.path(APP_ROOT, "app.R"), envir = environment)

shiny::testServer(environment$server, {
  session$setInputs(
    ddi_model_1 = "tacrolimus_woillard_ddi",
    ddi_route_1 = "Oral",
    ddi_dose_1 = 3,
    ddi_interval_1 = 12,
    ddi_infusion_1 = 0,
    ddi_target_parameter = "TVCL_TAC",
    ddi_interaction_type = "factor",
    ddi_factor = 0.5,
    ddi_model_2 = "voriconazole_vandenborn_ddi",
    ddi_route_2 = "Oral",
    ddi_dose_2 = 200,
    ddi_interval_2 = 12,
    ddi_infusion_2 = 0,
    ddi_start_day = 3,
    ddi_stop_day = 10,
    ddi_followup = 3,
    ddi_accept_disclaimer = TRUE,
    run_ddi = 1
  )
  session$flushReact()
  result <- ddi_store()
  stopifnot(
    !is.null(result),
    result$affected$id == "tacrolimus_woillard_ddi",
    result$driver$id == "voriconazole_vandenborn_ddi",
    result$config$target == "TVCL_TAC",
    result$metrics$auc_ratio > 1
  )

  session$setInputs(ddi_add_target = 1)
  session$flushReact()
  session$setInputs(
    ddi_target_fraction = 0.5,
    ddi_target_parameter_2 = "TVV1_TAC",
    ddi_target_fraction_2 = 0.25,
    run_ddi = 2
  )
  session$flushReact()
  result <- ddi_store()
  stopifnot(
    length(result$config$targets) == 2,
    identical(vapply(result$config$targets, `[[`, character(1), "parameter"), c("TVCL_TAC", "TVV1_TAC")),
    identical(vapply(result$config$targets, `[[`, numeric(1), "fraction"), c(0.5, 0.25)),
    length(ggplot2::ggplot_build(ddi_effect_figure(result))$data) >= 3
  )

  fit_1 <- list(
    id = result$affected$id,
    label = result$affected$label,
    model = result$affected$model,
    estimate = NULL,
    contract = list(route = result$affected$route, adm_cmt = result$affected$adm_cmt),
    split_lego = result$affected$split_lego,
    current_covariates = list(HT = 40),
    source_doses = data.frame(time = 0, amount = 4, interval = 12, infusion = 0)
  )
  analysis_store(list(fits = list(fit_1)))
  session$setInputs(ddi_use_tdm_1 = 1)
  session$flushReact()
  stopifnot(ddi_subject_1()$source == "tdm", ddi_subject_1()$id == result$affected$id)

  fit_2 <- list(
    id = result$driver$id,
    label = result$driver$label,
    model = result$driver$model,
    estimate = NULL,
    contract = list(route = result$driver$route, adm_cmt = result$driver$adm_cmt),
    split_lego = result$driver$split_lego,
    current_covariates = list(CRP = 75),
    source_doses = data.frame(time = 0, amount = 250, interval = 12, infusion = 0)
  )
  analysis_store(list(fits = list(fit_2)))
  session$setInputs(ddi_use_tdm_2 = 1)
  session$flushReact()
  stopifnot(
    ddi_subject_1()$source == "tdm",
    ddi_subject_2()$source == "tdm",
    ddi_subject_2()$id == result$driver$id
  )
})

cat("Generic DDI Shiny server test passed.\n")
