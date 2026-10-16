pd_builtin_iv_context <- function(volume, clearance, soloc, cache) {
  volume <- ddi_numeric(volume, "V", 1e-6, 1e6)
  clearance <- ddi_numeric(clearance, "CL", 1e-6, 1e6)
  # Compile one trusted Lego structure per session; changing CL/V only updates parameters.
  specification <- list(version = 3, nodes = list(list(id = 1, kind = "central", name = "CENT", vol = 20, dose = 100)),
    edges = list(list(from = 1, to = "OUT", kinetics = "first_order", eliminationParameterization = "clearance", cl = 4, k = .2)), covariates = list())
  code <- safe_lego_model_code(specification = specification)
  model <- compile_model(custom_code = code, allow_custom = FALSE, custom_soloc = soloc, custom_cache = cache)
  model <- safe_param(model, list(TV_v_L1_CENT = volume, TV_cl_L1_CENT = clearance))
  list(id = "custom-pk", label = "IV / 1 CMT", model = mrgsolve::zero_re(model), code = code,
    adm_cmt = 1L, route = "IV", source = "code", concentration = "CONC_L1_CENT", time_unit = "h", concentration_scale = 1,
    covariate_names = character())
}

pd_pk_units <- function(id, scale = 1, time_unit = "h") {
  if (id %in% MODEL_CATALOG$id) {
    record <- model_record(id)
    unit <- record$concentrationUnit[[1]]
    return(list(time_unit = record$timeUnit[[1]], concentration_scale = if (unit == "ng/mL") .001 else 1, unit = unit))
  }
  scale <- ddi_numeric(scale, "PK concentration unit", 1e-9, 1e9)
  list(time_unit = time_unit, concentration_scale = scale,
    unit = switch(as.character(scale), "1" = "mg/L", "0.001" = "ng/mL", "1000" = "mg/mL", "custom"))
}

pd_pk_tdm_context <- function(fit, scale = 1) {
  has_individual_fit <- !is.null(fit$estimate) ||
    (identical(fit$engine %||% "mapbayr", "pmetrics") && isTRUE(fit$pmetrics$mapped))
  if (!has_individual_fit) stop("A Bayesian PK fit is required. Run the TDM analysis first.")
  # PD uses the MAP estimate, never the experimental ML override.
  fit$ml_eta_override <- NULL
  value <- ddi_fit_context(fit)
  captures <- model_capture_names(value$model)
  preferred <- intersect(c("CP", "CONC", "IPRED", "CONC_PLASMA", "CONCENTRATION", "DV"), captures)
  if (!length(preferred)) stop("No PK concentration output found.")
  value$concentration <- preferred[[1]]
  utils::modifyList(value, pd_pk_units(value$id, scale))
}

pd_pk_ui <- function(id, t, allow_tdm = TRUE) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("source"), t("Source PK", "PK source"), c("MIPD" = "library", "mrgsolve / Lego" = "code", if (allow_tdm) c("TDM / MAP-BE" = "tdm")), inline = TRUE),
    conditionalPanel(sprintf("input['%s']=='library'", ns("source")),
      selectInput(ns("model"), t("Modele PK", "PK model"), catalog_choices_i18n(lang = t("fr", "en")), "vanco_pkjust"), uiOutput(ns("route_ui"))),
    conditionalPanel(sprintf("input['%s']=='code'", ns("source")),
      textAreaInput(ns("code"), "mrgsolve / C++", rows = 8),
      selectInput(ns("custom_route"), t("Voie", "Route"), c("IV", "Oral"))),
    conditionalPanel(sprintf("input['%s']!='tdm'", ns("source")),
      actionButton(ns("load"), t("Charger / compiler la PK", "Load / compile PK"), icon = icon("gears")),
      actionButton(ns("estimate"), t("Estimer la PK avec le TDM", "Estimate PK with TDM"), icon = icon("chart-line"))),
    conditionalPanel(sprintf("input['%s']=='tdm'", ns("source")), uiOutput(ns("tdm_models"))),
    uiOutput(ns("status")),
    conditionalPanel(sprintf("input['%s']=='code'", ns("source")),
      selectInput(ns("time_unit"), t("Unite de temps du code PK", "PK code time unit"), c("h" = "h", "day" = "day"))),
    uiOutput(ns("units")),
    uiOutput(ns("outputs")),
    uiOutput(ns("covariates")),
    conditionalPanel(sprintf("input['%s']!='tdm'", ns("source")),
      tags$details(tags$summary(t("Parametres populationnels avances", "Advanced population parameters")), uiOutput(ns("parameters"))))
  )
}

pd_pk_server <- function(id, soloc, cache, imported = reactive(NULL), analysis_store = reactive(NULL), open_tdm = NULL, allow_tdm = TRUE) {
  moduleServer(id, function(input, output, session) {
    t <- function(fr, en) app_t(app_language_from_query(session$clientData$url_search %||% ""), fr, en)
    loaded <- reactiveVal(NULL)
    load_error <- reactiveVal(NULL)
    generation <- reactiveVal(0L)
    pending <- reactiveVal(NULL)
    auto_import_key <- ""
    output$tdm_models <- renderUI({
      fits <- successful_fits(analysis_store()$fits %||% list())
      if (!length(fits)) return(p(t("Aucune estimation. Lancez une analyse TDM dans Analyse.", "No estimate. Run a TDM analysis in Analysis.")))
      selectInput(session$ns("tdm_model"), t("Estimation bayesienne", "Bayesian estimate"), setNames(names(fits), vapply(fits, `[[`, character(1), "label")))
    })
    observeEvent(input$estimate, tryCatch({
      if (is.null(open_tdm)) stop("TDM navigation is unavailable.")
      if (identical(input$source, "code") && identical(input$time_unit, "day")) stop(t("L'analyse TDM utilise des heures. Convertissez d'abord les constantes de temps du code PK en heures.", "TDM analysis uses hours. Convert PK code time constants to hours first."))
      open_tdm(list(source = input$source, id = input$model, route = if (input$source == "code") input$custom_route else input$route, code = input$code))
      if (allow_tdm) updateRadioButtons(session, "source", selected = "tdm")
    }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = 12)))
    output$route_ui <- renderUI({
      model <- input$model %||% DEFAULT_MODEL
      routes <- model_routes(model_record(model))
      selected <- if (identical(pending()$id, model)) pending()$route else isolate(input$route)
      if (is.null(selected) || !selected %in% routes) selected <- routes[[1]]
      selectInput(session$ns("route"), t("Voie", "Route"), routes, selected)
    })
    load_pk <- function(source, code = NULL, model_id = NULL, route = NULL) {
      if (!identical(source, "code")) return(ddi_library_context(model_id, route))
      if (!is.character(code) || !nzchar(trimws(code)) || nchar(code, type = "bytes") > 200000) stop("Supply a complete mrgsolve model (maximum 200 kB).")
      model <- compile_model(custom_code = code, allow_custom = ALLOW_CUSTOM_MODELS, custom_soloc = soloc, custom_cache = cache)
      model <- mrgsolve::zero_re(model)
      if (!length(model@cmtL) || !length(model_capture_names(model))) stop("PK requires compartments and a captured concentration.")
      adm <- tagged_compartment(model, "ADM")
      list(id = "custom-pk", label = "mrgsolve / Lego", model = model, code = code,
        adm_cmt = if (is.finite(adm)) adm else 1L, route = route %||% "IV", source = "code", covariate_names = parse_covariates(code)$name)
    }
    store_pk <- function(value, clear_pending = TRUE) {
      next_generation <- isolate(generation()) + 1L
      generation(next_generation)
      value$generation <- next_generation
      loaded(value)
      if (clear_pending) pending(NULL)
    }
    observeEvent(imported(), {
      spec <- imported(); if (is.null(spec)) return()
      key <- paste(spec$source %||% "", spec$id %||% "", spec$route %||% "", spec$code %||% "", sep = "\n")
      if (identical(auto_import_key, key)) return()
      auto_import_key <<- key
      pending(spec); loaded(NULL)
      load_error(NULL)
      updateRadioButtons(session, "source", selected = spec$source)
      if (spec$source == "code") {
        updateTextAreaInput(session, "code", value = spec$code)
        updateSelectInput(session, "custom_route", selected = spec$route %||% "IV")
      } else updateSelectInput(session, "model", selected = spec$id)
      updateSelectInput(session, "time_unit", selected = spec$time_unit %||% "h")
      updateSelectInput(session, "scale", selected = as.character(spec$concentration_scale %||% 1))
      loading_message <- t("Chargement PK", "Loading PK")
      session$onFlushed(function() tryCatch({
        value <- withProgress(message = loading_message, value = 0.3,
          load_pk(spec$source, spec$code, spec$id, spec$route))
        store_pk(value, clear_pending = FALSE)
        session$onFlushed(function() pending(NULL), once = TRUE)
      }, error = function(e) {
        auto_import_key <<- ""
        load_error(conditionMessage(e))
        showNotification(conditionMessage(e), type = "error", duration = 10)
      }), once = TRUE)
    })
    observeEvent(input$load, tryCatch({
      load_error(NULL)
      value <- withProgress(message = t("Chargement PK", "Loading PK"), value = 0.3,
        load_pk(input$source, input$code, input$model, if (identical(input$source, "code")) input$custom_route else input$route))
      store_pk(value)
    }, error = function(e) {
      load_error(conditionMessage(e))
      showNotification(conditionMessage(e), type = "error", duration = 10)
    }))
    valid_loaded <- reactive({
      if (identical(input$source, "tdm")) return(NULL)
      value <- loaded()
      if (is.null(value)) return(NULL)
      if (identical(input$source, "code")) {
        normalize_lines <- function(code) gsub("\r\n?", "\n", code %||% "")
        if (!identical(normalize_lines(value$code), normalize_lines(input$code)) || value$id != "custom-pk") return(NULL)
      } else if (!identical(value$id, input$model)) return(NULL)
      value
    })
    output$status <- renderUI({
      if (identical(input$source, "tdm")) return(p(t("PK issue de mapbayr; covariables de la derniere analyse. Modifiez les donnees dans Analyse puis relancez l'estimation.", "PK from mapbayr; covariates from the latest analysis. Edit data in Analysis and rerun estimation.")))
      if (!is.null(load_error())) return(p(class = "text-danger", load_error()))
      p(if (is.null(valid_loaded())) t("PK non chargee ou code modifie.", "PK not loaded or code modified.") else valid_loaded()$label)
    })
    output$units <- renderUI({
      id <- if (identical(input$source, "tdm")) input$tdm_model %||% "" else if (identical(input$source, "library")) input$model %||% "" else ""
      if (id %in% MODEL_CATALOG$id) return(p(paste(t("Unites PK :", "PK units:"), "h /", pd_pk_units(id)$unit)))
      selected <- as.character(pending()$concentration_scale %||% isolate(input$scale) %||% 1)
      choices <- c("mg/L (= ug/mL)" = "1", "ng/mL (= ug/L)" = "0.001", "mg/mL" = "1000")
      if (!selected %in% unname(choices)) choices <- c(choices, setNames(selected, t("Unite de l'atelier importe", "Imported workshop unit")))
      selectInput(session$ns("scale"), t("Unite des concentrations PK", "PK concentration unit"), choices, selected)
    })
    output$outputs <- renderUI({
      value <- valid_loaded(); shiny::req(value)
      captures <- model_capture_names(value$model)
      if (value$id %in% MODEL_CATALOG$id) return(NULL)
      preferred <- intersect(c("CP", "CONC", "IPRED", "CONC_PLASMA", "CONCENTRATION", "DV"), captures)
      tagList(selectInput(session$ns("concentration"), t("Sortie concentration", "Concentration output"), captures, if (length(preferred)) preferred[[1]] else captures[[1]]),
        if (identical(value$id, "custom-pk")) selectInput(session$ns("adm"), t("Compartiment d'administration", "Administration compartment"), setNames(seq_along(value$model@cmtL), value$model@cmtL), value$adm_cmt))
    })
    output$parameters <- renderUI({
      value <- valid_loaded(); shiny::req(value)
      p <- as.list(mrgsolve::param(value$model))
      names <- setdiff(names(p), c(value$covariate_names, grep("^ETA[0-9]+$", names(p), value = TRUE)))
      div(class = "pd-parameter-grid", lapply(names, function(name) numericInput(session$ns(paste0("param_", value$generation, "_", name)), name, p[[name]], step = NA)))
    })
    output$covariates <- renderUI({
      value <- valid_loaded(); shiny::req(value)
      p <- as.list(mrgsolve::param(value$model))
      covs <- parse_covariates(value$code)
      if (!nrow(covs)) return(NULL)
      tagList(h4(t("Covariables patient", "Patient covariates")), lapply(seq_len(nrow(covs)), function(i) {
        name <- covs$name[i]
        numericInput(session$ns(paste0("param_", value$generation, "_", name)), paste(name, covs$description[i], sep = " : "), p[[name]], step = NA)
      }))
    })
    context <- reactive({
      if (identical(input$source, "tdm")) {
        fit <- successful_fits(analysis_store()$fits %||% list())[[input$tdm_model %||% ""]]
        if (is.null(fit)) stop(t("Selectionnez une estimation TDM disponible.", "Select an available TDM estimate."))
        return(pd_pk_tdm_context(fit, input$scale %||% 1))
      }
      value <- valid_loaded()
      if (is.null(value)) stop(t("Chargez / compilez le modele PK selectionne.", "Load / compile the selected PK model."))
      if (identical(value$id, "custom-pk")) {
        value$route <- input$custom_route %||% "IV"
        value$adm_cmt <- as.integer(input$adm %||% value$adm_cmt)
      } else value <- ddi_context_with_route(value, input$route)
      if (!value$adm_cmt %in% seq_along(value$model@cmtL)) stop("Invalid administration compartment.")
      p <- as.list(mrgsolve::param(value$model))
      for (name in names(p)) p[[name]] <- ddi_numeric(input[[paste0("param_", value$generation, "_", name)]] %||% p[[name]], name)
      value$model <- mrgsolve::param(value$model, p)
      captures <- model_capture_names(value$model)
      preferred <- intersect(c("CP", "CONC", "IPRED", "CONC_PLASMA", "CONCENTRATION", "DV"), captures)
      value$concentration <- if (value$id %in% MODEL_CATALOG$id) preferred[[1]] else input$concentration %||% if (length(preferred)) preferred[[1]] else captures[[1]]
      if (!value$concentration %in% captures) stop("Select the PK concentration output.")
      utils::modifyList(value, pd_pk_units(value$id, input$scale %||% 1, input$time_unit %||% "h"))
    })
    session$onSessionEnded(function() { loaded(NULL); load_error(NULL); pending(NULL) })
    context
  })
}

# All input times use the caller's unit; infusion durations always use hours.
pd_pk_profile <- function(context, doses, horizon, hours_per_time = 1) {
  if (!context$time_unit %in% c("h", "day")) stop("PK time unit must be h or day.")
  if (context$route == "Oral" && any(doses$infusion > 0)) stop("Oral doses require infusion duration = 0.")
  if (nrow(doses) > 5000 || any(!is.finite(as.matrix(doses[c("time", "amount", "infusion")])))) stop("Invalid PK dose history.")
  scale_time <- hours_per_time / if (context$time_unit == "day") 24 else 1
  step <- 0.1 / hours_per_time
  knots <- c(doses$time, doses$time + doses$infusion / hours_per_time)
  grid <- sort(unique(c(seq(0, horizon, by = step), horizon, knots, knots - 1e-8 / hours_per_time)))
  grid <- grid[grid >= 0 & grid <= horizon]
  if (length(grid) > 200000) stop("PK preview grid limit exceeded.")
  obs <- data.frame(ID = 1, time = grid * scale_time, evid = 0, cmt = 0, amt = 0, rate = 0)
  events <- data.frame(ID = rep(1, nrow(doses)), time = doses$time * scale_time, evid = rep(1, nrow(doses)),
    cmt = rep(context$adm_cmt, nrow(doses)), amt = doses$amount,
    rate = ifelse(doses$infusion > 0, doses$amount / (doses$infusion / hours_per_time * scale_time), 0))
  data <- rbind(events, obs)
  data <- data[order(data$time, -data$evid), ]
  sim <- as.data.frame(mrgsolve::mrgsim_d(context$model, data = data, obsonly = TRUE, recsort = 3, nocb = FALSE))
  sim$time <- round(sim$time / scale_time, 10)
  sim <- sim[!duplicated(sim$time, fromLast = TRUE), ]
  c <- sim[[context$concentration]] * context$concentration_scale
  if (length(c) != nrow(sim) || any(!is.finite(c)) || any(c < -1e-7)) stop("Invalid PK concentration output.")
  data.frame(time = sim$time, concentration = pmax(0, c))
}

pd_pk_export_setup <- function(context) {
  if (is.null(context)) return("pk_context <- NULL")
  dump <- function(x) paste(capture.output(dput(x)), collapse = "\n")
  paste(c("# Requires mrgsolve and a compiler. Contains explicitly exported PK code and parameters.",
    "workdir <- tempfile('pd-pk-'); dir.create(workdir)",
    paste0("pk_context <- ", dump(context[intersect(names(context), c("label", "route", "adm_cmt", "concentration", "time_unit", "concentration_scale"))])),
    paste0("pk_context$model <- mrgsolve::mcode('pd_pk', ", dump(context$code), ", soloc=workdir, quiet=TRUE) |> mrgsolve::zero_re()"),
    paste0("pk_context$model <- mrgsolve::param(pk_context$model, ", dump(as.list(mrgsolve::param(context$model))), ")"),
    "# Unload the compiled model before deleting workdir after use."), collapse = "\n")
}
