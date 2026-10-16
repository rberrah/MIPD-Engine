pd_type_choices <- function(t) stats::setNames(PD_TYPES, c(
  t("Lineaire", "Linear"), "Emax", "Emax / Hill",
  t("Inhibition de la production", "Inhibit production"), t("Stimulation de la production", "Stimulate production"),
  t("Inhibition de la degradation", "Inhibit loss"), t("Stimulation de la degradation", "Stimulate loss")))

pd_general_panel <- function(lang = "fr") {
  t <- function(fr, en) app_t(lang, fr, en)
  ns <- shiny::NS("pd")
  nav_panel(t("PD generale", "General PD"), value = "general",
    div(class = "analysis-shell pd-shell",
      fileInput("pd_workshop_file", t("Importer un atelier (.json)", "Import workshop (.json)"), accept = ".json"),
      div(class = "safety-banner", tags$strong(t("Recherche et enseignement uniquement", "Research and teaching only")),
        span(t("Ajustement individuel par moindres carres; aucune recommandation clinique.", "Individual least-squares fitting; no clinical recommendation."))),
      div(class = "model-chain",
        div(class = "model-block exposure", icon("pills"), h3(t("Exposition", "Exposure")),
          selectInput(ns("exposure"), NULL, stats::setNames(c("exponential", "table", "pk"), c(t("Profil exponentiel", "Exponential profile"), t("Concentrations connues", "Known concentrations"), t("Modele PK libre", "Free PK model"))))),
        span(class = "model-arrow", icon("arrow-right")),
        div(class = "model-block delay", icon("clock"), h3(t("Delai", "Delay")), checkboxInput(ns("delay"), t("Compartiment d'effet", "Effect compartment"), FALSE)),
        span(class = "model-arrow", icon("arrow-right")),
        div(class = "model-block response", icon("chart-line"), h3(t("Reponse", "Response")), selectInput(ns("type"), NULL, pd_type_choices(t), "emax"))),
      tags$button(type = "button", class = "mobile-configure-button", onclick = "var r=this.closest('.pd-shell');var b=r.querySelector('.collapse-toggle');if(b&&b.getAttribute('aria-expanded')!=='true')b.click();", t("Parametres et donnees", "Parameters and data")),
      layout_sidebar(sidebar = sidebar(width = 360, open = "desktop",
        conditionalPanel(sprintf("input['%s']=='pk'", ns("exposure")), pd_pk_ui(ns("pk"), t),
          numericInput(ns("dose"), t("Dose (unite PK)", "Dose (PK unit)"), 100, min = 0.001),
          numericInput(ns("interval"), t("Intervalle (h)", "Interval (h)"), 24, min = 0.25),
          numericInput(ns("infusion"), t("Perfusion h (0 = oral / bolus)", "Infusion h (0 = oral / bolus)"), 0, min = 0)),
        h3(t("Parametres", "Parameters")), uiOutput(ns("parameters")),
        numericInput(ns("horizon"), t("Duree (h)", "Duration (h)"), 24, min = 0.01, max = 2400),
        conditionalPanel(sprintf("input['%s']=='exponential'", ns("exposure")),
          fluidRow(column(6, numericInput(ns("c0"), "C(0)", 10, min = 0)), column(6, numericInput(ns("kel"), "kel (1/h)", 0.1, min = 0)))),
        conditionalPanel(sprintf("input['%s']=='table'", ns("exposure")),
          textAreaInput(ns("profile"), t("Profil CSV : time,concentration", "CSV profile: time,concentration"), "time,concentration\n0,10\n6,5.49\n12,3.01\n24,0.91", rows = 6),
          actionButton(ns("use_tdm"), t("Reprendre l'exposition TDM", "Use TDM exposure"), icon = icon("file-import"))),
        actionButton(ns("simulate"), t("Simuler", "Simulate"), icon = icon("play"), class = "btn-primary"),
        tags$hr(), h3(t("Observations", "Observations")),
        radioButtons(ns("data_mode"), t("Type de donnees", "Data type"), stats::setNames(c("time", "concentration"), c(t("Temps-effet", "Time-effect"), t("Concentration-effet", "Concentration-effect")))),
        fileInput(ns("file"), t("Importer CSV", "Import CSV"), accept = ".csv"),
        pd_observations_ui(ns("observations"), t),
        actionButton(ns("example"), t("Donnees synthetiques", "Synthetic data"), icon = icon("flask")),
        uiOutput(ns("estimate_parameters")),
        checkboxInput(ns("accept"), t("Usage exploratoire, sous ma responsabilite", "Exploratory use, under my responsibility"), FALSE),
        actionButton(ns("fit"), t("Ajuster les parametres", "Fit parameters"), icon = icon("sliders"), class = "btn-primary")),
        div(class = "workspace pd-workspace", uiOutput(ns("status")),
          navset_tab(id = ns("tabs"),
            nav_panel(t("Simulation", "Simulation"), value = "simulation",
              plotOutput(ns("effect_plot"), height = "360px"), plotOutput(ns("hysteresis_plot"), height = "300px"),
              downloadButton(ns("download_profile"), t("Exporter les courbes", "Export curves"))),
            nav_panel(t("Ajustement", "Fit"), value = "fit",
              uiOutput(ns("fit_summary")), DTOutput(ns("fit_parameters")), plotOutput(ns("fit_plot"), height = "350px"), DTOutput(ns("fit_data")),
              downloadButton(ns("download_report"), t("Creer le rapport", "Create report"))),
            nav_panel(t("Modele", "Model"), value = "model", uiOutput(ns("equation")),
              tags$h4("mrgsolve / C++"), verbatimTextOutput(ns("code")), downloadButton(ns("download_code"), t("Exporter .cpp", "Export .cpp")),
              tags$h4("R / deSolve"), downloadButton(ns("download_r"), t("Exporter le script R", "Export R script"))),
            nav_panel(t("Methode et limites", "Methods and limitations"),
              tags$h3(t("Perimetre", "Scope")),
              tags$p(t("Les effets directs peuvent etre ajustes sur des couples concentration-effet. Un delai ou une reponse indirecte exige des temps et un historique d'exposition couvrant toute la simulation.", "Direct effects support concentration-effect pairs. A delay or indirect response requires times and an exposure history covering the simulation.")),
              tags$p(t("Les etats initiaux sont Ce(0)=0 et R(0)=E0. Pour les reponses indirectes, kin=E0*kout. L'inhibition est bornee entre 0 et 1. Les unites de concentration, EC50 et des observations doivent etre coherentes.", "Initial states are Ce(0)=0 and R(0)=E0. For indirect responses, kin=E0*kout. Inhibition is bounded between 0 and 1. Concentration, EC50 and observation units must be consistent.")),
              tags$p(t("La PK provient d'un modele populationnel avec covariables renseignees, ou de l'estimation bayesienne mapbayr realisee dans Analyse. Les unites de la bibliotheque sont connues; un code externe demande ses unites. Les concentrations sont exprimees en mg/L. La PK est ensuite maintenue fixe pendant l'ajustement PD.", "PK comes from a population model with entered covariates, or from the mapbayr Bayesian estimate in Analysis. Library units are known; external code requires its units. Concentrations use mg/L. PK is then held fixed during PD fitting.")),
              tags$p(t("Integration LSODA (deSolve). Ajustement PD par Levenberg-Marquardt borne, erreur additive : ce n'est pas une estimation bayesienne PD. L'incertitude PK n'est pas propagee. La convergence n'etablit pas l'identifiabilite ou la validite clinique.", "LSODA integration (deSolve). PD fitting uses bounded Levenberg-Marquardt with additive error, not Bayesian PD estimation. PK uncertainty is not propagated. Convergence does not establish identifiability or clinical validity.")),
              tags$p(t("Le fichier C++ utilise CP comme entree externe. Le script R exporte reproduit l'exposition et les equations utilisees ici. En PK libre, il contient le profil calcule, fige pour les conditions exportees, et ne recalcule pas la PK si vous modifiez la dose. Les donnees restent en session; seuls vos telechargements les conservent.", "The C++ file uses CP as an external input. The exported R script reproduces the exposure and equations used here. With free PK it contains the calculated profile, fixed for the exported conditions, and does not recalculate PK if you change the dose. Data remain in session; only your downloads retain them.")),
              tags$a(href = "https://desolve.r-forge.r-project.org/main.html", target = "_blank", rel = "noopener noreferrer", "deSolve"), " / ",
              tags$a(href = "https://mrgsolve.org/user-guide/specification.html", target = "_blank", rel = "noopener noreferrer", "mrgsolve")))))))
}

pd_server <- function(id, analysis_store, report_plot_uri, imported = reactive(NULL), soloc = tempdir(), cache = new.env(), open_tdm = NULL) {
  shiny::moduleServer(id, function(input, output, session) {
    t <- function(fr, en) app_t(app_language_from_query(session$clientData$url_search %||% ""), fr, en)
    fit_store <- reactiveVal(NULL)
    imported_parameters <- reactiveVal(NULL)
    pk <- pd_pk_server("pk", soloc, cache, reactive(if (identical(imported()$view, "pd")) imported()$models[[1]] else NULL), analysis_store, open_tdm)
    observations <- pd_observations_server("observations", reactive(c(if (identical(input$data_mode, "concentration")) "concentration" else "time", "effect")),
      data.frame(time = c(0, 2, 4, 8, 12, 18, 24), effect = c(140, 139, 136, 132, 127, 122, 115)))
    pk_profile <- reactive({
      shiny::req(identical(input$exposure, "pk"))
      regimen <- ddi_validate_regimen(list(dose = input$dose %||% 100, interval = input$interval %||% 24, infusion = input$infusion %||% 0), "PK")
      horizon <- input$horizon %||% 24
      time <- seq(0, horizon - 1e-8, by = regimen$interval)
      pd_pk_profile(pk(), data.frame(time = time, amount = regimen$dose, infusion = regimen$infusion), horizon)
    })
    names_active <- reactive(pd_parameters(input$type %||% "emax", isTRUE(input$delay)))
    initial_value <- function(name) if (name == "EMAX" && !grepl("_(in|out)$", input$type %||% "emax")) 50 else PD_DEFAULTS[[name]]
    output$parameters <- renderUI({
      div(class = "pd-parameter-grid", lapply(names_active(), function(name) {
        value <- isolate(input[[name]]) %||% (imported_parameters() %||% list())[[name]] %||% initial_value(name)
        if (name == "EMAX" && grepl("^inhibit", input$type %||% "emax") && value > 1) value <- 0.8
        numericInput(session$ns(name), name, value, step = if (name == "E0") 1 else 0.1)
      }))
    })
    output$estimate_parameters <- renderUI(checkboxGroupInput(session$ns("estimate"), t("Parametres a estimer", "Parameters to estimate"), names_active(), intersect(c("E0", "EMAX", "EC50", "SLOPE"), names_active())))
    configuration <- reactive({
      p <- setNames(lapply(names_active(), function(name) input[[name]] %||% (imported_parameters() %||% list())[[name]] %||% initial_value(name)), names_active())
      pd_config(list(type = input$type %||% "emax", delay = isTRUE(input$delay), parameters = p,
        exposure = input$exposure %||% "exponential", horizon = input$horizon %||% 24, c0 = input$c0 %||% 10, kel = input$kel %||% 0.1,
        profile = if (identical(input$exposure, "pk")) pk_profile() else if (identical(input$exposure, "table")) pd_read_table(input$profile, c("time", "concentration")) else NULL))
    })
    observeEvent(list(input$type, input$delay, input$exposure, input$horizon, input$c0, input$kel, input$profile, input$observations,
      input$data_mode, observations$data(), lapply(names(PD_DEFAULTS), function(name) input[[name]])), fit_store(NULL), ignoreInit = TRUE)
    observeEvent(if (identical(input$exposure, "pk")) list(tryCatch(pk(), error = function(e) NULL), input$dose, input$interval, input$infusion) else NULL, fit_store(NULL), ignoreInit = TRUE)
    observeEvent(input$simulate, fit_store(NULL))
    handle <- function(action) tryCatch(action(), error = function(e) showNotification(conditionMessage(e), type = "error", duration = 9))
    observeEvent(input$file, handle(function() {
      file <- input$file
      on.exit(unlink(file$datapath), add = TRUE)
      if (file$size > 200000) stop("CSV limit: 200 kB.")
      text <- paste(readLines(file$datapath, warn = FALSE), collapse = "\n")
      required <- c(if (identical(input$data_mode, "concentration")) "concentration" else "time", "effect")
      observations$replace(pd_read_table(text, required))
    }))
    observeEvent(input$example, handle(function() {
      config <- configuration()
      timed <- isTRUE(config$delay) || grepl("_(in|out)$", config$type)
      if (identical(input$data_mode, "concentration") && timed) {
        data <- pd_simulate(config, seq(0, config$horizon, length.out = 15))[, c("time", "effect")]
        observations$replace_for(data, c("time", "effect"))
        updateRadioButtons(session, "data_mode", selected = "time")
        showNotification(t("Le delai exige des observations temps-effet : le tableau a ete adapte.", "The delay requires time-effect observations: the table was adapted."), type = "message", duration = 5)
      } else if (identical(input$data_mode, "concentration")) {
        d <- data.frame(concentration = seq(0, config$c0, length.out = 15))
        observations$replace(transform(d, effect = pd_predict(config, d, "concentration")))
      } else {
        observations$replace(pd_simulate(config, seq(0, config$horizon, length.out = 15))[, c("time", "effect")])
      }
    }))
    observeEvent(input$use_tdm, handle(function() {
      result <- analysis_store()
      data <- result$fit_profiles$average
      if (is.null(data) || !nrow(data)) stop(t("Lancez un ajustement TDM dans cette session.", "Run a TDM fit in this session."))
      data <- data[data$time >= 0, c("time", "concentration"), drop = FALSE]
      data <- data[!duplicated(data$time, fromLast = TRUE), ]
      if (min(data$time) != 0) stop("The TDM profile must start at time 0.")
      if (nrow(data) > 2000) data <- data[unique(round(seq(1, nrow(data), length.out = 2000))), ]
      updateTextAreaInput(session, "profile", value = paste(capture.output(write.csv(data, row.names = FALSE)), collapse = "\n"))
      updateNumericInput(session, "horizon", value = max(data$time))
    }))
    observeEvent(input$fit, handle(function() {
      fit_store(NULL)
      if (!isTRUE(input$accept)) stop(t("Confirmez l'usage exploratoire.", "Confirm exploratory use."))
      mode <- input$data_mode %||% "time"
      data <- observations$data()
      fit <- withProgress(message = t("Ajustement PD", "PD fitting"), value = 0.3, pd_fit(configuration(), data, mode, input$estimate))
      fit_store(fit)
      bslib::nav_select("tabs", "fit", session = session)
    }))
    observeEvent(imported(), {
      payload <- imported(); if (!identical(payload$view, "pd")) return()
      tryCatch({
        config <- payload$config
        imported_parameters(config$parameters)
        updateSelectInput(session, "type", selected = config$type)
        updateCheckboxInput(session, "delay", value = config$delay)
        updateSelectInput(session, "exposure", selected = config$exposure)
        for (name in c("horizon", "c0", "kel")) updateNumericInput(session, name, value = config[[name]])
        if (identical(config$exposure, "pk")) for (name in c("dose", "interval", "infusion")) updateNumericInput(session, name, value = config$regimen[[name]])
        for (name in names(config$parameters)) updateNumericInput(session, name, value = config$parameters[[name]])
        fit_store(NULL)
        session$sendCustomMessage("workbench-ack", list(id = payload$id, ok = TRUE))
      }, error = function(e) session$sendCustomMessage("workbench-ack", list(id = payload$id, ok = FALSE, error = conditionMessage(e))))
    })
    session$onSessionEnded(function() { fit_store(NULL); imported_parameters(NULL) })
    current_config <- reactive(if (!is.null(fit_store())) fit_store()$config else configuration())
    curve <- reactive({
      tryCatch(pd_simulate(current_config()), error = function(e) shiny::validate(shiny::need(FALSE, conditionMessage(e))))
    })
    output$status <- renderUI({
      result <- fit_store()
      div(class = "pd-status", if (is.null(result)) t("Simulation avec les parametres saisis", "Simulation with entered parameters") else t("Parametres ajustes; conditions de l'ajustement conservees", "Fitted parameters; fitting conditions retained"))
    })
    effect_figure <- function() {
      data <- curve()
      scale <- max(1e-6, max(abs(data$effect))) / max(1e-6, max(data$concentration))
      ggplot(data, aes(time)) + geom_line(aes(y = effect, color = "E(t)", linetype = "E(t)"), linewidth = .9) +
        geom_line(aes(y = concentration * scale, color = "C(t)", linetype = "C(t)"), linewidth = .9) +
        scale_y_continuous(name = t("Effet E(t)", "Effect E(t)"), sec.axis = sec_axis(~ . / scale, name = "C(t) (mg/L)")) +
        scale_color_manual(values = c("E(t)" = "#a44430", "C(t)" = "#08796d")) +
        scale_linetype_manual(values = c("E(t)" = "solid", "C(t)" = "dashed")) +
        labs(x = t("Temps (h)", "Time (h)"), color = NULL, linetype = NULL) + theme_minimal(base_size = 13) +
        theme(legend.position = "bottom", axis.title.y.left = element_text(color = "#a44430"), axis.title.y.right = element_text(color = "#08796d"))
    }
    observe({
      payload <- imported(); if (!identical(payload$view, "pd")) return()
      data <- tryCatch(curve(), error = function(e) NULL)
      if (is.null(data)) {
        session$sendCustomMessage("workbench-result", list(id = payload$id, view = "pd", invalidated = TRUE)); return()
      }
      data <- data[unique(round(seq(1, nrow(data), length.out = min(1500, nrow(data))))), ]
      curves <- lapply(c("response", "concentration", "trajectory"), function(key) list(key = key, points = unname(lapply(seq_len(nrow(data)), function(i)
        list(x = if (key == "trajectory") data$concentration[i] else data$time[i], y = if (key == "concentration") data$concentration[i] else data$effect[i])))))
      session$sendCustomMessage("workbench-result", list(id = payload$id, view = "pd", curves = curves))
    })
    output$effect_plot <- renderPlot(effect_figure())
    output$hysteresis_plot <- renderPlot(ggplot(curve(), aes(concentration, effect, color = time)) + geom_path(linewidth = 1) +
      scale_color_gradient(low = "#c24c32", high = "#196f76") + labs(x = "Concentration", y = t("Effet", "Effect"), color = t("Temps", "Time")) + theme_minimal(base_size = 13))
    output$fit_summary <- renderUI({
      fit <- fit_store()
      if (is.null(fit)) return(p(t("Aucun ajustement disponible.", "No fit available.")))
      tagList(h3(paste("RMSE =", signif(fit$rmse, 4))), if (fit$ill_conditioned || fit$boundary) div(class = "ml-domain-warning",
        t("Parametres mal identifies ou proches d'une borne. Verifiez les valeurs initiales et les parametres fixes.", "Poorly identified parameters or a boundary solution. Review starting values and fixed parameters.")))
    })
    output$fit_parameters <- renderDT({ shiny::req(fit_store()); datatable(fit_store()$parameters, rownames = FALSE, options = list(dom = "t", scrollX = TRUE)) })
    output$fit_data <- renderDT({ shiny::req(fit_store()); datatable(fit_store()$data, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE)) })
    output$fit_plot <- renderPlot({
      fit <- fit_store(); shiny::req(fit)
      data <- fit$data
      data$x <- data[[if (fit$mode == "concentration") "concentration" else "time"]]
      ggplot(data, aes(x, effect)) + geom_point(color = "#b34930", size = 2) + geom_line(aes(y = prediction), color = "#196f76", linewidth = 1) +
        labs(x = if (fit$mode == "concentration") "Concentration" else t("Temps (h)", "Time (h)"), y = t("Effet", "Effect")) + theme_minimal(base_size = 13)
    })
    output$equation <- renderUI({
      config <- current_config()
      p(tags$code(switch(config$type, linear = "E = E0 + SLOPE*C", emax = "E = E0 + EMAX*C/(EC50+C)", hill = "E = E0 + EMAX*C^HILL/(EC50^HILL+C^HILL)",
        inhibit_in = "dR/dt = kin*(1-I(C)) - kout*R", stimulate_in = "dR/dt = kin*(1+S(C)) - kout*R", inhibit_out = "dR/dt = kin - kout*(1-I(C))*R", stimulate_out = "dR/dt = kin - kout*(1+S(C))*R")))
    })
    output$code <- renderText(pd_model_code(current_config()))
    output$download_code <- downloadHandler(filename = function() "pd-model.cpp", content = function(file) writeLines(pd_model_code(current_config()), file))
    output$download_r <- downloadHandler(filename = function() "pd-simulation.R", content = function(file) writeLines(pd_export_script(current_config()), file))
    output$download_profile <- downloadHandler(filename = function() "pd-profile.csv", content = function(file) write.csv(curve(), file, row.names = FALSE))
    output$download_report <- downloadHandler(filename = function() "pd-report.html", content = function(file) {
      config <- current_config(); fit <- fit_store()
      htmltools::save_html(tags$html(tags$head(tags$meta(charset = "utf-8"), tags$title("PD report")), tags$body(
        h1(t("Rapport pharmacodynamique", "Pharmacodynamic report")), p(t("Recherche uniquement. Ajustement individuel, pas de validation clinique.", "Research only. Individual fitting, not clinical validation.")),
        tags$pre(paste(capture.output(dput(config)), collapse = "\n")),
        if (!is.null(fit)) tagList(h2(paste("RMSE", signif(fit$rmse, 4))), tags$pre(paste(capture.output(print(fit$parameters, row.names = FALSE)), collapse = "\n")),
          p(paste("Boundary:", fit$boundary, "; ill-conditioned:", fit$ill_conditioned)),
          h2(t("Observations et residus", "Observations and residuals")), tags$pre(paste(capture.output(print(fit$data, row.names = FALSE)), collapse = "\n"))),
        tags$img(src = report_plot_uri(effect_figure()), alt = "PD profile"), tags$pre(pd_model_code(config)))), file = file)
    })
    list(configuration = configuration, fit = fit_store, curve = curve)
  })
}
