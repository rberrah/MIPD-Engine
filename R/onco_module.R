pd_panel <- function(lang = "fr") {
  t <- function(fr, en) app_t(lang, fr, en)
  nav_panel(t("Pharmacodynamie", "Pharmacodynamics"), value = "pd",
    navset_tab(id = "pd_workspace", pd_general_panel(lang), nav_panel(t("Oncologie", "Oncology"), value = "oncology", onco_ui(t)), infection_panel(lang)))
}

onco_parameter_labels <- function(t) c(V = t("Volume V (L)", "Volume V (L)"), CL = t("Clairance (L/jour)", "Clearance (L/day)"), T0 = t("Taille tumorale initiale (mm)", "Initial tumor size (mm)"), KG = t("Vitesse de croissance (1/jour)", "Growth rate (1/day)"), CAP = t("Taille limite (mm)", "Limiting size (mm)"),
  KILL = t("Effet maximal du traitement (1/jour)", "Maximum treatment effect (1/day)"), EC50 = "EC50 (mg/L)", RES = t("Perte de sensibilite (1/jour)", "Sensitivity loss (1/day)"), ANC0 = t("Neutrophiles initiaux (10^9/L)", "Initial neutrophils (10^9/L)"),
  MTT = t("MTT (jours)", "MTT (days)"), GAMMA = t("Retrocontrole gamma", "Feedback gamma"), SLOPE = "SLOPE (L/mg)")

onco_ui <- function(t) {
  ns <- NS("onco")
  div(class = "analysis-shell onco-shell",
    fileInput("onco_workshop_file", t("Importer un atelier (.json)", "Import workshop (.json)"), accept = ".json"),
    div(class = "safety-banner", tags$strong(t("Oncologie exploratoire", "Exploratory oncology")),
      span(t("Structures generiques, parametres illustratifs. Aucun protocole therapeutique valide.", "Generic structures, illustrative parameters. No validated treatment protocol."))),
    div(class = "model-chain",
      div(class = "model-block exposure", icon("syringe"), h3(t("Cycles et exposition", "Cycles and exposure")),
        selectInput(ns("pk_mode"), t("Structure PK", "PK structure"), stats::setNames(c("builtin", "free"), c(t("Exemple IV 1 compartiment", "One-compartment IV example"), t("Modele PK libre", "Free PK model"))))),
      span(class = "model-arrow", icon("arrow-right")),
      div(class = "model-block response", icon("chart-line"), h3(t("Reponse tumorale", "Tumor response")),
        selectInput(ns("growth"), t("Croissance", "Growth"), stats::setNames(c("exponential", "logistic", "gompertz"), c(t("Exponentielle", "Exponential"), t("Logistique", "Logistic"), "Gompertz")))),
      span(class = "model-arrow", icon("plus")),
      div(class = "model-block delay", icon("droplet"), h3(t("Toxicite optionnelle", "Optional toxicity")), checkboxInput(ns("toxicity"), t("Myelosuppression (Friberg)", "Myelosuppression (Friberg)"), FALSE))),
    tags$button(type = "button", class = "mobile-configure-button", onclick = "var r=this.closest('.onco-shell');var b=r.querySelector('.collapse-toggle');if(b&&b.getAttribute('aria-expanded')!=='true')b.click();", t("Parametres et donnees", "Parameters and data")),
    layout_sidebar(sidebar = sidebar(width = 370, open = "desktop",
      accordion(open = "cycles",
        accordion_panel(t("Cycles et historique", "Cycles and history"), value = "cycles",
          numericInput(ns("decision"), t("Prochaine dose / decision (jour)", "Next dose / decision (day)"), 42, min = 0),
          numericInput(ns("horizon"), t("Horizon (jours)", "Horizon (days)"), 84, min = 1, max = 730),
          tags$details(tags$summary(t("Doses passees", "Past doses")),
          pd_observations_ui(ns("history"), t),
          tags$small(t("Perfusion : 0 h = oral / bolus IV. Doses strictement anterieures au jour de decision.", "Infusion: 0 h = oral / IV bolus. Doses strictly before the decision day."))),
          numericInput(ns("dose"), t("Dose future de reference (unite PK)", "Reference future dose (PK unit)"), 100, min = 0),
          numericInput(ns("interval"), t("Intervalle de reference (jours)", "Reference interval (days)"), 21, min = 0.25),
          numericInput(ns("infusion"), t("Perfusion (h; 0 = oral / bolus IV)", "Infusion (h; 0 = oral / IV bolus)"), 1, min = 0)),
        accordion_panel(t("Modele PK libre", "Free PK model"), value = "free_pk",
          conditionalPanel(sprintf("input['%s']=='free'", ns("pk_mode")), pd_pk_ui(ns("pk"), t)),
          p(t("Les cycles utilisent des jours. La PK de la bibliotheque est automatiquement convertie en mg/L et synchronisee avec les cycles.", "Cycles use days. Library PK is automatically converted to mg/L and synchronized with cycles."))),
        accordion_panel(t("Parametres PK/PD", "PK/PD parameters"), value = "parameters", uiOutput(ns("parameters"))),
        accordion_panel(t("Observations et ajustement", "Observations and fitting"), value = "observations",
          fileInput(ns("file"), t("Importer CSV", "Import CSV"), accept = ".csv"),
          pd_observations_ui(ns("observations"), t),
          tags$small(t("Jour; tumor = somme des diametres (mm); anc = neutrophiles (10^9/L). Valeurs strictement positives, sans observation posterieure a la decision.", "Day; tumor = sum of diameters (mm); anc = neutrophils (10^9/L). Strictly positive values; no observation after the decision.")),
          actionButton(ns("example"), t("Generer des observations synthetiques", "Generate synthetic observations"), icon = icon("flask")),
          uiOutput(ns("estimate_ui")),
          tags$details(tags$summary(t("Erreur des observations", "Observation error")),
            numericInput(ns("sigma_tumor"), t("Ecart-type log : tumeur", "Log SD: tumor"), 0.15, min = 0.001),
            conditionalPanel(sprintf("input['%s']", ns("toxicity")), numericInput(ns("sigma_anc"), t("Ecart-type log : ANC", "Log SD: ANC"), 0.2, min = 0.001))),
          checkboxInput(ns("accept"), t("Usage de recherche uniquement", "Research use only"), FALSE),
          actionButton(ns("fit"), t("Ajuster sur les observations", "Fit observations"), icon = icon("sliders")),
          actionButton(ns("reset_fit"), t("Revenir aux parametres saisis", "Return to entered parameters"), icon = icon("rotate-left"))),
        accordion_panel(t("Cibles et grille de doses", "Targets and dose grid"), value = "future",
          fluidRow(column(4, numericInput(ns("dose_min"), "Dose min", 50, min = 0)),
            column(4, numericInput(ns("dose_max"), "Dose max", 150, min = 0)),
            column(4, numericInput(ns("dose_step"), t("Pas", "Step"), 25, min = .001))),
          selectizeInput(ns("intervals"), t("Intervalles proposes (jours)", "Candidate intervals (days)"), c(7, 14, 21, 28), selected = 21, multiple = TRUE, options = list(create = TRUE)),
          numericInput(ns("delay"), t("Report de la prochaine dose (jours)", "Delay next dose (days)"), 0, min = 0, max = 180),
          conditionalPanel(sprintf("input['%s']", ns("toxicity")), numericInput(ns("anc_floor"), t("Plancher ANC exploratoire (10^9/L)", "Exploratory ANC floor (10^9/L)"), 1, min = 0)),
          numericInput(ns("tumor_goal"), t("Plafond final tumeur / baseline", "Final tumor / baseline ceiling"), 0.8, min = 0.001))),
      actionButton(ns("compare"), t("Evaluer la grille de doses", "Evaluate dose grid"), icon = icon("play"), class = "btn-primary w-100")),
      div(class = "workspace", uiOutput(ns("status")), navset_tab(id = ns("tabs"),
        nav_panel(t("Reponse et toxicite", "Response and toxicity"), value = "simulation",
          uiOutput(ns("empty")), conditionalPanel(sprintf("output['%s']", ns("ready")),
            uiOutput(ns("candidate_ui")), DTOutput(ns("grid_metrics")),
            plotOutput(ns("tumor_plot"), height = "310px"),
            conditionalPanel(sprintf("input['%s']", ns("toxicity")), plotOutput(ns("anc_plot"), height = "280px")),
            plotOutput(ns("pk_plot"), height = "230px"), DTOutput(ns("metrics")), DTOutput(ns("schedule")),
            downloadButton(ns("download_curves"), t("Exporter les courbes", "Export curves")))),
        nav_panel(t("Ajustement", "Fit"), value = "fit", uiOutput(ns("fit_status")), DTOutput(ns("fit_parameters")), DTOutput(ns("fit_data"))),
        nav_panel(t("Code et rapport", "Code and report"), value = "code", verbatimTextOutput(ns("code")),
          downloadButton(ns("download_code"), "mrgsolve .cpp"), downloadButton(ns("download_r"), "R / deSolve"),
          conditionalPanel(sprintf("output['%s']", ns("ready")), downloadButton(ns("download_report"), t("Rapport HTML", "HTML report")))),
        nav_panel(t("Methodes et limites", "Methods and limitations"),
          p(t("PK IV d'exemple, modele populationnel avec covariables, ou estimation bayesienne mapbayr depuis Analyse. Cycles en jours, concentrations en mg/L; unites de la bibliotheque automatiques. La PK reste fixe pendant l'ajustement PD. Les scenarios maintenu et modifie partagent l'historique; sans traitement exclut toutes les doses depuis t=0.",
            "Built-in IV PK, population model with covariates, or mapbayr Bayesian estimate from Analysis. Cycles use days and concentrations mg/L; library units are automatic. PK remains fixed during PD fitting. Maintained and changed scenarios share history; untreated excludes all doses since t=0.")),
          p("dT/dt = growth(T) - KILL*C/(EC50+C)*exp(-RES*t)*T; ktr = 4/MTT."),
          p(t("La myelosuppression reprend les 5 compartiments et le retrocontrole de Friberg. Modification explicite : inhibition lineaire bornee a 1, avec valeurs initiales toutes egales a ANC0. Ce n'est pas une reproduction d'un modele medicament particulier. RES est une perte empirique de sensibilite avec le temps, pas une population cellulaire resistante.",
            "Myelosuppression uses Friberg's 5 compartments and feedback. Explicit modification: linear inhibition capped at 1, all initial counts equal ANC0. This is not a reproduction of a particular drug model. RES is empirical loss of sensitivity with time, not a resistant cell population.")),
          p(t("Ajustement PD par moindres carres ponderes sur log(valeur), pas de MAP-BE PD ni de propagation de l'incertitude PK. Classification deterministe des cibles sur une grille de 0,05 jour, pas une probabilite de succes. Aucun classement en recommandation. Sans toxicite active, aucune conclusion de securite n'est possible.",
            "PD fitting uses weighted least squares on log(value), not PD MAP-BE or PK uncertainty propagation. Deterministic target classification on a 0.05-day grid, not a success probability. No ranking as recommendations. Without active toxicity, no safety conclusion is possible.")),
          p(t("La taille est une somme des diametres (SLD), pas un volume. Aucun classement RECIST ni prediction de survie. Pas de G-CSF, de toxicite non hematologique, de combinaisons anticancereuses ni de parametres cliniques valides par molecule. Les donnees restent dans la session, sauf telechargement explicite.",
            "Size is sum of longest diameters (SLD), not volume. No RECIST classification or survival prediction. No G-CSF, nonhematological toxicity, anticancer combinations, or drug-specific validated clinical parameters. Data remain in session unless explicitly downloaded.")),
          tags$a(href = "https://doi.org/10.1200/JCO.2002.02.140", target = "_blank", rel = "noopener noreferrer", "Friberg et al., 2002"))))))
}

onco_server <- function(id, report_plot_uri, imported = reactive(NULL), soloc = tempdir(), cache = new.env(), analysis_store = reactive(NULL), open_tdm = NULL) {
  moduleServer(id, function(input, output, session) {
    t <- function(fr, en) app_t(app_language_from_query(session$clientData$url_search %||% ""), fr, en)
    fit_store <- reactiveVal(NULL); grid_results <- reactiveVal(NULL); imported_parameters <- reactiveVal(NULL)
    dose_grid <- reactive(list(min = input$dose_min, max = input$dose_max, step = input$dose_step, intervals = as.numeric(input$intervals)))
    comparison <- reactive({
      result <- grid_results(); if (is.null(result)) return(NULL)
      choices <- setdiff(result$regimens$scenario, "maintain")
      selected <- input$candidate
      if (is.null(selected) || !selected %in% choices) selected <- choices[[1]]
      onco_grid_selection(result, selected)
    })
    pk <- pd_pk_server("pk", soloc, cache, reactive(if (identical(imported()$view, "onco")) imported()$models[[1]] else NULL), analysis_store, open_tdm)
    pk_context <- reactive(if (identical(input$pk_mode, "free")) pk() else NULL)
    history <- pd_observations_server("history", reactive(c("time", "amount", "infusion")),
      data.frame(time = c(0, 21), amount = c(100, 100), infusion = c(1, 1)), dose_history = TRUE)
    observations <- pd_observations_server("observations", reactive(c("time", "endpoint", "value")),
      data.frame(time = c(0, 14, 28, 42, 0, 7, 14, 28), endpoint = rep(c("tumor", "anc"), each = 4), value = c(60, 52, 46, 44, 4, 2, 3, 2)), oncology = TRUE)
    active <- reactive(c("V", "CL", "T0", "KG", if ((input$growth %||% "exponential") != "exponential") "CAP", "KILL", "EC50", "RES", if (isTRUE(input$toxicity)) c("ANC0", "MTT", "GAMMA", "SLOPE")))
    output$parameters <- renderUI({
      groups <- list(initial = intersect(c("T0", "ANC0"), active()), pk = if (!identical(input$pk_mode, "free")) c("V", "CL"), pd = setdiff(active(), c("T0", "ANC0", "V", "CL")))
      labels <- c(initial = t("Conditions initiales", "Initial conditions"), pk = t("Parametres PK", "PK parameters"), pd = t("Parametres PD", "PD parameters"))
      control <- function(name) numericInput(session$ns(name), onco_parameter_labels(t)[[name]],
        isolate(input[[name]]) %||% (imported_parameters() %||% list())[[name]] %||% ONCO_DEFAULTS[[name]], min = 0)
      tagList(lapply(names(groups), function(group) if (length(groups[[group]])) tags$fieldset(tags$legend(labels[[group]]),
        div(class = "pd-parameter-grid", lapply(setdiff(groups[[group]], c("RES", "GAMMA")), control)),
        if (group == "pd") tags$details(tags$summary(t("Sensibilite et retrocontrole avances", "Advanced sensitivity and feedback")), lapply(intersect(groups[[group]], c("RES", "GAMMA")), control)))))
    })
    output$estimate_ui <- renderUI(checkboxGroupInput(session$ns("estimate"), t("Parametres a estimer", "Parameters to estimate"), setdiff(active(), c("V", "CL")), intersect(c("KILL", "SLOPE"), active())))
    configuration <- reactive({
      p <- ONCO_DEFAULTS
      for (name in active()) p[[name]] <- input[[name]] %||% (imported_parameters() %||% list())[[name]] %||% ONCO_DEFAULTS[[name]]
      values <- list(parameters = p, growth = input$growth %||% "exponential", toxicity = isTRUE(input$toxicity), free_pk = identical(input$pk_mode, "free"),
        history = history$data())
      for (name in c("decision", "horizon", "dose", "interval", "infusion", "anc_floor", "tumor_goal", "sigma_tumor", "sigma_anc")) values[[name]] <- input[[name]]
      onco_config(values)
    })
    handle <- function(fn) tryCatch(fn(), error = function(e) showNotification(conditionMessage(e), type = "error", duration = 10))
    observeEvent(list(input$growth, input$toxicity, history$data(), input$decision, input$horizon, input$dose, input$interval, input$infusion,
      input$pk_mode, observations$data(), input$sigma_tumor, input$sigma_anc, lapply(names(ONCO_DEFAULTS), function(name) input[[name]])), {
      fit_store(NULL); grid_results(NULL)
    }, ignoreInit = TRUE)
    observeEvent(tryCatch(pk_context(), error = function(e) NULL), { fit_store(NULL); grid_results(NULL) }, ignoreInit = TRUE, ignoreNULL = FALSE)
    observeEvent(list(dose_grid(), input$delay, input$anc_floor, input$tumor_goal), grid_results(NULL), ignoreInit = TRUE)
    observeEvent(input$reset_fit, { fit_store(NULL); grid_results(NULL) })
    current <- reactive({
      config <- if (is.null(fit_store())) configuration() else fit_store()$config
      config$anc_floor <- input$anc_floor %||% 1; config$tumor_goal <- input$tumor_goal %||% 0.8
      config
    })
    observeEvent(input$file, handle(function() {
      on.exit(unlink(input$file$datapath), add = TRUE)
      if (input$file$size > 200000) stop("CSV limit: 200 kB.")
      text <- paste(readLines(input$file$datapath, warn = FALSE), collapse = "\n")
      observations$replace(onco_read_observations(text))
    }))
    observeEvent(input$example, handle(function() {
      config <- configuration(); times <- seq(0, config$decision, length.out = 13)
      result <- onco_simulate(config, doses = config$history, times = times, pk_context = pk_context())
      data <- data.frame(time = times, endpoint = "tumor", value = stats::approx(result$time, result$TUMOR, times, rule = 2)$y)
      if (config$toxicity) data <- rbind(data, data.frame(time = times, endpoint = "anc", value = stats::approx(result$time, result$ANC, times, rule = 2)$y))
      observations$replace(data)
    }))
    observeEvent(input$fit, handle(function() {
      fit_store(NULL); grid_results(NULL)
      if (!isTRUE(input$accept)) stop(t("Confirmez l'usage de recherche.", "Confirm research use."))
      result <- withProgress(message = t("Ajustement oncologique", "Oncology fitting"), value = 0.2,
        onco_fit(configuration(), observations$data(), input$estimate, pk_context = pk_context()))
      fit_store(result); bslib::nav_select("tabs", "fit", session = session)
    }))
    observeEvent(input$compare, handle(function() {
      grid_results(NULL)
      result <- withProgress(message = t("Evaluation de la grille", "Evaluating dose grid"), value = 0.3,
        onco_compare_grid(current(), dose_grid(), input$delay, pk_context = pk_context()))
      grid_results(result); bslib::nav_select("tabs", "simulation", session = session)
    }))
    output$candidate_ui <- renderUI({
      result <- grid_results(); shiny::req(result)
      rows <- result$regimens[result$regimens$scenario != "maintain", ]
      selectInput(session$ns("candidate"), t("Posologie comparee", "Compared regimen"),
        setNames(rows$scenario, paste(rows$dose, "/", rows$interval, t("jours", "days"))), selected = isolate(input$candidate) %||% rows$scenario[[1]])
    })
    output$grid_metrics <- renderDT({
      result <- grid_results(); shiny::req(result)
      data <- result$metrics[, c("dose", "interval", "final_tumor_ratio", if (result$config$toxicity) "future_anc_nadir", "targets_met")]
      data$targets_met <- ifelse(data$targets_met, t("Oui", "Yes"), t("Non", "No"))
      names(data) <- c(t("Dose (unite PK)", "Dose (PK unit)"), t("Intervalle (jours)", "Interval (days)"), t("Ratio tumoral final", "Final tumor ratio"), if (result$config$toxicity) t("Nadir ANC futur", "Future ANC nadir"), t("Cibles atteintes", "Targets met"))
      data <- cbind(setNames(data.frame(ifelse(result$metrics$scenario == "maintain", t("Maintenir", "Maintain"), t("Candidate", "Candidate"))), t("Posologie", "Regimen")), data)
      datatable(data, rownames = FALSE, options = list(dom = "t", scrollX = TRUE, paging = FALSE)) |> DT::formatRound(columns = which(vapply(data, is.numeric, logical(1))), digits = 3)
    })
    observe({
      if (identical(imported()$view, "onco")) {
        result <- comparison()
        if (is.null(result)) { session$sendCustomMessage("workbench-result", list(id = imported()$id, view = "onco", invalidated = TRUE)); return() }
        curves <- lapply(c("untreated", "maintain"), function(name) {
          data <- result$curves[[name]]
          data <- data[unique(round(seq(1, nrow(data), length.out = min(1500, nrow(data))))), ]
          list(key = if (name == "maintain") "treated" else "untreated", points = unname(lapply(seq_len(nrow(data)), function(i) list(x = data$time[i], y = data$TUMOR[i]))))
        })
        session$sendCustomMessage("workbench-result", list(id = imported()$id, view = "onco", curves = curves))
      }
    })
    output$status <- renderUI(p(if (is.null(fit_store())) t("Parametres saisis, sans individualisation.", "Entered parameters, no individualization.") else
      t("Parametres ajustes sur l'historique disponible au jour de decision.", "Parameters fitted to history available at the decision day.")))
    output$ready <- reactive(!is.null(comparison()))
    outputOptions(output, "ready", suspendWhenHidden = FALSE)
    output$empty <- renderUI(if (is.null(comparison())) p(t("Aucune comparaison calculee pour cette configuration.", "No comparison calculated for this configuration.")))
    figure <- function(endpoint) {
      result <- comparison(); shiny::req(result)
      if (endpoint == "ANC") shiny::req(result$config$toxicity)
      scenarios <- if (endpoint == "TUMOR") names(result$curves) else setdiff(names(result$curves), "untreated")
      data <- do.call(rbind, lapply(scenarios, function(name) transform(result$curves[[name]], scenario = name)))
      data$value <- data[[endpoint]]
      labels <- c(untreated = t("Sans traitement", "Untreated"), maintain = t("Maintenir", "Maintain"), change = t("Modifier", "Change"))
      graph <- ggplot(data, aes(time, value, color = scenario, linetype = scenario)) + geom_line(linewidth = 0.9) +
        geom_vline(xintercept = result$config$decision, linetype = "dotted", color = "#444444") +
        scale_color_manual(values = c(untreated = "#62676c", maintain = "#a44430", change = "#08796d"), labels = labels) +
        scale_linetype_manual(values = c(untreated = "dashed", maintain = "solid", change = "solid"), labels = labels) +
        labs(x = t("Jour", "Day"), y = switch(endpoint, TUMOR = "SLD (mm)", ANC = "ANC (10^9/L)", concentration = "C (mg/L)"), color = NULL, linetype = NULL) +
        theme_minimal(base_size = 13) + theme(legend.position = "bottom")
      if (endpoint == "ANC") graph <- graph + geom_hline(yintercept = result$config$anc_floor, color = "#7c3e58", linetype = "dotdash")
      if (endpoint == "TUMOR") graph <- graph + geom_hline(yintercept = result$config$tumor_goal * result$config$parameters$T0, color = "#7c3e58", linetype = "dotdash")
      fit <- fit_store()
      if (!is.null(fit) && endpoint != "concentration") {
        observed <- fit$data[fit$data$endpoint == if (endpoint == "ANC") "anc" else "tumor", ]
        graph <- graph + geom_point(data = observed, aes(time, value), inherit.aes = FALSE, color = "#151515", size = 2)
      }
      graph
    }
    output$tumor_plot <- renderPlot(figure("TUMOR")); output$anc_plot <- renderPlot(figure("ANC")); output$pk_plot <- renderPlot(figure("concentration"))
    output$metrics <- renderDT({
      shiny::req(comparison()); data <- comparison()$metrics
      data$scenario <- ifelse(data$scenario == "maintain", t("Maintenir", "Maintain"), t("Modifier", "Change"))
      names(data) <- c(t("Scenario", "Scenario"), t("Ratio tumoral final", "Final tumor ratio"), t("Nadir ANC futur", "Future ANC nadir"),
        t("Jours sous plancher", "Days below floor"), t("Cible tumorale atteinte", "Tumor goal met"), t("Cible ANC atteinte", "ANC goal met"), t("Dose future totale (unite PK)", "Total future dose (PK unit)"))
      if (!comparison()$config$toxicity) data <- data[, -c(3, 4, 6), drop = FALSE]
      numeric_columns <- which(vapply(data, is.numeric, logical(1)))
      datatable(data, rownames = FALSE, options = list(dom = "t", scrollX = TRUE)) |> DT::formatRound(columns = numeric_columns, digits = 3)
    })
    output$schedule <- renderDT({ shiny::req(comparison()); datatable(do.call(rbind, lapply(names(comparison()$schedules), function(name) transform(comparison()$schedules[[name]], scenario = name))), rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE)) })
    output$fit_status <- renderUI({
      fit <- fit_store()
      if (is.null(fit)) return(p(t("Aucun ajustement.", "No fit.")))
      tagList(p(paste(t("Somme des residus log ponderes au carre :", "Sum of squared weighted log residuals:"), signif(fit$objective, 4))),
        if (fit$ill_conditioned || fit$boundary) div(class = "ml-domain-warning", t("Identifiabilite faible ou borne atteinte : interpretation prudente.", "Poor identifiability or boundary reached: interpret cautiously.")))
    })
    output$fit_parameters <- renderDT({ shiny::req(fit_store()); datatable(fit_store()$parameters, rownames = FALSE, options = list(dom = "t", scrollX = TRUE)) })
    output$fit_data <- renderDT({ shiny::req(fit_store()); datatable(fit_store()$data, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE)) })
    output$code <- renderText(onco_model_code(current()))
    output$download_code <- downloadHandler(filename = function() "oncology.cpp", content = function(file) writeLines(onco_model_code(current()), file))
    output$download_r <- downloadHandler(filename = function() "oncology.R", content = function(file) writeLines(onco_export_script(current(), delay = input$delay, pk_context = pk_context(), grid = dose_grid()), file))
    output$download_curves <- downloadHandler(filename = function() "oncology-curves.csv", content = function(file) {
      shiny::req(comparison()); write.csv(do.call(rbind, lapply(names(comparison()$curves), function(name) transform(comparison()$curves[[name]], scenario = name))), file, row.names = FALSE)
    })
    output$download_report <- downloadHandler(filename = function() "oncology-report.html", content = function(file) {
      result <- comparison(); shiny::req(result)
      htmltools::save_html(tags$html(tags$head(tags$meta(charset = "utf-8"), tags$title("Oncology report")), tags$body(
        h1(t("Oncologie : comparaison exploratoire", "Oncology: exploratory comparison")),
        p(t("Aucune recommandation clinique. Structures generiques, sans validation par molecule ni incertitude predictive.", "No clinical recommendation. Generic structures, no drug-specific validation or predictive uncertainty.")),
        tags$pre(paste(capture.output(dput(grid_results()[c("config", "pk", "grid", "delay", "metrics", "schedules")])), collapse = "\n")),
        tags$pre(paste(capture.output(dput(list(compared_regimen = result$change))), collapse = "\n")),
        tags$img(src = report_plot_uri(figure("TUMOR")), alt = "Tumor response"),
        if (result$config$toxicity) tags$img(src = report_plot_uri(figure("ANC")), alt = "ANC"),
        if (!is.null(fit_store())) tags$pre(paste(capture.output(dput(fit_store())), collapse = "\n")),
        tags$a(href = "https://doi.org/10.1200/JCO.2002.02.140", "Friberg 2002: structural basis"))), file = file)
    })
    observeEvent(imported(), {
      payload <- imported(); if (!identical(payload$view, "onco")) return()
      tryCatch({
        config <- onco_config(payload$config)
        imported_parameters(config$parameters)
        updateSelectInput(session, "growth", selected = config$growth)
        updateSelectInput(session, "pk_mode", selected = if (config$free_pk) "free" else "builtin")
        updateCheckboxInput(session, "toxicity", value = config$toxicity)
        for (name in c("decision", "horizon", "dose", "interval", "infusion", "anc_floor", "tumor_goal", "sigma_tumor", "sigma_anc")) updateNumericInput(session, name, value = config[[name]])
        for (name in names(ONCO_DEFAULTS)) updateNumericInput(session, name, value = config$parameters[[name]])
        updateNumericInput(session, "dose_min", value = config$dose * .5)
        updateNumericInput(session, "dose_max", value = config$dose * 1.5)
        updateNumericInput(session, "dose_step", value = max(.001, config$dose * .25))
        updateSelectizeInput(session, "intervals", choices = sort(unique(c(7, 14, 21, 28, config$interval))), selected = config$interval)
        history$replace(config$history)
        fit_store(NULL); grid_results(NULL)
        session$sendCustomMessage("workbench-ack", list(id = payload$id, ok = TRUE))
      }, error = function(e) session$sendCustomMessage("workbench-ack", list(id = payload$id, ok = FALSE, error = conditionMessage(e))))
    })
    session$onSessionEnded(function() { fit_store(NULL); grid_results(NULL); imported_parameters(NULL) })
    list(configuration = configuration, fit = fit_store, comparison = comparison)
  })
}
