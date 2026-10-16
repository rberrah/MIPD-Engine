infection_methods <- function(t) tagList(
  p(t("Exploration PK/PD, non validee pour une decision therapeutique. Cibles et fraction libre sont saisies par l'utilisateur, pas des recommandations. La concentration du modele doit etre convertie en mg/L et correspondre a la matrice de la cible.",
    "Exploratory PK/PD, not validated for treatment decisions. Targets and unbound fraction are user-defined, not recommendations. Model concentration must be converted to mg/L and match the target matrix.")),
  p(t("PTA populationnelle : tirages des ETA selon OMEGA, covariables fixes. Apres TDM : approximation normale de la distribution posterieure des ETA, covariance MAP requise. Model averaging : melange des distributions selon les poids TDM. Aucune erreur residuelle ajoutee. Les deux posologies utilisent les memes tirages aleatoires, a l'etat stationnaire.",
    "Population PTA: ETA draws from OMEGA, fixed covariates. After TDM: normal approximation to posterior ETA uncertainty, requiring MAP covariance. Model averaging mixes distributions using TDM weights. No residual measurement error. Both regimens use the same random draws at steady state.")),
  p(t("Indices : %T > k x CMI sur un intervalle de dose, AUC sur exactement 24 h / CMI, Cmax sur un intervalle / CMI. Fraction libre constante : C libre = fu x C totale. Pour une sortie deja libre, fu = 1. Pas de propagation de l'incertitude sur fu ou sur la CMI; la sensibilite CMI/2 et 2 x CMI est affichee. Les bandes de PTA sont des intervalles de Wilson a 95 % sur l'erreur Monte Carlo uniquement.",
    "Indices: %T > k x MIC over one dosing interval, AUC over exactly 24 h / MIC, Cmax over one interval / MIC. Constant unbound fraction: free C = fu x total C. For an already unbound output use fu = 1. No fu or MIC uncertainty propagation; MIC/2 and 2 x MIC sensitivity is shown. PTA bands are 95% Wilson intervals for Monte Carlo error only.")),
  p(t("L'historique TDM reste une estimation ponctuelle sur la fenetre reellement disponible, pas une PTA. Une AUC sur moins de 24 h n'est pas convertie en AUC24. Les schemas futurs sont stationnaires, sans transition ni dose de charge. Une grille fine limite, sans supprimer, l'erreur numerique sur Cmax et le temps au-dessus du seuil. Aucun modele de toxicite ou d'emergence de resistance.",
    "TDM history remains a point estimate over the actually available window, not a PTA. An AUC spanning less than 24 h is not converted into AUC24. Future regimens are at steady state, without a transition or loading dose. A fine time grid reduces but does not eliminate Cmax and time-above-threshold numerical error. No toxicity or resistance-emergence model.")),
  p(t("EUCAST : distributions aggregees de sources et periodes multiples, pas des taux locaux de resistance. ECOFF / TECOFF ne sont ni une CMI patient ni un seuil clinique S/I/R. Les effectifs des classes extremes peuvent etre censures : aucune couverture empirique ou CFR n'est calculee a partir de ces effectifs. Une CMI renseignee comme exacte ne doit pas etre une borne > ou <=.",
    "EUCAST: aggregated distributions from multiple sources and periods, not local resistance rates. ECOFF / TECOFF are neither patient MICs nor S/I/R clinical breakpoints. Extreme-bin counts may be censored: no empirical coverage or CFR is calculated from these counts. An MIC entered as exact must not be a > or <= bound.")),
  p(t("EUCAST n'est contacte que sur demande, avec le choix d'antimicrobien, jamais avec les donnees patient. Les resultats et modeles restent dans la session, sauf export explicite. Le connecteur lit les tableaux HTML publics; en cas d'indisponibilite, la saisie de CMI reste possible.",
    "EUCAST is contacted only on request, with the antimicrobial selection, never patient data. Results and models stay in session unless explicitly exported. The connector reads public HTML tables; manual MIC entry remains available during an outage.")),
  tags$a(href = "https://mic.eucast.org/", target = "_blank", rel = "noopener noreferrer", "EUCAST MIC distributions / ECOFF"), " | ",
  tags$a(href = "https://doi.org/10.1093/jac/dki079", target = "_blank", rel = "noopener noreferrer", "Mouton et al., 2005: PK/PD terminology"), " | ",
  tags$a(href = "https://doi.org/10.1128/AAC.46.3.913-916.2002", target = "_blank", rel = "noopener noreferrer", "Drusano et al., 2002: population PK / Monte Carlo"), " | ",
  tags$a(href = "https://doi.org/10.1086/383320", target = "_blank", rel = "noopener noreferrer", "Drusano et al., 2004: AUC/MIC / eradication"), " | ",
  tags$a(href = "https://abis.chu-limoges.fr/", target = "_blank", rel = "noopener noreferrer", "ABIS - CHU Limoges"))

infection_panel <- function(lang = "fr") {
  t <- function(fr, en) app_t(lang, fr, en); ns <- NS("infection")
  nav_panel(t("Infectiologie", "Infectiology"), value = "infection",
    div(class = "analysis-shell infection-shell",
      fileInput("infection_workshop_file", t("Importer un atelier (.json)", "Import workshop (.json)"), accept = ".json"),
      div(class = "safety-banner", strong(t("Infectiologie exploratoire", "Exploratory infectiology")),
        span(t("Cibles definies par l'utilisateur. Aucune recommandation clinique automatique.", "User-defined targets. No automatic clinical recommendation."))),
      div(class = "model-chain",
        div(class = "model-block exposure", icon("syringe"), h3(t("Exposition", "Exposure")),
          radioButtons(ns("source"), NULL, setNames(c("pk", "tdm"), c(t("Modele populationnel", "Population model"), "TDM / MAP-BE"))),
          conditionalPanel(sprintf("input['%s']=='pk'", ns("source")), selectInput(ns("pk_mode"), t("Structure PK", "PK structure"),
            setNames(c("builtin", "free"), c(t("IV 1 compartiment", "IV one compartment"), t("Modele PK libre", "Free PK model")))))),
        span(class = "model-arrow", icon("arrow-right")),
        div(class = "model-block delay", icon("microscope"), h3(t("CMI", "MIC")),
          numericInput(ns("mic"), t("CMI (mg/L)", "MIC (mg/L)"), 1, min = 1e-6)),
        span(class = "model-arrow", icon("arrow-right")),
        div(class = "model-block response", icon("chart-line"), h3(t("Cible PK/PD", "PK/PD target")),
          selectInput(ns("metric"), t("Indice", "Index"), setNames(c("time", "auc", "peak"), c(t("%T > k x CMI", "%T > k x MIC"), t("AUC24 / CMI (h)", "AUC24 / MIC (h)"), t("Cmax / CMI", "Cmax / MIC")))),
          numericInput(ns("target"), t("Temps au-dessus de k x CMI (%)", "Time above k x MIC (%)"), 100, min = .001, max = 100))),
      tags$button(type = "button", class = "mobile-configure-button", onclick = "var r=this.closest('.infection-shell');var b=r.querySelector('.collapse-toggle');if(b&&b.getAttribute('aria-expanded')!=='true')b.click();", t("Patient, CMI et doses", "Patient, MIC and doses")),
      layout_sidebar(sidebar = sidebar(width = 380, open = "desktop",
        accordion(open = "pk",
          accordion_panel("PK / TDM", value = "pk",
            conditionalPanel(sprintf("input['%s']=='pk' && input['%s']=='builtin'", ns("source"), ns("pk_mode")),
              numericInput(ns("v"), "V (L)", 20, min = 1e-6), numericInput(ns("cl"), "CL (L/h)", 4, min = 1e-6),
              p(t("Variances ETA independantes sur CL et V : 0,09 (CV 30,7 %). Exemple pedagogique, sans erreur residuelle dans la PTA.", "Independent CL and V ETA variances: 0.09 (CV 30.7%). Teaching example, no residual error in PTA."))),
            conditionalPanel(sprintf("input['%s']=='pk' && input['%s']=='free'", ns("source"), ns("pk_mode")), pd_pk_ui(ns("pk"), t, allow_tdm = FALSE)),
            conditionalPanel(sprintf("input['%s']=='tdm'", ns("source")),
              uiOutput(ns("tdm_source")), uiOutput(ns("tdm_units")),
              actionButton(ns("estimate"), t("Donnees patient et estimation mapbayr", "Patient data and mapbayr estimation"), icon = icon("chart-line")),
              actionButton(ns("last_regimen"), t("Reprendre la derniere posologie", "Use the latest regimen"), icon = icon("arrow-right")))) ,
          accordion_panel(t("CMI et cible PK/PD", "MIC and PK/PD target"), value = "target",
            radioButtons(ns("basis"), t("Concentration de la cible", "Target concentration basis"), setNames(c("free", "total"), c(t("Libre", "Unbound"), t("Totale", "Total"))), inline = TRUE),
            conditionalPanel(sprintf("input['%s']=='free'", ns("basis")), numericInput(ns("fu"), t("Fraction libre fu (0-1)", "Unbound fraction fu (0-1)"), 1, min = 1e-6, max = 1)),
            conditionalPanel(sprintf("input['%s']=='time'", ns("metric")), numericInput(ns("multiple"), "k", 1, min = .01)),
            p(t("Exemple numerique, pas une cible recommandee. fu = 1 si la sortie du modele est deja libre.", "Numerical example, not a recommended target. fu = 1 if the model output is already unbound."))),
          accordion_panel(t("Posologie actuelle", "Current regimen"), value = "regimens",
            numericInput(ns("dose"), t("Dose (unite PK)", "Dose (PK unit)"), 1000, min = 1e-6),
            numericInput(ns("interval"), t("Intervalle (h)", "Interval (h)"), 12, min = .25, max = 8760),
            numericInput(ns("infusion"), t("Perfusion (h; 0 = bolus IV / oral)", "Infusion (h; 0 = IV bolus / oral)"), 1, min = 0)),
          accordion_panel(t("Cible et grille de doses", "Target and dose grid"), value = "grid",
            numericInput(ns("pta_target"), t("Probabilite d'atteinte souhaitee (%)", "Desired target attainment probability (%)"), 90, min = .001, max = 100),
            fluidRow(column(4, numericInput(ns("dose_min"), "Dose min", 500, min = .001)),
              column(4, numericInput(ns("dose_max"), "Dose max", 1500, min = .001)),
              column(4, numericInput(ns("dose_step"), t("Pas", "Step"), 500, min = .001))),
            selectizeInput(ns("intervals"), t("Intervalles proposes (h)", "Candidate intervals (h)"), c(4, 6, 8, 12, 24, 48), c(8, 12), multiple = TRUE,
              options = list(create = TRUE, placeholder = t("Choisir ou saisir une valeur", "Select or enter a value"))),
            checkboxInput(ns("continuous"), t("Perfusion continue", "Continuous infusion"), FALSE),
            conditionalPanel(sprintf("!input['%s']", ns("continuous")), numericInput(ns("grid_infusion"), t("Perfusion (h; 0 = bolus IV / oral)", "Infusion (h; 0 = IV bolus / oral)"), 1, min = 0)))),
        tags$details(tags$summary(t("Simulation avancee", "Advanced simulation")),
            numericInput(ns("replicates"), t("Tirages Monte Carlo", "Monte Carlo draws"), 250, min = 50, max = 1000, step = 50),
            numericInput(ns("delta"), t("Pas temporel (h)", "Time step (h)"), .05, min = .01, max = .25, step = .01)),
        checkboxInput(ns("accept"), t("Unites, matrice, fraction libre et cible verifiees. Usage de recherche.", "Units, matrix, unbound fraction and target checked. Research use."), FALSE),
        conditionalPanel(sprintf("input['%s']=='tdm'", ns("source")), actionButton(ns("interpret"), t("Interpreter l'historique TDM", "Interpret TDM history"), icon = icon("magnifying-glass"))),
        actionButton(ns("simulate"), t("Evaluer la grille de doses", "Evaluate dose grid"), icon = icon("play"), class = "btn-primary w-100")),
        div(class = "workspace", uiOutput(ns("status")), navset_tab(
          nav_panel("PTA", uiOutput(ns("comparison")), plotOutput(ns("pta"), height = "320px"), DTOutput(ns("sensitivity")),
            plotOutput(ns("exposure"), height = "280px"), DTOutput(ns("regimens")),
            downloadButton(ns("csv"), t("Exporter PTA (CSV)", "Export PTA (CSV)")), downloadButton(ns("report"), t("Rapport HTML", "HTML report"))),
          nav_panel(t("Historique TDM", "TDM history"), uiOutput(ns("history_summary")), plotOutput(ns("history"), height = "300px")),
          nav_panel("EUCAST", p(t("Choisissez l'antimicrobien correspondant au modele PK. Ces donnees ne remplacent pas la CMI du patient.", "Select the antimicrobial corresponding to the PK model. These data do not replace the patient's MIC.")),
            actionButton(ns("connect"), t("Consulter EUCAST", "Connect to EUCAST"), icon = icon("globe")),
            selectInput(ns("antibiotic"), t("Antimicrobien EUCAST", "EUCAST antimicrobial"), choices = character()),
            actionButton(ns("fetch"), t("Charger les distributions", "Load distributions"), icon = icon("download")),
            selectInput(ns("species"), t("Espece", "Species"), choices = character()),
            uiOutput(ns("eucast_source")), plotOutput(ns("mic_plot"), height = "280px"), DTOutput(ns("mic_table"))),
          nav_panel(t("Methodes et limites", "Methods and limitations"), infection_methods(t)))))))
}

infection_server <- function(id, analysis_store, report_plot_uri, imported = reactive(NULL), soloc = tempdir(), cache = new.env(), open_tdm = NULL) {
  moduleServer(id, function(input, output, session) {
    t <- function(fr, en) app_t(app_language_from_query(session$clientData$url_search %||% ""), fr, en)
    begin_tdm <- function(spec = NULL) {
      if (is.null(open_tdm)) stop("TDM navigation unavailable.")
      open_tdm(spec)
      updateRadioButtons(session, "source", selected = "tdm")
    }
    pk <- pd_pk_server("pk", soloc, cache, reactive(if (identical(imported()$view, "infection") && length(imported()$models)) imported()$models[[1]] else NULL), analysis_store, begin_tdm, allow_tdm = FALSE)
    observeEvent(input$estimate, tryCatch(begin_tdm(), error = function(e) showNotification(conditionMessage(e), type = "error")))
    simulation <- reactiveVal(NULL); historical <- reactiveVal(NULL); remote <- reactiveVal(NULL); choices <- reactiveVal(character())
    config <- reactive({
      defaults <- infection_config()
      values <- setNames(lapply(names(defaults), function(name) input[[name]] %||% defaults[[name]]), names(defaults))
      values$exposure <- if (identical(input$pk_mode, "builtin")) "iv1" else "pk"
      if (!is.null(input$dose_min)) values$grid <- list(min = input$dose_min, max = input$dose_max, step = input$dose_step,
        intervals = as.numeric(input$intervals), infusion = input$grid_infusion, continuous = isTRUE(input$continuous))
      infection_config(values)
    })
    observeEvent(input$metric, {
      labels <- c(time = t("Temps au-dessus de k x CMI (%)", "Time above k x MIC (%)"), auc = t("AUC24 / CMI visee (h)", "Target AUC24 / MIC (h)"), peak = t("Rapport Cmax / CMI vise", "Target Cmax / MIC ratio"))
      updateNumericInput(session, "target", label = labels[[input$metric]], max = if (input$metric == "time") 100 else 1e6)
    })
    output$tdm_units <- renderUI({
      fits <- successful_fits(analysis_store()$fits %||% list())
      if (!length(fits) || all(names(fits) %in% MODEL_CATALOG$id)) return(NULL)
      selectInput(session$ns("tdm_scale"), t("Unite des concentrations PK", "PK concentration unit"), c("mg/L (= ug/mL)" = "1", "ng/mL (= ug/L)" = ".001", "mg/mL" = "1000"))
    })
    tdm_scale <- reactive({
      fits <- successful_fits(analysis_store()$fits %||% list())
      scales <- unique(vapply(fits, function(f) pd_pk_units(f$id, input$tdm_scale %||% 1)$concentration_scale, numeric(1)))
      if (!length(scales)) return(1)
      if (length(scales) != 1) stop("TDM model averaging requires consistent concentration units.")
      scales[[1]]
    })
    context <- reactive(if (identical(input$pk_mode, "builtin")) pd_builtin_iv_context(input$v, input$cl, soloc, cache) else pk())
    current_pk <- reactive(if (identical(input$source, "pk") && identical(input$pk_mode, "free")) tryCatch(pk(), error = function(e) NULL) else NULL)
    assumption_key <- reactive({
      x <- tryCatch(config(), error = function(e) NULL)
      context <- current_pk()
      model <- if (is.null(context)) NULL else list(id = context$id, route = context$route, adm_cmt = context$adm_cmt,
        concentration = context$concentration, time_unit = context$time_unit, concentration_scale = context$concentration_scale,
        parameters = as.list(mrgsolve::param(context$model)))
      jsonlite::toJSON(list(config = x, model = model, tdm_scale = input$tdm_scale, accept = isTRUE(input$accept)),
        auto_unbox = TRUE, null = "null", digits = 16)
    })
    observeEvent(assumption_key(), {
      if (!is.null(simulation()) && !identical(simulation()$assumption_key, assumption_key())) simulation(NULL)
      historical(NULL)
    }, ignoreInit = TRUE)
    observeEvent(analysis_store(), { simulation(NULL); historical(NULL) }, ignoreInit = TRUE)
    observeEvent(imported(), {
      value <- imported(); if (!identical(value$view, "infection")) return()
      updateSelectInput(session, "pk_mode", selected = if (identical(value$config$exposure, "iv1")) "builtin" else "free")
      for (name in setdiff(names(infection_config()), c("grid", "exposure"))) {
        if (name %in% c("source", "basis")) updateRadioButtons(session, name, selected = value$config[[name]])
        else if (name == "metric") updateSelectInput(session, name, selected = value$config[[name]])
        else updateNumericInput(session, name, value = value$config[[name]])
      }
      grid <- value$config$grid %||% list(min = value$config$dose2, max = value$config$dose2, step = 1, intervals = value$config$interval2, infusion = value$config$infusion2)
      for (key in c("min", "max", "step")) updateNumericInput(session, paste0("dose_", key), value = grid[[key]])
      intervals <- as.numeric(unlist(grid$intervals))
      updateSelectizeInput(session, "intervals", choices = sort(unique(c(4, 6, 8, 12, 24, 48, intervals))), selected = as.character(intervals))
      updateNumericInput(session, "grid_infusion", value = grid$infusion)
      updateCheckboxInput(session, "continuous", value = isTRUE(grid$continuous))
      updateCheckboxInput(session, "accept", value = FALSE)
      session$sendCustomMessage("workbench-ack", list(id = value$id, ok = TRUE))
    })
    output$tdm_source <- renderUI({
      result <- analysis_store()
      p(if (is.null(result)) t("Aucune analyse TDM dans cette session.", "No TDM analysis in this session.") else
        paste(paste(vapply(successful_fits(result$fits), `[[`, character(1), "label"), collapse = " / "), result$route))
    })
    observeEvent(input$last_regimen, tryCatch({
      result <- analysis_store(); if (is.null(result)) stop(t("Lancez d'abord une analyse TDM.", "Run a TDM analysis first."))
      reg <- result$current_exposure
      for (name in c("dose", "interval", "infusion")) updateNumericInput(session, name, value = reg[[name]])
    }, error = function(e) showNotification(conditionMessage(e), type = "error")))
    observeEvent(input$simulate, tryCatch({
      simulation(NULL)
      if (!isTRUE(input$accept)) stop(t("Verifiez les hypotheses et confirmez l'usage de recherche.", "Check the assumptions and confirm research use."))
      x <- config()
      result <- withProgress(message = "PTA / Monte Carlo", value = .2, {
        if (x$source == "tdm") {
          previous <- analysis_store(); if (is.null(previous)) stop(t("Aucune analyse TDM disponible.", "No TDM analysis available."))
          fits <- previous$fits
          scale <- tdm_scale()
          fits <- lapply(fits, function(f) { f$concentration_scale <- scale; f })
          infection_simulate(fits, previous$weights, x, previous$route)
        } else {
          selected <- context(); fit <- infection_population_fit(selected, soloc, cache)
          infection_simulate(setNames(list(fit), fit$id), setNames(1, fit$id), x, selected$route)
        }
      })
      result$assumption_key <- isolate(assumption_key())
      simulation(result)
    }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = 15)))
    observeEvent(input$interpret, tryCatch({
      historical(NULL)
      if (!isTRUE(input$accept) || !identical(input$source, "tdm")) stop(t("Selectionnez une analyse TDM et confirmez les hypotheses.", "Select a TDM analysis and confirm assumptions."))
      previous <- analysis_store(); if (is.null(previous)) stop(t("Aucune analyse TDM disponible.", "No TDM analysis available."))
      historical(infection_tdm_summary(previous, config(), tdm_scale()))
    }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = 15)))
    output$status <- renderUI({
      x <- simulation()
      if (is.null(x)) return(p(t("Aucun resultat PTA actuel. Toute modification des hypotheses invalide le calcul precedent.", "No current PTA result. Changing assumptions invalidates the previous calculation.")))
      p(strong(if (x$posterior) t("Probabilite posterieure d'atteinte de cible", "Posterior probability of target attainment") else t("PTA populationnelle", "Population PTA")),
        paste(" | ", paste(x$models, collapse = " / "), " | n =", x$config$replicates, " | ", x$route, " | SS"))
    })
    regimen_label <- function(g) paste0(g$dose, " / ", g$interval, " h; ", t("perfusion", "infusion"), " ", g$infusion, " h")
    output$comparison <- renderUI({
      x <- simulation(); shiny::req(x)
      grid <- infection_rank(x); grid <- grid[grid$regimen != x$regimens$regimen[1], ]
      selectInput(session$ns("compare_regimen"), t("Posologie comparee", "Compared regimen"), setNames(grid$regimen, regimen_label(grid)), grid$regimen[1])
    })
    compared <- reactive({
      x <- simulation(); shiny::req(x)
      available <- x$regimens$regimen[-1]
      selected <- input$compare_regimen
      if (is.null(selected) || !selected %in% available) selected <- infection_rank(x)$regimen[infection_rank(x)$regimen %in% available][1]
      c(x$regimens$regimen[1], selected)
    })
    compared_data <- function(data) {
      keys <- compared(); data <- data[data$regimen %in% keys, ]
      data$regimen <- factor(data$regimen, levels = keys, labels = c(t("Maintenir", "Maintain"), t("Comparer", "Compare")))
      data
    }
    observeEvent(list(simulation(), input$compare_regimen), {
      x <- simulation(); value <- imported()
      if (!identical(value$view, "infection")) return()
      if (is.null(x)) {
        session$sendCustomMessage("workbench-result", list(id = value$id, view = "infection", invalidated = TRUE))
        return()
      }
      keys <- compared()
      curves <- list()
      for (kind in c("exposure", "pta")) for (i in 1:2) {
        data <- if (kind == "exposure") x$exposure else x$curves
        data <- data[data$regimen == keys[i], ]
        data <- data[unique(round(seq(1, nrow(data), length.out = min(1500, nrow(data))))), ]
        xx <- if (kind == "exposure") data$time else data$mic
        yy <- if (kind == "exposure") data$median * if (x$config$basis == "free") x$config$fu else 1 else data$pta
        curves[[length(curves) + 1L]] <- list(key = paste(kind, c("current", "compare")[i], sep = "_"),
          points = unname(lapply(seq_along(xx), function(j) list(x = xx[j], y = yy[j]))))
      }
      session$sendCustomMessage("workbench-result", list(id = value$id, view = "infection", curves = curves))
    }, ignoreNULL = TRUE)
    palette <- c("#006d77", "#b63563")
    pta_plot <- reactive({
      x <- simulation(); shiny::req(x)
      ggplot(compared_data(x$curves), aes(mic, pta, color = regimen, fill = regimen)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .1, color = NA) + geom_line(aes(linetype = regimen), linewidth = 1) + scale_linetype_manual(values = c("dashed", "solid"), name = NULL) +
        geom_hline(yintercept = x$config$pta_target, linetype = "dotted", color = "#444444") +
        geom_vline(xintercept = x$config$mic, linetype = 2) + scale_x_log10(breaks = c(.01, .1, 1, 10, 100), labels = function(v) format(v, scientific = FALSE, trim = TRUE)) +
        scale_y_continuous(limits = c(0, 100)) + scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
        labs(x = t("CMI (mg/L)", "MIC (mg/L)"), y = "PTA (%)", color = NULL, fill = NULL) + theme_minimal(base_size = 12) + theme(legend.position = "bottom")
    })
    output$pta <- renderPlot(pta_plot())
    output$sensitivity <- renderDT({ x <- simulation(); shiny::req(x)
      d <- compared_data(x$curves[x$curves$mic %in% (x$config$mic * c(.5, 1, 2)), ])
      datatable(d, colnames = c(t("Posologie", "Regimen"), t("CMI (mg/L)", "MIC (mg/L)"), "PTA (%)", t("Borne basse MC", "MC lower bound"), t("Borne haute MC", "MC upper bound"), "n"),
        rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE)) |> formatRound(c("pta", "lower", "upper"), 1)
    })
    output$regimens <- renderDT({
      x <- simulation(); shiny::req(x); data <- infection_rank(x)
      data$regimen <- ifelse(data$regimen == x$regimens$regimen[1], t("Actuelle", "Current"), t("Candidate", "Candidate"))
      data$target_met <- ifelse(data$target_met, t("Oui", "Yes"), t("Non", "No"))
      datatable(data, colnames = c(t("Posologie", "Regimen"), "Dose", t("Intervalle (h)", "Interval (h)"), t("Perfusion (h)", "Infusion (h)"), t("Dose / jour", "Daily dose"), "PTA (%)", t("Seuil PTA atteint", "PTA threshold met")),
        rownames = FALSE, options = list(dom = "t", paging = FALSE, scrollX = TRUE, order = list())) |> formatRound("pta", 1)
    })
    output$exposure <- renderPlot({
      x <- simulation(); shiny::req(x); d <- compared_data(x$exposure)
      if (x$config$basis == "free") d[c("lower", "median", "upper")] <- d[c("lower", "median", "upper")] * x$config$fu
      ggplot(d, aes(time, median, color = regimen, fill = regimen)) + geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .12, color = NA) +
        geom_line(aes(linetype = regimen), linewidth = .8) + scale_linetype_manual(values = c("dashed", "solid"), name = NULL) + scale_color_manual(values = palette) + scale_fill_manual(values = palette) +
        geom_hline(yintercept = x$config$mic * if (x$config$metric == "time") x$config$multiple else 1, linetype = "dotted") +
        labs(x = "h", y = if (x$config$basis == "free") "C free (mg/L)" else "C total (mg/L)",
          subtitle = t("Mediane et percentiles 5-95, stationnaire", "Median and 5-95 percentiles, steady state"), color = NULL, fill = NULL) + theme_minimal(base_size = 12) + theme(legend.position = "bottom")
    })
    output$history_summary <- renderUI({
      x <- historical(); if (is.null(x)) return(p(t("Aucun historique interprete.", "No interpreted history.")))
      tagList(p(paste(t("Estimation ponctuelle sur la fenetre", "Point estimate over window"), paste(round(x$window, 2), collapse = " - "), "h")),
        p(paste(t("Indice PK/PD", "PK/PD index"), if (is.finite(x$value)) signif(x$value, 5) else t("Non calculable : AUC24 exige 24 h completes.", "Not calculable: AUC24 requires a full 24 h."))),
        p(paste("C(MAP) =", signif(x$concentration, 5), "mg/L")),
        p(t("Ni probabilite d'efficacite clinique, ni validation de posologie.", "Neither a clinical efficacy probability nor regimen validation.")))
    })
    output$history <- renderPlot({ x <- historical(); shiny::req(x)
      d <- x$profile; d <- d[d$time >= x$window[1] & d$time <= x$window[2], ]
      if (x$config$basis == "free") d$concentration <- d$concentration * x$config$fu
      ggplot(d, aes(time, concentration)) + geom_line(color = "#006d77", linewidth = 1) +
        geom_hline(yintercept = x$config$mic * if (x$config$metric == "time") x$config$multiple else 1, linetype = 2) +
        labs(x = "h", y = if (x$config$basis == "free") "C free (mg/L)" else "C total (mg/L)") + theme_minimal(base_size = 12)
    })
    observeEvent(input$connect, tryCatch({
      values <- eucast_choices(); choices(values); updateSelectInput(session, "antibiotic", choices = values)
    }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = 15)))
    observeEvent(input$antibiotic, remote(NULL), ignoreInit = TRUE)
    observeEvent(input$fetch, tryCatch({
      remote(NULL); ids <- choices(); if (!input$antibiotic %in% unname(ids)) stop("Select an EUCAST antimicrobial.")
      value <- withProgress(message = "EUCAST", value = .2, eucast_load(input$antibiotic, names(ids)[match(input$antibiotic, ids)]))
      remote(value); updateSelectInput(session, "species", choices = names(value$rows))
    }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = 15)))
    distribution <- reactive({ value <- remote(); if (is.null(value) || !input$species %in% names(value$rows)) return(NULL); value$rows[[input$species]] })
    output$eucast_source <- renderUI({ x <- distribution(); shiny::req(x); source <- remote()
      tagList(h4(paste(source$label, x$species, sep = " / ")), p(paste("(T)ECOFF:", x$ecoff, "mg/L | n =", x$observations, "|", t("Distributions", "Distributions"), x$distributions)),
        p(t("Reference epidemiologique, pas une CMI patient ni un seuil clinique. Aucun taux local de resistance ne peut en etre deduit.", "Epidemiological reference, not a patient MIC or clinical breakpoint. No local resistance rate can be inferred.")),
        p("EUCAST MIC and Zone diameter distribution website; last accessed ", source$accessed, ". ", tags$a(href = x$url, target = "_blank", rel = "noopener noreferrer", t("Source officielle", "Official source"))))
    })
    output$mic_plot <- renderPlot({ x <- distribution(); shiny::req(x)
      ggplot(x$data, aes(factor(mic, levels = mic), count)) + geom_col(fill = "#765b25") +
        labs(x = t("CMI (mg/L)", "MIC (mg/L)"), y = t("Effectif agrege", "Aggregated count")) + theme_minimal(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    })
    output$mic_table <- renderDT({ x <- distribution(); shiny::req(x); datatable(x$data, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) })
    output$csv <- downloadHandler(filename = function() "infection-pta.csv", content = function(file) { shiny::req(simulation()); write.csv(simulation()$curves, file, row.names = FALSE) })
    output$report <- downloadHandler(filename = function() "infection-report.html", content = function(file) {
      x <- simulation(); shiny::req(x)
      htmltools::save_html(tags$html(tags$head(tags$meta(charset = "utf-8"), tags$title("Infectiology / PTA")), tags$body(
        h1("Infectiology / PTA"), p(x$created), p(paste(x$models, collapse = " / ")), p(paste("Route:", x$route, "| Posterior:", x$posterior)),
        tags$pre(jsonlite::toJSON(list(config = x$config, weights = x$weights, parameters = x$parameters, seed = x$seed), auto_unbox = TRUE, pretty = TRUE)),
        tags$img(src = report_plot_uri(pta_plot()), alt = "PTA versus MIC"), infection_methods(t))), file = file)
    })
    session$onSessionEnded(function() { simulation(NULL); historical(NULL); remote(NULL); choices(character()) })
  })
}
