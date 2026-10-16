if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(mrgsolve)
  library(mapbayr)
  library(dplyr)
  library(ggplot2)
  library(DT)
  library(jsonlite)
})

`%||%` <- function(value, fallback) {
  if (is.null(value) || !length(value)) fallback else value
}

format_metric <- function(value, digits = 1, lang = "fr") {
  if (!length(value) || !is.finite(value[[1]])) return(app_t(lang, "Non calculable", "Not available"))
  format(round(as.numeric(value[[1]]), digits), nsmall = digits, trim = TRUE)
}

target_metric_label <- function(metric, lang = "fr") {
  labels <- if (identical(normalize_app_language(lang), "en")) {
    c(AUC24 = "AUC 0-24 h", Cmin = "Cmin", Cmax = "Cmax", TimeAbove = "% time above threshold", CminCmax = "Cmin + Cmax")
  } else {
    c(AUC24 = "AUC 0-24 h", Cmin = "Cmin", Cmax = "Cmax", TimeAbove = "% temps au-dessus du seuil", CminCmax = "Cmin + Cmax")
  }
  unname(labels[[metric]] %||% metric)
}

target_definition_label <- function(metric, low, high, options = list(), lang = "fr") {
  if (identical(metric, "TimeAbove")) {
    return(app_t(
      lang,
      paste0("%T > ", format_metric(options$concentration_threshold), " : ", format_metric(low), " - ", format_metric(high), " %"),
      paste0("%T > ", format_metric(options$concentration_threshold), ": ", format_metric(low), " - ", format_metric(high), "%")
    ))
  }
  if (identical(metric, "CminCmax")) {
    component_label <- function(name, operator, lower, upper) {
      if (identical(operator, "above")) return(paste(name, ">=", format_metric(lower)))
      if (identical(operator, "below")) return(paste(name, "<=", format_metric(upper)))
      paste(name, format_metric(lower), "-", format_metric(upper))
    }
    return(paste0(
      component_label("Cmin", options$cmin_operator %||% "range", low, high),
      " ; ",
      component_label("Cmax", options$cmax_operator %||% "range", options$cmax_low, options$cmax_high)
    ))
  }
  paste0(target_metric_label(metric, lang), " ", format_metric(low), " - ", format_metric(high))
}

target_component_plot_bounds <- function(operator, low, high) {
  if (identical(operator, "above")) return(c(low, Inf))
  if (identical(operator, "below")) return(c(-Inf, high))
  c(low, high)
}

target_value_label <- function(metric, row, lang = "fr") {
  if (identical(metric, "CminCmax")) {
    return(paste0("Cmin ", format_metric(row$cmin), " / Cmax ", format_metric(row$cmax)))
  }
  if (identical(metric, "TimeAbove")) return(paste0(format_metric(row$time_above_pct), " %"))
  format_metric(row$value)
}

APP_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "private_workspace.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "i18n.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "clinical_presets.R"), local = TRUE)
pmetrics_mode <- tolower(trimws(Sys.getenv("ENABLE_PMETRICS", "auto")))
if (!pmetrics_mode %in% c("auto", "true", "false")) stop("ENABLE_PMETRICS must be auto, true, or false.")
pmetrics_installed <- requireNamespace(paste0("P", "metrics"), quietly = TRUE)
if (identical(pmetrics_mode, "true") && !pmetrics_installed) stop("ENABLE_PMETRICS=true but Pmetrics is not installed.")
PMETRICS_ENABLED <- !identical(pmetrics_mode, "false") && pmetrics_installed
source(file.path(APP_ROOT, "R", if (PMETRICS_ENABLED) "pmetrics_engine.R" else "pmetrics_disabled.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "engine.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "ml_engine.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "ddi_engine.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "ddi_builder.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "pd_engine.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "pd_pk.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "pd_observations.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "pd_module.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "onco_engine.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "onco_module.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "infection_engine.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "eucast.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "infection_module.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "workbench_bridge.R"), local = TRUE)
source(file.path(APP_ROOT, "R", "teaching_lab.R"), local = TRUE)

ml_concordance <- function(result, lang = "fr") {
  map_auc24 <- suppressWarnings(as.numeric(result$current_exposure$steady_state_auc24 %||% NA_real_))
  ml_auc24 <- suppressWarnings(as.numeric(result$ml_status$auc24 %||% NA_real_))
  if (!is.finite(map_auc24) || map_auc24 <= 0 || !is.finite(ml_auc24) || ml_auc24 <= 0) return(NULL)
  relative_gap <- 100 * (ml_auc24 / map_auc24 - 1)
  absolute_gap <- abs(relative_gap)
  level <- if (absolute_gap <= 10) "close" else if (absolute_gap <= 20) "moderate" else "divergent"
  interpretation <- switch(level,
    close = app_t(lang,
      paste0("Les estimations sont proches (écart relatif ", format_metric(absolute_gap), " %). Cette concordance soutient la cohérence du résultat pour ce profil, sans constituer une validation clinique du ML."),
      paste0("The estimates are close (relative difference ", format_metric(absolute_gap), "%). This agreement supports consistency for this profile but is not a clinical validation of ML.")),
    moderate = app_t(lang,
      paste0("L'écart entre MAP-BE et ML est modéré (", format_metric(absolute_gap), " %). Vérifiez en priorité les horaires, les doses et les covariables avant d'interpréter l'estimation ML."),
      paste0("The MAP-BE and ML estimates differ moderately (", format_metric(absolute_gap), "%). Check times, doses, and covariates before interpreting the ML estimate.")),
    divergent = app_t(lang,
      paste0("Les deux méthodes divergent de ", format_metric(absolute_gap), " %. Cette discordance doit conduire à contrôler les données et à privilégier l'analyse pharmacométrique validée; elle ne permet pas de choisir la valeur la plus favorable."),
      paste0("The two methods differ by ", format_metric(absolute_gap), "%. Check the data and prioritize the validated pharmacometric analysis; do not select the more favorable value."))
  )
  list(
    map_auc24 = map_auc24,
    ml_auc24 = ml_auc24,
    unit = result$ml_status$unit %||% "mg.h/L",
    relative_gap = relative_gap,
    absolute_gap = absolute_gap,
    level = level,
    interpretation = interpretation
  )
}

build_ml_comparison_plot <- function(result, lang = "fr") {
  concordance <- ml_concordance(result, lang)
  if (is.null(concordance)) return(NULL)
  data <- data.frame(
    method = factor(c("MAP-BE", "ML"), levels = c("MAP-BE", "ML")),
    auc24 = c(concordance$map_auc24, concordance$ml_auc24)
  )
  color <- switch(concordance$level, close = "#287a4d", moderate = "#b16916", divergent = "#a43d35")
  limit_values <- data$auc24
  if (identical(result$target_metric, "AUC24")) {
    limit_values <- c(limit_values, result$target_low, result$target_high)
  }
  limit_span <- diff(range(limit_values, na.rm = TRUE))
  limit_padding <- max(25, 0.18 * if (is.finite(limit_span) && limit_span > 0) limit_span else max(limit_values))
  y_limits <- c(max(0, min(limit_values) - limit_padding), max(limit_values) + limit_padding)
  plot <- ggplot(data, aes(method, auc24, group = 1))
  if (identical(result$target_metric, "AUC24")) {
    plot <- plot + annotate(
      "rect",
      xmin = -Inf,
      xmax = Inf,
      ymin = result$target_low,
      ymax = result$target_high,
      fill = "#dcece5",
      alpha = 0.7
    )
  }
  plot +
    geom_line(color = color, linewidth = 1.2) +
    geom_point(aes(fill = method), shape = 21, size = 6, color = "white", stroke = 1.1) +
    geom_text(aes(label = sprintf("%.1f", auc24)), vjust = -1.15, fontface = "bold", color = "#17202a") +
    scale_fill_manual(values = c("MAP-BE" = "#176b70", "ML" = "#a4441f"), guide = "none") +
    coord_cartesian(ylim = y_limits, clip = "off") +
    labs(
      x = NULL,
      y = app_t(lang, paste0("AUC0-24 à l'état stationnaire (", concordance$unit, ")"), paste0("Steady-state AUC0-24 (", concordance$unit, ")")),
      subtitle = if (identical(result$target_metric, "AUC24")) app_t(lang, "La zone verte représente la cible définie dans les réglages.", "The green area represents the target defined in settings.") else NULL,
      caption = paste(strwrap(concordance$interpretation, width = 105), collapse = "\n")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(face = "bold", color = "#17202a"),
      plot.caption = element_text(hjust = 0, color = "#43515c", size = 10)
    )
}

ML_FEATURE_LABELS <- c(
  WT = "Poids",
  AGE = "Âge",
  CREAT = "Créatinine",
  SEX = "Sexe",
  DOSE = "Dose",
  INTERVAL = "Intervalle d'administration",
  INFUSION = "Durée de perfusion",
  POP_AUC24 = "AUC24 populationnelle",
  ML_AUC24 = "AUC24 prédite par ML",
  PREV_CONC = "Première concentration récente",
  PREV_TIME = "Horaire de la première concentration",
  LAST_CONC = "Dernière concentration",
  LAST_TIME = "Horaire de la dernière concentration",
  CONC_DIFF = "Variation entre les concentrations",
  TIME_DIFF = "Délai entre les concentrations",
  PREV_POP_CONC = "Première concentration populationnelle",
  LAST_POP_CONC = "Dernière concentration populationnelle",
  PREV_CONC_RATIO = "Rapport observé/prédit initial",
  LAST_CONC_RATIO = "Rapport observé/prédit final"
)

ML_FEATURE_LABELS_EN <- c(
  WT = "Weight", AGE = "Age", CREAT = "Creatinine", SEX = "Sex", DOSE = "Dose",
  INTERVAL = "Dosing interval", INFUSION = "Infusion duration",
  POP_AUC24 = "Population AUC24",
  ML_AUC24 = "ML-predicted AUC24",
  PREV_CONC = "First recent concentration", PREV_TIME = "First concentration time",
  LAST_CONC = "Last concentration", LAST_TIME = "Last concentration time",
  CONC_DIFF = "Concentration difference", TIME_DIFF = "Time between concentrations",
  PREV_POP_CONC = "First population concentration", LAST_POP_CONC = "Last population concentration",
  PREV_CONC_RATIO = "Initial observed/predicted ratio", LAST_CONC_RATIO = "Final observed/predicted ratio"
)

ml_feature_label <- function(name, lang = "fr") {
  labels <- if (identical(normalize_app_language(lang), "en")) ML_FEATURE_LABELS_EN else ML_FEATURE_LABELS
  label <- unname(labels[name])
  if (!length(label) || is.na(label) || !nzchar(label)) name else label
}

format_ml_feature_value <- function(name, value, lang = "fr", auc_unit = "mg.h/L") {
  value <- suppressWarnings(as.numeric(value))
  if (!is.finite(value)) return(app_t(lang, "valeur indisponible", "unavailable"))
  if (identical(name, "SEX")) return(if (value >= 0.5) app_t(lang, "femme (1)", "female (1)") else app_t(lang, "homme (0)", "male (0)"))
  concentration_unit <- sub("\\.h", "", auc_unit)
  unit <- switch(
    name,
    WT = " kg",
    AGE = app_t(lang, " ans", " years"),
    CREAT = " micromol/L",
    DOSE = " mg",
    INTERVAL = " h",
    INFUSION = " h",
    POP_AUC24 = paste0(" ", auc_unit),
    ML_AUC24 = paste0(" ", auc_unit),
    PREV_CONC = paste0(" ", concentration_unit),
    LAST_CONC = paste0(" ", concentration_unit),
    PREV_POP_CONC = paste0(" ", concentration_unit),
    LAST_POP_CONC = paste0(" ", concentration_unit),
    CONC_DIFF = paste0(" ", concentration_unit),
    PREV_TIME = app_t(lang, " h après la dose", " h after dose"),
    LAST_TIME = app_t(lang, " h après la dose", " h after dose"),
    TIME_DIFF = " h",
    ""
  )
  paste0(format(signif(value, 4), trim = TRUE), unit)
}

format_ml_domain_warning <- function(item, lang = "fr", auc_unit = "mg.h/L") {
  name <- as.character(item$feature %||% "")
  label <- ml_feature_label(name, lang)
  lower <- if (is.finite(item$min %||% NA_real_)) format_ml_feature_value(name, item$min, lang, auc_unit) else app_t(lang, "sans borne basse", "no lower bound")
  upper <- if (is.finite(item$max %||% NA_real_)) format_ml_feature_value(name, item$max, lang, auc_unit) else app_t(lang, "sans borne haute", "no upper bound")
  app_t(lang,
    paste0(label, " : ", format_ml_feature_value(name, item$value, lang, auc_unit), " ; plage d'entraînement [", lower, " ; ", upper, "]."),
    paste0(label, ": ", format_ml_feature_value(name, item$value, lang, auc_unit), "; training range [", lower, "; ", upper, "]."))
}

ml_domain_warning_block <- function(status, lang = "fr") {
  warnings <- status$domain_warnings %||% list()
  if (!length(warnings)) return(NULL)
  div(
    class = "ml-domain-warning",
    tags$strong(app_t(lang, "Avertissement : prédiction ML extrapolée", "Warning: extrapolated ML prediction")),
    tags$ul(lapply(warnings, function(item) tags$li(format_ml_domain_warning(item, lang, status$unit %||% "mg.h/L")))),
    span(app_t(lang, "Le résultat ML est affiché à titre expérimental, mais sa fiabilité peut être réduite. La recommandation reste calculée par MAP-BE.", "The ML result is shown for experimental use, but reliability may be reduced. The recommendation remains MAP-BE based."))
  )
}

build_ml_explanation_plot <- function(explanation, max_features = 8L, lang = "fr", unit = "mg.h/L") {
  if (!isTRUE((explanation %||% list())$available)) return(NULL)
  contributions <- as.data.frame(explanation$contributions, stringsAsFactors = FALSE)
  if (!nrow(contributions)) return(NULL)
  contributions <- contributions[order(-abs(contributions$contribution)), , drop = FALSE]
  labels <- vapply(seq_len(nrow(contributions)), function(index) {
    name <- contributions$variable[[index]]
    label <- ml_feature_label(name, lang)
    paste0(
      label,
      "\n",
      format_ml_feature_value(name, contributions$value[[index]], lang, unit)
    )
  }, character(1))
  effects <- contributions$contribution
  if (length(effects) > max_features) {
    labels <- c(labels[seq_len(max_features)], app_t(lang, "Autres informations\ncombinées", "Other combined\ninformation"))
    effects <- c(effects[seq_len(max_features)], sum(effects[(max_features + 1L):length(effects)]))
  }
  data <- data.frame(
    label = factor(labels, levels = rev(labels)),
    effect = effects,
    direction = ifelse(effects >= 0, "increase", "decrease")
  )
  range <- max(abs(data$effect), na.rm = TRUE)
  if (!is.finite(range) || range <= 0) range <- 1
  ggplot(data, aes(effect, label, fill = direction)) +
    geom_col(width = 0.66) +
    geom_vline(xintercept = 0, color = "#17202a", linewidth = 0.5) +
    geom_text(
      aes(label = sprintf("%+.1f", effect), hjust = ifelse(effect >= 0, -0.15, 1.15)),
      color = "#17202a",
      size = 3.4
    ) +
    scale_fill_manual(values = c(increase = "#287a4d", decrease = "#a4441f"), guide = "none") +
    scale_x_continuous(limits = c(-1.45 * range, 1.45 * range)) +
    labs(
      x = app_t(lang, paste0("Contribution à l'AUC24 ML (", unit, ")"), paste0("Contribution to ML AUC24 (", unit, ")")),
      y = NULL,
      subtitle = app_t(lang, "À droite, l'information augmente l'AUC prédite; à gauche, elle la diminue.", "Information to the right increases predicted AUC; information to the left decreases it."),
      caption = app_t(lang,
        paste0("Référence DALEX synthétique ", format_metric(explanation$baseline), " + contributions = AUC24 ML ", format_metric(explanation$prediction), " ", unit, " · référence fondée sur ", explanation$background_size, " profils simulés."),
        paste0("Synthetic DALEX baseline ", format_metric(explanation$baseline), " + contributions = ML AUC24 ", format_metric(explanation$prediction), " ", unit, " · baseline from ", explanation$background_size, " simulated profiles."))
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 9.5, lineheight = 1.05, color = "#17202a"),
      plot.subtitle = element_text(size = 10, color = "#43515c"),
      plot.caption = element_text(size = 10, color = "#43515c", hjust = 0)
    )
}

ml_status_message <- function(status, lang = "fr") {
  message <- as.character(status$message %||% "")
  if (!identical(normalize_app_language(lang), "en")) return(message)
  auc24 <- suppressWarnings(as.numeric(status$auc24 %||% NA_real_))
  if (is.finite(auc24)) {
    unit <- status$unit %||% "mg.h/L"
    suffix <- if (length(status$artifacts %||% character()) > 1L) " and model averaging." else "."
    warning <- if (length(status$domain_warnings %||% list())) paste0(" Warning: extrapolation outside the training domain for ", length(status$domain_warnings), " feature(s).") else ""
    validation_warning <- if (isFALSE(status$research %||% TRUE)) " One or more artifacts remain below the internal performance thresholds and are shown for experimental use only." else ""
    return(paste0("Experimental ML AUC24: ", round(auc24, 1), " ", unit, ". Dose projections remain MAP based", suffix, warning, validation_warning))
  }
  if (grepl("prédicteur\\(s\\).*disponible", message)) return(paste(status$available %||% 0L, "evaluated experimental direct AUC24 predictor(s) available."))
  if (grepl("non activée", message, fixed = TRUE)) return("Experimental ML AUC24 not enabled: MAP estimate retained.")
  if (grepl("Artefacts ML incomplets", message, fixed = TRUE)) return("Incomplete ML artifacts for the selected models: aggregate estimate not shown.")
  if (grepl("modèle de session", message, fixed = TRUE)) return("ML: unavailable for a session model.")
  if (grepl("Aucun prédicteur", message, fixed = TRUE)) return("No evaluated compatible direct AUC24 predictor: MAP estimate retained.")
  if (grepl("AUC24 ML non appliquée", message, fixed = TRUE)) {
    reason <- sub("^AUC24 ML non appliquée : ", "", message)
    reason <- sub("\\. Estimation MAP conservée\\.$", "", reason)
    reason <- gsub("le prédicteur AUC24 ML exige au moins", "the ML AUC24 predictor requires at least", reason, fixed = TRUE)
    reason <- gsub("concentrations dans un même intervalle posologique", "concentrations in the same dosing interval", reason, fixed = TRUE)
    reason <- gsub("le prédicteur AUC24 ML exige un C0 pré-dose suivi d'une concentration post-dose dans l'intervalle suivant", "the ML AUC24 predictor requires a pre-dose C0 followed by a post-dose concentration in the next dosing interval", reason, fixed = TRUE)
    reason <- gsub("une administration déclarée à l'état stationnaire", "an administration declared at steady state", reason, fixed = TRUE)
    return(paste0("ML AUC24 not applied: ", reason, ". MAP estimate retained."))
  }
  message
}

ALLOW_CUSTOM_MODELS <- identical(tolower(Sys.getenv("ALLOW_CUSTOM_MODELS", "false")), "true")
ANALYSIS_MODEL_IDS <- MODEL_CATALOG$id[vapply(seq_len(nrow(MODEL_CATALOG)), function(index) {
  model_analysis_eligible(MODEL_CATALOG[index, , drop = FALSE])
}, logical(1))]
DEFAULT_MODEL <- if ("vanco_roberts" %in% ANALYSIS_MODEL_IDS) "vanco_roberts" else ANALYSIS_MODEL_IDS[[1]]
DEFAULT_CODE <- read_library_code(DEFAULT_MODEL)
DEFAULT_ROUTE <- model_routes(model_record(DEFAULT_MODEL))[[1]]
DEFAULT_MODE <- model_administration_modes(model_record(DEFAULT_MODEL), DEFAULT_ROUTE)[[1]]
DEFAULT_TARGET_PRESET <- clinical_target_preset(model_record(DEFAULT_MODEL), DEFAULT_MODE)
PMETRICS_CATALOG <- pmetrics_artifact_catalog()
DEFAULT_PMETRICS_ARTIFACT <- if (nrow(PMETRICS_CATALOG)) PMETRICS_CATALOG$id[[1]] else ""
STANDARD_CANDIDATE_INTERVALS <- c(6, 8, 12, 24, 48)

APP_THEME <- bs_theme(
  version = 5,
  bg = "#f7f8fa",
  fg = "#17202a",
  primary = "#176b70",
  secondary = "#5c6670",
  success = "#287a4d",
  warning = "#b16916",
  danger = "#a43d35"
)

time_input_ui <- function(prefix, index, time = 0, days_ago = 0, clock = "08:00", calendar = "", lang = "fr") {
  localize_ui(tagList(
    conditionalPanel(
      condition = "input.time_entry_mode == 'relative_hours'",
      numericInput(paste0(prefix, "_time_", index), "Temps relatif (h)", time, step = 0.25)
    ),
    conditionalPanel(
      condition = "input.time_entry_mode == 'relative_days'",
      div(
        class = "compact-time-grid",
        numericInput(paste0(prefix, "_days_ago_", index), "Jours avant référence", days_ago, min = 0, step = 1),
        textInput(paste0(prefix, "_clock_", index), "Heure (HH:MM)", clock)
      )
    ),
    conditionalPanel(
      condition = "input.time_entry_mode == 'calendar'",
      textInput(
        paste0(prefix, "_calendar_", index),
        "Date/heure",
        calendar,
        placeholder = "JJ/MM HH:MM ou JJ/MM/AAAA HH:MM"
      )
    )
  ), lang)
}

dose_row_ui <- function(
  index,
  time = 0,
  amount = 1000,
  interval = 12,
  count = 4,
  infusion = 1,
  steady_state = FALSE,
  status = "administered",
  time_uncertainty = 0,
  lang = "fr"
) {
  localize_ui(div(
    id = paste0("dose_row_", index),
    class = "record-row dose-row",
    div(class = "record-index", paste0("A", index)),
    conditionalPanel(
      condition = sprintf("input.dose_ss_%d != true", index),
      div(class = "time-entry", time_input_ui("dose", index, time = time, lang = lang))
    ),
    conditionalPanel(
      condition = sprintf("input.dose_ss_%d == true", index),
      div(class = "route-readonly", tags$strong("Steady state"), span("Temps imposé à t = 0 h (ss = 1)."))
    ),
    numericInput(paste0("dose_amount_", index), "Dose (mg)", amount, min = 0, step = 50),
    selectInput(
      paste0("dose_status_", index),
      "Statut",
      choices = stats::setNames(
        c("administered", "uncertain", "missed"),
        if (identical(normalize_app_language(lang), "en")) c("Administered", "Uncertain time", "Missed / not taken") else c("Administr\u00e9e", "Horaire incertain", "Oubli\u00e9e / non prise")
      ),
      selected = status
    ),
    conditionalPanel(
      condition = sprintf("input.dose_status_%d == 'uncertain'", index),
      numericInput(paste0("dose_time_uncertainty_", index), "Incertitude horaire (± h)", time_uncertainty, min = 0, max = 24, step = 0.25)
    ),
    numericInput(paste0("dose_interval_", index), "Intervalle (h)", interval, min = 0, step = 1),
    div(
      class = "dose-repeat-control",
      checkboxInput(paste0("dose_ss_", index), "Steady state (ss = 1)", value = steady_state),
      conditionalPanel(
        condition = sprintf("input.dose_ss_%d != true", index),
        numericInput(paste0("dose_count_", index), "Nombre", count, min = 1, step = 1)
      )
    ),
    conditionalPanel(
      condition = "input.administration_route == 'IV'",
      numericInput(paste0("dose_infusion_", index), "Perfusion (h, 0 = bolus IV)", infusion, min = 0, step = 0.25)
    ),
    conditionalPanel(
      condition = "input.administration_route == 'Oral'",
      div(class = "route-readonly", tags$strong("Voie orale"), span("La perfusion est imposée à 0 h."))
    ),
    conditionalPanel(
      condition = "input.administration_route == 'IM'",
      div(class = "route-readonly", tags$strong("Voie intramusculaire"), span("La perfusion externe est imposée à 0 h; l'absorption est définie par le modèle."))
    ),
    div(
      class = "observation-covariates dose-covariates",
      uiOutput(paste0("dose_covariates_", index))
    ),
    div(class = "row-validation", uiOutput(paste0("dose_validation_", index)))
  ), lang)
}

observation_row_ui <- function(
  index,
  time = 47.5,
  concentration = 18,
  blq = FALSE,
  lloq = 0,
  matrix = "unspecified",
  time_uncertainty = 0,
  lang = "fr"
) {
  localize_ui(div(
    id = paste0("observation_row_", index),
    class = "record-row observation-row",
    div(class = "record-index", paste0("P", index)),
    div(class = "time-entry", time_input_ui("observation", index, time = time, lang = lang)),
    numericInput(paste0("observation_concentration_", index), "Concentration", concentration, min = 0, step = 0.1),
    checkboxInput(paste0("observation_blq_", index), "Sous la limite de quantification (BLQ)", blq),
    conditionalPanel(
      condition = sprintf("input.observation_blq_%d == true", index),
      numericInput(paste0("observation_lloq_", index), "LLOQ", lloq, min = 0, step = 0.1)
    ),
    selectInput(
      paste0("observation_matrix_", index),
      "Matrice biologique",
      choices = stats::setNames(
        c("unspecified", "plasma", "serum", "whole_blood"),
        if (identical(normalize_app_language(lang), "en")) c("Unspecified", "Plasma", "Serum", "Whole blood") else c("Non pr\u00e9cis\u00e9e", "Plasma", "S\u00e9rum", "Sang total")
      ),
      selected = matrix
    ),
    checkboxInput(paste0("observation_time_uncertain_", index), "Horaire incertain", is.finite(time_uncertainty) && time_uncertainty > 0),
    conditionalPanel(
      condition = sprintf("input.observation_time_uncertain_%d == true", index),
      numericInput(paste0("observation_time_uncertainty_", index), "Incertitude horaire (± h)", time_uncertainty, min = 0, max = 24, step = 0.25)
    ),
    div(
      class = "observation-covariates",
      uiOutput(paste0("observation_covariates_", index))
    ),
    div(class = "row-validation", uiOutput(paste0("observation_validation_", index)))
  ), lang)
}

ddi_panel <- function(lang = "fr") {
  t <- function(fr, en) app_t(lang, fr, en)
  choices <- catalog_choices_i18n(lang = lang)
  model_panel <- function(side, default, dose) {
    id <- function(name) paste0("ddi_", name, "_", side)
    accordion_panel(
      if (side == 1) t("1 - Molécule affectée", "1 - Affected drug") else t("3 - Molécule interagissante", "3 - Interacting drug"),
      value = paste0("ddi_model_", side, "_panel"),
      ddi_paste_ui(side, t),
      conditionalPanel(paste0("input.", id("source"), " != 'code'"),
        selectInput(id("model"), t("Modele", "Model"), choices, default), uiOutput(paste0("ddi_route_", side, "_ui"))),
      fluidRow(
        column(4, numericInput(id("dose"), t("Dose", "Dose"), dose, min = 0.001)),
        column(4, numericInput(id("interval"), t("Intervalle (h)", "Interval (h)"), 12, min = 0.25)),
        column(4, numericInput(id("infusion"), t("Perfusion (h)", "Infusion (h)"), 0, min = 0))),
      tags$small(class = "form-text", t("Perfusion = 0 pour oral ou bolus IV. Dose et concentration dans les unites du modele.", "Infusion = 0 for oral or IV bolus. Dose and concentration use model units.")),
      div(class = "ddi-patient-actions",
        actionButton(id("use_tdm"), t("Reprendre TDM / Lego", "Use TDM / Lego fit"), icon = icon("file-import"), class = "btn-outline-primary btn-sm"),
        actionButton(id("use_population"), t("Valeurs populationnelles", "Population values"), icon = icon("users"), class = "btn-outline-secondary btn-sm")),
      uiOutput(id("source_note")),
      if (side == 2) fluidRow(
        column(4, numericInput("ddi_start_day", t("Debut (jour)", "Start (day)"), 3, min = 0.25)),
        column(4, numericInput("ddi_stop_day", t("Arret (jour)", "Stop (day)"), 10, min = 0.5)),
        column(4, numericInput("ddi_followup", t("Suivi (jours)", "Follow-up (days)"), 3, min = 0))))
  }
  nav_panel(t("Interactions", "Interactions"), value = "ddi",
    div(class = "ddi-shell",
      fileInput("ddi_workshop_file", t("Importer un atelier (.json)", "Import workshop (.json)"), accept = ".json"),
      div(class = "safety-banner",
        div(tags$strong(t("Prototype de recherche et d'enseignement", "Research and teaching prototype")),
          span(t("Relation definie par l'utilisateur, sans validation clinique. Session uniquement.", "User-defined relationship, not clinically validated. Session only."))),
        span(class = "engine-badge", "mrgsolve / DDI")),
      tags$button(type = "button", class = "mobile-configure-button ddi-mobile-configure-button",
        onclick = "var r=this.closest('.ddi-shell');var b=r.querySelector('.collapse-toggle');if(b&&b.getAttribute('aria-expanded')!=='true')b.click();",
        t("Configurer l'interaction", "Configure interaction")),
      layout_sidebar(sidebar = sidebar(width = 410, open = "desktop",
        div(class = "sidebar-heading", h2(t("Atelier DDI", "DDI builder")), uiOutput("ddi_status")),
        accordion(id = "ddi_settings", open = c("ddi_model_1_panel", "ddi_interaction_panel"),
          model_panel(1, "tacrolimus_woillard_ddi", 3),
          accordion_panel(t("2 - Interaction", "2 - Interaction"), value = "ddi_interaction_panel",
            uiOutput("ddi_target_ui"),
            ddi_mechanism_ui(t)),
          model_panel(2, "voriconazole_vandenborn_ddi", 200)),
        div(class = "ddi-actions",
          actionButton("assemble_ddi", t("Generer le modele", "Generate model"), icon = icon("code"), class = "btn-outline-primary w-100"),
          checkboxInput("ddi_accept_disclaimer", t("Usage exploratoire, sous ma responsabilite", "Exploratory use, under my responsibility"), FALSE),
          actionButton("run_ddi", t("Simuler l'interaction", "Simulate interaction"), icon = icon("play"), class = "btn-primary w-100"))),
        div(class = "ddi-workspace",
          uiOutput("ddi_chain"),
          navset_tab(id = "ddi_tabs",
            nav_panel(t("Modele assemble", "Assembled model"), value = "model",
              uiOutput("ddi_code_status"), verbatimTextOutput("ddi_generated_code"),
              conditionalPanel("output.ddi_code_ready",
                div(class = "workbench-downloads",
                  downloadButton("download_ddi_code", t("Modele 1 + interaction (.cpp)", "Model 1 + interaction (.cpp)")),
                  downloadButton("download_ddi_r", t("Simulation couplee (.R)", "Coupled simulation (.R)"))))),
            nav_panel(t("Simulation", "Simulation"), value = "simulation",
              uiOutput("ddi_empty"), conditionalPanel("output.ddi_ready",
                uiOutput("ddi_metrics"), h3(t("Molecule affectee", "Affected drug")),
                plotOutput("ddi_concentration_plot", height = "420px"),
                h3(t("Parametre et exposition du modele 2", "Parameter and model 2 exposure")),
                plotOutput("ddi_effect_plot", height = "310px"), DTOutput("ddi_schedule"),
                downloadButton("download_ddi_report", t("Creer le rapport", "Create report")))),
            nav_panel(t("Methode et sources", "Methods and sources"),
              tags$p(t("Le modele 1 commence a l'etat stationnaire et sa posologie est maintenue. Le modele 2 commence sans exposition prealable au jour choisi, puis ses doses sont arretees. Les AUC et extrema comparent les 24 dernieres heures avant cet arret (ou la fenetre disponible si elle est plus courte).",
                "Model 1 starts at steady state and its regimen is maintained. Model 2 starts without prior exposure on the selected day, then its doses are stopped. AUC and extrema compare the final 24 hours before this stop (or the available window if shorter).")),
              tags$p(t("Couplage unidirectionnel sur une grille de 0,1 h. Les mecanismes dynamiques integrent une activite A(0)=1 avec renouvellement kdeg : inhibition dependante du temps (kinact/KI) ou induction de la synthese (Emax/EC50). La recuperation continue apres l'arret. La fraction affectee melange une part non affectee et une part modulee; elle est renseignee par l'utilisateur et ne constitue pas une estimation mecanistique de fm.",
                "Unidirectional coupling on a 0.1 h grid. Dynamic mechanisms integrate activity A(0)=1 with kdeg turnover: time-dependent inhibition (kinact/KI) or synthesis induction (Emax/EC50). Recovery continues after stopping treatment. The affected fraction mixes unaffected and modulated components; it is user supplied and is not a mechanistic fm estimate.")),
              tags$p(t("Pour les mecanismes dynamiques, le C++ attend DDI_ACTIVITY, calculee par le script R couple. Il ne contient pas lui-meme l'ODE enzymatique. Une inhibition enzymatique doit cibler un parametre compatible (p. ex. clairance), pas un volume sans justification mecanistique.",
                "For dynamic mechanisms, C++ expects DDI_ACTIVITY calculated by the coupled R script. It does not itself contain the enzyme ODE. Enzyme inhibition must target a compatible parameter (e.g. clearance), not a volume without mechanistic justification.")),
              tags$a(href = "https://www.fda.gov/regulatory-information/search-fda-guidance-documents/m12-drug-interaction-studies", target = "_blank", rel = "noopener noreferrer", "ICH M12: enzyme inhibition and induction"),
              tags$p(t("Le script R conserve les deux modeles et les parametres selectionnes. Le C++ expose DDI_CP (concentration externe) et DDI_ACTIVE (facteur actif); seul, il ne simule pas la molecule 2. Les valeurs du C++ sont populationnelles; le script R conserve aussi les instantanes TDM.",
                "The R script retains both models and selected parameters. C++ exposes DDI_CP (external concentration) and DDI_ACTIVE (active factor); alone it does not simulate drug 2. C++ defaults are population values; the R script also preserves TDM snapshots.")),
              uiOutput("ddi_method")))))))
}

app_ui <- function(request) {
  lang <- app_language_from_request(request)
  query <- parseQueryString(tryCatch(request$QUERY_STRING, error = function(error) "") %||% "")
  initial_view <- if ((query$view %||% "") %in% c("ddi", "pd")) query$view else "analysis"
  model_choices <- catalog_choices_i18n(lang = lang, analysis_only = TRUE)
  pmetrics_choices <- if (nrow(PMETRICS_CATALOG)) {
    stats::setNames(PMETRICS_CATALOG$id, if (identical(lang, "en")) PMETRICS_CATALOG$labelEn else PMETRICS_CATALOG$label)
  } else character()
  localize_ui(page_navbar(
  title = div(
    class = "brand",
    span(class = "brand-mark", "Pk"),
    span("MIPD Engine"),
    span(
      class = "app-language-switch",
      tags$button(type = "button", class = if (lang == "fr") "active" else NULL, onclick = "var u=new URL(location.href);u.searchParams.set('lang','fr');location.href=u;", "FR"),
      tags$button(type = "button", class = if (lang == "en") "active" else NULL, onclick = "var u=new URL(location.href);u.searchParams.set('lang','en');location.href=u;", "EN")
    )
  ),
  id = "main_navigation",
  selected = initial_view,
  theme = APP_THEME,
  fillable = FALSE,
  header = tags$head(
    tags$link(rel = "stylesheet", href = "app.css"),
    tags$script(src = "workbench.js"),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$script(HTML("(function () {
      var labRequests = new Map();
      function registerLabAck() {
        if (!window.Shiny || typeof Shiny.addCustomMessageHandler !== 'function') return;
        Shiny.addCustomMessageHandler('teaching-lab-ack', function (reply) {
          var request = labRequests.get(reply.id);
          if (request) { request.source.postMessage(Object.assign({ type: 'pk-lego-model-ack' }, reply), request.origin); labRequests.delete(reply.id); }
        });
      }
      document.addEventListener('shiny:connected', registerLabAck);
      document.addEventListener('DOMContentLoaded', registerLabAck);
      window.addEventListener('message', function (event) {
        var query = new URLSearchParams(window.location.search);
        var payload = event.data || {};
        if (query.get('bridge') !== 'lego' || event.source !== window.opener) return;
        if (query.has('origin') && event.origin !== query.get('origin')) return;
        if (payload.type !== 'pk-lego-model' || typeof payload.code !== 'string') return;
        if (!window.Shiny || typeof window.Shiny.setInputValue !== 'function' || !window.Shiny.shinyapp || !window.Shiny.shinyapp.config || !window.Shiny.shinyapp.config.sessionId) return;
        if (payload.teaching) {
          if (typeof payload.id !== 'string' || !/^[a-zA-Z0-9-]{1,80}$/.test(payload.id)) return;
          if (labRequests.has(payload.id)) return;
          labRequests.set(payload.id, { source: event.source, origin: event.origin });
          window.Shiny.setInputValue('teaching_lab_import', { id: payload.id, spec: payload.teaching }, { priority: 'event' });
          return;
        }
        window.Shiny.setInputValue('lego_model_import', {
          code: payload.code,
          spec: payload.spec || null,
          name: payload.name || 'lego_model',
          nonce: Date.now()
        }, { priority: 'event' });
        event.source.postMessage({ type: 'pk-lego-model-ack' }, event.origin);
      });
      function labelSidebarToggles() {
        document.querySelectorAll('.collapse-toggle').forEach(function (button) {
          var open = button.getAttribute('aria-expanded') === 'true';
          var english = new URLSearchParams(window.location.search).get('lang') === 'en';
          var label = open ? (english ? 'Close configuration' : 'Fermer la configuration') : (english ? 'Open configuration' : 'Ouvrir la configuration');
          button.setAttribute('title', label);
          button.setAttribute('aria-label', label);
        });
      }
      document.addEventListener('DOMContentLoaded', function () {
        labelSidebarToggles();
        new MutationObserver(labelSidebarToggles).observe(document.body, {
          subtree: true,
          attributes: true,
          attributeFilter: ['aria-expanded']
        });
      });
    })();"))
  ),
  nav_panel(
    "Analyse",
    value = "analysis",
    div(
      class = "analysis-shell",
      div(
        class = "safety-banner",
        div(
          tags$strong("Outil de recherche et d'enseignement"),
          span("Les résultats ne remplacent ni la validation locale du modèle ni le jugement clinique.")
        ),
        span(class = "engine-badge", if (PMETRICS_ENABLED) "mrgsolve + mapbayr / Pmetrics" else "mrgsolve + mapbayr")
      ),
      uiOutput("pd_tdm_return"),
      tags$button(
        type = "button",
        class = "mobile-configure-button",
        onclick = "var button=document.querySelector('.bslib-sidebar-layout .collapse-toggle');if(button&&button.getAttribute('aria-expanded')!=='true'){button.click();}",
        "Configurer et lancer"
      ),
      layout_sidebar(
        sidebar = sidebar(
          width = 390,
          open = "desktop",
          div(class = "sidebar-heading", h2("Configuration"), uiOutput("engine_status")),
          accordion(
            id = "settings_accordion",
            open = c("model_settings", "target_settings"),
            accordion_panel(
              "Modèle",
              value = "model_settings",
              radioButtons(
                "model_source",
                NULL,
                choices = stats::setNames(
                  c("library", if (PMETRICS_ENABLED) "pmetrics", "custom"),
                  c("Biblioth\u00e8que", if (PMETRICS_ENABLED) "Pmetrics", "Atelier Lego / C++")
                ),
                selected = "library",
                inline = TRUE
              ),
              conditionalPanel(
                "input.model_source == 'library'",
                selectInput("model_id", "Modèle principal", choices = model_choices, selected = DEFAULT_MODEL)
              ),
              conditionalPanel(
                "input.model_source == 'pmetrics'",
                selectInput("pmetrics_artifact_id", app_t(lang, "Artefact Pmetrics", "Pmetrics artifact"), choices = pmetrics_choices, selected = DEFAULT_PMETRICS_ARTIFACT),
                uiOutput("pmetrics_context_ui")
              ),
              conditionalPanel(
                "input.model_source == 'custom'",
                div(
                  class = "custom-note local",
                  if (ALLOW_CUSTOM_MODELS) {
                    "Mode local : le C++ libre est autorisé et s'exécute avec les droits du processus R."
                  } else {
                    "Serveur public : seuls les modèles portant la spécification contrôlée de l'Atelier Lego sont acceptés. Le C++ libre reste bloqué."
                  }
                )
              ),
              selectInput(
                "administration_route",
                "Voie d'administration",
                choices = "IV",
                selected = "IV"
              ),
              uiOutput("administration_mode_ui"),
              conditionalPanel(
                "input.model_source == 'library'",
                uiOutput("model_context_ui"),
                uiOutput("ml_status_ui"),
                checkboxInput("enable_averaging", "Activer le model averaging", FALSE),
                conditionalPanel(
                  "input.enable_averaging == true",
                  div(class = "custom-note", "Seuls les modèles de la même molécule, de la même voie et du même mode d'administration sont proposés."),
                  uiOutput("average_models_ui"),
                  selectInput(
                    "weighting_scheme",
                    "Pondération",
                    choices = stats::setNames(c("AIC", "LL"), if (lang == "en") c("Akaike criterion", "Log likelihood") else c("Crit\u00e8re d'Akaike", "Log-vraisemblance")),
                    selected = "AIC"
                  )
                )
              ),
              actionButton("validate_model", "Valider le modèle", class = "btn-outline-primary w-100")
            ),
            accordion_panel(
              "Cible et grille de doses",
              value = "target_settings",
              uiOutput("target_preset_ui"),
              selectInput(
                "target_metric",
                "Métrique",
                choices = c(
                  "AUC 0-24 h" = "AUC24",
                  "Concentration minimale" = "Cmin",
                  "Concentration maximale" = "Cmax",
                  "% du temps au-dessus d'une concentration" = "TimeAbove",
                  "Cmin et Cmax simultanément" = "CminCmax"
                ),
                selected = DEFAULT_TARGET_PRESET$target_metric %||% "AUC24"
              ),
              conditionalPanel(
                "input.target_metric == 'TimeAbove'",
                numericInput("target_concentration_threshold", "Seuil de concentration", 1, min = 0, step = 0.1),
                div(class = "custom-note", "Les bornes ci-dessous indiquent le pourcentage cible de l'intervalle posologique passé au-dessus du seuil.")
              ),
              conditionalPanel(
                "input.target_metric == 'CminCmax'",
                div(class = "custom-note", if (lang == "en") "Choose an independent rule for Cmin and Cmax. Both rules must be met." else "Choisissez une règle indépendante pour Cmin et Cmax. Les deux règles doivent être respectées."),
                selectInput(
                  "target_cmin_operator",
                  "Cmin",
                  choices = stats::setNames(c("range", "above", "below"), if (lang == "en") c("Within a range", "At or above a threshold", "At or below a threshold") else c("Dans un intervalle", "Au-dessus ou égal à un seuil", "En dessous ou égal à un seuil")),
                  selected = "range"
                ),
                conditionalPanel(
                  "input.target_cmin_operator == 'range'",
                  fluidRow(
                    column(6, numericInput("target_cmin_low", "Cmin basse", 5, min = 0, step = 0.1)),
                    column(6, numericInput("target_cmin_high", "Cmin haute", 15, min = 0, step = 0.1))
                  )
                ),
                conditionalPanel(
                  "input.target_cmin_operator != 'range'",
                  numericInput("target_cmin_threshold", if (lang == "en") "Cmin threshold" else "Seuil de Cmin", 10, min = 0, step = 0.1)
                ),
                selectInput(
                  "target_cmax_operator",
                  "Cmax",
                  choices = stats::setNames(c("range", "above", "below"), if (lang == "en") c("Within a range", "At or above a threshold", "At or below a threshold") else c("Dans un intervalle", "Au-dessus ou égal à un seuil", "En dessous ou égal à un seuil")),
                  selected = "range"
                ),
                conditionalPanel(
                  "input.target_cmax_operator == 'range'",
                  fluidRow(
                    column(6, numericInput("target_cmax_low", "Cmax basse", 20, min = 0, step = 0.1)),
                    column(6, numericInput("target_cmax_high", "Cmax haute", 40, min = 0, step = 0.1))
                  )
                ),
                conditionalPanel(
                  "input.target_cmax_operator != 'range'",
                  numericInput("target_cmax_threshold", if (lang == "en") "Cmax threshold" else "Seuil de Cmax", 30, min = 0, step = 0.1)
                )
              ),
              conditionalPanel(
                "input.target_metric != 'CminCmax'",
                fluidRow(
                  column(6, numericInput("target_low", "Borne basse", DEFAULT_TARGET_PRESET$target_low %||% 400, min = 0)),
                  column(6, numericInput("target_high", "Borne haute", DEFAULT_TARGET_PRESET$target_high %||% 600, min = 0))
                )
              ),
              fluidRow(
                column(4, numericInput("dose_min", "Dose min", DEFAULT_TARGET_PRESET$dose_min %||% 250, min = 0, step = 50)),
                column(4, numericInput("dose_max", "Dose max", DEFAULT_TARGET_PRESET$dose_max %||% 2000, min = 0, step = 50)),
                column(4, numericInput("dose_step", "Pas", DEFAULT_TARGET_PRESET$dose_step %||% 250, min = 1, step = 50))
              ),
              checkboxGroupInput(
                "candidate_intervals",
                "Intervalles testés",
                choices = stats::setNames(STANDARD_CANDIDATE_INTERVALS, paste(STANDARD_CANDIDATE_INTERVALS, "h")),
                selected = DEFAULT_TARGET_PRESET$intervals %||% c(8, 12, 24),
                inline = TRUE
              ),
              textInput(
                "candidate_intervals_custom",
                if (lang == "en") "Other intervals (h)" else "Autres intervalles (h)",
                value = "",
                placeholder = if (lang == "en") "e.g. 72, 168, 336" else "p. ex. 72, 168, 336"
              ),
              conditionalPanel(
                "input.administration_route == 'IV'",
                numericInput("future_infusion", "Durée de perfusion (h, 0 = bolus IV)", DEFAULT_TARGET_PRESET$infusion %||% 1, min = 0, step = 0.25)
              ),
              conditionalPanel(
                "input.administration_route == 'Oral'",
                div(class = "route-readonly compact", "Voie orale : durée de perfusion fixée à 0 h.")
              ),
              conditionalPanel(
                "input.administration_route == 'IM'",
                div(class = "route-readonly compact", "Voie intramusculaire : perfusion externe fixée à 0 h; absorption définie par le modèle.")
              )
            )
          ),
          div(
            class = "analysis-actions",
            checkboxInput(
              "accept_disclaimer",
              "J'ai lu le statut et j'accepte d'utiliser ce prototype sous ma responsabilité.",
              FALSE
            ),
            actionButton("run_analysis", "Lancer l'analyse", class = "btn-primary run-button w-100")
          )
        ),
        div(
          class = "workspace",
          navset_tab(
            id = "analysis_tabs",
            nav_panel(
              "Données",
              value = "data",
              div(
                class = "workspace-section",
                div(
                  class = "section-heading",
                  div(
                    h1("Historique thérapeutique"),
                    p("Saisissez des heures relatives, des jours avant une référence ou des dates civiles. Le moteur normalise ensuite l'origine sur la première administration.")
                  ),
                  div(
                    class = "patient-file-actions",
                    fileInput(
                      "upload_patient",
                      NULL,
                      accept = c("application/json", ".json"),
                      buttonLabel = "Importer JSON",
                      placeholder = "Aucun fichier"
                    ),
                    downloadButton("download_patient", "Exporter JSON", class = "btn-outline-secondary btn-sm")
                  )
                ),
                div(
                  class = "privacy-notice",
                  tags$strong("Session uniquement"),
                  span("Le fichier est lu, validé puis supprimé du stockage temporaire. Aucun identifiant, dossier patient ou code C++ n'est enregistré.")
                ),
                div(
                  class = "timeline-controls",
                  selectInput(
                    "time_entry_mode",
                    "Format temporel",
                    choices = stats::setNames(
                      c("relative_hours", "relative_days", "calendar"),
                      if (lang == "en") c("Relative hours", "Days before a reference date", "Dates and times") else c("Heures relatives", "Jours avant une date de r\u00e9f\u00e9rence", "Dates et heures")
                    ),
                    selected = "relative_hours"
                  ),
                  conditionalPanel(
                    "input.time_entry_mode != 'relative_hours'",
                    dateInput("reference_date", "Date de référence", value = Sys.Date(), format = if (lang == "en") "yyyy-mm-dd" else "dd/mm/yyyy", language = lang),
                    textInput("reference_clock", "Heure de référence", value = format(Sys.time(), "%H:%M"))
                  )
                ),
                h3("Administrations"),
                div(id = "dose_rows", class = "record-list", dose_row_ui(1L, lang = lang)),
                div(
                  class = "row-actions",
                  actionButton("add_dose", "Ajouter une administration", class = "btn-outline-secondary btn-sm"),
                  actionButton("remove_dose", "Retirer la dernière", class = "btn-outline-secondary btn-sm")
                ),
                hr(),
                h3("Concentrations observées"),
                uiOutput("observation_covariate_notice"),
                div(id = "observation_rows", class = "record-list", observation_row_ui(1L, lang = lang)),
                div(
                  class = "row-actions",
                  actionButton("add_observation", "Ajouter un prélèvement", class = "btn-outline-secondary btn-sm"),
                  actionButton("remove_observation", "Retirer le dernier", class = "btn-outline-secondary btn-sm")
                ),
                uiOutput("data_quality_preview"),
                plotOutput("therapy_timeline", height = "250px")
              )
            ),
            nav_panel(
              "Ajustement",
              value = "fit",
              div(
                class = "workspace-section",
                uiOutput("fit_empty"),
                div(
                  class = "table-toolbar result-toolbar",
                  h3("Ajustement bayésien"),
                  downloadButton("download_report", "Créer le rapport", class = "btn-outline-primary btn-sm")
                ),
                uiOutput("analysis_warnings"),
                uiOutput("current_exposure_summary"),
                plotOutput("fit_plot", height = "430px"),
                h3("Doses supplémentaires : poursuivre ou modifier"),
                uiOutput("future_comparison_summary"),
                plotOutput("future_comparison_plot", height = "430px"),
                uiOutput("ml_interpretation_panel"),
                h3("Sensibilité du model averaging"),
                uiOutput("averaging_sensitivity_summary"),
                DTOutput("averaging_sensitivity_table"),
                h3("Modèles et pondérations"),
                DTOutput("model_table"),
                conditionalPanel(
                  "input.model_source == 'pmetrics'",
                  h3(app_t(lang, "Distribution Pmetrics : points de support", "Pmetrics distribution: support points")),
                  p(app_t(lang, "La courbe centrale utilise les paramètres moyens. Les probabilités ci-dessous conservent la distribution jointe discrète, sans approximation normale.", "The central curve uses mean parameters. The probabilities below preserve the discrete joint distribution, without a normal approximation.")),
                  DTOutput("pmetrics_support_table")
                )
              )
            ),
            nav_panel(
              "Posologies",
              value = "dosing",
              div(
                class = "workspace-section",
                uiOutput("recommendation_summary"),
                plotOutput("recommended_profile", height = "410px"),
                h3("Distribution prédictive du meilleur scénario"),
                uiOutput("distribution_summary"),
                plotOutput("distribution_plot", height = "330px"),
                div(class = "table-toolbar", h3("Scénarios classés"), downloadButton("download_scenarios", "Exporter CSV", class = "btn-outline-secondary btn-sm")),
                DTOutput("recommendation_table")
              )
            ),
            nav_panel(
              "Modèle",
              value = "model",
              div(
                class = "workspace-section model-workspace",
                uiOutput("model_contract"),
                conditionalPanel(
                  "input.model_source == 'library'",
                  pre(class = "model-code", textOutput("library_code", container = span))
                ),
                conditionalPanel(
                  "input.model_source == 'custom'",
                  textAreaInput("custom_code", "Modèle mrgsolve complet", value = DEFAULT_CODE, rows = 34, width = "100%")
                ),
                conditionalPanel(
                  "input.model_source == 'pmetrics'",
                  pre(class = "model-code", textOutput("pmetrics_artifact_code", container = span))
                )
              )
            )
          )
        )
      )
    )
  ),
  ddi_panel(lang),
  pd_panel(lang),
  nav_panel(
    "Bibliothèque",
    value = "library",
    div(
      class = "plain-page",
      h1("Bibliothèque de modèles"),
      p("Chaque modèle est associé à l'article qui décrit ses paramètres et sa population. Toute nouvelle entrée est relue avant publication."),
      DTOutput("catalog_table")
    )
  ),
  nav_panel(
    "Réglages",
    value = "settings",
    div(
      class = "plain-page settings-page",
      h1("Réglages de simulation"),
      p("Ces paramètres s'appliquent à la prochaine analyse et sont inclus dans le rapport."),
      div(
        class = "settings-grid",
        div(
          class = "settings-group",
          h2("Résolution et projection"),
          selectInput(
            "simulation_delta",
            "Pas de simulation",
            choices = stats::setNames(c(0.05, 0.1, 0.25), if (lang == "en") c("High · 0.05 h", "Standard · 0.1 h", "Fast · 0.25 h") else c("Haute \u00b7 0,05 h", "Standard \u00b7 0,1 h", "Rapide \u00b7 0,25 h")),
            selected = 0.1
          ),
          numericInput("additional_doses", "Doses supplémentaires comparées", 6, min = 1, max = 30, step = 1),
          checkboxInput("show_component_profiles", "Afficher chaque modèle sur l'ajustement", TRUE)
        ),
        div(
          class = "settings-group",
          h2("Distribution prédictive"),
          numericInput("mc_replicates", "Réplications Monte Carlo", 250, min = 50, max = 1000, step = 50),
          selectInput("prediction_interval", "Intervalle prédictif", choices = c("80 %" = 80, "90 %" = 90, "95 %" = 95), selected = 90),
          checkboxGroupInput(
            "variability_components",
            "Composantes simulées",
            choices = stats::setNames(
              c("POSTERIOR", "RESID", "TIMING"),
              if (lang == "en") c("MAP posterior uncertainty", "Residual error", "Timing uncertainty") else c("Incertitude post\u00e9rieure MAP", "Erreur r\u00e9siduelle", "Incertitude des horaires")
            ),
            selected = "POSTERIOR"
          ),
          numericInput("posterior_replicates", "Réestimations avec horaires incertains", 20, min = 10, max = 40, step = 5),
          selectInput(
            "blq_method",
            "Concentrations BLQ",
            choices = if (lang == "en") c("Exclude from fit" = "exclude", "LLOQ / 2 (exploratory)" = "lloq_half") else c("Exclure de l'ajustement" = "exclude", "LLOQ / 2 (exploratoire)" = "lloq_half"),
            selected = "exclude"
          )
        ),
        div(
          class = "settings-group",
          h2("Erreur résiduelle MAP"),
          selectInput(
            "residual_error_mode",
            "Erreur utilisée pour l'ajustement",
            choices = c("Conserver l'erreur du modèle" = "model", "Fixer un CV proportionnel" = "fixed_cv"),
            selected = "model"
          ),
          conditionalPanel(
            "input.residual_error_mode == 'fixed_cv'",
            numericInput("fixed_residual_cv", "CV proportionnel fixé (%)", 1, min = 0.1, max = 100, step = 0.1)
          ),
          div(
            class = "custom-note",
            "Le mode fixé remplace les SIGMA de l'ajustement par le CV proportionnel choisi, ou son équivalent logarithmique pour un modèle exponentiel. Il modifie le poids des concentrations dans MAP-BE."
          )
        )
      ),
      div(class = "privacy-notice settings-privacy", tags$strong("Pas de lien partageable"), span("Les données patient ne sont pas encodées dans l'URL afin d'éviter leur présence dans l'historique du navigateur, les journaux réseau ou les outils d'analytique."))
    )
  ),
  nav_panel(
    "Méthode",
    value = "method",
    div(
      class = "plain-page method-page",
      h1("Méthode et limites"),
      h2("Estimation individuelle"),
      p("Les effets aléatoires individuels sont estimés par maximum a posteriori avec mapbayr. Les administrations, covariables datées et concentrations sont converties en événements NM-TRAN puis simulées avec mrgsolve."),
      p(app_t(lang,
        "Pour un artefact Pmetrics, PM_model$map() calcule les probabilités postérieures des points de support pondérés. La prédiction centrale utilise la moyenne postérieure des paramètres; la distribution prédictive échantillonne directement les points de support selon leurs probabilités postérieures. Les fichiers de calcul Pmetrics sont créés dans un dossier temporaire propre à l'appel puis supprimés.",
        "For a Pmetrics artifact, PM_model$map() computes posterior probabilities over the weighted support points. The central prediction uses the posterior parameter mean; the predictive distribution samples support points directly from their posterior probabilities. Pmetrics run files are created in a call-scoped temporary directory and then removed."
      )),
      h2("Erreur résiduelle"),
      p("Par défaut, l'ajustement conserve les variances résiduelles du modèle publié. Le réglage CV proportionnel fixé remplace chaque paire SIGMA par une variance proportionnelle commune, ou sa variance logarithmique équivalente si le modèle utilise une erreur exponentielle. Ce choix change le poids des observations dans l'estimation MAP et doit être documenté."),
      h2("Model averaging"),
      p("Chaque modèle analyse les mêmes données. Les prédictions sont moyennées avec des poids issus de la vraisemblance ou du critère d'Akaike. Le serveur refuse l'agrégation de modèles ne partageant pas la même molécule, la même voie et le même mode d'administration."),
      p("La robustesse est explorée avec des poids égaux, AIC, log-vraisemblance et des analyses laissant successivement de côté chaque modèle. Une divergence des doses proposées doit conduire à revoir l'applicabilité du model averaging."),
      h2("Exposition et scénarios"),
      p("L'exposition historique sur les dernières 24 heures est distinguée de l'exposition à l'état stationnaire du dernier schéma. La projection compare le maintien de la dernière posologie à l'application du scénario classé en tête, sur le nombre de doses supplémentaires choisi dans Réglages."),
      h2("Cibles pharmacocinétiques"),
      p("La cible temps au-dessus d'un seuil est calculée sur un intervalle posologique complet à l'état stationnaire et exprimée en pourcentage. Les franchissements du seuil sont interpolés entre deux points simulés. Pour la cible combinée, Cmin et Cmax suivent chacun une contrainte indépendante: intervalle, valeur minimale ou valeur maximale. Les deux contraintes doivent être respectées; aucune moyenne entre les deux critères n'est utilisée."),
      p("Une administration déclarée à l'état stationnaire est initialisée à t = 0 avec ss = 1. Les doses de la projection sont ensuite ajoutées explicitement avec ss = 0; la comparaison des doses supplémentaires continue donc au-delà de cet état initial."),
      p("Pour un modèle Lego qui répartit une dose entre une ou plusieurs voies d'entrée, l'état stationnaire est approché par 50 administrations explicites avant t = 0. Cette phase de préchauffage évite d'appliquer séparément ss = 1 à chaque voie, ce qui serait incorrect pour des cinétiques parallèles ou non linéaires."),
      p("Les scénarios sont d'abord classés sur leur prédiction moyenne. La probabilité d'atteindre la cible, de sous-exposition et de surexposition est ensuite simulée pour les douze scénarios moyens les plus proches, puis utilisée pour leur classement final."),
      h2("Qualité des données"),
      p("Les doses oubliées sont exclues. Les horaires incertains peuvent faire l'objet de réestimations de sensibilité. Les valeurs BLQ sont exclues par défaut; l'imputation LLOQ/2 est uniquement exploratoire. La matrice, les unités et le domaine exact des covariables ne sont jamais convertis ou déduits automatiquement."),
      h2("Distribution"),
      p("La distribution prédictive est exploratoire. Après un ajustement, elle repose sur un bootstrap paramétrique de l'estimation MAP, auquel peuvent s'ajouter l'erreur résiduelle, l'incertitude des horaires et l'incertitude entre modèles. Sans concentration exploitable, elle revient à une simulation populationnelle."),
      h2("Apprentissage automatique"),
      p("Le module expérimental adapte à chaque modèle de la bibliothèque la méthodologie de simulation publiée pour le tacrolimus (doi:10.1016/j.phrs.2021.105578). XGBoost apprend le logarithme du rapport entre l'AUC24 individuelle simulée et l'AUC24 populationnelle du même schéma. Il utilise les concentrations et horaires, leurs prédictions populationnelles, la dose, l'intervalle, la durée de perfusion et les covariables du modèle. Une décomposition locale DALEX explique la prédiction par rapport à un échantillon de référence entièrement synthétique."),
      p("Le prédicteur opérationnel exige une administration déclarée à l'état stationnaire (ss = 1), un C0 pré-dose et une concentration post-dose dans l'intervalle suivant. Le prélèvement post-dose le plus proche de la cible est utilisé : 1 h après la dose pour l'oral et le bolus, ou 1 h après la fin d'une perfusion intermittente. Le benchmark C0 seul est présenté sur le site comme comparateur, mais n'est pas utilisé par l'artefact opérationnel."),
      p("L'estimation ML est une AUC24 expérimentale. Elle est comparée à l'AUC24 MAP du schéma répété à l'état stationnaire, et non à l'intégrale MAP historique partielle lorsque moins de 24 h sont disponibles. En model averaging, chaque modèle produit sa propre AUC24 ML, puis les estimations sont agrégées avec les poids de l'analyse. Tous les modèles doivent disposer d'un artefact compatible et partager la même voie et le même mode d'administration."),
      p("Chaque artefact reste lié au même modèle PK, à la même voie, au même mode d'administration, au même schéma de variables et aux empreintes exactes du fichier mrgsolve et du booster. Le statut recherche exige le respect des seuils préspecifiés en validation croisée répétée et sur un test interne non touché; un artefact évalué mais sous ces seuils reste signalé comme expérimental. La transportabilité vers un autre modèle PopPK est rapportée séparément."),
      p("Une variable hors des bornes empiriques d'entraînement déclenche un avertissement sans masquer l'estimation. Un C0 ou un prélèvement post-dose manquant, l'absence d'état stationnaire, une perfusion continue, un protocole incompatible ou une AUC invalide restent bloquants. Ces entraînements synthétiques ne constituent pas une validation clinique pour les autres molécules. L'estimation ML ne remplace pas les projections de dose MAP tant qu'une validation favorable sur des patients réels indépendants n'est pas documentée."),
      h2("Sécurité"),
      p("Le serveur public ne compile jamais directement le C++ reçu. Pour un modèle Atelier Lego, il extrait une spécification JSON, la valide, régénère lui-même le code mrgsolve puis compile uniquement ce code contrôlé. Tout autre C++ reste refusé tant qu'il n'est pas exécuté dans un conteneur éphémère isolé."),
      p("Les imports JSON sont traités dans la session Shiny et leur fichier temporaire est supprimé immédiatement après lecture. Les exports sont produits à la demande sans base de données.")
    )
  ),
  nav_panel(
    "Statut",
    value = "status",
    div(
      class = "plain-page status-page",
      h1("Statut"),
      div(
        class = "status-disclaimer",
        tags$strong("Avertissement obligatoire"),
        p("Ce prototype est destiné à la recherche et à l'enseignement. Il n'est pas enregistré comme dispositif médical, ne garantit ni l'exactitude d'un résultat ni son applicabilité à un patient particulier et ne remplace pas le jugement clinique."),
        p("Toute décision de dose reste sous la responsabilité du professionnel de santé et exige la vérification de la voie, des horaires, des unités, de la population source, des covariables, des concentrations, de la fonction d'organe et des recommandations locales. Une validation indépendante et une gouvernance documentée sont nécessaires avant toute utilisation clinique."),
        p(
          "Si vous recherchez un véritable dispositif médical avec expertise, consultez ",
          tags$a(
            "ABIS du CHU de Limoges.",
            href = "https://abis.chu-limoges.fr/",
            target = "_blank",
            rel = "noopener noreferrer"
          )
        )
      )
    )
  ),
  footer = div(class = "app-footer", app_t(lang,
    paste0("PMx Explain · moteurs R mrgsolve/mapbayr", if (PMETRICS_ENABLED) " et Pmetrics" else "", " · aucun dossier patient n'est persisté"),
    paste0("PMx Explain · R engines mrgsolve/mapbayr", if (PMETRICS_ENABLED) " and Pmetrics" else "", " · no patient record is persisted")
  ))
), lang)
}

server <- function(input, output, session) {
  current_language <- reactive(app_language_from_query(session$clientData$url_search %||% ""))
  tx <- function(fr, en = NULL) app_t(current_language(), fr, en)
  localized <- function(value) localize_ui(value, current_language())
  current_model_id <- reactive({
    value <- input$model_id %||% ""
    if (length(value) && nzchar(value[[1]]) && value[[1]] %in% ANALYSIS_MODEL_IDS) value[[1]] else DEFAULT_MODEL
  })
  current_pmetrics_artifact <- reactive({
    id <- input$pmetrics_artifact_id %||% DEFAULT_PMETRICS_ARTIFACT
    shiny::req(nzchar(id))
    read_pmetrics_artifact(id)
  })
  current_administration_mode <- reactive({
    route <- input$administration_route %||% ""
    shiny::req(nzchar(route))
    if (identical(input$model_source %||% "library", "custom")) {
      return(if (identical(route, "Oral")) "ORAL" else if (identical(route, "IM")) "IM" else "IV_INTERMITTENT")
    }
    if (identical(input$model_source %||% "library", "pmetrics")) {
      artifact <- current_pmetrics_artifact()
      shiny::req(identical(route, artifact$route))
      return(artifact$mode)
    }
    record <- model_record(current_model_id())
    shiny::req(model_supports_route(record, route))
    modes <- model_administration_modes(record, route)
    shiny::req(length(modes))
    selected <- input$administration_mode %||% ""
    if (selected %in% modes) selected else modes[[1]]
  })
  dose_count <- reactiveVal(1L)
  observation_count <- reactiveVal(1L)
  analysis_store <- reactiveVal(NULL)
  ddi_store <- reactiveVal(NULL)
  ddi_assembly <- reactiveVal(NULL)
  ddi_custom_1 <- reactiveVal(NULL)
  ddi_custom_2 <- reactiveVal(NULL)
  ddi_subject_1 <- reactiveVal(NULL)
  ddi_subject_2 <- reactiveVal(NULL)
  ddi_target_count <- reactiveVal(1L)
  validation_store <- reactiveVal(NULL)
  query_applied <- reactiveVal(FALSE)
  pending_import_target <- reactiveVal(NULL)
  workbench_import <- reactiveVal(NULL)
  pending_ddi_target <- reactiveVal(NULL)
  pending_ddi_routes <- reactiveVal(list())
  last_target_preset_key <- reactiveVal(NULL)
  pending_import_covariates <- reactiveVal(NULL)
  pending_import_dose_covariates <- reactiveVal(NULL)
  pending_lego_covariates <- reactiveVal(NULL)
  session_model_dir <- new_private_workspace("pk-mipd-custom-session-")
  session_model_cache <- new.env(parent = emptyenv())

  numeric_input_value <- function(id, fallback = NA_real_) {
    value <- suppressWarnings(as.numeric(input[[id]] %||% fallback))
    if (length(value) != 1L || !is.finite(value)) fallback else value
  }

  candidate_interval_values <- function() {
    parse_dosing_intervals(
      isolate(input$candidate_intervals %||% character()),
      isolate(input$candidate_intervals_custom %||% "")
    )
  }

  update_candidate_intervals <- function(values) {
    values <- suppressWarnings(as.numeric(unlist(values %||% numeric())))
    values <- sort(unique(values[is.finite(values) & values > 0]))
    standard <- intersect(values, STANDARD_CANDIDATE_INTERVALS)
    custom <- setdiff(values, STANDARD_CANDIDATE_INTERVALS)
    updateCheckboxGroupInput(session, "candidate_intervals", selected = as.character(standard))
    updateTextInput(session, "candidate_intervals_custom", value = paste(custom, collapse = ", "))
  }

  update_target_inputs <- function(target) {
    if ((target$metric %||% "") %in% c("AUC24", "Cmin", "Cmax", "TimeAbove", "CminCmax")) {
      updateSelectInput(session, "target_metric", selected = target$metric)
    }
    cmin_operator <- target$cminOperator %||% "range"
    cmax_operator <- target$cmaxOperator %||% "range"
    if (cmin_operator %in% c("range", "above", "below")) updateSelectInput(session, "target_cmin_operator", selected = cmin_operator)
    if (cmax_operator %in% c("range", "above", "below")) updateSelectInput(session, "target_cmax_operator", selected = cmax_operator)
    target$cminLow <- target$cminLow %||% target$low
    target$cminHigh <- target$cminHigh %||% target$high
    target$cminThreshold <- target$cminThreshold %||% target$low
    target$cmaxThreshold <- target$cmaxThreshold %||% if (identical(cmax_operator, "below")) target$cmaxHigh else target$cmaxLow
    numeric_targets <- c(
      low = "target_low",
      high = "target_high",
      concentrationThreshold = "target_concentration_threshold",
      cminLow = "target_cmin_low",
      cminHigh = "target_cmin_high",
      cminThreshold = "target_cmin_threshold",
      cmaxLow = "target_cmax_low",
      cmaxHigh = "target_cmax_high",
      cmaxThreshold = "target_cmax_threshold",
      doseMin = "dose_min",
      doseMax = "dose_max",
      doseStep = "dose_step",
      infusion = "future_infusion"
    )
    for (name in names(numeric_targets)) {
      value <- suppressWarnings(as.numeric(target[[name]] %||% NA_real_))
      if (is.finite(value) && value >= 0) updateNumericInput(session, numeric_targets[[name]], value = value)
    }
    if (length(target$intervals %||% numeric())) update_candidate_intervals(target$intervals)
  }

  row_time_messages <- function(prefix, index) {
    mode <- input$time_entry_mode %||% "relative_hours"
    if (identical(mode, "relative_hours")) {
      value <- numeric_input_value(paste0(prefix, "_time_", index))
      return(if (!is.finite(value)) "Temps relatif manquant ou invalide." else character())
    }
    if (identical(mode, "relative_days")) {
      days <- numeric_input_value(paste0(prefix, "_days_ago_", index))
      clock <- trimws(as.character(input[[paste0(prefix, "_clock_", index)]] %||% ""))
      return(c(
        if (!is.finite(days) || days < 0) "Le nombre de jours doit être positif.",
        if (!grepl("^([01]?[0-9]|2[0-3]):[0-5][0-9]$", clock)) "L'heure doit respecter HH:MM."
      ))
    }
    value <- trimws(as.character(input[[paste0(prefix, "_calendar_", index)]] %||% ""))
    if (!grepl("^(\\d{1,2}/\\d{1,2}(/\\d{2,4})?|\\d{4}-\\d{2}-\\d{2})\\s+([01]?\\d|2[0-3]):[0-5]\\d$", value)) {
      "Utilisez JJ/MM HH:MM, JJ/MM/AAAA HH:MM ou AAAA-MM-JJ HH:MM."
    } else {
      character()
    }
  }

  dose_row_messages <- function(index) {
    status <- input[[paste0("dose_status_", index)]] %||% "administered"
    amount <- numeric_input_value(paste0("dose_amount_", index))
    interval <- numeric_input_value(paste0("dose_interval_", index))
    uncertainty <- if (identical(status, "uncertain")) numeric_input_value(paste0("dose_time_uncertainty_", index), 0) else 0
    steady_state <- isTRUE(input[[paste0("dose_ss_", index)]])
    count <- if (steady_state) 1 else numeric_input_value(paste0("dose_count_", index))
    infusion <- if (identical(input$administration_route, "IV")) numeric_input_value(paste0("dose_infusion_", index)) else 0
    administration_definition <- administration_covariate_definition()
    invalid_covariates <- if (nrow(administration_definition)) {
      administration_definition$name[!vapply(administration_definition$name, function(name) {
        is.finite(numeric_input_value(paste0("dose_cov_", name, "_", index)))
      }, logical(1))]
    } else character()
    c(
      if (!steady_state) row_time_messages("dose", index),
      if (!is.finite(amount) || amount <= 0) "La dose doit être strictement positive.",
      if (!is.finite(interval) || interval < 0) "L'intervalle doit être positif ou nul.",
      if (steady_state && (!is.finite(interval) || interval <= 0)) "Le steady state exige un intervalle strictement positif.",
      if (!steady_state && (!is.finite(count) || count < 1)) "Le nombre d'administrations doit être au moins 1.",
      if (!is.finite(infusion) || infusion < 0) "La perfusion doit être positive ou nulle.",
      if (is.finite(infusion) && is.finite(interval) && interval > 0 && infusion > interval) "La perfusion dépasse l'intervalle entre deux doses.",
      if (!is.finite(uncertainty) || uncertainty < 0) "L'incertitude horaire doit être positive ou nulle.",
      if (identical(status, "uncertain") && (!is.finite(uncertainty) || uncertainty <= 0)) "Précisez une incertitude horaire pour cette dose.",
      if (identical(status, "missed")) "Cette dose sera exclue des administrations reçues.",
      if (length(invalid_covariates)) paste0("Covariable(s) d'administration invalide(s) : ", paste(invalid_covariates, collapse = ", "), ".")
    )
  }

  observation_row_messages <- function(index) {
    blq <- isTRUE(input[[paste0("observation_blq_", index)]])
    concentration <- numeric_input_value(paste0("observation_concentration_", index))
    lloq <- if (blq) numeric_input_value(paste0("observation_lloq_", index)) else 0
    uncertain <- isTRUE(input[[paste0("observation_time_uncertain_", index)]])
    uncertainty <- if (uncertain) numeric_input_value(paste0("observation_time_uncertainty_", index), 0) else 0
    matrix <- input[[paste0("observation_matrix_", index)]] %||% "unspecified"
    definition <- covariate_definition()
    invalid_covariates <- if (nrow(definition)) {
      definition$name[!vapply(definition$name, function(name) {
        is.finite(numeric_input_value(paste0("observation_cov_", name, "_", index)))
      }, logical(1))]
    } else {
      character()
    }
    c(
      row_time_messages("observation", index),
      if (!blq && (!is.finite(concentration) || concentration < 0)) "La concentration doit être positive ou nulle.",
      if (blq && (!is.finite(lloq) || lloq <= 0)) "Une LLOQ strictement positive est requise pour une valeur BLQ.",
      if (!is.finite(uncertainty) || uncertainty < 0) "L'incertitude horaire doit être positive ou nulle.",
      if (uncertain && (!is.finite(uncertainty) || uncertainty <= 0)) "Précisez une incertitude horaire pour ce prélèvement.",
      if (identical(matrix, "unspecified")) "Matrice biologique non précisée : vérifiez sa compatibilité avec l'article.",
      if (length(invalid_covariates)) paste0("Covariable(s) invalide(s) : ", paste(invalid_covariates, collapse = ", "), ".")
    )
  }

  validation_message_ui <- function(messages) {
    messages <- unique(messages[nzchar(messages)])
    if (!length(messages)) return(span(class = "validation-ok", "Ligne valide"))
    tags$ul(class = "validation-messages", lapply(messages, tags$li))
  }

  register_dose_validation <- function(index) {
    local({
      row_index <- index
      output_id <- paste0("dose_validation_", row_index)
      output[[output_id]] <- renderUI(localized(validation_message_ui(dose_row_messages(row_index))))
      shiny::outputOptions(output, output_id, suspendWhenHidden = FALSE)
    })
  }

  register_observation_validation <- function(index) {
    local({
      row_index <- index
      output_id <- paste0("observation_validation_", row_index)
      output[[output_id]] <- renderUI(localized(validation_message_ui(observation_row_messages(row_index))))
      shiny::outputOptions(output, output_id, suspendWhenHidden = FALSE)
    })
  }

  session$onSessionEnded(function() {
    analysis_store(NULL)
    ddi_store(NULL)
    ddi_assembly(NULL)
    ddi_custom_1(NULL)
    ddi_custom_2(NULL)
    ddi_subject_1(NULL)
    ddi_subject_2(NULL)
    validation_store(NULL)
    cleanup_private_workspace(session_model_dir, session_model_cache)
  })

  insert_dose_row <- function(index, time = 0, amount = 1000, interval = 12, count = 4, infusion = 1, steady_state = FALSE, status = "administered", time_uncertainty = 0) {
    register_dose_covariates(index)
    register_dose_validation(index)
    insertUI(
      selector = "#dose_rows",
      where = "beforeEnd",
      ui = dose_row_ui(index, time, amount, interval, count, infusion, steady_state, status, time_uncertainty, lang = current_language())
    )
  }

  insert_observation_row <- function(index, time = 47.5, concentration = 18, blq = FALSE, lloq = 0, matrix = "unspecified", time_uncertainty = 0) {
    register_observation_covariates(index)
    register_observation_validation(index)
    insertUI(
      selector = "#observation_rows",
      where = "beforeEnd",
      ui = observation_row_ui(index, time, concentration, blq, lloq, matrix, time_uncertainty, lang = current_language())
    )
  }

  observeEvent(input$add_dose, {
    index <- dose_count() + 1L
    dose_count(index)
    insert_dose_row(index, time = 48, amount = 1000, interval = 12, count = 4, infusion = 1, steady_state = FALSE)
  })

  observeEvent(input$remove_dose, {
    index <- dose_count()
    if (index <= 1L) return(showNotification(tx("Au moins une administration est requise.", "At least one administration is required."), type = "warning"))
    removeUI(selector = paste0("#dose_row_", index))
    dose_count(index - 1L)
  })

  observeEvent(input$add_observation, {
    index <- observation_count() + 1L
    observation_count(index)
    insert_observation_row(index, time = 60, concentration = 15)
  })

  observeEvent(input$remove_observation, {
    index <- observation_count()
    if (index <= 1L) return(showNotification(tx("Conservez une ligne et saisissez 0 pour une analyse populationnelle.", "Keep one row and enter 0 for a population analysis."), type = "message"))
    removeUI(selector = paste0("#observation_row_", index))
    observation_count(index - 1L)
  })

  observe({
    if (query_applied()) return()
    query <- parseQueryString(session$clientData$url_search %||% "")
    requested <- query$model %||% ""
    if (requested %in% ANALYSIS_MODEL_IDS) updateSelectInput(session, "model_id", selected = requested)
    if (identical(query$source %||% "", "custom")) updateRadioButtons(session, "model_source", selected = "custom")
    if (identical(query$view %||% "", "ddi")) {
      session$onFlushed(function() bslib::nav_select("main_navigation", "ddi"), once = TRUE)
    }
    query_applied(TRUE)
  })

  observe({
    if (identical(input$model_source %||% "library", "custom")) {
      routes <- c("IV", "Oral", "IM")
    } else if (identical(input$model_source %||% "library", "pmetrics")) {
      routes <- current_pmetrics_artifact()$route
    } else {
      model_id <- current_model_id()
      routes <- model_routes(model_record(model_id))
    }
    current <- isolate(input$administration_route %||% "")
    selected <- if (current %in% routes) current else routes[[1]]
    updateSelectInput(
      session,
      "administration_route",
      choices = stats::setNames(routes, routes),
      selected = selected
    )
  })

  shiny::observeEvent(input$teaching_lab_import, {
    payload <- input$teaching_lab_import
    id <- payload$id
    if (!is.character(id) || length(id) != 1 || !grepl("^[a-zA-Z0-9-]{1,80}$", id)) return()
    tryCatch({
      imported <- teaching_lab_import(payload$spec)
      pending_lego_covariates(list(code = imported$code, definition = parse_covariates(imported$code)))
      updateRadioButtons(session, "model_source", selected = "custom")
      updateTextAreaInput(session, "custom_code", value = imported$code)
      updateSelectInput(session, "administration_route", selected = imported$route)
      updateSelectInput(session, "time_entry_mode", selected = "relative_hours")
      updateCheckboxInput(session, "enable_averaging", value = FALSE)
      updateCheckboxInput(session, "enable_experimental_ml", value = FALSE)
      replace_dose_rows(imported$doses)
      replace_observation_rows(data.frame(time = 0, concentration = NA_real_))
      validation_store(NULL)
      analysis_store(NULL)
      bslib::nav_select("analysis_tabs", "model")
      session$sendCustomMessage("teaching-lab-ack", list(id = id, ok = TRUE))
      mode <- switch(payload$spec$lab,
        absorption = tx("voie orale", "oral route"),
        infusion = tx("perfusion IV", "IV infusion"),
        tx("bolus IV", "IV bolus")
      )
      showNotification(paste0(
        tx("Experience synthetique importee : ", "Synthetic experiment imported: "), mode,
        tx(", mg, L, h. Aucune concentration observee. Variabilite illustrative uniquement.",
           ", mg, L, h. No observed concentration. Illustrative variability only.")
      ), duration = 12)
    }, error = function(error) session$sendCustomMessage("teaching-lab-ack", list(id = id, ok = FALSE, error = conditionMessage(error))))
  }, ignoreInit = TRUE)

  shiny::observeEvent(input$lego_model_import, {
    payload <- input$lego_model_import
    code <- payload$code %||% ""
    if (!is.character(code) || length(code) != 1L || !nzchar(trimws(code))) return()
    if (nchar(code, type = "bytes") > 200000) {
      return(showNotification(tx("Le modèle Lego dépasse la taille autorisée.", "The Lego model exceeds the size limit."), type = "error"))
    }
    imported_spec <- tryCatch(
      if (!is.null(payload$spec)) normalize_lego_spec(payload$spec) else lego_spec_from_code(code),
      error = function(error) NULL
    )
    safe_code <- tryCatch(
      safe_lego_model_code(code = code, specification = payload$spec %||% NULL),
      error = function(error) {
        showNotification(conditionMessage(error), type = "error", duration = 10)
        NULL
      }
    )
    if (is.null(safe_code)) return()
    pending_lego_covariates(list(code = safe_code, definition = parse_covariates(safe_code)))
    updateRadioButtons(session, "model_source", selected = "custom")
    updateTextAreaInput(session, "custom_code", value = safe_code)
    updateSelectInput(session, "residual_error_mode", selected = "model")
    validation_store(NULL)
    analysis_store(NULL)
    bslib::nav_select("analysis_tabs", "model")
    eta_count <- length(imported_spec$population$iivVariances %||% numeric())
    sigma_count <- if (is.null(imported_spec)) 0L else 2L
    showNotification(tx(
      paste0("Modèle Lego régénéré côté serveur · Ω : ", eta_count, " ETA · Σ : ", sigma_count, " composante(s). Le mode d'erreur résiduelle du modèle a été rétabli."),
      paste0("Lego model regenerated on the server · Ω: ", eta_count, " ETA · Σ: ", sigma_count, " component(s). Model residual-error mode was restored.")
    ), type = "message", duration = 8)
  }, ignoreNULL = TRUE)

  selected_model_ids <- reactive({
    if (!identical(input$model_source %||% "library", "library")) return(character())
    primary_id <- current_model_id()
    ids <- if (!isTRUE(input$enable_averaging)) primary_id else unique(c(primary_id, input$average_model_ids %||% character()))
    route <- input$administration_route %||% ""
    shiny::req(nzchar(route))
    primary_record <- model_record(primary_id)
    shiny::req(model_supports_route(primary_record, route))
    mode <- current_administration_mode()
    compatible <- vapply(ids, function(id) model_supports_administration_mode(model_record(id), route, mode), logical(1))
    shiny::validate(shiny::need(all(compatible), tx(
      "Tous les modèles moyennés doivent utiliser la même voie et le même mode d'administration.",
      "All averaged models must use the same administration route and mode."
    )))
    ids
  })

  output$average_models_ui <- renderUI({
    primary_id <- current_model_id()
    drug <- model_record(primary_id)$drug[[1]]
    route <- input$administration_route %||% ""
    shiny::req(nzchar(route))
    primary_record <- model_record(primary_id)
    shiny::req(model_supports_route(primary_record, route))
    mode <- current_administration_mode()
    choices <- catalog_choices_i18n(drug, route, mode, current_language(), analysis_only = TRUE)
    drug_label <- if (identical(current_language(), "en")) model_record(primary_id)$drugEn[[1]] else drug
    mode_label <- administration_mode_label(mode, current_language())
    localized(checkboxGroupInput(
      "average_model_ids",
      tx(paste0("Modèles de ", drug_label, " · ", mode_label), paste0(drug_label, " models · ", mode_label)),
      choices = choices,
      selected = primary_id
    ))
  })
  shiny::outputOptions(output, "average_models_ui", suspendWhenHidden = FALSE)

  output$model_context_ui <- renderUI({
    record <- model_record(current_model_id())
    english <- identical(current_language(), "en")
    tags_value <- if (english) record$populationTagsEn[[1]] %||% character() else record$populationTags[[1]] %||% character()
    categories <- if (english) record$administrationCategoriesEn[[1]] else record$administrationCategories[[1]]
    doi <- as.character(record$doi[[1]] %||% "")
    drug <- if (english) record$drugEn[[1]] else record$drug[[1]]
    population <- if (english) record$populationEn[[1]] else record$population[[1]]
    div(
      class = "model-context",
      div(class = "model-context-head", tags$strong(record$model[[1]]), span(paste0(drug, " · ", paste(categories, collapse = " + ")))),
      p(population),
      div(class = "population-tags", lapply(tags_value, function(value) span(value))),
      if (nzchar(doi) && !is.na(doi)) {
        tags$a(paste0("DOI ", doi), href = paste0("https://doi.org/", doi), target = "_blank", rel = "noopener noreferrer")
      },
      div(class = "domain-reminder", tx(
        "Vérifiez que le patient, la matrice et les covariables restent dans le domaine de l'article source.",
        "Check that the patient, biological matrix, and covariates remain within the source article's domain."
      ))
    )
  })
  shiny::outputOptions(output, "model_context_ui", suspendWhenHidden = FALSE)

  output$pmetrics_context_ui <- renderUI({
    artifact <- current_pmetrics_artifact()
    english <- identical(current_language(), "en")
    div(
      class = "model-context",
      div(class = "model-context-head", tags$strong(if (english) artifact$labelEn else artifact$label), span("Pmetrics · POSTPROB")),
      p(if (english) artifact$descriptionEn else artifact$description),
      div(class = "population-tags",
        span(tx("Points de support pondérés", "Weighted support points")),
        span(tx("Artefact de démonstration", "Demonstration artifact"))
      ),
      div(class = "domain-reminder", tx(
        "Distribution synthétique : ne pas utiliser pour une décision clinique.",
        "Synthetic distribution: do not use for a clinical decision."
      ))
    )
  })
  shiny::outputOptions(output, "pmetrics_context_ui", suspendWhenHidden = FALSE)

  output$ml_status_ui <- renderUI({
    if (!identical(input$model_source %||% "library", "library")) {
      return(div(class = "ml-status", tx("ML : indisponible pour un modèle de session.", "ML: unavailable for a session model.")))
    }
    status <- ml_status_summary(selected_model_ids(), input$administration_route %||% "", current_administration_mode())
    localized(tagList(
      div(class = if (status$available) "ml-status available" else "ml-status", ml_status_message(status, current_language())),
      if (status$available) {
        checkboxInput(
          "enable_experimental_ml",
          tx("Afficher l'estimation AUC24 ML expérimentale", "Show the experimental ML AUC24 estimate"),
          value = isTRUE(isolate(input$enable_experimental_ml))
        )
      }
    ))
  })
  shiny::outputOptions(output, "ml_status_ui", suspendWhenHidden = FALSE)

  covariate_definition <- reactive({
    if (identical(input$model_source, "custom")) return(parse_covariates(input$custom_code %||% DEFAULT_CODE))
    if (identical(input$model_source, "pmetrics")) {
      return(data.frame(name = character(), value = numeric(), description = character(), scope = character()))
    }
    primary_id <- current_model_id()
    ids <- if (isTRUE(input$enable_averaging)) {
      unique(c(primary_id, input$average_model_ids %||% character()))
    } else {
      primary_id
    }
    ids <- intersect(ids, MODEL_CATALOG$id)
    if (!length(ids)) ids <- DEFAULT_MODEL
    definitions <- lapply(ids, function(id) parse_covariates(read_library_code(id)))
    definitions <- Filter(nrow, definitions)
    if (!length(definitions)) return(data.frame(name = character(), value = numeric(), description = character(), scope = character()))
    combined <- do.call(rbind, definitions)
    combined[!duplicated(combined$name), , drop = FALSE]
  })

  patient_covariate_definition <- reactive({
    definition <- covariate_definition()
    definition[definition$scope != "administration", , drop = FALSE]
  })

  administration_covariate_definition <- reactive({
    definition <- covariate_definition()
    definition[definition$scope == "administration", , drop = FALSE]
  })

  register_dose_covariates <- function(index) {
    local({
      row_index <- index
      output_id <- paste0("dose_covariates_", row_index)
      output[[output_id]] <- renderUI({
        definition <- administration_covariate_definition()
        if (!nrow(definition)) return(NULL)
        tagList(lapply(seq_len(nrow(definition)), function(definition_index) {
          name <- definition$name[[definition_index]]
          input_id <- paste0("dose_cov_", name, "_", row_index)
          pending <- pending_lego_covariates()
          current <- if (!is.null(pending) && identical(input$custom_code, pending$code)) NULL else isolate(input[[input_id]])
          value <- suppressWarnings(as.numeric(current))
          if (length(value) != 1L || !is.finite(value)) value <- definition$value[[definition_index]]
          numericInput(input_id, paste0(name, " · ", tx("à cette administration", "at this administration")), value = value)
        }))
      })
      shiny::outputOptions(output, output_id, suspendWhenHidden = FALSE)
    })
  }

  register_observation_covariates <- function(index) {
    local({
      row_index <- index
      output_id <- paste0("observation_covariates_", row_index)
      output[[output_id]] <- renderUI({
        definition <- patient_covariate_definition()
        if (!nrow(definition)) {
          return(div(class = "empty-state compact", tx("Aucune covariable pour ce modèle.", "This model has no covariates.")))
        }

        tagList(lapply(seq_len(nrow(definition)), function(definition_index) {
          name <- definition$name[[definition_index]]
          input_id <- paste0("observation_cov_", name, "_", row_index)
          pending <- pending_lego_covariates()
          current <- if (!is.null(pending) && identical(input$custom_code, pending$code)) NULL else isolate(input[[input_id]])
          value <- suppressWarnings(as.numeric(current))
          if (length(value) != 1L || !is.finite(value)) value <- definition$value[[definition_index]]
          numericInput(
            input_id,
            if (nzchar(definition$description[[definition_index]])) {
              paste0(name, " · ", definition$description[[definition_index]])
            } else {
              name
            },
            value = value
          )
        }))
      })
      shiny::outputOptions(output, output_id, suspendWhenHidden = FALSE)
    })
  }

  register_dose_validation(1L)
  register_dose_covariates(1L)
  register_observation_covariates(1L)
  register_observation_validation(1L)

  observe({
    pending <- pending_lego_covariates()
    shiny::req(pending)
    shiny::req(identical(input$model_source, "custom"))
    shiny::req(identical(input$custom_code, pending$code))
    definition <- pending$definition
    patient_definition <- definition[definition$scope != "administration", , drop = FALSE]
    administration_definition <- definition[definition$scope == "administration", , drop = FALSE]
    input_ids <- c(
      unlist(lapply(seq_len(observation_count()), function(index) paste0("observation_cov_", patient_definition$name, "_", index))),
      unlist(lapply(seq_len(dose_count()), function(index) paste0("dose_cov_", administration_definition$name, "_", index)))
    )
    if (length(input_ids)) shiny::req(all(vapply(input_ids, function(id) !is.null(input[[id]]), logical(1))))

    for (index in seq_len(observation_count())) {
      for (definition_index in seq_len(nrow(patient_definition))) {
        updateNumericInput(
          session,
          paste0("observation_cov_", patient_definition$name[[definition_index]], "_", index),
          value = patient_definition$value[[definition_index]]
        )
      }
    }
    for (index in seq_len(dose_count())) {
      for (definition_index in seq_len(nrow(administration_definition))) {
        updateNumericInput(session, paste0("dose_cov_", administration_definition$name[[definition_index]], "_", index), value = administration_definition$value[[definition_index]])
      }
    }
    pending_lego_covariates(NULL)
  })

  output$observation_covariate_notice <- renderUI({
    definition <- patient_covariate_definition()
    if (!nrow(definition)) {
      return(div(class = "covariate-notice empty", tx("Ce modèle ne déclare aucune covariable.", "This model declares no covariates.")))
    }
    div(
      class = "covariate-notice",
      tags$strong(tx("Covariables au moment du prélèvement", "Covariates at sampling time")),
      span(tx(
        "Renseignez leur valeur sur chaque ligne. La dernière valeur connue est conservée jusqu'au prélèvement suivant.",
        "Enter their values on each row. The latest known value is carried forward until the next sample."
      ))
    )
  })
  shiny::outputOptions(output, "observation_covariate_notice", suspendWhenHidden = FALSE)

  parse_clock_minutes <- function(value, label) {
    value <- trimws(as.character(value %||% ""))
    match <- regexec("^([01]?[0-9]|2[0-3]):([0-5][0-9])$", value)
    parts <- regmatches(value, match)[[1]]
    shiny::validate(shiny::need(length(parts) == 3L, paste0(label, tx(" doit respecter HH:MM.", " must use HH:MM."))))
    as.numeric(parts[[2]]) * 60 + as.numeric(parts[[3]])
  }

  reference_datetime <- function() {
    date <- as.Date(input$reference_date %||% Sys.Date())
    minutes <- parse_clock_minutes(input$reference_clock %||% "00:00", tx("L'heure de référence", "Reference time"))
    as.POSIXct(date, tz = "UTC") + minutes * 60
  }

  parse_calendar_datetime <- function(value, reference, label) {
    value <- trimws(as.character(value %||% ""))
    has_year <- grepl("^\\d{1,2}/\\d{1,2}/\\d{2,4}", value) || grepl("^\\d{4}-", value)
    year <- format(reference, "%Y", tz = "UTC")
    candidates <- if (has_year) value else paste(value, year)
    formats <- if (has_year) {
      c("%d/%m/%Y %H:%M", "%d/%m/%y %H:%M", "%Y-%m-%d %H:%M")
    } else {
      "%d/%m %H:%M %Y"
    }
    parsed <- NA_real_
    for (format in formats) {
      candidate <- suppressWarnings(as.POSIXct(strptime(candidates, format = format, tz = "UTC")))
      if (!is.na(candidate)) {
        parsed <- candidate
        break
      }
    }
    shiny::validate(shiny::need(!is.na(parsed), paste0(label, tx(
      " est invalide. Utilisez JJ/MM HH:MM ou JJ/MM/AAAA HH:MM.",
      " is invalid. Use DD/MM HH:MM, DD/MM/YYYY HH:MM, or YYYY-MM-DD HH:MM."
    ))))
    if (!has_year && parsed > reference + 24 * 3600) {
      previous <- paste(value, as.integer(year) - 1L)
      parsed <- as.POSIXct(strptime(previous, format = "%d/%m %H:%M %Y", tz = "UTC"))
    }
    parsed
  }

  read_record_time <- function(prefix, index) {
    mode <- isolate(input$time_entry_mode %||% "relative_hours")
    if (identical(mode, "relative_hours")) {
      return(as.numeric(isolate(input[[paste0(prefix, "_time_", index)]])))
    }
    reference <- reference_datetime()
    if (identical(mode, "relative_days")) {
      days <- as.numeric(isolate(input[[paste0(prefix, "_days_ago_", index)]]))
      clock <- parse_clock_minutes(isolate(input[[paste0(prefix, "_clock_", index)]]), tx("L'heure saisie", "Entered time"))
      reference_clock <- as.numeric(format(reference, "%H", tz = "UTC")) * 60 + as.numeric(format(reference, "%M", tz = "UTC"))
      return(-days * 24 + (clock - reference_clock) / 60)
    }
    parsed <- parse_calendar_datetime(
      isolate(input[[paste0(prefix, "_calendar_", index)]]),
      reference,
      tx(paste0("La date de la ligne ", index), paste0("Date on row ", index))
    )
    as.numeric(difftime(parsed, reference, units = "hours"))
  }

  normalize_timeline <- function(doses, observation_records) {
    origin <- min(doses$time, na.rm = TRUE)
    doses$time <- doses$time - origin
    observation_records$raw$time <- observation_records$raw$time - origin
    observation_records$covariate_history$time <- observation_records$covariate_history$time - origin
    observation_records$observations$time <- observation_records$observations$time - origin
    shiny::validate(shiny::need(all(observation_records$raw$time >= 0), tx(
      "Une concentration ou covariable est datée avant la première administration.",
      "A concentration or covariate is dated before the first administration."
    )))
    list(doses = doses, observations = observation_records, origin = origin)
  }

  read_doses <- function(include_missed = FALSE) {
    administration_definition <- isolate(administration_covariate_definition())
    rows <- lapply(seq_len(isolate(dose_count())), function(index) {
      steady_state <- isTRUE(isolate(input[[paste0("dose_ss_", index)]]))
      status <- as.character(isolate(input[[paste0("dose_status_", index)]]) %||% "administered")
      count <- if (steady_state) 1L else as.integer(isolate(input[[paste0("dose_count_", index)]]))
      row <- data.frame(
        time = if (steady_state) 0 else read_record_time("dose", index),
        amount = as.numeric(isolate(input[[paste0("dose_amount_", index)]])),
        interval = as.numeric(isolate(input[[paste0("dose_interval_", index)]])),
        count = count,
        infusion = if (identical(isolate(input$administration_route), "IV")) as.numeric(isolate(input[[paste0("dose_infusion_", index)]])) else 0,
        ss = as.integer(steady_state),
        status = status,
        time_uncertainty = if (identical(status, "uncertain")) as.numeric(isolate(input[[paste0("dose_time_uncertainty_", index)]]) %||% 0) else 0,
        stringsAsFactors = FALSE
      )
      for (definition_index in seq_len(nrow(administration_definition))) {
        name <- administration_definition$name[[definition_index]]
        candidate <- suppressWarnings(as.numeric(isolate(input[[paste0("dose_cov_", name, "_", index)]])))
        if (length(candidate) != 1L || !is.finite(candidate)) candidate <- administration_definition$value[[definition_index]]
        row[[name]] <- candidate
      }
      row
    })
    data <- do.call(rbind, rows)
    numeric_columns <- c("time", "amount", "interval", "count", "infusion", "ss", "time_uncertainty")
    shiny::validate(shiny::need(all(is.finite(as.matrix(data[, numeric_columns, drop = FALSE]))), tx("Toutes les administrations doivent être numériques.", "All administration values must be numeric.")))
    valid <- data$amount > 0 & data$interval >= 0 & data$count >= 1 & data$infusion >= 0 &
      data$ss %in% c(0L, 1L) & (data$ss == 0L | data$interval > 0) &
      data$status %in% c("administered", "uncertain", "missed") & data$time_uncertainty >= 0
    shiny::validate(shiny::need(all(valid), tx("Administration invalide. Un steady state nécessite un intervalle strictement positif.", "Invalid administration. Steady state requires a strictly positive interval.")))
    if (isTRUE(include_missed)) return(data)
    data <- data[data$status != "missed", , drop = FALSE]
    shiny::validate(shiny::need(nrow(data) > 0, tx("Au moins une dose administrée ou incertaine est requise.", "At least one administered or uncertain dose is required.")))
    data
  }

  read_observation_records <- function() {
    definition <- isolate(patient_covariate_definition())
    rows <- lapply(seq_len(isolate(observation_count())), function(index) {
      blq <- isTRUE(isolate(input[[paste0("observation_blq_", index)]]))
      row <- data.frame(
        time = read_record_time("observation", index),
        concentration = as.numeric(isolate(input[[paste0("observation_concentration_", index)]])),
        blq = blq,
        lloq = if (blq) as.numeric(isolate(input[[paste0("observation_lloq_", index)]]) %||% 0) else 0,
        matrix = as.character(isolate(input[[paste0("observation_matrix_", index)]]) %||% "unspecified"),
        time_uncertainty = if (isTRUE(isolate(input[[paste0("observation_time_uncertain_", index)]]))) as.numeric(isolate(input[[paste0("observation_time_uncertainty_", index)]]) %||% 0) else 0,
        stringsAsFactors = FALSE
      )
      for (definition_index in seq_len(nrow(definition))) {
        name <- definition$name[[definition_index]]
        candidate <- suppressWarnings(as.numeric(isolate(input[[paste0("observation_cov_", name, "_", index)]])))
        if (length(candidate) != 1L || !is.finite(candidate)) candidate <- definition$value[[definition_index]]
        row[[name]] <- candidate
      }
      row
    })
    data <- do.call(rbind, rows)
    shiny::validate(shiny::need(all(is.finite(data$time)), tx("Le temps de prélèvement doit être valide.", "Sample time must be valid.")))
    concentration_valid <- (!data$blq & is.finite(data$concentration) & data$concentration >= 0) |
      (data$blq & is.finite(data$lloq) & data$lloq > 0)
    shiny::validate(shiny::need(all(concentration_valid), tx("Une concentration BLQ exige une LLOQ positive; les autres concentrations doivent être positives.", "A BLQ concentration requires a positive LLOQ; other concentrations must be non-negative.")))
    shiny::validate(shiny::need(all(data$matrix %in% c("unspecified", "plasma", "serum", "whole_blood")), tx("Matrice biologique invalide.", "Invalid biological matrix.")))
    shiny::validate(shiny::need(all(is.finite(data$time_uncertainty) & data$time_uncertainty >= 0), tx("Incertitude horaire invalide.", "Invalid time uncertainty.")))
    data <- data[order(data$time), , drop = FALSE]

    covariate_names <- definition$name
    covariate_history <- data[, c("time", covariate_names), drop = FALSE]
    baseline <- if (length(covariate_names)) {
      as.list(data[1, covariate_names, drop = FALSE])
    } else {
      list()
    }
    usable <- !data$blq & data$concentration > 0
    observations <- data[usable, c("time", "concentration", "time_uncertainty"), drop = FALSE]
    blq_method <- isolate(input$blq_method %||% "exclude")
    if (identical(blq_method, "lloq_half") && any(data$blq)) {
      blq_rows <- data[data$blq, c("time", "lloq", "time_uncertainty"), drop = FALSE]
      names(blq_rows)[names(blq_rows) == "lloq"] <- "concentration"
      blq_rows$concentration <- blq_rows$concentration / 2
      observations <- rbind(observations, blq_rows)
      observations <- observations[order(observations$time), , drop = FALSE]
    }
    list(observations = observations, covariate_history = covariate_history, baseline = baseline, raw = data)
  }

  output$data_quality_preview <- renderUI({
    reactiveValuesToList(input)
    dose_messages <- unlist(lapply(seq_len(dose_count()), dose_row_messages), use.names = FALSE)
    observation_messages <- unlist(lapply(seq_len(observation_count()), observation_row_messages), use.names = FALSE)
    messages <- unique(c(dose_messages, observation_messages))
    messages <- messages[nzchar(messages)]
    if (!length(messages)) {
      return(localized(div(class = "quality-summary ok", tags$strong("Saisie cohérente"), span("Aucune anomalie de forme détectée avant l'analyse."))))
    }
    localized(div(
      class = "quality-summary warning",
      tags$strong(tx(paste(length(messages), "point(s) à vérifier"), paste(length(messages), "item(s) to check"))),
      tags$ul(lapply(head(messages, 8), tags$li)),
      if (length(messages) > 8) span(tx(paste0("+ ", length(messages) - 8, " autre(s) message(s) dans les lignes."), paste0("+ ", length(messages) - 8, " other row message(s).")))
    ))
  })
  shiny::outputOptions(output, "data_quality_preview", suspendWhenHidden = TRUE)

  output$therapy_timeline <- renderPlot({
    reactiveValuesToList(input)
    preview <- tryCatch({
      doses <- read_doses(include_missed = TRUE)
      records <- read_observation_records()
      origin_candidates <- doses$time[doses$status != "missed"]
      if (!length(origin_candidates)) stop("Ajoutez une administration reçue.")
      origin <- min(origin_candidates, na.rm = TRUE)
      doses$time <- doses$time - origin
      records$raw$time <- records$raw$time - origin
      list(doses = doses, observations = records$raw)
    }, error = function(error) NULL)
    shiny::validate(shiny::need(!is.null(preview), tx("La frise apparaîtra lorsque les horaires seront cohérents.", "The timeline will appear once the times are consistent.")))

    dose_events <- do.call(rbind, lapply(seq_len(nrow(preview$doses)), function(index) {
      row <- preview$doses[index, , drop = FALSE]
      repetitions <- if (row$ss[[1]] == 1L || row$status[[1]] == "missed") 0 else max(0L, row$count[[1]] - 1L)
      data.frame(
        time = row$time[[1]] + seq.int(0L, repetitions) * row$interval[[1]],
        lane = tx("Administrations", "Administrations"),
        type = row$status[[1]],
        label = if (row$ss[[1]] == 1L) paste0(row$amount[[1]], " mg · SS") else paste0(row$amount[[1]], " mg"),
        stringsAsFactors = FALSE
      )
    }))
    samples <- data.frame(
      time = preview$observations$time,
      lane = tx("Prélèvements", "Samples"),
      type = ifelse(preview$observations$blq, "blq", "sample"),
      label = ifelse(preview$observations$blq, paste0("BLQ < ", preview$observations$lloq), as.character(preview$observations$concentration)),
      stringsAsFactors = FALSE
    )
    events <- rbind(dose_events, samples)
    ggplot(events, aes(time, lane, color = type, shape = type)) +
      geom_hline(yintercept = c(1, 2), color = "#dfe4e8", linewidth = 0.5) +
      geom_point(size = 3) +
      geom_text(aes(label = label), nudge_y = 0.16, check_overlap = TRUE, size = 3, color = "#17202a") +
      scale_color_manual(values = c(
        administered = "#176b70", uncertain = "#b16916", missed = "#a43d35",
        sample = "#a4441f", blq = "#66717c"
      ), drop = FALSE) +
      scale_shape_manual(values = c(administered = 16, uncertain = 17, missed = 4, sample = 15, blq = 1), drop = FALSE) +
      labs(x = tx("Temps depuis la première administration reçue (h)", "Time since first received administration (h)"), y = NULL, title = tx("Frise thérapeutique", "Treatment timeline"), color = NULL, shape = NULL) +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank(), legend.position = "bottom", plot.title = element_text(face = "bold", size = 12))
  })
  shiny::outputOptions(output, "therapy_timeline", suspendWhenHidden = TRUE)

  build_quality_diagnostics <- function(dose_records, observation_records, model_ids, blq_method) {
    english <- identical(current_language(), "en")
    messages <- c(
      if (any(dose_records$status == "missed")) tx(paste(sum(dose_records$status == "missed"), "dose(s) oubliée(s) exclue(s) du calcul."), paste(sum(dose_records$status == "missed"), "missed dose(s) excluded from calculations.")),
      if (any(dose_records$status == "uncertain")) tx(paste(sum(dose_records$status == "uncertain"), "dose(s) avec horaire incertain."), paste(sum(dose_records$status == "uncertain"), "dose(s) with uncertain time.")),
      if (any(observation_records$time_uncertainty > 0)) tx(paste(sum(observation_records$time_uncertainty > 0), "prélèvement(s) avec incertitude horaire."), paste(sum(observation_records$time_uncertainty > 0), "sample(s) with timing uncertainty.")),
      if (any(observation_records$blq)) {
        tx(
          paste(sum(observation_records$blq), "valeur(s) BLQ", if (identical(blq_method, "exclude")) "exclue(s) de l'ajustement." else "imputée(s) à LLOQ/2 (méthode exploratoire)."),
          paste(sum(observation_records$blq), "BLQ value(s)", if (identical(blq_method, "exclude")) "excluded from the fit." else "imputed as LLOQ/2 (exploratory method).")
        )
      },
      if (any(observation_records$matrix == "unspecified")) tx("Matrice biologique non précisée pour au moins un prélèvement.", "Biological matrix is unspecified for at least one sample."),
      if (length(unique(observation_records$matrix)) > 1L) tx("Plusieurs matrices biologiques sont mélangées sans conversion.", "Several biological matrices are mixed without conversion.")
    )

    bounds <- list(
      WT = c(2, 300), WEIGHT = c(2, 300), AGE = c(0, 110), BMI = c(10, 80),
      CRCL = c(0, 250), CLCR = c(0, 250), EGFR = c(0, 250), CREAT = c(5, 2000),
      SCR = c(5, 2000), ALB = c(5, 60), ALBUMIN = c(5, 60), HT = c(40, 230), HEIGHT = c(40, 230)
    )
    covariate_names <- intersect(names(observation_records), names(bounds))
    for (name in covariate_names) {
      limits <- bounds[[name]]
      values <- suppressWarnings(as.numeric(observation_records[[name]]))
      if (any(is.finite(values) & (values < limits[[1]] | values > limits[[2]]))) {
        messages <- c(messages, tx(
          paste0(name, " sort de la plage physiologique de contrôle ", limits[[1]], "–", limits[[2]], "."),
          paste0(name, " is outside the physiological screening range ", limits[[1]], "–", limits[[2]], ".")
        ))
      }
    }

    populations <- lapply(model_ids, function(id) {
      record <- model_record(id)
      list(
        model = paste(if (english) record$drugEn[[1]] else record$drug[[1]], record$model[[1]], sep = " - "),
        population = if (english) record$populationEn[[1]] else record$population[[1]],
        tags = if (english) record$populationTagsEn[[1]] %||% character() else record$populationTags[[1]] %||% character(),
        doi = record$doi[[1]] %||% "",
        model_type = if (english) record$modelTypeEn[[1]] %||% "" else record$modelType[[1]] %||% "",
        note = if (english) record$noteEn[[1]] %||% "" else record$note[[1]] %||% ""
      )
    })
    engineering <- Filter(function(item) grepl("déterministe|ingénierie|adapté|deterministic|engineering|adapted", item$model_type, ignore.case = TRUE), populations)
    if (length(engineering)) {
      messages <- c(messages, tx("Au moins un modèle est adapté ou déterministe : relisez sa note avant interprétation.", "At least one model is adapted or deterministic: review its note before interpretation."))
    }
    list(messages = unique(messages[nzchar(messages)]), populations = populations)
  }

  build_fit_diagnostics <- function(fits) {
    lapply(fits, function(fit) {
      diagnostic <- fit_numeric_diagnostics(fit)
      messages <- character()
      if (identical(diagnostic$status, "error")) {
        messages <- tx("Ajustement non disponible.", "Fit unavailable.")
      } else if (identical(diagnostic$status, "population")) {
        messages <- tx("Aucune concentration exploitable : prédiction populationnelle, sans estimation individuelle ni incertitude postérieure.", "No usable concentration: population prediction without individual estimation or posterior uncertainty.")
      } else if (identical(diagnostic$engine, "pmetrics")) {
        messages <- if (is.finite(diagnostic$effective_support)) tx(
          paste0("Posterior Pmetrics calculé sur les points de support pondérés; nombre effectif de points : ", format_metric(diagnostic$effective_support), "."),
          paste0("Pmetrics posterior computed over weighted support points; effective number of points: ", format_metric(diagnostic$effective_support), ".")
        ) else tx("Posterior Pmetrics calculé; dispersion des points de support indisponible.", "Pmetrics posterior computed; support-point dispersion unavailable.")
      } else {
        if (isFALSE(diagnostic$converged)) messages <- c(messages, tx("L'optimiseur n'a pas signalé une convergence normale.", "The optimizer did not report normal convergence."))
        if (is.na(diagnostic$converged)) messages <- c(messages, tx("Statut de convergence de l'optimiseur indisponible.", "Optimizer convergence status unavailable."))
        if (isTRUE(diagnostic$sparse)) messages <- c(messages, tx(
          paste0(diagnostic$observations, " concentration(s) pour ", diagnostic$eta_count, " ETA estimés : les données individuelles sont peu informatives par rapport à la dimension du posterior."),
          paste0(diagnostic$observations, " concentration(s) for ", diagnostic$eta_count, " estimated ETAs: individual data are sparse relative to posterior dimension.")
        ))
        if (isTRUE(diagnostic$at_boundary)) messages <- c(messages, tx("Au moins un ETA est proche d'une borne numérique de l'optimisation mapbayr.", "At least one ETA is near a numerical mapbayr optimization bound."))
        if (!isTRUE(diagnostic$uncertainty_available)) messages <- c(messages, tx("Covariance postérieure indisponible : aucun intervalle d'incertitude individuel ne doit être interprété.", "Posterior covariance unavailable: no individual uncertainty interval should be interpreted."))
        if (isFALSE(diagnostic$uncertainty_stable)) messages <- c(messages, tx("Covariance postérieure non définie positive ou numériquement mal conditionnée : incertitude individuelle instable.", "Posterior covariance is not positive definite or is numerically ill-conditioned: individual uncertainty is unstable."))
        if (!length(messages)) messages <- tx("Convergence rapportée, aucun ETA proche d'une borne et covariance postérieure numériquement exploitable.", "Convergence reported, no ETA near a bound, and posterior covariance numerically usable.")
      }
      list(label = diagnostic$label, messages = messages)
    })
  }

  fit_diagnostics_ui <- function(items) {
    rows <- lapply(items, function(item) tags$li(tags$strong(item$label), tags$ul(lapply(item$messages, tags$li))))
    div(
      class = "analysis-diagnostics",
      tags$strong(tx("Diagnostics numériques de l'ajustement", "Numerical fit diagnostics")),
      tags$ul(rows),
      p(tx("Ces signaux décrivent l'optimisation et l'information disponible; ils ne démontrent ni l'identifiabilité du modèle ni sa validité clinique.", "These signals describe optimization and available information; they demonstrate neither model identifiability nor clinical validity."))
    )
  }

  analysis_provenance <- function(specifications) {
    package_version <- function(name) {
      tryCatch(as.character(utils::packageVersion(name)), error = function(error) "non disponible")
    }
    models <- lapply(specifications, function(specification) {
      if (identical(specification$engine %||% "mapbayr", "pmetrics")) {
        path <- pmetrics_artifact_path(specification$artifact_id %||% specification$id)
        hash <- if (requireNamespace("digest", quietly = TRUE)) digest::digest(file = path, algo = "sha256") else short_hash(paste(readLines(path, warn = FALSE), collapse = "\n"))
        list(id = specification$id, engine = "pmetrics", sha256 = hash)
      } else if (specification$id %in% MODEL_CATALOG$id) {
        list(id = specification$id, sha256 = model_sha256(specification$id))
      } else {
        code <- specification$code %||% ""
        hash <- if (requireNamespace("digest", quietly = TRUE)) digest::digest(code, algo = "sha256", serialize = FALSE) else short_hash(code)
        list(id = "custom", sha256 = hash)
      }
    })
    list(
      generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
      r = R.version.string,
      packages = list(mrgsolve = package_version("mrgsolve"), mapbayr = package_version("mapbayr"),
        Pmetrics = package_version("Pmetrics"), shiny = package_version("shiny")),
      models = models,
      assumptions = if (identical(current_language(), "en")) c(
        "Concentration units match the model; no automatic conversion.",
        "Covariates are carried forward to the next dated value.",
        "BLQ values are handled according to session settings.",
        "Exact domain ranges must be checked in each source article."
      ) else c(
        "Unités de concentration identiques à celles du modèle; aucune conversion automatique.",
        "Covariables reportées jusqu'à la valeur datée suivante.",
        "Valeurs BLQ traitées selon le réglage de session.",
        "Les plages de domaine exactes doivent être vérifiées dans chaque article source."
      )
    )
  }

  replace_dose_rows <- function(data) {
    removeUI(selector = "#dose_rows .dose-row", multiple = TRUE, immediate = TRUE)
    dose_count(nrow(data))
    for (index in seq_len(nrow(data))) {
      insert_dose_row(
        index,
        time = data$time[[index]],
        amount = data$amount[[index]],
        interval = data$interval[[index]],
        count = data$count[[index]],
        infusion = data$infusion[[index]],
        steady_state = isTRUE(data$ss[[index]] == 1),
        status = as.character(data$status[[index]] %||% "administered"),
        time_uncertainty = as.numeric(data$time_uncertainty[[index]] %||% 0)
      )
    }
    pending_import_dose_covariates(data)
  }

  replace_observation_rows <- function(data, expected_model_ids = character()) {
    removeUI(selector = "#observation_rows .observation-row", multiple = TRUE, immediate = TRUE)
    observation_count(nrow(data))
    for (index in seq_len(nrow(data))) {
      insert_observation_row(
        index,
        data$time[[index]],
        data$concentration[[index]],
        blq = isTRUE(as.logical(data$blq[[index]] %||% FALSE)),
        lloq = as.numeric(data$lloq[[index]] %||% NA_real_),
        matrix = as.character(data$matrix[[index]] %||% "unspecified"),
        time_uncertainty = as.numeric(data$time_uncertainty[[index]] %||% 0)
      )
    }
    pending_import_covariates(list(data = data, expected_model_ids = expected_model_ids))
  }

  observe({
    data <- pending_import_dose_covariates()
    shiny::req(data)
    definition <- administration_covariate_definition()
    names_to_apply <- intersect(definition$name, names(data))
    input_ids <- unlist(lapply(seq_len(nrow(data)), function(index) paste0("dose_cov_", names_to_apply, "_", index)))
    if (length(input_ids)) shiny::req(all(vapply(input_ids, function(id) !is.null(input[[id]]), logical(1))))
    for (index in seq_len(nrow(data))) {
      for (name in names_to_apply) {
        value <- suppressWarnings(as.numeric(data[[name]][[index]]))
        if (is.finite(value)) updateNumericInput(session, paste0("dose_cov_", name, "_", index), value = value)
      }
    }
    pending_import_dose_covariates(NULL)
  })

  observe({
    pending <- pending_import_covariates()
    shiny::req(pending)
    expected_model_ids <- pending$expected_model_ids %||% character()
    if (length(expected_model_ids)) {
      shiny::req(all(expected_model_ids %in% selected_model_ids()))
    }
    data <- pending$data
    definition <- patient_covariate_definition()
    observation_fields <- c("time", "concentration", "blq", "lloq", "matrix", "time_uncertainty")
    covariate_names <- intersect(setdiff(names(data), observation_fields), definition$name)
    input_ids <- unlist(lapply(seq_len(nrow(data)), function(index) {
      paste0("observation_cov_", covariate_names, "_", index)
    }))
    if (length(input_ids)) shiny::req(all(vapply(input_ids, function(id) !is.null(input[[id]]), logical(1))))

    for (index in seq_len(nrow(data))) {
      for (name in covariate_names) {
        value <- suppressWarnings(as.numeric(data[[name]][[index]]))
        if (is.finite(value)) {
          updateNumericInput(session, paste0("observation_cov_", name, "_", index), value = value)
        }
      }
    }
    pending_import_covariates(NULL)
    showNotification(tx("Données patient importées dans cette session.", "Patient data imported into this session."), type = "message", duration = 4)
  })

  validate_patient_table <- function(value, required, label, maximum_rows = 100L, character_columns = character()) {
    if (is.null(value)) stop(tx(paste0("Section ", label, " absente du fichier."), paste0("Missing ", label, " section in the file.")))
    data <- as.data.frame(value, stringsAsFactors = FALSE, check.names = FALSE)
    missing <- setdiff(required, names(data))
    if (length(missing)) stop(tx(
      paste0("Colonnes manquantes dans ", label, " : ", paste(missing, collapse = ", ")),
      paste0("Missing columns in ", label, ": ", paste(missing, collapse = ", "))
    ))
    if (!nrow(data) || nrow(data) > maximum_rows) stop(tx(
      paste0(label, " doit contenir entre 1 et ", maximum_rows, " lignes."),
      paste0(label, " must contain between 1 and ", maximum_rows, " rows.")
    ))
    if (ncol(data) > 50L || any(!grepl("^[A-Za-z][A-Za-z0-9_]*$", names(data)))) {
      stop(tx(paste0("Noms ou nombre de colonnes invalides dans ", label, "."), paste0("Invalid column names or count in ", label, ".")))
    }
    identity_pattern <- "(^id$|patient|name|nom|prenom|birth|naissance|email|phone|telephone|address|adresse|mrn|ipp|nir)"
    if (any(grepl(identity_pattern, names(data), ignore.case = TRUE))) {
      stop(tx("Le fichier ne doit contenir aucun identifiant patient.", "The file must not contain any patient identifier."))
    }
    numeric_columns <- setdiff(names(data), character_columns)
    for (name in numeric_columns) data[[name]] <- suppressWarnings(as.numeric(data[[name]]))
    if (length(numeric_columns) && !all(is.finite(as.matrix(data[, numeric_columns, drop = FALSE])))) {
      stop(tx(paste0("Toutes les valeurs numériques de ", label, " doivent être valides."), paste0("All numeric values in ", label, " must be valid.")))
    }
    data
  }

  output$download_patient <- downloadHandler(
    filename = function() paste0("tdm-patient-", Sys.Date(), ".json"),
    content = function(file) {
      timeline <- normalize_timeline(read_doses(), read_observation_records())
      doses <- read_doses(include_missed = TRUE)
      doses$time <- doses$time - timeline$origin
      observations <- timeline$observations$raw
      model_source <- isolate(input$model_source %||% "library")
      document <- list(
        schema = "pk-mipd-patient",
        version = 6L,
        privacy = list(containsIdentity = FALSE, customCodeIncluded = FALSE),
        timeline = list(mode = "relative_hours", origin = "first_administration"),
        model = list(
          source = model_source,
          id = if (identical(model_source, "library")) isolate(input$model_id) else if (identical(model_source, "pmetrics")) isolate(input$pmetrics_artifact_id) else NULL,
          route = isolate(input$administration_route),
          mode = isolate(current_administration_mode()),
          averaging = list(
            enabled = identical(model_source, "library") && isTRUE(isolate(input$enable_averaging)),
            ids = if (identical(model_source, "library")) isolate(input$average_model_ids %||% character()) else character(),
            scheme = isolate(input$weighting_scheme %||% "AIC")
          )
        ),
        doses = doses,
        observations = observations,
        target = list(
          metric = isolate(input$target_metric),
          low = isolate(input$target_low),
          high = isolate(input$target_high),
          concentrationThreshold = isolate(input$target_concentration_threshold),
          cminOperator = isolate(input$target_cmin_operator %||% "range"),
          cminLow = isolate(input$target_cmin_low),
          cminHigh = isolate(input$target_cmin_high),
          cminThreshold = isolate(input$target_cmin_threshold),
          cmaxOperator = isolate(input$target_cmax_operator %||% "range"),
          cmaxLow = isolate(input$target_cmax_low),
          cmaxHigh = isolate(input$target_cmax_high),
          cmaxThreshold = isolate(input$target_cmax_threshold),
          doseMin = isolate(input$dose_min),
          doseMax = isolate(input$dose_max),
          doseStep = isolate(input$dose_step),
          intervals = candidate_interval_values(),
          infusion = if (identical(isolate(input$administration_route), "IV")) isolate(input$future_infusion) else 0
        ),
        settings = list(
          delta = as.numeric(isolate(input$simulation_delta)),
          additionalDoses = as.integer(isolate(input$additional_doses)),
          monteCarloReplicates = as.integer(isolate(input$mc_replicates)),
          predictionInterval = as.numeric(isolate(input$prediction_interval)),
          variability = isolate(input$variability_components %||% character()),
          posteriorReplicates = as.integer(isolate(input$posterior_replicates %||% 20)),
          blqMethod = isolate(input$blq_method %||% "exclude"),
          residualErrorMode = isolate(input$residual_error_mode %||% "model"),
          fixedResidualCv = as.numeric(isolate(input$fixed_residual_cv %||% 1)),
          experimentalML = isTRUE(isolate(input$enable_experimental_ml))
        )
      )
      jsonlite::write_json(document, file, pretty = TRUE, auto_unbox = TRUE, dataframe = "rows", null = "null")
    }
  )

  observeEvent(input$upload_patient, {
    file_info <- input$upload_patient
    if (is.null(file_info)) return()
    tryCatch({
      if (file_info$size[[1]] > 1024 * 1024) stop(tx("Le fichier patient est limité à 1 Mo.", "The patient file is limited to 1 MB."))
      document <- jsonlite::fromJSON(file_info$datapath[[1]], simplifyDataFrame = TRUE)
      unlink(file_info$datapath[[1]], force = TRUE)
      version <- as.integer(document$version %||% 0)
      if (!identical(document$schema %||% "", "pk-mipd-patient") || !version %in% c(1L, 2L, 3L, 4L, 5L, 6L)) {
        stop(tx("Format de fichier patient non reconnu.", "Unrecognized patient file format."))
      }

      doses <- validate_patient_table(
        document$doses,
        c("time", "amount", "interval", "count", "infusion"),
        "administrations",
        character_columns = "status"
      )
      observations <- validate_patient_table(
        document$observations,
        c("time", "concentration"),
        "observations",
        character_columns = "matrix"
      )
      if (!"ss" %in% names(doses)) doses$ss <- 0
      if (!"status" %in% names(doses)) doses$status <- "administered"
      if (!"time_uncertainty" %in% names(doses)) doses$time_uncertainty <- 0
      if (!"blq" %in% names(observations)) observations$blq <- 0
      if (!"lloq" %in% names(observations)) observations$lloq <- 0
      if (!"matrix" %in% names(observations)) observations$matrix <- "unspecified"
      if (!"time_uncertainty" %in% names(observations)) observations$time_uncertainty <- 0
      invalid_doses <- doses$time < 0 | doses$amount <= 0 | doses$interval < 0 | doses$count < 1 | doses$infusion < 0 |
        !doses$ss %in% c(0, 1) | (doses$ss == 1 & doses$interval <= 0) |
        !doses$status %in% c("administered", "uncertain", "missed") | doses$time_uncertainty < 0
      if (any(invalid_doses)) stop(tx(
        "Administration invalide dans le fichier. Un steady state nécessite ss = 1 et un intervalle positif.",
        "Invalid administration in the file. Steady state requires ss = 1 and a positive interval."
      ))
      invalid_observations <- observations$time < 0 | observations$concentration < 0 |
        !observations$matrix %in% c("unspecified", "plasma", "serum", "whole_blood") |
        observations$time_uncertainty < 0 | (as.logical(observations$blq) & observations$lloq <= 0)
      if (any(invalid_observations)) stop(tx("Observation invalide dans le fichier.", "Invalid observation in the file."))
      doses$count <- as.integer(doses$count)
      doses$ss <- as.integer(doses$ss)

      target <- document$target %||% list()
      model <- document$model %||% list()
      expected_model_ids <- character()
      if (identical(model$source %||% "", "library") && (model$id %||% "") %in% MODEL_CATALOG$id) {
        updateRadioButtons(session, "model_source", selected = "library")
        updateSelectInput(session, "model_id", selected = model$id)
        record <- model_record(model$id)
        routes <- model_routes(record)
        route <- as.character(model$route %||% routes[[1]])
        if (!route %in% routes) route <- routes[[1]]
        modes <- model_administration_modes(record, route)
        mode <- as.character(model$mode %||% modes[[1]])
        if (!mode %in% modes) mode <- modes[[1]]
        imported_target_key <- paste("library", model$id, route, mode, sep = "::")
        if (!identical(isolate(target_preset_key()), imported_target_key)) {
          pending_import_target(list(key = imported_target_key, target = target))
        } else {
          pending_import_target(NULL)
        }
        averaging <- model$averaging %||% list()
        updateCheckboxInput(session, "enable_averaging", value = isTRUE(averaging$enabled))
        if ((averaging$scheme %||% "AIC") %in% c("AIC", "LL")) {
          updateSelectInput(session, "weighting_scheme", selected = averaging$scheme)
        }
        model_drug <- record$drug[[1]]
        route_compatible <- vapply(seq_len(nrow(MODEL_CATALOG)), function(index) {
          MODEL_CATALOG$drug[[index]] == model_drug && model_supports_administration_mode(MODEL_CATALOG[index, , drop = FALSE], route, mode)
        }, logical(1))
        allowed_ids <- MODEL_CATALOG$id[route_compatible]
        selected_ids <- intersect(averaging$ids %||% character(), allowed_ids)
        expected_model_ids <- if (isTRUE(averaging$enabled)) unique(c(model$id, selected_ids)) else model$id
        session$onFlushed(function() {
          updateSelectInput(session, "administration_route", selected = route)
          session$onFlushed(function() {
            updateSelectInput(session, "administration_mode", selected = mode)
            updateCheckboxGroupInput(session, "average_model_ids", selected = selected_ids)
          }, once = TRUE)
        }, once = TRUE)
      } else if (identical(model$source %||% "", "pmetrics") && (model$id %||% "") %in% PMETRICS_CATALOG$id) {
        artifact <- read_pmetrics_artifact(model$id)
        updateRadioButtons(session, "model_source", selected = "pmetrics")
        updateSelectInput(session, "pmetrics_artifact_id", selected = artifact$id)
        pending_import_target(NULL)
        session$onFlushed(function() {
          updateSelectInput(session, "administration_route", selected = artifact$route)
        }, once = TRUE)
      } else if (identical(model$source %||% "", "custom")) {
        showNotification(tx(
          "Le code C++ personnalisé n'est jamais inclus dans l'export; recollez-le avant l'analyse.",
          "Custom C++ code is never included in the export; paste it again before analysis."
        ), type = "warning", duration = 8)
      }

      update_target_inputs(target)

      settings <- document$settings %||% list()
      if ((settings$delta %||% "") %in% c(0.05, 0.1, 0.25)) updateSelectInput(session, "simulation_delta", selected = as.character(settings$delta))
      additional_doses <- suppressWarnings(as.integer(settings$additionalDoses %||% NA_integer_))
      if (is.finite(additional_doses) && additional_doses >= 1 && additional_doses <= 30) updateNumericInput(session, "additional_doses", value = additional_doses)
      mc_replicates <- suppressWarnings(as.integer(settings$monteCarloReplicates %||% NA_integer_))
      if (is.finite(mc_replicates) && mc_replicates >= 50 && mc_replicates <= 1000) updateNumericInput(session, "mc_replicates", value = mc_replicates)
      prediction_interval <- suppressWarnings(as.numeric(settings$predictionInterval %||% NA_real_))
      if (prediction_interval %in% c(80, 90, 95)) updateSelectInput(session, "prediction_interval", selected = as.character(prediction_interval))
      variability <- intersect(settings$variability %||% character(), c("POSTERIOR", "RESID", "TIMING", "IIV"))
      variability[variability == "IIV"] <- "POSTERIOR"
      if (length(variability)) updateCheckboxGroupInput(session, "variability_components", selected = variability)
      posterior_replicates <- suppressWarnings(as.integer(settings$posteriorReplicates %||% NA_integer_))
      if (is.finite(posterior_replicates) && posterior_replicates >= 10 && posterior_replicates <= 40) {
        updateNumericInput(session, "posterior_replicates", value = posterior_replicates)
      }
      if ((settings$blqMethod %||% "") %in% c("exclude", "lloq_half")) {
        updateSelectInput(session, "blq_method", selected = settings$blqMethod)
      }
      if ((settings$residualErrorMode %||% "") %in% c("model", "fixed_cv")) {
        updateSelectInput(session, "residual_error_mode", selected = settings$residualErrorMode)
      }
      fixed_residual_cv <- suppressWarnings(as.numeric(settings$fixedResidualCv %||% NA_real_))
      if (is.finite(fixed_residual_cv) && fixed_residual_cv > 0 && fixed_residual_cv <= 100) {
        updateNumericInput(session, "fixed_residual_cv", value = fixed_residual_cv)
      }
      updateCheckboxInput(session, "enable_experimental_ml", value = isTRUE(settings$experimentalML))

      updateSelectInput(session, "time_entry_mode", selected = "relative_hours")
      replace_dose_rows(doses)
      replace_observation_rows(observations, expected_model_ids)
      analysis_store(NULL)
      validation_store(NULL)
    }, error = function(error) {
      pending_import_target(NULL)
      if (file.exists(file_info$datapath[[1]])) unlink(file_info$datapath[[1]], force = TRUE)
      showNotification(conditionMessage(error), type = "error", duration = 10)
    })
  }, ignoreInit = TRUE)

  model_specifications <- function() {
    route <- isolate(input$administration_route %||% "")
    if (identical(isolate(input$model_source), "custom")) {
      if (!nzchar(route)) route <- "IV"
      mode <- if (identical(route, "Oral")) "ORAL" else if (identical(route, "IM")) "IM" else "IV_INTERMITTENT"
      return(list(list(id = "custom", label = "Modèle personnalisé", code = isolate(input$custom_code), route = route, mode = mode, adm_cmt_name = NULL)))
    }
    if (identical(isolate(input$model_source), "pmetrics")) {
      artifact <- read_pmetrics_artifact(isolate(input$pmetrics_artifact_id %||% DEFAULT_PMETRICS_ARTIFACT))
      return(list(list(
        id = artifact$id,
        artifact_id = artifact$id,
        engine = "pmetrics",
        label = if (identical(current_language(), "en")) artifact$labelEn else artifact$label,
        code = NULL,
        route = artifact$route,
        mode = artifact$mode,
        adm_cmt_name = NULL
      )))
    }
    primary_id <- isolate(input$model_id)
    primary_record <- model_record(primary_id)
    if (!nzchar(route)) route <- model_routes(primary_record)[[1]]
    ids <- if (!isTRUE(isolate(input$enable_averaging))) {
      primary_id
    } else {
      unique(c(primary_id, isolate(input$average_model_ids %||% character())))
    }
    mode <- isolate(current_administration_mode())
    compatible <- vapply(ids, function(id) model_supports_administration_mode(model_record(id), route, mode), logical(1))
    shiny::validate(shiny::need(all(compatible), tx("Tous les modèles moyennés doivent utiliser la même voie et le même mode d'administration.", "All averaged models must use the same administration route and mode.")))
    lapply(ids, function(id) {
      record <- model_record(id)
      list(
        id = id,
        label = record$label[[1]],
        code = NULL,
        route = route,
        mode = mode,
        adm_cmt_name = model_administration_cmt(record, route)
      )
    })
  }

  validate_primary_model <- function() {
    if (identical(isolate(input$model_source), "pmetrics")) {
      artifact <- read_pmetrics_artifact(isolate(input$pmetrics_artifact_id %||% DEFAULT_PMETRICS_ARTIFACT))
      build_pmetrics_model(artifact, pmetrics_error_spec(
        artifact,
        isolate(input$residual_error_mode %||% "model"),
        isolate(input$fixed_residual_cv %||% 1)
      ))
      compile_pmetrics_bridge(artifact, session_model_dir, session_model_cache)
      compartments <- pmetrics_compartment_contract(artifact)
      return(list(
        ok = TRUE, errors = character(), warnings = character(),
        engine = "pmetrics", route = artifact$route, mode = artifact$mode,
        adm_cmt = compartments$adm_cmt, obs_cmt = compartments$obs_cmt, n_eta = 0L,
        n_sigma = length(artifact$error$coefficients), n_support = nrow(artifact$support)
      ))
    }
    if (identical(isolate(input$model_source), "custom")) {
      model <- compile_model(
        custom_code = isolate(input$custom_code),
        allow_custom = ALLOW_CUSTOM_MODELS,
        custom_soloc = session_model_dir,
        custom_cache = session_model_cache
      )
    } else {
      model <- compile_model(model_id = isolate(input$model_id))
    }
    model <- apply_residual_error_setting(
      model,
      isolate(input$residual_error_mode %||% "model"),
      isolate(input$fixed_residual_cv %||% 1)
    )
    contract <- validate_model_contract(model)
    specification <- model_specifications()[[1]]
    if (isTRUE(contract$ok)) contract$adm_cmt <- resolve_administration_cmt(model, contract, specification)
    contract$route <- specification$route
    contract
  }

  observeEvent(input$validate_model, {
    withProgress(message = tx("Compilation du modèle", "Compiling model"), value = 0.3, {
      result <- tryCatch(validate_primary_model(), error = function(error) list(ok = FALSE, errors = conditionMessage(error), warnings = character()))
      validation_store(result)
      incProgress(0.7)
    })
  })

  output$engine_status <- renderUI({
    result <- validation_store()
    if (is.null(result)) return(span(class = "status-pill idle", tx("Non validé", "Not validated")))
    if (isTRUE(result$ok)) span(class = "status-pill ok", tx("Modèle valide", "Valid model")) else span(class = "status-pill error", tx("Erreur modèle", "Model error"))
  })

  output$library_code <- renderText({
    read_library_code(current_model_id())
  })

  output$pmetrics_artifact_code <- renderText({
    path <- pmetrics_artifact_path(current_pmetrics_artifact()$id)
    paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  })

  output$model_contract <- renderUI({
    result <- validation_store()
    if (is.null(result)) return(div(class = "contract-note", tx(
      "Cliquez sur « Valider le modèle » pour vérifier les tags, OMEGA, SIGMA et sorties capturées.",
      "Select ‘Validate model’ to check tags, OMEGA, SIGMA, and captured outputs."
    )))
    if (!isTRUE(result$ok)) {
      return(div(class = "contract-error", tags$strong(tx("Modèle incompatible", "Incompatible model")), tags$ul(lapply(result$errors, tags$li))))
    }
    if (identical(result$engine %||% "mapbayr", "pmetrics")) {
      return(div(
        class = "contract-ok",
        tags$strong(tx("Contrat Pmetrics valide", "Valid Pmetrics contract")),
        span(tx(
          paste0("Voie ", result$route, " · ", result$n_support, " points de support · erreur polynomiale à ", result$n_sigma, " coefficient(s)"),
          paste0("Route ", result$route, " · ", result$n_support, " support points · polynomial error with ", result$n_sigma, " coefficient(s)")
        ))
      ))
    }
    div(
      class = "contract-ok",
      tags$strong(tx("Contrat mapbayr valide", "Valid mapbayr contract")),
      span(tx(
        paste0("Voie ", result$route, " · Administration CMT ", result$adm_cmt, " · Observation CMT ", result$obs_cmt, " · ", result$n_eta, " ETA · ", result$n_sigma, " erreurs résiduelles"),
        paste0("Route ", result$route, " · Administration CMT ", result$adm_cmt, " · Observation CMT ", result$obs_cmt, " · ", result$n_eta, " ETA · ", result$n_sigma, " residual errors")
      )),
      if (length(result$warnings)) tags$ul(lapply(result$warnings, tags$li))
    )
  })

  observeEvent(input$run_analysis, {
    analysis_stage <- tx("validation initiale", "initial validation")
    tryCatch({
      shiny::validate(shiny::need(isTRUE(input$accept_disclaimer), tx("Lisez et acceptez l'avertissement avant de lancer l'analyse.", "Read and accept the warning before running the analysis.")))
      analysis_stage <- tx("lecture de toutes les administrations", "reading all administrations")
      all_dose_records <- read_doses(include_missed = TRUE)
      analysis_stage <- tx("lecture des administrations reçues", "reading received administrations")
      administered_doses <- read_doses()
      analysis_stage <- tx("lecture des concentrations et covariables", "reading concentrations and covariates")
      entered_observations <- read_observation_records()
      analysis_stage <- tx("normalisation de la chronologie", "normalizing the timeline")
      timeline <- normalize_timeline(administered_doses, entered_observations)
      doses <- timeline$doses
      all_dose_records$time <- all_dose_records$time - timeline$origin
      observation_records <- timeline$observations
      observations <- observation_records$observations
      covariates <- observation_records$baseline
      administration_names <- intersect(administration_covariate_definition()$name, names(doses))
      if (length(administration_names)) {
        covariates <- c(covariates, as.list(doses[nrow(doses), administration_names, drop = FALSE]))
      }
      specifications <- model_specifications()
      library_model_ids <- intersect(vapply(specifications, `[[`, character(1), "id"), MODEL_CATALOG$id)
      intervals <- candidate_interval_values()
      route <- isolate(input$administration_route)
      infusion <- if (identical(route, "IV")) as.numeric(isolate(input$future_infusion)) else 0
      delta <- as.numeric(isolate(input$simulation_delta %||% 0.1))
      variability <- isolate(input$variability_components %||% character())
      residual_error_mode <- isolate(input$residual_error_mode %||% "model")
      fixed_residual_cv <- as.numeric(isolate(input$fixed_residual_cv %||% 1))
      target_metric <- isolate(input$target_metric %||% "AUC24")
      cmin_operator <- isolate(input$target_cmin_operator %||% "range")
      cmax_operator <- isolate(input$target_cmax_operator %||% "range")
      target_low <- if (identical(target_metric, "CminCmax")) {
        if (identical(cmin_operator, "range")) as.numeric(isolate(input$target_cmin_low)) else as.numeric(isolate(input$target_cmin_threshold))
      } else {
        as.numeric(isolate(input$target_low))
      }
      target_high <- if (identical(target_metric, "CminCmax")) {
        if (identical(cmin_operator, "range")) as.numeric(isolate(input$target_cmin_high)) else target_low
      } else {
        as.numeric(isolate(input$target_high))
      }
      cmax_low <- if (identical(cmax_operator, "range")) as.numeric(isolate(input$target_cmax_low)) else as.numeric(isolate(input$target_cmax_threshold))
      cmax_high <- if (identical(cmax_operator, "range")) as.numeric(isolate(input$target_cmax_high)) else cmax_low
      target_options <- list(
        concentration_threshold = as.numeric(isolate(input$target_concentration_threshold %||% 1)),
        cmin_operator = cmin_operator,
        cmax_operator = cmax_operator,
        cmax_low = cmax_low,
        cmax_high = cmax_high
      )

      analysis_stage <- tx("validation des réglages", "validating settings")
      shiny::validate(shiny::need(length(intervals), tx("Sélectionnez au moins un intervalle de dose.", "Select at least one dosing interval.")))
      shiny::validate(shiny::need(target_metric %in% c("AUC24", "Cmin", "Cmax", "TimeAbove", "CminCmax"), tx("La métrique cible est invalide.", "The target metric is invalid.")))
      shiny::validate(shiny::need(is.finite(target_low) && is.finite(target_high) && (
        identical(target_metric, "CminCmax") || target_high > target_low
      ), tx("La borne haute doit dépasser la borne basse.", "The upper bound must exceed the lower bound.")))
      shiny::validate(shiny::need(!identical(target_metric, "TimeAbove") || (
        is.finite(target_options$concentration_threshold) && target_options$concentration_threshold >= 0 && target_low >= 0 && target_high <= 100
      ), tx("La cible temps > seuil exige un seuil positif ou nul et des bornes comprises entre 0 et 100 %.", "The time-above-threshold target requires a non-negative threshold and bounds between 0 and 100%.")))
      shiny::validate(shiny::need(!identical(target_metric, "CminCmax") || (
        cmin_operator %in% c("range", "above", "below") && cmax_operator %in% c("range", "above", "below") &&
          target_low >= 0 && target_high >= 0 && is.finite(cmax_low) && is.finite(cmax_high) && cmax_low >= 0 && cmax_high >= 0 &&
          (!identical(cmin_operator, "range") || target_high > target_low) &&
          (!identical(cmax_operator, "range") || cmax_high > cmax_low)
      ), tx("Les contraintes combinées Cmin/Cmax sont invalides.", "The combined Cmin/Cmax constraints are invalid.")))
      shiny::validate(shiny::need(input$dose_max >= input$dose_min && input$dose_step > 0, tx("La grille de doses est invalide.", "The dose grid is invalid.")))
      shiny::validate(shiny::need(residual_error_mode %in% c("model", "fixed_cv"), tx("Le réglage de l'erreur résiduelle est invalide.", "The residual error setting is invalid.")))
      shiny::validate(shiny::need(!identical(residual_error_mode, "fixed_cv") || (is.finite(fixed_residual_cv) && fixed_residual_cv > 0 && fixed_residual_cv <= 100), tx("Le CV résiduel fixé doit être compris entre 0 et 100 %.", "The fixed residual CV must be between 0 and 100%.")))

      result <- withProgress(message = tx("Analyse pharmacométrique", "Pharmacometric analysis"), value = 0, {
        analysis_stage <- tx("ajustement MAP", "MAP fitting")
        incProgress(0.12, detail = tx("Compilation et validation", "Compilation and validation"))
        fits <- fit_model_set(
          specifications,
          doses,
          observations,
          covariates,
          ALLOW_CUSTOM_MODELS,
          covariate_history = observation_records$covariate_history,
          custom_soloc = session_model_dir,
          custom_cache = session_model_cache,
          residual_error_mode = residual_error_mode,
          fixed_residual_cv = fixed_residual_cv
        )
        fits <- apply_hybrid_ml_to_fits(
          fits,
          route,
          enabled = isTRUE(isolate(input$enable_experimental_ml))
        )
        valid <- successful_fits(fits)
        if (!length(valid)) {
          messages <- vapply(fits, function(item) item$message %||% "Unknown error", character(1))
          stop(paste(messages, collapse = " | "))
        }

        analysis_stage <- tx("exposition actuelle", "current exposure")
        incProgress(0.32, detail = tx("Estimation MAP", "MAP estimation"))
        weights <- compute_model_weights(fits, isolate(input$weighting_scheme %||% "AIC"))
        dose_ss <- if ("ss" %in% names(doses)) as.integer(doses$ss) else rep(0L, nrow(doses))
        known_dose_times <- doses$time + ifelse(dose_ss == 1L, 0, doses$interval * pmax(0, doses$count - 1))
        end_time <- max(c(known_dose_times, observations$time, 24), na.rm = TRUE) + 48
        profiles <- fit_profiles(fits, weights, end_time, delta = max(delta, 0.1))
        current_exposure <- current_regimen_exposure(fits, weights, doses, observations)

        analysis_stage <- tx("grille de posologies", "dosage grid")
        incProgress(0.58, detail = tx("Exploration des posologies", "Exploring dosage regimens"))
        recommendations <- recommend_regimens(
          fits = fits,
          weights = weights,
          dose_min = isolate(input$dose_min),
          dose_max = isolate(input$dose_max),
          dose_step = isolate(input$dose_step),
          intervals = intervals,
          infusion = infusion,
          metric = target_metric,
          target_low = target_low,
          target_high = target_high,
          target_options = target_options,
          delta = delta
        )
        analysis_stage <- tx("probabilité d'atteindre la cible", "target attainment probability")
        recommendations <- rank_regimens_by_pta(
          recommendations = recommendations,
          fits = fits,
          weights = weights,
          metric = target_metric,
          target_low = target_low,
          target_high = target_high,
          target_options = target_options,
          replicates = min(150L, as.integer(isolate(input$mc_replicates %||% 250))),
          delta = delta,
          include_posterior = "POSTERIOR" %in% variability,
          include_residual = "RESID" %in% variability
        )
        best <- recommendations[1, , drop = FALSE]
        best_profile <- simulate_averaged_regimen(
          fits,
          weights,
          dose = best$dose[[1]],
          interval = best$interval[[1]],
          infusion = infusion,
          delta = delta,
          concentration_threshold = target_options$concentration_threshold
        )$profile

        analysis_stage <- tx("sensibilité du model averaging", "model averaging sensitivity")
        averaging_sensitivity <- model_averaging_sensitivity(
          fits = fits,
          selected_weights = weights,
          dose_min = isolate(input$dose_min),
          dose_max = isolate(input$dose_max),
          dose_step = isolate(input$dose_step),
          intervals = intervals,
          infusion = infusion,
          metric = target_metric,
          target_low = target_low,
          target_high = target_high,
          target_options = target_options,
          delta = delta
        )

        analysis_stage <- tx("comparaison des doses futures", "future dose comparison")
        incProgress(0.76, detail = "Comparaison des doses futures")
        future_comparison <- compare_future_regimens(
          fits,
          weights,
          doses,
          observations,
          recommended = best,
          additional_doses = isolate(input$additional_doses %||% 6),
          delta = delta
        )

        analysis_stage <- tx("distribution prédictive", "predictive distribution")
        incProgress(0.86, detail = tx("Distribution prédictive", "Predictive distribution"))
        distribution <- simulate_regimen_distribution(
          fits,
          weights,
          dose = best$dose[[1]],
          interval = best$interval[[1]],
          infusion = infusion,
          metric = target_metric,
          replicates = isolate(input$mc_replicates %||% 250),
          interval_level = as.numeric(isolate(input$prediction_interval %||% 90)),
          delta = delta,
          include_posterior = "POSTERIOR" %in% variability,
          include_residual = "RESID" %in% variability,
          include_timing = "TIMING" %in% variability,
          timing_refits = as.integer(isolate(input$posterior_replicates %||% 20)),
          target_low = target_low,
          target_high = target_high,
          target_options = target_options
        )

        analysis_stage <- tx("qualité et traçabilité", "quality and traceability")
        incProgress(0.9, detail = tx("Préparation des résultats", "Preparing results"))
        quality <- build_quality_diagnostics(
          all_dose_records,
          observation_records$raw,
          library_model_ids,
          isolate(input$blq_method %||% "exclude")
        )
        list(
          fits = fits,
          weights = weights,
          model_summary = model_summary(fits, weights),
          fit_profiles = profiles,
          current_exposure = current_exposure,
          recommendations = recommendations,
          averaging_sensitivity = averaging_sensitivity,
          best = best,
          best_profile = best_profile,
          future_comparison = future_comparison,
          distribution = distribution,
          doses = doses,
          all_dose_records = all_dose_records,
          observations = observations,
          observation_records = observation_records$raw,
          route = route,
          infusion = infusion,
          target_metric = target_metric,
          target_low = target_low,
          target_high = target_high,
          target_options = target_options,
          target_preset = applied_target_preset(),
          weighting_scheme = isolate(input$weighting_scheme %||% "AIC"),
          quality = quality,
          fit_diagnostics = build_fit_diagnostics(fits),
          provenance = analysis_provenance(specifications),
          ml_status = ml_application_summary(fits, weights),
          settings = list(
            delta = delta,
            additional_doses = as.integer(isolate(input$additional_doses %||% 6)),
            mc_replicates = as.integer(isolate(input$mc_replicates %||% 250)),
            posterior_replicates = as.integer(isolate(input$posterior_replicates %||% 20)),
            prediction_interval = as.numeric(isolate(input$prediction_interval %||% 90)),
            variability = variability,
            blq_method = isolate(input$blq_method %||% "exclude"),
            residual_error_mode = residual_error_mode,
            fixed_residual_cv = if (identical(residual_error_mode, "fixed_cv")) fixed_residual_cv else NA_real_,
            show_component_profiles = isTRUE(isolate(input$show_component_profiles))
          )
        )
      })
      analysis_store(result)
      bslib::nav_select("analysis_tabs", "fit")
      showNotification(tx("Analyse terminée.", "Analysis completed."), type = "message", duration = 3)
    }, error = function(error) {
      detail <- conditionMessage(error)
      if (!nzchar(detail)) detail <- tx("Une entrée réactive requise n'est pas disponible.", "A required reactive input is unavailable.")
      message("TDM analysis error [", analysis_stage, "] class=", paste(class(error), collapse = "/"))
      showNotification(paste0(analysis_stage, " : ", detail), type = "error", duration = 10)
    })
  })

  output$fit_empty <- renderUI({
    if (is.null(analysis_store())) div(class = "empty-results", tx("Renseignez les données puis lancez l'analyse.", "Enter the data, then run the analysis."))
  })

  output$analysis_warnings <- renderUI({
    result <- analysis_store()
    if (is.null(result)) return(NULL)
    population_rows <- lapply(result$quality$populations, function(item) {
      tags$li(tags$strong(item$model), paste0(" : ", item$population))
    })
    tagList(
      div(
        class = "analysis-diagnostics",
        tags$strong(tx("Applicabilité et qualité des données", "Applicability and data quality")),
        if (length(result$quality$messages)) tags$ul(lapply(result$quality$messages, tags$li)) else span(tx("Aucune anomalie générique détectée.", "No generic issue detected.")),
        if (length(population_rows)) tagList(span(tx("Populations sources à comparer au patient :", "Source populations to compare with the patient:")), tags$ul(population_rows)),
        if (identical(result$settings$residual_error_mode, "fixed_cv")) {
          div(
            class = "residual-setting-warning",
            tx(
              paste0("Ajustement MAP avec un CV résiduel proportionnel fixé à ", format_metric(result$settings$fixed_residual_cv), " %; représentation SIGMA adaptée à l'échelle d'erreur du modèle."),
              paste0("MAP fit with proportional residual CV fixed at ", format_metric(result$settings$fixed_residual_cv), "%; SIGMA representation adapted to the model error scale.")
            )
          )
        },
        span(class = "ml-result-status", ml_status_message(result$ml_status, current_language())),
        ml_domain_warning_block(result$ml_status, current_language())
      ),
      fit_diagnostics_ui(result$fit_diagnostics)
    )
  })

  output$current_exposure_summary <- renderUI({
    result <- analysis_store()
    if (is.null(result)) return(NULL)
    exposure <- result$current_exposure
    regimen_label <- if (isTRUE(exposure$steady_state)) {
      tx(paste0(format_metric(exposure$dose), " mg / ", format_metric(exposure$interval), " h · État stationnaire"), paste0(format_metric(exposure$dose), " mg / ", format_metric(exposure$interval), " h · Steady state"))
    } else if (isTRUE(exposure$single_dose)) {
      tx(paste0("Dose unique · ", format_metric(exposure$dose), " mg"), paste0("Single dose · ", format_metric(exposure$dose), " mg"))
    } else {
      paste0(format_metric(exposure$dose), " mg / ", format_metric(exposure$interval), " h")
    }
    c0_detail <- if (isTRUE(exposure$single_dose)) {
      tx("À 24 h après la dose", "At 24 h after dose")
    } else {
      tx(paste0("À ", format_metric(exposure$interval), " h, avant la dose suivante"), paste0("At ", format_metric(exposure$interval), " h, before the next dose"))
    }
    pmetrics_mapped <- any(vapply(successful_fits(result$fits), function(fit) {
      identical(fit$engine %||% "mapbayr", "pmetrics") && isTRUE(fit$pmetrics$mapped)
    }, logical(1)))
    method <- if (pmetrics_mapped) {
      tx("MAP non paramétrique · Pmetrics", "Nonparametric MAP · Pmetrics")
    } else if (any(vapply(successful_fits(result$fits), function(fit) !is.null(fit$estimate), logical(1)))) {
      tx("MAP bayésienne · mapbayr", "Bayesian MAP · mapbayr")
    } else {
      tx("Prédiction populationnelle", "Population prediction")
    }
    partial_auc_detail <- if (exposure$historical_coverage_hours < 24 - 1e-8) {
      tx("Intégrale MAP partielle, non directement comparable à une AUC24", "Partial MAP integral, not directly comparable with an AUC24")
    } else {
      tx("Fenêtre MAP complète de 24 h", "Complete 24-hour MAP window")
    }
    div(
      class = "exposure-strip",
      div(
        span(tx("AUC24 MAP · état stationnaire", "MAP AUC24 · steady state")),
        strong(format_metric(exposure$steady_state_auc24)),
        tags$small(tx("Horizon de 24 h comparable au ML", "24-hour horizon comparable with ML"))
      ),
      if (is.finite(result$ml_status$auc24 %||% NA_real_)) {
        div(
          span(tx("AUC24 ML expérimentale", "Experimental ML AUC24")),
          strong(format_metric(result$ml_status$auc24)),
          tags$small(tx("Estimation directe issue des prélèvements récents · non utilisée pour la recommandation", "Direct estimate from recent samples · not used for the recommendation"))
        )
      },
      div(
        span(tx(paste0("AUC MAP · fenêtre ", format_metric(exposure$historical_coverage_hours), " h"), paste0("MAP AUC · ", format_metric(exposure$historical_coverage_hours), "-hour window"))),
        strong(format_metric(exposure$historical_auc24)),
        tags$small(partial_auc_detail)
      ),
      div(
        span(tx("C0 actuelle", "Current C0")),
        strong(format_metric(exposure$historical_c0)),
        tags$small(tx(paste0("MAP juste avant la prochaine dose à h ", format_metric(exposure$historical_c0_time)), paste0("MAP immediately before the next dose at h ", format_metric(exposure$historical_c0_time))))
      ),
      div(span(tx("C0 à l'état stationnaire", "Steady-state C0")), strong(format_metric(exposure$steady_state_c0)), tags$small(c0_detail)),
      div(
        span(tx("Schéma actuel", "Current regimen")),
        strong(regimen_label),
        tags$small(if (exposure$infusion > 0) {
          tx(paste0("Perfusion IV : ", format_metric(exposure$infusion), " h"), paste0("IV infusion: ", format_metric(exposure$infusion), " h"))
        } else if (identical(result$route, "Oral")) {
          tx("Administration orale · perfusion = 0 h", "Oral administration · infusion = 0 h")
        } else if (identical(result$route, "IM")) {
          tx("Administration intramusculaire · absorption définie par le modèle", "Intramuscular administration · absorption defined by the model")
        } else {
          tx("Bolus IV · perfusion = 0 h", "IV bolus · infusion = 0 h")
        })
      ),
      div(span(tx("Estimation", "Estimate")), strong(method), tags$small(if (length(result$weights) > 1) tx("Prédiction pondérée par model averaging", "Prediction weighted by model averaging") else tx("Modèle individuel sélectionné", "Selected individual model")))
    )
  })

  output$fit_plot <- renderPlot({
    result <- analysis_store()
    shiny::req(result)
    profiles <- result$fit_profiles
    shiny::validate(shiny::need(nrow(profiles$average), tx("Ajoutez au moins une concentration observée pour produire l'ajustement MAP.", "Add at least one observed concentration to produce the MAP fit.")))

    plot <- ggplot()
    if (isTRUE(result$settings$show_component_profiles)) {
      plot <- plot + geom_line(
          data = profiles$per_model,
          aes(time, concentration, group = model),
          color = "#8c969f",
          linewidth = 0.6,
          alpha = 0.55
        )
    }
    plot +
      geom_line(data = profiles$average, aes(time, concentration), color = "#176b70", linewidth = 1.2) +
      geom_point(data = result$observations, aes(time, concentration), color = "#a4441f", size = 3) +
      labs(x = tx("Temps (h)", "Time (h)"), y = tx("Concentration", "Concentration"), title = tx("Prédiction individuelle et observations", "Individual prediction and observations")) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  })

  output$future_comparison_summary <- renderUI({
    result <- analysis_store()
    if (is.null(result)) return(div(class = "empty-results", tx("La comparaison sera calculée avec l'analyse.", "The comparison will be calculated with the analysis.")))
    comparison <- result$future_comparison
    current <- comparison$scenarios[[SCENARIO_MAINTAIN]]
    recommended <- comparison$scenarios[[SCENARIO_RECOMMENDED]]
    div(
      class = "comparison-strip",
      div(span(tx("Prochaine dose simulée", "Next simulated dose")), strong(paste0("h ", format_metric(comparison$future_start)))),
      div(span(tx("Poursuite", "Maintain")), strong(paste0(current$dose, " mg / ", current$interval, " h"))),
      div(span(tx("Recommandation", "Recommendation")), strong(paste0(recommended$dose, " mg / ", recommended$interval, " h"))),
      div(span(tx("Projection", "Projection")), strong(tx(paste0(comparison$additional_doses, " doses supplémentaires"), paste0(comparison$additional_doses, " additional doses"))))
    )
  })

  output$future_comparison_plot <- renderPlot({
    result <- analysis_store()
    shiny::req(result)
    comparison <- result$future_comparison
    window_start <- max(0, comparison$decision_time - 24)
    ggplot(comparison$profiles, aes(time, concentration, color = scenario, linetype = scenario)) +
      geom_vline(xintercept = comparison$decision_time, color = "#66717c", linetype = "dashed") +
      geom_vline(xintercept = comparison$future_start, color = "#17202a", linetype = "dotted") +
      geom_line(linewidth = 1.2) +
      coord_cartesian(xlim = c(window_start, comparison$end_time)) +
      scale_color_manual(values = stats::setNames(c("#3f4a54", "#176b70"), c(SCENARIO_MAINTAIN, SCENARIO_RECOMMENDED))) +
      scale_linetype_manual(values = stats::setNames(c("dashed", "solid"), c(SCENARIO_MAINTAIN, SCENARIO_RECOMMENDED))) +
      labs(x = tx("Temps depuis la première administration (h)", "Time since first administration (h)"), y = "Concentration", color = NULL, linetype = NULL) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), legend.position = "top")
  })

  output$ml_interpretation_panel <- renderUI({
    result <- analysis_store()
    if (is.null(result)) return(NULL)
    concordance <- ml_concordance(result, current_language())
    if (is.null(concordance)) return(NULL)
    explanation <- result$ml_status$explanation %||% list(available = FALSE, reason = "Explication DALEX indisponible.")
    agreement_label <- switch(
      concordance$level,
      close = tx("Estimations proches", "Close estimates"),
      moderate = tx("Écart modéré", "Moderate difference"),
      divergent = tx("Divergence importante", "Substantial divergence")
    )
    tagList(
      h3(tx("Concordance MAP-BE et ML", "MAP-BE and ML agreement")),
      ml_domain_warning_block(result$ml_status, current_language()),
      div(
        class = "comparison-strip ml-comparison-strip",
        div(span("AUC24 MAP-BE"), strong(format_metric(concordance$map_auc24)), tags$small("mrgsolve + mapbayr")),
        div(span("AUC24 ML"), strong(format_metric(concordance$ml_auc24)), tags$small("XGBoost direct")),
        div(
          span(tx("Écart relatif", "Relative difference")),
          strong(paste0(if (concordance$relative_gap >= 0) "+" else "", format_metric(concordance$relative_gap), " %")),
          tags$small(tx("ML par rapport au MAP-BE", "ML relative to MAP-BE"))
        ),
        div(span(tx("Lecture", "Interpretation")), strong(agreement_label), tags$small(concordance$interpretation))
      ),
      plotOutput("ml_comparison_plot", height = "330px"),
      h3(tx("Explication locale DALEX de l'AUC24 ML", "Local DALEX explanation of ML AUC24")),
      if (isTRUE(explanation$available)) {
        tagList(
          plotOutput("ml_explanation_plot", height = "470px"),
          div(
            class = "analysis-diagnostics ml-explanation-help",
            tags$strong(tx("Comment lire cette figure", "How to read this figure")),
            tags$ul(
              tags$li(tx("Chaque barre montre comment une information de ce patient déplace la prédiction par rapport à la référence synthétique DALEX.", "Each bar shows how one item of patient information shifts the prediction from the synthetic DALEX baseline.")),
              tags$li(tx("Les contributions expliquent le calcul de XGBoost; elles ne démontrent pas une relation causale.", "Contributions explain the XGBoost calculation; they do not establish causality.")),
              tags$li(tx("Une concordance MAP-BE/ML est un contrôle de cohérence, pas une validation clinique indépendante.", "MAP-BE/ML agreement is a consistency check, not an independent clinical validation."))
            )
          )
        )
      } else {
        div(class = "empty-state compact", tx(paste0("Explication DALEX indisponible : ", explanation$reason %||% "raison inconnue"), paste0("DALEX explanation unavailable: ", explanation$reason %||% "unknown reason")))
      }
    )
  })

  output$administration_mode_ui <- renderUI({
    if (identical(input$model_source %||% "library", "custom")) return(NULL)
    if (identical(input$model_source %||% "library", "pmetrics")) {
      return(div(class = "custom-note", paste(
        tx("Mode", "Mode"),
        administration_mode_label(current_pmetrics_artifact()$mode, current_language())
      )))
    }
    record <- model_record(current_model_id())
    route <- input$administration_route %||% ""
    shiny::req(model_supports_route(record, route))
    modes <- model_administration_modes(record, route)
    selected <- isolate(input$administration_mode %||% "")
    if (!selected %in% modes) selected <- modes[[1]]
    selectInput(
      "administration_mode",
      tx("Mode d'administration", "Administration mode"),
      choices = stats::setNames(modes, vapply(modes, administration_mode_label, character(1), lang = current_language())),
      selected = selected
    )
  })

  active_target_preset <- reactive({
    if (!identical(input$model_source %||% "library", "library")) return(NULL)
    record <- model_record(current_model_id())
    clinical_target_preset(record, current_administration_mode())
  })

  applied_target_preset <- reactive({
    preset <- active_target_preset()
    if (is.null(preset)) return(NULL)
    if (!identical(input$target_metric %||% "", preset$target_metric)) return(NULL)
    low <- suppressWarnings(as.numeric(input$target_low %||% NA_real_))
    high <- suppressWarnings(as.numeric(input$target_high %||% NA_real_))
    if (!isTRUE(all.equal(low, preset$target_low)) || !isTRUE(all.equal(high, preset$target_high))) return(NULL)
    preset
  })

  target_preset_key <- reactive({
    paste(
      input$model_source %||% "library",
      if (identical(input$model_source %||% "library", "pmetrics")) current_pmetrics_artifact()$id else current_model_id(),
      input$administration_route %||% "",
      current_administration_mode(),
      sep = "::"
    )
  })

  observeEvent(target_preset_key(), {
    key <- target_preset_key()
    pending <- pending_import_target()
    if (!is.null(pending)) {
      if (!identical(key, pending$key)) return()
      last_target_preset_key(key)
      update_target_inputs(pending$target)
      pending_import_target(NULL)
      return()
    }
    if (is.null(last_target_preset_key())) {
      last_target_preset_key(key)
      return()
    }
    if (identical(last_target_preset_key(), key)) return()
    last_target_preset_key(key)
    preset <- active_target_preset()
    if (is.null(preset)) return()
    updateSelectInput(session, "target_metric", selected = preset$target_metric)
    updateNumericInput(session, "target_low", value = preset$target_low)
    updateNumericInput(session, "target_high", value = preset$target_high)
    updateNumericInput(session, "dose_min", value = preset$dose_min)
    updateNumericInput(session, "dose_max", value = preset$dose_max)
    updateNumericInput(session, "dose_step", value = preset$dose_step)
    update_candidate_intervals(preset$intervals)
    updateNumericInput(session, "future_infusion", value = preset$infusion)
  }, ignoreInit = FALSE)

  output$target_preset_ui <- renderUI({
    preset <- active_target_preset()
    if (is.null(preset)) {
      return(div(
        class = "target-preset unavailable",
        tags$strong(tx("Aucun préréglage clinique universel", "No universal clinical preset")),
        span(tx(
          "La cible et la grille affichées doivent être définies selon l'indication, le protocole local et le patient; elles ne constituent pas une recommandation.",
          "The displayed target and grid must be defined for the indication, local protocol, and patient; they are not a recommendation."
        ))
      ))
    }
    if (is.null(applied_target_preset())) {
      return(div(
        class = "target-preset unavailable",
        tags$strong(tx("Cible personnalisée", "Custom target")),
        span(tx(
          "Les valeurs affichées diffèrent du préréglage documenté de cette molécule; elles seront utilisées telles quelles.",
          "The displayed values differ from this drug's documented preset and will be used as entered."
        ))
      ))
    }
    scope <- if (identical(current_language(), "en")) preset$scope_en else preset$scope_fr
    grid_note <- if (identical(current_language(), "en")) preset$grid_note_en else preset$grid_note_fr
    div(
      class = "target-preset",
      tags$strong(tx("Préréglage documenté appliqué", "Documented preset applied")),
      span(paste0(scope, " ", grid_note)),
      tags$a(
        paste0(preset$citation, " DOI ", preset$doi),
        href = paste0("https://doi.org/", preset$doi),
        target = "_blank",
        rel = "noopener noreferrer"
      )
    )
  })
  shiny::outputOptions(output, "target_preset_ui", suspendWhenHidden = FALSE)

  output$ml_comparison_plot <- renderPlot({
    result <- analysis_store()
    shiny::req(result)
    plot <- build_ml_comparison_plot(result, current_language())
    shiny::validate(shiny::need(!is.null(plot), tx("La comparaison MAP-BE/ML n'est pas calculable.", "The MAP-BE/ML comparison cannot be calculated.")))
    plot
  })

  output$ml_explanation_plot <- renderPlot({
    result <- analysis_store()
    shiny::req(result)
    plot <- build_ml_explanation_plot(result$ml_status$explanation, lang = current_language(), unit = result$ml_status$unit %||% "mg.h/L")
    shiny::validate(shiny::need(!is.null(plot), tx("L'explication DALEX n'est pas disponible.", "The DALEX explanation is unavailable.")))
    plot
  })

  output$averaging_sensitivity_summary <- renderUI({
    result <- analysis_store()
    if (is.null(result)) return(NULL)
    table <- result$averaging_sensitivity
    if (!nrow(table) || nrow(table) == 1L) {
      return(div(class = "empty-state compact", tx("Un seul jeu de poids distinct : aucune sensibilité inter-modèles calculable.", "Only one distinct weight set: no between-model sensitivity can be calculated.")))
    }
    target_range <- paste0(format_metric(min(table$target_value)), "–", format_metric(max(table$target_value)))
    target_label <- tx("Valeur cible", "Target value")
    if (identical(result$target_metric, "TimeAbove")) {
      target_range <- paste0(format_metric(min(table$time_above_pct)), "–", format_metric(max(table$time_above_pct)), " %")
      target_label <- target_metric_label(result$target_metric, current_language())
    } else if (identical(result$target_metric, "CminCmax")) {
      target_label <- tx("Écart maximal aux deux cibles", "Maximum distance from both targets")
    }
    div(
      class = "sensitivity-strip",
      div(span(tx("Dose proposée", "Proposed dose")), strong(paste0(min(table$dose), "–", max(table$dose), " mg"))),
      div(span(tx("Intervalles proposés", "Proposed intervals")), strong(paste(sort(unique(table$interval)), collapse = ", "), " h")),
      div(span(target_label), strong(target_range))
    )
  })

  output$averaging_sensitivity_table <- renderDT({
    result <- analysis_store()
    shiny::req(result)
    table <- result$averaging_sensitivity
    if (identical(result$target_metric, "TimeAbove")) {
      table$target_value <- NULL
    } else {
      table$time_above_pct <- NULL
    }
    rounding_columns <- intersect(c("auc24", "cmin", "cmax", "time_above_pct", "target_value"), names(table))
    table[rounding_columns] <- lapply(table[rounding_columns], round, 2)
    labels <- if (identical(current_language(), "en")) {
      c(scenario = "Analysis", dose = "Dose", interval = "Interval", auc24 = "AUC24", cmin = "Cmin", cmax = "Cmax", time_above_pct = "%T > threshold", target_value = "Target value", in_target = "In target", weights = "Weights")
    } else {
      c(scenario = "Analyse", dose = "Dose", interval = "Intervalle", auc24 = "AUC24", cmin = "Cmin", cmax = "Cmax", time_above_pct = "%T > seuil", target_value = "Valeur cible", in_target = "Dans cible", weights = "Poids")
    }
    if (identical(result$target_metric, "CminCmax")) labels[["target_value"]] <- tx("Écart cible max.", "Maximum target distance")
    datatable(
      table,
      rownames = FALSE,
      options = list(dom = "t", pageLength = nrow(table), scrollX = TRUE),
      colnames = unname(labels[names(table)])
    )
  })

  output$pmetrics_support_table <- renderDT({
    result <- analysis_store()
    shiny::req(result)
    fits <- Filter(function(fit) identical(fit$engine %||% "mapbayr", "pmetrics"), successful_fits(result$fits))
    shiny::req(length(fits) == 1L)
    fit <- fits[[1]]
    parameter_names <- names(fit$posterior_parameters)
    table <- fit$posterior_points[, c(parameter_names, "prob"), drop = FALSE]
    parameter_labels <- c(ka = "ka (1/h)", ke = "ke (1/h)", k12 = "k12 (1/h)", k21 = "k21 (1/h)", v = "V (L)")
    labels <- unname(ifelse(parameter_names %in% names(parameter_labels), parameter_labels[parameter_names], parameter_names))
    datatable(table, rownames = FALSE,
      colnames = c(labels, if (isTRUE(fit$pmetrics$mapped)) tx("Probabilité postérieure", "Posterior probability") else tx("Probabilité a priori", "Prior probability")),
      options = list(dom = "t", scrollX = TRUE, pageLength = nrow(table)))
  })

  output$model_table <- renderDT({
    result <- analysis_store()
    shiny::req(result)
    table <- result$model_summary
    table$weight <- round(table$weight, 4)
    table$clearance <- round(table$clearance, 3)
    datatable(
      table,
      rownames = FALSE,
      options = list(dom = "t", pageLength = nrow(table), scrollX = TRUE),
      colnames = app_t(current_language(), c("Modèle", "Statut", "Poids", "CL", "ETA", "AUC24 ML", "Détail"))
    )
  })

  output$recommendation_summary <- renderUI({
    result <- analysis_store()
    if (is.null(result)) return(div(class = "empty-results", tx("Aucun scénario calculé.", "No scenario calculated.")))
    best <- result$best
    div(
      class = "recommendation-strip",
      div(span(tx("Scénario le plus proche de la cible", "Scenario closest to target")), strong(tx(paste0(best$dose, " mg toutes les ", best$interval, " h"), paste0(best$dose, " mg every ", best$interval, " h")))),
      div(span(target_metric_label(result$target_metric, current_language())), strong(target_value_label(result$target_metric, best, current_language()))),
      div(span(tx("Probabilité d'atteindre la cible", "Target attainment probability")), strong(if (is.finite(best$p_target)) paste0(round(100 * best$p_target), " %") else tx("Non calculée", "Not calculated"))),
      div(span("AUC 0-24 h"), strong(round(best$auc24, 1))),
      div(span("Cmin / Cmax"), strong(paste(round(best$cmin, 1), round(best$cmax, 1), sep = " / ")))
    )
  })

  output$recommended_profile <- renderPlot({
    result <- analysis_store()
    shiny::req(result)
    metric <- result$target_metric
    profile <- result$best_profile
    if (!identical(metric, "AUC24")) profile <- profile[profile$time <= result$best$interval[[1]] + 1e-8, , drop = FALSE]
    plot <- ggplot(profile, aes(time, concentration))
    if (metric %in% c("Cmin", "Cmax")) {
      plot <- plot + annotate("rect", xmin = -Inf, xmax = Inf, ymin = result$target_low, ymax = result$target_high, fill = "#dcece5", alpha = 0.65)
    } else if (identical(metric, "CminCmax")) {
      cmin_bounds <- target_component_plot_bounds(result$target_options$cmin_operator %||% "range", result$target_low, result$target_high)
      cmax_bounds <- target_component_plot_bounds(result$target_options$cmax_operator %||% "range", result$target_options$cmax_low, result$target_options$cmax_high)
      plot <- plot +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = cmin_bounds[[1]], ymax = cmin_bounds[[2]], fill = "#dcece5", alpha = 0.55) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = cmax_bounds[[1]], ymax = cmax_bounds[[2]], fill = "#f3dfc5", alpha = 0.45)
    } else if (identical(metric, "TimeAbove")) {
      plot <- plot + geom_hline(yintercept = result$target_options$concentration_threshold, color = "#a4441f", linetype = "dashed", linewidth = 0.9)
    }
    plot +
      geom_line(color = "#176b70", linewidth = 1.2) +
      labs(
        x = tx("Temps depuis la dose (h)", "Time since dose (h)"),
        y = "Concentration",
        title = tx("Profil moyen pondéré", "Weighted mean profile"),
        subtitle = target_definition_label(metric, result$target_low, result$target_high, result$target_options, current_language())
      ) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  })

  output$distribution_summary <- renderUI({
    result <- analysis_store()
    if (is.null(result)) return(div(class = "empty-results", tx("La distribution sera calculée avec l'analyse.", "The distribution will be calculated with the analysis.")))
    distribution <- result$distribution
    components <- c(
      if (isTRUE(distribution$include_posterior)) tx("incertitude postérieure MAP", "MAP posterior uncertainty"),
      if (isTRUE(distribution$include_residual)) {
        if (identical(distribution$residual_error_mode, "fixed_cv")) {
          tx(
            paste0("erreur résiduelle fixée à ", format_metric(distribution$fixed_residual_cv), " % CV"),
            paste0("residual error fixed at ", format_metric(distribution$fixed_residual_cv), "% CV")
          )
        } else {
          tx("erreur résiduelle du modèle", "model residual error")
        }
      },
      if (isTRUE(distribution$include_timing) && isTRUE(distribution$timing_available)) tx("incertitude des horaires", "timing uncertainty"),
      if (length(result$weights) > 1L) tx("incertitude entre modèles", "between-model uncertainty")
    )
    if (identical(distribution$metric, "CminCmax")) {
      cmin <- distribution$component_quantiles$cmin
      cmax <- distribution$component_quantiles$cmax
      return(div(
        class = "distribution-strip",
        div(span(tx("Cmin médiane", "Median Cmin")), strong(format_metric(cmin[["median"]]))),
        div(span(tx(paste0("Cmin · intervalle ", distribution$interval_level, " %"), paste0("Cmin · ", distribution$interval_level, "% interval"))), strong(paste0(format_metric(cmin[["lower"]]), " – ", format_metric(cmin[["upper"]])))),
        div(span(tx("Cmax médiane", "Median Cmax")), strong(format_metric(cmax[["median"]]))),
        div(span(tx(paste0("Cmax · intervalle ", distribution$interval_level, " %"), paste0("Cmax · ", distribution$interval_level, "% interval"))), strong(paste0(format_metric(cmax[["lower"]]), " – ", format_metric(cmax[["upper"]])))),
        div(span(tx("Composantes", "Components")), strong(if (length(components)) paste(components, collapse = ", ") else tx("Aucune variabilité aléatoire", "No random variability")))
      ))
    }
    div(
      class = "distribution-strip",
      div(span(tx("Médiane", "Median")), strong(paste0(format_metric(distribution$median), if (identical(distribution$metric, "TimeAbove")) " %" else ""))),
      div(span(tx(paste0("Intervalle prédictif ", distribution$interval_level, " %"), paste0("Prediction interval ", distribution$interval_level, "%"))), strong(paste0(format_metric(distribution$lower), " – ", format_metric(distribution$upper)))),
      div(span(tx("Composantes", "Components")), strong(if (length(components)) paste(components, collapse = ", ") else tx("Aucune variabilité aléatoire", "No random variability")))
    )
  })

  output$distribution_plot <- renderPlot({
    result <- analysis_store()
    shiny::req(result)
    distribution <- result$distribution
    title <- tx(paste0("Distribution du scénario ", result$best$dose, " mg / ", result$best$interval, " h"), paste0("Scenario distribution ", result$best$dose, " mg / ", result$best$interval, " h"))
    if (identical(distribution$metric, "CminCmax")) {
      cmin_bounds <- target_component_plot_bounds(result$target_options$cmin_operator %||% "range", result$target_low, result$target_high)
      cmax_bounds <- target_component_plot_bounds(result$target_options$cmax_operator %||% "range", result$target_options$cmax_low, result$target_options$cmax_high)
      return(ggplot(distribution$data, aes(cmin, cmax)) +
        annotate("rect", xmin = cmin_bounds[[1]], xmax = cmin_bounds[[2]], ymin = cmax_bounds[[1]], ymax = cmax_bounds[[2]], fill = "#dcece5", color = "#287a4d", alpha = 0.55) +
        geom_point(aes(color = in_target), alpha = 0.55, size = 1.8) +
        scale_color_manual(values = c(`TRUE` = "#287a4d", `FALSE` = "#a4441f"), guide = "none") +
        labs(x = "Cmin", y = "Cmax", title = title, subtitle = target_definition_label(distribution$metric, result$target_low, result$target_high, result$target_options, current_language())) +
        theme_minimal(base_size = 13) +
        theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")))
    }
    metric_label <- target_metric_label(distribution$metric, current_language())
    ggplot(distribution$data, aes(value)) +
      annotate("rect", xmin = result$target_low, xmax = result$target_high, ymin = -Inf, ymax = Inf, fill = "#dcece5", alpha = 0.55) +
      geom_histogram(bins = 30, fill = "#9fc9c8", color = "#176b70", linewidth = 0.25) +
      geom_vline(xintercept = c(distribution$lower, distribution$upper), color = "#a4441f", linetype = "dashed") +
      geom_vline(xintercept = distribution$median, color = "#176b70", linewidth = 1) +
      labs(x = metric_label, y = tx("Réplications", "Replicates"), title = title) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  })

  output$recommendation_table <- renderDT({
    result <- analysis_store()
    shiny::req(result)
    table <- head(result$recommendations, 50)
    table$auc24 <- round(table$auc24, 1)
    table$cmin <- round(table$cmin, 2)
    table$cmax <- round(table$cmax, 2)
    table$time_above_pct <- round(table$time_above_pct, 2)
    table$value <- round(table$value, 2)
    table$p_under <- round(100 * table$p_under, 1)
    table$p_target <- round(100 * table$p_target, 1)
    table$p_over <- round(100 * table$p_over, 1)
    table$distance <- round(table$distance, 3)
    table$center_distance <- NULL
    table$pta_evaluated <- NULL
    if (identical(result$target_metric, "TimeAbove")) {
      table$value <- NULL
    } else {
      table$time_above_pct <- NULL
    }
    if (identical(result$target_metric, "CminCmax")) {
      table$value <- NULL
      table$p_under <- NULL
      table$p_over <- NULL
    }
    labels <- if (identical(current_language(), "en")) {
      c(dose = "Dose", interval = "Interval", infusion = "Infusion", auc24 = "AUC24", cmin = "Cmin", cmax = "Cmax", time_above_pct = "%T > threshold", value = "Target value", in_target = "In target", distance = "Distance", p_under = "Below target (%)", p_target = "In target (%)", p_over = "Above target (%)")
    } else {
      c(dose = "Dose", interval = "Intervalle", infusion = "Perfusion", auc24 = "AUC24", cmin = "Cmin", cmax = "Cmax", time_above_pct = "%T > seuil", value = "Valeur cible", in_target = "Dans cible", distance = "Distance", p_under = "Sous-cible (%)", p_target = "Dans cible (%)", p_over = "Sur-cible (%)")
    }
    if (identical(result$target_metric, "CminCmax")) labels[["distance"]] <- tx("Écart cible max.", "Maximum target distance")
    datatable(
      table,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 12, scrollX = TRUE),
      colnames = unname(labels[names(table)])
    ) |>
      formatStyle("in_target", target = "row", backgroundColor = styleEqual(c(TRUE, FALSE), c("#e8f4ed", "transparent")))
  })

  report_table <- function(data) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    tags$table(
      tags$thead(tags$tr(lapply(names(data), tags$th))),
      tags$tbody(lapply(seq_len(nrow(data)), function(index) {
        tags$tr(lapply(data[index, , drop = FALSE], function(value) tags$td(as.character(value))))
      }))
    )
  }

  report_plot_uri <- function(plot, width = 9, height = 4.5) {
    path <- tempfile(fileext = ".png")
    on.exit(unlink(path, force = TRUE), add = TRUE)
    ggsave(path, plot = plot, width = width, height = height, dpi = 150, bg = "white")
    base64enc::dataURI(file = path, mime = "image/png")
  }

  output$download_report <- downloadHandler(
    filename = function() paste0(if (identical(current_language(), "en")) "tdm-report-" else "rapport-tdm-", Sys.Date(), ".html"),
    content = function(file) {
      result <- analysis_store()
      shiny::req(result)
      exposure <- result$current_exposure
      best <- result$best
      partial_auc_detail <- if (exposure$historical_coverage_hours < 24 - 1e-8) {
        tx("Intégrale MAP partielle, non directement comparable à une AUC24", "Partial MAP integral, not directly comparable with an AUC24")
      } else {
        tx("Fenêtre MAP complète de 24 h", "Complete 24-hour MAP window")
      }
      fit_plot <- ggplot() +
        geom_line(data = result$fit_profiles$average, aes(time, concentration), color = "#176b70", linewidth = 1.1) +
        geom_point(data = result$observations, aes(time, concentration), color = "#a4441f", size = 2.5) +
        labs(x = tx("Temps (h)", "Time (h)"), y = "Concentration") +
        theme_minimal(base_size = 11)
      comparison <- result$future_comparison
      comparison_plot <- ggplot(comparison$profiles, aes(time, concentration, color = scenario, linetype = scenario)) +
        geom_vline(xintercept = comparison$decision_time, color = "#66717c", linetype = "dashed") +
        geom_line(linewidth = 1.15) +
        scale_color_manual(values = stats::setNames(c("#3f4a54", "#176b70"), c(SCENARIO_MAINTAIN, SCENARIO_RECOMMENDED))) +
        scale_linetype_manual(values = stats::setNames(c("dashed", "solid"), c(SCENARIO_MAINTAIN, SCENARIO_RECOMMENDED))) +
        labs(x = tx("Temps (h)", "Time (h)"), y = "Concentration", color = NULL, linetype = NULL) +
        theme_minimal(base_size = 11) +
        theme(legend.position = "top")
      concordance <- ml_concordance(result, current_language())
      ml_comparison_plot <- build_ml_comparison_plot(result, current_language())
      ml_explanation <- result$ml_status$explanation %||% list(available = FALSE)
      ml_explanation_plot <- build_ml_explanation_plot(ml_explanation, lang = current_language(), unit = result$ml_status$unit %||% "mg.h/L")
      target_definition <- target_definition_label(result$target_metric, result$target_low, result$target_high, result$target_options, current_language())
      distribution_summary_text <- if (identical(result$distribution$metric, "CminCmax")) {
        cmin <- result$distribution$component_quantiles$cmin
        cmax <- result$distribution$component_quantiles$cmax
        tx(
          paste0("Cmin médiane ", format_metric(cmin[["median"]]), " (intervalle prédictif ", result$distribution$interval_level, " % : ", format_metric(cmin[["lower"]]), " - ", format_metric(cmin[["upper"]]), ") ; Cmax médiane ", format_metric(cmax[["median"]]), " (", format_metric(cmax[["lower"]]), " - ", format_metric(cmax[["upper"]]), ")."),
          paste0("Median Cmin ", format_metric(cmin[["median"]]), " (", result$distribution$interval_level, "% prediction interval: ", format_metric(cmin[["lower"]]), " - ", format_metric(cmin[["upper"]]), "); median Cmax ", format_metric(cmax[["median"]]), " (", format_metric(cmax[["lower"]]), " - ", format_metric(cmax[["upper"]]), ").")
        )
      } else {
        suffix <- if (identical(result$distribution$metric, "TimeAbove")) " %" else ""
        tx(
          paste0("Médiane ", format_metric(result$distribution$median), suffix, " · intervalle prédictif ", result$distribution$interval_level, " % : ", format_metric(result$distribution$lower), " - ", format_metric(result$distribution$upper), suffix, "."),
          paste0("Median ", format_metric(result$distribution$median), suffix, " · ", result$distribution$interval_level, "% prediction interval: ", format_metric(result$distribution$lower), " - ", format_metric(result$distribution$upper), suffix, ".")
        )
      }

      model_table <- result$model_summary
      model_table$weight <- round(model_table$weight, 4)
      model_table$clearance <- round(model_table$clearance, 3)
      dose_table <- result$all_dose_records
      dose_names <- c(
        time = "Temps (h)", amount = "Dose (mg)", interval = "Intervalle (h)", count = "Nombre",
        infusion = "Perfusion (h)", ss = "SS", status = "Statut", time_uncertainty = "Incertitude horaire (± h)"
      )
      if (identical(current_language(), "en")) dose_names <- c(
        time = "Time (h)", amount = "Dose (mg)", interval = "Interval (h)", count = "Count",
        infusion = "Infusion (h)", ss = "SS", status = "Status", time_uncertainty = "Time uncertainty (± h)"
      )
      renamed_dose_columns <- unname(dose_names[names(dose_table)])
      renamed_dose_columns[is.na(renamed_dose_columns)] <- names(dose_table)[is.na(renamed_dose_columns)]
      names(dose_table) <- renamed_dose_columns
      observation_table <- result$observation_records
      names(observation_table)[names(observation_table) == "time"] <- "Temps (h)"
      names(observation_table)[names(observation_table) == "concentration"] <- "Concentration"
      if (identical(current_language(), "en")) names(observation_table)[names(observation_table) == "Temps (h)"] <- "Time (h)"
      references <- lapply(names(successful_fits(result$fits)), function(id) {
        fit <- successful_fits(result$fits)[[id]]
        if (identical(fit$engine %||% "mapbayr", "pmetrics")) {
          return(tags$li(tx("Artefact Pmetrics de démonstration · distribution synthétique", "Pmetrics demonstration artifact · synthetic distribution")))
        }
        if (!id %in% MODEL_CATALOG$id) return(tags$li(tx("Modèle personnalisé de session", "Custom session model")))
        record <- model_record(id)
        tags$li(record$citation[[1]], " · ", tags$a(paste0("DOI ", record$doi[[1]]), href = paste0("https://doi.org/", record$doi[[1]])))
      })
      if (!is.null(result$target_preset)) {
        references <- c(references, list(tags$li(
          tx("Préréglage de cible : ", "Target preset: "),
          result$target_preset$citation,
          " · ",
          tags$a(paste0("DOI ", result$target_preset$doi), href = paste0("https://doi.org/", result$target_preset$doi))
        )))
      }
      sensitivity_table <- result$averaging_sensitivity
      if (identical(result$target_metric, "TimeAbove")) {
        sensitivity_table$target_value <- NULL
      } else {
        sensitivity_table$time_above_pct <- NULL
      }
      rounding_columns <- intersect(c("auc24", "cmin", "cmax", "time_above_pct", "target_value"), names(sensitivity_table))
      sensitivity_table[rounding_columns] <- lapply(sensitivity_table[rounding_columns], round, 2)
      sensitivity_names <- if (identical(current_language(), "en")) {
        c(scenario = "Analysis", dose = "Dose", interval = "Interval", auc24 = "AUC24", cmin = "Cmin", cmax = "Cmax", time_above_pct = "%T > threshold", target_value = "Target value", in_target = "In target", weights = "Weights")
      } else {
        c(scenario = "Analyse", dose = "Dose", interval = "Intervalle", auc24 = "AUC24", cmin = "Cmin", cmax = "Cmax", time_above_pct = "%T > seuil", target_value = "Valeur cible", in_target = "Dans cible", weights = "Poids")
      }
      if (identical(result$target_metric, "CminCmax")) sensitivity_names[["target_value"]] <- tx("Écart cible max.", "Maximum target distance")
      names(sensitivity_table) <- unname(sensitivity_names[names(sensitivity_table)])
      provenance <- result$provenance
      model_fingerprints <- lapply(provenance$models, function(item) {
        tags$li(paste0(item$id, " : ", item$sha256 %||% tx("empreinte indisponible", "fingerprint unavailable")))
      })

      document <- localize_ui(tags$html(
        tags$head(
          tags$meta(charset = "utf-8"),
          tags$title("Rapport TDM"),
          tags$style(HTML("body{font-family:Arial,sans-serif;color:#17202a;max-width:1050px;margin:32px auto;padding:0 24px;line-height:1.45}h1{border-bottom:3px solid #176b70;padding-bottom:12px}h2{margin-top:30px;border-bottom:1px solid #dfe4e8;padding-bottom:6px}.summary{display:grid;grid-template-columns:repeat(3,1fr);border:1px solid #dfe4e8}.summary div{padding:12px;border-right:1px solid #dfe4e8;border-bottom:1px solid #dfe4e8}.summary span{display:block;color:#66717c;font-size:11px;text-transform:uppercase}.summary strong{display:block;margin-top:4px}.summary small{display:block;margin-top:3px;color:#66717c}table{width:100%;border-collapse:collapse;font-size:12px}th,td{border:1px solid #dfe4e8;padding:7px;text-align:left}th{background:#f1f4f5}img{width:100%;height:auto}.warning{border-left:4px solid #b16916;background:#fffaf0;padding:12px 15px;margin-top:28px;color:#5e4b23}.meta{color:#66717c;font-size:12px}.provenance{font-family:Consolas,monospace;font-size:11px;overflow-wrap:anywhere}@media print{body{margin:0}.no-print{display:none}h2{break-after:avoid}table,img{break-inside:avoid}}"))
        ),
        tags$body(
          h1("Rapport d'analyse TDM"),
          p(class = "meta", tx(
            paste0("Généré le ", format(Sys.time(), "%d/%m/%Y à %H:%M"), " · voie ", result$route, " · aucune donnée conservée après la session"),
            paste0("Generated on ", format(Sys.time(), "%Y-%m-%d at %H:%M"), " · route ", result$route, " · no data retained after the session")
          )),
          div(
            class = "summary",
            div(span(tx("AUC24 MAP · état stationnaire", "MAP AUC24 · steady state")), strong(format_metric(exposure$steady_state_auc24)), tags$small(tx("Horizon de 24 h comparable au ML", "24-hour horizon comparable with ML"))),
            if (is.finite(result$ml_status$auc24 %||% NA_real_)) {
              div(span(tx("AUC24 ML expérimentale", "Experimental ML AUC24")), strong(format_metric(result$ml_status$auc24)), tags$small(tx("Prédiction directe sur 24 h", "Direct 24-hour prediction")))
            },
            div(span(tx(paste0("AUC MAP · fenêtre ", format_metric(exposure$historical_coverage_hours), " h"), paste0("MAP AUC · ", format_metric(exposure$historical_coverage_hours), "-hour window"))), strong(format_metric(exposure$historical_auc24)), tags$small(partial_auc_detail)),
            div(span("C0 actuelle (MAP)"), strong(format_metric(exposure$historical_c0))),
            div(span("C0 à l'état stationnaire"), strong(format_metric(exposure$steady_state_c0))),
            div(span("Recommandation"), strong(paste0(best$dose, " mg / ", best$interval, " h"))),
            div(span("PTA / cible"), strong(paste0(if (is.finite(best$p_target)) paste0(round(100 * best$p_target), " %") else tx("Non calculée", "Not calculated"), " · ", target_definition)))
          ),
          h2("Applicabilité et qualité"),
          if (length(result$quality$messages)) tags$ul(lapply(result$quality$messages, tags$li)) else p("Aucune anomalie générique détectée."),
          tags$ul(lapply(result$quality$populations, function(item) tags$li(tags$strong(item$model), paste0(" : ", item$population)))),
          fit_diagnostics_ui(result$fit_diagnostics),
          h2("Réglages appliqués"),
          p(tx(
            if (identical(result$settings$residual_error_mode, "fixed_cv")) {
              paste0("Erreur résiduelle MAP : CV proportionnel fixé à ", format_metric(result$settings$fixed_residual_cv), " %, avec représentation SIGMA adaptée à l'échelle d'erreur du modèle.")
            } else {
              "Erreur résiduelle MAP : valeurs du modèle."
            },
            if (identical(result$settings$residual_error_mode, "fixed_cv")) {
              paste0("MAP residual error: proportional CV fixed at ", format_metric(result$settings$fixed_residual_cv), "%, with SIGMA representation adapted to the model error scale.")
            } else {
              "MAP residual error: model values."
            }
          )),
          p(tx("Cible analysée : ", "Analyzed target: "), tags$strong(target_definition)),
          if (!is.null(result$target_preset)) p(tx("La cible a été initialisée depuis le préréglage documenté cité dans les références; les valeurs du rapport sont celles utilisées lors de l'analyse.", "The target was initialized from the documented preset cited in the references; report values are those used for the analysis.")),
          h2("Ajustement"),
          tags$img(src = report_plot_uri(fit_plot), alt = "Ajustement pharmacocinétique"),
          h2("Comparaison des doses supplémentaires"),
          p(tx(paste0(comparison$additional_doses, " doses simulées à partir de h ", format_metric(comparison$future_start), "."), paste0(comparison$additional_doses, " doses simulated from h ", format_metric(comparison$future_start), "."))),
          tags$img(src = report_plot_uri(comparison_plot), alt = "Comparaison de posologies"),
          if (!is.null(concordance) && !is.null(ml_comparison_plot)) {
            tagList(
              h2("Concordance MAP-BE et ML"),
              p(concordance$interpretation),
              tags$img(src = report_plot_uri(ml_comparison_plot, height = 4), alt = "Comparaison des AUC24 MAP-BE et ML")
            )
          },
          if (length(result$ml_status$domain_warnings %||% list())) {
            div(
              class = "warning",
              tags$strong("Avertissement : prédiction ML extrapolée"),
              tags$ul(lapply(
                result$ml_status$domain_warnings,
                function(item) tags$li(format_ml_domain_warning(item, current_language(), result$ml_status$unit %||% "mg.h/L"))
              )),
              p("La recommandation de dose reste calculée par MAP-BE.")
            )
          },
          if (isTRUE(ml_explanation$available) && !is.null(ml_explanation_plot)) {
            tagList(
              h2("Explication locale DALEX de l'AUC24 ML"),
              p("Les contributions décrivent le calcul de XGBoost par rapport à une référence synthétique; elles ne sont pas causales."),
              tags$img(src = report_plot_uri(ml_explanation_plot, height = 6), alt = "Décomposition locale DALEX de l'AUC24 ML")
            )
          },
          h2("Distribution prédictive"),
          p(distribution_summary_text),
          h2("Sensibilité du model averaging"),
          report_table(sensitivity_table),
          h2("Modèles et pondérations"),
          report_table(model_table),
          h2("Administrations"),
          report_table(dose_table),
          h2("Concentrations et covariables"),
          report_table(observation_table),
          h2("Références"),
          tags$ul(references),
          h2("Traçabilité"),
          p(class = "provenance", provenance$r),
          p(class = "provenance", paste0("mrgsolve ", provenance$packages$mrgsolve, " · mapbayr ", provenance$packages$mapbayr,
            " · Pmetrics ", provenance$packages$Pmetrics, " · shiny ", provenance$packages$shiny)),
          tags$strong("Empreintes SHA-256 des modèles"),
          tags$ul(class = "provenance", model_fingerprints),
          tags$strong("Hypothèses"),
          tags$ul(lapply(provenance$assumptions, tags$li)),
          p(class = "provenance", ml_status_message(result$ml_status, current_language())),
          div(
            class = "warning",
            tags$strong("Prototype de recherche et d'enseignement"),
            p("Ce rapport ne constitue pas une prescription ni une recommandation clinique validée. Vérifiez les données, unités, horaires, voie d'administration, population source et limites du modèle. Toute décision reste sous la responsabilité du professionnel de santé.")
          )
        )
      ), current_language())
      htmltools::save_html(document, file = file)
    }
  )

  pd_return <- reactiveVal(NULL)
  pd_open_tdm <- function(view) function(spec = NULL) {
    if (!is.null(spec)) {
      updateRadioButtons(session, "model_source", selected = if (spec$source == "code") "custom" else "library")
      if (spec$source == "code") updateTextAreaInput(session, "custom_code", value = spec$code)
      else updateSelectInput(session, "model_id", selected = spec$id)
      updateCheckboxInput(session, "enable_averaging", value = FALSE)
      session$onFlushed(function() updateSelectInput(session, "administration_route", selected = spec$route), once = TRUE)
    }
    pd_return(view)
    bslib::nav_select("main_navigation", "analysis", session = session)
    bslib::nav_select("analysis_tabs", "data", session = session)
    showNotification(tx("Renseignez doses, concentrations et covariables, puis lancez l'analyse mapbayr. Les donnees existantes ne sont pas remplacees.", "Enter doses, concentrations and covariates, then run the mapbayr analysis. Existing data are not replaced."), duration = 15)
  }
  output$pd_tdm_return <- renderUI({
    shiny::req(pd_return())
    actionButton("return_to_pd", tx("Retour a la pharmacodynamie", "Return to pharmacodynamics"), icon = icon("arrow-left"))
  })
  observeEvent(input$return_to_pd, {
    bslib::nav_select("main_navigation", "pd", session = session)
    bslib::nav_select("pd_workspace", pd_return(), session = session)
  })
  pd_server("pd", analysis_store, report_plot_uri, reactive(workbench_import()), session_model_dir, session_model_cache, pd_open_tdm("general"))
  onco_server("onco", report_plot_uri, reactive(workbench_import()), session_model_dir, session_model_cache, analysis_store, pd_open_tdm("oncology"))
  infection_server("infection", analysis_store, report_plot_uri, reactive(workbench_import()), session_model_dir, session_model_cache, pd_open_tdm("infection"))

  apply_workbench <- function(raw) {
    payload <- workbench_payload(raw)
    if (payload$view == "ddi") {
      config <- payload$config
      pending_ddi_target(config$targets)
      ddi_target_count(length(config$targets))
      pending_ddi_routes(stats::setNames(lapply(payload$models, function(model) list(model = model$id, route = model$route)), c("1", "2")))
      ddi_subject_1(NULL); ddi_subject_2(NULL); ddi_custom_1(NULL); ddi_custom_2(NULL)
      for (side in 1:2) {
        model <- payload$models[[side]]
        updateRadioButtons(session, paste0("ddi_source_", side), selected = model$source)
        if (model$source == "code") updateTextAreaInput(session, paste0("ddi_code_", side), value = model$code)
        else {
          updateSelectInput(session, paste0("ddi_model_", side), selected = model$id)
          updateSelectInput(session, paste0("ddi_route_", side), selected = model$route)
        }
        regimen <- if (side == 1) config$affected else config$driver
        for (name in c("dose", "interval", "infusion")) updateNumericInput(session, paste0("ddi_", name, "_", side), value = regimen[[name]])
      }
      updateSelectInput(session, "ddi_interaction_type", selected = config$type)
      for (name in c("factor", "strength", "c50", "hill", "kdeg", "kinact", "start_day", "stop_day")) updateNumericInput(session, paste0("ddi_", name), value = config[[name]])
      updateNumericInput(session, "ddi_followup", value = config$followup_days)
      ddi_store(NULL); ddi_assembly(NULL)
      bslib::nav_select("main_navigation", "ddi", session = session)
      session$sendCustomMessage("workbench-ack", list(id = payload$id, ok = TRUE))
    } else {
      bslib::nav_select("main_navigation", "pd", session = session)
      bslib::nav_select("pd_workspace", switch(payload$view, onco = "oncology", infection = "infection", "general"), session = session)
    }
    workbench_import(payload)
  }
  observeEvent(input$workbench_payload, tryCatch(apply_workbench(input$workbench_payload),
    error = function(e) session$sendCustomMessage("workbench-ack", list(id = input$workbench_payload$id, ok = FALSE, error = conditionMessage(e)))))
  for (file_id in c("ddi_workshop_file", "pd_workshop_file", "onco_workshop_file", "infection_workshop_file")) local({
    id <- file_id
    observeEvent(input[[id]], tryCatch({
      file <- input[[id]]
      on.exit(unlink(file$datapath), add = TRUE)
      if (file$size > 450000) stop("Workshop JSON limit: 450 kB.")
      payload <- jsonlite::fromJSON(file$datapath, simplifyVector = FALSE)
      if (!is.list(payload) || !identical(payload$type, "pk-workbench")) stop("Invalid workshop JSON.")
      payload$id <- paste0("file-", gsub("[^a-zA-Z0-9]", "", basename(tempfile())))
      apply_workbench(payload)
      showNotification(tx("Atelier importe. Le code personnel doit etre compile explicitement.", "Workshop imported. Custom code must be compiled explicitly."), type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = 10)))
  })
  session$onSessionEnded(function() { workbench_import(NULL); pending_ddi_target(NULL); pending_ddi_routes(list()) })

  ddi_choices_for <- function(side) {
    choices <- catalog_choices_i18n(lang = current_language())
    snapshot <- if (side == 1L) ddi_subject_1() else ddi_subject_2()
    if (!is.null(snapshot) && !snapshot$id %in% unname(choices)) {
      choices <- c(choices, stats::setNames(snapshot$id, paste0(snapshot$label, " · TDM/Lego")))
    }
    choices
  }

  ddi_route_ui <- function(side) {
    model_id <- input[[paste0("ddi_model_", side)]] %||% ""
    snapshot <- if (side == 1L) ddi_subject_1() else ddi_subject_2()
    routes <- if (!is.null(snapshot) && identical(snapshot$id, model_id) && !model_id %in% MODEL_CATALOG$id) {
      snapshot$route
    } else if (model_id %in% MODEL_CATALOG$id) {
      model_routes(model_record(model_id))
    } else character()
    shiny::req(length(routes))
    pending <- pending_ddi_routes()[[as.character(side)]]
    selected <- if (!is.null(pending) && identical(pending$model, model_id)) pending$route else isolate(input[[paste0("ddi_route_", side)]]) %||% routes[[1]]
    if (!selected %in% routes) selected <- routes[[1]]
    labels <- vapply(routes, function(route) {
      if (identical(route, "Oral")) tx("Orale", "Oral") else if (identical(route, "IV")) "IV" else route
    }, character(1))
    selectInput(paste0("ddi_route_", side), tx("Voie d'administration", "Administration route"), choices = stats::setNames(routes, labels), selected = selected)
  }

  output$ddi_route_1_ui <- renderUI(ddi_route_ui(1L))
  output$ddi_route_2_ui <- renderUI(ddi_route_ui(2L))
  observeEvent(list(input$ddi_route_1, input$ddi_route_2, input$ddi_model_1, input$ddi_model_2), {
    pending <- pending_ddi_routes()
    for (side in names(pending)) {
      value <- pending[[side]]
      if (identical(value$model, input[[paste0("ddi_model_", side)]]) && identical(value$route, input[[paste0("ddi_route_", side)]])) pending[[side]] <- NULL
    }
    pending_ddi_routes(pending)
  })

  ddi_context_for <- function(side) {
    if (identical(input[[paste0("ddi_source_", side)]], "code")) {
      context <- if (side == 1L) ddi_custom_1() else ddi_custom_2()
      if (is.null(context) || !identical(context$code, input[[paste0("ddi_code_", side)]])) stop(tx("Compilez d'abord le code actuel.", "Compile the current code first."))
      context$adm_cmt <- as.integer(input[[paste0("ddi_adm_", side)]] %||% 1L)
      if (!context$adm_cmt %in% seq_along(context$model@cmtL)) stop("Invalid administration compartment.")
      return(context)
    }
    model_id <- input[[paste0("ddi_model_", side)]] %||% ""
    route <- input[[paste0("ddi_route_", side)]] %||% NULL
    snapshot <- if (side == 1L) ddi_subject_1() else ddi_subject_2()
    context <- if (!is.null(snapshot) && identical(snapshot$id, model_id)) snapshot else ddi_library_context(model_id, route)
    ddi_context_with_route(context, route)
  }

  ddi_target_parameters <- reactive({
    model_id <- input$ddi_model_1 %||% ""
    snapshot <- ddi_subject_1()
    tryCatch({
      if (identical(input$ddi_source_1, "code")) {
        context <- ddi_custom_1()
        if (!is.null(context) && identical(context$code, input$ddi_code_1)) ddi_parameter_table(context) else data.frame()
      } else if (!is.null(snapshot) && identical(snapshot$id, model_id)) ddi_parameter_table(snapshot)
      else if (model_id %in% MODEL_CATALOG$id) ddi_code_parameter_table(read_library_code(model_id))
      else data.frame()
    }, error = function(error) data.frame())
  })

  ddi_target_input_id <- function(index, field = "parameter") {
    stem <- if (field == "fraction") "ddi_target_fraction" else "ddi_target_parameter"
    if (index == 1L) stem else paste0(stem, "_", index)
  }

  ddi_targets_from_inputs <- function() lapply(seq_len(ddi_target_count()), function(index) list(
    parameter = input[[ddi_target_input_id(index)]] %||% "",
    fraction = numeric_input_value(ddi_target_input_id(index, "fraction"), 1)
  ))

  output$ddi_target_ui <- renderUI({
    parameters <- ddi_target_parameters()
    if (!nrow(parameters)) return(div(class = "ddi-source-note", tx("Aucun paramètre structurel positif n'est disponible.", "No positive structural parameter is available.")))
    choices <- stats::setNames(parameters$name, paste0(parameters$name, " (", signif(parameters$value, 5), ")"))
    count <- min(ddi_target_count(), nrow(parameters), 6L)
    pending <- pending_ddi_target()
    used <- character()
    rows <- lapply(seq_len(count), function(index) {
      imported <- if (length(pending) >= index) pending[[index]] else NULL
      selected <- imported$parameter %||% isolate(input[[ddi_target_input_id(index)]]) %||% ""
      if (!selected %in% parameters$name || selected %in% used) {
        clearance <- setdiff(grep("(^CL$|CL_|_CL$|TVCL)", parameters$name, value = TRUE, ignore.case = TRUE), used)
        selected <- if (length(clearance)) clearance[[1]] else setdiff(parameters$name, used)[[1]]
      }
      used <<- c(used, selected)
      fraction <- imported$fraction %||% isolate(input[[ddi_target_input_id(index, "fraction")]]) %||% 1
      div(class = "ddi-target-row",
        selectInput(ddi_target_input_id(index), paste0(tx("Paramètre cible ", "Target parameter "), index), choices = choices, selected = selected),
        numericInput(ddi_target_input_id(index, "fraction"), tx("Fraction affectée (analogue fm)", "Affected fraction (fm analogue)"), fraction, min = 0, max = 1, step = 0.05))
    })
    if (!is.null(pending)) session$onFlushed(function() pending_ddi_target(NULL), once = TRUE)
    tagList(
      rows,
      div(class = "ddi-target-actions",
        actionButton("ddi_add_target", tx("Ajouter une cible", "Add target"), icon = icon("plus"), disabled = if (count >= min(6L, nrow(parameters))) "disabled" else NULL),
        actionButton("ddi_remove_target", tx("Retirer la dernière", "Remove last"), icon = icon("minus"), disabled = if (count <= 1L) "disabled" else NULL)),
      tags$small(class = "form-text", tx("Pour chaque cible : P/P0 = (1-f) + f×M(C). Les covariables, ETA et erreurs résiduelles sont exclues; les cibles doivent être distinctes.", "For each target: P/P0 = (1-f) + f×M(C). Covariates, ETAs, and residual errors are excluded; targets must be distinct."))
    )
  })

  observeEvent(input$ddi_add_target, {
    maximum <- min(6L, nrow(ddi_target_parameters()))
    ddi_target_count(min(ddi_target_count() + 1L, maximum))
  })
  observeEvent(input$ddi_remove_target, ddi_target_count(max(1L, ddi_target_count() - 1L)))
  observe({
    maximum <- max(1L, min(6L, nrow(ddi_target_parameters())))
    if (ddi_target_count() > maximum) ddi_target_count(maximum)
  })

  observeEvent(input$ddi_model_1, {
    snapshot <- ddi_subject_1()
    if (!is.null(snapshot) && !identical(snapshot$id, input$ddi_model_1)) ddi_subject_1(NULL)
    ddi_store(NULL)
  }, ignoreInit = TRUE)
  observeEvent(input$ddi_model_2, {
    snapshot <- ddi_subject_2()
    if (!is.null(snapshot) && !identical(snapshot$id, input$ddi_model_2)) ddi_subject_2(NULL)
    ddi_store(NULL)
  }, ignoreInit = TRUE)

  ddi_snapshot_tdm <- function(side) {
    result <- analysis_store()
    if (is.null(result)) {
      return(showNotification(tx("Lancez d'abord un ajustement dans l'onglet Analyse.", "First run a fit in the Analysis tab."), type = "warning", duration = 6))
    }
    fits <- successful_fits(result$fits)
    if (!length(fits)) return(showNotification(tx("Aucun ajustement exploitable n'est disponible.", "No usable fit is available."), type = "warning", duration = 6))
    selected <- input[[paste0("ddi_model_", side)]] %||% ""
    indexes <- which(vapply(fits, function(fit) identical(fit$id, selected), logical(1)))
    if (!length(indexes) && length(fits) == 1L) indexes <- 1L
    if (!length(indexes)) {
      return(showNotification(tx("Sélectionnez dans ce bloc un modèle présent dans le dernier ajustement TDM.", "In this block, select a model included in the latest TDM fit."), type = "warning", duration = 7))
    }
    fit <- fits[[indexes[[1]]]]
    context <- ddi_fit_context(fit)
    updateRadioButtons(session, paste0("ddi_source_", side), selected = "library")
    if (side == 1L) ddi_subject_1(context) else ddi_subject_2(context)
    updateSelectInput(session, paste0("ddi_model_", side), choices = ddi_choices_for(side), selected = context$id)

    doses <- fit$source_doses %||% data.frame()
    if (nrow(doses)) {
      dose <- doses[order(doses$time), , drop = FALSE][nrow(doses), , drop = FALSE]
      if (is.finite(as.numeric(dose$amount))) updateNumericInput(session, paste0("ddi_dose_", side), value = as.numeric(dose$amount))
      if (is.finite(as.numeric(dose$interval)) && as.numeric(dose$interval) > 0) updateNumericInput(session, paste0("ddi_interval_", side), value = as.numeric(dose$interval))
      if (is.finite(as.numeric(dose$infusion))) updateNumericInput(session, paste0("ddi_infusion_", side), value = as.numeric(dose$infusion))
    }
    ddi_store(NULL)
    showNotification(tx(paste0("Ajustement TDM repris pour le modèle ", side, "."), paste0("TDM fit loaded for model ", side, ".")), type = "message", duration = 4)
  }

  observeEvent(input$ddi_use_tdm_1, ddi_snapshot_tdm(1L))
  observeEvent(input$ddi_use_tdm_2, ddi_snapshot_tdm(2L))

  ddi_use_population <- function(side) {
    if (identical(input[[paste0("ddi_source_", side)]], "code")) return(ddi_load_custom(side))
    model_id <- input[[paste0("ddi_model_", side)]] %||% ""
    if (side == 1L) ddi_subject_1(NULL) else ddi_subject_2(NULL)
    choices <- catalog_choices_i18n(lang = current_language())
    if (!model_id %in% unname(choices)) model_id <- if (side == 1L) DEFAULT_MODEL else unname(choices)[[min(2, length(choices))]]
    updateSelectInput(session, paste0("ddi_model_", side), choices = choices, selected = model_id)
    ddi_store(NULL)
  }
  observeEvent(input$ddi_use_population_1, ddi_use_population(1L))
  observeEvent(input$ddi_use_population_2, ddi_use_population(2L))

  ddi_source_ui <- function(side) {
    if (identical(input[[paste0("ddi_source_", side)]], "code")) return(NULL)
    snapshot <- if (side == 1L) ddi_subject_1() else ddi_subject_2()
    model_id <- input[[paste0("ddi_model_", side)]] %||% ""
    if (!is.null(snapshot) && identical(snapshot$id, model_id)) {
      div(class = "ddi-source-note personalized", tx("Paramètres du dernier ajustement TDM ou modèle Lego, conservés uniquement dans cette session.", "Parameters from the latest TDM fit or Lego model, retained only in this session."))
    } else {
      div(class = "ddi-source-note", tx("Paramètres populationnels du modèle sélectionné.", "Population parameters from the selected model."))
    }
  }
  output$ddi_source_note_1 <- renderUI(ddi_source_ui(1L))
  output$ddi_source_note_2 <- renderUI(ddi_source_ui(2L))

  ddi_load_custom <- function(side) tryCatch({
    context <- withProgress(message = tx("Compilation mrgsolve", "Compiling mrgsolve"), value = 0.2,
      ddi_custom_context(input[[paste0("ddi_code_", side)]], side, session_model_dir, session_model_cache, ALLOW_CUSTOM_MODELS))
    if (side == 1L) ddi_custom_1(context) else ddi_custom_2(context)
    ddi_store(NULL); ddi_assembly(NULL)
  }, error = function(error) showNotification(conditionMessage(error), type = "error", duration = 10))
  observeEvent(input$ddi_load_code_1, ddi_load_custom(1L))
  observeEvent(input$ddi_load_code_2, ddi_load_custom(2L))
  for (side in 1:2) local({
    side <- side
    context <- reactive(if (side == 1L) ddi_custom_1() else ddi_custom_2())
    output[[paste0("ddi_custom_status_", side)]] <- renderUI({
      value <- context()
      p(class = "ddi-source-note", if (is.null(value) || !identical(value$code, input[[paste0("ddi_code_", side)]]))
        tx("Code non compile ou modifie.", "Code not compiled or modified.") else tx("Modele compile pour cette session.", "Model compiled for this session."))
    })
    output[[paste0("ddi_custom_adm_", side)]] <- renderUI({
      value <- context(); shiny::req(value)
      selectInput(paste0("ddi_adm_", side), tx("Compartiment d'administration", "Administration compartment"),
        stats::setNames(seq_along(value$model@cmtL), value$model@cmtL), value$adm_cmt)
    })
  })

  observeEvent(input$ddi_interaction_type, {
    inhibition <- input$ddi_interaction_type %in% c("inhibition", "hill_inhibition")
    updateNumericInput(session, "ddi_strength", max = if (inhibition) 1 else 20,
      value = if (inhibition) min(input$ddi_strength %||% 1, 1) else input$ddi_strength %||% 1)
  })

  ddi_config_from_inputs <- function() {
    targets <- ddi_targets_from_inputs()
    list(
      affected = list(dose = numeric_input_value("ddi_dose_1"), interval = numeric_input_value("ddi_interval_1"), infusion = numeric_input_value("ddi_infusion_1", 0)),
      driver = list(dose = numeric_input_value("ddi_dose_2"), interval = numeric_input_value("ddi_interval_2"), infusion = numeric_input_value("ddi_infusion_2", 0)),
      target = targets[[1]]$parameter,
      target_fraction = targets[[1]]$fraction,
      targets = targets,
      type = input$ddi_interaction_type %||% "factor",
      factor = numeric_input_value("ddi_factor", 0.5),
      strength = numeric_input_value("ddi_strength", 1),
      c50 = numeric_input_value("ddi_c50", 1),
      hill = numeric_input_value("ddi_hill", 2),
      kdeg = numeric_input_value("ddi_kdeg", 0.02),
      kinact = numeric_input_value("ddi_kinact", 0.1),
      start_day = numeric_input_value("ddi_start_day", 3),
      stop_day = numeric_input_value("ddi_stop_day", 10),
      followup_days = numeric_input_value("ddi_followup", 3)
    )
  }

  output$ddi_status <- renderUI({
    if (is.null(ddi_store())) span(class = "status-pill idle", tx("Non simulé", "Not simulated"))
    else span(class = "status-pill ok", tx("Simulation prête", "Simulation ready"))
  })

  observeEvent(list(ddi_config_from_inputs(), input$ddi_model_1, input$ddi_model_2,
    input$ddi_source_1, input$ddi_source_2, input$ddi_code_1, input$ddi_code_2, input$ddi_adm_1, input$ddi_adm_2,
    input$ddi_route_1, input$ddi_route_2, ddi_subject_1(), ddi_subject_2()), {
    ddi_store(NULL); ddi_assembly(NULL)
  }, ignoreInit = TRUE)
  output$ddi_chain <- renderUI({
    label <- function(side) {
      if (identical(input[[paste0("ddi_source_", side)]], "code")) return(paste("mrgsolve", side))
      id <- input[[paste0("ddi_model_", side)]] %||% ""
      if (id %in% MODEL_CATALOG$id) model_record(id)$label[[1]] else paste("TDM / Lego", side)
    }
    div(class = "model-chain",
      div(class = "model-block exposure", icon("pills"), h3(label(1)), tags$small(tx("Molecule affectee", "Affected drug"))),
      span(class = "model-arrow", icon("arrow-right")),
      div(class = "model-block delay", icon("link"), h3(paste(vapply(ddi_targets_from_inputs(), function(target) paste0(target$parameter, " (f=", format_metric(target$fraction, 2), ")"), character(1)), collapse = " + ")),
        tags$code(ddi_mechanism_equation(input$ddi_interaction_type %||% "factor"))),
      span(class = "model-arrow", icon("arrow-left")),
      div(class = "model-block response", icon("pills"), h3(label(2)), tags$small(tx("Molecule interagissante", "Interacting drug"))))
  })
  observeEvent(input$assemble_ddi, tryCatch({
    ddi_assembly(NULL)
    contexts <- withProgress(message = tx("Assemblage", "Assembly"), value = 0.2, list(ddi_context_for(1L), ddi_context_for(2L)))
    config <- ddi_config_from_inputs()
    code <- ddi_model_code(contexts[[1]], config)
    # Compile the generated model before offering a C++ download.
    compile_model(custom_code = code, allow_custom = ALLOW_CUSTOM_MODELS, custom_soloc = session_model_dir, custom_cache = session_model_cache)
    ddi_assembly(list(code = code, script = ddi_export_script(config, contexts[[1]], contexts[[2]])))
    bslib::nav_select("ddi_tabs", "model", session = session)
  }, error = function(error) showNotification(conditionMessage(error), type = "error", duration = 10)))
  output$ddi_code_ready <- reactive(!is.null(ddi_assembly()))
  outputOptions(output, "ddi_code_ready", suspendWhenHidden = FALSE)
  output$ddi_code_status <- renderUI(p(if (is.null(ddi_assembly()))
    tx("Aucun modele assemble pour les valeurs actuelles.", "No model assembled for the current values.") else
    tx("C++ compile. Entrees externes : DDI_CP, DDI_ACTIVE ou DDI_ACTIVITY selon le mecanisme. Le script R calcule le couplage dynamique et simule les deux modeles.", "C++ compiled. External inputs: DDI_CP, DDI_ACTIVE or DDI_ACTIVITY depending on mechanism. The R script computes dynamic coupling and simulates both models.")))
  output$ddi_generated_code <- renderText({ shiny::req(ddi_assembly()); ddi_assembly()$code })
  output$download_ddi_code <- downloadHandler(filename = function() "pk-interaction.cpp", content = function(file) { shiny::req(ddi_assembly()); writeLines(ddi_assembly()$code, file) })
  output$download_ddi_r <- downloadHandler(filename = function() "pk-interaction.R", content = function(file) { shiny::req(ddi_assembly()); writeLines(ddi_assembly()$script, file) })

  observeEvent(input$run_ddi, {
    ddi_store(NULL)
    tryCatch({
      shiny::validate(shiny::need(isTRUE(input$ddi_accept_disclaimer), tx("Lisez et acceptez l'avertissement avant de lancer la simulation.", "Read and accept the warning before running the simulation.")))
      result <- withProgress(message = tx("Simulation de l'interaction", "Simulating interaction"), value = 0.1, {
        context_1 <- ddi_context_for(1L)
        incProgress(0.3, detail = tx("Préparation du modèle 1", "Preparing model 1"))
        context_2 <- ddi_context_for(2L)
        incProgress(0.3, detail = tx("Préparation du modèle 2", "Preparing model 2"))
        ddi_simulate(ddi_config_from_inputs(), context_1, context_2, delta = 0.1)
      })
      ddi_store(result)
      bslib::nav_select("ddi_tabs", "simulation", session = session)
      showNotification(tx("Simulation DDI terminée.", "DDI simulation completed."), type = "message", duration = 3)
    }, error = function(error) {
      showNotification(paste0(tx("Simulation impossible : ", "Simulation failed: "), conditionMessage(error)), type = "error", duration = 9)
    })
  })

  output$ddi_ready <- reactive(!is.null(ddi_store()))
  shiny::outputOptions(output, "ddi_ready", suspendWhenHidden = FALSE)
  output$ddi_empty <- renderUI({
    if (!is.null(ddi_store())) return(NULL)
    div(class = "empty-results ddi-empty", tags$strong(tx("Assemblez puis simulez l'interaction.", "Assemble and simulate the interaction.")), span(tx("Le modèle 2 pilote une ou plusieurs cibles structurelles du modèle 1, avec une fraction affectée propre à chaque cible.", "Model 2 drives one or more structural targets in model 1, each with its own affected fraction.")))
  })

  ddi_targets_label <- function(result) paste(vapply(result$config$targets, function(target) paste0(target$parameter, " (f=", format_metric(target$fraction, 2), ")"), character(1)), collapse = " + ")

  ddi_concentration_figure <- function(result) {
    data <- result$affected_profile
    labels <- c(baseline = tx("Sans interaction", "Without interaction"), interaction = tx("Avec interaction", "With interaction"))
    data$scenario <- factor(data$scenario, levels = names(labels), labels = unname(labels))
    ggplot(data, aes(day, concentration, color = scenario, linetype = scenario)) +
      annotate("rect", xmin = result$config$start_day, xmax = result$config$stop_day, ymin = -Inf, ymax = Inf, fill = "#f4e4dc", alpha = 0.6) +
      geom_vline(xintercept = c(result$config$start_day, result$config$stop_day), color = "#8b756d", linetype = "dotted") +
      geom_line(linewidth = 1) +
      scale_color_manual(values = stats::setNames(c("#59636c", "#176b70"), unname(labels))) +
      scale_linetype_manual(values = stats::setNames(c("dashed", "solid"), unname(labels))) +
      labs(x = tx("Jour", "Day"), y = tx("Concentration (unité du modèle)", "Concentration (model unit)"), color = NULL, linetype = NULL, caption = paste0(result$affected$label, " · ", ddi_targets_label(result))) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), legend.position = "top", plot.caption = element_text(hjust = 0, color = "#66717c"))
  }

  ddi_effect_figure <- function(result) {
    effect <- result$effect
    maximum <- max(effect$driver_concentration, na.rm = TRUE)
    driver <- unique(effect[, c("day", "driver_concentration"), drop = FALSE])
    driver_relative <- if (is.finite(maximum) && maximum > 0) 100 * driver$driver_concentration / maximum else 0
    target_metrics <- paste0(effect$target, " (f=", vapply(effect$fraction, format_metric, character(1), digits = 2), ")")
    data <- rbind(
      data.frame(day = effect$day, value = 100 * effect$modifier, metric = target_metrics),
      data.frame(day = driver$day, value = driver_relative, metric = tx("Exposition du modèle 2 (normalisée)", "Model 2 exposure (normalized)"))
    )
    metrics <- unique(data$metric)
    colors <- stats::setNames(c(grDevices::hcl.colors(length(metrics) - 1L, "Dark 3"), "#a4441f"), metrics)
    ggplot(data, aes(day, value, color = metric, linetype = metric)) +
      annotate("rect", xmin = result$config$start_day, xmax = result$config$stop_day, ymin = -Inf, ymax = Inf, fill = "#f4e4dc", alpha = 0.6) +
      geom_hline(yintercept = 100, color = "#8c969f", linewidth = 0.5) +
      geom_vline(xintercept = c(result$config$start_day, result$config$stop_day), color = "#8b756d", linetype = "dotted") +
      geom_line(linewidth = 1) +
      scale_color_manual(values = colors) +
      labs(x = tx("Jour", "Day"), y = tx("Valeur relative (%)", "Relative value (%)"), color = NULL, linetype = NULL) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), legend.position = "top")
  }

  output$ddi_metrics <- renderUI({
    result <- ddi_store()
    shiny::req(result)
    metrics <- result$metrics
    target_cards <- lapply(result$config$targets, function(target) div(
      span(tx("Paramètre initial", "Initial parameter")),
      strong(format_metric(metrics$baseline_parameters[[target$parameter]], 3)),
      tags$small(paste0(target$parameter, " · f=", format_metric(target$fraction, 2)))
    ))
    div(
      class = "exposure-strip ddi-metrics",
      div(span(tx("AUC24 sans interaction", "AUC24 without interaction")), strong(format_metric(metrics$baseline$auc)), tags$small(tx("unité concentration·h", "concentration unit·h"))),
      div(span(tx("AUC24 avec interaction", "AUC24 with interaction")), strong(format_metric(metrics$interaction$auc)), tags$small(paste0("×", format_metric(metrics$auc_ratio, 2)))),
      div(span(tx("Cmin avec / sans", "Cmin with / without")), strong(paste0("×", format_metric(metrics$cmin_ratio, 2))), tags$small(paste(format_metric(metrics$baseline$cmin), "→", format_metric(metrics$interaction$cmin)))),
      div(span(tx("Cmax avec / sans", "Cmax with / without")), strong(paste0("×", format_metric(metrics$cmax_ratio, 2))), tags$small(paste(format_metric(metrics$baseline$cmax), "→", format_metric(metrics$interaction$cmax)))),
      target_cards,
      div(span(tx("Facteur min · max", "Minimum · maximum factor")), strong(paste(format_metric(metrics$minimum_modifier, 2), "·", format_metric(metrics$maximum_modifier, 2))), tags$small(tx("fractions affectées incluses", "including affected fractions")))
    )
  })

  output$ddi_concentration_plot <- renderPlot({ shiny::req(ddi_store()); ddi_concentration_figure(ddi_store()) })
  output$ddi_effect_plot <- renderPlot({ shiny::req(ddi_store()); ddi_effect_figure(ddi_store()) })

  ddi_schedule_table <- function(result) {
    table <- result$schedule
    data.frame(
      `Modèle` = table$model,
      Source = ifelse(table$source == "tdm", "TDM", tx("Population", "Population")),
      Voie = table$route,
      `Début (jour)` = table$start_day,
      `Fin (jour)` = table$stop_day,
      Dose = table$dose,
      `Intervalle (h)` = table$interval,
      `Perfusion (h)` = table$infusion,
      check.names = FALSE
    )
  }

  output$ddi_schedule <- renderDT({
    result <- ddi_store()
    shiny::req(result)
    table <- ddi_schedule_table(result)
    if (identical(current_language(), "en")) names(table) <- c("Model", "Source", "Route", "Start (day)", "Stop (day)", "Dose", "Interval (h)", "Infusion (h)")
    datatable(table, rownames = FALSE, options = list(dom = "t", pageLength = 2, scrollX = TRUE))
  })

  ddi_interaction_description <- function(result) {
    targets <- ddi_targets_label(result)
    if (result$config$type %in% c("reversible", "hill_inhibition", "tdi", "turnover_induction")) {
      fields <- c("c50", "hill", "kdeg", "kinact", "strength")
      return(paste(targets, "|", ddi_mechanism_equation(result$config$type), "|", paste(fields, unlist(result$config[fields]), sep = "=", collapse = "; ")))
    }
    switch(result$config$type,
      factor = tx(paste0(targets, " · M=", format_metric(result$config$factor, 2), " pendant le traitement du modèle 2."), paste0(targets, " · M=", format_metric(result$config$factor, 2), " during model 2 treatment.")),
      inhibition = paste0(targets, " · M(C2) = 1 - ", format_metric(result$config$strength, 2), "·C2/(", format_metric(result$config$c50, 3), "+C2)"),
      induction = paste0(targets, " · M(C2) = 1 + ", format_metric(result$config$strength, 2), "·C2/(", format_metric(result$config$c50, 3), "+C2)")
    )
  }

  ddi_reference_tag <- function(context) {
    tagList(
      tags$strong(context$label),
      p(context$citation),
      if (nzchar(context$doi)) tags$a(paste0("DOI ", context$doi), href = paste0("https://doi.org/", context$doi), target = "_blank", rel = "noopener noreferrer")
    )
  }

  output$ddi_method <- renderUI({
    result <- ddi_store()
    shiny::req(result)
    tagList(
      div(class = "ddi-method-grid", div(ddi_reference_tag(result$affected)), div(ddi_reference_tag(result$driver)), div(tags$strong(tx("Relation simulée", "Simulated relationship")), p(ddi_interaction_description(result)))),
      div(class = "analysis-diagnostics", tags$strong(tx("Limites d'interprétation", "Interpretation limits")), p(tx("Le lien entre les deux modèles est une relation utilisateur appliquée aux cibles choisies. La fraction affectée est un analogue structurel de fm renseigné par l'utilisateur, pas une fraction métabolisée estimée par le modèle. Cette construction ne démontre ni un mécanisme biologique ni une interaction clinique validée. Pour les relations Emax, la valeur IC50/EC50 doit employer exactement l'unité de concentration capturée par le modèle 2.", "The link between models is a user-defined relationship applied to the selected targets. The affected fraction is a user-supplied structural analogue of fm, not a fraction metabolized estimated by the model. This construction demonstrates neither a biological mechanism nor a validated clinical interaction. For Emax relationships, IC50/EC50 must use exactly the concentration unit captured by model 2.")), p(tx("Les valeurs populationnelles conservent les covariables par défaut du fichier de modèle. L'option TDM reprend un instantané des paramètres individuels de chaque molécule; elle ne fusionne pas les dossiers ni les données brutes.", "Population values retain the model file's default covariates. The TDM option uses a snapshot of each drug's individual parameters; it does not merge records or raw data."))),
      div(class = "privacy-notice ddi-privacy", tags$strong(tx("Session uniquement", "Session only")), span(tx("Aucun modèle importé, paramètre individuel ou résultat de simulation n'est conservé par ce module.", "No imported model, individual parameter, or simulation result is retained by this module.")))
    )
  })

  output$download_ddi_report <- downloadHandler(
    filename = function() paste0(if (identical(current_language(), "en")) "ddi-report-" else "rapport-ddi-", Sys.Date(), ".html"),
    content = function(file) {
      result <- ddi_store()
      shiny::req(result)
      metrics <- result$metrics
      schedule <- ddi_schedule_table(result)
      if (identical(current_language(), "en")) names(schedule) <- c("Model", "Source", "Route", "Start (day)", "Stop (day)", "Dose", "Interval (h)", "Infusion (h)")
      document <- tags$html(
        tags$head(tags$meta(charset = "utf-8"), tags$title(tx("Rapport de simulation DDI", "DDI simulation report")), tags$style(HTML("body{font-family:Arial,sans-serif;color:#17202a;max-width:1050px;margin:32px auto;padding:0 20px;line-height:1.5}img{max-width:100%;height:auto}table{width:100%;border-collapse:collapse}th,td{padding:8px;border:1px solid #dfe4e8;text-align:left}.warning{border-left:4px solid #b16916;background:#fffaf0;padding:12px;margin-top:24px}"))),
        tags$body(
          h1(tx("Simulation d'interaction médicamenteuse", "Drug-interaction simulation")),
          p(paste(result$affected$label, "+", result$driver$label)),
          p(ddi_interaction_description(result)),
          h2(tx("Exposition comparée", "Exposure comparison")),
          p(paste0("AUC24: ", format_metric(metrics$baseline$auc), " → ", format_metric(metrics$interaction$auc), " (×", format_metric(metrics$auc_ratio, 2), ")")),
          tags$img(src = report_plot_uri(ddi_concentration_figure(result), height = 5), alt = tx("Concentrations comparées", "Compared concentrations")),
          h2(tx("Relation d'interaction", "Interaction relationship")),
          tags$img(src = report_plot_uri(ddi_effect_figure(result), height = 4), alt = tx("Paramètre et exposition relative", "Relative parameter and exposure")),
          h2(tx("Scénario", "Scenario")),
          report_table(schedule),
          h2(tx("Sources des modèles", "Model sources")),
          tags$ul(tags$li(result$affected$citation, if (nzchar(result$affected$doi)) paste0(" DOI ", result$affected$doi)), tags$li(result$driver$citation, if (nzchar(result$driver$doi)) paste0(" DOI ", result$driver$doi))),
          div(class = "warning", tags$b(tx("Prototype de recherche et d'enseignement", "Research and teaching prototype")), p(tx("La relation d'interaction est une hypothèse définie par l'utilisateur. Ce rapport n'est ni une prescription ni un dispositif médical.", "The interaction relationship is a user-defined assumption. This report is neither a prescription nor a medical device.")))
        )
      )
      htmltools::save_html(document, file = file)
    }
  )

  output$download_scenarios <- downloadHandler(
    filename = function() paste0("tdm-scenarios-", Sys.Date(), ".csv"),
    content = function(file) {
      shiny::req(analysis_store())
      utils::write.csv(analysis_store()$recommendations, file, row.names = FALSE, na = "")
    }
  )

  output$catalog_table <- renderDT({
    english <- identical(current_language(), "en")
    table <- MODEL_CATALOG[, c(if (english) "drugEn" else "drug", "model", if (english) "administrationCategoriesEn" else "administrationCategories", "citation", if (english) "populationEn" else "population", "doi", if (english) "modelTypeEn" else "modelType", if (english) "implementationStatusLabelEn" else "implementationStatusLabel", if (english) "analysisEligibilityLabelEn" else "analysisEligibilityLabel")]
    administration_column <- if (english) "administrationCategoriesEn" else "administrationCategories"
    table[[administration_column]] <- vapply(table[[administration_column]], function(values) paste(values, collapse = " + "), character(1))
    has_doi <- !is.na(table$doi) & nzchar(table$doi)
    table$doi <- ifelse(
      has_doi,
      paste0('<a href="https://doi.org/', table$doi, '" target="_blank" rel="noopener noreferrer">', table$doi, "</a>"),
      tx("Source à confirmer", "Source to be confirmed")
    )
    datatable(
      table,
      rownames = FALSE,
      filter = "top",
      escape = FALSE,
      options = list(pageLength = 20, scrollX = TRUE),
      colnames = app_t(current_language(), c("Molécule", "Modèle", "Administration", "Article source", "Population de l'article", "DOI", "Type", "Implémentation", "Usage TDM"))
    )
  })
}

shinyApp(ui = app_ui, server = server)
