ddi_mechanism_choices <- function(t) stats::setNames(
  c("factor", "inhibition", "induction", "reversible", "hill_inhibition", "tdi", "turnover_induction"),
  c(t("Facteur constant", "Constant factor"), t("Inhibition Emax", "Emax inhibition"), t("Stimulation Emax instantanee", "Instantaneous Emax stimulation"),
    t("Inhibition reversible (Ki)", "Reversible inhibition (Ki)"), t("Inhibition sigmoide (Hill)", "Sigmoid inhibition (Hill)"),
    t("Inhibition dependante du temps", "Time-dependent inhibition"), t("Induction avec renouvellement", "Induction with turnover")))

ddi_mechanism_equation <- function(type) switch(type,
  factor = "P = P0 * factor (active treatment)", inhibition = "P/P0 = max(0.01, 1-Imax*C/(IC50+C))",
  induction = "P/P0 = 1+Emax*C/(EC50+C)", reversible = "P/P0 = max(0.01, 1/(1+C/Ki))",
  hill_inhibition = "P/P0 = max(0.01, 1-Imax*C^h/(IC50^h+C^h))",
  tdi = "dA/dt = kdeg*(1-A)-kinact*C/(KI+C)*A; P/P0=max(0.01,A)",
  turnover_induction = "dA/dt = kdeg*(1+Emax*C/(EC50+C)-A); P/P0=A")

ddi_mechanism_ui <- function(t) tagList(
  selectInput("ddi_interaction_type", t("Mecanisme", "Mechanism"), ddi_mechanism_choices(t)),
  conditionalPanel("input.ddi_interaction_type == 'factor'", numericInput("ddi_factor", t("Facteur multiplicatif", "Multiplicative factor"), 0.5, min = 0.01, max = 20)),
  conditionalPanel("input.ddi_interaction_type != 'factor'",
    numericInput("ddi_c50", "IC50 / EC50 / Ki / KI", 1, min = 0.000001),
    tags$small(t("Meme unite que la concentration du modele 2; pas de conversion libre/totale automatique.", "Same unit as model 2 concentration; no automatic unbound/total conversion."))),
  conditionalPanel("['inhibition','induction','hill_inhibition','turnover_induction'].includes(input.ddi_interaction_type)", numericInput("ddi_strength", "Imax / Emax", 1, min = 0, max = 20)),
  conditionalPanel("input.ddi_interaction_type == 'hill_inhibition'", numericInput("ddi_hill", "Hill", 2, min = 0.1, max = 10)),
  conditionalPanel("['tdi','turnover_induction'].includes(input.ddi_interaction_type)", numericInput("ddi_kdeg", "kdeg (1/h)", 0.02, min = 0.000001, max = 10, step = 0.01)),
  conditionalPanel("input.ddi_interaction_type == 'tdi'", numericInput("ddi_kinact", "kinact (1/h)", 0.1, min = 0, max = 10, step = 0.01)))

ddi_paste_ui <- function(side, t) {
  id <- function(name) paste0("ddi_", name, "_", side)
  tagList(
    radioButtons(id("source"), NULL, stats::setNames(c("library", "code"), c(t("Bibliotheque / TDM", "Library / TDM"), "mrgsolve / C++")), inline = TRUE),
    conditionalPanel(paste0("input.", id("source"), " == 'code'"),
      textAreaInput(id("code"), t("Code mrgsolve (bibliotheque, Lego ou personnel)", "mrgsolve code (library, Lego or custom)"), rows = 10, width = "100%"),
      actionButton(id("load_code"), t("Compiler le modele", "Compile model"), icon = icon("gears")),
      uiOutput(id("custom_status")), uiOutput(id("custom_adm"))))
}

ddi_custom_context <- function(code, side, soloc, cache, allow_custom = FALSE) {
  if (!is.character(code) || length(code) != 1 || !nzchar(trimws(code)) || nchar(code, type = "bytes") > 200000) stop("Supply a complete mrgsolve model (maximum 200 kB).")
  model <- compile_model(custom_code = code, allow_custom = allow_custom, custom_soloc = soloc, custom_cache = cache)
  model <- mrgsolve::zero_re(model)
  if (!length(model@cmtL)) stop("The PK model must have at least one compartment.")
  if (!any(model_capture_names(model) %in% c("DV", "CP", "CONC", "CONC_PLASMA", "CONCENTRATION"))) stop("Capture a concentration named DV, CP or CONC in the PK model.")
  adm <- tagged_compartment(model, "ADM")
  list(id = paste0("custom_", side), label = paste("Custom PK", side), drug = "Custom", source = "custom",
    code = code, model = model, route = "Custom", adm_cmt = if (is.finite(adm)) adm else 1L, split_lego = is_split_lego_model(model),
    covariate_names = parse_covariates(code)$name, citation = "User-supplied model", doi = "")
}

ddi_model_code <- function(context, config) {
  config <- ddi_validate_mechanism(config)
  config <- ddi_validate_targets(config, ddi_parameter_table(context)$name)
  targets <- vapply(config$targets, `[[`, character(1), "parameter")
  fractions <- vapply(config$targets, `[[`, numeric(1), "fraction")
  ddi_numeric(config$factor, "Factor", 0.01, 20)
  ddi_numeric(config$strength, "Imax / Emax", 0, if (config$type == "inhibition") 1 else 20)
  ddi_numeric(config$c50, "IC50 / EC50", 1e-9, 1e9)
  code <- context$code
  if (is.null(code) || !nzchar(code)) stop("Original model code unavailable.")
  if (grepl("\\bDDI_", code, perl = TRUE)) stop("Reserved DDI names already exist in this model.")
  # Edit only a parameter declaration; never replace arbitrary C++ expressions.
  lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  lines <- lines[!grepl("PK_LEGO_SPEC", lines, fixed = TRUE)]
  block <- ""
  replaced <- stats::setNames(integer(length(targets)), targets)
  for (i in seq_along(lines)) {
    header <- regmatches(lines[[i]], regexec("^\\s*(?:\\$([A-Za-z]+)|\\[([A-Za-z]+)\\])", lines[[i]], perl = TRUE))[[1]]
    if (length(header)) block <- toupper(paste0(header[[2]], header[[3]]))
    if (block != "PARAM") next
    for (target in targets) {
      pattern <- paste0("(?<![A-Za-z0-9_])", target, "(?=\\s*[:=])")
      if (grepl(pattern, sub("//.*", "", lines[[i]]), perl = TRUE)) {
        lines[[i]] <- sub(pattern, paste0("DDI_BASE_", target), lines[[i]], perl = TRUE)
        replaced[[target]] <- replaced[[target]] + 1L
      }
    }
  }
  if (any(replaced != 1L)) stop("Each target must have one explicit numeric $PARAM declaration.")
  multiplier <- switch(config$type,
    factor = "(DDI_ACTIVE > 0.5 ? DDI_FACTOR : 1.0)",
    inhibition = "fmax(0.01, 1.0-DDI_STRENGTH*fmax(0.0,DDI_CP)/(DDI_C50+fmax(0.0,DDI_CP)))",
    induction = "(1.0+DDI_STRENGTH*fmax(0.0,DDI_CP)/(DDI_C50+fmax(0.0,DDI_CP)))",
    reversible = "fmax(0.01,1.0/(1.0+fmax(0.0,DDI_CP)/DDI_C50))",
    hill_inhibition = "fmax(0.01,1.0-DDI_STRENGTH/(1.0+pow(DDI_C50/fmax(1e-300,DDI_CP),DDI_HILL)))",
    tdi = "fmax(0.01,DDI_ACTIVITY)",
    turnover_induction = "fmax(0.01,DDI_ACTIVITY)",
    stop("Unknown interaction type."))
  declarations <- vapply(seq_along(targets), function(i) paste0(
    "double ", targets[[i]], " = DDI_BASE_", targets[[i]],
    " * ((1.0-DDI_FRACTION_", i, ")+DDI_FRACTION_", i, "*", multiplier, ");"
  ), character(1))
  declaration <- paste(declarations, collapse = "\n")
  main <- grep("^\\s*(\\$(MAIN|PK)\\b|\\[(MAIN|PK)\\])", lines, perl = TRUE, ignore.case = TRUE)
  if (length(main)) {
    i <- main[[1]]
    header <- regmatches(lines[[i]], regexpr("^\\s*(\\$(MAIN|PK)\\b|\\[(MAIN|PK)\\])", lines[[i]], perl = TRUE, ignore.case = TRUE))
    rest <- substring(lines[[i]], nchar(header) + 1L)
    lines[[i]] <- paste(header, declaration, rest, sep = "\n")
  } else lines <- c(lines, "$MAIN", declaration)
  fraction_parameters <- paste0("DDI_FRACTION_", seq_along(fractions), "=", fractions)
  paste(c("// Research only. Supply model 2 concentrations in DDI_CP at each time step.",
    "// For a constant factor use DDI_ACTIVE=1 only during interaction treatment.",
    "// Each target uses P/P0 = (1-fraction) + fraction*DDI_modifier.",
    if (config$type %in% c("tdi", "turnover_induction")) c(
      "// Dynamic mechanism: DDI_ACTIVITY is an EXTERNAL input, initialized at 1.",
      "// Use the coupled R export to integrate activity; DDI_CP alone is insufficient.",
      if (config$type == "tdi") "// dA/dt = kdeg*(1-A) - kinact*C/(KI+C)*A" else "// dA/dt = kdeg*(1+Emax*C/(EC50+C)-A)",
      paste0("// kdeg=", config$kdeg, " /h; kinact=", config$kinact, " /h")),
    lines, "$PARAM", paste(c(paste0("DDI_CP=0, DDI_ACTIVE=0, DDI_ACTIVITY=1, DDI_HILL=", config$hill,
      ", DDI_FACTOR=", config$factor, ", DDI_STRENGTH=", config$strength, ", DDI_C50=", config$c50), fraction_parameters), collapse = ", ")), collapse = "\n")
}

ddi_export_script <- function(config, affected, driver) {
  dump <- function(x) paste(capture.output(dput(x)), collapse = "\n")
  context_code <- function(context, name) {
    metadata <- context[setdiff(names(context), c("model", "code"))]
    paste0(name, " <- ", dump(metadata), "\n", name, "$model <- mrgsolve::mcode('", name, "', ", dump(context$code),
      ", soloc=workdir, quiet=TRUE) |> mrgsolve::zero_re()\n", name, "$model <- mrgsolve::param(", name, "$model, ", dump(as.list(mrgsolve::param(context$model))), ")")
  }
  helpers <- c("model_param_names", "pick_concentration_column", "trap_auc", "ddi_numeric", "ddi_filter_parameters", "ddi_parameter_table",
    "ddi_validate_regimen", "ddi_validate_targets", "ddi_validate_mechanism", "ddi_dose_data", "ddi_simulate_profile", "ddi_modifier", "ddi_profile_window", "ddi_simulate")
  paste(c("# Coupled PK simulation, research only. Requires mrgsolve and a C++ compiler.",
    "# Contains the selected models and parameter snapshots, not raw TDM observations.",
    "# Model 2 drives one or more parameters of model 1 on a 0.1 h grid; no reciprocal interaction.",
    "library(mrgsolve)", "`%||%` <- function(x,y) if (is.null(x) || !length(x)) y else x",
    "LEGO_STEADY_STATE_WARMUP_DOSES <- 50L",
    vapply(helpers, function(name) paste0(name, " <- ", paste(deparse(get(name, mode = "function")), collapse = "\n")), character(1)),
    "workdir <- tempfile('ddi-'); dir.create(workdir)",
    context_code(affected, "affected"), context_code(driver, "driver"), paste0("config <- ", dump(config)),
    "result <- ddi_simulate(config, affected, driver)",
    "with(subset(result$affected_profile, scenario=='baseline'), plot(day, concentration, type='l', lty=2, xlab='Day', ylab='Concentration'))",
    "with(subset(result$affected_profile, scenario=='interaction'), lines(day, concentration, col='#196f76'))",
    "print(result$metrics)",
    "# After use, unload compiled models before deleting workdir (especially on Windows)."), collapse = "\n\n")
}
