workbench_pk_model <- function(model) {
  if (!is.list(model)) stop("Invalid PK model specification.")
  if (identical(model$source, "code")) {
    if (!is.character(model$code) || length(model$code) != 1 || !nzchar(trimws(model$code)) || nchar(model$code, type = "bytes") > 200000) stop("Invalid mrgsolve code (maximum 200 kB).")
    if (!is.null(model$route) && !model$route %in% c("IV", "Oral")) stop("Custom PK route must be IV or Oral.")
  } else if (!identical(model$source, "library") || !model$id %in% MODEL_CATALOG$id || !model$route %in% model_routes(model_record(model$id))) stop("Unknown library model or route.")
  if (!is.null(model$time_unit) && !model$time_unit %in% c("h", "day")) stop("PK time unit must be h or day.")
  if (!is.null(model$concentration_scale)) ddi_numeric(model$concentration_scale, "Concentration scale", 1e-9, 1e9)
  model
}

workbench_payload <- function(payload) {
  if (!is.list(payload) || !identical(payload$version, 1L) && !identical(payload$version, 1) ||
      !is.character(payload$id) || length(payload$id) != 1 || !grepl("^[a-zA-Z0-9-]{1,80}$", payload$id) ||
      !is.character(payload$view) || length(payload$view) != 1 || !payload$view %in% c("ddi", "pd", "onco", "infection")) stop("Unsupported workshop payload.")
  if (nchar(jsonlite::toJSON(payload, auto_unbox = TRUE), type = "bytes") > 450000) stop("Workshop payload too large.")
  if (!is.list(payload$config)) stop("Missing workshop configuration.")
  if (!is.null(payload$models)) payload$models <- lapply(payload$models, workbench_pk_model)
  if (payload$view == "onco") {
    history <- payload$config$history
    if (is.list(history) && !is.data.frame(history)) {
      payload$config$history <- if (!length(history)) data.frame(time = numeric(), amount = numeric(), infusion = numeric()) else
        do.call(rbind, lapply(history, function(row) as.data.frame(row, stringsAsFactors = FALSE)))
    }
    payload$config <- onco_config(payload$config)
  } else if (payload$view == "infection") {
    payload$config <- infection_config(payload$config)
    expected <- if (payload$config$source == "pk") 1L else 0L
    if (length(payload$models) != expected) stop("Infectiology requires one PK model, or the current engine TDM session.")
    if (expected == 1L && payload$models[[1]]$route == "Oral" && (payload$config$infusion > 0 || payload$config$infusion2 > 0)) stop("Oral dosing requires infusion = 0.")
  } else if (payload$view == "pd") {
    exposure <- payload$config$exposure
    if (!exposure %in% c("exponential", "pk")) stop("Import concentration tables inside the engine.")
    if (exposure == "pk") {
      payload$config$regimen <- ddi_validate_regimen(payload$config$regimen, "PK")
      payload$config$exposure <- "exponential"
    }
    payload$config <- pd_config(payload$config)
    payload$config$exposure <- exposure
  } else {
    payload$config <- ddi_validate_mechanism(payload$config)
    if (!is.character(payload$config$target) || length(payload$config$target) != 1 || !grepl("^[A-Za-z_][A-Za-z0-9_]*$", payload$config$target)) stop("Invalid interaction target.")
    payload$config$affected <- ddi_validate_regimen(payload$config$affected, "Model 1")
    payload$config$driver <- ddi_validate_regimen(payload$config$driver, "Model 2")
    payload$config$start_day <- ddi_numeric(payload$config$start_day, "Start", 0.25, 180)
    payload$config$stop_day <- ddi_numeric(payload$config$stop_day, "Stop", payload$config$start_day + 0.25, 365)
    payload$config$followup_days <- ddi_numeric(payload$config$followup_days, "Follow-up", 0, 180)
    if (!is.list(payload$models) || length(payload$models) != 2) stop("Two PK models are required.")
  }
  if ((payload$view == "onco" && payload$config$free_pk || payload$view == "pd" && payload$config$exposure == "pk") && length(payload$models) != 1) stop("One PK model is required.")
  payload
}
