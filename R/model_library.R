if (!exists("APP_ROOT", inherits = TRUE)) {
  APP_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

MODEL_ROOT <- file.path(APP_ROOT, "models")
MODEL_CACHE_DIR <- file.path(tempdir(), "pk-mipd-model-cache")
dir.create(MODEL_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

catalog_path <- file.path(MODEL_ROOT, "catalog.json")
if (!file.exists(catalog_path)) {
  stop("Missing model catalog. Run `npm run tdm:index` from the repository root.")
}

MODEL_CATALOG <- jsonlite::fromJSON(catalog_path, simplifyDataFrame = TRUE)$models
MODEL_CATALOG$label <- paste(MODEL_CATALOG$drug, MODEL_CATALOG$model, sep = " - ")

.model_cache <- new.env(parent = emptyenv())

short_hash <- function(text) {
  hash <- 2166136261 %% 2147483647
  for (byte in as.integer(charToRaw(enc2utf8(text)))) {
    hash <- (hash * 31 + byte) %% 2147483647
  }
  sprintf("%08x", as.integer(hash))
}

safe_model_id <- function(value) {
  value <- gsub("[^A-Za-z0-9_]", "_", value)
  substr(value, 1, 48)
}

model_routes <- function(record) {
  values <- record$routes
  if (is.list(values)) values <- values[[1]]
  unique(as.character(values[!is.na(values) & nzchar(values)]))
}

model_supports_route <- function(record, route) {
  route %in% model_routes(record)
}

model_administration_modes <- function(record, route = NULL) {
  values <- record$administrationModes
  if (is.list(values)) values <- values[[1]]
  values <- unique(as.character(values[!is.na(values) & nzchar(values)]))
  if (is.null(route)) return(values)
  values[
    if (identical(route, "Oral")) values == "ORAL"
    else if (identical(route, "IM")) values == "IM"
    else startsWith(values, "IV_")
  ]
}

model_administration_mode <- function(record, route, selected = NULL) {
  modes <- model_administration_modes(record, route)
  selected <- as.character(selected %||% "")
  if (length(selected) && selected[[1]] %in% modes) return(selected[[1]])
  if (length(modes) != 1L) stop("The model must declare exactly one administration mode for route ", route, ".")
  modes[[1]]
}

model_supports_administration_mode <- function(record, route, mode) {
  model_supports_route(record, route) && mode %in% model_administration_modes(record, route)
}

administration_mode_label <- function(mode, lang = "fr") {
  labels <- if (identical(lang, "en")) {
    c(ORAL = "Oral", IM = "Intramuscular", IV_INTERMITTENT = "Intermittent IV", IV_CONTINUOUS = "Continuous IV")
  } else {
    c(ORAL = "Orale", IM = "Intramusculaire", IV_INTERMITTENT = "IV intermittente", IV_CONTINUOUS = "IV continue")
  }
  unname(labels[[mode]] %||% mode)
}

model_administration_cmt <- function(record, route) {
  column <- if (identical(route, "Oral")) "oralCmt" else if (identical(route, "IM")) "imCmt" else "ivCmt"
  value <- record[[column]]
  if (is.list(value)) value <- value[[1]]
  value <- as.character(value[[1]] %||% "")
  if (is.na(value) || !nzchar(value)) stop("No ", route, " administration compartment is declared for this model.")
  value
}

model_analysis_eligible <- function(record) {
  value <- record$analysisEligible
  if (is.list(value)) value <- value[[1]]
  isTRUE(as.logical(value[[1]] %||% FALSE))
}

catalog_choices <- function(drug = NULL, route = NULL, mode = NULL, analysis_only = FALSE) {
  rows <- MODEL_CATALOG
  if (isTRUE(analysis_only)) {
    keep <- vapply(seq_len(nrow(rows)), function(index) model_analysis_eligible(rows[index, , drop = FALSE]), logical(1))
    rows <- rows[keep, , drop = FALSE]
  }
  if (!is.null(drug)) rows <- rows[rows$drug == drug, , drop = FALSE]
  if (!is.null(route)) {
    keep <- vapply(seq_len(nrow(rows)), function(index) model_supports_route(rows[index, , drop = FALSE], route), logical(1))
    rows <- rows[keep, , drop = FALSE]
  }
  if (!is.null(mode)) {
    keep <- vapply(seq_len(nrow(rows)), function(index) model_supports_administration_mode(rows[index, , drop = FALSE], route, mode), logical(1))
    rows <- rows[keep, , drop = FALSE]
  }
  stats::setNames(rows$id, rows$label)
}

model_record <- function(model_id) {
  row <- MODEL_CATALOG[MODEL_CATALOG$id == model_id, , drop = FALSE]
  if (nrow(row) != 1) stop("Unknown library model: ", model_id)
  row
}

read_library_code <- function(model_id) {
  row <- model_record(model_id)
  path <- file.path(MODEL_ROOT, row$file[[1]])
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

parse_covariates <- function(code) {
  lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  starts <- which(grepl("^\\s*(?:\\$PARAM|\\[PARAM\\]).*@covariates", lines, ignore.case = TRUE, perl = TRUE))
  if (!length(starts)) {
    return(data.frame(name = character(), value = numeric(), description = character(), scope = character()))
  }

  parsed <- list()
  index <- 1L
  for (start in starts) {
    following <- if (start < length(lines)) seq.int(start + 1L, length(lines)) else integer()
    stop_at <- following[grepl("^\\s*(?:\\$[A-Za-z]|\\[[A-Za-z]+\\])", lines[following], perl = TRUE)][1]
    if (is.na(stop_at)) stop_at <- length(lines) + 1L
    block <- lines[seq.int(start + 1L, stop_at - 1L)]

    for (line in block) {
      clean <- trimws(sub("//.*$", "", line))
      if (!nzchar(clean)) next
      parts <- trimws(strsplit(clean, ":", fixed = TRUE)[[1]])
      if (length(parts) < 2 || !grepl("^[A-Za-z_][A-Za-z0-9_]*$", parts[[1]])) next
      value <- suppressWarnings(as.numeric(parts[[2]]))
      if (!is.finite(value)) next
      description <- if (length(parts) > 2) paste(parts[-c(1, 2)], collapse = ": ") else ""
      scope <- if (grepl("\\[administration\\]", description, ignore.case = TRUE)) "administration" else "patient"
      description <- trimws(gsub("\\s*\\[administration\\]\\s*", "", description, ignore.case = TRUE))
      parsed[[index]] <- data.frame(
        name = parts[[1]],
        value = value,
        description = description,
        scope = scope,
        stringsAsFactors = FALSE
      )
      index <- index + 1L
    }
  }

  if (!length(parsed)) {
    return(data.frame(name = character(), value = numeric(), description = character(), scope = character()))
  }
  unique(do.call(rbind, parsed))
}

LEGO_SPEC_PREFIX <- "// PK_LEGO_SPEC_V1:"
LEGO_NODE_KINDS <- c("depot", "transit", "central", "periph", "metab", "effect", "response", "tumor", "interaction")
LEGO_VOLUME_KINDS <- c("central", "periph", "metab")
LEGO_PD_KINDS <- c("effect", "response", "tumor", "interaction")
LEGO_RESERVED_COVARIATES <- c(
  "ID", "TIME", "CMT", "AMT", "EVID", "RATE", "II", "ADDL", "SS", "DV", "MDV", "IPRED", "ETA", "EPS"
)

lego_record_list <- function(value, field) {
  if (is.null(value)) return(list())
  if (is.data.frame(value)) {
    return(lapply(seq_len(nrow(value)), function(index) as.list(value[index, , drop = FALSE])))
  }
  if (!is.list(value)) stop("Lego field `", field, "` must be an array.")
  value
}

lego_scalar <- function(value, field, default = NULL) {
  if (is.null(value) || !length(value)) {
    if (!is.null(default)) return(default)
    stop("Missing Lego field `", field, "`.")
  }
  if (is.list(value)) value <- unlist(value, recursive = TRUE, use.names = FALSE)
  if (length(value) != 1L) stop("Lego field `", field, "` must contain one value.")
  value[[1]]
}

lego_number <- function(value, field, minimum, maximum, default = NULL, integer = FALSE) {
  raw <- lego_scalar(value, field, default)
  number <- suppressWarnings(as.numeric(raw))
  if (length(number) != 1L || !is.finite(number) || number < minimum || number > maximum) {
    stop("Invalid Lego numeric field `", field, "`.")
  }
  if (isTRUE(integer) && abs(number - round(number)) > .Machine$double.eps^0.5) {
    stop("Lego field `", field, "` must be an integer.")
  }
  if (isTRUE(integer)) as.integer(round(number)) else number
}

lego_text <- function(value, field, pattern, maximum = 32L) {
  text <- as.character(lego_scalar(value, field))
  if (nchar(text, type = "bytes") > maximum || !grepl(pattern, text, perl = TRUE)) {
    stop("Invalid Lego text field `", field, "`.")
  }
  text
}


lego_advanced_fields <- function(node) {
  if (identical(node$kind, "tumor")) return(c("t0", "kg", if (node$growth != "exponential") "cap", "kill", "ec50", "res"))
  if (!identical(node$kind, "interaction")) return(character())
  switch(node$mechanism, factor = "factor", reversible = "c50",
    inhibition = c("strength", "c50"), hill_inhibition = c("strength", "c50", "hill"),
    induction = c("strength", "c50"), tdi = c("c50", "kinact", "kdeg"),
    turnover_induction = c("strength", "c50", "kdeg"))
}

# Only validated graph fields reach these templates, never user-supplied C++.
lego_advanced_expression <- function(node, cp, state, p) {
  c <- paste0("fmax(0.0,", cp, ")")
  if (node$kind == "tumor") {
    size <- paste0("fmax(0.0,", state, ")")
    growth <- switch(node$growth,
      exponential = "", logistic = paste0("*(1-", size, "/", p("cap"), ")"),
      gompertz = paste0("*log(", p("cap"), "/fmax(1e-12,", size, "))"))
    return(list(initial = p("t0"), factor = "1", derivative = paste0(
      p("kg"), "*", size, growth, "-", p("kill"), "*", c, "/(", p("ec50"), "+", c,
      ")*exp(-", p("res"), "*SOLVERTIME)*", size)))
  }
  raised <- function(x) paste0("pow(", x, ",", p("hill"), ")")
  factor <- switch(node$mechanism,
    factor = p("factor"),
    reversible = paste0("1/(1+", c, "/", p("c50"), ")"),
    inhibition = paste0("1-", p("strength"), "*", c, "/(", p("c50"), "+", c, ")"),
    hill_inhibition = paste0("1-", p("strength"), "*", raised(c), "/(", raised(p("c50")), "+", raised(c), ")"),
    induction = paste0("1+", p("strength"), "*", c, "/(", p("c50"), "+", c, ")"),
    tdi = state, turnover_induction = state)
  derivative <- switch(node$mechanism,
    tdi = paste0(p("kdeg"), "*(1-", state, ")-", p("kinact"), "*", c, "/(", p("c50"), "+", c, ")*", state),
    turnover_induction = paste0(p("kdeg"), "*(1+", p("strength"), "*", c, "/(", p("c50"), "+", c, ")-", state, ")"),
    "0")
  if (!node$mechanism %in% c("factor", "induction")) factor <- paste0("fmax(0.01,", factor, ")")
  list(initial = "1", derivative = derivative, factor = paste0("(", factor, ")"))
}

normalize_lego_spec <- function(specification) {
  if (!is.list(specification)) stop("The Lego specification must be a JSON object.")
  version <- lego_number(specification$version, "version", 1, 4, integer = TRUE)
  input_nodes <- lego_record_list(specification$nodes, "nodes")
  input_edges <- lego_record_list(specification$edges, "edges")
  input_covariates <- lego_record_list(specification$covariates, "covariates")

  if (!length(input_nodes) || length(input_nodes) > 20L) {
    stop("A Lego model must contain between 1 and 20 compartments.")
  }
  if (length(input_edges) > 60L) stop("A Lego model cannot contain more than 60 transfers.")
  if (length(input_covariates) > 50L) stop("A Lego model cannot contain more than 50 covariates.")

  nodes <- lapply(seq_along(input_nodes), function(index) {
    input <- input_nodes[[index]]
    if (!is.list(input)) stop("Each Lego compartment must be an object.")
    id <- lego_number(input$id, paste0("nodes[", index, "].id"), 1, 1000000, integer = TRUE)
    kind <- lego_text(input$kind, paste0("nodes[", index, "].kind"), "^[a-z]+$", 16L)
    if (!kind %in% LEGO_NODE_KINDS) stop("Unsupported Lego compartment kind: ", kind)
    name <- lego_text(input$name, paste0("nodes[", index, "].name"), "^[A-Za-z_][A-Za-z0-9_]*$", 32L)
    node <- list(
      id = id,
      kind = kind,
      name = name,
      dose = lego_number(input$dose, paste0("nodes[", index, "].dose"), 0, 1e9, default = 0)
    )
    if (!kind %in% LEGO_PD_KINDS) {
      node$inputType <- lego_text(input$inputType %||% "bolus", paste0("nodes[", index, "].inputType"), "^(bolus|zero_order)$", 16L)
      node$inputDuration <- lego_number(input$inputDuration, paste0("nodes[", index, "].inputDuration"), 1e-6, 1e6, default = 1)
      node$tlag <- lego_number(input$tlag, paste0("nodes[", index, "].tlag"), 0, 1e6, default = 0)
      node$doseFraction <- lego_number(input$doseFraction, paste0("nodes[", index, "].doseFraction"), 1e-6, 100, default = 100)
      node$inputDurationTlagOf <- if (is.null(input$inputDurationTlagOf)) NULL else lego_number(input$inputDurationTlagOf, paste0("nodes[", index, "].inputDurationTlagOf"), 1, 1000000, integer = TRUE)
      node$fractionComplementOf <- if (is.null(input$fractionComplementOf)) NULL else lego_number(input$fractionComplementOf, paste0("nodes[", index, "].fractionComplementOf"), 1, 1000000, integer = TRUE)
    }
    if (kind %in% LEGO_VOLUME_KINDS) {
      node$vol <- lego_number(input$vol, paste0("nodes[", index, "].vol"), 1e-6, 1e9)
    }
    if (identical(kind, "effect")) {
      node$ke0 <- lego_number(input$ke0, paste0("nodes[", index, "].ke0"), 0, 1e6, default = 0.4)
      node$source <- lego_number(input$source, paste0("nodes[", index, "].source"), 1, 1000000, integer = TRUE)
    }
    if (identical(kind, "response")) {
      node$kin <- lego_number(input$kin, paste0("nodes[", index, "].kin"), 0, 1e12, default = 10)
      node$kout <- lego_number(input$kout, paste0("nodes[", index, "].kout"), 1e-12, 1e6, default = 0.15)
      node$smax <- lego_number(input$smax, paste0("nodes[", index, "].smax"), -1e6, 1e6, default = 3)
      node$sc50 <- lego_number(input$sc50, paste0("nodes[", index, "].sc50"), 1e-12, 1e12, default = 3)
      node$source <- lego_number(input$source, paste0("nodes[", index, "].source"), 1, 1000000, integer = TRUE)
    }
    if (kind %in% c("tumor", "interaction")) {
      node$source <- lego_number(input$source, "source", 1, 1000000, integer = TRUE)
      if (kind == "tumor") node$growth <- lego_text(input$growth, "growth", "^(exponential|logistic|gompertz)$", 16L)
      if (kind == "interaction") {
        node$mechanism <- lego_text(input$mechanism, "mechanism", "^(factor|reversible|inhibition|hill_inhibition|induction|tdi|turnover_induction)$", 24L)
        node$targetFrom <- lego_number(input$targetFrom, "targetFrom", 1, 1000000, integer = TRUE)
        node$targetTo <- if (identical(input$targetTo, "OUT")) "OUT" else lego_number(input$targetTo, "targetTo", 1, 1000000, integer = TRUE)
      }
      for (key in lego_advanced_fields(node)) {
        lower <- if (key %in% c("kg", "kill", "res", "strength", "kinact")) 0 else 1e-6
        upper <- if (key == "strength" && grepl("inhibition", node$mechanism %||% "")) 1 else 1e6
        node[[key]] <- lego_number(input[[key]], key, lower, upper)
      }
    }
    if (kind %in% LEGO_PD_KINDS && node$dose != 0) stop("A PD or interaction block cannot receive drug mass.")
    node
  })

  node_ids <- vapply(nodes, `[[`, integer(1), "id")
  node_names <- vapply(nodes, `[[`, character(1), "name")
  if (anyDuplicated(node_ids)) stop("Lego compartment identifiers must be unique.")
  if (anyDuplicated(node_names)) stop("Lego compartment names must be unique.")

  for (node in nodes) {
    if (node$kind %in% LEGO_PD_KINDS && !node$source %in% node_ids) {
      stop("The source of Lego compartment `", node$name, "` does not exist.")
    }
  }
  for (node in Filter(function(n) n$kind %in% LEGO_PD_KINDS, nodes)) {
    source <- nodes[[match(node$source, node_ids)]]
    if (source$id == node$id || !source$kind %in% c(LEGO_VOLUME_KINDS, "effect")) stop("A block source must be a concentration.")
    visited <- node$id
    while (identical(source$kind, "effect")) {
      if (source$id %in% visited) stop("Cyclic Lego effect sources are not supported.")
      visited <- c(visited, source$id)
      source <- nodes[[match(source$source, node_ids)]]
    }
  }
  dosed_nodes <- Filter(function(node) !node$kind %in% LEGO_PD_KINDS && node$dose > 0, nodes)
  dosed_ids <- vapply(dosed_nodes, `[[`, integer(1), "id")
  for (node in dosed_nodes) {
    if (!is.null(node$fractionComplementOf)) {
      if (identical(node$id, node$fractionComplementOf) || !node$fractionComplementOf %in% dosed_ids) {
        stop("A complementary Lego fraction must reference another dosed compartment.")
      }
      source <- nodes[[match(node$fractionComplementOf, node_ids)]]
      if (!is.null(source$fractionComplementOf)) stop("Chained complementary Lego fractions are not supported.")
    }
    if (!is.null(node$inputDurationTlagOf)) {
      if (identical(node$id, node$inputDurationTlagOf) || !node$inputDurationTlagOf %in% dosed_ids) {
        stop("A linked Lego duration must reference another dosed compartment.")
      }
      source <- nodes[[match(node$inputDurationTlagOf, node_ids)]]
      if (source$tlag <= 0) stop("A linked Lego duration requires a strictly positive source Tlag.")
    }
  }
  if (length(dosed_nodes) > 1L) {
    fraction_value <- function(node) {
      if (is.null(node$fractionComplementOf)) return(node$doseFraction)
      100 - nodes[[match(node$fractionComplementOf, node_ids)]]$doseFraction
    }
    fractions <- vapply(dosed_nodes, fraction_value, numeric(1))
    if (any(fractions <= 0)) stop("Parallel Lego dose fractions must be strictly positive.")
    fraction_total <- sum(fractions)
    if (abs(fraction_total - 100) > 1e-3) stop("Parallel Lego dose fractions must total 100%.")
  }

  edges <- lapply(seq_along(input_edges), function(index) {
    input <- input_edges[[index]]
    if (!is.list(input)) stop("Each Lego transfer must be an object.")
    from <- lego_number(input$from, paste0("edges[", index, "].from"), 1, 1000000, integer = TRUE)
    if (!from %in% node_ids) stop("A Lego transfer references an unknown source compartment.")
    raw_to <- lego_scalar(input$to, paste0("edges[", index, "].to"))
    to <- if (identical(as.character(raw_to), "OUT")) {
      "OUT"
    } else {
      lego_number(raw_to, paste0("edges[", index, "].to"), 1, 1000000, integer = TRUE)
    }
    if (!identical(to, "OUT") && !to %in% node_ids) {
      stop("A Lego transfer references an unknown destination compartment.")
    }
    if (!identical(to, "OUT") && identical(from, to)) stop("A Lego transfer cannot loop to itself.")
    kinetics <- lego_text(input$kinetics %||% "first_order", paste0("edges[", index, "].kinetics"), "^(first_order|michaelis_menten|hill)$", 24L)
    parameterization <- lego_text(input$eliminationParameterization %||% "rate", paste0("edges[", index, "].eliminationParameterization"), "^(rate|clearance)$", 16L)
    transfer_parameterization <- lego_text(input$transferParameterization %||% "rate", paste0("edges[", index, "].transferParameterization"), "^(rate|clearance)$", 16L)
    source_node <- nodes[[match(from, node_ids)]]
    if (source_node$kind %in% LEGO_PD_KINDS || (!identical(to, "OUT") && nodes[[match(to, node_ids)]]$kind %in% LEGO_PD_KINDS)) stop("Mass transfers cannot enter or leave PD/interaction blocks.")
    if (identical(parameterization, "clearance") && (!identical(to, "OUT") || !source_node$kind %in% LEGO_VOLUME_KINDS)) {
      stop("Clearance parameterization requires elimination from a volume compartment.")
    }
    if (identical(transfer_parameterization, "clearance")) {
      target_node <- if (identical(to, "OUT")) NULL else nodes[[match(to, node_ids)]]
      if (is.null(target_node) || !source_node$kind %in% LEGO_VOLUME_KINDS || !target_node$kind %in% LEGO_VOLUME_KINDS) {
        stop("Intercompartmental clearance requires a transfer between two volume compartments.")
      }
    }
    list(
      from = from,
      to = to,
      kinetics = kinetics,
      k = lego_number(input[["k"]], paste0("edges[", index, "].k"), 0, 1e6, default = 0.2),
      vmax = lego_number(input$vmax, paste0("edges[", index, "].vmax"), 1e-12, 1e12, default = 10),
      km = lego_number(input$km, paste0("edges[", index, "].km"), 1e-12, 1e12, default = 10),
      gamma = lego_number(input$gamma, paste0("edges[", index, "].gamma"), 1e-6, 100, default = 1),
      eliminationParameterization = parameterization,
      cl = lego_number(input$cl, paste0("edges[", index, "].cl"), 1e-12, 1e12, default = 5),
      transferParameterization = transfer_parameterization,
      q = lego_number(input$q, paste0("edges[", index, "].q"), 1e-12, 1e12, default = 5)
    )
  })

  edge_keys <- vapply(edges, function(edge) paste(edge$from, edge$to, sep = "->"), character(1))
  if (anyDuplicated(edge_keys)) stop("Duplicate Lego transfers are not supported.")
  q_edges <- Filter(function(edge) identical(edge$transferParameterization, "clearance"), edges)
  if (length(q_edges)) {
    q_keys <- vapply(q_edges, function(edge) paste(sort(c(edge$from, edge$to)), collapse = "::"), character(1))
    for (key in unique(q_keys)) {
      values <- vapply(q_edges[q_keys == key], `[[`, numeric(1), "q")
      if (length(values) > 1L && max(values) - min(values) > 1e-10) stop("Reciprocal Lego transfers must share the same intercompartmental clearance Q.")
    }
  }
  for (node in Filter(function(n) n$kind == "interaction", nodes)) {
    if (!paste(node$targetFrom, node$targetTo, sep = "->") %in% edge_keys) stop("The interaction target flux does not exist.")
  }
  if (!any(vapply(nodes, function(node) node$kind %in% LEGO_VOLUME_KINDS, logical(1)))) {
    stop("A Lego TDM model requires a central, peripheral or metabolite compartment.")
  }

  covariates <- lapply(seq_along(input_covariates), function(index) {
    input <- input_covariates[[index]]
    if (!is.list(input)) stop("Each Lego covariate must be an object.")
    name <- lego_text(input$name, paste0("covariates[", index, "].name"), "^[A-Z][A-Z0-9_]*$", 24L)
    if (name %in% LEGO_RESERVED_COVARIATES || grepl("^(ETA|EPS)[0-9]*$", name)) {
      stop("Reserved Lego covariate name: ", name)
    }
    type <- if (is.null(input$type)) {
      "continuous"
    } else {
      lego_text(input$type, paste0("covariates[", index, "].type"), "^(continuous|categorical)$", 16L)
    }
    form <- if (identical(type, "continuous") && identical(input$form, "linear-additive")) "linear-additive" else "power"
    reference <- lego_number(
      input$reference,
      paste0("covariates[", index, "].reference"),
      if (identical(type, "continuous") && identical(form, "power")) 1e-12 else -1e12,
      1e12
    )
    comparison_default <- if (identical(type, "continuous")) reference * 1.25 else if (reference == 0) 1 else 0
    comparison <- lego_number(
      input$comparison,
      paste0("covariates[", index, "].comparison"),
      if (identical(type, "continuous")) 1e-12 else -1e12,
      1e12,
      default = comparison_default
    )
    if (identical(type, "categorical") && comparison == reference) {
      stop("Categorical Lego covariate reference and comparison must differ: ", name)
    }
    list(
      name = name,
      type = type,
      form = form,
      scope = lego_text(input$scope %||% "patient", paste0("covariates[", index, "].scope"), "^(patient|administration)$", 16L),
      target = lego_text(input$target, paste0("covariates[", index, "].target"), "^[A-Za-z_][A-Za-z0-9_]*$", 96L),
      reference = reference,
      comparison = comparison,
      maximum = if (is.null(input$maximum)) NULL else lego_number(input$maximum, paste0("covariates[", index, "].maximum"), reference + 1e-12, 1e12),
      beta = lego_number(input$beta, paste0("covariates[", index, "].beta"), -10, 10)
    )
  })
  covariate_names <- vapply(covariates, `[[`, character(1), "name")
  effect_keys <- vapply(covariates, function(item) paste(item$name, item$target, sep = "::"), character(1))
  if (anyDuplicated(effect_keys)) stop("Lego covariate effects must target a parameter only once.")
  for (name in unique(covariate_names)) {
    definitions <- Filter(function(item) identical(item$name, name), covariates)
    signatures <- vapply(definitions, function(item) paste(item$type, item$scope, item$reference, item$comparison, sep = "::"), character(1))
    if (length(unique(signatures)) > 1L) stop("Repeated Lego covariates must share type, scope, reference and comparison: ", name)
  }

  population_input <- specification$population %||% list()
  if (!is.list(population_input)) stop("The Lego population model must be an object.")
  named_numbers <- function(values, field, lower, upper, key_pattern) {
    if (is.null(values)) return(NULL)
    if (is.atomic(values) && !is.null(names(values))) values <- as.list(values)
    if (!is.list(values) || is.null(names(values)) || length(values) > 210L) stop(field, " must be a named object.")
    if (any(!grepl(key_pattern, names(values), perl = TRUE)) || anyDuplicated(names(values))) stop("Invalid ", field, " parameter name.")
    stats::setNames(vapply(seq_along(values), function(index) lego_number(values[[index]], paste0(field, ".", names(values)[index]), lower, upper), numeric(1)), names(values))
  }
  named_distributions <- function(values) {
    if (is.null(values)) return(NULL)
    if (is.atomic(values) && !is.null(names(values))) values <- as.list(values)
    if (!is.list(values) || is.null(names(values)) || length(values) > 210L) stop("iivDistributions must be a named object.")
    if (any(!grepl("^[A-Za-z_][A-Za-z0-9_]*$", names(values))) || anyDuplicated(names(values))) stop("Invalid iivDistributions parameter name.")
    stats::setNames(vapply(seq_along(values), function(index) {
      lego_text(values[[index]], paste0("iivDistributions.", names(values)[index]), "^(normal|lognormal|logit)$", 16L)
    }, character(1)), names(values))
  }
  residual_input <- population_input$residualError %||% list(type = "combined", additive = 0.1, proportional = 0.2)
  if (!is.list(residual_input)) stop("Residual error must be an object.")
  residual <- list(
    type = lego_text(residual_input$type %||% "combined", "residualError.type", "^(additive|proportional|combined)$", 16L),
    additive = lego_number(residual_input$additive, "residualError.additive", 1e-9, 1e6, default = 0.1),
    proportional = lego_number(residual_input$proportional, "residualError.proportional", 1e-9, 1e3, default = 0.2)
  )
  population <- list(
    iivVariances = named_numbers(population_input$iivVariances, "iivVariances", 1e-9, 1e6, "^[A-Za-z_][A-Za-z0-9_]*$"),
    iivDistributions = named_distributions(population_input$iivDistributions),
    iivCovariances = named_numbers(population_input$iivCovariances, "iivCovariances", -1e6, 1e6, "^[A-Za-z_][A-Za-z0-9_]*::[A-Za-z_][A-Za-z0-9_]*$"),
    residualError = residual
  )

  list(version = version, nodes = nodes, edges = edges, covariates = covariates, population = population)
}

lego_spec_from_code <- function(code) {
  if (!is.character(code) || length(code) != 1L || nchar(code, type = "bytes") > 200000L) {
    stop("The pasted Lego model is empty or too large.")
  }
  lines <- strsplit(code, "\n", fixed = TRUE)[[1]]
  marker <- lines[startsWith(trimws(lines), LEGO_SPEC_PREFIX)][1]
  if (is.na(marker) || !nzchar(marker)) {
    stop("Public compilation accepts only models generated by Atelier Lego. The Lego specification marker is missing.")
  }
  encoded <- substring(trimws(marker), nchar(LEGO_SPEC_PREFIX) + 1L)
  if (!nzchar(encoded) || nchar(encoded, type = "bytes") > 100000L) stop("Invalid Lego specification marker.")
  decoded <- utils::URLdecode(encoded)
  specification <- tryCatch(
    jsonlite::fromJSON(decoded, simplifyVector = FALSE),
    error = function(error) stop("Invalid Lego specification JSON: ", conditionMessage(error))
  )
  normalize_lego_spec(specification)
}

lego_format_number <- function(value) {
  formatted <- format(as.numeric(value), digits = 15, scientific = FALSE, trim = TRUE, nsmall = 0)
  sub("\\.$", "", formatted)
}

lego_model_code <- function(specification) {
  spec <- normalize_lego_spec(specification)
  nodes <- spec$nodes
  edges <- spec$edges
  ids <- vapply(nodes, `[[`, integer(1), "id")
  internal <- stats::setNames(
    vapply(nodes, function(node) paste0("L", node$id, "_", node$name), character(1)),
    as.character(ids)
  )
  name_for <- function(id) unname(internal[[as.character(id)]])
  client_name_for <- function(id) nodes[[match(id, ids)]]$name
  q_name_for <- function(edge, client = FALSE) {
    resolver <- if (isTRUE(client)) client_name_for else name_for
    paste0("q_", paste(sort(c(resolver(edge$from), resolver(edge$to))), collapse = "_"))
  }
  volume_nodes <- Filter(function(node) node$kind %in% LEGO_VOLUME_KINDS, nodes)
  central_nodes <- Filter(function(node) identical(node$kind, "central"), nodes)
  observed <- if (length(central_nodes)) central_nodes[[1]] else volume_nodes[[1]]
  dosed <- Filter(function(node) is.finite(node$dose) && node$dose > 0, nodes)
  mass_nodes <- Filter(function(node) !node$kind %in% LEGO_PD_KINDS, nodes)
  adm <- if (length(dosed)) dosed[[1]] else if (length(mass_nodes)) mass_nodes[[1]] else nodes[[1]]
  node_for <- function(id) nodes[[match(id, ids)]]
  fraction_parameter <- function(node) {
    source <- if (is.null(node$fractionComplementOf)) node else node_for(node$fractionComplementOf)
    paste0("f_", name_for(source$id))
  }
  fraction_expression <- function(node) {
    parameter <- fraction_parameter(node)
    if (is.null(node$fractionComplementOf)) parameter else paste0("(1-", parameter, ")")
  }
  duration_parameter <- function(node) {
    if (is.null(node$inputDurationTlagOf)) paste0("tk0_", name_for(node$id)) else paste0("tlag_", name_for(node$inputDurationTlagOf))
  }
  route_doses <- length(dosed) > 1L || any(vapply(dosed, function(node) identical(node$inputType, "zero_order") || node$tlag > 0, logical(1)))

  parameters <- list()
  for (edge in edges) {
    from <- name_for(edge$from)
    to <- if (identical(edge$to, "OUT")) "e" else name_for(edge$to)
    client_suffix <- paste0(client_name_for(edge$from), "_", if (identical(edge$to, "OUT")) "e" else client_name_for(edge$to))
    if (edge$kinetics %in% c("michaelis_menten", "hill")) {
      parameters <- c(parameters, list(
        list(name = paste0("vmax_", from, "_", to), client_name = paste0("vmax_", client_suffix), value = edge$vmax, note = paste0("maximum transfer rate from ", from, " to ", to), iiv = FALSE),
        list(name = paste0("km_", from, "_", to), client_name = paste0("km_", client_suffix), value = edge$km, note = paste0("Michaelis constant from ", from, " to ", to), iiv = FALSE)
      ))
      if (identical(edge$kinetics, "hill")) {
        parameters[[length(parameters) + 1L]] <- list(
          name = paste0("gamma_", from, "_", to), client_name = paste0("gamma_", client_suffix), value = edge$gamma,
          note = paste0("Hill exponent from ", from, " to ", to), iiv = FALSE
        )
      }
    } else if (identical(edge$eliminationParameterization, "clearance")) {
      parameters[[length(parameters) + 1L]] <- list(
        name = paste0("cl_", from), client_name = paste0("cl_", client_name_for(edge$from)), value = edge$cl,
        note = paste0("clearance from ", from), iiv = TRUE
      )
    } else if (identical(edge$transferParameterization, "clearance")) {
      client_name <- q_name_for(edge, client = TRUE)
      existing <- if (length(parameters)) vapply(parameters, `[[`, character(1), "client_name") else character()
      if (!client_name %in% existing) {
        parameters[[length(parameters) + 1L]] <- list(
          name = q_name_for(edge), client_name = client_name, value = edge$q,
          note = paste0("intercompartmental clearance between ", from, " and ", to), iiv = FALSE
        )
      }
    } else {
      parameters[[length(parameters) + 1L]] <- list(
        name = paste0("k_", from, "_", to),
        client_name = paste0("k_", client_suffix),
        value = edge$k,
        note = if (identical(edge$to, "OUT")) paste0("elimination from ", from) else paste0("transfer ", from, " to ", to),
        iiv = identical(edge$to, "OUT")
      )
    }
  }
  for (node in volume_nodes) {
    name <- name_for(node$id)
    parameters[[length(parameters) + 1L]] <- list(
      name = paste0("v_", name), client_name = paste0("v_", node$name), value = node$vol, note = paste0("volume of ", name), iiv = identical(node$kind, "central")
    )
  }
  for (node in nodes) {
    name <- name_for(node$id)
    for (key in lego_advanced_fields(node)) parameters[[length(parameters) + 1L]] <- list(
      name = paste0(key, "_", name), client_name = paste0(key, "_", node$name),
      value = node[[key]], note = paste(key, name), iiv = FALSE
    )
    if (node$dose > 0) {
      if (node$tlag > 0) parameters[[length(parameters) + 1L]] <- list(name = paste0("tlag_", name), client_name = paste0("tlag_", node$name), value = node$tlag, note = paste0("administration lag to ", name), iiv = FALSE)
      if (identical(node$inputType, "zero_order") && is.null(node$inputDurationTlagOf)) parameters[[length(parameters) + 1L]] <- list(name = paste0("tk0_", name), client_name = paste0("tk0_", node$name), value = node$inputDuration, note = paste0("administration duration to ", name), iiv = FALSE)
      if (length(dosed) > 1L && is.null(node$fractionComplementOf)) parameters[[length(parameters) + 1L]] <- list(name = paste0("f_", name), client_name = paste0("f_", node$name), value = node$doseFraction / 100, note = paste0("dose fraction to ", name), iiv = FALSE)
    }
    if (identical(node$kind, "effect")) {
      parameters[[length(parameters) + 1L]] <- list(name = paste0("ke0_", name), client_name = paste0("ke0_", node$name), value = node$ke0, note = paste0("effect equilibration for ", name), iiv = FALSE)
    }
    if (identical(node$kind, "response")) {
      parameters <- c(parameters, list(
        list(name = paste0("kin_", name), client_name = paste0("kin_", node$name), value = node$kin, note = paste0("production of ", name), iiv = FALSE),
        list(name = paste0("kout_", name), client_name = paste0("kout_", node$name), value = node$kout, note = paste0("degradation of ", name), iiv = FALSE),
        list(name = paste0("smax_", name), client_name = paste0("smax_", node$name), value = node$smax, note = paste0("maximum effect on ", name), iiv = FALSE),
        list(name = paste0("sc50_", name), client_name = paste0("sc50_", node$name), value = node$sc50, note = paste0("half effect concentration for ", name), iiv = FALSE)
      ))
    }
  }
  if (!length(parameters)) stop("The Lego model does not expose any pharmacometric parameter.")

  client_parameter_names <- vapply(parameters, `[[`, character(1), "client_name")
  if (anyDuplicated(client_parameter_names)) stop("Lego parameter targets must be unique.")
  covariate_effects <- lapply(seq_along(spec$covariates), function(index) {
    covariate <- spec$covariates[[index]]
    target_index <- match(covariate$target, client_parameter_names)
    if (is.na(target_index)) stop("Unknown Lego covariate target: ", covariate$target)
    c(covariate, list(
      target_name = parameters[[target_index]]$name,
      beta_name = paste0("BETA_", covariate$name, "_", index)
    ))
  })
  data_covariates <- covariate_effects[!duplicated(vapply(covariate_effects, `[[`, character(1), "name"))]

  variance_input <- spec$population$iivVariances
  random_parameters <- if (is.null(variance_input)) {
    Filter(function(parameter) isTRUE(parameter$iiv), parameters)
  } else {
    unknown <- setdiff(names(variance_input), client_parameter_names)
    if (length(unknown)) stop("Unknown Lego random-effect parameter: ", unknown[[1]])
    parameters[client_parameter_names %in% names(variance_input)]
  }
  variance_for <- function(parameter) if (is.null(variance_input)) 0.09 else unname(variance_input[parameter$client_name])
  distribution_input <- spec$population$iivDistributions
  if (!is.null(distribution_input)) {
    unknown <- setdiff(names(distribution_input), client_parameter_names)
    if (length(unknown)) stop("Unknown Lego random-effect distribution parameter: ", unknown[[1]])
  }
  distribution_for <- function(parameter) {
    distribution <- if (!is.null(distribution_input) && parameter$client_name %in% names(distribution_input)) {
      unname(distribution_input[[parameter$client_name]])
    } else if (parameter$value > 0) "lognormal" else "normal"
    if (identical(distribution, "logit") && (parameter$value <= 0 || parameter$value >= 1)) {
      stop("A logit-normal Lego parameter must have a typical value strictly between 0 and 1: ", parameter$client_name)
    }
    distribution
  }
  covariance_input <- spec$population$iivCovariances %||% numeric()
  covariance_for <- function(left, right) {
    key <- paste(sort(c(left$client_name, right$client_name)), collapse = "::")
    if (key %in% names(covariance_input)) unname(covariance_input[[key]]) else 0
  }
  if (length(random_parameters)) {
    omega <- matrix(0, length(random_parameters), length(random_parameters))
    diag(omega) <- vapply(random_parameters, variance_for, numeric(1))
    for (row in seq_along(random_parameters)) for (column in seq_len(row - 1L)) {
      omega[row, column] <- omega[column, row] <- covariance_for(random_parameters[[row]], random_parameters[[column]])
    }
    if (any(!is.finite(omega)) || any(eigen(omega, symmetric = TRUE, only.values = TRUE)$values <= 1e-12)) stop("The Lego Omega matrix must be positive definite.")
  }
  eta_index <- stats::setNames(seq_along(random_parameters), vapply(random_parameters, `[[`, character(1), "name"))
  width <- max(8L, nchar(paste0("TV_", vapply(parameters, `[[`, character(1), "name"))))
  pad <- function(text) sprintf("%-*s", width, text)
  marker_spec <- spec
  for (field in c("iivVariances", "iivDistributions", "iivCovariances")) {
    if (!is.null(marker_spec$population[[field]])) marker_spec$population[[field]] <- as.list(marker_spec$population[[field]])
  }
  spec_json <- jsonlite::toJSON(marker_spec, auto_unbox = TRUE, null = "null", digits = 10)
  marker <- paste0(LEGO_SPEC_PREFIX, utils::URLencode(spec_json, reserved = TRUE))
  lines <- c(marker, if (route_doses) "$PLUGIN evtools", "$PARAM @annotated")

  for (parameter in parameters) {
    lines <- c(lines, paste0(pad(paste0("TV_", parameter$name)), " : ", lego_format_number(parameter$value), " : ", parameter$note))
  }
  for (covariate in covariate_effects) {
    effect_label <- if (identical(covariate$type, "categorical")) "categorical effect of " else if (identical(covariate$form, "linear-additive")) "additive linear contribution of " else "power effect of "
    lines <- c(lines, paste0(covariate$beta_name, " : ", lego_format_number(covariate$beta), " : ", effect_label, covariate$name, " on ", covariate$target))
  }
  for (index in seq_along(random_parameters)) {
    parameter <- random_parameters[[index]]
    lines <- c(lines, paste0(pad(paste0("ETA", index)), " : 0 : individual effect on ", parameter$name))
  }
  if (length(covariate_effects)) {
    lines <- c(lines, "", "$PARAM @covariates @annotated")
    for (covariate in data_covariates) {
      description <- if (identical(covariate$type, "categorical")) {
        paste0("categorical covariate, reference ", lego_format_number(covariate$reference), ", affected category ", lego_format_number(covariate$comparison))
      } else {
        "continuous covariate, reference value"
      }
      scope <- if (identical(covariate$scope, "administration")) " [administration]" else ""
      lines <- c(lines, paste0(covariate$name, " : ", lego_format_number(covariate$reference), " : ", description, scope))
    }
  }

  lines <- c(lines, "", "$OMEGA @block")
  if (length(random_parameters)) {
    for (row in seq_along(random_parameters)) {
      values <- vapply(seq_len(row), function(column) if (row == column) variance_for(random_parameters[[row]]) else covariance_for(random_parameters[[row]], random_parameters[[column]]), numeric(1))
      lines <- c(lines, paste(lego_format_number(values), collapse = " "))
    }
  } else lines <- c(lines, "0 FIX")
  residual <- spec$population$residualError
  lines <- c(
    lines,
    "",
    "$SIGMA @annotated"
  )
  lines <- c(
    lines,
    paste0("PROP : ", lego_format_number(if (residual$type %in% c("proportional", "combined")) residual$proportional^2 else 0), " : proportional residual variance"),
    paste0("ADD : ", lego_format_number(if (residual$type %in% c("additive", "combined")) residual$additive^2 else 0), " : additive residual variance")
  )
  lines <- c(lines, "", "$CMT @annotated")
  if (route_doses) {
    lines <- c(lines, paste0(pad("LEGO_INPUT"), " : split dose input [ADM]"))
  }
  for (node in nodes) {
    name <- name_for(node$id)
    tags <- c(if (!route_doses && identical(node$id, adm$id)) "ADM", if (identical(node$id, observed$id)) "OBS")
    tag_text <- if (length(tags)) paste0(" [", paste(tags, collapse = ", "), "]") else ""
    lines <- c(lines, paste0(pad(name), " : model state ", name, tag_text))
  }

  lines <- c(lines, "", "$MAIN")
  if (route_doses) lines <- c(lines, "F_LEGO_INPUT = 0;")
  for (parameter in parameters) {
    index <- unname(eta_index[parameter$name])
    effects <- Filter(function(covariate) identical(covariate$target_name, parameter$name), covariate_effects)
    multiplicative_effects <- Filter(function(covariate) !identical(covariate$form, "linear-additive"), effects)
    additive_effects <- Filter(function(covariate) identical(covariate$form, "linear-additive"), effects)
    effect_code <- paste0(vapply(
      multiplicative_effects,
      function(covariate) {
        if (identical(covariate$type, "categorical")) {
          paste0(" * exp(", covariate$beta_name, " * (", covariate$name, " == ", lego_format_number(covariate$comparison), "))")
        } else {
          paste0(" * pow(", covariate$name, "/", lego_format_number(covariate$reference), ", ", covariate$beta_name, ")")
        }
      },
      character(1)
    ), collapse = "")
    additive_code <- paste0(vapply(additive_effects, function(covariate) {
      source <- if (is.null(covariate$maximum)) covariate$name else paste0("fmin(", covariate$name, ", ", lego_format_number(covariate$maximum), ")")
      paste0(" + ", covariate$beta_name, " * (", source, " - ", lego_format_number(covariate$reference), ")")
    }, character(1)), collapse = "")
    typical <- if (length(additive_effects)) paste0("(TV_", parameter$name, effect_code, additive_code, ")") else paste0("TV_", parameter$name, effect_code)
    if (length(index) && is.finite(index)) {
      random <- paste0("ETA", index, " + ETA(", index, ")")
      distribution <- distribution_for(parameter)
      if (identical(distribution, "normal")) {
        individual <- paste0(typical, " + ", random)
      } else if (identical(distribution, "logit")) {
        link_effects <- paste0(vapply(multiplicative_effects, function(covariate) {
          if (identical(covariate$type, "categorical")) {
            paste0(" + ", covariate$beta_name, " * (", covariate$name, " == ", lego_format_number(covariate$comparison), ")")
          } else {
            paste0(" + ", covariate$beta_name, " * log(", covariate$name, "/", lego_format_number(covariate$reference), ")")
          }
        }, character(1)), collapse = "")
        individual <- paste0("1.0/(1.0 + exp(-(log(TV_", parameter$name, "/(1.0-TV_", parameter$name, "))", link_effects, " + ", random, ")))")
      } else {
        individual <- paste0(typical, " * exp(", random, ")")
      }
    } else {
      individual <- typical
    }
    lines <- c(lines, paste0("double ", parameter$name, " = ", individual, ";"))
  }
  response_nodes <- Filter(function(node) identical(node$kind, "response"), nodes)
  for (node in response_nodes) {
    name <- name_for(node$id)
    lines <- c(lines, paste0(name, "_0 = kin_", name, "/kout_", name, ";"))
  }

  driver_concentration <- function(node) {
    source <- nodes[[match(node$source, ids)]]
    source_name <- name_for(source$id)
    if (source$kind %in% LEGO_VOLUME_KINDS) paste0("(", source_name, "/v_", source_name, ")") else source_name
  }
  special <- function(node) lego_advanced_expression(node, driver_concentration(node), name_for(node$id), function(key) paste0(key, "_", name_for(node$id)))
  for (node in Filter(function(n) n$kind %in% c("tumor", "interaction"), nodes)) {
    lines <- c(lines, paste0(name_for(node$id), "_0 = ", special(node)$initial, ";"))
  }
  mass_terms <- function(node) {
    incoming <- Filter(function(edge) !identical(edge$to, "OUT") && identical(edge$to, node$id), edges)
    outgoing <- Filter(function(edge) identical(edge$from, node$id), edges)
    base_flux <- function(edge) {
      from <- name_for(edge$from)
      to <- if (identical(edge$to, "OUT")) "e" else name_for(edge$to)
      amount <- from
      if (identical(edge$kinetics, "hill")) {
        gamma <- paste0("gamma_", from, "_", to)
        return(paste0("vmax_", from, "_", to, "*pow(", amount, ", ", gamma, ")/(pow(km_", from, "_", to, ", ", gamma, ") + pow(", amount, ", ", gamma, "))"))
      }
      if (identical(edge$kinetics, "michaelis_menten")) {
        return(paste0("vmax_", from, "_", to, "*", amount, "/(km_", from, "_", to, " + ", amount, ")"))
      }
      if (identical(edge$eliminationParameterization, "clearance")) {
        return(paste0("cl_", from, "*", amount, "/v_", from))
      }
      if (identical(edge$transferParameterization, "clearance")) {
        return(paste0(q_name_for(edge), "*", amount, "/v_", from))
      }
      paste0("k_", from, "_", to, "*", amount)
    }
    flux <- function(edge) {
      modifiers <- Filter(function(n) n$kind == "interaction" && n$targetFrom == edge$from && identical(n$targetTo, edge$to), nodes)
      paste0(base_flux(edge), paste0(vapply(modifiers, function(n) paste0("*(", special(n)$factor, ")"), character(1)), collapse = ""))
    }
    terms <- c(
      vapply(incoming, function(edge) paste0("+ ", flux(edge)), character(1)),
      vapply(outgoing, function(edge) paste0("- ", flux(edge)), character(1))
    )
    if (length(terms)) paste(terms, collapse = " ") else "0"
  }

  lines <- c(lines, "", "$ODE")
  if (route_doses) lines <- c(lines, "dxdt_LEGO_INPUT = 0;")
  for (node in nodes) {
    name <- name_for(node$id)
    if (node$kind %in% c("tumor", "interaction")) {
      lines <- c(lines, paste0("dxdt_", name, " = ", special(node)$derivative, ";"))
    } else if (identical(node$kind, "effect")) {
      lines <- c(lines, paste0("dxdt_", name, " = ke0_", name, "*", "(", driver_concentration(node), " - ", name, ");"))
    } else if (identical(node$kind, "response")) {
      driver <- driver_concentration(node)
      lines <- c(lines, paste0("dxdt_", name, " = kin_", name, "*(1 + smax_", name, "*", driver, "/(sc50_", name, " + ", driver, ")) - kout_", name, "*", name, ";"))
    } else {
      lines <- c(lines, paste0("dxdt_", name, " = ", mass_terms(node), ";"))
    }
  }

  if (route_doses) {
    lines <- c(lines, "", "$EVENT", "if ((EVID == 1 || EVID == 4) && CMT == 1) {")
    for (index in seq_along(dosed)) {
      node <- dosed[[index]]
      name <- name_for(node$id)
      amount <- if (length(dosed) > 1L) paste0("AMT*", fraction_expression(node)) else "AMT"
      cmt <- match(node$id, ids) + 1L
      event_name <- paste0("route_", index)
      if (identical(node$inputType, "zero_order")) {
        lines <- c(lines, paste0("  evt::ev ", event_name, " = evt::infuse(", amount, ", ", cmt, ", ", amount, "/", duration_parameter(node), ");"))
      } else {
        lines <- c(lines, paste0("  evt::ev ", event_name, " = evt::bolus(", amount, ", ", cmt, ");"))
      }
      if (node$tlag > 0) lines <- c(lines, paste0("  evt::retime(", event_name, ", TIME + tlag_", name, ");"))
      lines <- c(lines, paste0("  self.push(", event_name, ");"))
    }
    lines <- c(lines, "}")
  }

  lines <- c(lines, "", "$TABLE")
  for (node in volume_nodes) {
    name <- name_for(node$id)
    lines <- c(lines, paste0("double CONC_", name, " = ", name, "/v_", name, ";"))
  }
  observed_name <- name_for(observed$id)
  for (node in Filter(function(n) n$kind == "interaction", nodes)) {
    lines <- c(lines, paste0("double MOD_", name_for(node$id), " = ", special(node)$factor, ";"))
  }
  lines <- c(
    lines,
    paste0("double IPRED = CONC_", observed_name, ";"),
    switch(residual$type,
      additive = "double DV = IPRED + EPS(2);",
      proportional = "double DV = IPRED * (1 + EPS(1));",
      "double DV = IPRED * (1 + EPS(1)) + EPS(2);"
    ),
    "if (DV < 0) DV = 0;",
    "",
    "$CAPTURE @annotated",
    "DV : simulated concentration with residual error"
  )
  for (node in volume_nodes) {
    name <- name_for(node$id)
    lines <- c(lines, paste0("CONC_", name, " : concentration in ", name))
  }
  for (node in Filter(function(n) n$kind == "interaction", nodes)) {
    lines <- c(lines, paste0("MOD_", name_for(node$id), " : flux divided by baseline flux"))
  }
  paste(lines, collapse = "\n")
}

safe_lego_model_code <- function(code = NULL, specification = NULL) {
  spec <- if (!is.null(specification)) normalize_lego_spec(specification) else lego_spec_from_code(code)
  lego_model_code(spec)
}

compile_model <- function(
  model_id = NULL,
  custom_code = NULL,
  allow_custom = FALSE,
  custom_soloc = NULL,
  custom_cache = NULL
) {
  if (!is.null(custom_code)) {
    if (!nzchar(trimws(custom_code))) stop("Paste a complete mrgsolve model before compiling.")
    safe_lego <- !isTRUE(allow_custom)
    if (safe_lego) custom_code <- safe_lego_model_code(code = custom_code)
    if (is.null(custom_soloc) || !dir.exists(custom_soloc)) {
      stop("A session-scoped compilation directory is required for custom models.")
    }
    key <- safe_model_id(paste0("custom_", short_hash(custom_code)))
    if (!is.null(custom_cache) && !is.null(custom_cache[[key]])) return(custom_cache[[key]])
    if (safe_lego && !is.null(custom_cache) && length(ls(custom_cache, all.names = TRUE)) >= 5L) {
      stop("This session has reached the limit of five distinct Lego model compilations.")
    }
    model <- mrgsolve::mcode(key, custom_code, soloc = custom_soloc, quiet = TRUE)
    if (!is.null(custom_cache)) custom_cache[[key]] <- model
    return(model)
  }

  row <- model_record(model_id)
  key <- row$id[[1]]
  if (!is.null(.model_cache[[key]])) return(.model_cache[[key]])
  model <- mrgsolve::mread(
    model = key,
    project = MODEL_ROOT,
    file = row$file[[1]],
    soloc = MODEL_CACHE_DIR,
    quiet = TRUE
  )
  .model_cache[[key]] <- model
  model
}

model_param_names <- function(model) names(as.list(mrgsolve::param(model)))

safe_param <- function(model, values) {
  if (is.null(values) || !length(values)) return(model)
  values <- as.list(values)
  values <- values[names(values) %in% model_param_names(model)]
  if (!length(values)) return(model)
  mrgsolve::param(model, values)
}

tagged_compartment <- function(model, tag) {
  annotations <- tryCatch(model@annot$data, error = function(error) data.frame())
  if (!nrow(annotations) || !all(c("block", "name", "options") %in% names(annotations))) return(NA_integer_)
  rows <- annotations[
    annotations$block == "CMT" & grepl(paste0("(^|,\\s*)", tag, "(\\s*,|$)"), annotations$options),
    ,
    drop = FALSE
  ]
  if (!nrow(rows)) return(NA_integer_)
  match(rows$name[[1]], model@cmtL)
}

model_capture_names <- function(model) {
  simulation <- model |>
    mrgsolve::zero_re() |>
    mrgsolve::mrgsim(end = 0, delta = 1) |>
    as.data.frame()
  setdiff(names(simulation), c("ID", "time"))
}

validate_model_contract <- function(model) {
  errors <- character()
  warnings <- character()
  adm <- tagged_compartment(model, "ADM")
  obs <- tagged_compartment(model, "OBS")
  captures <- tryCatch(model_capture_names(model), error = function(error) character())
  omega <- tryCatch(as.matrix(mrgsolve::omat(model)), error = function(error) matrix(numeric(), 0, 0))
  sigma <- tryCatch(as.matrix(mrgsolve::smat(model)), error = function(error) matrix(numeric(), 0, 0))

  if (!is.finite(adm)) errors <- c(errors, "No [ADM] compartment tag was found.")
  if (!is.finite(obs)) errors <- c(errors, "No [OBS] compartment tag was found.")
  if (!length(captures)) errors <- c(errors, "The model does not expose any captured prediction.")
  if (!"DV" %in% captures) warnings <- c(warnings, "DV is not captured; another concentration column will be used.")
  if (!nrow(omega)) errors <- c(errors, "OMEGA is missing; MAP estimation requires random effects.")
  if (!nrow(sigma)) errors <- c(errors, "SIGMA is missing; MAP estimation requires a residual error model.")
  if (!all(paste0("ETA", seq_len(nrow(omega))) %in% model_param_names(model))) {
    errors <- c(errors, "Each OMEGA term must have a matching ETA1, ETA2, ... parameter for mapbayr.")
  }

  list(
    ok = !length(errors),
    errors = unique(errors),
    warnings = unique(warnings),
    adm_cmt = adm,
    obs_cmt = obs,
    captures = captures,
    n_eta = nrow(omega),
    n_sigma = nrow(sigma)
  )
}

pick_concentration_column <- function(data) {
  preferred <- c("DV", "CP", "CONC", "CONC_PLASMA", "CONCENTRATION")
  found <- preferred[preferred %in% names(data)][1]
  if (is.na(found)) stop("No concentration output found. Capture DV, CP or CONC in the model.")
  found
}
