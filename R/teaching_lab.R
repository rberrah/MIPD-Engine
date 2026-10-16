# Only bounded synthetic lab parameters cross this bridge; no observations or code.
teaching_lab_import <- function(spec) {
  labs <- c("distribution", "accumulation", "absorption", "infusion")
  if (!is.list(spec) || !identical(spec$type, "pk-teaching-lab") ||
      !isTRUE(spec$version == 1) || length(spec$lab) != 1 || !spec$lab %in% labs) stop("Invalid teaching laboratory.")
  ranges <- list(
    dose = c(1, 1000), cl = c(0.1, 30), vc = c(5, 100), q = c(0, 50), vp = c(5, 200),
    tau = c(1, 48), count = c(1, 20), loading = c(1, 4), end = c(12, 336),
    ka = c(0.01, 10), f = c(0, 1), duration = c(0.1, 72)
  )
  fields <- switch(spec$lab,
    distribution = c("dose", "cl", "vc", "q", "vp", "tau", "count", "loading", "end"),
    accumulation = c("dose", "cl", "vc", "q", "vp", "tau", "count", "loading", "end"),
    absorption = c("dose", "cl", "vc", "ka", "f", "end"),
    infusion = c("dose", "cl", "vc", "duration", "end")
  )
  p <- spec$parameters
  if (!is.list(p) || !setequal(names(p), fields) || anyDuplicated(names(p))) stop("Invalid teaching parameters.")
  for (key in fields) {
    v <- p[[key]]
    if (!is.numeric(v) || length(v) != 1 || !is.finite(v) || v < ranges[[key]][1] || v > ranges[[key]][2]) stop(paste("Invalid teaching parameter:", key))
  }
  if (identical(spec$lab, "accumulation") && p$count != floor(p$count)) stop("Dose count must be an integer.")
  two <- identical(spec$lab, "distribution")
  repeated <- identical(spec$lab, "accumulation")
  doses <- data.frame(time = if (repeated) (seq_len(p$count) - 1) * p$tau else 0,
    amount = p$dose, interval = if (repeated) p$tau else p$end, count = 1L,
    infusion = if (identical(spec$lab, "infusion")) p$duration else 0, ss = 0L,
    status = "administered", time_uncertainty = 0)
  if (repeated) doses$amount[1] <- p$dose * p$loading
  central_dose <- if (spec$lab %in% c("distribution", "accumulation", "infusion")) doses$amount[1] else 0
  nodes <- list(list(id = 1, kind = "central", name = "CENT", vol = p$vc, dose = central_dose))
  edges <- list(list(from = 1, to = "OUT", k = p$cl / p$vc, kinetics = "first_order", eliminationParameterization = "clearance", cl = p$cl))
  if (two) {
    nodes[[2]] <- list(id = 2, kind = "periph", name = "PERI", vol = p$vp, dose = 0)
    edges <- c(edges, list(list(from = 1, to = 2, k = p$q / p$vc), list(from = 2, to = 1, k = p$q / p$vp)))
  } else if (identical(spec$lab, "absorption")) {
    nodes[[2]] <- list(id = 2, kind = "depot", name = "DEPOT", dose = p$dose)
    edges <- c(edges, list(list(from = 2, to = 1, k = p$ka * p$f), list(from = 2, to = "OUT", k = p$ka * (1 - p$f))))
  } else if (identical(spec$lab, "infusion")) {
    nodes[[1]]$inputType <- "zero_order"
    nodes[[1]]$inputDuration <- p$duration
  }
  code <- safe_lego_model_code(specification = list(version = 3, nodes = nodes, edges = edges, covariates = list()))
  list(code = code, doses = doses, route = if (identical(spec$lab, "absorption")) "Oral" else "IV")
}
