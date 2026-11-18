# Original, session-only A/B comparison of *already-ranked* dose candidates.
# This does not change MAP fitting, model weights, recommendation ranking or
# source models. It reuses the same simulation engine as the dosing workflow.
regimen_ab_candidates <- function(result, max_candidates = 20L, max_interval_h = 336) {
  if (is.null(result) || !is.data.frame(result$recommendations)) return(data.frame())
  candidates <- result$recommendations
  required <- c("dose", "interval", "infusion")
  if (!all(required %in% names(candidates))) return(data.frame())
  usable <- is.finite(candidates$dose) & candidates$dose > 0 &
    is.finite(candidates$interval) & candidates$interval > 0 &
    candidates$interval <= max_interval_h &
    is.finite(candidates$infusion) & candidates$infusion >= 0 &
    candidates$infusion <= candidates$interval
  candidates <- candidates[usable, required, drop = FALSE]
  candidates <- unique(candidates)
  candidates <- head(candidates, as.integer(max_candidates))
  if (!nrow(candidates)) return(data.frame())
  candidates$rank <- seq_len(nrow(candidates))
  rownames(candidates) <- NULL
  candidates
}

regimen_ab_compare <- function(result, a = 1L, b = 2L) {
  candidates <- regimen_ab_candidates(result)
  a <- suppressWarnings(as.integer(a))
  b <- suppressWarnings(as.integer(b))
  if (nrow(candidates) < 2L || length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || a == b ||
      a < 1L || b < 1L || a > nrow(candidates) || b > nrow(candidates)) return(NULL)

  selected <- candidates[c(a, b), , drop = FALSE]
  horizon <- max(24, selected$interval)
  # A cap keeps the optional display from multiplying the parent fitting cost.
  # This is a *display* grid; clinical candidate ranking uses its own settings.
  grid_step <- 0.25
  profiles <- vector("list", 2L)
  metrics <- vector("list", 2L)
  for (i in seq_len(2L)) {
    regimen <- selected[i, , drop = FALSE]
    output <- simulate_averaged_regimen(
      fits = result$fits, weights = result$weights,
      dose = regimen$dose[[1]], interval = regimen$interval[[1]],
      infusion = regimen$infusion[[1]], horizon = horizon,
      delta = grid_step,
      concentration_threshold = result$target_options$concentration_threshold %||% NA_real_
    )
    profiles[[i]] <- transform(output$profile, scenario = c("A", "B")[[i]])
    metrics[[i]] <- data.frame(
      scenario = c("A", "B")[[i]],
      dose_mg = regimen$dose[[1]],
      interval_h = regimen$interval[[1]],
      infusion_h = regimen$infusion[[1]],
      auc0_24 = output$metrics$auc24[[1]],
      cmin_interval = output$metrics$cmin[[1]],
      cmax_interval = output$metrics$cmax[[1]],
      time_above_pct = output$metrics$time_above_pct[[1]]
    )
  }
  list(profiles = do.call(rbind, profiles), metrics = do.call(rbind, metrics),
       horizon = horizon, grid_step = grid_step, candidates = selected)
}
