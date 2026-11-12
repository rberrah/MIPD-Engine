# MAP-BE diagnostic positioning, research/teaching only.
# Never infer an unreported absolute PK parameter distribution from an ETA.
# The reference is the *marginal population prior*, not the posterior nor an
# empirical distribution of observed patients.

map_be_eta_position <- function(eta, omega, labels = NULL, selected = NULL,
                                posterior_covariance = NULL, condition_limit = 1e8) {
  omega <- as.matrix(omega)
  if (!is.numeric(omega) || length(dim(omega)) != 2L ||
      nrow(omega) != ncol(omega) || !nrow(omega) ||
      any(!is.finite(omega))) return(NULL)
  n <- nrow(omega)
  if (any(abs(omega - t(omega)) > 1e-7)) return(NULL)
  eta_values <- suppressWarnings(as.numeric(eta))
  selected <- if (is.null(selected)) seq_len(n) else suppressWarnings(as.integer(selected))
  if (!length(selected) || anyNA(selected) || any(selected < 1L | selected > n) ||
      anyDuplicated(selected)) return(NULL)
  if (length(eta_values) == n) {
    full_eta <- eta_values
  } else if (length(eta_values) == length(selected)) {
    full_eta <- rep(NA_real_, n)
    full_eta[selected] <- eta_values
  } else return(NULL)
  if (is.null(labels) || length(labels) != n) labels <- paste0("ETA", seq_len(n))
  labels <- as.character(labels)
  labels[is.na(labels) | !nzchar(labels)] <- paste0("ETA", which(is.na(labels) | !nzchar(labels)))
  # A zero-variance or unselected ETA is not an estimable stochastic parameter.
  retained <- selected[is.finite(full_eta[selected]) & diag(omega)[selected] > 0]
  if (!length(retained)) return(NULL)
  sd_prior <- sqrt(diag(omega)[retained])
  z <- full_eta[retained] / sd_prior
  q <- stats::qnorm(0.95)
  output <- data.frame(
    index = retained,
    parameter = make.unique(labels[retained]),
    eta = full_eta[retained],
    prior_sd = sd_prior,
    z = z,
    percentile = 100 * stats::pnorm(z),
    prior_5 = -q * sd_prior,
    prior_95 = q * sd_prior,
    posterior_sd = rep(NA_real_, length(retained)),
    posterior_5 = rep(NA_real_, length(retained)),
    posterior_95 = rep(NA_real_, length(retained)),
    stringsAsFactors = FALSE
  )
  post <- tryCatch(as.matrix(posterior_covariance), error = function(e) matrix(numeric(), 0, 0))
  valid_cov <- is.numeric(post) && nrow(post) == ncol(post) && nrow(post) > 0L &&
    all(is.finite(post)) && max(abs(post - t(post))) < 1e-7
  if (valid_cov) {
    eigenvalues <- tryCatch(eigen((post + t(post)) / 2, symmetric = TRUE,
      only.values = TRUE)$values, error = function(e) numeric())
    valid_cov <- length(eigenvalues) == nrow(post) && all(is.finite(eigenvalues)) &&
      min(eigenvalues) > 0 && max(eigenvalues) / min(eigenvalues) <= condition_limit
  }
  if (valid_cov) {
    if (nrow(post) == n) {
      post_sd <- sqrt(diag(post)[retained])
    } else if (nrow(post) == length(selected)) {
      post_sd <- sqrt(diag(post)[match(retained, selected)])
    } else {
      post_sd <- numeric()
    }
    if (length(post_sd) == nrow(output) && all(is.finite(post_sd))) {
      output$posterior_sd <- post_sd
      output$posterior_5 <- output$eta - q * post_sd
      output$posterior_95 <- output$eta + q * post_sd
    }
  }
  output
}

map_be_weighted_quantile <- function(values, weights, p) {
  keep <- is.finite(values) & is.finite(weights) & weights >= 0
  values <- values[keep]; weights <- weights[keep]
  if (!length(values) || sum(weights) <= 0) return(NA_real_)
  order <- order(values); values <- values[order]; weights <- weights[order] / sum(weights)
  values[which(cumsum(weights) >= p)[1]]
}

map_be_pmetrics_position <- function(prior, posterior, parameter_names) {
  if (!is.data.frame(prior) || !is.data.frame(posterior) ||
      !length(parameter_names) || !all(c(parameter_names, "prob") %in% names(prior)) ||
      !all(c(parameter_names, "prob") %in% names(posterior))) return(NULL)
  wp <- as.numeric(prior$prob); wq <- as.numeric(posterior$prob)
  if (any(!is.finite(wp) | wp < 0) || any(!is.finite(wq) | wq < 0) ||
      sum(wp) <= 0 || sum(wq) <= 0) return(NULL)
  wp <- wp / sum(wp); wq <- wq / sum(wq)
  rows <- lapply(parameter_names, function(name) {
    xp <- suppressWarnings(as.numeric(prior[[name]]))
    xq <- suppressWarnings(as.numeric(posterior[[name]]))
    if (any(!is.finite(xp)) || any(!is.finite(xq))) return(NULL)
    patient <- sum(xq * wq)
    data.frame(
      parameter = name, prior_mean = sum(xp * wp), patient = patient,
      percentile = 100 * sum(wp[xp <= patient]),
      prior_5 = map_be_weighted_quantile(xp, wp, 0.05),
      prior_95 = map_be_weighted_quantile(xp, wp, 0.95),
      stringsAsFactors = FALSE
    )
  })
  if (any(vapply(rows, is.null, logical(1)))) return(NULL)
  table <- do.call(rbind, rows)
  points <- do.call(rbind, lapply(parameter_names, function(name) {
    rbind(
      data.frame(parameter = name, value = as.numeric(prior[[name]]),
        prob = wp, distribution = "Prior", stringsAsFactors = FALSE),
      data.frame(parameter = name, value = as.numeric(posterior[[name]]),
        prob = wq, distribution = "Posterior", stringsAsFactors = FALSE)
    )
  }))
  list(table = table, points = points)
}

map_be_diagnostic <- function(fit) {
  if (is.null(fit) || inherits(fit, "tdm_fit_error")) {
    return(list(type = "unavailable", reason = "fit_error"))
  }
  if (identical(fit$engine %||% "mapbayr", "pmetrics")) {
    if (!isTRUE(fit$pmetrics$mapped)) return(list(type = "unavailable", reason = "population"))
    if (!exists("read_pmetrics_artifact", mode = "function", inherits = TRUE))
      return(list(type = "unavailable", reason = "pmetrics_backend"))
    artifact <- tryCatch(read_pmetrics_artifact(fit$pmetrics$artifact_id),
      error = function(e) NULL)
    if (is.null(artifact)) return(list(type = "unavailable", reason = "pmetrics_backend"))
    prior <- artifact$support
    posterior <- fit$posterior_points %||% NULL
    parameter_names <- as.character(artifact$parameters$name)
    result <- map_be_pmetrics_position(prior, posterior, parameter_names)
    if (is.null(result)) return(list(type = "unavailable", reason = "incompatible"))
    return(c(list(type = "pmetrics", label = fit$label %||% fit$id,
                  observations = nrow(fit$source_observations %||% data.frame())),
             result))
  }
  if (is.null(fit$estimate)) return(list(type = "unavailable", reason = "population"))
  omega <- tryCatch(as.matrix(mrgsolve::omat(fit$model)), error = function(e) NULL)
  if (is.null(omega) || !nrow(omega))
    return(list(type = "unavailable", reason = "no_omega"))
  eta <- fit$estimate$final_eta[[1]] %||% NULL
  selected <- fit$estimate$arg.optim$select_eta %||% NULL
  covariance <- fit$estimate$covariance[[1]] %||% NULL
  labels <- rownames(omega)
  # mrgsolve keeps @annotated IIV labels on omat(), otherwise show ETA index.
  if (is.null(labels)) {
    labels <- tryCatch(unlist(mrgsolve::labels(mrgsolve::omat(fit$model)),
                              use.names = FALSE), error = function(e) NULL)
  }
  rows <- map_be_eta_position(eta, omega, labels, selected, covariance)
  if (is.null(rows)) return(list(type = "unavailable", reason = "incompatible"))
  list(type = "mapbayr", label = fit$label %||% fit$id,
       table = rows,
       observations = nrow(fit$source_observations %||% data.frame()),
       eta_count = nrow(rows),
       covariance_available = all(is.finite(rows$posterior_sd)),
       nonestimated = nrow(omega) - nrow(rows))
}

map_be_diagnostic_plot <- function(data, english = FALSE) {
  if (identical(data$type, "mapbayr")) {
    rows <- data$table
    curves <- do.call(rbind, lapply(seq_len(nrow(rows)), function(i) {
      row <- rows[i, , drop = FALSE]
      limit <- max(3.8 * row$prior_sd, abs(row$eta) * 1.12)
      x <- seq(-limit, limit, length.out = 180L)
      data.frame(parameter = row$parameter, eta = x,
        density = stats::dnorm(x, sd = row$prior_sd))
    }))
    return(ggplot2::ggplot(curves, ggplot2::aes(x = eta, y = density)) +
      ggplot2::geom_line(color = "#176b70", linewidth = 1) +
      ggplot2::geom_vline(data = rows,
        ggplot2::aes(xintercept = eta), inherit.aes = FALSE,
        color = "#b24c2c", linewidth = 1, linetype = "dashed") +
      ggplot2::facet_wrap(~parameter, scales = "free", ncol = 3) +
      ggplot2::labs(
        x = if (english) "ETA (dimensionless)" else "ETA (sans unité)",
        y = if (english) "Marginal prior density" else "Densité a priori marginale",
        subtitle = if (english) "Green: model prior N(0, omega); dashed: individual MAP estimate"
                   else "Vert : prior du modèle N(0, omega) ; tirets : estimation MAP individuelle") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank()))
  }
  if (identical(data$type, "pmetrics")) {
    return(ggplot2::ggplot(data$points,
         ggplot2::aes(x = value, y = prob, color = distribution)) +
      ggplot2::geom_segment(ggplot2::aes(xend = value, yend = 0),
        linewidth = 0.65, alpha = 0.75) +
      ggplot2::geom_point(size = 2.6) +
      ggplot2::geom_vline(data = data$table,
        ggplot2::aes(xintercept = patient), inherit.aes = FALSE,
        color = "#17202a", linewidth = .8, linetype = "dashed") +
      ggplot2::scale_color_manual(values = c(Prior = "#176b70", Posterior = "#b24c2c")) +
      ggplot2::facet_wrap(~parameter, scales = "free_x", ncol = 3) +
      ggplot2::labs(
        x = if (english) "Parameter (model units)" else "Paramètre (unités du modèle)",
        y = if (english) "Support-point probability" else "Probabilité du point de support",
        color = NULL,
        subtitle = if (english) "Weighted support points; dashed: posterior mean, not a Gaussian MAP"
                   else "Points de support pondérés ; tirets : moyenne postérieure, pas de MAP gaussien") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                     legend.position = "top"))
  }
  NULL
}
