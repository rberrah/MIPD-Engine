pmetrics_unavailable <- function(...) stop("Pmetrics is unavailable on this server.")

pmetrics_artifact_catalog <- function() data.frame(
  id = character(), label = character(), labelEn = character(),
  stringsAsFactors = FALSE
)

pmetrics_artifact_path <- pmetrics_unavailable
read_pmetrics_artifact <- pmetrics_unavailable
build_pmetrics_model <- pmetrics_unavailable
pmetrics_error_spec <- pmetrics_unavailable
compile_pmetrics_bridge <- pmetrics_unavailable
pmetrics_compartment_contract <- pmetrics_unavailable
fit_one_pmetrics_model <- pmetrics_unavailable
pmetrics_residual_values <- pmetrics_unavailable
