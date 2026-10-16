arguments <- commandArgs(trailingOnly = FALSE)
path <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
source(file.path(dirname(path), "..", "R", "ddi_engine.R"))
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
time <- seq(0, 240, by = 0.1)
c <- ifelse(time < 120, 2, 0)
for (type in c("tdi", "turnover_induction")) {
  values <- ddi_modifier(type, c, time, 0, 120, 0.5, 3, 1, kdeg = 0.02, kinact = 0.1)
  if (type == "tdi") {
    loss <- 0.02 + 0.1 * 2 / 3
    equilibrium <- 0.02 / loss
    expected <- equilibrium + (1 - equilibrium) * exp(-loss * time[time <= 120])
    stopifnot(values[1201] < values[length(values)], values[1201] < 0.3)
  } else {
    expected <- 3 + (1 - 3) * exp(-0.02 * time[time <= 120])
    stopifnot(values[1201] > values[length(values)], values[1201] > 2.5)
  }
  stopifnot(max(abs(values[time <= 120] - expected)) < 1e-10)
  stopifnot(all(ddi_modifier(type, rep(0,length(time)), time, 0, 120, 1, 3, 1) == 1))
}
stopifnot(abs(ddi_modifier("reversible", 1, 0, 0, 1, 1, 1, 1) - 0.5) < 1e-12)
stopifnot(abs(ddi_modifier("hill_inhibition", 1, 0, 0, 1, 1, 0.8, 1, hill = 3) - 0.6) < 1e-12)
cat("DDI mechanisms: analytical turnover, recovery, zero exposure, Ki and Hill passed.\n")
