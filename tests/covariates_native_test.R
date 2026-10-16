suppressPackageStartupMessages(library(mrgsolve))
root <- commandArgs(trailingOnly = TRUE)[[1]]
fixtures <- jsonlite::fromJSON(file.path(root, "fixtures.json"))
for (index in seq_len(nrow(fixtures))) {
  case <- fixtures[index, ]
  model <- mrgsolve::mread(case$name, project = root, soloc = root, quiet = TRUE)
  model <- mrgsolve::zero_re(model)
  out <- model |> mrgsolve::ev(amt = 100, cmt = 1) |> mrgsolve::mrgsim(end = 24, delta = 0.2)
  stopifnot(max(abs(out$CL - case$cl)) < 1e-7)
  expected <- 100 / case$volume * exp(-case$cl / case$volume * out$time)
  after_dose <- out$time > 0
  stopifnot(max(abs(out$IPRED[after_dose] - expected[after_dose])) < 1e-5)
  if (case$name == "lesson") {
    out <- model |> mrgsolve::param(WT = 35, GENO = 0) |>
      mrgsolve::ev(amt = 100, cmt = 1) |> mrgsolve::mrgsim(end = 24, delta = 0.1)
    stopifnot(abs(out$CL[[1]] - 4 * (35/70)^0.75) < 1e-8, abs(out$V[[1]] - 15) < 1e-8)
  }
}
cat("Native mrgsolve covariate exports and course example match analytic IV predictions.\n")
