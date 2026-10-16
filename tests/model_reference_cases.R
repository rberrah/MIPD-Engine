suppressPackageStartupMessages({
  library(mrgsolve)
  library(jsonlite)
})

arguments <- commandArgs(trailingOnly = FALSE)
file_argument <- sub("^--file=", "", grep("^--file=", arguments, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file_argument), ".."), winslash = "/", mustWork = TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local = TRUE)

# These are deterministic software-regression cases, not independent
# reproductions of the source articles. Population covariates stay at the model
# defaults, ETAs and residual errors are zero, and concentrations are sampled at
# 1, 4, 12 and 24 h after a single administration.
cases <- list(
  list(id = "amox_carlier", dose = 1000, rate = 0,
    expected = c(21.84032290331795, 7.400625190145801, 0.630142232944522, 0.015659561827004)),
  list(id = "vanco_pkjust", dose = 1000, rate = 0,
    expected = c(6.503660611652944, 5.924335680997364, 4.619458756670076, 3.180667326606865)),
  list(id = "vanco_goti", dose = 1000, rate = 0,
    expected = c(14.48030909761354, 9.75027207278402, 5.62221726420194, 3.523574723211349)),
  list(id = "vanco_roberts", dose = 1000, rate = 1000 / 24,
    expected = c(0.357550684630095, 1.355530487122397, 3.54177729758866, 5.832931321652911)),
  list(id = "tacrolimus_woillard_ddi", dose = 5, rate = 0,
    expected = c(20.63937886394410, 9.551734546808834, 5.723997609759323, 3.301936934475965)),
  list(id = "dapto_dvorchik_IV_adults", dose = 500, rate = 0,
    expected = c(51.41574065140168, 27.96362497629674, 14.44638508584297, 5.431270811207022)),
  list(id = "cefepime_jonckheere", dose = 2000, rate = 4000,
    expected = c(76.06237029265395, 37.45561193565089, 14.29672286563470, 3.518472907650924)),
  list(id = "rifampicine_jing", dose = 450, rate = 0,
    expected = c(5.963798392839188, 6.147878990875527, 3.531830992973112, 1.532981243434618))
)

times <- c(1, 4, 12, 24)
recording <- any(vapply(cases, function(case) is.null(case$expected), logical(1)))

for (case in cases) {
  record <- model_record(case$id)
  route <- model_routes(record)[[1]]
  model <- mrgsolve::zero_re(compile_model(model_id = case$id))
  compartment <- match(model_administration_cmt(record, route), model@cmtL)
  event <- mrgsolve::ev(amt = case$dose, cmt = compartment, rate = case$rate)
  output <- as.data.frame(mrgsolve::mrgsim(
    model, event, tgrid = times, obsonly = TRUE, recsort = 3,
    atol = 1e-12, rtol = 1e-12
  ))
  column <- pick_concentration_column(output)
  actual <- unname(output[[column]])

  if (is.null(case$expected)) {
    cat(case$id, " = c(", paste(format(actual, digits = 16, scientific = FALSE), collapse = ", "), ")\n", sep = "")
  } else if (!isTRUE(all.equal(actual, case$expected, tolerance = 1e-7, check.attributes = FALSE))) {
    stop(case$id, " deterministic reference profile changed: ",
      paste(format(actual, digits = 12), collapse = ", "))
  }
}

if (recording) quit(status = 2L)
cat("Deterministic model reference cases passed: ", length(cases), ".\n", sep = "")
