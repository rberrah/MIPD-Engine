suppressPackageStartupMessages({ library(mrgsolve); library(jsonlite) })
APP_ROOT <- normalizePath("tdm-engine", winslash = "/", mustWork = TRUE)
source("tdm-engine/R/model_library.R")
source("tdm-engine/R/teaching_lab.R")
run <- function() {
  directory <- tempfile("pk-lab-test-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  fixture <- file.path(directory, "synthetic.json")
  status <- system2("node", c("scripts/test_labs.mjs", shQuote(fixture)))
  stopifnot(status == 0)
  fixtures <- jsonlite::fromJSON(fixture, simplifyVector = FALSE)
  stopifnot(abs(as.numeric(lego_format_number(1e-9)) - 1e-9) < 1e-20)
  for (index in seq_along(fixtures)) {
    f <- fixtures[[index]]
    imported <- teaching_lab_import(list(type = "pk-teaching-lab", version = 1, lab = f$lab, parameters = f$parameters))
    expected <- do.call(rbind, lapply(f$points, as.data.frame))
    dose_table <- do.call(rbind, lapply(f$doses, as.data.frame))
    stopifnot(
      isTRUE(all.equal(imported$doses$time, dose_table$time)),
      isTRUE(all.equal(imported$doses$amount, dose_table$amount)),
      isTRUE(all.equal(imported$doses$infusion, dose_table$infusion)),
      identical(imported$route, if (f$lab == "absorption") "Oral" else "IV")
    )
    for (variant in c("export", "server")) {
      code <- if (variant == "export") f$code else imported$code
      mod <- mrgsolve::mcode(paste0("lab", index, variant), code, soloc = directory, quiet = TRUE)
      cmt <- if (f$lab == "absorption" && variant == "server") 2L else 1L
      rate <- ifelse(dose_table$infusion > 0, dose_table$amount / dose_table$infusion, 0)
      events <- mrgsolve::as.ev(data.frame(time = dose_table$time, amt = dose_table$amount, cmt = cmt, evid = 1, rate = rate))
      result <- as.data.frame(mrgsolve::mrgsim(mrgsolve::ev(mrgsolve::zero_re(mod), events), tgrid = expected$t, obsonly = TRUE, recsort = 3, atol = 1e-10, rtol = 1e-10))
      stopifnot(nrow(result) == nrow(expected))
      concentration <- if ("IPRED" %in% names(result)) result$IPRED else result$DV
      stopifnot(max(abs(concentration - expected$c)) < 1e-5)
      if (f$lab == "distribution" && max(abs(result$L2_PERI - expected$peripheral)) >= 1e-5) {
        print(data.frame(time = result$time, actual = result$L2_PERI, expected = expected$peripheral))
        stop("Peripheral mismatch in scenario ", index, ", ", variant)
      }
    }
  }
  invalid <- list(type = "pk-teaching-lab", version = 1, lab = "accumulation", parameters = fixtures[[6]]$parameters)
  for (key in c("dose", "count", "cl")) {
    bad <- invalid; bad$parameters[[key]] <- -1
    stopifnot(inherits(try(teaching_lab_import(bad), silent = TRUE), "try-error"))
  }
  message("Teaching laboratories: 32 native/regenerated mrgsolve comparisons and input validation passed.")
}
run()
