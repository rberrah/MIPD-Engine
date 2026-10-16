# Run scripts/audit_lego_roundtrip.mjs first. Only generated preset data is used.
suppressPackageStartupMessages({ library(jsonlite); library(mrgsolve) })
cases <- fromJSON("test-results/lego-export-audit.json", simplifyVector = FALSE)
work <- tempfile("lego-native-")
dir.create(work)
failures <- character()
for (case in cases) {
  tryCatch({
    if (case$format == "mrgsolve") {
      model <- mcode(case$id, case$code, soloc = work, quiet = TRUE)
      variants <- Filter(function(x) identical(x$preset, case$preset), cases)
      for (variant in variants) {
        other <- mcode(paste0(variant$id, "_import"), variant$reimportedMrgsolve, soloc = work, quiet = TRUE)
        run <- function(mod, spec, comparison) {
          names <- names(mod@init)
          dosed <- Filter(function(n) n$dose > 0, spec$nodes)
          cmt <- if ("LEGO_INPUT" %in% names) match("LEGO_INPUT", names) else match(dosed[[1]]$name, names)
          # Dose and horizon are data, not part of a plain structural model export.
          horizon <- case$spec$simulation$horizon
          if (is.null(horizon)) horizon <- 24
          doses <- ev(amt = case$spec$nodes[[which(vapply(case$spec$nodes, function(n) n$dose > 0, logical(1)))[1]]]$dose,
                      cmt = cmt, ii = horizon / 2, addl = 1)
          if (comparison && length(spec$covariates)) {
            # Comparison settings are simulation data and are not in plain code.
            cov <- setNames(lapply(case$spec$covariates, function(x) x$comparison), vapply(case$spec$covariates, function(x) x$name, character(1)))
            mod <- param(mod, cov[!duplicated(names(cov))])
          }
          as.data.frame(mrgsim(zero_re(mod), events = doses, end = horizon, delta = horizon / 500, recsort = 3))
        }
        for (comparison in c(FALSE, TRUE)) {
          a <- run(model, case$spec, comparison)
          b <- run(other, variant$reimportedSpec, comparison)
          stopifnot(identical(a$time, b$time))
          # Match states case-insensitively; names are upper-case in NM-TRAN.
          states <- c(intersect(tolower(names(model@init)), tolower(names(other@init))), "dv",
                      tolower(grep("^CONC_", names(a), value = TRUE)))
          for (state in states) {
            x <- a[[match(state, tolower(names(a)))]]
            y <- b[[match(state, tolower(names(b)))]]
            stopifnot(length(x) == length(y), all(is.finite(x)), max(abs(x - y)) < 1e-6 * max(1, abs(x)))
          }
        }
        message("PASS native mrgsolve ", variant$id)
      }
    }
  }, error = function(e) {
    failures <<- c(failures, paste(case$id, conditionMessage(e)))
    message("FAIL ", tail(failures, 1))
  })
}
unlink(work, recursive = TRUE)
if (length(failures)) stop(paste(failures, collapse = "\n"))
