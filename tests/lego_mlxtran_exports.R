if (!requireNamespace("lixoftConnectors", quietly = TRUE)) {
  cat("MONOLIX_NATIVE_STATUS=not-executed; lixoftConnectors is unavailable.\n")
  quit(status = 0L)
}
suppressPackageStartupMessages({ library(jsonlite); library(lixoftConnectors); library(mrgsolve) })
initialized <- tryCatch({
  isTRUE(initializeLixoftConnectors("simulx", path = Sys.getenv("MONOLIX_HOME", "C:/Program Files/Lixoft/MonolixSuite2024R1")))
}, error = function(error) {
  cat("MONOLIX_NATIVE_STATUS=not-executed; ", conditionMessage(error), "\n", sep = "")
  FALSE
})
if (!initialized) {
  cat("MONOLIX_NATIVE_STATUS=not-executed; Simulx could not be initialized with a usable installation and licence.\n")
  quit(status = 0L)
}
fixture_path <- "test-results/lego-export-audit.json"
if (!file.exists(fixture_path)) stop("Run node scripts/audit_lego_roundtrip.mjs before the native Monolix harness.")
cases <- fromJSON(fixture_path, simplifyVector = FALSE)
work <- tempfile("lego-mlxtran-")
dir.create(work)
failures <- character()
for (case in Filter(function(x) x$format == "mlxtran", cases)) {
  tryCatch({
    path <- file.path(work, paste0(case$id, ".txt"))
    writeLines(case$code, path)
    newProject(modelFile = path)
    pop <- getPopulationElements()[[1]]$data
    for (line in strsplit(case$code, "\n", fixed = TRUE)[[1]]) {
      match <- regmatches(line, regexec("^;[[:space:]]+([A-Za-z_][A-Za-z_0-9]*) = ([-0-9.e]+)", line))[[1]]
      if (length(match) == 3 && match[[2]] %in% names(pop)) pop[[match[[2]]]] <- as.numeric(match[[3]])
    }
    pop[grepl("^omega_", names(pop))] <- 0
    if ("b" %in% names(pop)) pop$b <- 0.2
    definePopulationElement("auditParameters", pop)
    group <- getGroups()[[1]]$name
    setGroupElement(group, "auditParameters")
    setGroupSize(group, 1)
    if (length(case$spec$covariates)) {
      cov <- setNames(lapply(case$spec$covariates, function(x) x$reference), vapply(case$spec$covariates, function(x) x$name, character(1)))
      defineCovariateElement("auditCovariates", as.data.frame(cov[!duplicated(names(cov))]))
      setGroupElement(group, "auditCovariates")
    }
    horizon <- case$spec$simulation$horizon
    dosed <- Filter(function(n) n$dose > 0, case$spec$nodes)
    defineTreatmentElement("auditDoses", list(admID = 1, data = data.frame(time = c(0, horizon/2), amount = dosed[[1]]$dose)))
    central <- Filter(function(n) n$kind == "central", case$spec$nodes)[[1]]
    times <- seq(horizon/500, horizon, length.out = 500)
    defineOutputElement("auditOutput", list(data = data.frame(time = times), output = paste0("C_", central$name)))
    setGroupElement(group, c("auditDoses", "auditOutput"))
    runSimulation()
    result <- getSimulationResults()
    prediction <- result$res[[paste0("C_", central$name)]]
    if (is.null(prediction) || nrow(prediction) != length(times)) stop("Simulx returned no requested concentration profile")
    cpp <- Filter(function(x) identical(x$preset, case$preset) && x$format == "mrgsolve", cases)[[1]]
    mod <- mcode(paste0(case$id, "_reference"), cpp$code, soloc = work, quiet = TRUE)
    cmt <- if ("LEGO_INPUT" %in% names(mod@init)) 1 else match(dosed[[1]]$name, names(mod@init))
    reference <- as.data.frame(mrgsim(zero_re(mod), events = ev(amt = dosed[[1]]$dose, cmt = cmt, ii = horizon/2, addl = 1), tgrid = times, recsort = 3, obsonly = TRUE))
    x <- prediction[[paste0("C_", central$name)]]
    y <- reference[[paste0("CONC_", central$name)]]
    stopifnot(length(x) == length(y), all(is.finite(x)), max(x) > 0)
    error <- max(abs(x - y)) / max(1, abs(y))
    if (error > 1e-4) stop("Simulx/mrgsolve discrepancy: ", error)
    message("PASS native MLXTRAN ", case$id)
  }, error = function(e) {
    failures <<- c(failures, paste(case$id, conditionMessage(e)))
    message("FAIL ", tail(failures, 1))
  })
}
unlink(work, recursive = TRUE)
if (length(failures)) stop(paste(failures, collapse = "\n"))
cat("MONOLIX_NATIVE_STATUS=passed; native Simulx comparisons completed.\n")
