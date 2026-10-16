# Run from the repository root; add --deploy after reviewing the manifest.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
with_pmetrics <- "--with-pmetrics" %in% args
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (with_pmetrics) options(repos = c(LAPKB = "https://lapkb.r-universe.dev", CRAN = getOption("repos")[["CRAN"]]))

tracked <- system2("git", c("ls-files"), stdout = TRUE)
files <- tracked[grepl("^(app\\.R|R/[^/]+\\.R|www/[^/]+\\.(css|js)|models/[^/]+\\.(cpp|json)|pmetrics/(artifacts|schema)/[^/]+\\.json|ml/registry\\.json|ml/artifacts/[^/]+\\.rds)$", tracked)]
if (!with_pmetrics) {
  files <- files[files != "R/pmetrics_engine.R" & !grepl("^pmetrics/", files)]
  Sys.setenv(ENABLE_PMETRICS = "false")
} else {
  files <- files[files != "R/pmetrics_disabled.R"]
  Sys.setenv(ENABLE_PMETRICS = "true")
}
stopifnot("app.R" %in% files, "models/catalog.json" %in% files,
          all(file.exists(file.path(root, files))))
deps <- rsconnect::appDependencies(root, appFiles = files)
required <- c("deSolve", "minpack.lm", "mrgsolve", "mapbayr")
if (with_pmetrics) required <- c(required, "Pmetrics")
stopifnot(all(required %in% deps$Package))
if (!with_pmetrics) stopifnot(!"Pmetrics" %in% deps$Package)
if ("--deploy" %in% args) {
  Sys.setenv(ALLOW_CUSTOM_MODELS = "false")
  rsconnect::deployApp(appDir = root, appFiles = files, appName = "MIPD_Engine", launch.browser = FALSE)
}
