# Standalone public Shiny deployment to the EXISTING tdmhub app.
# Normal invocation validates dependencies without touching a live service.
# --verify-bundle tests production packaging offline; --deploy uploads it.
# --deploy bundles an immutable, public release fingerprint so post-deploy
# monitoring can distinguish a tested Git SHA from whatever is actually live.
run_mipd_deployment <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  root <- normalizePath(".", winslash = "/", mustWork = TRUE)
  with_pmetrics <- "--with-pmetrics" %in% args
  do_deploy <- "--deploy" %in% args
  verify_bundle <- "--verify-bundle" %in% args
  bundle_release <- do_deploy || verify_bundle
  if (bundle_release && with_pmetrics) {
    stop("Pmetrics deployment is unsupported on the current Shinyapps image.")
  }
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  if (with_pmetrics) {
    options(repos = c(LAPKB = "https://lapkb.r-universe.dev",
                     CRAN = getOption("repos")[["CRAN"]]))
  }

  tracked <- system2("git", c("ls-files"), stdout = TRUE)
  files <- tracked[grepl(
    "^(app\\.R|R/[^/]+\\.R|www/[^/]+\\.(css|js)|models/[^/]+\\.(cpp|json)|integration/(tdm-guide-(recommendations|clinical-presets)\\.v1\\.json|tdm-drugs/[a-z][a-z0-9_-]*\\.json)|pmetrics/(artifacts|schema)/[^/]+\\.json|ml/registry\\.json|ml/artifacts/[^/]+\\.rds)$",
    tracked
  )]
  if (!with_pmetrics) {
    files <- files[files != "R/pmetrics_engine.R" & !grepl("^pmetrics/", files)]
    Sys.setenv(ENABLE_PMETRICS = "false")
  } else {
    files <- files[files != "R/pmetrics_disabled.R"]
    Sys.setenv(ENABLE_PMETRICS = "true")
  }
  stopifnot("app.R" %in% files, "R/deployment_profile.R" %in% files,
            "models/catalog.json" %in% files,
            all(file.exists(file.path(root, files))))
  # The public Shiny runtime reads these data files from APP_ROOT on startup.
  # Verify the *actual bundle manifest*, not only the working checkout:
  # without these JSON artifacts a green CI could still deploy a broken app.
  catalogue <- jsonlite::fromJSON(file.path(root, "models", "catalog.json"))$models
  drug_keys <- sort(unique(as.character(catalogue$drugKey)))
  expected_drug_files <- file.path("integration", "tdm-drugs",
                                   paste0(drug_keys, ".json"))
  included_drug_files <- grep(
    "^integration/tdm-drugs/[a-z][a-z0-9_-]*\\.json$", files, value = TRUE
  )
  required_guide_files <- c(
    "integration/tdm-guide-recommendations.v1.json",
    "integration/tdm-guide-clinical-presets.v1.json"
  )
  stopifnot(
    length(drug_keys) == 23L,
    identical(sort(included_drug_files), sort(expected_drug_files)),
    all(required_guide_files %in% files),
    all(file.exists(file.path(root, c(required_guide_files, expected_drug_files))))
  )

  release_path <- file.path(root, "www", "deployment-release.json")
  if (bundle_release) {
    # Fail closed if the local checkout contains modifications that would
    # make the published release SHA inconsistent with the deployed source.
    source(file.path(root, "R", "deployment_profile.R"), local = TRUE)
    assert_clean_mipd_release_checkout(root)
    if (file.exists(release_path)) {
      stop("Refusing to overwrite pre-existing local deployment-release.json.")
    }
    sha <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
    if (length(sha) != 1L || !grepl("^[a-f0-9]{40}$", sha)) {
      stop("Cannot identify an immutable MIPD commit for release.")
    }
    release <- list(
      schema = "mipd-public-release.v1",
      target = "tdmhub/MIPD_Engine",
      appId = "17776160",
      sourceSha = sha[[1]],
      modelCount = 46L,
      allowCustomModels = FALSE,
      enablePmetrics = FALSE
    )
    source(file.path(root, "R", "deployment_profile.R"), local = TRUE)
    validate_mipd_release(release)
    jsonlite::write_json(release, release_path, auto_unbox = TRUE, pretty = TRUE)
    on.exit(unlink(release_path), add = TRUE)
    stopifnot(file.exists(release_path))
    files <- c(files, "www/deployment-release.json")
    Sys.setenv(ALLOW_CUSTOM_MODELS = "false", ENABLE_PMETRICS = "false")
    cat("Bundled public MIPD release commit: ", sha[[1]], "\n", sep = "")
  } else {
    # A bare manifest check never bundles or creates a release fingerprint.
    if (file.exists(release_path)) {
      stop("Unexpected deployment-release.json in manifest-only validation.")
    }
  }

  deps <- rsconnect::appDependencies(root, appFiles = files)
  required <- c("deSolve", "minpack.lm", "mrgsolve", "mapbayr")
  if (with_pmetrics) required <- c(required, "Pmetrics")
  stopifnot(all(required %in% deps$Package))
  if (!with_pmetrics) stopifnot(!"Pmetrics" %in% deps$Package)

  if (verify_bundle) cat("Dry-run public release bundle PASS: ", sha[[1]], "\n", sep = "")

  if (do_deploy) {
    account <- Sys.getenv("SHINYAPPS_ACCOUNT")
    if (!nzchar(account)) stop("SHINYAPPS_ACCOUNT is required for deployment.")
    if (!identical(account, "tdmhub")) {
      stop("Refusing deployment to unexpected Shinyapps account: ", account)
    }
    rsconnect::deployApp(
      appDir = root, appFiles = files,
      appId = "17776160", appName = "MIPD_Engine", account = account,
      server = "shinyapps.io", launch.browser = FALSE
    )
    cat("Published MIPD commit: ", sha[[1]], "\n", sep = "")
  }
  invisible(files)
}
run_mipd_deployment()
