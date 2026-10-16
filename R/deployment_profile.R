# A deployment fingerprint is generated only for the fixed public Shinyapps target.
# The marker is shipped inside the app bundle and cannot be overridden by
# environment variables supplied by the hosting process.
validate_mipd_release <- function(profile) {
  stopifnot(is.list(profile))
  if (!identical(profile$schema, "mipd-public-release.v1") ||
      !identical(profile$target, "tdmhub/MIPD_Engine") ||
      !identical(as.character(profile$appId), "17776160") ||
      !identical(profile$modelCount, 46L) ||
      !isFALSE(profile$allowCustomModels) ||
      !isFALSE(profile$enablePmetrics) ||
      !is.character(profile$sourceSha) ||
      length(profile$sourceSha) != 1L ||
      !grepl("^[0-9a-f]{40}$", profile$sourceSha)) {
    stop("Invalid or unsafe MIPD public deployment fingerprint.")
  }
  invisible(TRUE)
}

enforce_mipd_deployment_profile <- function(root) {
  release_path <- file.path(root, "www", "deployment-release.json")
  if (!file.exists(release_path)) return(invisible(FALSE))
  release <- jsonlite::fromJSON(release_path, simplifyVector = TRUE)
  validate_mipd_release(release)
  # An explicit public deployment marker ALWAYS dominates host environment.
  Sys.setenv(ALLOW_CUSTOM_MODELS = "false", ENABLE_PMETRICS = "false")
  invisible(TRUE)
}
