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

# The SHA embedded in a production release must describe ALL shipped tracked
# content. An uncommitted local edit would make a truthful commit fingerprint
# impossible. GitHub Actions checkouts are clean; manual local deployment must
# satisfy the same rule. Untracked files are not included in the rsconnect
# manifest, which is built exclusively from git ls-files.
assert_clean_mipd_release_checkout <- function(root) {
  original_dir <- getwd()
  on.exit(setwd(original_dir), add = TRUE)
  setwd(root)
  status <- suppressWarnings(system2(
    "git", c("status", "--porcelain", "--untracked-files=no"),
    stdout = TRUE, stderr = TRUE
  ))
  exit_status <- attr(status, "status")
  if (!is.null(exit_status) && exit_status != 0L) {
    stop("Cannot check the Git state before publishing a release.")
  }
  if (length(status)) {
    stop("Cannot publish an immutable Git SHA with modified tracked files: ",
         paste(status, collapse = "; "))
  }
  invisible(TRUE)
}
