args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[[1]] else "core"
contract <- jsonlite::fromJSON("runtime/reference-runtime.json", simplifyVector = TRUE)

stopifnot(identical(as.character(getRversion()), contract$rVersion))

check_versions <- function(expected) {
  for (pkg in names(expected)) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing required package: ", pkg)
    actual <- as.character(utils::packageVersion(pkg))
    if (!identical(actual, unname(expected[[pkg]]))) {
      stop("Runtime drift for ", pkg, ": expected ", expected[[pkg]], ", got ", actual)
    }
    cat(pkg, actual, "OK\n")
  }
}

check_lockfile <- function(path = "renv.lock") {
  if (!file.exists(path)) stop("Missing runtime lockfile: ", path)
  lock <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (!identical(lock$R$Version, contract$rVersion)) {
    stop("Lockfile R drift: expected ", contract$rVersion, ", got ", lock$R$Version)
  }

  repository_urls <- vapply(lock$R$Repositories, `[[`, character(1), "URL")
  portable_snapshot <- sub("/__linux__/noble/", "/", contract$coreCranSnapshot, fixed = TRUE)
  if (!any(repository_urls %in% c(contract$coreCranSnapshot, portable_snapshot))) {
    stop("Lockfile does not use the dated CRAN snapshot.")
  }

  locked <- lock$Packages
  for (pkg in names(contract$directPackages)) {
    if (is.null(locked[[pkg]])) stop("Package missing from lockfile: ", pkg)
    expected <- package_version(contract$directPackages[[pkg]])
    actual <- package_version(locked[[pkg]]$Version)
    if (actual != expected) {
      stop("Lockfile drift for ", pkg, ": expected ", expected, ", got ", actual)
    }
  }
  if (is.null(locked$rsconnect)) stop("Package missing from lockfile: rsconnect")
  if (!is.null(locked$Pmetrics)) stop("Pmetrics must remain outside the public lockfile.")
  cat("Lockfile contract PASS:", length(locked), "packages\n")
}

if (identical(mode, "core")) {
  check_lockfile()
  check_versions(contract$directPackages)
} else if (identical(mode, "pmetrics")) {
  check_versions(c(
    Pmetrics = contract$pmetrics$version,
    mrgsolve = contract$directPackages[["mrgsolve"]],
    mapbayr = contract$directPackages[["mapbayr"]],
    jsonlite = contract$directPackages[["jsonlite"]]
  ))
} else {
  stop("Unknown runtime check mode: ", mode)
}

cat("Runtime contract PASS:", mode, "R", as.character(getRversion()), "\n")
