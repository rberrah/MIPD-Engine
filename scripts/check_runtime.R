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

if (identical(mode, "core")) {
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
