`%||%` <- function(value, fallback) if (is.null(value) || !length(value)) fallback else value

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][[1]]
script_path <- sub("^--file=", "", script_arg)
APP_ROOT <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
if (!dir.exists(file.path(APP_ROOT, "R"))) APP_ROOT <- normalizePath("tdm-engine", winslash = "/", mustWork = TRUE)
source(file.path(APP_ROOT, "R", "private_workspace.R"), local = TRUE)

sentinels <- c(
  concentration = "PRIVATE_CONCENTRATION_7f61b3",
  covariate = "PRIVATE_COVARIATE_7f61b3",
  patient = "PRIVATE_PATIENT_7f61b3",
  model = "PRIVATE_MODEL_7f61b3"
)

session_a <- new_private_workspace("pk-mipd-session-a-")
session_b <- new_private_workspace("pk-mipd-session-b-")
on.exit(cleanup_private_workspace(session_a), add = TRUE)
on.exit(cleanup_private_workspace(session_b), add = TRUE)
writeLines(sentinels, file.path(session_a, "private-inputs.txt"), useBytes = TRUE)
stopifnot(!any(vapply(list.files(session_b, full.names = TRUE), function(file) {
  any(grepl(paste(sentinels, collapse = "|"), readLines(file, warn = FALSE)))
}, logical(1))))
stopifnot(isTRUE(cleanup_private_workspace(session_a)), !dir.exists(session_a), dir.exists(session_b))

error_path <- NA_character_
error <- tryCatch(
  with_private_workspace("pk-mipd-error-", function(path) {
    error_path <<- path
    writeLines(sentinels, file.path(path, "interrupted-calculation.txt"), useBytes = TRUE)
    stop("INTERRUPTED_CALCULATION")
  }),
  error = identity
)
stopifnot(inherits(error, "error"), grepl("INTERRUPTED_CALCULATION", conditionMessage(error), fixed = TRUE))
stopifnot(is.character(error_path), !dir.exists(error_path))

stopifnot(isTRUE(cleanup_private_workspace(session_b)), !dir.exists(session_b))
cat("Private workspace isolation and interrupted-calculation cleanup PASS\n")
