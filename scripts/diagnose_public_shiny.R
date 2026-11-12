# Fetch ONLY technical startup errors from the existing public Shiny application.
# Public GitHub Actions logs must never contain unfiltered server logs, secret
# values, request payloads, concentrations, or other patient-shaped content.
suppressPackageStartupMessages(library(rsconnect))

account <- Sys.getenv("SHINYAPPS_ACCOUNT", "")
token <- Sys.getenv("SHINYAPPS_TOKEN", "")
secret <- Sys.getenv("SHINYAPPS_SECRET", "")
if (!identical(account, "tdmhub") || !nzchar(token) || !nzchar(secret)) {
  stop("Shinyapps production environment secrets not configured")
}
rsconnect::setAccountInfo(name = account, token = token, secret = secret)

logs <- tryCatch(rsconnect::getLogs(
  appName = "MIPD_Engine", account = account,
  server = "shinyapps.io", entries = 250
), error = function(e) {
  stop("Could not retrieve Shinyapps logs (no application data printed): ",
       class(e)[[1]], call. = FALSE)
})
if (is.null(logs) || !nrow(logs) || !"message" %in% names(logs)) {
  stop("Shinyapps returned no structured log lines")
}
messages <- as.character(logs$message)
# Filter out ordinary session traffic. Only pre-session R startup/errors matter.
has_error <- grepl(
  "Error in |Error:|ERROR|Execution halted|failed to|cannot open|not found|not available|invalid|fatal|terminat|exited",
  messages, ignore.case = TRUE
)
# R/shiny emits its actual exception *after* the generic
# "Error in ..stacktraceon.." prefix. Retain one following line as context.
# Never print the raw lines: everything is sanitized below.
error_indices <- tail(which(has_error), 15L)
selected <- sort(unique(unlist(lapply(error_indices, function(i) {
  seq.int(max(1L, i), min(length(messages), i + 1L))
}), use.names = FALSE)))
errors <- tail(messages[selected], 30L)
if (!length(errors)) {
  cat("SHINY_STARTUP_DIAGNOSTIC: no error-pattern lines in recent server logs; ",
      nrow(logs), " entries inspected.\n", sep = "")
} else {
  sanitize <- function(x) {
    # Hide identifiers, file path details, all quoted strings and numeric
    # values; retain only standardized R error categories and operations.
    x <- gsub("https?://[^[:space:]]+", "<url>", x, perl = TRUE)
    x <- gsub("[[:alnum:]._%+-]+@[[:alnum:].-]+", "<email>", x, perl = TRUE)
    x <- gsub('"[^"]*"', '"<quoted>"', x, perl = TRUE)
    x <- gsub("'[^']*'", "'<quoted>'", x, perl = TRUE)
    x <- gsub("(?i)(token|secret|password)[[:space:]]*[:=][^[:space:]]+",
              "<credential>", x, perl = TRUE)
    x <- gsub("[0-9a-fA-F]{20,}", "<hex>", x, perl = TRUE)
    x <- gsub("(?<![[:alnum:]_])/[[:graph:]]+", "<path>", x, perl = TRUE)
    x <- gsub("[[:digit:]]+", "#", x, perl = TRUE)
    substr(trimws(x), 1L, 190L)
  }
  cat("SHINY_STARTUP_DIAGNOSTIC: ", length(errors),
      " filtered technical errors from recent server logs:\n", sep = "")
  cat(paste(sanitize(errors), collapse = "\n"), "\n")
}
