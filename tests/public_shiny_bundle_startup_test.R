# Reproduce a public Shinyapps start from ONLY the files in rsconnect appFiles.
# Existing model tests can pass while the shipped runtime cannot boot.
# No patient inputs, private keys, or network access are used here.
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
Sys.setenv(ALLOW_CUSTOM_MODELS = "false", ENABLE_PMETRICS = "false")
Sys.unsetenv("TDM_GUIDE_PREVIEW_JSON")

# Call the production manifest builder without --deploy / --verify-bundle.
# With a fresh Rscript invocation it returns the exact appFiles selection.
manifest <- source(file.path(root, "scripts", "deploy_shiny.R"),
                   local = TRUE)$value
stopifnot(is.character(manifest), length(manifest) > 0L, "app.R" %in% manifest)
sandbox <- tempfile("mipd-public-appfiles-")
dir.create(sandbox, recursive = TRUE)
on.exit(unlink(sandbox, recursive = TRUE), add = TRUE)
for (relative in manifest) {
  destination <- file.path(sandbox, relative)
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(file.path(root, relative), destination, overwrite = FALSE)) {
    stop("Cannot reproduce production appFiles bundle for: ", relative)
  }
}

# Simulate the immutable release marker exactly as deploy_shiny.R creates it.
sha <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE)
stopifnot(length(sha) == 1L, grepl("^[0-9a-f]{40}$", sha))
release <- list(
  schema = "mipd-public-release.v1", target = "tdmhub/MIPD_Engine",
  appId = "17776160", sourceSha = sha[[1]], modelCount = 46L,
  allowCustomModels = FALSE, enablePmetrics = FALSE
)
jsonlite::write_json(release, file.path(sandbox, "www", "deployment-release.json"),
                     auto_unbox = TRUE, pretty = TRUE)

# A Shiny app must load successfully when the checkout's other files are absent.
old <- setwd(sandbox)
on.exit(setwd(old), add = TRUE)
cat("Starting isolated Shiny production bundle with ", length(manifest),
    " tracked files.\n", sep = "")
stopifnot("R/_disable_autoload.R" %in% manifest)
stopifnot(file.exists(file.path(sandbox, "R", "_disable_autoload.R")))

# This is the real deployment loading path. Unlike sys.source("app.R"),
# shinyAppDir() first calls shiny::loadSupport(), which tries to autoload
# every R/ file in alphabetical order unless _disable_autoload.R is shipped.
# The previous 500 escaped CI because this behavior was not exercised.
old_autoload <- getOption("shiny.autoload.r")
on.exit(options(shiny.autoload.r = old_autoload), add = TRUE)
options(shiny.autoload.r = TRUE)
shiny_bundle <- tryCatch(shiny::shinyAppDir(sandbox), error = function(e) {
  stop("Shinyapps autoload startup failed: ", conditionMessage(e), call. = FALSE)
})
stopifnot(inherits(shiny_bundle, "shiny.appobj"))
cat("Actual Shiny appDir autoload startup PASS\n")

app_env <- new.env(parent = globalenv())
tryCatch(sys.source("app.R", envir = app_env), error = function(e) {
  stop("Public appFiles startup failed: ", conditionMessage(e), call. = FALSE)
})
stopifnot(is.function(app_env$app_ui), is.function(app_env$server))
app <- shiny::shinyApp(ui = app_env$app_ui, server = app_env$server)
stopifnot(inherits(app, "shiny.appobj"))
cat("Public appFiles startup PASS: ", sha[[1]], "\n", sep = "")

# Exercise the actual UI factory in both supported languages, not just data load.
for (language in c("fr", "en")) {
  request <- list(QUERY_STRING = paste0("?lang=", language))
  html <- tryCatch({
    ui <- app_env$app_ui(request)
    htmltools::renderTags(ui)$html
  }, error = function(e) {
    stop("Public ", language, " UI failed: ", conditionMessage(e), call. = FALSE)
  })
  stopifnot(is.character(html), length(html) == 1L, nzchar(html))
  cat("Public ", toupper(language), " UI PASS\n", sep = "")
}
