suppressPackageStartupMessages(library(mrgsolve))
args <- commandArgs(trailingOnly = FALSE)
file <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(file), ".."), winslash = "/", mustWork = TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
source(file.path(APP_ROOT, "R", "model_library.R"))
fixtures <- file.path(APP_ROOT, "..", "test-results", "advanced-native")
files <- list.files(fixtures, pattern = "\\.cpp$", full.names = TRUE)
stopifnot(length(files) >= 9)
cache <- tempfile("advanced-lego-")
dir.create(cache)
for (file in files) {
  code <- paste(readLines(file, warn = FALSE), collapse = "\n")
  spec <- lego_spec_from_code(code)
  name <- gsub("[^A-Za-z0-9_]", "_", tools::file_path_sans_ext(basename(file)))
  native <- zero_re(mcode(paste0("raw_", name), code, soloc = cache, quiet = TRUE))
  safe <- zero_re(mcode(paste0("safe_", name), lego_model_code(spec), soloc = cache, quiet = TRUE))
  raw <- as.data.frame(mrgsim(native, events = ev(amt = 100, cmt = 1), end = 168, delta = 0.21, recsort = 3))
  regenerated <- as.data.frame(mrgsim(safe, events = ev(amt = 100, cmt = 1), end = 168, delta = 0.21, recsort = 3))
  for (node in spec$nodes) {
    a <- raw[[node$name]]
    b <- regenerated[[paste0("L", node$id, "_", node$name)]]
    stopifnot(length(a) == length(b), all(is.finite(a)), all(is.finite(b)), max(abs(a-b)) < 1e-5)
    if (node$kind == "interaction") stopifnot(max(abs(raw[[paste0("MOD_",node$name)]]-regenerated[[paste0("MOD_L",node$id,"_",node$name)]])) < 1e-5)
  }
  if (name == "composed") {
    preview <- jsonlite::fromJSON(file.path(fixtures, "composed-preview.json"))
    preview <- preview[!duplicated(preview$id), ]
    for (i in seq_len(nrow(preview))) {
      node <- Filter(function(n) n$id == preview$id[i], spec$nodes)[[1]]
      column <- if (node$kind == "interaction") paste0("MOD_",node$name) else node$name
      expected <- tail(raw[[column]], 1) / (node$vol %||% 1)
      stopifnot(abs(expected-preview$value[i]) < 1e-4)
    }
  }
  tumor <- Filter(function(n) n$kind == "tumor", spec$nodes)
  if (length(tumor)) {
    node <- tumor[[1]]
    zero_kill <- param(safe, setNames(list(0), paste0("TV_kill_L", node$id, "_", node$name)))
    out <- as.data.frame(mrgsim(zero_kill, end = 168, delta = 1))
    expected <- if (node$growth == "gompertz") node$cap*exp(log(node$t0/node$cap)*exp(-node$kg*out$time)) else if (node$growth == "logistic") node$cap/(1+(node$cap/node$t0-1)*exp(-node$kg*out$time)) else node$t0*exp(node$kg*out$time)
    stopifnot(max(abs(out[[paste0("L",node$id,"_",node$name)]]-expected)) < 1e-4)
  }
  cat(basename(file), ": native and safe trajectories identical\n")
}
bad <- spec
bad$nodes[[length(bad$nodes)+1]] <- list(id=99,kind="interaction",name="bad",dose=0,source=1,mechanism="factor",factor=1,targetFrom=1,targetTo=999)
error <- tryCatch({normalize_lego_spec(bad);""},error=conditionMessage)
stopifnot(grepl("target flux",error))
bad$nodes[[length(bad$nodes)]]$mechanism <- "system_call"
stopifnot(inherits(try(normalize_lego_spec(bad),silent=TRUE),"try-error"))
cat("Advanced safe compilation, analytic no-treatment limit and invalid input rejection OK.\n")
