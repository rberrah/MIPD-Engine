new_private_workspace <- function(prefix = "pk-mipd-session-") {
  path <- tempfile(prefix)
  if (!dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create the private session workspace.")
  }
  path
}

cleanup_private_workspace <- function(path, cache = NULL) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) return(invisible(FALSE))
  session_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  loaded <- getLoadedDLLs()
  dll_paths <- vapply(loaded, function(dll) dll[["path"]] %||% "", character(1))
  dll_paths <- dll_paths[nzchar(dll_paths)]
  normalized <- normalizePath(dll_paths, winslash = "/", mustWork = FALSE)
  for (dll_path in unique(dll_paths[startsWith(normalized, session_path)])) {
    try(dyn.unload(dll_path), silent = TRUE)
  }
  if (is.environment(cache)) rm(list = ls(cache, all.names = TRUE), envir = cache)
  invisible(gc())
  unlink(path, recursive = TRUE, force = TRUE)
  invisible(!dir.exists(path))
}

with_private_workspace <- function(prefix, action) {
  path <- new_private_workspace(prefix)
  on.exit(cleanup_private_workspace(path), add = TRUE)
  action(path)
}
