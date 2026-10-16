suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(FALSE)
file <- sub("^--file=", "", grep("^--file=", args, value=TRUE)[[1]])
root <- normalizePath(file.path(dirname(file), ".."), mustWork=TRUE)
source(file.path(root, "R", "deployment_profile.R"), local=TRUE)

tmp <- tempfile("mipd-release-contract-")
dir.create(file.path(tmp, "www"), recursive=TRUE)
on.exit(unlink(tmp, recursive=TRUE), add=TRUE)
release <- list(schema="mipd-public-release.v1",
                target="tdmhub/MIPD_Engine", appId="17776160",
                sourceSha=paste(rep("a",40),collapse=""),
                modelCount=46L,allowCustomModels=FALSE,enablePmetrics=FALSE)
validate_mipd_release(release)
old_custom <- Sys.getenv("ALLOW_CUSTOM_MODELS", unset=NA_character_)
old_pmetrics <- Sys.getenv("ENABLE_PMETRICS", unset=NA_character_)
on.exit({
  if (is.na(old_custom)) Sys.unsetenv("ALLOW_CUSTOM_MODELS") else Sys.setenv(ALLOW_CUSTOM_MODELS=old_custom)
  if (is.na(old_pmetrics)) Sys.unsetenv("ENABLE_PMETRICS") else Sys.setenv(ENABLE_PMETRICS=old_pmetrics)
}, add=TRUE)

Sys.setenv(ALLOW_CUSTOM_MODELS="true",ENABLE_PMETRICS="true")
stopifnot(!enforce_mipd_deployment_profile(tmp))
stopifnot(identical(Sys.getenv("ALLOW_CUSTOM_MODELS"),"true"))
write_json(release,file.path(tmp,"www","deployment-release.json"),auto_unbox=TRUE,pretty=TRUE)
stopifnot(isTRUE(enforce_mipd_deployment_profile(tmp)))
stopifnot(identical(Sys.getenv("ALLOW_CUSTOM_MODELS"),"false"),
          identical(Sys.getenv("ENABLE_PMETRICS"),"false"))

for (broken in list(
  within(release, allowCustomModels <- TRUE),
  within(release, enablePmetrics <- TRUE),
  within(release, sourceSha <- "main"),
  within(release, appId <- "000000")
)) {
  write_json(broken,file.path(tmp,"www","deployment-release.json"),auto_unbox=TRUE)
  err <- tryCatch(enforce_mipd_deployment_profile(tmp),error=identity)
  stopifnot(inherits(err,"error"))
}
cat("MIPD production release fingerprint and host-policy override: PASS.\n")
