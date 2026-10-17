# All 13 experimental MAP priors: numerical covariance and synthetic
# structural-contract checks. This is NOT independent primary-article
# reproduction; the only independent source-equation/MAP oracles currently
# are documented separately for Roberts, Revilla and Goti.
suppressPackageStartupMessages({library(mrgsolve); library(jsonlite)})
arguments <- commandArgs(FALSE)
script <- sub("^--file=", "", grep("^--file=", arguments,value=TRUE)[[1]])
APP_ROOT <- normalizePath(file.path(dirname(script),".."),mustWork=TRUE)
source(file.path(APP_ROOT,"R","model_library.R"),local=TRUE)
records <- MODEL_CATALOG[MODEL_CATALOG$analysisEligible, , drop=FALSE]
stopifnot(nrow(records)==13L)
source_evidence <- jsonlite::fromJSON(
  file.path(APP_ROOT,"scientific-validation","source-evidence.json"),
  simplifyVector=FALSE)$entries
partial <- vapply(source_evidence,`[[`,character(1),"modelId")
stopifnot(setequal(partial,c("vanco_roberts","vanco_pkjust","vanco_goti")))
issues <- character()
rows <- list()
for (id in records$id) tryCatch({
  source_model <- compile_model(model_id=id)
  omega <- as.matrix(mrgsolve::omat(source_model))
  sigma <- as.matrix(mrgsolve::smat(source_model))
  model <- mrgsolve::zero_re(source_model)
  assert_covariance <- function(mat,name) {
    if (length(mat)==0 || nrow(mat)==0 || nrow(mat)!=ncol(mat) ||
        any(!is.finite(mat)) || max(abs(mat-t(mat)))>1e-6)
      stop(name," must be a finite symmetric matrix")
    eig <- eigen((mat+t(mat))/2,symmetric=TRUE,only.values=TRUE)$values
    if (min(eig) < -1e-7) stop(name," has a negative eigenvalue")
    if (name=="OMEGA" && max(eig)<=0) stop("MAP prior must have nonzero variability")
    list(dim=nrow(mat),minEigenvalue=min(eig))
  }
  om <- assert_covariance(omega,"OMEGA")
  sm <- assert_covariance(sigma,"SIGMA")
  routes <- model_routes(model_record(id))
  if (length(routes)<1L) stop("No specified route")
  route <- routes[[1]]
  cmt <- match(model_administration_cmt(model_record(id),route),model@cmtL)
  if (is.na(cmt)) stop("Unrecognized administration compartment")
  dose <- if (model_record(id)$concentrationUnit[[1]]=="ng/mL") 1 else 100
  out <- as.data.frame(mrgsolve::mrgsim(
    model,mrgsolve::ev(amt=dose,cmt=cmt,rate=if(route=="IV") dose else 0),
    tgrid=c(1.4,3,8,16,24),obsonly=TRUE,recsort=3,
    atol=1e-10,rtol=1e-10))
  exposure <- as.numeric(out[[pick_concentration_column(out)]])
  if (length(exposure)!=5L || any(!is.finite(exposure)) ||
      !any(exposure>0)) stop("MAP prior has no finite positive synthetic profile")
  rows[[length(rows)+1L]] <- list(
    modelId=id,sourceEquationEvidence=if(id %in% partial)
      "PARTIAL_SOURCE_EQUATIONS" else "NOT_INDEPENDENTLY_REPRODUCED",
    omegaDimension=om$dim,omegaMinEigenvalue=om$minEigenvalue,
    sigmaDimension=sm$dim,sigmaMinEigenvalue=sm$minEigenvalue,
    syntheticConcentrations=signif(exposure,digits=10),
    softwareContract="PASS",externalClinicalValidation=FALSE)
},error=function(e){issues<<-c(issues,paste(id,conditionMessage(e),sep=": "))})
report <- list(schema="mipd-map-prior-numeric-matrix.v1",
  scope="13 MAP-eligible models: positive-semidefinite OMEGA/SIGMA and finite synthetic profiles; not an independent primary-source/MAP optimizer reproduction",
  expectedPriors=13,examined=length(rows),partialIndependentSourceModels=partial,
  completeSourceReproductions=0,externalClinicalValidations=0,
  failures=issues,models=rows)
path <- file.path(APP_ROOT,"scientific-validation","map-prior-matrix.generated.json")
jsonlite::write_json(report,path,pretty=TRUE,auto_unbox=TRUE,digits=10)
if(length(issues)||length(rows)!=13L)
  stop("MAP prior structural checks incomplete: ",paste(issues,collapse=" | "))
cat("MAP prior numerical matrix PASS: 13/13 software priors, ",
    length(partial)," partial independent source models, 0 full source reproductions.\n",sep="")
