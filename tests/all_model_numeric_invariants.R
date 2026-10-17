# All 46 models: numerical integrator stability and physical shape checks.
# This is NOT source-paper reproduction and is NOT a frozen numerical
# reference profile. The output is an auditable CANDIDATE fixture report to
# review independently before any reference value is promoted.
suppressPackageStartupMessages({library(mrgsolve); library(jsonlite)})
args <- commandArgs(FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value=TRUE)[[1]])
APP_ROOT <- normalizePath(file.path(dirname(script), ".."), mustWork=TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local=TRUE)

stopifnot(nrow(MODEL_CATALOG) == 46L, !anyDuplicated(MODEL_CATALOG$id))
times <- c(0.25, 0.6, 1.3, 2, 4, 8, 12, 24, 48)
cases <- list()
errors <- character()
for (id in MODEL_CATALOG$id) {
  record <- model_record(id)
  model <- mrgsolve::zero_re(compile_model(model_id=id))
  for (route in model_routes(record)) {
    modes <- model_administration_modes(record,route)
    for (mode in modes) for (scenario in 1:2) {
      dose <- if (record$concentrationUnit[[1]] == "ng/mL") 1 else 100
      dose <- dose * if (scenario == 1L) 1 else 1.7
      duration <- if (mode == "IV_CONTINUOUS") 24 else 1
      rate <- if (route == "IV") dose / duration else 0
      cmt <- match(model_administration_cmt(record,route), model@cmtL)
      key <- paste(id,route,mode,scenario,sep="/")
      tryCatch({
        stopifnot(length(cmt)==1L,!is.na(cmt),is.finite(dose),dose>0)
        event <- mrgsolve::ev(amt=dose,cmt=cmt,rate=rate)
        simulate <- function(tol) as.data.frame(mrgsolve::mrgsim(
          model,event,tgrid=times,obsonly=TRUE,recsort=3,
          atol=tol,rtol=tol
        ))
        high <- simulate(1e-11)
        ordinary <- simulate(1e-8)
        stopifnot(nrow(high)==length(times),nrow(ordinary)==length(times))
        col <- pick_concentration_column(high)
        a <- as.numeric(high[[col]])
        b <- as.numeric(ordinary[[col]])
        stopifnot(length(a)==length(times),length(b)==length(times),
                  all(is.finite(a)),all(is.finite(b)))
        shape_tolerance <- 1e-6 + 1e-4*max(abs(a))
        if (min(a) < -shape_tolerance || max(a) <= 0)
          stop("negative concentrations or absent positive exposure")
        relative_solver_discrepancy <- max(abs(a-b)/(1e-6+0.015*pmax(abs(a),abs(b))))
        if (!is.finite(relative_solver_discrepancy) || relative_solver_discrepancy > 1)
          stop("solver discrepancy exceeds 1.5% relative + 1e-6 absolute")
        auc <- sum(diff(times)*(head(a,-1)+tail(a,-1))/2)
        stopifnot(is.finite(auc),auc>0)
        cases[[length(cases)+1L]] <- list(
          modelId=id,route=route,administrationMode=mode,scenario=scenario,
          syntheticDose=dose,syntheticInfusionRate=rate,
          concentrations=round(a,signif=10),
          times=times,auc_0_25_to_48=auc,
          maxSolverNormalizedError=relative_solver_discrepancy,
          status="NUMERIC_INVARIANT_PASS")
      },error=function(e) {
        errors <<- c(errors,paste(key,conditionMessage(e),sep=": "))
      })
    }
  }
}
report <- list(
  schema="mipd-all-model-numerical-invariants.v1",
  scope="Two synthetic dosing scenarios per administration mode, two solver tolerances; NOT source reproduction, clinical validation or a reviewed frozen oracle",
  activeModels=nrow(MODEL_CATALOG),
  testedModels=length(unique(vapply(cases,`[[`,character(1),"modelId"))),
  caseCount=length(cases),failures=errors,cases=cases
)
report_file <- file.path(APP_ROOT,"scientific-validation","all-model-invariants.generated.json")
jsonlite::write_json(report,report_file,pretty=TRUE,auto_unbox=TRUE,digits=10)
if (length(errors) || report$testedModels != 46L || length(cases) < 92L)
  stop("All-model numerical qualification incomplete: ",
       paste(head(errors,20),collapse=" | "),
       "; tested=",report$testedModels,", cases=",length(cases))
cat("All 46 model numerical invariant checks passed across ",length(cases),
    " synthetic dosing scenarios. Unreviewed candidate trajectories stored in ",
    report_file,".\n",sep="")
