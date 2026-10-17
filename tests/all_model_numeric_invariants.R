# All 46 models: numerical integrator stability and physical shape checks.
# This test additionally checks a FROZEN 46-MODEL SOFTWARE REGRESSION
# reference derived from the qualified, deployed 2026-10-17 MIPD source SHA.
# Its AUC, sampled peak and 24h concentration are synthetic outputs of the
# model code, NOT independent primary-paper or clinical validation.
suppressPackageStartupMessages({library(mrgsolve); library(jsonlite)})
args <- commandArgs(FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value=TRUE)[[1]])
APP_ROOT <- normalizePath(file.path(dirname(script), ".."), mustWork=TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local=TRUE)

stopifnot(nrow(MODEL_CATALOG) == 46L, !anyDuplicated(MODEL_CATALOG$id))
times <- c(0.25, 0.6, 1.3, 2, 4, 8, 12, 24, 48)
baseline <- jsonlite::fromJSON(
  file.path(APP_ROOT,"tests","fixtures","all_model_software_baselines_v1.json"),
  simplifyVector=FALSE
)
stopifnot(
  identical(baseline$schema,"mipd-46-model-software-reference.v1"),
  identical(baseline$sourceSha,"fec16cbdf8e33ed6c791226ccf236761333f8de8"),
  isTRUE(all.equal(unlist(baseline$timesHours),times)),
  length(baseline$metrics)==49L
)
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
        # Three independent summary dimensions per synthetic dose:
        # AUC0.25-48, sampled peak, and 24h concentration. Two doses per
        # route/mode give six FROZEN values. This catches source-code drift
        # that positive/finite shape invariants alone would not detect.
        group <- paste(id,route,mode,sep="/")
        expected <- unlist(baseline$metrics[[group]],use.names=FALSE)
        if (length(expected)!=6L || any(!is.finite(expected))) {
          stop("Missing or invalid frozen 46-model baseline for ",group)
        }
        index <- c(scenario,scenario+2L,scenario+4L)
        observed <- c(auc,max(a),a[[8L]])
        allowed <- 1e-6+0.002*abs(expected[index])
        if (any(!is.finite(observed)) || any(abs(observed-expected[index])>allowed)) {
          stop(sprintf("FROZEN SOFTWARE REGRESSION DRIFT %s: observed AUC/peak/C24=(%.9g, %.9g, %.9g); frozen=(%.9g, %.9g, %.9g)",
            key,observed[[1]],observed[[2]],observed[[3]],
            expected[index][[1]],expected[index][[2]],expected[index][[3]]))
        }
        cases[[length(cases)+1L]] <- list(
          modelId=id,route=route,administrationMode=mode,scenario=scenario,
          syntheticDose=dose,syntheticInfusionRate=rate,
          concentrations=signif(a,digits=10),
          times=times,auc_0_25_to_48=auc,
          maxSolverNormalizedError=relative_solver_discrepancy,
          status="NUMERIC_INVARIANT_PASS")
      },error=function(e) {
        errors <<- c(errors,paste(key,conditionMessage(e),sep=": "))
      })
    }
  }
}
actual_groups <- unique(vapply(cases,function(entry)
  paste(entry$modelId,entry$route,entry$administrationMode,sep="/"),character(1)))
if (!setequal(actual_groups,names(baseline$metrics)))
  errors <- c(errors,"Frozen regression route/model keys differ from all 46 current models")
report <- list(
  frozenSoftwareBaseline="tests/fixtures/all_model_software_baselines_v1.json",
  frozenSourceSha=baseline$sourceSha,
  frozenRegressionGroups=length(actual_groups),
  frozenSoftwareRegressions=length(cases),
  schema="mipd-all-model-numerical-invariants.v1",
  scope="Two synthetic scenarios/mode: frozen software-baseline AUC/peak/C24 + two solver tolerances; NOT an independently published-source oracle and NOT external clinical validation",
  activeModels=nrow(MODEL_CATALOG),
  testedModels=length(unique(vapply(cases,`[[`,character(1),"modelId"))),
  caseCount=length(cases),failures=errors,cases=cases
)
report_file <- file.path(APP_ROOT,"scientific-validation","all-model-invariants.generated.json")
jsonlite::write_json(report,report_file,pretty=TRUE,auto_unbox=TRUE,digits=10)
if (length(errors) || report$testedModels != 46L || length(cases) < 92L)
  stop("All-model numerical/frozen regression qualification incomplete: ",
       paste(head(errors,20),collapse=" | "),
       "; tested=",report$testedModels,", cases=",length(cases))
cat("All 46 frozen SOFTWARE reference and numerical invariant checks passed across ",length(cases),
    " synthetic dosing scenarios. Additional candidate trajectories stored in ",
    report_file,".\n",sep="")
