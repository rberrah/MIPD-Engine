# Independent MAP objective oracles - Roberts 2011, Revilla 2010 and Goti 2018.
# Sources: doi:10.1128/AAC.01708-10 (Table 2);
#          doi:10.1111/j.1365-2125.2010.03679.x (Table 4);
#          doi:10.1097/FTD.0000000000000490 (Methods and Table 2).
#
# Test data are SYNTHETIC. Reference C(t), random-effect prior and Gaussian
# observation likelihood are written here independently from the model C++
# and from mapbayr. The final comparison uses mapbayr only as the software
# under test. This is NOT reproduction of the authors' original EBEs.
suppressPackageStartupMessages({
  library(mrgsolve)
  library(mapbayr)
  library(jsonlite)
})

args <- commandArgs(FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value=TRUE)[[1]])
APP_ROOT <- normalizePath(file.path(dirname(script), ".."), mustWork=TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local=TRUE)

near <- function(a,b,limit,label) {
  if (length(a)!=length(b) || any(!is.finite(a)) || any(!is.finite(b)) ||
      max(abs(a-b)) > limit) {
    stop(sprintf("%s: max absolute discrepancy %.8g exceeds %.8g",
                 label, max(abs(a-b)), limit))
  }
}

cases <- list(
  list(
    id="vanco_roberts", paper="Roberts 2011 Table 2", dose=1000,
    rate=1000/24, obs_times=c(4,12,20),
    cov=list(CRCL=90.7,WT=74.8),
    cv=c(0.389,0.374), sd_add=2.4, cv_residual=0.199,
    source_omega=c(0.1407,0.1308), source_sigma=c(0.0396,5.76),
    true_eta=c(0.18,-0.11), noise=c(0.20,-0.25,0.33)
  ),
  list(
    id="vanco_pkjust", paper="Revilla 2010 Table 4", dose=1000,
    rate=0, obs_times=c(1,4,12),
    cov=list(AGE=61.1,CRCL=74.7,CREAT=123.8,WT=73),
    cv=c(0.3013,0.2283), sd_add=4.23, cv_residual=0,
    source_omega=c(0.09,0.053), source_sigma=c(0,17.8929),
    true_eta=c(0.16,-0.09), noise=c(-0.20,0.32,-0.28)
  )
)

# Independent analytical 1-compartment predictions, single infusion/bolus.
analytic <- function(case, eta, times) {
  p <- case$cov
  if (case$id == "vanco_roberts") {
    cl <- 4.58*p$CRCL/100*exp(eta[[1]])
    vol <- 1.53*p$WT*exp(eta[[2]])
  } else {
    cl <- (0.67*(p$CRCL/p$WT)+p$AGE^(-0.24))*p$WT*0.06*exp(eta[[1]])
    vol <- 0.82*2.49^(as.integer(p$CREAT > 88.4))*p$WT*exp(eta[[2]])
  }
  stopifnot(cl>0,vol>0)
  if (case$rate > 0) {
    duration <- case$dose/case$rate
    prediction <- ifelse(times<=duration,
      case$rate/cl*(-expm1(-cl*times/vol)),
      case$rate/cl*(-expm1(-cl*duration/vol))*exp(-cl*(times-duration)/vol))
  } else {
    prediction <- case$dose/vol*exp(-cl*times/vol)
  }
  as.numeric(prediction)
}

for (case in cases) {
  model <- compile_model(model_id=case$id)
  omega <- diag(as.matrix(mrgsolve::omat(model)))
  sigma <- diag(as.matrix(mrgsolve::smat(model)))
  stopifnot(length(omega)==2L,length(sigma)==2L)
  near(omega,case$source_omega,1e-8,paste(case$id,"OMEGA against encoded paper values"))
  near(sigma,case$source_sigma,1e-8,paste(case$id,"SIGMA against encoded paper values"))
  near(sigma[[2]],case$sd_add^2,1e-7,paste(case$id,"additive SD to variance"))
  # Published CV to log-normal variance conversion is approximate for these
  # source transcriptions. Report its magnitude rather than claiming equality.
  cv_variance <- log1p(case$cv^2)
  delta <- omega-cv_variance
  stopifnot(max(abs(delta))<0.005)
  cat(sprintf("IIV conversion audit %s: stored-source log-variance differences %.6g / %.6g\n",
              case$id,delta[[1]],delta[[2]]))
  stopifnot(abs(sqrt(sigma[[1]])-case$cv_residual)<0.001)

  observed <- analytic(case,case$true_eta,case$obs_times)+case$noise
  stopifnot(all(observed>0), length(observed)==length(case$obs_times))
  covariates <- case$cov
  model <- do.call(mrgsolve::param,c(list(model),covariates))
  cmt <- match(model_administration_cmt(model_record(case$id),"IV"),model@cmtL)
  stopifnot(length(cmt)==1L,!is.na(cmt))
  n <- length(observed)
  fit_data <- data.frame(
    ID=rep(1L,n+1L),
    time=c(0,case$obs_times),
    evid=c(1L,rep(0L,n)),
    amt=c(case$dose,rep(0,n)),
    cmt=rep(cmt,n+1L),
    rate=c(case$rate,rep(0,n)),
    DV=c(NA_real_,observed),
    mdv=c(1L,rep(0L,n))
  )
  # Independent -2 log posterior up to eta-independent constants:
  # log(sigma_prop*C^2+sigma_add)+(y-C)^2/variance + eta' Omega^-1 eta.
  objective <- function(eta) {
    pred <- analytic(case,eta,case$obs_times)
    var <- sigma[[1]]*pred^2 + sigma[[2]]
    stopifnot(all(var>0))
    sum(log(var)+(observed-pred)^2/var)+sum(eta^2/omega)
  }
  external <- stats::optim(c(0,0),objective,method="BFGS",
                           control=list(reltol=1e-12,maxit=2000))
  stopifnot(external$convergence==0,all(is.finite(external$par)))
  # Check the optimum is locally stable independently of optimizer choice.
  for (start in list(c(-0.25,0.20),c(0.30,-0.20))) {
    rerun <- stats::optim(start,objective,method="BFGS",
                          control=list(reltol=1e-12,maxit=2000))
    stopifnot(rerun$convergence==0)
    near(rerun$par,external$par,0.0005,paste(case$id,"independent optimizer stability"))
  }
  fitted <- mapbayr::mapbayest(x=model,data=fit_data,
                               hessian=stats::optimHess,verbose=FALSE,
                               progress=FALSE)
  eta_map <- as.numeric(mapbayr::get_eta(fitted))
  stopifnot(length(eta_map)==2L,all(is.finite(eta_map)))
  near(eta_map,external$par,0.035,paste(case$id,"MAP ETAs versus independent optimizer"))
  stopifnot(objective(eta_map) <= objective(external$par)+0.025)
  stopifnot(objective(eta_map) < objective(c(0,0)))
  cat(sprintf("Independent synthetic MAP PASS %s: eta=(%.5f,%.5f); reference=(%.5f,%.5f); OFV difference=%.6g\n",
      case$id,eta_map[[1]],eta_map[[2]],external$par[[1]],external$par[[2]],
      objective(eta_map)-objective(external$par)))
}
cat("Independent Roberts/Revilla IIV, SIGMA and two-dimensional MAP test: PASS on synthetic cases.\n")

# Independent three-dimensional Goti oracle. This synthetic case invokes both
# source renal preprocessing rules: creatinine floor for age >60 and CrCl cap.
goti <- list(
  id="vanco_goti", dose=1000, obs_times=c(0.5,2,8,24),
  cov=list(AGE=61,WT=140,CREAT=50,SEX=0,DIAL=0),
  source_omega=log1p(c(0.398,0.816,0.571)^2),
  source_sigma=c(0.227^2,3.4^2),
  true_eta=c(0.12,-0.08,0.10), noise=c(0.16,-0.22,0.11,-0.14)
)

goti_analytic <- function(case,eta,times) {
  p <- case$cov
  creat_used <- if (p$AGE>60 && p$CREAT<88.4) 88.4 else p$CREAT
  sex_factor <- if (p$SEX>=0.5) 1.04 else 1.23
  crcl <- min(150,sex_factor*p$WT*(140-p$AGE)/creat_used)
  dialysis <- as.integer(p$DIAL>=0.5)
  cl <- 4.5*(crcl/120)^0.8*0.7^dialysis*exp(eta[[1]])
  v1 <- 58.4*(p$WT/70)*0.5^dialysis*exp(eta[[2]])
  v2 <- 38.4*exp(eta[[3]])
  q <- 6.5
  k10 <- cl/v1
  k12 <- q/v1
  k21 <- q/v2
  root <- sqrt((k10+k12+k21)^2-4*k10*k21)
  alpha <- ((k10+k12+k21)+root)/2
  beta <- ((k10+k12+k21)-root)/2
  as.numeric(case$dose/v1*(
    (alpha-k21)/(alpha-beta)*exp(-alpha*times)+
    (k21-beta)/(alpha-beta)*exp(-beta*times)
  ))
}

goti_model <- compile_model(model_id=goti$id)
goti_omega <- diag(as.matrix(mrgsolve::omat(goti_model)))
goti_sigma <- diag(as.matrix(mrgsolve::smat(goti_model)))
near(goti_omega,goti$source_omega,1e-9,"vanco_goti OMEGA against paper CV values")
near(goti_sigma,goti$source_sigma,1e-9,"vanco_goti SIGMA against paper values")

goti_observed <- goti_analytic(goti,goti$true_eta,goti$obs_times)+goti$noise
stopifnot(all(goti_observed>0))
goti_model <- do.call(mrgsolve::param,c(list(goti_model),goti$cov))
goti_cmt <- match(model_administration_cmt(model_record(goti$id),"IV"),goti_model@cmtL)
stopifnot(length(goti_cmt)==1L,!is.na(goti_cmt))
goti_n <- length(goti_observed)
goti_data <- data.frame(
  ID=rep(1L,goti_n+1L),
  time=c(0,goti$obs_times),
  evid=c(1L,rep(0L,goti_n)),
  amt=c(goti$dose,rep(0,goti_n)),
  cmt=rep(goti_cmt,goti_n+1L),
  rate=rep(0,goti_n+1L),
  DV=c(NA_real_,goti_observed),
  mdv=c(1L,rep(0L,goti_n))
)
goti_objective <- function(eta) {
  pred <- goti_analytic(goti,eta,goti$obs_times)
  variance <- goti_sigma[[1]]*pred^2+goti_sigma[[2]]
  sum(log(variance)+(goti_observed-pred)^2/variance)+sum(eta^2/goti_omega)
}
goti_external <- stats::optim(c(0,0,0),goti_objective,method="BFGS",
                              control=list(reltol=1e-12,maxit=3000))
stopifnot(goti_external$convergence==0,all(is.finite(goti_external$par)))
for (start in list(c(-0.2,0.2,-0.1),c(0.25,-0.15,0.2))) {
  rerun <- stats::optim(start,goti_objective,method="BFGS",
                        control=list(reltol=1e-12,maxit=3000))
  stopifnot(rerun$convergence==0)
  near(rerun$par,goti_external$par,0.0005,"vanco_goti independent optimizer stability")
}
goti_fit <- mapbayr::mapbayest(x=goti_model,data=goti_data,
                               hessian=stats::optimHess,verbose=FALSE,
                               progress=FALSE)
goti_eta_map <- as.numeric(mapbayr::get_eta(goti_fit))
stopifnot(length(goti_eta_map)==3L,all(is.finite(goti_eta_map)))
near(goti_eta_map,goti_external$par,0.045,"vanco_goti MAP ETAs versus independent optimizer")
stopifnot(goti_objective(goti_eta_map)<=goti_objective(goti_external$par)+0.04)
stopifnot(goti_objective(goti_eta_map)<goti_objective(c(0,0,0)))
cat(sprintf("Independent synthetic MAP PASS vanco_goti: eta=(%.5f,%.5f,%.5f); reference=(%.5f,%.5f,%.5f); OFV difference=%.6g\n",
    goti_eta_map[[1]],goti_eta_map[[2]],goti_eta_map[[3]],
    goti_external$par[[1]],goti_external$par[[2]],goti_external$par[[3]],
    goti_objective(goti_eta_map)-goti_objective(goti_external$par)))
cat("Independent Goti IIV, SIGMA and three-dimensional MAP test: PASS on a synthetic renal-boundary case.\n")
