# Amoxicillin Carlier 2013 — independent linear-system invariants.
# Bibliography: doi:10.1093/jac/dkt240 (Carlier et al., JAC 2013).
#
# This test independently reconstructs the TWO-COMPARTMENT LINEAR IV solution
# using the current catalogue's explicitly audited parameterization.
# It is a SOFTWARE/PHYSICAL-INVARIANT TEST, NOT primary-article source
# reproduction: the complete paper's executable model and study data are not
# independently reproduced here. No clinical use implied.
suppressPackageStartupMessages({
  library(mrgsolve)
  library(jsonlite)
})
args <- commandArgs(FALSE)
script <- sub("^--file=","",grep("^--file=",args,value=TRUE)[[1]])
APP_ROOT <- normalizePath(file.path(dirname(script),".."),mustWork=TRUE)
source(file.path(APP_ROOT,"R","model_library.R"),local=TRUE)

# Closed-form central impulse and rectangular-infusion convolution derived from
# mass-balance ODEs; does not import predictions or model-generated fixtures.
amox_oracle <- function(t,dose,duration,crcl) {
  cl <- 10*(crcl/102)
  vc <- 13.7
  vp <- 13.7
  q <- 15.6
  k10 <- cl/vc
  k12 <- q/vc
  k21 <- q/vp
  d <- sqrt((k10+k12+k21)^2-4*k10*k21)
  alpha <- (k10+k12+k21+d)/2
  beta  <- (k10+k12+k21-d)/2
  stopifnot(alpha>beta,beta>0)
  aa <- (alpha-k21)/(alpha-beta)
  bb <- (k21-beta)/(alpha-beta)
  u <- pmin(t,duration)
  rate <- dose/duration
  concentration <- rate/vc*(
    aa*(-expm1(-alpha*u))/alpha*exp(-alpha*(t-u))+
    bb*(-expm1(-beta*u))/beta*exp(-beta*(t-u))
  )
  list(C=as.numeric(concentration), CL=cl, VC=vc, VP=vp, Q=q)
}

near <- function(actual,pred,label,rtol=3e-5,atol=1e-6) {
  stopifnot(length(actual)==length(pred),all(is.finite(actual)),all(is.finite(pred)))
  if (any(abs(actual-pred)>atol+rtol*abs(pred))) {
    i <- which.max(abs(actual-pred)/(atol+rtol*abs(pred)))
    stop(sprintf("%s: at row %d got %.10g expected %.10g",label,i,actual[[i]],pred[[i]]))
  }
}

model <- mrgsolve::zero_re(compile_model(model_id="amox_carlier"))
route <- model_record("amox_carlier")
cmt <- match(model_administration_cmt(route,"IV"),model@cmtL)
stopifnot(length(cmt)==1L,!is.na(cmt))
test_times <- c(0.25,0.5,1,2,4,8,12,24,48)
auc_times <- seq(0,96,by=0.05)
scenarios <- list(
  list(label="lower renal function, 30 min infusion",crcl=60,duration=0.5,dose=1000),
  list(label="reference renal function, 1 h infusion",crcl=102,duration=1,dose=1000),
  list(label="augmented renal function, 4 h infusion",crcl=160,duration=4,dose=2000)
)
auc_values <- numeric(length(scenarios))
for (i in seq_along(scenarios)) {
  s <- scenarios[[i]]
  m <- mrgsolve::param(model,CRCL=s$crcl)
  event <- mrgsolve::ev(amt=s$dose,cmt=cmt,rate=s$dose/s$duration)
  output <- as.data.frame(mrgsolve::mrgsim(
    m,event,tgrid=test_times,obsonly=TRUE,recsort=3,
    atol=1e-12,rtol=1e-12
  ))
  expected <- amox_oracle(test_times,s$dose,s$duration,s$crcl)
  stopifnot(nrow(output)==length(test_times))
  for (name in c("CL","VC","VP","Q")) {
    near(output[[name]],rep(expected[[name]],length(test_times)),
         paste(s$label,name),rtol=1e-9)
  }
  near(output$DV,expected$C,paste(s$label,"C(t)"))
  stopifnot(all(output$DV>=0))
  dense <- as.data.frame(mrgsolve::mrgsim(
    m,event,tgrid=auc_times,obsonly=TRUE,recsort=3,
    atol=1e-11,rtol=1e-11
  ))
  auc <- sum(diff(dense$time)*(head(dense$DV,-1)+tail(dense$DV,-1))/2)
  target_auc <- s$dose/expected$CL
  if (abs(auc/target_auc-1)>0.003) {
    stop(sprintf("%s: AUC numerical integration %.7g does not satisfy dose/CL=%.7g",
                 s$label,auc,target_auc))
  }
  auc_values[[i]] <- auc
  cat(sprintf("AMOX INVARIANT PASS %s: CL=%.4f L/h, AUC=%.5f, dose/CL=%.5f\n",
              s$label,expected$CL,auc,target_auc))
}
# Dose linearity checked at unchanged renal function and duration.
s <- scenarios[[2]]
single <- as.data.frame(mrgsolve::mrgsim(
  mrgsolve::param(model,CRCL=s$crcl),
  mrgsolve::ev(amt=1000,cmt=cmt,rate=1000),
  tgrid=test_times,obsonly=TRUE,recsort=3))
double <- as.data.frame(mrgsolve::mrgsim(
  mrgsolve::param(model,CRCL=s$crcl),
  mrgsolve::ev(amt=2000,cmt=cmt,rate=2000),
  tgrid=test_times,obsonly=TRUE,recsort=3))
near(double$DV,2*single$DV,"dose-proportional concentration",rtol=2e-5)
cat("Carlier independent two-compartment IV analytic/physical invariants: PASS (3 regimens, dose linearity, mass-exposure law).\n")
