# An et al., J Antimicrob Chemother 2023;78:1460-1470.
# DOI 10.1093/jac/dkad106, article equation (2), Table 2.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10474939/
#
# INDEPENDENT primary-publication fixed-effect structural cross-check.
# Concentrations are computed from the published TWO-compartment parameters
# by an analytic zero-order-infusion impulse response, not by mrgsolve,
# the engine's regression fixtures or an exported simulation.
# Synthetic subjects only: NO primary-cohort, IIV-distribution, clinical,
# full-original-NONMEM or MAP posterior reproduction claim.
#
# IMPORTANT LIMITATION: the article describes CLCR based on lean body weight
# but does not specify the exact LBW estimator or any LBW floor in its HTML
# Methods. The case-wise Boer estimator is therefore an explicit *local
# implementation assumption*; its fidelity to the article needs original
# supplements/author source. No patient requires an LBW clamp here.

suppressPackageStartupMessages(library(mrgsolve))
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value=TRUE)[[1]])
APP_ROOT <- normalizePath(file.path(dirname(script), ".."), mustWork=TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local=TRUE)

times <- c(0.25, 0.5, 0.75, 1, 2, 3, 4, 8, 12, 24, 36, 48)
check <- function(x,y,label,rtol=3e-5,atol=1e-7) {
  if(length(x)!=length(y) || any(!is.finite(x)) || any(!is.finite(y)))
    stop("Invalid independent cefepime comparison: ",label)
  differences <- abs(as.numeric(x)-as.numeric(y))
  tolerances <- atol + rtol*abs(as.numeric(y))
  if(any(differences>tolerances)) {
    i <- which.max(differences/tolerances)
    stop(sprintf("%s: at %d observed %.10g expected %.10g, delta %.6g > %.6g",
      label,i,x[[i]],y[[i]],differences[[i]],tolerances[[i]]))
  }
}
# Two-compartment IV zero-order input integrated over 0 <= t <= infusion
# duration. k10 = CL/Vc, k12 = Q/Vc, k21 = Q/Vp.
# Impulse coefficients follow algebraically from the 2x2 linear ODE.
infusion_analytic <- function(t,dose,duration,cl,vc,vp,q) {
  stopifnot(duration>0,cl>0,vc>0,vp>0,q>0)
  k10 <- cl/vc
  k12 <- q/vc
  k21 <- q/vp
  discriminant <- sqrt((k10+k12+k21)^2 - 4*k10*k21)
  alpha <- (k10+k12+k21+discriminant)/2
  beta <- (k10+k12+k21-discriminant)/2
  stopifnot(alpha>beta,beta>0)
  a <- (alpha-k21)/(alpha-beta)
  b <- (k21-beta)/(alpha-beta)
  exposure <- pmin(t,duration)
  elapsed_after_stop <- pmax(t-duration,0)
  rate <- dose/duration
  rate/vc*(
    a*(-expm1(-alpha*exposure))/alpha*exp(-alpha*elapsed_after_stop)+
    b*(-expm1(-beta*exposure))/beta*exp(-beta*elapsed_after_stop)
  )
}

# Exact article Table 2 central fixed-effects and reported percentages;
# paper: CL_R=2.00 L/h at CLCR(LBW)=54 mL/min, CL_NR=0.526,
# CL_CRRT=1.64, V1=13.4 L at TBW 90 kg, V2=7.52 L, Q=12 L/h.
# Both combined residual terms are checked in variance units.
reference_model <- compile_model(model_id="cefepime_an")
omega <- diag(as.matrix(mrgsolve::omat(reference_model)))
sigma <- diag(as.matrix(mrgsolve::smat(reference_model)))
check(omega,c(.299^2,.654^2,.173^2),"An 2023 Table 2 IIV",rtol=0,atol=1e-10)
check(sigma,c(.224^2,6.44^2),"An 2023 Table 2 residuals",rtol=0,atol=1e-10)
model <- mrgsolve::zero_re(reference_model)
record <- model_record("cefepime_an")
cmt <- match(model_administration_cmt(record,"IV"),model@cmtL)
stopifnot(length(cmt)==1L,!is.na(cmt))

# Within-range synthetic profiles, changing renal clearance, gender, total
# weight and CRRT status. Actual article patient data are not used.
cases <- list(
 list(label="nonCRRT typical",AGE=65,WT=90,HT=170,CREAT=70,SEX=0,CRRT=0,duration=0.5,dose=2000),
 list(label="nonCRRT renal impairment",AGE=76,WT=72,HT=169,CREAT=170,SEX=0,CRRT=0,duration=3,dose=1000),
 list(label="nonCRRT augmented clearance",AGE=46,WT=112,HT=185,CREAT=55,SEX=0,CRRT=0,duration=0.5,dose=2000),
 list(label="nonCRRT female",AGE=67,WT=74,HT=165,CREAT=86,SEX=1,CRRT=0,duration=3,dose=2000),
 list(label="CRRT male continuous infusion",AGE=65,WT=90,HT=170,CREAT=70,SEX=0,CRRT=1,duration=24,dose=3000),
 list(label="CRRT female short infusion",AGE=70,WT=85,HT=166,CREAT=155,SEX=1,CRRT=1,duration=0.5,dose=1000)
)
for(case in cases) {
  # Case conversion from age/TBW/height/sex/creatinine to CLCR(LBW)
  # is a documented implementation convention (NOT independently
  # established as the exact source's LBW estimator).
  lbw <- if(case$SEX==0) {
    .407*case$WT + .267*case$HT - 19.2
  } else {
    .252*case$WT + .473*case$HT - 48.3
  }
  stopifnot(lbw>30)
  crcl_lbw <- (if(case$SEX==1) .85 else 1)*1.23*(140-case$AGE)*lbw/case$CREAT
  cl <- 2.00*(crcl_lbw/54)*(1-case$CRRT)+.526+1.64*case$CRRT
  vc <- 13.4*(case$WT/90)
  vp <- 7.52
  q <- 12
  oracle <- infusion_analytic(times,case$dose,case$duration,cl,vc,vp,q)
  p <- case[c("AGE","WT","HT","CREAT","SEX","CRRT")]
  configured <- do.call(mrgsolve::param,c(list(model),p))
  result <- as.data.frame(mrgsolve::mrgsim(
    configured,mrgsolve::ev(amt=case$dose,cmt=cmt,rate=case$dose/case$duration),
    tgrid=times,obsonly=TRUE,recsort=3,atol=1e-12,rtol=1e-12
  ))
  stopifnot(nrow(result)==length(times))
  check(result$CL,rep(cl,length(times)),paste(case$label,"published total CL"),rtol=1e-9)
  check(result$VC,rep(vc,length(times)),paste(case$label,"published Vc"),rtol=1e-9)
  check(result$LBW,rep(lbw,length(times)),paste(case$label,"LBW convention"),rtol=1e-9)
  check(result$clcr_calc,rep(crcl_lbw,length(times)),paste(case$label,"CLCR conversion"),rtol=1e-9)
  check(result$on_crrt,rep(case$CRRT,length(times)),paste(case$label,"CRRT switch"),rtol=0)
  check(result$CP,oracle,paste(case$label,"independent two-compartment C(t)"))
  # Independent physical check: the analytical concentration integrates to
  # the IV mass-balance identity AUC(0-infinity) = dose / total clearance.
  # Split at infusion completion to avoid a quadrature kink.
  kernel <- function(x) infusion_analytic(x,case$dose,case$duration,cl,vc,vp,q)
  auc <- stats::integrate(kernel,0,case$duration,rel.tol=1e-8)$value +
    stats::integrate(kernel,case$duration,Inf,rel.tol=1e-8)$value
  check(auc,case$dose/cl,paste(case$label,"independent AUC mass balance"),rtol=1e-5)
  cat("AN 2023 SOURCE-STRUCTURE PASS: ",case$label,"\n",sep="")
}
cat(sprintf("An 2023 source fixed-effects and two-compartment infusion oracle: %d synthetic scenarios PASS (not original full-model reproduction).\n",length(cases)))
