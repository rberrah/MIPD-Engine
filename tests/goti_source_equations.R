# Goti 2018 independent fixed-effect IV two-compartment vancomycin checks.
# Source: Goti et al. Ther Drug Monit 2018;40:212-221.
# DOI: 10.1097/FTD.0000000000000490
# Published final effects independently cross-checked in
# Antibiotics 2023;12:301, Table 2 (doi:10.3390/antibiotics12020301).
#
# The oracle derives the analytical IV bolus impulse response from the
# published two-compartment rate matrix, without using the C++ model
# parameters, imported trajectories, or historical regression fixtures.
# Coverage is limited to covariates not invoking ambiguous upper bounds.
# OMEGA, SIGMA, MAP, dialysis over time, clinical validation: NOT tested.

suppressPackageStartupMessages({
  library(mrgsolve)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1]])
APP_ROOT <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local = TRUE)

times <- c(0.25, 1, 3, 8, 12, 24)
dose <- 1000  # synthetic IV bolus in mg

# Independent typical-value equations, expressed with units.
goti_oracle <- function(age, wt, creat_umol, female, dial) {
  stopifnot(age > 0, age <= 65 || creat_umol >= 88.4, wt > 0,
            creat_umol > 0, female %in% c(0, 1), dial %in% c(0, 1))
  sex_factor <- if (female == 1) 1.04 else 1.25
  crcl <- sex_factor * wt * (140 - age) / creat_umol
  stopifnot(crcl > 0, crcl < 150) # avoids disputed source renal caps
  cl <- 4.5 * (crcl / 120)^0.8 * (0.7^dial)    # L/h
  v1 <- 58.4 * (wt / 70) * (0.5^dial)         # L
  v2 <- 38.4                                  # L
  q <- 6.5                                    # L/h
  k10 <- cl / v1
  k12 <- q / v1
  k21 <- q / v2
  root <- sqrt((k10 + k12 + k21)^2 - 4*k10*k21)
  alpha <- ((k10 + k12 + k21) + root) / 2
  beta <- ((k10 + k12 + k21) - root) / 2
  stopifnot(alpha > beta, beta > 0)
  # Independent impulse response for central concentration after IV bolus:
  # C(t)=D/V1*[(alpha-k21)/(alpha-beta)*e^(-alpha*t)
  #            +(k21-beta)/(alpha-beta)*e^(-beta*t)].
  predicted <- dose / v1 * (
    ((alpha - k21) / (alpha - beta)) * exp(-alpha * times) +
    ((k21 - beta) / (alpha - beta)) * exp(-beta * times)
  )
  list(CRCL=crcl, CL=cl, V1=v1, V2=v2, Q=q, C=predicted)
}

check_close <- function(actual, expected, label, rtol=2e-5, atol=1e-7) {
  stopifnot(length(actual) == length(expected), all(is.finite(actual)),
            all(is.finite(expected)))
  error <- abs(actual - expected)
  if (any(error > atol + rtol * abs(expected))) {
    idx <- which.max(error / (atol + rtol * abs(expected)))
    stop(sprintf("%s mismatch at %d: got %.9g, expected %.9g (error %.4g)",
                 label, idx, actual[[idx]], expected[[idx]], error[[idx]]))
  }
}

cases <- list(
  list(label="male without dialysis", AGE=50, WT=70, CREAT=80, SEX=0, DIAL=0),
  list(label="female with dialysis", AGE=55, WT=85, CREAT=125, SEX=1, DIAL=1),
  list(label="older male adequate serum creatinine", AGE=70, WT=75, CREAT=110, SEX=0, DIAL=0)
)

record <- model_record("vanco_goti")
model <- mrgsolve::zero_re(compile_model(model_id="vanco_goti"))
cmt <- match(model_administration_cmt(record, "IV"), model@cmtL)
stopifnot(length(cmt)==1L, !is.na(cmt))
for (case in cases) {
  expected <- with(case, goti_oracle(AGE, WT, CREAT, SEX, DIAL))
  args <- case[c("AGE","WT","CREAT","SEX","DIAL")]
  compiled <- do.call(mrgsolve::param, c(list(model), args))
  observed <- as.data.frame(mrgsolve::mrgsim(
    compiled, mrgsolve::ev(amt=dose, cmt=cmt, rate=0),
    tgrid=times, obsonly=TRUE, recsort=3, atol=1e-12, rtol=1e-12
  ))
  stopifnot(nrow(observed)==length(times))
  for (name in c("CRCL","CL","V1","V2","Q")) {
    check_close(observed[[name]], rep(expected[[name]], length(times)),
                paste(case$label, name), rtol=1e-9)
  }
  check_close(observed$DV, expected$C, paste(case$label, "C(t)"))
  cat("GOTI ANALYTICAL PASS:", case$label, "\n")
}

# The secondary source documents renal boundary handling not yet accounted
# for in the MIPD implementation. Audit these branches without silently
# changing the frozen model, or calling the source model fully reproduced.
published_crcl <- function(age, wt, creat_umol, female) {
  scr <- if (age > 65) max(creat_umol, 88.4) else creat_umol
  min(150, (if (female == 1) 1.04 else 1.25) * wt * (140-age) / scr)
}
boundary_cases <- list(
  list(label="augmented renal clearance upper bound", AGE=40, WT=100,
       CREAT=50, SEX=0, DIAL=0),
  list(label="older person with low creatinine", AGE=75, WT=70,
       CREAT=50, SEX=0, DIAL=0)
)
for (case in boundary_cases) {
  compiled <- do.call(mrgsolve::param,
                      c(list(model), case[c("AGE","WT","CREAT","SEX","DIAL")]))
  obs <- as.data.frame(mrgsolve::mrgsim(
    compiled, mrgsolve::ev(amt=dose, cmt=cmt, rate=0),
    tgrid=c(1), obsonly=TRUE, recsort=3
  ))
  expected_crcl <- with(case, published_crcl(AGE, WT, CREAT, SEX))
  actual_crcl <- obs$CRCL[[1]]
  # Explicit nonfatal SCIENCE GAP, not a falsely passing source oracle.
  cat(sprintf("GOTI BOUNDARY SCIENCE GAP %s: model CrCl=%.4f, published secondary rule CrCl=%.4f mL/min\n",
              case$label, actual_crcl, expected_crcl))
  stopifnot(is.finite(actual_crcl), is.finite(expected_crcl))
}
cat("Goti 2018: 3 in-domain independent fixed-effect cases passed; 2 boundary discrepancies require primary-source adjudication.\n")
