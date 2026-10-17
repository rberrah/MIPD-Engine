# Goti 2018 independent fixed-effect IV two-compartment vancomycin checks.
# Primary source: Goti et al. Ther Drug Monit. 2018;40:212-221.
# DOI: 10.1097/FTD.0000000000000490
# Source locations: Methods (renal covariate preprocessing) and Table 2.
#
# The source oracle is independent from the C++ implementation. It covers the
# published creatinine floor for age >60 years, CrCl cap, typical parameters,
# random-effect and residual variances, concentration-time response, and AUC.
# Test patients and observations are synthetic. This is not external validation.

suppressPackageStartupMessages({
  library(mrgsolve)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[[1]])
APP_ROOT <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
source(file.path(APP_ROOT, "R", "model_library.R"), local = TRUE)

times <- c(0.25, 1, 3, 8, 12, 24)
dose <- 1000

check_close <- function(actual, expected, label, rtol = 2e-5, atol = 1e-7) {
  stopifnot(length(actual) == length(expected), all(is.finite(actual)),
            all(is.finite(expected)))
  error <- abs(actual - expected)
  if (any(error > atol + rtol * abs(expected))) {
    idx <- which.max(error / (atol + rtol * abs(expected)))
    stop(sprintf("%s mismatch at %d: got %.9g, expected %.9g (error %.4g)",
                 label, idx, actual[[idx]], expected[[idx]], error[[idx]]))
  }
}

goti_oracle <- function(age, wt, creat_umol, female, dial) {
  stopifnot(age > 0, age < 140, wt > 0, creat_umol > 0,
            female %in% c(0, 1), dial %in% c(0, 1))
  creat_used <- if (age > 60 && creat_umol < 88.4) 88.4 else creat_umol
  sex_factor <- if (female == 1) 1.04 else 1.23
  crcl_raw <- sex_factor * wt * (140 - age) / creat_used
  crcl <- min(crcl_raw, 150)
  cl <- 4.5 * (crcl / 120)^0.8 * 0.7^dial
  v1 <- 58.4 * (wt / 70) * 0.5^dial
  v2 <- 38.4
  q <- 6.5
  k10 <- cl / v1
  k12 <- q / v1
  k21 <- q / v2
  root <- sqrt((k10 + k12 + k21)^2 - 4 * k10 * k21)
  alpha <- ((k10 + k12 + k21) + root) / 2
  beta <- ((k10 + k12 + k21) - root) / 2
  coef_alpha <- (alpha - k21) / (alpha - beta)
  coef_beta <- (k21 - beta) / (alpha - beta)
  concentration <- dose / v1 * (
    coef_alpha * exp(-alpha * times) + coef_beta * exp(-beta * times)
  )
  auc_coefficients <- dose / v1 * (coef_alpha / alpha + coef_beta / beta)
  list(
    CREAT_USED = creat_used, CRCL_RAW = crcl_raw, CRCL = crcl,
    CL = cl, V1 = v1, V2 = v2, Q = q, C = concentration,
    AUC = dose / cl, AUC_COEFFICIENTS = auc_coefficients
  )
}

cases <- list(
  list(label = "male without dialysis", AGE = 50, WT = 70, CREAT = 80, SEX = 0, DIAL = 0),
  list(label = "female with dialysis", AGE = 55, WT = 85, CREAT = 125, SEX = 1, DIAL = 1),
  list(label = "older male above creatinine floor", AGE = 70, WT = 75, CREAT = 110, SEX = 0, DIAL = 0),
  list(label = "age 60 does not invoke floor", AGE = 60, WT = 50, CREAT = 80, SEX = 0, DIAL = 0),
  list(label = "age 61 invokes floor", AGE = 61, WT = 50, CREAT = 80, SEX = 0, DIAL = 0),
  list(label = "creatinine exactly 1 mg/dL", AGE = 70, WT = 70, CREAT = 88.4, SEX = 0, DIAL = 0),
  list(label = "augmented renal clearance is capped", AGE = 40, WT = 100, CREAT = 50, SEX = 0, DIAL = 0),
  list(label = "raw clearance exactly cap", AGE = 40, WT = 100, CREAT = 82, SEX = 0, DIAL = 0)
)

source_model <- compile_model(model_id = "vanco_goti")
omega <- diag(as.matrix(mrgsolve::omat(source_model)))
sigma <- diag(as.matrix(mrgsolve::smat(source_model)))
check_close(omega, log1p(c(0.398, 0.816, 0.571)^2), "published IIV CV conversion", rtol = 0, atol = 5e-10)
check_close(sigma, c(0.227^2, 3.4^2), "published residual variance conversion", rtol = 0, atol = 1e-10)
model <- mrgsolve::zero_re(source_model)

record <- model_record("vanco_goti")
cmt <- match(model_administration_cmt(record, "IV"), model@cmtL)
stopifnot(length(cmt) == 1L, !is.na(cmt))

for (case in cases) {
  expected <- with(case, goti_oracle(AGE, WT, CREAT, SEX, DIAL))
  parameters <- case[c("AGE", "WT", "CREAT", "SEX", "DIAL")]
  compiled <- do.call(mrgsolve::param, c(list(model), parameters))
  observed <- as.data.frame(mrgsolve::mrgsim(
    compiled, mrgsolve::ev(amt = dose, cmt = cmt, rate = 0),
    tgrid = times, obsonly = TRUE, recsort = 3, atol = 1e-12, rtol = 1e-12
  ))
  stopifnot(nrow(observed) == length(times))
  for (name in c("CREAT_USED", "CRCL_RAW", "CRCL", "CL", "V1", "V2", "Q")) {
    check_close(observed[[name]], rep(expected[[name]], length(times)),
                paste(case$label, name), rtol = 1e-9)
  }
  check_close(observed$CP, expected$C, paste(case$label, "C(t)"))
  check_close(expected$AUC_COEFFICIENTS, expected$AUC,
              paste(case$label, "analytical AUC identity"), rtol = 1e-10)
  check_close(dose / observed$CL[[1]], expected$AUC,
              paste(case$label, "implementation AUC"), rtol = 1e-9)
  cat("GOTI SOURCE PASS:", case$label, "\n")
}

cat(sprintf("Goti 2018 source-equation test: PASS (%d synthetic cases, including renal boundaries).\n",
            length(cases)))
