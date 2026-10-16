$PROB
# Vancomycin PopPK -- Roberts et al. (2011)
# Antimicrob Agents Chemother. 2011;55:2704-2709.
# DOI: 10.1128/AAC.01708-10
#
# One-compartment model developed in 206 septic adult ICU patients receiving
# continuous infusion. The published renal covariate is measured 24-hour urinary
# creatinine clearance normalized to 1.73 m2; it is therefore entered directly as
# CRCL and is not replaced by Cockcroft-Gault.

$PARAM @annotated
TVCL : 4.58 : Clearance at CRCL 100 mL/min/1.73 m2 (L/h)
TVV  : 1.53 : Volume of distribution per kg total body weight (L/kg)
ETA1 : 0    : Individual ETA on clearance (mapbayr)
ETA2 : 0    : Individual ETA on volume (mapbayr)

$PARAM @annotated @covariates
CRCL : 90.7 : Measured 24-hour urinary creatinine clearance (mL/min/1.73 m2)
WT   : 74.8 : Total body weight (kg)

$OMEGA @annotated @diagonal
ETA_CL : 0.1407 : IIV on CL, from 38.9 %CV
ETA_V  : 0.1308 : IIV on V, from 37.4 %CV

$SIGMA @annotated @diagonal
PROP : 0.0396 : Proportional residual variance, from 19.9 %CV
ADD  : 5.76   : Additive residual variance, from SD 2.4 mg/L

$CMT @annotated
CENT : Central compartment amount (mg) [ADM, OBS]

$MAIN
double safe_CRCL = CRCL;
if(safe_CRCL < 0.0) safe_CRCL = 0.0;

double safe_WT = WT;
if(safe_WT < 1.0) safe_WT = 1.0;

double CL = TVCL * (safe_CRCL / 100.0) * exp(ETA1 + ETA(1));
double V  = TVV * safe_WT * exp(ETA2 + ETA(2));

$ODE
dxdt_CENT = -(CL / V) * CENT;

$TABLE
double IPRED = CENT / V;
double DV = IPRED * (1.0 + EPS(1)) + EPS(2);

int i = 0;
while(DV < 0.0 && i < 100) {
  simeps();
  DV = IPRED * (1.0 + EPS(1)) + EPS(2);
  ++i;
}

$CAPTURE @annotated
DV    : Simulated concentration with residual error (mg/L)
CL    : Individual clearance (L/h)
V     : Individual volume of distribution (L)
