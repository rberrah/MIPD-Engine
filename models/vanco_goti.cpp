$PROB
// Vancomycin population PK model (Goti et al., 2018).
// Population: hospitalized adults with and without intermittent hemodialysis.
// DIAL is the study hemodialysis-status covariate, not a time-varying dialysis clearance.
// Reference: doi:10.1097/FTD.0000000000000490.

$PARAM @annotated
TVCL     : 4.5   : Typical clearance at CRCL 120 mL/min without dialysis (L/h)
TVV1     : 58.4  : Typical central volume at 70 kg without dialysis (L)
TVV2     : 38.4  : Typical peripheral volume (L)
TVQ      : 6.5   : Typical intercompartmental clearance (L/h)
CRCL_CL  : 0.8   : Cockcroft-Gault clearance exponent on CL
DIAL_CL  : 0.7   : Hemodialysis-status factor on CL
DIAL_V1  : 0.5   : Hemodialysis-status factor on V1
REF_CRCL : 120.0 : Reference creatinine clearance (mL/min)
REF_WT   : 70.0  : Reference body weight (kg)
ETA1 : 0.0 : Posterior ETA offset on CL
ETA2 : 0.0 : Posterior ETA offset on V1
ETA3 : 0.0 : Posterior ETA offset on V2

$PARAM @covariates @annotated
AGE   : 50.0 : Age (years)
WT    : 70.0 : Body weight (kg)
CREAT : 80.0 : Serum creatinine (micromol/L)
SEX   : 0.0  : Sex (0 = male, 1 = female)
DIAL  : 0.0  : Hemodialysis status (0 = no, 1 = yes)

$CMT @annotated
CENT   : Central compartment (mg) [ADM, OBS]
PERIPH : Peripheral compartment (mg)

$OMEGA @annotated @diagonal
ETA_CL : 0.158 : Between-subject variance on CL
ETA_V1 : 0.510 : Between-subject variance on V1
ETA_V2 : 0.280 : Between-subject variance on V2

$SIGMA @annotated @diagonal
PROP : 0.051529 : Proportional residual variance (CV 22.7%)
ADD  : 11.56    : Additive residual variance (SD 3.4 mg/L)

$MAIN
double age_safe = (AGE > 0.0 && AGE < 140.0) ? AGE : 50.0;
double wt_safe = (WT > 0.0) ? WT : REF_WT;
double creat_safe = (CREAT > 0.0) ? CREAT : 80.0;
double sex_coefficient = (SEX >= 0.5) ? 1.04 : 1.25;
double CRCL = sex_coefficient * wt_safe * (140.0 - age_safe) / creat_safe;
if (CRCL < 0.1) CRCL = 0.1;
double dialysis = (DIAL >= 0.5) ? 1.0 : 0.0;

double CL = TVCL * pow(CRCL / REF_CRCL, CRCL_CL) * pow(DIAL_CL, dialysis) * exp(ETA(1) + ETA1);
double V1 = TVV1 * (wt_safe / REF_WT) * pow(DIAL_V1, dialysis) * exp(ETA(2) + ETA2);
double V2 = TVV2 * exp(ETA(3) + ETA3);
double Q = TVQ;

$ODE
double CP = CENT / V1;
double CPP = PERIPH / V2;
dxdt_CENT = -CL * CP - Q * (CP - CPP);
dxdt_PERIPH = Q * (CP - CPP);

$TABLE
double DV = CP * (1.0 + EPS(1)) + EPS(2);

$CAPTURE @annotated
DV   : Simulated plasma concentration (mg/L)
CP   : Predicted plasma concentration (mg/L)
CL   : Individual clearance (L/h)
V1   : Individual central volume (L)
V2   : Individual peripheral volume (L)
Q    : Intercompartmental clearance (L/h)
CRCL : Cockcroft-Gault creatinine clearance (mL/min)
