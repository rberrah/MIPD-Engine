$PROB
// Amoxicillin population PK model (Carlier et al., 2013).
// Population: critically ill adults receiving intravenous amoxicillin.
// Structure: two compartments, linear elimination and zero-order IV input.
// Renal covariate: measured 24-hour urinary creatinine clearance (mL/min).
// Reference: doi:10.1093/jac/dkt240.

$PARAM @annotated
TVCL     : 10.0  : Typical clearance at CRCL 102 mL/min (L/h)
TVVC     : 13.7  : Typical central volume (L)
TVVP     : 13.7  : Typical peripheral volume (L)
TVQ      : 15.6  : Typical intercompartmental clearance (L/h)
REF_CRCL : 102.0 : Median 24-hour urinary creatinine clearance (mL/min)

// Posterior ETA offsets used by mapbayr.
ETA1 : 0.0 : Posterior ETA on CL
ETA2 : 0.0 : Posterior ETA on VC

$PARAM @covariates @annotated
CRCL : 102.0 : Measured 24-hour urinary creatinine clearance (mL/min)

$OMEGA @annotated @diagonal
ETA_CL : 0.148 : Between-subject variance on CL
ETA_VC : 0.140 : Between-subject variance on VC

$SIGMA @annotated @diagonal
PROP : 0.0473 : Exponential proportional residual variance
ADD  : 0.0    : Additive residual variance fixed to zero

$CMT @annotated
CENT   : Central amount (mg) [ADM, OBS]
PERIPH : Peripheral amount (mg)

$MAIN
// The publication used measured urinary CRCL, not an estimated GFR equation.
double CRCL_USED = (CRCL > 0.0) ? CRCL : REF_CRCL;

double CL = TVCL * (CRCL_USED / REF_CRCL) * exp(ETA(1) + ETA1);
double VC = TVVC * exp(ETA(2) + ETA2);
double VP = TVVP;
double Q  = TVQ;

double K10 = CL / VC;
double K12 = Q  / VC;
double K21 = Q  / VP;

$ODE
dxdt_CENT   = K21 * PERIPH - (K10 + K12) * CENT;
dxdt_PERIPH = K12 * CENT - K21 * PERIPH;

$TABLE
double CP = CENT / VC;
double DV = CP * exp(EPS(1)) + EPS(2);

$CAPTURE @annotated
DV        : Simulated total plasma concentration (mg/L)
CP        : Prediction without residual error (mg/L)
CL        : Individual clearance (L/h)
VC        : Individual central volume (L)
Q         : Intercompartmental clearance (L/h)
VP        : Peripheral volume (L)
K10       : Elimination rate constant (1/h)
K12       : Central-to-peripheral rate constant (1/h)
K21       : Peripheral-to-central rate constant (1/h)
CRCL_USED : Measured 24-hour urinary creatinine clearance used by the model (mL/min)
