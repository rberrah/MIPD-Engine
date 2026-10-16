$PROB
// Cefepime population PK model (Jonckheere et al., 2016).
// Population: critically ill adults, with and without intermittent hemodialysis.
// Structure: two compartments with separate renal, non-renal and dialysis clearance.
// CRCL is the Cockcroft-Gault creatinine clearance in mL/min.
// Reference: doi:10.1093/jac/dkw171.

$PARAM @annotated
TVV1      : 18.3  : Typical central volume (L)
TVV2      : 11.1  : Typical peripheral volume (L)
TVQ       : 6.63  : Typical intercompartmental clearance (L/h)
THETA1    : 2.88  : Renal clearance at Cockcroft-Gault 3.40 L/h (L/h)
THETA2    : 0.368 : Linear Cockcroft-Gault effect on renal clearance (h/L)
REF_CG_LH : 3.40  : Centering Cockcroft-Gault clearance (L/h)
CLDIAL    : 5.74  : Clearance during intermittent hemodialysis (L/h)
CLOTHER   : 0.87  : Non-renal clearance (L/h)
ETA1      : 0.0   : Posterior ETA offset on renal clearance
ETA2      : 0.0   : Posterior ETA offset on central volume

$PARAM @covariates @annotated
CRCL : 56.67 : Cockcroft-Gault creatinine clearance (mL/min)
IHD  : 0.0   : Intermittent hemodialysis at this time (0 = no, 1 = yes)

$CMT @annotated
CENT   : Central compartment (mg) [ADM, OBS]
PERIPH : Peripheral compartment (mg)

$OMEGA @annotated @diagonal
ETA_CLREN : 0.2191 : Between-subject variance on renal clearance (CV 49.5%)
ETA_V1    : 0.1463 : Between-subject variance on central volume (CV 39.7%)

$SIGMA @annotated @diagonal
PROP : 0.110224 : Proportional plasma residual variance (CV 33.2%)
ADD  : 0.0      : Additive residual variance fixed to zero

$MAIN
double crcl_safe = (CRCL > 0.0) ? CRCL : 0.0;
double cg_lh = crcl_safe * 0.06;
double renal_effect = 1.0 + THETA2 * (cg_lh - REF_CG_LH);
if (renal_effect < 0.0) renal_effect = 0.0;

double CLREN = (IHD >= 0.5) ? 0.0 : THETA1 * renal_effect * exp(ETA(1) + ETA1);
double CLHD = (IHD >= 0.5) ? CLDIAL : 0.0;
double CL = CLREN + CLOTHER + CLHD;
double V1 = TVV1 * exp(ETA(2) + ETA2);
double V2 = TVV2;
double Q = TVQ;

$ODE
double CP = CENT / V1;
double CPP = PERIPH / V2;
dxdt_CENT = -CL * CP - Q * (CP - CPP);
dxdt_PERIPH = Q * (CP - CPP);

$TABLE
double DV = CP * (1.0 + EPS(1)) + EPS(2);

$CAPTURE @annotated
DV    : Simulated plasma concentration (mg/L)
CP    : Predicted plasma concentration (mg/L)
CL    : Total individual clearance (L/h)
CLREN : Individual renal clearance (L/h)
CLHD  : Dialysis clearance active at this time (L/h)
V1    : Individual central volume (L)
V2    : Peripheral volume (L)
Q     : Intercompartmental clearance (L/h)
cg_lh : Cockcroft-Gault creatinine clearance converted to L/h
