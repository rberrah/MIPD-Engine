$PROB
// Oral levofloxacin population PK model (Canoui et al., 2022).
// Population: adults treated for bone and joint infection.
// Structure: one compartment, first-order absorption and elimination.
// Reference: doi:10.1093/jac/dkac031.

$PARAM @annotated
TVCL     : 5.57   : Typical apparent clearance (L/h)
TVV      : 96.3   : Typical apparent volume (L)
KA       : 1.6    : Fixed absorption rate constant (1/h)
CLCR_CL  : 0.684  : Creatinine clearance exponent on CL/F
AGE_CL   : -0.312 : Age exponent on CL/F
REF_CLCR : 90.0   : Reference Cockcroft-Gault clearance (mL/min)
REF_AGE  : 61.5   : Reference age (years)
ETA1     : 0.0    : Posterior ETA offset on V/F
ETA2     : 0.0    : Posterior ETA offset on CL/F

$PARAM @covariates @annotated
CLCR : 90.0 : Cockcroft-Gault creatinine clearance (mL/min)
AGE  : 61.5 : Age (years)

$CMT @annotated
GUT  : Oral dosing compartment (mg) [ADM]
CENT : Central compartment (mg) [OBS]

$OMEGA @annotated @diagonal
ETA_V  : 0.080089 : Between-subject variance on V/F (SD 0.283)
ETA_CL : 0.134689 : Between-subject variance on CL/F (SD 0.367)

$SIGMA @annotated @diagonal
PROP : 0.274576 : Multiplicative residual variance (SD 0.524)
ADD  : 0.0      : Additive residual variance fixed to zero

$MAIN
double clcr_safe = (CLCR > 0.0) ? CLCR : REF_CLCR;
double age_safe = (AGE > 0.0) ? AGE : REF_AGE;
double V = TVV * exp(ETA(1) + ETA1);
double CL = TVCL * pow(clcr_safe / REF_CLCR, CLCR_CL) * pow(age_safe / REF_AGE, AGE_CL) * exp(ETA(2) + ETA2);

$ODE
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - CL * (CENT / V);

$TABLE
double CP = CENT / V;
double DV = CP * (1.0 + EPS(1)) + EPS(2);

$CAPTURE @annotated
DV : Simulated plasma concentration (mg/L)
CP : Predicted plasma concentration (mg/L)
CL : Individual apparent clearance (L/h)
V  : Individual apparent volume (L)
KA : Absorption rate constant (1/h)
