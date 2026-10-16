$PROB
# Revilla et al. Br J Clin Pharmacol. 2010;70:201-212.
# doi:10.1111/j.1365-2125.2010.03679.x
# One-compartment IV model for adult ICU patients.
# CRCL is the measured 24 h creatinine clearance used in the article.

$PARAM @annotated
TVCLCR : 0.67  : Slope of measured CRCL on clearance
TVV1   : 0.82  : Central volume (L/kg)
AGE_CL : -0.24 : Exponent of age on clearance
CREAT_V: 2.49  : Volume multiplier when serum creatinine is above 1 mg/dL

ETA1: 0 : Clearance  (L/h)
ETA2: 0 : Central volume (L)

$PARAM @annotated @covariates
AGE   : 61.1 : Age (years)
CRCL  : 74.7 : Measured 24 h creatinine clearance (mL/min)
CREAT : 123.8 : Serum creatinine (micromol/L)
WT    : 73.0 : Total body weight (kg)

$OMEGA 0.09 0.053

$SIGMA @annotated
PROP : 0       : Proportional residual variance (not used in the article)
ADD  : 17.8929 : Additive residual variance (SD 4.23 mg/L)

$CMT @annotated
CENT  : Central compartment (mg/L)[ADM, OBS]

$TABLE
double DV = (CENT/V1) * (1 + EPS(1)) + EPS(2);

int i = 0;
while(DV<0 && i <100) {
simeps();
DV = (CENT/V1) * (1 + EPS(1)) + EPS(2);
++i;
}

$MAIN
double A = (CREAT / 88.4 > 1 ? 1 : 0);
double clcr_per_kg = CRCL / WT;
double CL = ((TVCLCR * clcr_per_kg + pow(AGE, AGE_CL)) * WT * 0.06) * exp(ETA1 + ETA(1));
double V1 = (TVV1 * pow(CREAT_V, A) * WT) * exp(ETA2 + ETA(2));

$ODE
dxdt_CENT   =  - (CL / V1)  * CENT ;

$CAPTURE DV CL V1 clcr_per_kg
