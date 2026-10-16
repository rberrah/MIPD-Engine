$PLUGIN tad
$PROB
// Isavuconazole population PK model (Desai et al.)
// Article: Desai A et al. Antimicrob Agents Chemother. 2016;60:5483-5491.
// DOI: 10.1128/AAC.02819-15
// Population: Healthy adults and adults with invasive fungal infections in phase 1 and SECURE phase 3 trials.
// Implementation: deterministic PK module adapted for standalone TDM.
// The published time-dependent Weibull absorption is retained; the residual prior is an engineering adaptation for mapbayr.

$PARAM @annotated
TVCL : 2.36 : Typical clearance (L/h)
TVV1 : 49.1 : Typical central volume (L)
TVV2 : 417 : Typical peripheral volume (L)
TVQ : 26.6 : Typical intercompartmental clearance (L/h)
KAMAX : 1.08 : Maximum absorption rate (1/h)
RA : 0.72 : Weibull scale parameter
GAM1 : 4.88 : Weibull shape parameter
EFF_BMI_V2 : 0.060 : Linear BMI effect on peripheral volume
BMI_REF : 24.8 : Reference BMI (kg/m2)
ETA1 : 0 : MAP random effect on clearance

$PARAM @annotated @covariates
BMI : 23.6 : Body mass index (kg/m2)

$OMEGA @annotated
IIV_CL : 0.09 : Engineering prior variance on clearance

$SIGMA @annotated
PROP : 0.04 : Engineering prior proportional residual variance
ADD  : 0.0  : Additive residual variance fixed to zero

$CMT @annotated
GUT : Oral depot [ADM]
CENT : Central compartment [OBS]
PERI : Peripheral compartment

$MAIN
double TDOSE = TIME - self.tad();
double CL = TVCL * exp(ETA1 + ETA(1));
double V2 = TVV2 * (1.0 + EFF_BMI_V2 * (BMI - BMI_REF));

$ODE
double TAD = SOLVERTIME - TDOSE;
double KA = 0.0;
if (TAD > 0.0) KA = KAMAX * (1.0 - exp(-pow(RA * TAD, GAM1)));
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - ((CL + TVQ) / TVV1) * CENT + (TVQ / V2) * PERI;
dxdt_PERI = (TVQ / TVV1) * CENT - (TVQ / V2) * PERI;

$TABLE
double IPRED = CENT / TVV1;
double DV = IPRED * (1.0 + EPS(1));

$CAPTURE DV CL
