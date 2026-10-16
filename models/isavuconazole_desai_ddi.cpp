$PLUGIN tad
$PROB
// Isavuconazole population PK model (Desai et al.)
// Article: Desai A et al. Antimicrob Agents Chemother. 2016;60:5483-5491.
// DOI: 10.1128/AAC.02819-15
// Population: Healthy adults and adults with invasive fungal infections in phase 1 and SECURE phase 3 trials.
// Implementation: deterministic PK module adapted for standalone TDM.
// Published Weibull absorption, race and patient-status effects, BMI effect and reported IIV are retained. Patient CL variability is the default and the healthy/patient residual branches are not separated, so this remains simulation-only.

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
ETA2 : 0 : Random effect on Q
ETA3 : 0 : Random effect on peripheral volume
ETA4 : 0 : Random effect on Weibull RA
ETA5 : 0 : Random effect on Weibull GAM1

$PARAM @annotated @covariates
BMI : 23.6 : Body mass index (kg/m2)
ASIAN : 0 : Asian race indicator
PATIENT : 1 : SECURE patient (1) or healthy phase-1 participant (0)

$OMEGA @annotated
IIV_CL : 0.3290 : Patient CL IIV, log(1 + 0.6244^2)
IIV_Q : 0.2158 : Q IIV, log(1 + 0.4909^2)
IIV_V2 : 0.0962 : Vp IIV, log(1 + 0.3178^2)
IIV_RA : 0.1502 : RA IIV, log(1 + 0.4024^2)
IIV_GAM1 : 0.1900 : GAM1 IIV, log(1 + 0.4571^2)

$SIGMA @annotated
PROP : 0.201960 : Published exponential residual SD 0.4494 squared
ADD  : 0.0  : Additive residual variance fixed to zero

$CMT @annotated
GUT : Oral depot [ADM]
CENT : Central compartment [OBS]
PERI : Peripheral compartment

$MAIN
double TDOSE = TIME - self.tad();
double CL = (ASIAN > 0.5 ? 1.51 : TVCL) * exp(ETA1 + ETA(1));
double Q = TVQ * exp(ETA2 + ETA(2));
double V2_BASE = PATIENT > 0.5 ? TVV2 : 260.0;
double V2 = V2_BASE * (1.0 + EFF_BMI_V2 * (BMI - BMI_REF)) * exp(ETA3 + ETA(3));
double RA_IND = RA * exp(ETA4 + ETA(4));
double GAM1_IND = GAM1 * exp(ETA5 + ETA(5));
if (V2 < 1.0) V2 = 1.0;

$ODE
double TAD = SOLVERTIME - TDOSE;
double KA = 0.0;
if (TAD > 0.0) KA = KAMAX * (1.0 - exp(-pow(RA_IND * TAD, GAM1_IND)));
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - ((CL + Q) / TVV1) * CENT + (Q / V2) * PERI;
dxdt_PERI = (Q / TVV1) * CENT - (Q / V2) * PERI;

$TABLE
double IPRED = CENT / TVV1;
double DV = IPRED * (1.0 + EPS(1));

$CAPTURE DV CL
