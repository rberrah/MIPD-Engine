$PROB
// Daptomycin population PK model (Garreau et al., 2021).
// Population: adults treated for bone and joint infection.
// Structure: two compartments and linear elimination after IV administration.
// The publication estimated inter-occasion variability on CL and V1; this
// session model retains between-subject variability only because occasions are
// not represented by independent random effects in the current MAP workflow.
// Reference: doi:10.1093/jac/dkab006.

$PARAM @annotated
TVCL : 0.365 : Baseline clearance before covariate effects (L/h)
TVV1 : 3.59  : Baseline central volume before covariate effects (L)
TVQ  : 0.752 : Typical intercompartmental clearance (L/h)
TVV2 : 4.71  : Typical peripheral volume (L)
CLCR_CL : 0.430  : Creatinine clearance effect on CL
MALE_CL : 0.232  : Male effect on CL
AGE_V1  : 0.263  : Age effect on V1
WT_V1   : 0.603  : Body-weight effect on V1
RIF_V1  : -0.121 : Rifampicin co-administration effect on V1
MALE_V1 : 0.117  : Male effect on V1
REF_CLCR : 109.0 : Reference creatinine clearance (mL/min)
REF_AGE  : 60.4  : Reference age (years)
REF_WT   : 79.2  : Reference body weight (kg)
ETA1 : 0.0 : Posterior ETA offset on CL
ETA2 : 0.0 : Posterior ETA offset on V1
ETA3 : 0.0 : Posterior ETA offset on Q
ETA4 : 0.0 : Posterior ETA offset on V2

$PARAM @covariates @annotated
CLCR : 109.0 : Cockcroft-Gault creatinine clearance (mL/min)
AGE  : 60.4  : Age (years)
WT   : 79.2  : Body weight (kg)
SEX  : 0.0   : Sex (0 = male, 1 = female)
RIF  : 0.0   : Rifampicin co-administration (0 = no, 1 = yes)

$CMT @annotated
CENT : Central compartment (mg) [ADM, OBS]
PERI : Peripheral compartment (mg)

$OMEGA @annotated @diagonal
ETA_CL : 0.038416 : Between-subject variance on CL (SD 0.196)
ETA_V1 : 0.015129 : Between-subject variance on V1 (SD 0.123)
ETA_Q  : 1.276900 : Between-subject variance on Q (SD 1.13)
ETA_V2 : 0.084100 : Between-subject variance on V2 (SD 0.29)

$SIGMA @annotated @diagonal
PROP : 0.0496 : Proportional residual variance
ADD  : 2.55   : Additive residual variance

$MAIN
double clcr_safe = (CLCR > 0.0) ? CLCR : REF_CLCR;
double age_safe = (AGE > 0.0) ? AGE : REF_AGE;
double wt_safe = (WT > 0.0) ? WT : REF_WT;
double male = (SEX < 0.5) ? 1.0 : 0.0;
double rifampicin = (RIF >= 0.5) ? 1.0 : 0.0;

double CL = TVCL * exp(CLCR_CL * (clcr_safe / REF_CLCR) + MALE_CL * male) * exp(ETA(1) + ETA1);
double V1 = TVV1 * exp(AGE_V1 * (age_safe / REF_AGE) + WT_V1 * (wt_safe / REF_WT) + RIF_V1 * rifampicin + MALE_V1 * male) * exp(ETA(2) + ETA2);
double Q = TVQ * exp(ETA(3) + ETA3);
double V2 = TVV2 * exp(ETA(4) + ETA4);

$ODE
double CP = CENT / V1;
double CPP = PERI / V2;
dxdt_CENT = -CL * CP - Q * (CP - CPP);
dxdt_PERI = Q * (CP - CPP);

$TABLE
double DV = CP * (1.0 + EPS(1)) + EPS(2);

$CAPTURE @annotated
DV : Simulated plasma concentration (mg/L)
CP : Predicted plasma concentration (mg/L)
CL : Individual clearance (L/h)
V1 : Individual central volume (L)
Q  : Individual intercompartmental clearance (L/h)
V2 : Individual peripheral volume (L)
