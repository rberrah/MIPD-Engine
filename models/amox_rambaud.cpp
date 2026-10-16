$PROB
// Amoxicillin PK model (Rambaud et al., 2020).
// Population: adults receiving continuous IV amoxicillin for infective endocarditis.
// Structure: two compartments parameterized with microconstants.
// The source was estimated with the non-parametric NPAG algorithm. The OMEGA
// terms below are an approximate log-normal MAP prior and do not reproduce the
// published discrete support-point distribution.
// Reference: doi:10.1093/jac/dkaa232.

$PARAM @annotated
TVKE1      : 1.688  : Baseline elimination rate constant (1/h)
TVKE2      : 0.946  : Absolute GFR exponent on elimination
TVKCP      : 0.247  : Central-to-peripheral rate constant (1/h)
TVKPC      : 11.427 : Peripheral-to-central rate constant (1/h)
TVV        : 5.698  : Central volume (L)
REF_AGFR   : 64.92  : Reference absolute GFR (mL/min)
ERR_C0     : 2.5    : Published assay SD intercept (mg/L)
ERR_C1     : 0.15   : Published assay SD slope
ERR_LAMBDA : 5.65   : Published final lambda process-noise term
ETA1 : 0.0 : Posterior ETA offset on KE1
ETA2 : 0.0 : Posterior ETA offset on V
ETA3 : 0.0 : Posterior ETA offset on KCP
ETA4 : 0.0 : Posterior ETA offset on KPC
ETA5 : 0.0 : Posterior ETA offset on KE2

$PARAM @covariates @annotated
CREAT : 88.4 : Serum creatinine (micromol/L)
SEX   : 0.0  : Sex (0 = male, 1 = female)
BLACK : 0.0  : Black ethnicity factor in the published 2009 CKD-EPI equation (0 = no, 1 = yes)
AGE   : 60.0 : Age (years)
BSA   : 1.73 : Body surface area (m2)

$CMT @annotated
CENT   : Central compartment (mg) [ADM, OBS]
PERIPH : Peripheral compartment (mg)

$OMEGA @annotated @diagonal
ETA_KE1 : 0.0367 : Approximate MAP variance on KE1
ETA_V   : 0.0615 : Approximate MAP variance on V
ETA_KCP : 1.1100 : Approximate MAP variance on KCP
ETA_KPC : 0.1940 : Approximate MAP variance on KPC
ETA_KE2 : 0.0964 : Approximate MAP variance on KE2

$SIGMA @annotated @diagonal
ASSAY   : 1.0 : Unit variance scaled by the published SD polynomial
PROCESS : 1.0 : Unit variance scaled by the published lambda term

$MAIN
double creat_mgdl = (CREAT > 0.0 ? CREAT : 88.4) / 88.4;
double age_safe = (AGE > 0.0) ? AGE : 60.0;
double bsa_safe = (BSA > 0.0) ? BSA : 1.73;
double female = (SEX >= 0.5) ? 1.0 : 0.0;
double black = (BLACK >= 0.5) ? 1.0 : 0.0;
double kappa = (female > 0.5) ? 0.7 : 0.9;
double alpha = (female > 0.5) ? -0.329 : -0.411;
double ratio = creat_mgdl / kappa;
double ratio_low = (ratio < 1.0) ? ratio : 1.0;
double ratio_high = (ratio > 1.0) ? ratio : 1.0;
double sex_factor = (female > 0.5) ? 1.018 : 1.0;
double ethnicity_factor = (black > 0.5) ? 1.159 : 1.0;

double EGFR = 141.0 * pow(ratio_low, alpha) * pow(ratio_high, -1.209) *
              pow(0.993, age_safe) * sex_factor * ethnicity_factor;
double AGFR = EGFR * bsa_safe / 1.73;
if (AGFR < 0.1) AGFR = 0.1;

double KE1 = TVKE1 * exp(ETA(1) + ETA1);
double KE2 = TVKE2 * exp(ETA(5) + ETA5);
double K10 = KE1 * pow(AGFR / REF_AGFR, KE2);
double V = TVV * exp(ETA(2) + ETA2);
double KCP = TVKCP * exp(ETA(3) + ETA3);
double KPC = TVKPC * exp(ETA(4) + ETA4);
double CL = K10 * V;
double Q = KCP * V;
double VP = Q / KPC;

$ODE
dxdt_CENT = KPC * PERIPH - (K10 + KCP) * CENT;
dxdt_PERIPH = KCP * CENT - KPC * PERIPH;

$TABLE
double CP = CENT / V;
// Pmetrics weighted observations by sqrt((C0 + C1*obs)^2 + lambda^2).
// Prediction replaces the unknown observation in this mrgsolve residual approximation.
double DV = CP + (ERR_C0 + ERR_C1 * CP) * EPS(1) + ERR_LAMBDA * EPS(2);

$CAPTURE @annotated
DV   : Simulated total plasma concentration (mg/L)
CP   : Predicted total plasma concentration (mg/L)
CL   : Individual clearance (L/h)
V    : Individual central volume (L)
Q    : Individual intercompartmental clearance (L/h)
VP   : Individual peripheral volume (L)
K10  : Individual elimination rate constant (1/h)
KCP  : Central-to-peripheral rate constant (1/h)
KPC  : Peripheral-to-central rate constant (1/h)
EGFR : 2009 CKD-EPI GFR indexed to 1.73 m2
AGFR : Absolute GFR used by the published model (mL/min)
