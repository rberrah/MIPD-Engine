$PROB
Cefazolin PopPK Model (Lanoiselee et al. 2021)
Reference: Lanoiselee J et al. Sci Rep. 2021;11(1):19763.
Population: Adults (24-91 years), total hip arthroplasty
Structure: 2 compartments
$PARAM @annotated
TVCL       : 2.86  : Typical Clearance (L/h)
THETA_CLCR : 0.79  : Power exponent of CrCl on CL
TVVC       : 5.20  : Typical Central Volume (L)
TVVP       : 4.56  : Typical Peripheral Volume (L)
TVQ        : 10.9  : Typical Inter-compartmental Clearance (L/h)
REF_CLCR   : 80.0  : Reference CrCl (mL/min/1.73m2)
FU         : 0.20  : Unbound fraction
ETA1 : 0 : Mapbayr ETA on CL
ETA2 : 0 : Mapbayr ETA on VC
ETA3 : 0 : Mapbayr ETA on Q
ETA4 : 0 : Mapbayr ETA on VP
$PARAM @annotated @covariates
AGE  : 67    : Age (years)
WT   : 76.0  : Total Body Weight (kg)
HT   : 168.0 : Height (cm)
CREAT  : 80.0  : Serum Creatinine (µmol/L)
SEX  : 0     : Sex (0=Male, 1=Female)
$CMT @annotated
CENT   : Central Compartment (mg) [ADM, OBS]
PERIPH : Peripheral Compartment (mg)
$OMEGA @block @annotated
IIV_CL  : 0.097719          : Variance on CL
IIV_VC  : 0.137660 0.281483 : Covariance CL-VC / Variance on VC
$OMEGA @annotated
IIV_Q   : 0.361541 : Variance on Q
IIV_VP  : 0.009950 : Variance on VP
$SIGMA @annotated
PROP : 0.014298 : Proportional error variance
ADD  : 0.000000 : Additive error variance
$MAIN
double wt_safe  = (WT  < 1.0)  ? 1.0  : WT;
double ht_safe  = (HT  < 50.0) ? 50.0 : HT;
double scr_safe = (CREAT < 10.0) ? 10.0 : CREAT;
double scr_mgdl = scr_safe / 88.4;
double CLCR = REF_CLCR;
double kappa;
double alpha;
double ratio;
if (SEX == 1) {
    kappa = 0.7;
    alpha = (scr_mgdl <= kappa) ? -0.329 : -1.209;
    ratio = scr_mgdl / kappa;
    CLCR = 144.0 * pow(ratio, alpha) * pow(0.993, AGE);
} else {
    kappa = 0.9;
    alpha = (scr_mgdl <= kappa) ? -0.411 : -1.209;
    ratio = scr_mgdl / kappa;
    CLCR = 141.0 * pow(ratio, alpha) * pow(0.993, AGE);
}
CLCR = (CLCR > 200.0) ? 200.0 : CLCR;
CLCR = (CLCR < 1.0)   ? 1.0   : CLCR;
double CL = TVCL * pow(CLCR / REF_CLCR, THETA_CLCR) * exp(ETA(1) + ETA1);
double VC = TVVC * exp(ETA(2) + ETA2);
double Q  = TVQ  * exp(ETA(3) + ETA3);
double VP = TVVP * exp(ETA(4) + ETA4);
$ODE
dxdt_CENT   = -(CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERIPH;
dxdt_PERIPH =  (Q/VC)*CENT - (Q/VP)*PERIPH;
$TABLE
double CP      = CENT / VC;
double CP_FREE = CP * FU;
double DV      = CP * (1.0 + EPS(1)) + EPS(2);
$CAPTURE @annotated
CP       : Predicted Total Concentration (mg/L)
CP_FREE  : Predicted Unbound Concentration (mg/L)
CL       : Individual Clearance (L/h)
VC       : Individual Central Volume (L)
Q        : Individual Inter-compartmental Clearance (L/h)
VP       : Individual Peripheral Volume (L)
CLCR     : Calculated CrCl (mL/min/1.73m2)
DV       : Simulated Total Concentration (mg/L)
