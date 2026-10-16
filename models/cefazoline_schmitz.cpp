$PROB
Cefazolin PopPK Model (Schmitz et al. 2015)
Reference: Schmitz ML et al. (2015) [PMCID: PMC4468723]
Population: Pediatric surgical patients (10-12y) and Adults
Structure: 2 compartments
Adaptation: the published CL and V1 inter-occasion effects are not represented.
The published residual SD is 0.845 + sqrt(0.0112) * prediction (combined1);
the two-EPS combined2 form below is retained for engine compatibility, so this
file is simulation-only and is not exposed as an article-faithful MAP prior.
$PARAM @annotated
TVCL_NR  : 0.153 : Typical Non-Renal Clearance (L/h/70kg)
TVCL_R   : 3.63  : Typical Renal Clearance (L/h/70kg)
TVVC     : 4.72  : Typical Central Volume (L/70kg)
TVVP     : 3.70  : Typical Peripheral Volume (L/70kg)
TVQ      : 7.89  : Typical Inter-compartmental Clearance (L/h/70kg)
POW_RCC  : 0.337 : Power of CLCR on Renal CL
SHFT_PED : 0.335 : Pediatric reduction fraction on Renal CL
ALLO_CL  : 0.75  : Allometric exponent for CL and Q
ALLO_V   : 1.00  : Allometric exponent for Vc and Vp
ETA1 : 0.0 : Mapbayr ETA on CL
ETA2 : 0.0 : Mapbayr ETA on VC
$PARAM @annotated @covariates
WT       : 45.6  : Weight (kg)
AGE      : 12    : Age (years)
CLCR     : 90.0  : Creatinine clearance used by the published model (mL/min/1.73m2)

$CMT @annotated
CENT   : Central Compartment (mg) [ADM, OBS]
PERIPH : Peripheral Compartment (mg)
$OMEGA @annotated
IIV_CL : 0.0122 : Variance on CL
IIV_VC : 0.0351 : Variance on VC
$SIGMA @annotated
PROP : 0.0112 : Proportional error variance
ADD  : 0.714025 : Squared published additive component (0.845 mg/L; combined2 adaptation)
$MAIN
double wt_safe    = (WT < 1.0)     ? 1.0  : WT;
double PEDS = (AGE < 18.0) ? 1 : 0;

double clcr_safe = (CLCR > 0.0) ? CLCR : 90.0;
double allo_cl   = pow(wt_safe / 70.0, ALLO_CL);
double allo_v    = pow(wt_safe / 70.0, ALLO_V);
double renal_fac = pow(clcr_safe / 90.0, POW_RCC);
double peds_fac  = (1.0 - SHFT_PED * PEDS);
double cl_ind = TVCL_NR + (TVCL_R * renal_fac * peds_fac);
double CL = cl_ind * allo_cl * exp(ETA(1) + ETA1);
double VC = TVVC   * allo_v  * exp(ETA(2) + ETA2);
double Q  = TVQ    * allo_cl;
double VP = TVVP   * allo_v;
$ODE
dxdt_CENT   = -(CL/VC)*CENT - (Q/VC)*CENT + (Q/VP)*PERIPH;
dxdt_PERIPH =  (Q/VC)*CENT - (Q/VP)*PERIPH;
$TABLE
double CP = CENT / VC;
double DV = CP * (1.0 + EPS(1)) + EPS(2);
$CAPTURE @annotated
CP       : Predicted Total Concentration (mg/L)
CL       : Individual Clearance (L/h)
VC       : Individual Central Volume (L)
Q        : Individual Inter-compartmental Clearance (L/h)
VP       : Individual Peripheral Volume (L)
clcr_safe: Creatinine clearance used by the model (mL/min/1.73m2)
DV       : Simulated Total Concentration (mg/L)
