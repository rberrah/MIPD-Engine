// Standalone PK adaptation for TDM.
// The DDI mechanism placeholders were removed: this file describes baseline PK only.
// No patient data or pasted model is persisted by the application.
$PLUGIN tad
$SET end = 168, delta = 0.1

$PARAM @annotated
// --- CICLOSPORINE BASE (Press et al. 2010, Renal Transplant) ---
TVCL_CSA : 15.0  : CL app CsA (L/h)
TVV1_CSA : 56.0  : V1 app CsA (L)
TVQ_CSA  : 14.0  : Q CsA (L/h)
TVV2_CSA : 125.0 : V2 CsA (L)
TVKA_CSA : 2.0   : Ka CsA (1/h)
TVF_CSA  : 0.5   : F CsA
TVTLAG   : 1.0   : Transit lag time (h)

// --- EFFET DES COVARIABLES ---
WTCL  : 0.75  : Allometric WT on CL
WTV   : 1.0   : Allometric WT on V
WTREF : 76.0  : Published reference body weight (kg)
PREDNI_KA : 0.45 : Ka multiplier with prednisolone >= 20 mg/day
PREDNI_F  : 0.78 : F multiplier with prednisolone >= 20 mg/day

// Fraction of ciclosporine apparent oral clearance that CYP3A actually drives.
// Ciclosporine is extensively CYP3A-metabolised in vitro, but the OBSERVED
// interactions are far smaller than for the other immunosuppressants: the
// voriconazole label reports a 1.7-fold rise in ciclosporine AUC, against
// 3.2-fold for tacrolimus and 11-fold for sirolimus. This model previously
// reused the tacrolimus value of 0.95 and predicted a 7.6-fold rise, which
// would have advised cutting the dose to a seventh where the label says half.
// FM_CYP3A is therefore calibrated on that clinical figure. Raise it to make
// every CYP3A interaction proportionally stronger.
FM_CYP3A : 0.50 : CYP3A-dependent fraction of CL/F

// --- ETAs---
ETA1:0:CL
ETA2:0:V1
ETA3:0:Ka

// Covariates MUST live in an @covariates block: mapbayr refuses a data set whose
// columns match a plain $PARAM entry. The perpetrator and mechanism parameters
// are declared here too, because a perpetrator may bring a covariate of its own.
$PARAM @annotated @covariates
WT    : 76    : Body weight (kg)
PREDNI : 0    : Prednisolone >= 20 mg/day (0 = no, 1 = yes)


$OMEGA 0.03 0.12 0.09
$SIGMA 0.07 0.0

$CMT @annotated
DEPOT_CSA : CsA depot [ADM]
TR1_CSA   : Transit 1
CENT_CSA  : CsA central [OBS]
PERI_CSA  : CsA peripheral

$MAIN
// Initialisation des pools d'enzymes si necessaire

// Calculs de base de la Ciclosporine
double CL0 = TVCL_CSA * pow(WT/WTREF, WTCL) * exp(ETA1 + ETA(1));
double V1  = TVV1_CSA * pow(WT/WTREF, WTV) * exp(ETA2 + ETA(2));
double Q   = TVQ_CSA;
double V2  = TVV2_CSA;
double KA  = TVKA_CSA * (PREDNI > 0.5 ? PREDNI_KA : 1.0) * exp(ETA3 + ETA(3));
double KTR = 1.0 / TVTLAG;

// Interaction-only clearance split. Press et al. did not retain a CYP3A5
// genotype effect, so baseline PK and this split do not depend on genotype.
double f3A5   = 0.0;
double f3A4   = FM_CYP3A;
double fm_sys = 1.0 - FM_CYP3A;

F_DEPOT_CSA = TVF_CSA * (PREDNI > 0.5 ? PREDNI_F : 1.0);


$ODE
// Variables de communication standardisees (Par defaut = Pas d'effet)
double ACT_CYP3A4 = 1.0;
double ACT_CYP3A5 = 1.0;
double ACT_PGP    = 1.0;
double F_BOOST    = 1.0;
double CP_PERP    = 0.0;

// 1. Calcul de la Pharmacocinetique du Perpetrateur

// 2. Calcul des Mecanismes d'interaction (Modifient les ACT_)

// 3. Application a la clairance de la Ciclosporine
double CL_CSA_DDI = CL0 * (f3A4 * ACT_CYP3A4 + f3A5 * ACT_CYP3A5 + fm_sys);

dxdt_DEPOT_CSA = -KTR * DEPOT_CSA;
dxdt_TR1_CSA   = KTR * DEPOT_CSA - KA * TR1_CSA;
dxdt_CENT_CSA  = KA * TR1_CSA * F_BOOST - (CL_CSA_DDI + Q)*(CENT_CSA/V1) + Q*(PERI_CSA/V2);
dxdt_PERI_CSA  = Q*(CENT_CSA/V1) - Q*(PERI_CSA/V2);

$TABLE
double CONC = CENT_CSA / V1 * 1000.0; // ng/mL
double DV = CONC * (1.0 + EPS(1)) + EPS(2);
double ratio = (CL_CSA_DDI / CL0) / F_BOOST;
double CL_OUT = CL_CSA_DDI;

$CAPTURE DV CONC ratio CL_CSA_DDI CP_PERP CL_OUT
