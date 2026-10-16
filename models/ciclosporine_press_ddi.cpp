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
CYPCL : 1.3   : Ratio Genotype CYP3A5

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
WT    : 70    : Poids (kg)
CYP   : 0     : CYP3A5 Status (0=non-expr, 1=expr)


$OMEGA 0.03 0.12 0.09
$SIGMA 0.07 25.0

$CMT @annotated
DEPOT_CSA : CsA depot [ADM]
TR1_CSA   : Transit 1
CENT_CSA  : CsA central [OBS]
PERI_CSA  : CsA peripheral

$MAIN
// Initialisation des pools d'enzymes si necessaire

// Calculs de base de la Ciclosporine
double CL0 = TVCL_CSA * pow(WT/70.0, WTCL) * pow(CYPCL, CYP) * exp(ETA1 + ETA(1));
double V1  = TVV1_CSA * pow(WT/70.0, WTV) * exp(ETA2 + ETA(2));
double Q   = TVQ_CSA;
double V2  = TVV2_CSA;
double KA  = TVKA_CSA * exp(ETA3 + ETA(3));
double KTR = 1.0 / TVTLAG;

// Fractions metaboliques de la Ciclosporine. The CYP3A-dependent part is split
// between CYP3A4 and CYP3A5 according to genotype; everything else is left
// untouched by a CYP3A inhibitor, so the three fractions always sum to 1.
double f3A5   = FM_CYP3A * ((CYPCL - 1.0) / CYPCL) * CYP;
double f3A4   = FM_CYP3A - f3A5;
double fm_sys = 1.0 - FM_CYP3A;

F_DEPOT_CSA = TVF_CSA;


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
