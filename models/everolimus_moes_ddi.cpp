// Standalone PK adaptation for TDM.
// The DDI mechanism placeholders were removed: this file describes baseline PK only.
// No patient data or pasted model is persisted by the application.
$PLUGIN tad
$SET end = 168, delta = 0.1

$PARAM @annotated
// --- EVEROLIMUS BASE ----------------------------------------------------------
// Moes DJAR, Press RR, den Hartigh J, van der Straaten T, de Fijter JW,
// Guchelaar HJ. Population pharmacokinetics and pharmacogenetics of everolimus
// in renal transplant patients. Clin Pharmacokinet 2012;51:467-480.
//
// Two-compartment model with first-order absorption and a lag time, fitted on
// 783 samples from 53 renal transplant patients switched to a calcineurin-
// inhibitor-free regimen. Every value below is the FINAL MODEL column of the
// paper's table V (the abstract quotes a mix of base and final estimates).
TVCL_EVE  : 17.9  : CL/F (L/h)
TVV1_EVE  : 148   : V1/F (L)
TVQ_EVE   : 55.7  : Q/F (L/h)
TVV2_EVE  : 498   : V2/F (L)
TVKA_EVE  : 7.55  : Absorption rate constant (1/h)
TVALAG    : 0.709 : Absorption lag time (h)

// Ideal body weight on V1/F, standardised to the median of the dataset:
//     V1/F = theta * (65.75 / IBW)^(-1.41)
// The exponent is NEGATIVE, so the central volume RISES with ideal body weight -
// which the authors explain by everolimus being >75% partitioned into red cells
// and largely protein-bound, with length and sex entering the IBW formula.
// This is the only covariate retained: none of the ABCB1, CYP3A5, CYP2C8 or PXR
// polymorphisms had a significant effect, hence no CYP3A5 route in this model.
IBW_REF   : 65.75 : Median ideal body weight of the source population (kg)
IBW_V1    : -1.41 : Exponent of (IBW_REF / IBW) on V1/F

// The paper also reports CL/F = theta * (DOSE/2.25)^0.532, but states that the
// relationship "appeared to be caused by strict TDM": patients with a high
// clearance are moved to a higher dose, which creates an apparent
// dose-clearance link. Reproducing it here would make a simulated dose increase
// raise the clearance - exactly the artefact - so it is deliberately left out
// and the clearance is dose-independent.

// --- ETAs ---
ETA1:0:CL
ETA2:0:V1
ETA3:0:Q
ETA4:0:V2
ETA5:0:KA

// Covariates MUST live in an @covariates block: mapbayr refuses a data set whose
// columns match a plain $PARAM entry. The perpetrator and mechanism parameters
// are declared here too, because a perpetrator may bring a covariate of its own.
$PARAM @annotated @covariates
IBW   : 65.75 : Ideal body weight (kg)
WT    : 70    : Body weight (kg)


// Inter-individual variability, final model of table V, converted from CV% to
// the log scale as omega2 = ln(1 + CV^2):
//   CL/F 26.2% -> 0.066 | V1/F 27.7% -> 0.074 | ka 108.6% -> 0.778
// Moes estimated no random effect on Q or V2; the values kept here are
// pragmatic and only serve to let the Bayesian step move them slightly.
$OMEGA 0.066 0.074 0.20 0.20 0.778
$SIGMA 0.04 0.25

$CMT @annotated
DEPOT_EVE : Everolimus depot [ADM]
CENT_EVE  : Everolimus central [OBS]
PERI_EVE  : Everolimus peripheral

$MAIN

double CL0 = TVCL_EVE * exp(ETA1 + ETA(1));
double V1  = TVV1_EVE * pow(IBW_REF / IBW, IBW_V1) * exp(ETA2 + ETA(2));
double Q   = TVQ_EVE  * exp(ETA3 + ETA(3));
double V2  = TVV2_EVE * exp(ETA4 + ETA(4));
double KA  = TVKA_EVE * exp(ETA5 + ETA(5));

ALAG_DEPOT_EVE = TVALAG;

// Everolimus is cleared essentially by CYP3A4 and is a P-gp substrate. The
// source model found NO significant effect of CYP3A5 genotype, so there is no
// CYP3A5 route here - unlike tacrolimus, sirolimus or ciclosporine.
double fm_sys = 0.10;
double f3A4 = 1.0 - fm_sys;
double f3A5 = 0.0;


$ODE
// Standardised communication variables (default = no effect)
double ACT_CYP3A4 = 1.0;
double ACT_CYP3A5 = 1.0;
double ACT_PGP    = 1.0;
double F_BOOST    = 1.0;
double CP_PERP    = 0.0;

// 1. Perpetrator pharmacokinetics

// 2. Interaction mechanisms (they modify the ACT_ variables)

// 3. Application to everolimus clearance
double CL_EVE_DDI = CL0 * (f3A4 * ACT_CYP3A4 + f3A5 * ACT_CYP3A5 + fm_sys);

dxdt_DEPOT_EVE = -KA * DEPOT_EVE;
dxdt_CENT_EVE  =  KA * DEPOT_EVE * F_BOOST - (CL_EVE_DDI + Q) * (CENT_EVE / V1)
                  + Q * (PERI_EVE / V2);
dxdt_PERI_EVE  =  Q * (CENT_EVE / V1) - Q * (PERI_EVE / V2);

$TABLE
double CONC = CENT_EVE / (V1 / 1000.0);      // ng/mL
double DV = CONC * (1.0 + EPS(1)) + EPS(2);
double ratio = (CL_EVE_DDI / CL0) / F_BOOST;
double CL_OUT = CL_EVE_DDI;

$CAPTURE DV CONC ratio CL_EVE_DDI CP_PERP CL_OUT
