// Standalone PK adaptation for TDM.
// The DDI mechanism placeholders were removed: this file describes baseline PK only.
// No patient data or pasted model is persisted by the application.
$PLUGIN tad

$SET end = 168, delta = 0.1

$PARAM @annotated
// --- TACROLIMUS BASE ---
TVCL  : 19.6 : Typical value of clearance (L/h at 70 kg)
TVV1  : 123  : Typical central volume of distribution (L per 70 kg)
TVQ   : 74.9 : Typical intercompartmental clearance (L/h)
TVV2  : 500  : Typical peripheral volume of distribution (Fixed) (L per 70 kg)
TVKTR : 0.752 : Typical transit rate constant (1/h)
LAG   : 2.29  : Lag time (h)

// --- EFFET DES COVARIABLES ---
CYPCL : 1.625 : effect of CYP on cl

// --- ETAs---
ETA1  : 0 : ETA on clearance
ETA2  : 0 : ETA on V1
ETA3  : 0 : ETA on Q
ETA4  : 0 : ETA on V2
ETA5  : 0 : ETA on KTR

// Covariates MUST live in an @covariates block: mapbayr refuses a data set whose
// columns match a plain $PARAM entry. The perpetrator and mechanism parameters
// belong here too, since a perpetrator may bring a covariate of its own (CRP).
// This model has no haematocrit term - that is intentional, and the application
// warns the user that Envarsus ignores the haematocrit.
$PARAM @annotated @covariates
CYP  :    0 : expressor (1) non exp(0)


$OMEGA 0.08 0.10 0.29 0.36 0.06
$SIGMA 0.208 0.307

$CMT @annotated
DEPOT_TAC   : Tac depot [ADM]
TRANS1_TAC  : Transit compartment 1 (mg)
TRANS2_TAC  : Transit compartment 2 (mg)
CENT_TAC    : Tac cent [OBS]
PERI_TAC    : Tac peri

// --- COMPARTIMENTS INTERACTION ---

$MAIN
// Initialisation des pools d'enzymes si nécessaire

// Calculs Tacrolimus
double CL0 = TVCL * pow(CYPCL,CYP) * exp(ETA1 + ETA(1));
double V1  = TVV1 * exp(ETA2 + ETA(2));
double Q   = TVQ * exp(ETA3 + ETA(3));
double V2  = TVV2 * exp(ETA4 + ETA(4));
double KTR = TVKTR * exp(ETA5 + ETA(5));
double ALAG = LAG;

ALAG_DEPOT_TAC = ALAG;

// Fractions métaboliques du Tacrolimus
double fm_sys = 0.05;
double f3A4 = 0.95;
double f3A5 = 0.0;
if(CYP == 1) {
    f3A4 = (1.0 - fm_sys) / CYPCL;
    f3A5 = (CYPCL - 1.0 - fm_sys) / CYPCL;
}


$ODE
// Variables de communication standardisées (Par défaut = Pas d'effet)
double ACT_CYP3A4 = 1.0;
double ACT_CYP3A5 = 1.0;
double ACT_PGP    = 1.0;
double F_BOOST    = 1.0;
double CP_PERP    = 0.0;

// 1. Calcul de la Pharmacocinétique du Perpétrateur

// 2. Calcul des Mécanismes d'interaction (Modifient les ACT_)

// 3. Application à la clairance du Tacrolimus
double CL_TAC_DDI = CL0 * (f3A4 * ACT_CYP3A4 + f3A5 * ACT_CYP3A5 + fm_sys);


// ODE TACROLIMUS (Standardisée)
dxdt_DEPOT_TAC    = -KTR * DEPOT_TAC;
dxdt_TRANS1_TAC   = KTR * DEPOT_TAC - KTR * TRANS1_TAC;
dxdt_TRANS2_TAC   = KTR * TRANS1_TAC   - KTR * TRANS2_TAC;

// Application CL_TAC_DDI et F_BOOST
dxdt_CENT_TAC  = KTR * TRANS2_TAC * F_BOOST - (CL_TAC_DDI + Q)*(CENT_TAC/V1) + Q*(PERI_TAC/V2);
dxdt_PERI_TAC  = Q*(CENT_TAC/V1) - Q*(PERI_TAC/V2);

$TABLE
double CONC = CENT_TAC / (V1/1000.0);
double DV = CONC * (1.0 + EPS(1)) + EPS(2);
double ratio =  (CL_TAC_DDI / CL0) / F_BOOST;
double CL_OUT = CL_TAC_DDI;
double V1_OUT = V1;
double F_OUT = F_BOOST;

$CAPTURE DV CONC ratio CL_OUT V1_OUT F_OUT CP_PERP