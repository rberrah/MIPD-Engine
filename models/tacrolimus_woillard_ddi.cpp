// Standalone PK adaptation for TDM.
// The DDI mechanism placeholders were removed: this file describes baseline PK only.
// No patient data or pasted model is persisted by the application.
$PLUGIN tad
$SET end = 168, delta = 0.1

$PARAM @annotated
// --- TACROLIMUS BASE ---
TVCL_TAC : 21.2 : CL app Tacro (L/h)
TVV1_TAC : 486  : V1 app Tacro (L)
TVQ_TAC  : 79   : Q Tacro (L/h)
TVV2_TAC : 271  : V2 Tacro (L)
TVKTR_TAC: 3.34 : KTR Tacro (1/h)
HTCL  : -1.14 : Effet HT
STV1  : 0.29  : Effet ST
STKTR : 1.53  : Effet KTR
CYPCL : 2.00  : Ratio Genotype
FORME_OTH   : 1    : Non-inhibitable clearance: 0 = proportional, 1 = absolute
CL_OTHER_LH : 2.12 : Non-inhibitable clearance in L/h, used when FORME_OTH = 1
FRACT_OTH   : 0.10 : Non-inhibitable fraction of CL0, used when FORME_OTH = 0

// --- ETAs---
ETA1:0:CL
ETA2:0:V1
ETA3:0:Q
ETA4:0:V2
ETA5:0:KTR

$PARAM @annotated @covariates
HT    : 35    : Hematocrite
ST    : 1     : Galenique
CYP   : 0     : CYP Status

// --- PARAMETRES PERPETRATEUR & MECANISMES ---

$OMEGA 0.08 0.10 0.29 0.36 0.06
$SIGMA 0.012 0.1

$CMT @annotated
DEPOT_TAC : Tac depot [ADM]
TR1_TAC   : TR1
TR2_TAC   : TR2
TR3_TAC   : TR3
CENT_TAC  : Tac cent [OBS]
PERI_TAC  : Tac peri

$MAIN
// Initialisation des pools d'enzymes si nécessaire

// Calculs de base du Tacrolimus
double CL0 = TVCL_TAC * pow(HT/35.0, HTCL) * pow(CYPCL, CYP) * exp(ETA1 + ETA(1));
double V1  = TVV1_TAC * pow(STV1, ST) * exp(ETA2 + ETA(2));
double Q   = TVQ_TAC * exp(ETA3 + ETA(3));
double V2  = TVV2_TAC * exp(ETA4 + ETA(4));
double KTR = TVKTR_TAC * pow(STKTR, ST) * exp(ETA5 + ETA(5));

// Clearance that CYP3A inhibition cannot reach.
//
// Two forms, selected by FORME_OTH. In BOTH the individual baseline is preserved
// - CL_base = CL_other + CL_indiv - and the inhibitor scales only the second term.
//   FORME_OTH = 0  proportional: CL_other = FRACT_OTH x CL0, the same share in
//                  everyone, so the ceiling is 1 / FRACT_OTH for every patient.
//   FORME_OTH = 1  absolute: CL_other = CL_OTHER_LH, in L/h. The non-inhibitable
//                  share is then larger the lower the patient's clearance, which
//                  protects the slow metaboliser - the behaviour actually seen in
//                  patients, where slow metabolisers show a smaller interaction.
//
// The absolute form is the one externally validated on 122 tacrolimus
// concentrations from 21 voriconazole-treated recipients; 2.12 L/h is 10% of the
// 21.2 L/h population clearance of Woillard 2011, the value the source manuscript
// used. It is the default here so the application and the published evaluation
// run the same structure.
//
// Splitting rather than adding also fixes an arithmetic slip: writing
// CL0 x (...) + CL_other makes CL = CL0 + CL_other with no inhibitor present,
// which inflates the estimated baseline. Here the two terms partition CL0.
double CL_OTH = (FORME_OTH > 0.5) ? CL_OTHER_LH : (FRACT_OTH * CL0);
double CL_IND = (CL0 > CL_OTH + 0.01) ? (CL0 - CL_OTH) : 0.01;

// Share of the inhibitable clearance carried by each isoform. These sum to 1 by
// construction, for any CYPCL; the previous form subtracted the non-inhibitable
// fraction inside the split as well as outside it, and was exact only at
// CYPCL = 2 - which is the fixed value, so nothing was wrong in practice.
double f3A4 = 1.0;
double f3A5 = 0.0;
if(CYP == 1) {
    f3A4 = 1.0 / CYPCL;
    f3A5 = (CYPCL - 1.0) / CYPCL;
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
double CL_TAC_DDI = CL_IND * (f3A4 * ACT_CYP3A4 + f3A5 * ACT_CYP3A5) + CL_OTH;

dxdt_DEPOT_TAC = -KTR * DEPOT_TAC;
dxdt_TR1_TAC   = KTR * DEPOT_TAC - KTR * TR1_TAC;
dxdt_TR2_TAC   = KTR * TR1_TAC   - KTR * TR2_TAC;
dxdt_TR3_TAC   = KTR * TR2_TAC   - KTR * TR3_TAC;
dxdt_CENT_TAC  = KTR * TR3_TAC * F_BOOST - (CL_TAC_DDI + Q)*(CENT_TAC/V1) + Q*(PERI_TAC/V2);
dxdt_PERI_TAC  = Q*(CENT_TAC/V1) - Q*(PERI_TAC/V2);

$TABLE
double CONC = CENT_TAC / (V1/1000.0);
double DV = CONC * (1.0 + EPS(1)) + EPS(2);
double ratio = (CL_TAC_DDI / CL0) / F_BOOST;
double CL_OUT = CL_TAC_DDI;

$CAPTURE DV CONC ratio CL_TAC_DDI CP_PERP CL_OUT
