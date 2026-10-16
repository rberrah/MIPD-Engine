// Standalone PK adaptation for TDM.
// The DDI mechanism placeholders were removed: this file describes baseline PK only.
// No patient data or pasted model is persisted by the application.
$PLUGIN tad
$SET end = 168, delta = 0.1

$PARAM @annotated
// --- SIROLIMUS BASE (Djebli et al., Clin Pharmacokinet 2006;45:1135-48) -------
// Kidney transplant recipients on sirolimus + mycophenolate, no calcineurin
// inhibitor. Two-compartment model with an Erlang absorption chain.
//   CL/F = 14.1 + 14.2 x CYP3A5   (14.1 L/h for *3/*3, 28.3 L/h for *1 carriers)
//   Vc/F 218 L | Vp/F 292 L | Q/F 38.7 L/h | ktr 5.25 1/h
TVCL_SIR  : 14.1  : CL/F in CYP3A5 non-expressers (L/h)
TVCLC_SIR : 14.2  : Additional CL/F per CYP3A5 expresser status (L/h)
TVV1_SIR  : 218   : V1/F (L)
TVQ_SIR   : 38.7  : Q/F (L/h)
TVV2_SIR  : 292   : V2/F (L)
TVKTR_SIR : 5.25  : Transit rate constant (1/h)

// --- ETAs ---
ETA1:0:CL
ETA2:0:V1
ETA3:0:Q
ETA4:0:V2
ETA5:0:KTR

// Covariates MUST live in an @covariates block: mapbayr refuses a data set whose
// columns match a plain $PARAM entry ("Variables found both in the model and in
// the data"). The perpetrator and mechanism parameters are declared here too,
// because a perpetrator may contribute a covariate of its own - voriconazole
// brings CRP.
$PARAM @annotated @covariates
CYP   : 0     : CYP3A5 status (0 = non-expresser, 1 = expresser)
WT    : 70    : Body weight (kg)


// Inter-individual variability. Djebli reports 49.3% CV on CL/F, which is
// 0.219 on the variance scale; the remaining variances are pragmatic values of
// the same order, the paper reporting no usable estimate for them.
$OMEGA 0.219 0.15 0.25 0.25 0.20
$SIGMA 0.04 1.0

$CMT @annotated
DEPOT_SIR : Sirolimus depot [ADM]
TR1_SIR   : Transit 1
TR2_SIR   : Transit 2
TR3_SIR   : Transit 3
CENT_SIR  : Sirolimus central [OBS]
PERI_SIR  : Sirolimus peripheral

$MAIN

// Djebli's CYP3A5 effect is ADDITIVE on clearance, so it is written exactly as
// published rather than converted into a ratio.
double CL_NS  = TVCL_SIR;               // CYP3A4 + non-CYP3A routes
double CL_3A5 = TVCLC_SIR * CYP;        // CYP3A5 route, expressers only
double CL0 = (CL_NS + CL_3A5) * exp(ETA1 + ETA(1));
double V1  = TVV1_SIR * exp(ETA2 + ETA(2));
double Q   = TVQ_SIR  * exp(ETA3 + ETA(3));
double V2  = TVV2_SIR * exp(ETA4 + ETA(4));
double KTR = TVKTR_SIR * exp(ETA5 + ETA(5));

// Metabolic fractions used by the interaction machinery. Sirolimus is cleared
// almost entirely by CYP3A; a small non-CYP3A fraction is kept so that a strong
// inhibitor cannot drive clearance to zero.
double fm_sys = 0.05;
double f3A5 = CL_3A5 / (CL_NS + CL_3A5);
double f3A4 = 1.0 - f3A5 - fm_sys;


$ODE
// Standardised communication variables (default = no effect)
double ACT_CYP3A4 = 1.0;
double ACT_CYP3A5 = 1.0;
double ACT_PGP    = 1.0;
double F_BOOST    = 1.0;
double CP_PERP    = 0.0;

// 1. Perpetrator pharmacokinetics

// 2. Interaction mechanisms (they modify the ACT_ variables)

// 3. Application to sirolimus clearance
double CL_SIR_DDI = CL0 * (f3A4 * ACT_CYP3A4 + f3A5 * ACT_CYP3A5 + fm_sys);

dxdt_DEPOT_SIR = -KTR * DEPOT_SIR;
dxdt_TR1_SIR   =  KTR * DEPOT_SIR - KTR * TR1_SIR;
dxdt_TR2_SIR   =  KTR * TR1_SIR   - KTR * TR2_SIR;
dxdt_TR3_SIR   =  KTR * TR2_SIR   - KTR * TR3_SIR;
dxdt_CENT_SIR  =  KTR * TR3_SIR * F_BOOST - (CL_SIR_DDI + Q) * (CENT_SIR / V1)
                  + Q * (PERI_SIR / V2);
dxdt_PERI_SIR  =  Q * (CENT_SIR / V1) - Q * (PERI_SIR / V2);

$TABLE
double CONC = CENT_SIR / (V1 / 1000.0);      // ng/mL
double DV = CONC * (1.0 + EPS(1)) + EPS(2);
double ratio = (CL_SIR_DDI / CL0) / F_BOOST;
double CL_OUT = CL_SIR_DDI;

$CAPTURE DV CONC ratio CL_SIR_DDI CP_PERP CL_OUT
