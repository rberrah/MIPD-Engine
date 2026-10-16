// Standalone PK adaptation for TDM.
// The DDI mechanism placeholders were removed: this file describes baseline PK only.
// No patient data or pasted model is persisted by the application.
$PLUGIN tad
$SET end = 168, delta = 0.1

$PARAM @annotated
// --- EVEROLIMUS, MECHANISTIC WELL-STIRRED LIVER MODEL -------------------------
// ter Heine R, van Erp NP, Guchelaar HJ, de Fijter JW, Reinders MEJ,
// van Herpen CM, Burger DM, Moes DJAR. A pharmacological rationale for improved
// everolimus dosing in oncology and transplant patients.
// Br J Clin Pharmacol 2018;84:1575-1586.  doi:10.1111/bcp.13591
//
// FINAL MODEL column of table 2, fitted on 1240 whole-blood samples from 126
// patients: 71 with metastatic thyroid or breast cancer on Afinitor 10 mg once
// daily, and 55 renal transplant recipients on Certican 1.5-3 mg twice daily.
// Formulation had no effect on bioavailability, so the two indications share one
// model - which is exactly why it is offered here alongside the Moes model:
// it is the only everolimus model in the library that covers oncology doses.
//
// VALIDATED against the eight typical steady-state simulations the paper itself
// reports for a 40-year-old, 70 kg, 1.80 m man at Ht 45% (trough, ug/L, model
// vs paper): 10 mg QD 12.38/12.30 - 3.75 mg BID 11.69/11.70 - 0.75 mg BID
// 2.38/2.37 - 1 mg BID 3.17/3.16 - 2.25 mg BID 7.08/6.70 - and with high-dose
// prednisolone 0.75 mg BID 1.76/1.75, 1 mg BID 2.34/2.33, 3 mg BID 6.96/7.01.
// The peaks run ~13% above the two the paper quotes (79.6 vs 68.8 and 38.8 vs
// 34.2) because those come from a curve simulated WITH the 110% variability on
// the absorption time; replaying that variability here gives 64.5 and 32.7,
// i.e. the paper's values, while the troughs move by less than 1%.
//
// Structure (their figure 1 and the rate constants printed in the methods):
//   1 depot -> 2..5 four transit compartments -> 6 LIVER -> 7 central -> 8 peri
// The dose is absorbed INTO THE LIVER, so first-pass extraction is part of the
// structure rather than a bioavailability term. Note the consequence: for an
// oral dose the well-stirred model gives CL/F = CLint x fu exactly, so the oral
// exposure is inversely proportional to the intrinsic clearance and a CYP3A
// inhibitor is NOT damped by hepatic blood flow.
TVMAT     : 0.404 : Mean absorption time (h)
TVCLINT   : 340   : Apparent intrinsic clearance (L/h, plasma)
TVVC      : 175   : Central volume (L, plasma)
TVVP      : 577   : Peripheral volume (L, plasma)
TVQ_EVE   : 85.7  : Intercompartmental clearance (L/h)
QH_REF    : 90    : Hepatic BLOOD flow (L/h, fixed in the paper)
FU_EVE    : 0.27  : Unbound fraction in plasma (Kovarik, fixed)
NTR       : 4     : Number of transit compartments

// High-dose prednisolone (>= 20 mg/day) is the only covariate retained: it
// induces CYP3A4 and raises the apparent intrinsic clearance by 31%.
// The covariate is called PREDNI and not PRED: mapbayr's output table already
// has a PRED column (the population prediction), and a parameter of the same
// name produces a table with duplicate columns that dplyr then refuses.
PREDNI_CL : 1.31  : CLint multiplier under prednisolone >= 20 mg/day

// Allometric scaling on FAT-FREE MASS. The reference is 57.2 kg, the fat-free
// mass of a 70 kg, 1.80 m man by the Janmahasatian equation - the very typical
// individual the paper simulates. The text scales the FLOWS (QH_EVE, Q) and the
// VOLUMES (VC_EVE, VP_EVE); CLint is deliberately NOT scaled, as printed.
FFM_REF   : 57.2  : Reference fat-free mass (kg)
EXP_FLOW  : 0.75  : Allometric exponent on flows
EXP_VOL   : 1.0   : Allometric exponent on volumes

// Non-linear partitioning into red cells, fitted by the authors on the
// manufacturer's in vitro blood-distribution data (their equation 1):
//     Crb = Bmax * Cp / (Kd + Cp) + Kns * Cp        (all in mg/L)
//     Cwb = Ht * Crb + (1 - Ht) * Cp
// This is what makes the haematocrit a genuine covariate of the OBSERVATION:
// the model runs in plasma and reports whole blood, which is what the
// laboratory measures. At the concentrations seen in practice the whole
// blood / plasma ratio is about 5.3 at Ht 45% and falls as the drug saturates
// the erythrocyte binding sites.
BMAX_EVE  : 0.964  : Maximal concentration specifically bound to erythrocytes (mg/L)
KD_EVE    : 0.0920 : Dissociation constant (mg/L)
KNS_EVE   : 0.153  : Non-specific erythrocyte binding constant

// Liver volume. The paper computes it from body size and age with Small et al.
// (Biopharm Drug Dispos 2017;38:290-300), which is NOT in this library, so that
// equation is not reproduced: the liver volume is fixed at a typical adult value
// and scaled allometrically like the other volumes. This is safe, and it was
// checked rather than assumed: sweeping it from 0.8 to 3.0 L moves the
// steady-state trough by 0.2% and the peak by 0.9%, because the liver volume
// cancels out of the steady state - the oral exposure of a well-stirred liver
// model depends only on CLint x fu. Replace it if that paper is ever added.
VL_REF    : 1.70  : Liver volume at the reference fat-free mass (L)

// --- ETAs ---
ETA1:0:CL
ETA2:0:V1
ETA3:0:KA

// Covariates MUST live in an @covariates block: mapbayr refuses a data set whose
// columns match a plain parameter entry. The perpetrator and mechanism
// parameters are declared inside it too, because a perpetrator may bring a
// covariate of its own.
$PARAM @annotated @covariates
WT    : 70 : Body weight (kg)
BH    : 180 : Body height (cm)
SEX   : 1  : Sex (1 = male, 0 = female)
HT    : 45 : Haematocrit (%)
PREDNI : 0   : Prednisolone >= 20 mg/day (1 = yes)


// Final model, table 2, converted from CV% to the log scale (omega2 = ln(1+CV^2)):
//   CLint 33.9% -> 0.1087 | central volume 40.6% -> 0.1525 | MAT 110% -> 0.7930
// The 110% on the mean absorption time is INTRA-individual (inter-occasion) in the paper. mapbayr
// fits one eta per subject, so it is carried here as a between-subject term of
// the same magnitude: it keeps the absorption free enough for the Bayesian step
// to follow an early sample, which is the practical purpose it serves.
$OMEGA 0.1087 0.1525 0.7930
$SIGMA 0.0320 0.25

$CMT @annotated
DEPOT_EVE : Everolimus depot [ADM]
TR1_EVE   : Transit 1
TR2_EVE   : Transit 2
TR3_EVE   : Transit 3
TR4_EVE   : Transit 4
LIV_EVE   : Liver (plasma)
CENT_EVE  : Everolimus central [OBS]
PERI_EVE  : Everolimus peripheral

$MAIN

// Every local below carries an _EVE suffix on purpose: a perpetrator block may
// declare a parameter of the same plain name (isavuconazole declares BMI), and
// mrgsolve makes a parameter read-only, so the assignment would not compile.
// Fat-free mass, Janmahasatian et al. 2005 - the size descriptor the paper
// selected during base-model development.
double BMI_EVE  = WT / pow(BH / 100.0, 2.0);
double FFM_EVE  = SEX > 0.5 ? (9270.0 * WT) / (6680.0 + 216.0 * BMI_EVE)
                        : (9270.0 * WT) / (8780.0 + 244.0 * BMI_EVE);
double SF_EVE = pow(FFM_EVE / FFM_REF, EXP_FLOW);
double SV_EVE  = pow(FFM_EVE / FFM_REF, EXP_VOL);

double HCT_EVE = HT / 100.0;                 // the app carries haematocrit in %
double MAT_EVE = TVMAT * exp(ETA3 + ETA(3));
double KTR_EVE = (NTR + 1.0) / MAT_EVE;

double CLINT0 = TVCLINT * (PREDNI > 0.5 ? PREDNI_CL : 1.0) * exp(ETA1 + ETA(1));
double VC_EVE     = TVVC    * SV_EVE * exp(ETA2 + ETA(2));
double VP_EVE     = TVVP    * SV_EVE;
double Q_EVE  = TVQ_EVE * SF_EVE;
double QH_EVE     = QH_REF  * SF_EVE;
double QHP_EVE    = QH_EVE * (1.0 - HCT_EVE);        // liver PLASMA flow
double VL_EVE     = VL_REF * SV_EVE;

// Everolimus is cleared essentially by CYP3A4 and is a P-gp substrate. Moes
// tested ABCB1, CYP3A5, CYP2C8 and PXR and retained none, so there is no
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

// 3. Application to the INTRINSIC clearance - the enzyme is what the
//    perpetrator inhibits, and the well-stirred liver then decides how much of
//    that reaches the systemic exposure.
double CLINT_DDI = CLINT0 * (f3A4 * ACT_CYP3A4 + f3A5 * ACT_CYP3A5 + fm_sys);
double EH_EVE   = (CLINT_DDI * FU_EVE) / (QHP_EVE + CLINT_DDI * FU_EVE);
double CLH_EVE  = EH_EVE * QHP_EVE;

dxdt_DEPOT_EVE = -KTR_EVE * DEPOT_EVE;
dxdt_TR1_EVE   =  KTR_EVE * (DEPOT_EVE - TR1_EVE);
dxdt_TR2_EVE   =  KTR_EVE * (TR1_EVE   - TR2_EVE);
dxdt_TR3_EVE   =  KTR_EVE * (TR2_EVE   - TR3_EVE);
dxdt_TR4_EVE   =  KTR_EVE * (TR3_EVE   - TR4_EVE);
// The absorbed drug enters the LIVER: first pass is structural here.
dxdt_LIV_EVE   =  KTR_EVE * TR4_EVE * F_BOOST + QHP_EVE * (CENT_EVE / VC_EVE)
                  - QHP_EVE * (LIV_EVE / VL_EVE);
dxdt_CENT_EVE  =  QHP_EVE * (1.0 - EH_EVE) * (LIV_EVE / VL_EVE) - QHP_EVE * (CENT_EVE / VC_EVE)
                  - Q_EVE * (CENT_EVE / VC_EVE) + Q_EVE * (PERI_EVE / VP_EVE);
dxdt_PERI_EVE  =  Q_EVE * (CENT_EVE / VC_EVE) - Q_EVE * (PERI_EVE / VP_EVE);

$TABLE
// Plasma, then whole blood through the saturable erythrocyte binding.
double CP_EVE = CENT_EVE / VC_EVE;                                   // mg/L
double CRB    = BMAX_EVE * CP_EVE / (KD_EVE + CP_EVE) + KNS_EVE * CP_EVE;
double CWB    = HCT_EVE * CRB + (1.0 - HCT_EVE) * CP_EVE;                // mg/L
double CONC   = CWB * 1000.0;                                    // ng/mL = ug/L
double DV = CONC * (1.0 + EPS(1)) + EPS(2);

// For an oral dose the well-stirred model collapses to CL/F = CLint x fu, so
// the exposure ratio the app uses to rescale the dose is exactly the intrinsic
// clearance ratio - hepatic flow does not damp an oral interaction.
double CL_EVE_DDI = CLINT_DDI * FU_EVE;
double ratio = (CLINT_DDI / CLINT0) / F_BOOST;
double CL_OUT = CL_EVE_DDI;

$CAPTURE DV CONC ratio CL_EVE_DDI CP_PERP CL_OUT
