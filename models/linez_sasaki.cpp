// =========================================================================
// linez_sasaki.cpp -- Linezolid, adult Japanese inpatients (Sasaki 2011)
//
// SOURCE
//   Sasaki T, Takane H, Ogawa K, Isagawa S, Hirota T, Higuchi S, Horii T,
//   Otsubo K, Ieiri I. "Population Pharmacokinetic and Pharmacodynamic
//   Analysis of Linezolid and a Hematologic Side Effect, Thrombocytopenia,
//   in Japanese Patients." Antimicrob Agents Chemother 2011;55(5):1867-1873.
//   doi:10.1128/AAC.01185-10
//   PDF: popPK linezolide Sasaki.pdf
//   Read here with `pdftotext -layout` (/mingw64/bin/pdftotext, 513 lines,
//   all 7 pages = journal pp. 1867-1873). The final-model equations were read
//   off the Results paragraph on p.1870 verbatim; Table 2 was cross-checked
//   line by line. `pdftoppm` and `pdfinfo` are NOT installed on this machine,
//   so no page image was rendered by me -- the rendered-image confirmations
//   quoted in the upstream verification (Table 2 at 5-6x, the negative
//   Box-Cox signs, the "ka (1iter/h)" typo) are that analyst's work, not mine,
//   and I am relying on them for those three points only.
//
//   ORIGINAL MODEL AUTHOR = Sasaki. The PopPK model is original to this
//   paper. Only the fixed ka value is borrowed (Abe et al. 2009,
//   J Clin Pharmacol 49:1071-1078), which does not make it Abe's model.
//
// POPULATION
//   50 Japanese adult inpatients at Tottori University Hospital (Yonago),
//   prospective open-label, June 2006 - July 2010. 135 STEADY-STATE plasma
//   concentrations (HPLC, LLOQ 0.1 ug/mL), 1-5 samples per patient.
//     Age    69.1 y (SD 12.8; range 32-92)   -- ADULT, minimum 32 y
//     BW     57.3 kg (SD 12.1; range 38.4-100)
//     Sex    36 M / 14 F
//     SCr    1.10 mg/dL (SD 0.86; range 0.20-4.24) == approx 97 umol/L
//     CLcr   74.0 mL/min (SD 54.5; range 9.43-330), Cockcroft-Gault
//     Route  39 i.v. only / 8 p.o. only / 3 both; 300-600 mg twice daily;
//            i.v. given as a constant infusion "within 1 to 2 h"
//     Duration of therapy 11.1 days (SD 4.50; range 4-23)
//   Indication is described by the paper only as "Japanese patients with
//   infectious disease". NO pathogen breakdown is printed. NOT prosthetic
//   joint infection; no orthopaedic subgroup anywhere in the paper.
//   Reference weight 57.9 kg -- adult, nowhere near the 6.8 kg paediatric
//   trip-wire that has already bitten this project.
//
// STRUCTURE
//   ONE compartment, first-order absorption (ka FIXED, no IIV), linear
//   first-order elimination, parameterised as CL and V (not micro-constants).
//   NONMEM VI level 1.0, FOCE-INTER; PsN bootstrap 1000 resamples; case
//   deletion diagnostics (all parameters within 20% of original).
//   Deliberately LINEAR: the authors fitted Plock's autoinhibition model to
//   these data and REJECTED it (no OFV/GOF gain, poor convergence).
//
// PARAMETERS (Results p.1870, printed verbatim; Table 2 for RSE/bootstrap)
//   ka = 0.583 h^-1                   FIXED from Abe 2009, %RSE NE, no IIV
//   CL (L/h) = 2.85 * (CLCR/60.9)^0.618 * 0.472^CIR
//        2.85   %RSE 5.93   bootstrap 2.85 [2.55; 3.26]
//        0.618  %RSE 15.1   bootstrap 0.624 [0.408; 0.815]
//        0.472  %RSE 14.2   bootstrap 0.470 [0.333; 0.673]
//   V  (L)   = 33.6 * (BW/57.9)
//        33.6   %RSE 4.82   bootstrap 33.6 [30.7; 37.1]
//        The BW exponent is IMPLICIT 1 -- it is not an estimated parameter
//        and does not appear in Table 2. This is NOT allometric 0.75.
//   IIV, exponential, Table 2 column "IIV (%; %RSE)":
//        CL 35.2% (%RSE 30.6), eta-shrinkage 2.89%
//        V  30.8% (%RSE 35.8), eta-shrinkage 30.8%
//   Residual error: ADDITIVE ONLY, SD 1.43 ug/mL (%RSE 12.5), bootstrap
//        median 1.42, eps-shrinkage 30.5%.
//
//   OMEGA DERIVATION -- WHICH RULE WAS USED AND WHY
//     The paper states an exponential IIV model and prints a bare "IIV (%)".
//     It never disambiguates 100*sqrt(omega^2) from 100*sqrt(exp(omega^2)-1).
//     I used the FIRST reading: the printed % IS omega on the log scale, so
//     omega^2 = (%/100)^2. This is NOT the ln(1+CV^2) rule.
//       CL: omega^2 = 0.352^2 = 0.123904
//       V : omega^2 = 0.308^2 = 0.094864
//     Justification, not just convention: with omega = 0.352 on CL the
//     published Table 3 Monte-Carlo attainment percentages reproduce to
//     within ~2 points across all six CLCR/cirrhosis groups. The ln(1+CV^2)
//     reading would give 0.1174 / 0.0904 -- under 6% different in variance
//     and immaterial here, but the reading above is the one the paper's own
//     simulations are consistent with. NO off-diagonal is published; a
//     DIAGONAL omega is assumed.
//
//   SIGMA DERIVATION
//     Additive SD 1.43 ug/mL is on the natural concentration scale, so
//     $SIGMA holds its VARIANCE: 1.43^2 = 2.0449. There is NO proportional
//     component in the PK model. The 19.8% proportional / 15.7 x10^3/uL
//     additive pair in the LOWER half of Table 2 belongs to the PLATELET PD
//     model and must not be transplanted here.
//
// -------------------------------------------------------------------------
// LIMITS -- every reason this model could mislead in an adult PJI cohort
//
//   1. BIOAVAILABILITY WAS NEVER REPORTED. F appears nowhere in Table 2, was
//      neither estimated nor stated to be fixed. I grepped the full extracted
//      text for "bioavail"/"F1"; the only hits are the Introduction's
//      background sentence "LZD has 100% oral bioavailability" and reference
//      4 (Beringer 2005). Because 39/50 PK subjects were i.v.-only, CL = 2.85
//      L/h is a TOTAL clearance, not CL/F. TVF below is hard-coded to 1 as an
//      ASSUMPTION of mine, on the strength of that background sentence alone.
//      Any oral simulation from this file inherits that assumption. If F is
//      in truth below 1 in a post-operative PJI patient (ileus, opioids,
//      enteral feeding), oral AUC here is OPTIMISTIC.
//
//   2. CL CARRIES NO EXPLICIT WEIGHT TERM, BUT AUC IS NOT WEIGHT-INDEPENDENT.
//      This is the trap. Weight enters V explicitly and CL only through
//      Cockcroft-Gault CLcr, which is itself proportional to WT, so
//      CL scales as WT^0.618 in this implementation. Worked example, 70 y
//      male, SCr 1.0 mg/dL: 50 kg -> CLcr 48.6, CL 2.48 L/h, AUC24 484
//      mg.h/L at 1200 mg/day; 100 kg -> CLcr 97.2, CL 3.80 L/h, AUC24 315.
//      A 1.54-fold swing, not zero. The real extrapolation risk in a heavy
//      European PJI patient is UNDER-exposure, and it is entirely mediated by
//      a covariate equation the paper did not print (see LIMIT 3), not by a
//      fitted weight effect. If the harness ever supplies CLCR directly
//      instead of letting this file compute it, the weight dependence
//      DISAPPEARS -- that is a real behavioural difference between two
//      defensible wirings of the same published model.
//
//   3. THE COCKCROFT-GAULT VARIANT IS NOT SPECIFIED. The paper says only
//      "Creatinine clearance (CLCR) was calculated by the Cockcroft-Gault
//      equation (8)". It does NOT say whether the 0.85 female multiplier was
//      applied, and 14 of 50 development subjects were women. Applying it
//      swings CLcr by 15% in women, which propagates to about 9% on CL
//      through the 0.618 power. THIS FILE APPLIES IT (CGFEM = 0.85, exposed
//      as a parameter so it can be set to 1.0 without editing equations).
//      That choice is mine and unverifiable from the source. Total body
//      weight is used, not IBW or adjusted BW -- also unstated in the paper.
//      SEX is used ONLY here, and only with harness coding 0 = male,
//      1 = female. Sex is NOT a covariate of the published model, so there is
//      no sign convention in the paper to transpose; but note that sex WAS
//      GAM-screened onto CL and dropped at the NONMEM step, so it re-enters
//      implicitly and undocumented through this 0.85 factor.
//
//   4. UNIT TRAP ON CREATININE. The paper works in mg/dL (mean 1.10 mg/dL);
//      the harness supplies CREAT in umol/L. CREAT/88.4 is applied in $MAIN.
//      Getting this wrong scales CLcr by 88 and, through the 0.618 power,
//      changes CL by roughly 16-fold. The harness covariate is NOT renamed or
//      re-united; the conversion happens inside $MAIN only.
//
//   5. THE PROVENANCE OF THE 60.9 mL/min CENTERING CONSTANT IS UNKNOWN. The
//      paper prints 60.9 in the equation and reports the cohort MEAN CLcr as
//      74.0 mL/min (SD 54.5) in Table 1. It reports no median and never says
//      what 60.9 is. Do not "correct" it to 74.0. Likewise V is centred on
//      57.9 kg while the cohort mean BW is 57.3 kg -- use the printed 57.9.
//
//   6. WRONG POPULATION SHAPE FOR EUROPEAN PJI. Elderly (mean 69.1 y),
//      low-body-weight (mean 57.3 kg, only the very top of the 38.4-100 kg
//      range overlapping a typical European PJI cohort), Japanese, and
//      hospitalised with unspecified general infectious disease. Nothing
//      orthopaedic. Reference 30 (Rao & Hamilton, linezolid in Gram-positive
//      orthopaedic infection) is cited in the Discussion but contributes no
//      data to the fit.
//
//   7. EXPOSURE IS HIGH RELATIVE TO WESTERN LINEZOLID EXPERIENCE. CL = 2.85
//      L/h at the reference CLcr gives AUC24 = 421 mg.h/L at 600 mg q12h.
//      Western adults are usually quoted nearer 7 L/h -- but that ~7 L/h
//      figure is NOT in this paper and is NOT cited to any source here; it is
//      recollection, and must be labelled as an unverified cross-model
//      comparison, not a finding. Two honest qualifications: (a) a
//      normal-renal PJI patient (CLcr 100 mL/min) gets CL 3.87 L/h and AUC24
//      ~310 mg.h/L, which is unremarkable; (b) 1.4x-4.2x disagreement between
//      published models of one drug is the measured norm in this library, so
//      this is not disqualifying. The DIRECTION is what matters: against a
//      Western comparator this model will make AUC/MIC attainment look
//      optimistic and thrombocytopenia risk look pessimistic.
//
//   8. THE SAMPLING WINDOW CANNOT IDENTIFY A TERMINAL PHASE. Samples were
//      taken pre-dose and at 1.5, 2, 4, 5, 8 and 9 h only, at steady state
//      only (>= 3 days), 1-5 per patient, 135 concentrations from 50
//      patients. A one-compartment model was IMPOSED (the paper says so); no
//      two-compartment alternative is reported as having been tested. Any
//      distribution phase and any true terminal phase beyond ~9 h are
//      unidentified. The derived t1/2 of 8.2 h (ke = 2.85/33.6 = 0.0848 h^-1)
//      is a fitted mono-exponential slope, not an observed terminal slope.
//
//   9. SHRINKAGE. eta-shrinkage on V is 30.8% and eps-shrinkage 30.5% (both
//      moderately high), so IIV on V is plausibly UNDERSTATED and individual
//      V estimates are partly shrunk to the typical value. eta-shrinkage on
//      CL is fine (2.89%). MAP estimation of V from sparse PJI TDM will be
//      correspondingly weak.
//
//  10. ka IS FIXED AND ORAL PEAK SHAPE IS UNCHARACTERISED. ka = 0.583 h^-1
//      was imported from Abe 2009 because of "a lack of sampling points"
//      early after dosing, carries no IIV, and only 11 of 50 patients took
//      any oral dose. No absorption lag time is reported (assumed zero here).
//      Use this file for oral AUC and trough; DO NOT trust its oral Cmax or
//      Tmax. Table 2 prints the ka unit as "(1iter/h)" -- a journal
//      typesetting error; a first-order rate constant cannot be liter/h and
//      the value is h^-1.
//
//  11. THE MODEL IS LINEAR AND LINEZOLID IS NOT. Autoinhibition of linezolid
//      metabolism (Plock 2007) is real and the authors acknowledge it; they
//      simply could not fit it to sparse data. Meagher's alternative (linear
//      renal + Michaelis-Menten nonrenal CL) is discussed and not refitted.
//      Dose-escalation or heavy-accumulation simulations from this file are
//      extrapolating beyond a structure the authors knew was incomplete.
//
//  12. THE RENAL COVARIATE IS A STEADY-STATE PHENOMENON. The authors note
//      their CLcr effect contradicts Brier et al., who found none after a
//      SINGLE dose, and attribute the difference to saturation of the
//      nonrenal pathway at steady state (approximately 30% of linezolid is
//      eliminated unchanged renally). Do not use this model to predict
//      first-dose kinetics.
//
//  13. CIRRHOSIS: n = 4. The 0.472 multiplier on CL rests on four Child-Pugh
//      grade C patients (troughs 32.5, 36.4, 40.8, 45.4 ug/mL), bootstrap
//      [0.333; 0.673]. CIR is NOT a harness covariate (the harness has AGE,
//      WT, SEX, CREAT, CREAT2, HT, ALB, DIAL, IHD, DM, TEMP, FUS, PERIOD --
//      no cirrhosis flag). It is exposed below as a plain $PARAM switch
//      DEFAULTING TO 0, which is the safe default. Setting CIR = 1 roughly
//      doubles AUC and is supported by four patients.
//
//  14. RIFAMPIN-TREATED PATIENTS WERE EXCLUDED BY DESIGN. Methods: patients
//      "taking drugs known to interact with LZD (including rifampin) were
//      excluded". Rifampicin plus linezolid is a standard staphylococcal PJI
//      combination, this library already holds four rifampicin models, and
//      the paper itself cites evidence that rifampin lowers linezolid AUC.
//      Co-simulating this file alongside a rifampicin model reproduces
//      neither the interaction nor the population it was fitted in. This is
//      the single largest reason the model does not describe the commonest
//      real linezolid regimen in PJI.
//
//  15. NO BONE, SYNOVIAL OR TISSUE DATA WHATSOEVER. The word "bone" occurs in
//      this paper only as "bone marrow". It contributes nothing to the bone
//      penetration category. Bone exposure must come from elsewhere (e.g.
//      Stolle's linezolid microdialysis paper, already in Reference/).
//
//  16. NOTHING HERE SOURCES THE BIOFILM MIC MULTIPLIER. The paper's only PD
//      target is the PLANKTONIC AUC24/MIC >= 100 at MIC 2.0 mg/L from a
//      murine model (Andes 2002). The presence of a numeric PK/PD target in
//      this paper must NOT be read as support for the simulator's unsourced
//      8-256x biofilm MIC factor.
//
//  17. DURATION TOXICITY, WHICH PJI MAKES WORSE. The paper's companion
//      Friberg-type platelet model (NOT implemented in this file -- it is a
//      PD model, and its feedback term (Circ0/Circ)^gamma is not printed in
//      this paper) predicts >30% thrombocytopenia by DAY 14 at 1200 mg/day
//      when CLcr < 30 mL/min or Child-Pugh C. PJI courses run 6-12 weeks,
//      far past that horizon and far past the 4-23 day exposure actually
//      observed. Long-course simulations from this PK file are outside the
//      duration the source population ever experienced.
//
//  18. INTERNAL VALIDATION ONLY. The PK model was checked by bootstrap and
//      case deletion; it was never externally validated. Only the PLATELET
//      model was validated externally (60 further patients). Supplemental
//      Figs S1-S3 (GOF diagnostics) are not in this PDF and I did not fetch
//      them from the AAC website, so anything they say about misspecification
//      is unchecked here.
//
//  19. THE INFUSION DURATION IS A CHOICE. The dataset used 1-2 h and it is
//      not recoverable per record; only the paper's own Table 3 simulations
//      pin it at 1 h. Use 1 h in the dosing record for consistency with
//      Table 3. Nothing in this file forces it -- it comes from the data set.
//
// SMOKE EXPECTATION
//   Standard adult IV regimen 600 mg q12h, 1 h infusion, 80 kg / 70 y male /
//   CREAT 80 umol/L (CLcr 85.9 mL/min by the CG variant coded here):
//     CL 3.53 L/h, V 46.4 L, t1/2 9.1 h
//     steady-state AUC24 = 340 mg.h/L      -- IN the 200-400 target band
//     steady-state peak (end of infusion) ~21 mg/L
//     steady-state trough (12 h)          ~9.0 mg/L
//   THE TROUGH IS ABOVE THE 2-8 mg/L BAND GIVEN IN THE BRIEF, by roughly 1
//   mg/L, and this is the model, not a coding error: Sasaki's CL of 2.85 L/h
//   is low relative to Western linezolid experience (LIMIT 7), so troughs run
//   high while AUC24 still lands in range. Parameters have NOT been adjusted
//   to make the trough fit. A 100 kg 60 y male at CREAT 80 gives AUC24 273
//   mg.h/L, peak 16.7, trough 7.3 mg/L -- fully in band; a 50 kg 80 y female
//   at CREAT 90 gives AUC24 595 mg.h/L, trough 16 mg/L, which is the paper's
//   own message that low-weight elderly renally-impaired patients are
//   over-exposed at 1200 mg/day.
// =========================================================================

[SET] end=72, delta=0.1


$PARAM @annotated
TVCL   : 2.85    : Typical clearance at CLcr 60.9 mL/min, no cirrhosis (L/h)
TVV    : 33.6    : Typical volume of distribution at BW 57.9 kg (L)
KA     : 0.583   : Absorption rate constant (1/h) -- FIXED, Abe 2009, no IIV
CLCRREF: 60.9    : CLcr centering constant (mL/min) -- printed, provenance unknown
BWREF  : 57.9    : Body weight centering constant (kg) -- printed
PWRCLCR: 0.618   : Power exponent on CLcr for CL (dimensionless)
EFFCIR : 0.472   : Multiplier on CL in Child-Pugh grade C cirrhosis
// CIR is NOT a harness covariate and defaults to 0 -- see LIMITS 13
CIR    : 0       : Cirrhosis flag, 0 = none, 1 = Child-Pugh grade C
TVF    : 1.0     : Oral bioavailability -- ASSUMED 1, never reported, LIMITS 1
CGFEM  : 0.85    : Female factor in Cockcroft-Gault -- ASSUMED applied, LIMITS 3

ETA1 : 0 : IIV on clearance
ETA2 : 0 : IIV on volume of distribution


$PARAM @annotated @covariates
WT    : 57.9 : total body weight (kg)
AGE   : 69.1 : age (years)
SEX   : 0    : 0 = male, 1 = female
CREAT : 97.0 : serum creatinine (umol/L)


// Exponential IIV, Table 2 "IIV (%)". The printed % is read as omega on the
// LOG scale, so omega^2 = (%/100)^2 -- NOT ln(1+CV^2). See the OMEGA
// DERIVATION note in the header for why this reading was chosen.
//   ETA1 CL : 35.2% -> 0.352^2 = 0.123904   (eta-shrinkage 2.89%)
//   ETA2 V  : 30.8% -> 0.308^2 = 0.094864   (eta-shrinkage 30.8%, understated)
// No off-diagonal is published; diagonal assumed.
$OMEGA 0.123904 0.094864


// PUBLISHED residual error: ADDITIVE ONLY, SD 1.43 ug/mL -> variance 2.0449.
// This is NOT the house 0.01/0.0001 default -- the paper prints it. There is
// no proportional component. Do not import the platelet-model error pair.
$SIGMA 0 2.0449


$CMT @annotated
GUT  : Oral absorption compartment (mg) [ADM]
CENT : Central compartment (mg) [ADM, OBS]


$MAIN
// Cockcroft-Gault. CREAT arrives from the harness in umol/L; the paper works
// in mg/dL (1 mg/dL = 88.4 umol/L). The harness covariate is not renamed or
// re-united -- the conversion is local. See LIMITS 4.
double CRMGDL = CREAT / 88.4;
double CLCR   = (140.0 - AGE) * WT / (72.0 * CRMGDL) * ((SEX == 0) ? 1.0 : CGFEM);
if (CLCR < 1.0) CLCR = 1.0;    // guard: keep the power term finite

double CL = TVCL * pow(CLCR/CLCRREF, PWRCLCR) * pow(EFFCIR, CIR)
                 * exp(ETA1 + ETA(1));
double V  = TVV  * (WT/BWREF) * exp(ETA2 + ETA(2));

// Bioavailability stated explicitly rather than left implicit, because it was
// never published. F = 1 is an ASSUMPTION, not a fitted value. IV doses go
// into CENT and bypass this.
F_GUT = TVF;


$ODE
dxdt_GUT  = - KA * GUT;
dxdt_CENT =   KA * GUT - (CL / V) * CENT;


$TABLE
// Additive residual error only, SD 1.43 ug/mL (variance in $SIGMA).
double CP = CENT / V;
double DV = CP + EPS(2);
int i = 0;
while(DV<0 && i <100) {
simeps();
DV = CP + EPS(2);
++i;
}


$CAPTURE DV CP CL V CLCR
