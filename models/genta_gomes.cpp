// =========================================================================
// genta_gomes.cpp -- Gentamicin, adult infective endocarditis (Gomes 2017)
//
// SOURCE
//   Gomes A, van der Wijk L, Proost JH, Sinha B, Touw DJ. "Pharmacokinetic
//   modeling of gentamicin in treatment of infective endocarditis: model
//   development and validation of existing models." PLoS ONE 2017;12(5):
//   e0177324.
//   PDF: Pharmacokinetic modeling of gentamicin in
//        treatment of infective endocarditis ... Gomes.pdf
//
// POPULATION
//   65 adults with infective endocarditis, The Hague, 2011-2013, 221 serum
//   samples. Mean age 69.3 y (32-92), mean weight 76.2 kg (46-121).
//   Cockcroft-Gault CLcr mean 64.3 mL/min, range 8.7-157.5 -- a WIDE and
//   CONTINUOUS renal range, which is what makes this family usable where the
//   haemodialysis model (genta_franck) is not. DIALYSIS PATIENTS EXCLUDED.
//   Once-daily 3 mg/kg starting dose, combined with a beta-lactam.
//
// SHARED STRUCTURE (Gomes 2017 Table 3 -- three parameter sets, one structure)
//   One compartment, first-order elimination, IV infusion (~30 min).
//     CL = CLm * (BW/70) + fr * CLcr        fr is a DIMENSIONLESS fraction of
//                                           creatinine clearance
//     V  = Vd_coefficient * LBMc            LBMc = lean body mass corrected
//                                           for fat distribution
//   CLm is a non-renal floor, so CL never falls to zero as renal function
//   declines. That is a real advantage over models whose clearance is strictly
//   proportional to CLcr.
//
// -------------------------------------------------------------------------
// LIMITS SHARED BY ALL THREE MODELS OF THIS FAMILY
//
//   1. THE LEAN-BODY-MASS EQUATION IS NOT PRINTED IN THE SOURCE. Gomes cites
//      Chennavasin (ref 27) for LBM and Schwartz (ref 28) for the fat
//      correction but reproduces neither, so V cannot be computed from the
//      paper alone. Implemented here as Devine ideal body weight plus a 0.4
//      fat correction -- the MW\Pharm convention, and the one that reproduces
//      the ~69.6 kg cross-check for a 174 cm male implied by the cohort
//      description. THIS IS AN ASSUMPTION AND MUST BE CONFIRMED against the
//      two original references before the model is trusted quantitatively.
//      The correction factor is exposed as FATC so it can be changed without
//      editing the equations.
//
//   2. THE COCKCROFT-GAULT VARIANT IS NOT SPECIFIED. The paper tabulates both
//      raw and BSA-normalised CLcr, which differ by about 10% in its cohort,
//      and states neither the equation nor which form feeds CL. Standard
//      Cockcroft-Gault on TOTAL body weight with the 0.85 female factor is
//      used here.
//
//   3. UNIT TRAP ON CLcr. fr is dimensionless and CL is in L/h, so CLcr must
//      be converted to L/h. The paper reports it in mL/min. Using mL/min
//      directly would inflate clearance 16.7-fold.
//
//   4. THE THREE MODELS ARE NOT INDEPENDENT. They share this structure and
//      this covariate algebra, and the endocarditis model is a Bayesian
//      re-estimation SEEDED from the Evers ICU model. Holding one out and
//      training on another is a weak external validation. The structurally
//      independent partner in this library is `genta_debord`, which drives
//      elimination from a non-steady-state Jelliffe creatinine clearance
//      built on two successive creatinines.
//
//   5. NO PUBLISHED RESIDUAL ERROR MODEL. Fitting used the EMIT/Architect
//      assay error function SD = 0.0766 + 0.0006*C + 0.0064*C^2 (mg/L). That
//      polynomial is implemented directly in $TABLE rather than approximated
//      by a proportional term, because its quadratic term dominates at peak
//      concentrations and no proportional/additive pair matches it across the
//      0.5-20 mg/L range.
//
//   6. IIV IS REPORTED AS AN SD IN NATURAL UNITS, not as a variance and not
//      as a CV%. Table 3 columns are headed "mean" and "SD", and the paper
//      states the random effects are lognormal. $OMEGA below therefore holds
//      ln(1 + (SD/mean)^2); the arithmetic is shown per parameter.
//
//   MODEL-SPECIFIC NOTES
//     - CLm was FIXED to 0.277 with NO between-subject variability, because
//       fewer than 10% of the cohort had renal function below 15 mL/min and
//       the non-renal component could not be identified. $OMEGA for ETA1 is
//       therefore 0: that is the paper's estimate, not an omission here.
//     - fr and Vd carry nonparametric bootstrap 95% CIs (1000 replicates):
//       fr [0.610; 0.794], Vd [0.292; 0.331].
//     - Endocarditis, not bone infection. Transporting it to a PJI cohort is
//       an extrapolation, though a far shorter one than from dialysis or from
//       febrile neutropenia.
//
//   SMOKE EXPECTATION
//     80 kg / 175 cm / 65 y male, CREAT 80 umol/L, 400 mg over 30 min:
//     CL 4.17 L/h, V 23.2 L, peak 17.3 mg/L, t1/2 3.85 h, C24 0.23 mg/L.
// =========================================================================

[SET] end=72, delta=0.1


$PARAM @annotated
TVCLM : 0.277 : Non-renal clearance at 70 kg (L/h)
TVFR  : 0.698 : Renal clearance as a fraction of CLcr (dimensionless)
TVVD  : 0.312 : Volume coefficient (L/kg of corrected lean body mass)
FATC  : 0.4 : Fat correction on ideal body weight -- see LIMITS 1

ETA1 : 0 : IIV on non-renal clearance
ETA2 : 0 : IIV on the renal clearance fraction
ETA3 : 0 : IIV on the volume coefficient


$PARAM @annotated @covariates
WT    : 76.0 : total body weight (kg)
HT    : 174.0 : height (cm)
AGE   : 69.0 : age (years)
SEX   : 0 : 0 = male, 1 = female
CREAT : 86.0 : serum creatinine (umol/L)


// ETA1 CLm : SD 0     (fixed, no IIV)      -> omega^2 0
// ETA2 fr  : SD 0.358 / mean 0.698 = 51.3% -> ln(1+CV^2) = 0.2335
// ETA3 Vd  : SD 0.076 / mean 0.312 = 24.4% -> ln(1+CV^2) = 0.0576
$OMEGA 0.0000 0.2335 0.0576


$SIGMA 1   // unit variance; the assay error polynomial scales it in $TABLE


$CMT @annotated
CENT : Central compartment (mg) [ADM, OBS]


$MAIN
// Corrected lean body mass. Devine ideal body weight plus a fat correction --
// a SUBSTITUTE for the Chennavasin/Schwartz equations the paper omits.
double IBW = (SEX == 0) ? 50.0 + 2.3*(HT/2.54 - 60.0)
                        : 45.5 + 2.3*(HT/2.54 - 60.0);
if (IBW < 1.0) IBW = WT;              // guard: very short stature
double LBMC = IBW + FATC*(WT - IBW);
if (LBMC > WT) LBMC = WT;             // corrected lean mass cannot exceed total

// Cockcroft-Gault. CREAT arrives in umol/L; 1 mg/dL = 88.4 umol/L.
double CRMGDL = CREAT / 88.4;
double CCR_MLMIN = (140.0 - AGE) * WT / (72.0 * CRMGDL) * ((SEX == 0) ? 1.0 : 0.85);
double CCR = CCR_MLMIN * 0.06;        // mL/min -> L/h; fr is dimensionless

double CLM = TVCLM * exp(ETA1 + ETA(1));
double FR  = TVFR  * exp(ETA2 + ETA(2));
double CL  = CLM * (WT/70.0) + FR * CCR;
double V   = TVVD * exp(ETA3 + ETA(3)) * LBMC;


$ODE
dxdt_CENT = - (CL / V) * CENT;


$TABLE
// EMIT/Architect assay error function, implemented exactly rather than
// approximated: SD = 0.0766 + 0.0006*C + 0.0064*C^2 (mg/L). See LIMITS 5.
double CP = CENT / V;
double DV = CP + (0.0766 + 0.0006*CP + 0.0064*CP*CP) * EPS(1);
int i = 0;
while(DV<0 && i <100) {
simeps();
DV = CP + (0.0766 + 0.0006*CP + 0.0064*CP*CP) * EPS(1);
++i;
}


$CAPTURE DV CP CL V CCR LBMC
