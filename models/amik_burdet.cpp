// =========================================================================
// amik_burdet.cpp -- Amikacin, adult ICU (Burdet 2015, BASIC covariate-free
//                    model, 60 patients)
//
// SOURCE
//   Burdet C, Pajot O, Couffignal C, Armand-Lefevre L, Foucrier A,
//   Laouenan C, Wolff M, Massias L, Mentre F. "Population pharmacokinetics
//   of single-dose amikacin in critically ill patients with suspected
//   ventilator-associated pneumonia." Eur J Clin Pharmacol 2015;71(1):75-83.
//   doi:10.1007/s00228-014-1766-y. IMPACT trial, NCT00950222.
//   PDF: popPK single dose amikacin Burdet.pdf
//   Parameters below are Table 2, "Basic model (60 patients)" column.
//
// POPULATION
//   60 critically ill ADULTS, three French ICUs (Bichat, Paris; V. Dupouy,
//   Argenteuil), 2009-2011. Mechanically ventilated >48 h with clinically
//   suspected Gram-negative-bacilli ventilator-associated pneumonia; all
//   received concomitant imipenem. Median age 61.5 y (28-84), median total
//   body weight 78 kg (45-126), 47 male / 13 female. Under-18s were an
//   explicit exclusion criterion, as were patients needing renal
//   replacement therapy and patients already on imipenem or amikacin.
//   Median 4-h MEASURED urinary creatinine clearance 82 mL/min (4-412).
//   Median albumin 19 g/L (10-44); SOFA 7; SAPS-II 42; septic shock 43%.
//   291 concentrations, median 5 per patient (3-5), at 0.5, 1, 8, 16, 24 h
//   after a SINGLE 30-min infusion of a median 20 mg/kg (11-28).
//   Estimation by SAEM in Monolix 4.2; model choice by BIC.
//
// STRUCTURE
//   Two compartments, linear, first-order elimination from the central
//   compartment, drug in as a zero-order IV infusion (30 min in the study).
//   Parameterised CL / V1 / Q / V2. A one-compartment model was tested and
//   rejected by BIC. CL is a TRUE systemic clearance, not CL/F: the route is
//   IV only and there are no oral data. Units are L/h and L, self-checked
//   twice by the paper ("4.3 L/h (72 mL/min)"; V1+V2 = 37.3 L against
//   "approximately 37 L"). There is NO allometric or per-70-kg scaling.
//
// PARAMETERS (Table 2, basic model; RSE in brackets)
//   CL 4.0 L/h [7%]   V1 15.3 L [5%]   Q 12.2 L/h [6%]   V2 22.1 L [7%]
//   Derived: Vss 37.4 L, t1/2-alpha 0.46 h, t1/2-beta 7.28 h.
//
//   OMEGA DERIVATION. Methods p.77: "We used an exponential random effects
//   model for each pharmacokinetic parameter. We assumed the random effects
//   to have a normal distribution with a mean of 0 and a variance of
//   omega^2." So Table 2's omega column is the SD OF THE LOG-SCALE RANDOM
//   EFFECT. An SD already on the log scale is SQUARED (house rule case 3):
//     omega_CL 0.6 -> 0.36    omega_V1 0.3 -> 0.09
//     omega_Q  0.3 -> 0.09    omega_V2 0.5 -> 0.25
//   ln(1+CV^2) is NOT used here and would be wrong -- see LIMITS 7.
//
//   Off-diagonals from the printed correlations (basic model column):
//     rho(CL,V1) 0.3  rho(CL,Q) -0.4  rho(V1,Q) 0.4
//     rho(CL,V2) 0.3  rho(V1,V2) 0.2  rho(Q,V2) 0.2
//   cov(i,j) = rho * omega_i * omega_j, e.g. cov(CL,Q) = -0.4*0.6*0.3
//   = -0.072. See LIMITS 8 for why the block is implemented rather than
//   discarded, and for the positive-definiteness check.
//
//   RESIDUAL ERROR (published): a = 0.2 mg/L [17%] additive SD,
//   b = 0.1 i.e. 10% [12%] proportional SD. $SIGMA holds variances:
//   b^2 = 0.01 and a^2 = 0.04. NOT the house default -- these are the
//   paper's own numbers, which happen to coincide with the default
//   proportional term.
//
// -------------------------------------------------------------------------
// LIMITS -- every reason this model could mislead in an adult PJI cohort
//
//   1. THE SOURCE PHYSIOLOGY IS NOT PJI PHYSIOLOGY, AND THE VOLUME IS
//      INFLATED BY DESIGN. Every patient was mechanically ventilated with
//      suspected VAP; 43% were in septic shock and median albumin was
//      19 g/L, i.e. profoundly hypoalbuminaemic and third-spacing. Vss is
//      37.4 L = 0.48 L/kg at the 78 kg reference, against roughly
//      0.25-0.30 L/kg for a stable ward adult. The authors state the
//      finding plainly: clearance "decreased" and volume "increased" versus
//      published norms. Applied to an ambulatory or ward PJI patient this
//      model will UNDER-PREDICT the peak, by roughly the ratio of the true
//      volume to 37 L -- potentially by a third or more. The paper concedes
//      exactly this: "the results may not apply to patients who do not
//      require mechanical ventilation."
//
//   2. IT IS COVARIATE-FREE BY DELIBERATE CHOICE, SO EVERY PATIENT LOOKS
//      IDENTICAL A PRIORI. The published FINAL model has covariates, but
//      they cannot be supplied by this harness (see the companion note in
//      `not_drafted`). The consequence of using the basic model instead is
//      that a 45 kg anuric patient and a 126 kg patient with augmented
//      renal clearance receive the SAME typical CL 4.0 L/h and V1 15.3 L,
//      and all of that real, explainable difference is dumped into the
//      random effects. omega_CL is 0.6, i.e. a 60% spread on clearance;
//      the covariate model would have halved it (the paper reports "a 50 %
//      maximal decrease for amikacin CL"). THIS MODEL IS THEREFORE ONLY
//      HONEST AS A POPULATION PRIOR FOR MAP/TDM FEEDBACK. It should not be
//      used for a priori individualised dosing.
//
//   3. NO WEIGHT SCALING WHATSOEVER. V1 and V2 are fixed volumes in litres,
//      valid at the reference TBW of 78 kg. There is no /70 term to rescale
//      and none is implied. Using it at 50 kg or 120 kg is unscaled
//      extrapolation. For reference, the final model did find TBW^0.9 on
//      V1 -- but it explicitly did NOT find weight on V2 ("The BIC was not
//      improved by adding total body weight as a covariate for V2"), so
//      even the covariate model would not rescale the peripheral volume.
//      Observed weight range was 45-126 kg and no lean or adjusted body
//      weight was collected at all, so there is no obesity information
//      behind any of this.
//
//   4. NO RENAL COVARIATE, AND NO DIALYSIS INFORMATION AT ALL. Patients
//      requiring renal replacement therapy were excluded outright. Because
//      the basic model carries no renal term, a patient in renal failure is
//      assigned the typical CL of 4.0 L/h, which OVER-ESTIMATES clearance
//      and under-predicts accumulation and nephrotoxic exposure -- the
//      dangerous direction for an aminoglycoside. Conversely, augmented
//      renal clearance (the cohort reached 412 mL/min) is under-predicted.
//      Do not use this model for dialysis, and do not use it unsupervised
//      at the extremes of renal function.
//
//   5. SINGLE DOSE, 24-HOUR OBSERVATION WINDOW. Every patient received
//      exactly ONE infusion and was followed for 24 h with 5 samples. There
//      is no steady-state and no accumulation information anywhere in the
//      data. Simulating a repeated q24h regimen from this model is an
//      assumption of linearity and time-invariance, not an observation. The
//      terminal half-life is 7.28 h, so the 24 h window spans only about
//      3.3 terminal half-lives -- enough to see the beta phase but not to
//      exclude a slower one. The authors themselves flag a compartment the
//      design cannot see: total urine recovery of amikacin is incomplete at
//      24 h and aminoglycosides accumulate in the kidney. Multi-dose
//      accumulation is therefore likely to be UNDER-PREDICTED.
//
//   6. THE PERIPHERAL COMPARTMENT IS THE WEAKEST PART OF THE MODEL. V2
//      carries the largest random effect (omega 0.5) and rests on only
//      three late samples (8, 16, 24 h). Two covariates were tested on it
//      and both failed: creatinine clearance was dropped (p=0.4) and total
//      body weight did not improve BIC. The authors' own reading is that
//      "V2 might represent a weight-independent compartment in which
//      amikacin accumulates." Any conclusion that depends on the terminal
//      phase -- 24 h trough, accumulation on repeated dosing, tissue-side
//      exposure -- inherits that weakness.
//
//   7. OMEGA SCALE TRAP -- DO NOT RE-CONVERT. Table 2's omega is the SD on
//      the log scale and is squared above. The paper ALSO quotes omega*100
//      as a percentage in prose (basic model: "ranged from 30 % for Q to
//      60 % for CL", against printed 0.3 and 0.6; abstract: 31/22/27/47%
//      for the final model against printed 0.3/0.2/0.3/0.5). That is the
//      authors' labelling habit, not a CV in the
//      CV = sqrt(exp(omega^2)-1) sense -- checked arithmetically: under
//      that convention omega 0.47 would be quoted as 50%, not 47%, and
//      omega 0.31 as 32%, not 31%. Anyone who reads those percentages as
//      CV% and applies ln(1+CV^2) will get 0.30 instead of 0.36 for CL.
//
//   8. THE CORRELATION BLOCK IS IMPLEMENTED AS PRINTED, AND THAT IS A
//      DELIBERATE DEVIATION FROM THE HOUSE DIAGONAL $OMEGA. The basic
//      model's six correlations are all >= |0.2|, consistent with the
//      paper's own stated retention rule, and rho(CL,Q) = -0.4 and
//      rho(V1,Q) = 0.4 are large enough that discarding them would
//      misstate the simulated spread of peak and AUC. I verified positive
//      definiteness rather than assuming it: eigenvalues of the covariance
//      block are 0.0264 / 0.1086 / 0.2277 / 0.4272, all positive, so it
//      will Cholesky-decompose. CAVEATS: the correlations are printed to
//      ONE decimal place and the full variance-covariance matrix is
//      deferred to Online Resource 3, which is NOT in this 9-page PDF;
//      and several basic-model correlations have RSEs of 52-100%, so they
//      are directionally real but numerically soft. To fall back to
//      independent etas, replace the block with
//      `$OMEGA 0.36 0.09 0.09 0.25` -- nothing else needs to change.
//
//   9. THE COMBINED RESIDUAL ERROR FORM IS NOT SPECIFIED BY THE PAPER. It
//      defines only "a being the standard deviation of the additive
//      component and b the standard deviation of the proportional
//      component", which does not distinguish Monolix combined1
//      (SD = a + b*f) from combined2 (SD = sqrt(a^2 + (b*f)^2)). The
//      $TABLE idiom below implements combined2. With a = 0.2 mg/L against
//      concentrations of 40-60 mg/L the choice is immaterial except near
//      the 0.5 mg/L LOQ.
//
//  10. ASSAY AND CENSORING. Concentrations were fluorescence polarization
//      immunoassay (Innofluor), not LC-MS/MS; intra- and inter-assay CV
//      4.1% and 5.8%. Fourteen of 291 values (4.8%) were below the 0.5
//      mg/L LOQ and handled inside the SAEM likelihood -- good practice,
//      but it means the terminal-phase parameters rest partly on censored
//      data.
//
//  11. NO BONE, TISSUE, OR BIOFILM DATA OF ANY KIND. Plasma only. This
//      model says nothing about amikacin bone penetration; the library's
//      bone-penetration category for amikacin must be sourced elsewhere
//      and must NOT be inferred from anything here.
//
//  12. CLASS REDUNDANCY. Amikacin is an aminoglycoside, the same class as
//      the four gentamicin models already in the library. This file does
//      not touch the system's worst measured failure mode (unseen drug
//      class, R2 0.43); linezolid remains the higher-value gap.
//
//  13. TYPICAL-VALUE CURVE RUNS ABOVE THE PAPER'S OWN REPORTED MEDIAN C1h.
//      At the reference 78 kg and a nominal 20 mg/kg (1560 mg) this model
//      predicts C1h = 53 mg/L, against the reported observed median C1h of
//      45 mg/L (22-87). The most likely explanation is that actual
//      administered milligram doses were vial-rounded below 20 mg/kg x
//      78 kg (the dose range really was 11-28 mg/kg, and the paper states
//      it modelled "the total dose administered to each patient"), so the
//      cohort's median mg dose was below 1560 mg. Do not treat agreement
//      with the reported median C1h as a validation criterion for this
//      file.
//
// -------------------------------------------------------------------------
// PD TARGETS AND MIC SENSITIVITY (this paper, Results p.79 and Fig. 4 --
// recorded because it is in the article itself, not in the missing Online
// Resources, and because it bears directly on the biofilm MIC multiplier)
//   Targets used, both literature-referenced: C1h/MIC >= 10, and
//   AUC/MIC >= 90. EUCAST susceptibility breakpoint 8 mg/L for
//   Enterobacteriaceae and P. aeruginosa. Monte Carlo over doses
//   20/25/30/35/40 mg/kg and MICs 0.25-64 mg/L.
//     C1h/MIC >= 10:  20 mg/kg -> 100% at MIC <=2, 80% at MIC 4, 4% at
//                     MIC 8;  25 mg/kg -> 96% at MIC 4, 20% at MIC 8;
//                     40 mg/kg -> 80% at MIC 8.
//     AUC/MIC >= 90:  20 mg/kg -> 90% at MIC <=2, <20% at MIC 8;
//                     25 mg/kg -> 69% at MIC 4, 22% at MIC 8;
//                     40 mg/kg -> 90% at MIC 4, 52% at MIC 8.
//   NOTE FOR THE BIOFILM MULTIPLIER: a FOURFOLD rise in MIC (2 -> 8 mg/L)
//   collapses attainment from 100% to 4% at an unchanged 20 mg/kg dose.
//   The simulator's 8-256x biofilm multiplier is far larger than fourfold,
//   so under that multiplier no amikacin dose in this range attains target.
//   That is a real published data point against which the unsourced
//   multiplier's consequences can be sanity-checked -- it does not source
//   the multiplier, and this paper contains no biofilm data.
//
// SMOKE EXPECTATION
//   Amikacin 15 mg/kg q24h, 80 kg = 1200 mg, 30-min IV infusion, typical
//   patient (all ETA = 0):
//     end-of-infusion peak (t = 0.5 h)  61.9 mg/L   -- IN the 45-70 target
//     C1h (this paper's own PD index)   41.0 mg/L
//     24 h trough                        2.6 mg/L   -- ABOVE the usual
//                                                      aminoglycoside <1
//     AUC(0-24) 272 mg*h/L; AUCinf = dose/CL = 300 mg*h/L
//     Vss 37.4 L; t1/2-beta 7.28 h
//   Population spread with the omega block (2x10^5 draws): peak median
//   61.0, 5th-95th 39-94 mg/L; 24 h trough median 2.6, 5th-95th 0.2-11.5.
//   THE TROUGH DOES NOT MEET THE <1 mg/L AMINOGLYCOSIDE EXPECTATION, and
//   that is a genuine property of this ICU source population -- Vss 37.4 L
//   and a 7.3 h terminal half-life -- NOT a transcription error and NOT
//   something to tune away. Read it as this model's honest statement that
//   its patients did not clear a 24 h dosing interval. Parameters are as
//   printed in Table 2 and have not been adjusted.
// =========================================================================

[SET] end=72, delta=0.1


$PARAM @annotated
TVCL : 4.0  : Population clearance (L/h) -- true systemic CL, not CL/F
TVV1 : 15.3 : Population central volume (L) -- NOT per kg, see LIMITS 3
TVQ  : 12.2 : Population intercompartmental clearance (L/h)
TVV2 : 22.1 : Population peripheral volume (L) -- NOT per kg, see LIMITS 3

ETA1 : 0 : IIV on clearance
ETA2 : 0 : IIV on central volume
ETA3 : 0 : IIV on intercompartmental clearance
ETA4 : 0 : IIV on peripheral volume


$PARAM @annotated @covariates
// This model has NO covariate relations -- that is the point of the basic
// model (LIMITS 2). WT is declared only so the harness's WT column binds
// somewhere and so the reference weight is on the record. It does NOT
// enter any equation below. Dose in mg is supplied by the harness.
WT : 78 : reference total body weight (kg) -- DECLARED BUT UNUSED


// Log-scale variances on the diagonal, cov(i,j) = rho*omega_i*omega_j off
// it. Table 2 omegas are SDs on the log scale, so they are SQUARED, not
// passed through ln(1+CV^2) -- see PARAMETERS and LIMITS 7.
//   diag: 0.6^2=0.36  0.3^2=0.09  0.3^2=0.09  0.5^2=0.25
//   CL-V1:  0.3*0.6*0.3 =  0.054      CL-Q :  -0.4*0.6*0.3 = -0.072
//   V1-Q :  0.4*0.3*0.3 =  0.036      CL-V2:   0.3*0.6*0.5 =  0.090
//   V1-V2:  0.2*0.3*0.5 =  0.030      Q-V2 :   0.2*0.3*0.5 =  0.030
// Positive definite: eigenvalues 0.0264 0.1086 0.2277 0.4272 (LIMITS 8).
// Order is CL, V1, Q, V2 -> ETA(1..4).
$OMEGA @block
0.36
0.054  0.09
-0.072 0.036 0.09
0.090  0.030 0.030 0.25


// Published combined error, as variances: b = 0.1 -> 0.01 proportional,
// a = 0.2 mg/L -> 0.04 additive. These are the paper's values, not the
// house defaults. Combined form is combined2 -- see LIMITS 9.
$SIGMA
0.01
0.04


$CMT @annotated
CENT   : Central compartment (mg) [ADM, OBS]
PERIPH : Peripheral compartment (mg)


$MAIN
double CL = TVCL * exp(ETA1 + ETA(1));
double V1 = TVV1 * exp(ETA2 + ETA(2));
double Q  = TVQ  * exp(ETA3 + ETA(3));
double V2 = TVV2 * exp(ETA4 + ETA(4));


$ODE
dxdt_CENT   =  Q*(PERIPH/V2) - Q*(CENT/V1) - CL*(CENT/V1);
dxdt_PERIPH =  Q*(CENT/V1)   - Q*(PERIPH/V2);


$TABLE
double CP = CENT / V1;
double DV = CP * (1 + EPS(1)) + EPS(2);
int i = 0;
while(DV<0 && i <100) {
simeps();
DV = CP * (1 + EPS(1)) + EPS(2);
++i;
}


$CAPTURE DV CP CL V1 Q V2
