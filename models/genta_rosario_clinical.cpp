// ====================================================================
// genta_rosario_clinical.cpp -- Gentamicin, adult cancer patients,
//                    CLINICAL FACTORS parameterisation
//                    (Rosario 1998, Table 5 RIGHT column = model 19b)
//
// SOURCE
//   Rosario MC, Thomson AH, Jodrell DI, Sharp CA, Elliott HL.
//   "Population pharmacokinetics of gentamicin in patients with cancer."
//   Br J Clin Pharmacol 1998; 46: 229-236. DOI 10.1046/j.1365-2125.1998.00779.x
//   Clinical Pharmacokinetics and Biometrics Unit, University Department of
//   Medicine and Therapeutics, and Beatson Oncology Centre, West Glasgow
//   Hospitals University NHS Trust, Glasgow, UK.
//   PDF: Rosario - Population pharmacokinetics of
//        gentamicin in patients with cancer.pdf
//   Read here as pdftotext -layout (643 lines, all 8 pages) plus 300-400 dpi
//   renders of Table 1 (p.232) and Table 5 (p.234), whose columns the text
//   layer scrambles. The Read tool refuses this PDF with "pdftoppm is not
//   installed"; that failure is real, was worked around, and blocked nothing.
//
// POPULATION
//   Identical to genta_rosario_clcr.cpp -- SAME 210 adult human cancer
//   patients, same merged index (140) + test (70) data, same refit. HUMAN,
//   ADULT. Median age 50 y (15-81) index / 53 y (14-77) test; median weight
//   66 / 65 kg (38-117); median height 170 / 165 cm; median BSA 1.7 m2
//   (1.3-2.3); median albumin 34.0 g/L (17-53) against a 36-50 g/L reference
//   range; median serum creatinine 71 umol/L (26-258); Cockcroft-Gault CLcr
//   median 89.7 mL/min (22.9-319.7). Routine TDM, Glasgow, Jan 1993-Aug 1996,
//   IV gentamicin 40-300 mg titrated to a 1 h peak of 8-12 mg/L and trough
//   <2 mg/L. Patients with rapidly changing renal function excluded; care
//   setting never stated and dialysis never mentioned.
//
//   WHY THIS FILE EXISTS ALONGSIDE genta_rosario_clcr.cpp: the paper reports
//   TWO CO-EQUAL final models. "Comparison of the two best models (model 8a
//   and 19b - Table 4) by examination of the plots and by the error associated
//   with the parameter estimates allowed no firm conclusion to be drawn"; the
//   CLcr version was adopted "due to its simplicity". For an mrgsolve port
//   THIS one is the better target, because it uses raw creatinine, age and BSA
//   -- covariates the harness holds directly -- and therefore sidesteps the
//   unstated Cockcroft-Gault body-weight descriptor entirely. It needs no SEX
//   at all. See LIMITS 21 for what it costs.
//
// STRUCTURE
//   Two-compartment linear disposition, IV input. NONMEM version IV,
//   ADVAN3 / TRANS4, parameterised CL / V1 / Q / V2. CL, NOT CL/F -- IV route
//   throughout, F never estimated. Same structure as the CLcr model, different
//   covariate parameterisation of clearance. Random effects lognormal:
//   Pj = P * exp(eta_pj), eta ~ N(0, omega^2) (Methods, verbatim).
//
// PARAMETERS as printed in Table 5, right column (header reads
//   CL = th1*(CR70/71)^th5 * (age/46)^th6 * BSA,
//   V1 = th2*BSA*(alb/34)^th7):
//     th1  CL coefficient on BSA               2.46   L/h/m2
//     th5  creatinine exponent on CL          -0.78   dimensionless
//     th6  age exponent on CL                 -0.33   dimensionless
//     th2  V1 coefficient on BSA               8.52   L/m2
//     th7  albumin exponent on V1             -0.40   dimensionless
//     th3  Q                                   1.34   L/h
//     th4  V2                                  9.05   L
//   Derived, printed: CL 4.1 L/h, Vss 0.38 L/kg, t1/2,alpha 1.7 h (Discussion
//   gives the individual range for THIS 1.7 h figure as 0.8-3.8 h),
//   t1/2,z 7.5 h. Reproduced from the thetas: at BSA 1.7 m2, CR70 = 71 umol/L
//   and age 46 y, CL = 2.46 * 1 * 1 * 1.7 = 4.18 L/h vs printed 4.1.
//
//   CENTERING CONSTANTS -- three of them, and the first two look alike:
//     creatinine  FLOORED at 70 umol/L, then DIVIDED BY 71 umol/L
//                 (71 is the index-set median). Floor 70, divisor 71.
//                 The CLcr model of the same paper floors at 60 and does not
//                 centre at all. Mixing the two is a silent error.
//     age         divided by 46 years
//     albumin     divided by 34 g/L
//     BSA         MULTIPLICATIVE AND UNCENTRED in both CL and V1 -- it enters
//                 as a bare proportional factor with no exponent, so the
//                 thetas 2.46 and 8.52 are per m2, not typical values.
//
//   OMEGA DERIVATION -- spelled out because this is a failure mode the project
//   has already been bitten by. Table 5 labels the rows "CV omega_CL (%)" =
//   19.6 and "CV omega_Q (%)" = 30.1. The Table 5 FOOTNOTE glosses omega as
//   "variance of interpatient variability", contradicting the row heading.
//   The ABSTRACT resolves it for the companion model in the same units
//   ("Coefficient of variation was 18.5% on clearance and 28.2% on Q"), so
//   these are CV PERCENTAGES. With the lognormal form from Methods,
//   omega^2 = ln(1 + CV^2):
//       CL : ln(1 + 0.196^2) = ln(1.038416) = 0.037696   (omega = 0.1942)
//       Q  : ln(1 + 0.301^2) = ln(1.090601) = 0.086729   (omega = 0.2945)
//   Coding OMEGA = 19.6 from the footnote would inflate the SD ~5.3-fold.
//   NO IIV IS PRINTED FOR V1 OR V2 -- the "Interpatient variability" block of
//   Table 5 has exactly two rows. ETA3 on V1 carries OMEGA = 0 below; that is
//   the published estimate, not an omission. See LIMITS 6.
//
//   RESIDUAL ERROR as printed: combined, Y = C*exp(e1) + e2, variances
//   error1 = 0.025 and error2 = 0.102 (mg/L)^2. The pdftotext layer flattens
//   the superscript footnote marker into the digits and renders these as
//   "0.0252" and "0.1022"; the page render confirms three decimals plus a
//   marker, and the arithmetic corroborates it:
//       sqrt(1*0.025 + 0.102)  = 0.356 mg/L vs the paper's footnote 0.36
//       sqrt(64*0.025 + 0.102) = 1.305 mg/L vs the paper's footnote 1.30
//   (the four-decimal reading 0.0252/0.1022 gives 1.310 at 8 mg/L, a worse
//   match to the paper's own footnote -- so the three-decimal reading is
//   confirmed arithmetically as well as visually). The proportional term is
//   EXPONENTIAL in the paper, not (1+eps), and is coded that way in $TABLE.
//
// --------------------------------------------------------------------
// LIMITS -- every reason this model could mislead in an adult PJI cohort
//
//  1. THE COHORT IS PROFOUNDLY MYELOSUPPRESSED, AND THE PAPER DISPUTES THE
//     OBVIOUS CONCLUSION. Table 1 (400 dpi render): median white cell count
//     0.8 x10^9/L index and 1.1 test against a printed reference range of
//     4.1-11; median platelets 52 x10^9/L (reference 150-400); median
//     haemoglobin 10.4 g/dL (reference 13-18); 22.2% of index samples drawn at
//     >=38 C. The paper never calls its patients neutropenic and reports no
//     absolute neutrophil count, so "febrile-neutropenia haemato-oncology
//     cohort" is an INFERENCE from the WCC distribution. Vss is 0.38 L/kg.
//     BUT the AUTHORS EXPLICITLY REJECT the inflation reading: Discussion,
//     "In our study the value of 0.38 l kg-1, equivalent to 24.6 l, is
//     consistent with the values reported by both Bertino et al. and MacGowan
//     et al. in similar patients. However, the differences in volume of
//     distribution in the published literature may simply reflect different
//     sampling protocols and the application of different pharmacokinetic
//     models"; Conclusions, "Volume of distribution was similar to values
//     estimated by other authors". The 0.25-0.31 L/kg general-adult figure is
//     cited by the paper inside a paragraph that presents the question as
//     unresolved, immediately followed by Bertino's 235 cancer patients at
//     0.35 L/kg against 0.34 L/kg in controls. SAFEBONE POSITION: the worry
//     that an oncology TDM volume does not transport to an arthroplasty or PJI
//     patient is legitimate and belongs in the catalogue as OUR inference; if
//     V1 is inflated, simulated Cmax is depressed, which is the wrong
//     direction for an aminoglycoside judged on Cmax/MIC. Do not cite the
//     paper as the source of that claim.
//
//  2. THE HAEMATOLOGICAL COVARIATES WERE NEVER FITTED IN NONMEM. Only pyrexia,
//     obesity and gender were formally tested ("No significant improvement in
//     fit was obtained with the addition of pyrexia, obesity or gender to
//     either of the full models"). White cell count, haemoglobin and platelet
//     count were dropped at the graphical screening stage on scatter plots of
//     POSTHOC Bayesian estimates: "Haematological indices (haemoglobin, white
//     cell count, platelet count) showed no obvious trends". So the paper
//     provides no likelihood-ratio test of whether myelosuppression drives
//     volume, and there is no covariate here that can be switched off to
//     de-adjust the model back toward general-adult physiology.
//
//  3. RENAL RANGE AND SHAPE. Unlike the CLcr model, clearance here is a POWER
//     function of creatinine with exponent -0.78, so it does at least go to
//     zero as creatinine goes to infinity rather than resting on a 0.88 L/h
//     intercept. That does NOT make it safe at the low end: the authors' own
//     warning applies to the whole analysis -- "few patients had creatinine
//     clearance estimates <30 ml min-1. It is therefore important that these
//     results are not extrapolated to patients with very poor renal function
//     in whom gentamicin clearance might be overestimated." Serum creatinine
//     ranged 26-258 umol/L; nothing above 258 was observed, and an exponent of
//     -0.78 fitted over 70-258 umol/L extrapolated to 500 umol/L is
//     unsupported. Dialysis is nowhere in this data set -- use genta_franck
//     for that, not this.
//
//  4. NO COCKCROFT-GAULT, WHICH IS THIS MODEL'S MAIN ADVANTAGE. Because CL is
//     driven by raw creatinine, age and BSA, this file does NOT inherit the
//     unstated CG body-weight descriptor that afflicts the CLcr version (where
//     18% of concentrations came from obese patients, actual weight >120% of
//     ideal, median actual 66 kg vs median IBW 59.2 kg). It also needs no SEX
//     at all: gender was significant as a single covariate but NOT in the
//     final clinical-factors model, and the authors attribute that to BSA
//     already absorbing it. If gender matters in the PJI cohort, this model
//     cannot express it.
//
//  5. THE CREATININE FLOOR IS LOAD-BEARING AND ITS CONSTANTS ARE CONFUSABLE.
//     Creatinine <=70 umol/L is raised to 70, then divided by 71. The 70 floor
//     gave the optimal fit for THIS model (objective function improvement 35).
//     The companion CLcr model floors at 60 and centres at nothing. Floor 70,
//     divisor 71, and they are different numbers doing different jobs.
//     Omitting the floor is a misimplementation, not a simplification: without
//     it a PJI patient with a low creatinine from low muscle mass -- common
//     after a long orthopaedic illness -- gets a clearance the model was never
//     fitted to support.
//
//  6. NO BETWEEN-SUBJECT VARIABILITY ON V1 OR V2. This is the largest
//     practical obstacle to Monte Carlo target attainment from this file,
//     because V1 sets the peak and the peak is the aminoglycoside index. Run
//     AS PUBLISHED, a virtual population from this file UNDERSTATES the spread
//     of Cmax. The only donor anywhere in the paper is the covariate-free base
//     model's 15.5% CV on V1 -- ln(1+0.155^2) = 0.023741 -- which belongs to a
//     different model fitted to a different 140-patient subset with no
//     covariates. It is deliberately NOT applied here; the ETA3 slot exists
//     with OMEGA 0 so that an analyst who needs it substitutes it knowingly
//     and records it as a kludge.
//
//  7. NO STANDARD ERRORS, RSEs OR CONFIDENCE INTERVALS ANYWHERE IN THE PAPER.
//     Table 5 has exactly two columns, "Parameter" and "Estimate" (confirmed
//     on the render). Parameter uncertainty cannot be propagated. No
//     off-diagonal omega elements are printed and the paper never states the
//     omega matrix was diagonal; diagonal is ASSUMED here.
//
//  8. THE TERMINAL PHASE IS UNDER-IDENTIFIED BY DESIGN, AND THE AUTHORS SAY
//     SO. Sampling ran 1-26 h post-dose; 71% of samples were within 12 h and
//     only two patients contributed a sample beyond 25 h. t1/2,z 7.5 h is far
//     shorter than the deep tissue-release phase (order of 100 h) reported
//     elsewhere for gentamicin. Over the multi-week course a PJI simulation
//     needs, this model will UNDER-PREDICT accumulation, therefore under-
//     predict trough rise and nephrotoxicity exposure.
//
//  9. THE INFUSION DURATION CODED IN THE NONMEM DATA SET IS NOT STATED. The
//     Methods say only "either by a short infusion over 10-30 min or as a slow
//     bolus over 2-3 min", mixed in the same data set. With an alpha half-life
//     of 1.7 h the assumed duration moves Cmax appreciably. The harness
//     supplies it on the dose record; no duration in this file can be said to
//     match the fit.
//
// 10. FITTED ALMOST ENTIRELY TO LOW CONCENTRATIONS UNDER CONVENTIONAL DOSING.
//     57% of the measured concentrations were TROUGHS; the whole range was
//     0.1-13.5 mg/L with a MEDIAN OF 1.6 mg/L; doses were 40-300 mg titrated
//     to a 1 h peak of 8-12 mg/L. A modern 5-7 mg/kg once-daily simulation
//     with a 20+ mg/L peak lies ABOVE the entire observed concentration range,
//     where linearity is assumed rather than demonstrated.
//
// 11. AGE IS A COVARIATE HERE, AND ITS RANGE MATTERS. (age/46)^-0.33 was
//     fitted over 14-81 y. The minimum age of 14-15 y means a few adolescents
//     are pooled in -- this is NOT a paediatric model (reference weight 66 kg,
//     not 6.8 kg) but it is not purely adult either. At the elderly end the
//     PJI population is older than this cohort: an 80 y patient gets
//     (80/46)^-0.33 = 0.83, an 18% clearance reduction driven by an exponent
//     estimated where the data thin out.
//
// 12. THE EXTERNAL EVALUATION PASSED, BUT WITH A NON-STANDARD DENOMINATOR AND
//     IN A DIRECTION WORTH KNOWING. Mean %PE was -6.6% for this model on 70
//     held-out patients / 132 concentrations, no significant bias, and 90% of
//     measured concentrations fell inside the 95% interval of 1000 simulations
//     against a nominal 95%. But the paper defines %PE = (C_pred - Y_meas)*100
//     / C_PRED -- normalised to the PREDICTED concentration, not the measured
//     one -- so it is not comparable to the usual metric, and the negative
//     sign means the model UNDER-PREDICTS observed concentrations by ~7%.
//
// 13. NO BONE, TISSUE, MICRODIALYSIS, HOMOGENATE, MIC OR BIOFILM DATA OF ANY
//     KIND. Plasma-only TDM paper. It contributes nothing to the bone-
//     penetration-ratio problem and nothing to the unsourced 8-256x biofilm
//     MIC multiplier. Any bone concentration from this file is a plasma
//     concentration times a number that came from somewhere else.
//
// 14. INDICATION MISMATCH. Gentamicin was given for prophylaxis and treatment
//     of gram-negative infection in cancer patients. No prosthesis, no bone
//     infection, no surgical physiology, no long course in the data.
//
// 15. 1998 VINTAGE. Fluorescence polarization immunoassay (TDx, LOQ 0.1 mg/L,
//     interassay CV 6.3% at 1 mg/L, 3.7% at 4, 4.3% at 8). NONMEM "version IV"
//     and nothing more -- the ESTIMATION METHOD (FO vs FOCE vs FOCE-I) IS
//     NEVER STATED. FO-era variance components from sparse TDM data are
//     routinely biased; a modern reanalysis would likely give different
//     omegas.
//
// 16. NO INTER-OCCASION VARIABILITY, although three patients contributed
//     several separate courses (14, 17 and 22 measurements) pooled without an
//     occasion random effect. Part of what the residual sigma below absorbs is
//     really between-occasion variability, inflating sigma and deflating the
//     CL omega.
//
// 17. "NON-ICU, NON-DIALYSIS" IS AN ARGUMENT FROM SILENCE. The only stated
//     exclusions are missing clinical/dosage/sampling data and rapidly
//     changing renal function. The care setting is never characterised and
//     dialysis is never mentioned in eight pages. CLcr does reach 15.8-22.9
//     mL/min at the low end. Catalogue it as "setting not stated".
//
// 18. ALBUMIN CENTERED BELOW NORMAL. The centering value 34 g/L is the study
//     median and is BELOW the paper's own 36-50 g/L reference range. A
//     normo-albuminaemic PJI patient at 40 g/L therefore gets V1 scaled DOWN
//     by (40/34)^-0.40 = 0.937, about 6%. The paper quantifies the covariate's
//     reach: Vss ran 0.55 down to 0.35 L/kg across the observed albumin range
//     of 17-53 g/L. Outside 17-53 g/L this power function is extrapolation.
//     The harness ALB unit is not declared; this model needs g/L (see the
//     guard in $MAIN).
//
// 19. THIS MODEL RECEIVED LESS DOWNSTREAM SCRUTINY WITHIN THE PAPER THAN ITS
//     TWIN. Only the CLcr model (8a) was carried into the paper's nomogram
//     simulations and Figure 2. Nothing in the paper plots or stress-tests
//     THIS parameterisation beyond the Table 5 estimates and the -6.6% MPE, so
//     there is no published concentration-time figure to regress this file
//     against. Use the Figure 2 cross-check on genta_rosario_clcr.cpp for
//     that, and treat agreement between the two files as a self-consistency
//     check only.
//
// 20. STRATEGIC. This is a gentamicin model, and the library already codes the
//     aminoglycoside class (gomes, evers_standard, evers_icu, franck). Per the
//     project's own measurement that the worst failure mode is an unseen drug
//     CLASS (R2 0.43 +/- 0.13), this is worth less than the still-missing
//     linezolid oxazolidinone model. Its genuine contribution is structural:
//     it is a TWO-compartment model with albumin on central volume, so it is
//     not another re-estimation of the one-compartment Gomes/Evers structure.
//
// 21. THIS MODEL AND genta_rosario_clcr.cpp ARE THE SAME DATA, NOT INDEPENDENT
//     EVIDENCE. Two co-equal covariate parameterisations refitted on the same
//     210 patients. "Comparison of the two best models (model 8a and 19b) ...
//     allowed no firm conclusion to be drawn." Never treat their agreement as
//     corroboration and never hold one out to validate the other. Wire ONE of
//     them into any given analysis.
//
// SMOKE EXPECTATION
//   Standard adult regimen: gentamicin 5 mg/kg q24h = 400 mg over 30 min,
//   80 kg, 175 cm, 65 y, CREAT 80 umol/L, ALB 40 g/L (SEX is not used).
//     -> BSA (Du Bois) 1.96 m2; creatinine 80 is above the 70 umol/L floor so
//        no truncation applies; (80/71)^-0.78 = 0.911; (65/46)^-0.33 = 0.892
//     -> CL 3.91 L/h, V1 15.62 L, Q 1.34 L/h, V2 9.05 L, Vss 24.7 L
//     -> t1/2,alpha 1.78 h, t1/2,z 7.27 h
//     -> FIRST-DOSE PEAK at end of infusion (t = 0.5 h) 23.6 mg/L
//     -> 24 h TROUGH 0.48 mg/L (steady state ~0.55 mg/L; the q24h accumulation
//        factor on the beta phase is only 1.12)
//   IN RANGE: the target band is peak 15-25 mg/L and 24 h trough <1 mg/L, and
//   this model gives 23.6 mg/L and 0.48 mg/L. As with its twin it sits in the
//   UPPER half of the peak band, because V1 scales on BSA rather than weight,
//   so V1/kg falls as the patient gets heavier than the 66 kg cohort median.
//   Expect systematically HIGHER peaks than the one-compartment Gomes/Evers
//   gentamicin files at the same dose (genta_gomes gives 17.3 mg/L for a
//   comparable patient) -- a real between-model disagreement of the size the
//   project has already measured (1.4x to 4.2x on AUC24), not a defect.
//   Cross-check that a correct build must reproduce: at BSA 1.7 m2, CREAT
//   71 umol/L and age 46 y, CL = 2.46 * 1.7 = 4.18 L/h against the Table 5
//   derived value of 4.1 L/h; at BSA 1.7 and ALB 34 g/L, V1 = 14.48 L and
//   Vss = 14.48 + 9.05 = 23.5 L. NOTE that 23.5 L over the 66 kg cohort median
//   is 0.356 L/kg, whereas Table 5 prints Vss = 0.38 L/kg for BOTH models. The
//   CLcr model reconstructs to 24.4 L = 0.37 L/kg, so the 0.38 figure is
//   evidently a rounded value carried across both columns rather than a
//   per-model recomputation. Do not treat a 0.36 L/kg build as a bug; treat a
//   build that does not give V1 = 14.48 L at BSA 1.7 / ALB 34 as a bug.
// ====================================================================

[SET] end=72, delta=0.1


$PARAM @annotated
TVCLBSA : 2.46  : CL coefficient on body surface area, theta1 (L/h per m2)
CREXP   : -0.78 : exponent of serum creatinine on CL, theta5 (dimensionless)
AGEEXP  : -0.33 : exponent of age on CL, theta6 (dimensionless)
TVV1BSA : 8.52  : V1 coefficient on body surface area, theta2 (L/m2)
ALBEXP  : -0.40 : exponent of albumin on V1, theta7 (dimensionless)
TVQ     : 1.34  : intercompartmental clearance, theta3 (L/h)
TVV2    : 9.05  : peripheral volume, theta4 (L)

CRFLOOR : 70.0  : creatinine floor (umol/L) -- values <=70 are set to 70
CRREF   : 71.0  : creatinine centering value (umol/L) -- index-set median
AGEREF  : 46.0  : age centering value (years)
ALBREF  : 34.0  : albumin centering value (g/L) -- study median, below normal

ETA1 : 0 : IIV on clearance
ETA2 : 0 : IIV on intercompartmental clearance
ETA3 : 0 : IIV on central volume -- NOT PUBLISHED, omega fixed to 0, LIMITS 6

// FLOOR 70, DIVISOR 71. Two similar numbers doing two different jobs. The
// companion CLcr model floors at 60 and centres at nothing. See LIMITS 5.


$PARAM @annotated @covariates
WT    : 70.0  : total body weight (kg) -- paper's own reference patient
HT    : 170.0 : height (cm)
AGE   : 50.0  : age (years)
CREAT : 71.0  : serum creatinine (umol/L) -- cohort median
ALB   : 40.0  : serum albumin (g/L) -- reference patient value
SEX   : 0     : 0 = MALE, 1 = FEMALE -- NOT USED by this model

// Defaults above are the paper's OWN simulated reference patient (Dosage
// guidelines / Figure 2: 50 y, 70 kg, 170 cm, albumin 40 g/L), not the cohort
// medians. SEX is declared only so the harness can pass its standard covariate
// set unchanged: gender was significant as a single covariate but NOT in the
// final clinical factors model, and this model contains no gender term at all.


// OMEGA. Table 5 rows are "CV omega_CL (%)" = 19.6 and "CV omega_Q (%)" =
// 30.1. These are COEFFICIENTS OF VARIATION, not variances: the row headings
// say CV(%), the Abstract confirms the convention for the companion model
// ("Coefficient of variation was 18.5% on clearance and 28.2% on Q"), and only
// the Table 5 footnote dissents by calling omega a variance. Lognormal random
// effects (Methods: Pj = P*Exp(eta_pj)), so omega^2 = ln(1 + CV^2):
//   ETA1 CL : ln(1 + 0.196^2) = ln(1.038416) = 0.037696
//   ETA2 Q  : ln(1 + 0.301^2) = ln(1.090601) = 0.086729
//   ETA3 V1 : NOT REPORTED IN THE FINAL MODEL -> 0. The paper's covariate-free
//             base model offers 15.5% CV, i.e. ln(1+0.155^2) = 0.023741, but
//             that is a different model on a different subset. Substitute it
//             only deliberately and record it as a kludge. See LIMITS 6.
// Omega assumed DIAGONAL; the paper prints no off-diagonals and never says.
$OMEGA 0.037696 0.086729 0.000000


// Combined residual error as published: Y = C*exp(e1) + e2, variances
// error1 = 0.025 (dimensionless, exponential) and error2 = 0.102 (mg/L)^2.
// Three-decimal printed values; the trailing digit in the PDF text layer is a
// superscript footnote marker, not a fourth decimal.
$SIGMA 0.025 0.102


$CMT @annotated
CENT   : Central compartment (mg) [ADM, OBS]
PERIPH : Peripheral compartment (mg)


$MAIN
// --- body size -------------------------------------------------------
// Du Bois body surface area (the paper's reference 16 is Du Bois & Du Bois,
// "Clinical calorimetry"), WT in kg and HT in cm. BSA drives BOTH CL and V1
// here, multiplicatively and uncentred.
double BSA = 0.007184 * pow(WT, 0.425) * pow(HT, 0.725);
if(BSA < 0.5) BSA = 0.5;      // guard only; observed range was 1.3-2.3 m2
if(BSA > 3.5) BSA = 3.5;

// --- albumin ---------------------------------------------------------
// The harness carries ALB without a declared unit; this model needs g/L
// (centred on 34). A value below 10 can only be g/dL, so it is converted here
// explicitly rather than silently mis-scaled -- 4.0 g/dL taken as 4 g/L would
// inflate V1 by (34/4)^0.40 = 2.4-fold. The value used is captured as ALBGL.
double ALBGL = (ALB < 10.0) ? ALB * 10.0 : ALB;
if(ALBGL < 10.0) ALBGL = 10.0;   // observed albumin range was 17-53 g/L
if(ALBGL > 60.0) ALBGL = 60.0;

// --- renal function --------------------------------------------------
// "All creatinine concentrations less than or equal to 70 umol/L were set to
// 70 umol/L", then centred on 71 umol/L. CREAT is umol/L in the harness and in
// the paper (the PDF text layer renders micro as "m"; it is NOT mmol/L).
// No Cockcroft-Gault and no SEX are involved -- that is this model's point.
double CREAT70 = (CREAT < CRFLOOR) ? CRFLOOR : CREAT;

// --- parameters ------------------------------------------------------
double CL = TVCLBSA * pow(CREAT70 / CRREF, CREXP)
                    * pow(AGE / AGEREF, AGEEXP)
                    * BSA * exp(ETA1 + ETA(1));
double Q  = TVQ * exp(ETA2 + ETA(2));
double V1 = TVV1BSA * BSA * pow(ALBGL / ALBREF, ALBEXP) * exp(ETA3 + ETA(3));
double V2 = TVV2;


$ODE
dxdt_CENT   = - (CL / V1) * CENT - (Q / V1) * CENT + (Q / V2) * PERIPH;
dxdt_PERIPH =   (Q / V1) * CENT - (Q / V2) * PERIPH;


$TABLE
// Combined error exactly as published: Y = C * exp(e1) + e2.
// Exponential proportional term is the paper's form, not the house
// (1+EPS(1)) default.
double CP = CENT / V1;
double DV = CP * exp(EPS(1)) + EPS(2);
int i = 0;
while(DV<0 && i <100) {
simeps();
DV = CP * exp(EPS(1)) + EPS(2);
++i;
}


$CAPTURE DV CP CL V1 V2 Q BSA CREAT70 ALBGL
