// ====================================================================
// genta_rosario_clcr.cpp -- Gentamicin, adult cancer patients
//                           (Rosario 1998, Table 5 LEFT column = model 8a)
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
//   page renders of Table 1 (p.232) and Table 5 (p.234), because the text
//   layer scrambles the column alignment of both tables. The Read tool itself
//   refuses this PDF with "pdftoppm is not installed"; that failure is real
//   and was worked around, not used as an excuse. Nothing was unreadable.
//
// POPULATION
//   HUMAN. ADULT. 210 patients with cancer (solid and haematological; the
//   split is not reported) treated at West Glasgow Hospitals, Jan 1993 -
//   Aug 1996, from ROUTINE THERAPEUTIC DRUG MONITORING. 140 patients / 378
//   concentrations in the index set, 70 patients / 132 concentrations in the
//   held-out evaluation set; the Table 5 parameters below come from a refit
//   on the merged 210 (the concentration count for that fit is not stated).
//     Age      median 50 y (15-81) index, 53 y (14-77) test. Model B of the
//              same paper centres age on 46 y. NOT paediatric; a handful of
//              14-15 y adolescents are pooled in (see LIMITS 12).
//     Weight   median 66 kg (38-117) index, 65 kg (38-95) test.
//     Height   median 170 cm (146-188) / 165 cm (147-188).
//     BSA      median 1.7 m2 (1.3-2.3). Median IBW 59.2 kg (38.0-82.2).
//     Albumin  median 34.0 g/L (17-53), i.e. the typical patient is MILDLY
//              HYPOALBUMINAEMIC against the paper's own 36-50 g/L range, and
//              the model is centred there.
//     Renal    Cockcroft-Gault CLcr median 89.7 mL/min (22.9-319.7); with the
//              60 umol/L creatinine floor applied, 86.1 mL/min (22.9-183.7).
//              Serum creatinine median 71 umol/L (26-258), reference 60-110.
//              Patients with RAPIDLY CHANGING renal function were EXCLUDED.
//     Setting  NEVER STATED. The paper does not say ward, oncology unit or
//              ICU, and never mentions dialysis or renal replacement anywhere
//              in eight pages. Do not record this as a "non-ICU, non-dialysis"
//              model; record it as "setting not stated, unstable renal
//              function excluded, dialysis never mentioned". See LIMITS 18.
//     Dosing   IV, individual doses 40-300 mg, initial dose by a local renal
//              nomogram, titrated to a 1 h peak of 8-12 mg/L and trough
//              <2 mg/L -- i.e. CONVENTIONAL multiple-daily dosing targets,
//              not the modern 5-7 mg/kg once-daily peak. See LIMITS 11.
//
// STRUCTURE
//   Two-compartment linear disposition, IV input, first-order elimination.
//   NONMEM version IV, ADVAN3 / TRANS4, parameterised CL / V1 / Q / V2.
//   CL, NOT CL/F -- the route is intravenous throughout and F was never
//   estimated. Two compartments beat one on a log-likelihood difference of
//   86.3 and removed a time trend in the weighted residuals.
//   Random effects lognormal: Pj = P * exp(eta_pj), eta ~ N(0, omega^2)
//   (Methods, verbatim).
//
// PARAMETERS as printed in Table 5, left column (header reads
//   CL = th1*(1 + th5*CLcr60), V1 = th2*BSA*(alb/34)^th6):
//     th1  CL intercept                        0.88   L/h
//     th5  slope of CL on CLcr60               0.043  per (mL/min)
//     th2  V1 coefficient on BSA               8.59   L/m2
//     th6  albumin exponent on V1             -0.39   dimensionless
//     th3  Q                                   1.30   L/h
//     th4  V2                                  9.79   L (abstract rounds 9.80)
//   Derived, printed: CL 4.2 L/h, Vss 0.38 L/kg (24.6 L), t1/2,alpha 1.8 h,
//   t1/2,z 8.0 h (Discussion gives the individual terminal range 3.7-17.8 h;
//   the initial-phase range 0.8-3.8 h quoted in the Discussion belongs to the
//   1.7 h figure of the CLINICAL FACTORS model, not to this 1.8 h column).
//   Reproduced here from the thetas: at the truncated median CLcr 86.1 mL/min,
//   CL = 0.88*(1+0.043*86.1) = 4.14 L/h vs printed 4.2; at BSA 1.7 and
//   ALB 34, V1 = 14.6 L, Vss = 14.6+9.79 = 24.4 L vs printed 24.6 L.
//
//   OMEGA DERIVATION -- this is the trap this project has already been bitten
//   by, so it is spelled out. Table 5 labels the rows "CV omega_CL (%)" = 18.5
//   and "CV omega_Q (%)" = 28.2. The Table 5 FOOTNOTE glosses the symbol omega
//   as "variance of interpatient variability", which contradicts the row
//   heading. The ABSTRACT settles it, verbatim: "Coefficient of variation was
//   18.5% on clearance and 28.2% on Q." So 18.5 and 28.2 are CV PERCENTAGES.
//   With the lognormal form stated in Methods, omega^2 = ln(1 + CV^2):
//       CL : ln(1 + 0.185^2) = ln(1.034225) = 0.033652   (omega = 0.1835)
//       Q  : ln(1 + 0.282^2) = ln(1.079524) = 0.076520   (omega = 0.2766)
//   An implementer who followed the footnote and coded OMEGA = 18.5 would
//   inflate the SD about 5.4-fold. (The naive CV^2 reading, 0.0342 and 0.0795,
//   is what a careless conversion gives; the difference from ln(1+CV^2) is
//   immaterial here but the correct form is used.)
//   NO IIV IS PRINTED FOR V1 OR V2. The "Interpatient variability" block of
//   Table 5 has exactly two rows. ETA3 on V1 therefore carries OMEGA = 0
//   below -- that is the PAPER'S estimate as published, not an omission here.
//   See LIMITS 6 for the only available donor.
//
//   RESIDUAL ERROR as printed: combined, Y = C*exp(e1) + e2, with variances
//   error1 = 0.025 (exponential/proportional) and error2 = 0.107 (additive,
//   (mg/L)^2). The pdftotext layer flattens the SUPERSCRIPT FOOTNOTE MARKERS
//   into the numbers and renders these as "0.0251" and "0.1071"; the page
//   render confirms three decimals plus a marker. Arithmetic corroborates the
//   three-decimal reading against the paper's own footnote:
//       sqrt(1*0.025 + 0.107)  = 0.363 mg/L vs footnote 0.36 at 1 mg/L
//       sqrt(64*0.025 + 0.107) = 1.307 mg/L vs footnote 1.32 at 8 mg/L
//   The proportional term is EXPONENTIAL in the paper, not (1+eps), and is
//   implemented that way in $TABLE.
//
// --------------------------------------------------------------------
// LIMITS -- every reason this model could mislead in an adult PJI cohort
//
//  1. THE COHORT IS PROFOUNDLY MYELOSUPPRESSED, AND THE PAPER DOES NOT AGREE
//     WITH THE OBVIOUS CONCLUSION. Table 1 (read from a 400 dpi render):
//     median white cell count 0.8 x10^9/L index and 1.1 test, against a
//     printed reference range of 4.1-11; median platelets 52 x10^9/L
//     (reference 150-400); median haemoglobin 10.4 g/dL (reference 13-18);
//     22.2% of index-set samples drawn at temperature >=38 C. The paper never
//     uses the word neutropenic about its own patients and never reports an
//     absolute neutrophil count, so "febrile neutropenia cohort" is an
//     INFERENCE from the WCC distribution, not a stated classification.
//     Vss is 0.38 L/kg. HOWEVER -- and this must not be laundered -- the
//     AUTHORS EXPLICITLY REJECT the reading that this is inflated: Discussion,
//     "In our study the value of 0.38 l kg-1, equivalent to 24.6 l, is
//     consistent with the values reported by both Bertino et al. and MacGowan
//     et al. in similar patients. However, the differences in volume of
//     distribution in the published literature may simply reflect different
//     sampling protocols and the application of different pharmacokinetic
//     models"; Conclusions, "Volume of distribution was similar to values
//     estimated by other authors"; Introduction, "the reason for the apparent
//     increase in volume is not known". The 0.25-0.31 L/kg general-adult range
//     is cited by the paper inside a paragraph presenting the question as
//     UNRESOLVED, immediately followed by Bertino's 235 cancer patients at
//     0.35 L/kg vs 0.34 L/kg in controls -- i.e. no increase in cancer.
//     PRACTICAL POSITION FOR SAFEBONE: an oncology TDM cohort's central volume
//     may still not transport to an elective-arthroplasty or PJI patient, and
//     if it is inflated the bias depresses simulated Cmax -- the wrong
//     direction for an aminoglycoside judged on Cmax/MIC. Carry that as OUR
//     concern, labelled as ours. Do not cite the paper for it.
//
//  2. THE HAEMATOLOGICAL COVARIATES WERE NEVER FITTED. Only pyrexia, obesity
//     and gender were entered into NONMEM ("No significant improvement in fit
//     was obtained with the addition of pyrexia, obesity or gender to either
//     of the full models"). White cell count, haemoglobin and platelet count
//     were dropped at the GRAPHICAL screening stage -- "Haematological indices
//     (haemoglobin, white cell count, platelet count) showed no obvious trends
//     nor did the presence or absence of pyrexia" -- which refers to scatter
//     plots of POSTHOC Bayesian estimates, not to a likelihood-ratio test.
//     So the paper offers NO formal test of whether myelosuppression drives
//     volume, and there is no covariate in this model that could be switched
//     off to de-adjust it back toward a general-adult Vd. Absence of evidence,
//     from an eyeballed plot, is weak evidence of absence.
//
//  3. RENAL RANGE AND AN UNSAFE EXTRAPOLATION THE AUTHORS THEMSELVES FLAG.
//     CL = 0.88*(1 + 0.043*CLcr) is ADDITIVE-LINEAR with a NON-ZERO INTERCEPT
//     of 0.88 L/h at zero renal function. The paper: "few patients had
//     creatinine clearance estimates <30 ml min-1. It is therefore important
//     that these results are not extrapolated to patients with very poor renal
//     function in whom gentamicin clearance might be overestimated." A PJI
//     patient with CLcr below ~30 mL/min simulated with this model gets too
//     much clearance, therefore too little accumulation and too low a
//     predicted nephrotoxicity risk over a multi-week course. That is the
//     unsafe direction.
//
//  4. THE COCKCROFT-GAULT VARIANT IS NOT SPECIFIED. The paper never states
//     whether CLcr was computed on total, ideal or adjusted body weight. This
//     implementation uses TOTAL body weight with the arbitrary 0.85 female
//     factor (the paper itself notes the female factor "was chosen
//     arbitrarily"). 18% of concentrations came from OBESE patients (actual
//     weight >120% of ideal; median actual 66 kg against median IBW 59.2 kg),
//     so the choice moves CL materially -- in an obese PJI patient TBW-based
//     CG can overstate CLcr by tens of percent. The companion file
//     genta_rosario_clinical.cpp does NOT have this problem and is the better
//     choice where it applies.
//
//  5. UNIT TRAP ON CLcr. The slope 0.043 is per mL/min and CLcr enters RAW and
//     UNCENTRED. A harness that supplied renal function in L/h would collapse
//     CL to about 0.9 L/h regardless of renal function. CREAT is umol/L
//     throughout the paper (the text layer renders every micro sign as "m", so
//     the extracted text falsely reads "mmol l-1"; the render, the 60-110
//     reference range and the paper's own "88.4 umol l-1 (1 mg dl-1)"
//     equivalence all confirm micromol). Note that the Table 1 UREA row is
//     genuinely mmol/L -- applying the micro-correction globally to the text
//     layer would corrupt urea 1000-fold. Urea is not in this model.
//
//  6. THE CREATININE FLOOR IS LOAD-BEARING, AND THE TWO ROSARIO MODELS USE
//     DIFFERENT ONES. Here: creatinine <=60 umol/L is raised to 60 BEFORE the
//     Cockcroft-Gault calculation, and CLcr is then NOT centred. Truncating
//     improved the objective function by 42 (floor 60) and 45 (floor 70)
//     against only 12.9 for the commonly recommended 88.4 umol/L floor.
//     Implementing this model without the floor is a misimplementation. The
//     clinical-factors model floors at 70 AND divides by 71 -- two similar
//     numbers that are easy to conflate silently. Note also that the 70 floor
//     scored marginally BETTER here (45 vs 42); the authors adopted 60 anyway
//     and 60 is what Table 5 documents, so 60 is what is coded.
//
//  7. NO BETWEEN-SUBJECT VARIABILITY ON V1 OR V2. This is the single largest
//     obstacle to using this file for Monte Carlo target attainment, because
//     V1 is what sets the peak and the peak is what an aminoglycoside is
//     judged on. A virtual population run from this file AS PUBLISHED will
//     understate the spread of Cmax. The only donor in the paper is the
//     covariate-free base model's 15.5% CV on V1, which belongs to a DIFFERENT
//     model fitted to a DIFFERENT 140-patient subset without covariates;
//     ln(1 + 0.155^2) = 0.023741 would be the substitution. It is a documented
//     kludge, not a published parameter, and it is deliberately NOT applied
//     here -- the ETA3 slot exists with OMEGA 0 so an analyst can make that
//     substitution consciously and record it.
//
//  8. NO STANDARD ERRORS, RSEs OR CONFIDENCE INTERVALS ANYWHERE. Table 5 has
//     exactly two columns, "Parameter" and "Estimate" (confirmed on the page
//     render). Parameter uncertainty cannot be propagated from this source.
//     There is also no off-diagonal omega and no statement that the omega
//     matrix was diagonal; diagonal is ASSUMED here.
//
//  9. THE TERMINAL PHASE IS UNDER-IDENTIFIED BY DESIGN, AND THE AUTHORS SAY
//     SO. Sampling times ran 1-26 h post-dose; 71% of samples were within 12 h
//     and only two patients contributed a sample beyond 25 h. The estimated
//     t1/2,z of 8.0 h is far shorter than the deep tissue-release phase (of
//     the order of 100 h) reported elsewhere for gentamicin. For a PJI use
//     case involving a multi-day or multi-week course, this model will
//     UNDER-PREDICT accumulation and therefore under-predict trough rise and
//     nephrotoxicity exposure. It also means the model carries no information
//     about the deep compartment that any tissue-distribution argument would
//     need.
//
// 10. THE INFUSION DURATION CODED IN THE NONMEM DATA SET IS NOT STATED. The
//     Methods say only that gentamicin was given "either by a short infusion
//     over 10-30 min or as a slow bolus over 2-3 min", mixed within the same
//     data set. With an alpha half-life of 1.8 h the assumed infusion duration
//     moves simulated Cmax appreciably. The harness must supply the duration
//     on the dose record; nothing in this file fixes it, and no duration in
//     this file can be said to match the fit.
//
// 11. THE MODEL WAS FITTED ALMOST ENTIRELY TO LOW CONCENTRATIONS UNDER A
//     CONVENTIONAL DOSING TARGET. 57% of the measured concentrations were
//     TROUGHS; the whole concentration range was 0.1-13.5 mg/L with a MEDIAN
//     OF 1.6 mg/L; doses were 40-300 mg titrated to a 1 h peak of 8-12 mg/L.
//     Simulating a modern 5-7 mg/kg once-daily regimen with a 20+ mg/L peak
//     puts the simulation ABOVE the entire observed concentration range. The
//     linear disposition assumption is untested there.
//
// 12. AGE RANGE INCLUDES ADOLESCENTS. Minimum age 14-15 y across the two sets.
//     The reference patient is a 46-50 y adult and this is emphatically NOT a
//     paediatric model (no 6.8 kg reference weight problem here -- reference
//     weight is 66 kg), but the pooled data are not purely adult.
//
// 13. THE EXTERNAL EVALUATION PASSED, BUT IN A DIRECTION WORTH KNOWING, AND
//     WITH A NON-STANDARD DENOMINATOR. Mean %PE was -7.2% on 70 held-out
//     patients / 132 concentrations, with no significant bias, and 90% of
//     measured concentrations fell inside the 95% interval of 1000 simulations
//     against a nominal 95%. But the paper defines %PE = (C_pred - Y_meas)*100
//     / C_PRED -- normalised to the PREDICTED concentration, not the measured
//     one, so it is not comparable to the conventional metric. The negative
//     sign means the model UNDER-PREDICTS observed concentrations by about 7%.
//     For a Cmax/MIC argument, an under-prediction is conservative on efficacy
//     and anti-conservative on toxicity.
//
// 14. NO BONE, TISSUE, MICRODIALYSIS, HOMOGENATE, MIC OR BIOFILM DATA OF ANY
//     KIND. This is a plasma-only TDM paper. It contributes nothing to the
//     bone-penetration-ratio problem and nothing to the unsourced 8-256x
//     biofilm MIC multiplier. Any bone concentration derived from this file is
//     a plasma concentration multiplied by a number that came from somewhere
//     else.
//
// 15. INDICATION MISMATCH. Gentamicin here was given for PROPHYLAXIS AND
//     TREATMENT OF GRAM-NEGATIVE INFECTION in cancer patients. It is not a
//     bone-infection, implant, or orthopaedic population, and there is no
//     prosthesis, no surgical washout physiology and no long course in the
//     data.
//
// 16. 1998 VINTAGE. Fluorescence polarization immunoassay (TDx, LOQ 0.1 mg/L,
//     interassay CV 6.3% at 1 mg/L, 3.7% at 4, 4.3% at 8). NONMEM version IV
//     -- the paper says "version IV" and nothing more; the ESTIMATION METHOD
//     (FO vs FOCE vs FOCE-I) IS NEVER STATED, and FO-era omegas from sparse
//     TDM data are routinely biased. A modern reanalysis of the same data
//     would very likely give different variance components.
//
// 17. NO INTER-OCCASION VARIABILITY, although three patients contributed
//     several separate courses of treatment (14, 17 and 22 measurements each)
//     pooled without an occasion random effect. Some of what is attributed to
//     residual error is really between-occasion variability, which means the
//     residual sigma below is inflated relative to true assay-plus-model
//     error, and the CL omega is correspondingly deflated.
//
// 18. "NON-ICU, NON-DIALYSIS" IS AN ARGUMENT FROM SILENCE. The paper states
//     only two exclusions: missing clinical/dosage/sampling data, and rapidly
//     changing renal function. It never characterises the care setting and
//     never mentions dialysis. The exclusion of unstable renal function is
//     what actually rules out most AKI/ICU physiology, and CLcr does reach
//     15.8-22.9 mL/min at the low end. Catalogue it honestly.
//
// 19. GENDER IS NOT A COVARIATE IN THIS MODEL. It enters ONLY inside the
//     Cockcroft-Gault equation, via the 0.85 female factor that the paper
//     itself calls arbitrary. Harness convention SEX 0 = MALE, 1 = FEMALE is
//     used below; the 0.85 applies to SEX == 1. A transposed SEX coding here
//     changes CLcr by 15% in both directions.
//
// 20. STRATEGIC. This is the FIFTH gentamicin file in a library that already
//     codes the aminoglycoside class (gomes, evers_standard, evers_icu,
//     franck). Per the project's own measurement that the worst failure mode
//     is an unseen drug CLASS (R2 0.43 +/- 0.13), this adds less than a
//     linezolid oxazolidinone model would. Its value is that it is a
//     structurally different (two-compartment, albumin-on-volume) partner to
//     the one-compartment Gomes/Evers family, so it is a genuine independent
//     structure rather than another re-estimation of the same one.
//
// 21. THIS MODEL AND genta_rosario_clinical.cpp ARE THE SAME DATA. They are
//     co-equal alternative covariate parameterisations refitted on the same
//     210 patients, not independent evidence. The paper: "Comparison of the
//     two best models (model 8a and 19b) ... allowed no firm conclusion to be
//     drawn"; the CLcr version was picked "due to its simplicity". Never treat
//     agreement between them as external corroboration, and never hold one out
//     to validate the other.
//
// SMOKE EXPECTATION
//   Standard adult regimen: gentamicin 5 mg/kg q24h = 400 mg over 30 min,
//   80 kg, 175 cm, 65 y, MALE, CREAT 80 umol/L, ALB 40 g/L.
//     -> BSA (Du Bois) 1.96 m2; CLcr (CG, TBW) 92.1 mL/min (creatinine 80 is
//        above the 60 umol/L floor, so no truncation applies)
//     -> CL 4.36 L/h, V1 15.77 L, Q 1.30 L/h, V2 9.79 L, Vss 25.6 L
//     -> t1/2,alpha 1.73 h, t1/2,z 7.55 h
//     -> FIRST-DOSE PEAK at end of infusion (t = 0.5 h) 23.2 mg/L
//     -> 24 h TROUGH 0.38 mg/L (steady state ~0.43 mg/L; the q24h
//        accumulation factor on the beta phase is only 1.12)
//   IN RANGE: the target band is peak 15-25 mg/L and 24 h trough <1 mg/L, and
//   this model gives 23.2 mg/L and 0.38 mg/L. It sits in the UPPER half of the
//   peak band because Vss 0.38 L/kg applies to a 66 kg cohort median while the
//   smoke patient is 80 kg -- V1 scales on BSA, not on weight, so V1/kg falls
//   as weight rises and the peak is pushed up. Expect this file to give HIGHER
//   peaks than the one-compartment Gomes/Evers gentamicin files at the same
//   dose (genta_gomes gives 17.3 mg/L for a comparable patient), which is a
//   real structural disagreement between published models, not a defect.
//   Cross-checks against the paper's own numbers, which should be reproduced
//   exactly by a correct build:
//     - at the truncated cohort median CLcr 86.1 mL/min, CL = 4.14 L/h
//       (Table 5 derived value 4.2 L/h)
//     - at the paper's own simulated reference patient (Dosage guidelines /
//       Figure 2: 50 y, 70 kg, 170 cm, albumin 40 g/L) BSA = 1.81 m2,
//       V1 = 14.59 L, Vss = 14.59 + 9.79 = 24.4 L (paper 24.6 L, 0.38 L/kg)
//     - at CLcr 105 mL/min for that reference patient, CL = 4.85 L/h; the
//       paper simulated exactly CLcr 25 and 105 mL/min for Figure 2, so that
//       figure is a ready-made visual regression test for this file.
// ====================================================================

[SET] end=72, delta=0.1


$PARAM @annotated
TVCLINT : 0.88  : CL intercept, theta1 (L/h)
TVCLSLP : 0.043 : slope of CL on Cockcroft-Gault CLcr60, theta5 (per mL/min)
TVV1BSA : 8.59  : V1 coefficient on body surface area, theta2 (L/m2)
ALBEXP  : -0.39 : exponent of albumin on V1, theta6 (dimensionless)
TVQ     : 1.30  : intercompartmental clearance, theta3 (L/h)
TVV2    : 9.79  : peripheral volume, theta4 (L)

ALBREF  : 34.0  : albumin centering value (g/L) -- study median, below normal
CRFLOOR : 60.0  : creatinine floor before Cockcroft-Gault (umol/L) -- 60 HERE

ETA1 : 0 : IIV on clearance
ETA2 : 0 : IIV on intercompartmental clearance
ETA3 : 0 : IIV on central volume -- NOT PUBLISHED, omega fixed to 0, LIMITS 7

// CRFLOOR is 60 in THIS model. The clinical-factors model floors at 70 AND
// divides by 71. Do not mix the two sets of constants. See LIMITS 6.


$PARAM @annotated @covariates
WT    : 70.0  : total body weight (kg) -- paper's own reference patient
HT    : 170.0 : height (cm)
AGE   : 50.0  : age (years)
SEX   : 0     : 0 = MALE, 1 = FEMALE -- used ONLY inside Cockcroft-Gault
CREAT : 71.0  : serum creatinine (umol/L) -- cohort median
ALB   : 40.0  : serum albumin (g/L) -- reference patient value

// Defaults above are the paper's OWN simulated reference patient (Dosage
// guidelines / Figure 2: 50 y, 70 kg, 170 cm, albumin 40 g/L), not the cohort
// medians (66 kg, 170 cm, 50 y, ALB 34 g/L, CREAT 71 umol/L), so that a bare
// run of this file reproduces a figure that exists in the source. SEX is not
// a covariate of the model; it enters only via the arbitrary 0.85 female
// factor of Cockcroft-Gault.


// OMEGA. Table 5 rows are "CV omega_CL (%)" = 18.5 and "CV omega_Q (%)" =
// 28.2, and the Abstract confirms these are COEFFICIENTS OF VARIATION, not
// variances, against a Table 5 footnote that says otherwise. Lognormal
// random effects (Methods: Pj = P*Exp(eta_pj)), so omega^2 = ln(1 + CV^2):
//   ETA1 CL : ln(1 + 0.185^2) = ln(1.034225) = 0.033652
//   ETA2 Q  : ln(1 + 0.282^2) = ln(1.079524) = 0.076520
//   ETA3 V1 : NOT REPORTED IN THE FINAL MODEL -> 0. The paper's covariate-free
//             base model offers 15.5% CV on V1, i.e. ln(1+0.155^2) = 0.023741,
//             but that belongs to a different model on a different subset.
//             Substitute it only deliberately and record it. See LIMITS 7.
// Omega assumed DIAGONAL; the paper prints no off-diagonals and does not say.
$OMEGA 0.033652 0.076520 0.000000


// Combined residual error as published: Y = C*exp(e1) + e2, variances
// error1 = 0.025 (dimensionless, exponential) and error2 = 0.107 (mg/L)^2.
// These are the printed three-decimal values; the trailing digit the PDF text
// layer appends to each is a superscript footnote marker.
$SIGMA 0.025 0.107


$CMT @annotated
CENT   : Central compartment (mg) [ADM, OBS]
PERIPH : Peripheral compartment (mg)


$MAIN
// --- body size -------------------------------------------------------
// Du Bois body surface area (the paper's reference 16 is Du Bois & Du Bois,
// "Clinical calorimetry"), WT in kg and HT in cm.
double BSA = 0.007184 * pow(WT, 0.425) * pow(HT, 0.725);
if(BSA < 0.5) BSA = 0.5;      // guard only; observed range was 1.3-2.3 m2
if(BSA > 3.5) BSA = 3.5;

// --- albumin ---------------------------------------------------------
// The harness carries ALB without a declared unit. This model needs g/L
// (centred on 34). A value below 10 can only be g/dL, so it is converted here
// rather than silently mis-scaled -- a g/dL value of 4.0 taken as g/L would
// inflate V1 by (34/4)^0.39 = 2.3-fold. The value actually used is captured.
double ALBGL = (ALB < 10.0) ? ALB * 10.0 : ALB;
if(ALBGL < 10.0) ALBGL = 10.0;   // observed albumin range was 17-53 g/L
if(ALBGL > 60.0) ALBGL = 60.0;

// --- renal function --------------------------------------------------
// Creatinine floor: "All creatinine concentrations less than or equal to
// 60 umol/L were set to 60 umol/L". Applied BEFORE Cockcroft-Gault.
double CREAT60 = (CREAT < CRFLOOR) ? CRFLOOR : CREAT;

// Cockcroft-Gault. CREAT is umol/L in the harness and in the paper;
// 1 mg/dL = 88.4 umol/L. The body-weight descriptor is NOT stated by the
// paper -- TOTAL body weight is assumed here. See LIMITS 4.
// CLcr stays in mL/MIN because the slope 0.043 is per mL/min. See LIMITS 5.
double CRMGDL = CREAT60 / 88.4;
double CLCR = (140.0 - AGE) * WT / (72.0 * CRMGDL) * ((SEX == 0) ? 1.0 : 0.85);

// --- parameters ------------------------------------------------------
double CL = TVCLINT * (1.0 + TVCLSLP * CLCR) * exp(ETA1 + ETA(1));
double Q  = TVQ * exp(ETA2 + ETA(2));
double V1 = TVV1BSA * BSA * pow(ALBGL / ALBREF, ALBEXP) * exp(ETA3 + ETA(3));
double V2 = TVV2;


$ODE
dxdt_CENT   = - (CL / V1) * CENT - (Q / V1) * CENT + (Q / V2) * PERIPH;
dxdt_PERIPH =   (Q / V1) * CENT - (Q / V2) * PERIPH;


$TABLE
// Combined error exactly as published: Y = C * exp(e1) + e2.
// The exponential form is the paper's, not the house (1+EPS(1)) default.
double CP = CENT / V1;
double DV = CP * exp(EPS(1)) + EPS(2);
int i = 0;
while(DV<0 && i <100) {
simeps();
DV = CP * exp(EPS(1)) + EPS(2);
++i;
}


$CAPTURE DV CP CL V1 V2 Q BSA CLCR ALBGL
