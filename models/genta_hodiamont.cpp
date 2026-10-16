// =========================================================================
// genta_hodiamont.cpp -- Gentamicin, adult critically ill (Hodiamont 2017)
//
// SOURCE
//   Hodiamont CJ, Juffermans NP, Bouman CSC, de Jong MD, Mathot RAA,
//   van Hest RM. "Determinants of gentamicin concentrations in critically ill
//   patients: a population pharmacokinetic analysis."
//   Int J Antimicrob Agents 2017;49:204-211.
//   PDF: Determinants of gentamicin concentrations in
//        critically ill patients Hodiamont 2017.pdf
//
// POPULATION
//   44 adult ICU patients, Academic Medical Center Amsterdam, Jan-Jun 2011.
//   303 concentrations, 174 doses. Age median 61 y (20-78). TBW median 70.5 kg
//   (42-116). IBW median 68.2 kg (55.6-87.5). Serum albumin median 21.5 g/L
//   (range 10-36) -- a HYPOALBUMINAEMIC population. Measured 6 h urinary
//   creatinine clearance median 48.3 mL/min (0-130). 5/44 patients on CVVH.
//   Starting dose 4 mg/kg total body weight over 30 min.
//
// WHY THIS MODEL IS IN THE LIBRARY
//   It is the SECOND STRUCTURE for gentamicin. genta_gomes and the two Evers
//   variants share one structure and one covariate algebra, so none of them
//   validates another. This model is two-compartment, allometric on IDEAL body
//   weight, and carries albumin on V1 -- structurally independent of all three.
//   It is what lets gentamicin enter leave-one-model-out external validation.
//
// STRUCTURE
//   Two compartments, first-order elimination, 30 min IV infusion.
//   All disposition scaled to 70 kg IDEAL body weight.
//     CL(off CVVH) = 1.15 * (IBW/70)^0.75
//                          * (1 + 0.0131*(CLCR - 30)*FLAG) * 1.39^(1-FLAG)
//     CL(on CVVH)  = 2.13 * (IBW/70)^0.75
//     V1           = 21.2 * (IBW/70)^1 * (ALB/22)^-0.833
//     Q            = 1.96 * (IBW/70)^0.75
//     V2           = 18.4 * (IBW/70)^1
//
// PARAMETERS as printed (Table 2), with RSE
//   CLnoCVVH 1.15 (14.9%)   CLCVVH 2.13 (7.9%)   V1 21.2 (5.4%)
//   Q 1.96 (20.7%)          V2 18.4 (9.5%)
//   CLCR slope 0.0131/mL/min centred at 30      ALB exponent -0.833 (21%)
//   FLAG multiplier 1.39 (16%)
//   Table 2 prints the CLCR slope as 0.0132; 0.0131 is the value used in
//   Eq. 4 and in the running text, and it is the one that reproduces the
//   published worked examples, so it is the one implemented.
//
//   Four independent arithmetic checks reproduce from these values:
//     V1 = 16.37 vs printed 16.3 L at albumin 30 g/L
//     V1 = 29.17 vs printed 29.2 L at albumin 15 g/L
//     CL = 1.903 vs printed 1.90 L/h at CLCR 80 mL/min
//     CL = 0.924 vs printed 0.92 L/h at CLCR 15 mL/min
//
// IIV reported as CV%, lognormal. The Table 2 footnote defines CV% as
//   sqrt(exp(omega^2)-1)*100 (the published PDF lost the superscript 2).
//   That reading is confirmed independently: the paper states the covariates
//   explained 36% of IIV in V1 and 64% in CLnoCVVH, and this transform
//   reproduces 36.1% and 64.0%.
//     CLnoCVVH 42.5 CV% -> omega^2 0.1660
//     CLCVVH   29.5 CV% -> omega^2 0.0834
//     V1       17.2 CV% -> omega^2 0.0292
//     corr(CLnoCVVH, V1) = 0.54 -> covariance 0.0376
//
// -------------------------------------------------------------------------
// LIMITS
//
//   1. THE RENAL COVARIATE IS NOT COCKCROFT-GAULT. CLCR here is a MEASURED
//      creatinine clearance from a timed 6 h urine collection
//      (00:00-06:00), not an estimate from serum creatinine. The harness
//      supplies only serum CREAT, so Cockcroft-Gault is substituted in $MAIN.
//      That substitution is a documented approximation, and the paper tests it
//      directly: in its own visual predictive checks, serum creatinine and
//      Cockcroft-Gault both OVERESTIMATE clearance and so UNDERESTIMATE
//      troughs -- but ONLY beyond 36 h post-dose (observed 1.47 vs simulated
//      0.79 mg/L with SCr; 1.19 vs 0.80 with Cockcroft-Gault). In the 23-36 h
//      window all four renal measures are essentially unbiased (SCr 1.91 vs
//      1.87; Cockcroft-Gault 1.98 vs 1.89; measured 1.91 vs 2.09). A once-daily
//      regimen in a patient with preserved renal function lives in the second
//      regime, so the substitution is acceptable there and NOT acceptable for
//      extended intervals in severe renal impairment. CLCR is exposed as a
//      covariate so a measured value can be supplied instead when one exists.
//
//   2. ALBUMIN IS EXTRAPOLATED FOR MOST PJI PATIENTS. The model was built at
//      median albumin 21.5 g/L, observed range 10-36. Only 42.7% of the
//      simulated SAFEBONE cohort falls inside that range; the simulated median
//      is 37 g/L, just above the top of it. A power term with exponent -0.833
//      applied well outside its data is exactly how a model produces confident
//      nonsense, so ALB is CLAMPED to the observed range by default. The bounds
//      are exposed as ALBLO and ALBHI; set ALBHI high to disable the clamp and
//      reproduce the published equation verbatim. Clamping at 36 rather than
//      extrapolating to 45 raises V1 by about 20%.
//
//   3. THE IDEAL BODY WEIGHT EQUATION IS NOT PRINTED. The paper cites only
//      Devine (ref 21). Devine is implemented here, which also keeps this file
//      consistent with genta_gomes and the Evers variants. Note that the
//      reported IBW range (55.6-87.5 kg) does not reconcile with the reported
//      height range (154-195 cm) under Devine at the lower end, though that is
//      explicable if the two ranges are marginal rather than paired.
//
//   4. Q AND V2 EXPONENTS ARE NOT PRINTED. Only the Table 2 units ("L/h/70 kg",
//      "L/70 kg") and the statement that parameters were allometrically scaled
//      imply them. 0.75 for Q and 1 for V2 are assumed by analogy with CL and
//      V1 and with the cited Anderson and Holford framework. ASSUMPTION.
//
//   5. THE CVVH BRANCH IS NOT DRIVEN BY THE HARNESS. CVVH is continuous
//      veno-venous haemofiltration, which is neither DIAL (chronic dialysis)
//      nor IHD (intermittent haemodialysis session). The harness supplies no
//      CVVH indicator, so this model always runs the off-CVVH branch -- which
//      is the correct branch for a prosthetic-joint-infection ward cohort. The
//      authors warn in any case that CLCVVH is specific to their filter and
//      flow settings and is likely not applicable to other modalities.
//
//   6. NO INTEROCCASION VARIABILITY. IOV on CLnoCVVH was estimated in the basic
//      model (38.6 CV%) and REJECTED from the final one. Adding it back would
//      introduce variability the authors deliberately removed. The $OMEGA block
//      below contains none.
//
//   7. THE CL-V1 CORRELATION IS POORLY IDENTIFIED. r = 0.54 carries RSE 91%,
//      above the 80% acceptance threshold the authors set for random effects,
//      and its bootstrap confidence interval spans -1 to 0.69. It is
//      implemented because it is part of the published final model, not
//      because it is well estimated.
//
//   8. CRITICALLY ILL, NOT ORTHOPAEDIC. Capillary leak and augmented renal
//      clearance are ICU phenomena. Applied to a stable ward patient with
//      preserved renal function this model gives a clearance around 2.1 L/h
//      where the Gomes family gives 4.2-5.3 L/h -- a two-fold disagreement.
//      That is within the range this library shows between published models of
//      the same drug (rifampicin spans 5.0-fold, vancomycin 3.2-fold), but the
//      direction is systematic: the linear renal term, centred at 30 mL/min in
//      a population whose median was 48, under-extrapolates at 90 mL/min.
//
//   SMOKE EXPECTATION
//     80 kg / 175 cm / 65 y male, CREAT 80 umol/L, ALB 40 (clamped to 36),
//     400 mg over 30 min q24h. MEASURED at steady state:
//       CL 2.10 L/h, V1 14.2 L, Vss 32.7 L, peak 29.2 mg/L, 24 h trough
//       2.95 mg/L, AUC24 191 mg*h/L.
//     Compare the same patient and regimen under the other three gentamicin
//     models: peak 15.3-18.8, trough 0.10-0.24, AUC24 76-96. The 2.0-2.5 fold
//     exposure gap and the trough sitting ABOVE the 1 mg/L re-dosing threshold
//     are both consequences of limit 8 and are expected, not a transcription
//     error. What WOULD indicate a transcription error is this model landing
//     in the 15-20 mg/L peak band alongside the others.
// =========================================================================

[SET] end=72, delta=0.1


$PARAM @annotated
TVCL0 : 1.15 : Clearance off CVVH at 70 kg IBW and CLCR 30 mL/min (L/h)
TVCLH : 2.13 : Clearance on CVVH at 70 kg IBW (L/h)
TVV1  : 21.2 : Central volume at 70 kg IBW and albumin 22 g/L (L)
TVQ   : 1.96 : Intercompartmental clearance at 70 kg IBW (L/h)
TVV2  : 18.4 : Peripheral volume at 70 kg IBW (L)
SLCLCR : 0.0131 : Linear effect of CLCR on clearance (per mL/min)
CTCLCR : 30.0 : Centring value for CLCR (mL/min)
EXALB : -0.833 : Power exponent of albumin on V1
CTALB : 22.0 : Centring value for albumin (g/L)
FMISS : 1.39 : Multiplier applied when CLCR is unavailable (FLAG = 0)
EXCL  : 0.75 : Allometric exponent on IBW for CL and Q
EXV   : 1.0 : Allometric exponent on IBW for V1 and V2
ALBLO : 10.0 : Lower clamp on albumin, the study minimum -- see LIMITS 2
ALBHI : 36.0 : Upper clamp on albumin, the study maximum -- see LIMITS 2

ETA1 : 0 : IIV on clearance off CVVH
ETA2 : 0 : IIV on clearance on CVVH
ETA3 : 0 : IIV on central volume


$PARAM @annotated @covariates
WT    : 70.5 : total body weight (kg)
HT    : 170.0 : height (cm)
AGE   : 61.0 : age (years)
SEX   : 0 : 0 = male, 1 = female
CREAT : 115.0 : serum creatinine (umol/L)
ALB   : 21.5 : serum albumin (g/L)
CVVH  : 0 : continuous veno-venous haemofiltration in progress -- see LIMITS 5
FLAG  : 1 : 1 = CLCR known, 0 = unknown -- see LIMITS 1


// Lower triangle of the omega block. CV% converts as omega^2 = ln(1+CV^2).
//   ETA1 CLnoCVVH : 42.5 CV% -> 0.1660
//   ETA2 CLCVVH   : 29.5 CV% -> 0.0834   (uncorrelated with the others)
//   ETA3 V1       : 17.2 CV% -> 0.0292
//   cov(ETA1,ETA3) = 0.54 * sqrt(0.1660*0.0292) = 0.0376   -- see LIMITS 7
$OMEGA @block
0.1660
0.0000 0.0834
0.0376 0.0000 0.0292


$SIGMA
0.1142   // proportional, 33.8% CV as published; sigma^2 = 0.338^2
0        // additive residual variance fixed to zero


$CMT @annotated
CENT   : Central compartment (mg) [ADM, OBS]
PERIPH : Peripheral compartment (mg)


$MAIN
// Ideal body weight, Devine. Cited by the paper but never printed -- LIMITS 3.
double IBW = (SEX == 0) ? 50.0 + 2.3*(HT/2.54 - 60.0)
                        : 45.5 + 2.3*(HT/2.54 - 60.0);
if (IBW < 1.0) IBW = WT;                 // guard: very short stature
double FSZCL = pow(IBW/70.0, EXCL);
double FSZV  = pow(IBW/70.0, EXV);

// Cockcroft-Gault standing in for the measured 6 h urinary clearance -- LIMITS 1.
double CRMGDL = CREAT / 88.4;            // 1 mg/dL = 88.4 umol/L
double CLCR = (140.0 - AGE) * WT / (72.0 * CRMGDL) * ((SEX == 0) ? 1.0 : 0.85);

// Albumin held inside the range the model was estimated over -- LIMITS 2.
double ALBC = ALB;
if (ALBC < ALBLO) ALBC = ALBLO;
if (ALBC > ALBHI) ALBC = ALBHI;

double CLREN = 1.0 + SLCLCR*(CLCR - CTCLCR)*FLAG;
if (CLREN < 0.05) CLREN = 0.05;          // guard: the linear term can go negative
double CL0 = TVCL0 * FSZCL * CLREN * pow(FMISS, 1.0 - FLAG) * exp(ETA1 + ETA(1));
double CLH = TVCLH * FSZCL * exp(ETA2 + ETA(2));
double CL  = (CVVH > 0.5) ? CLH : CL0;

double V1 = TVV1 * FSZV * pow(ALBC/CTALB, EXALB) * exp(ETA3 + ETA(3));
double Q  = TVQ  * FSZCL;
double V2 = TVV2 * FSZV;


$ODE
dxdt_CENT   = - (CL/V1)*CENT - (Q/V1)*CENT + (Q/V2)*PERIPH;
dxdt_PERIPH =                  (Q/V1)*CENT - (Q/V2)*PERIPH;


$TABLE
double CP = CENT / V1;
double DV = CP * (1 + EPS(1));
int i = 0;
while(DV<0 && i <100) {
simeps();
DV = CP * (1 + EPS(1));
++i;
}


$CAPTURE DV CP CL V1 Q V2 CLCR IBW ALBC
