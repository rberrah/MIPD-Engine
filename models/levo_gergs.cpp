// ===========================================================================
// levo_gergs.cpp -- Levofloxacin, adult hip/knee arthroplasty patients
//
// SOURCE
//   Gergs U. et al., "Population Pharmacokinetics of Levofloxacin in Plasma
//   and Bone of Patients Undergoing Hip or Knee Surgery".
//   PDF: Gergs - Population Pharmacokinetics of
//        Levofloxacin in Plasma and Bone of Patients.pdf
//
// POPULATION  42 adult patients undergoing elective hip or knee arthroplasty,
//   age 46-84 y. Single 500 mg IV infusion (60 min) as surgical prophylaxis.
//   This is the ONLY model in the library besides rifampicine_marsot that was
//   developed in an orthopaedic-surgical population.
//
// STRUCTURE   One compartment, IV infusion input, linear elimination.
//   CL and V are reported PER KILOGRAM, so weight enters as a simple linear
//   scaling -- not allometric, no reference weight, no centering constant.
//
// PARAMETERS as printed (Table 2)
//   CL 0.164 L/h/kg (10% RSE)   interpatient CV 32% (32% RSE)
//   V  0.835 L/kg   (12% RSE)   interpatient CV 51% (18% RSE)
//   Variability is reported as % CV on a lognormal scale, so
//   omega^2 = ln(1 + CV^2):  CL -> ln(1.1024) = 0.0975
//                            V  -> ln(1.2601) = 0.2312
//
// ---------------------------------------------------------------------------
// LIMITS -- read before using this model for anything but external validation
//
//   1. NO RENAL COVARIATE. Elevated creatinine was an EXCLUSION CRITERION, so
//      the study contains no renal impairment and the model has no mechanism
//      for it. CREAT/CREAT2 supplied by the harness bind to nothing here. Do
//      not apply to impaired renal function.
//
//   2. THE FITTED HALF-LIFE IS TOO SHORT. ln2 * V/CL = 0.693 * 0.835/0.164
//      = 3.5 h, against the 8-13 h terminal half-life the paper itself cites
//      for levofloxacin. Sampling stopped before 5 h and a one-compartment
//      model was forced onto it. Any AUC24 or fAUC/MIC derived from this model
//      is therefore biased LOW. This is the dominant quantitative hazard and
//      the reason the model is annotated rather than silently trusted.
//
//   3. RESIDUAL ERROR IS ASSUMED, NOT PUBLISHED. The paper describes a
//      combined error model but never prints sigma0 or sigma1. The $SIGMA
//      block below is an ASSUMPTION (10% proportional, small additive),
//      chosen to be typical of the reverse-phase HPLC assay described in the
//      Methods. It is NOT a published estimate. Do not quote it as one.
//
//   4. NO BONE COMPARTMENT IMPLEMENTED. The paper fits an irreversible bone
//      uptake compartment (penetration t1/2 4.2 h cortical, 5.4 h trabecular).
//      It is deliberately omitted: an irreversible-uptake compartment cannot
//      represent a decline phase and cannot be used for repeated dosing, so
//      wiring it in would produce bone concentrations that only ever rise.
//      The published bone:plasma ratios (3.8 at 1 h rising to 20.9 at 3 h) are
//      a consequence of that structure, not a measurement, and must not be
//      used as penetration factors.
//
//   5. Single dose, no steady state, uninfected prophylaxis cohort
//      (osteomyelitis was excluded). Mean body weight is never reported: the
//      abstract's 14.0 L/h implies ~85 kg while its 77 L implies ~92 kg.
//
//   PAIRING NOTE. levo_canoui is an ORAL model with an absorption constant;
//   this one is IV-only. The leave-one-model-out comparison between them is a
//   disposition comparison, and the two must not be cross-labelled by route.
// ===========================================================================

[SET] end=48, delta=0.1


$PARAM @annotated
TVCL : 0.164 : Clearance per kg (L/h/kg)
TVV  : 0.835 : Volume of distribution per kg (L/kg)

ETA1 : 0 : IIV on clearance
ETA2 : 0 : IIV on volume of distribution


$PARAM @annotated @covariates
WT : 70 : body weight (kg)


$OMEGA 0.0975 0.2312   // ln(1+CV^2) for CV = 32% and 51%


$SIGMA
0.01   // proportional -- ASSUMED, not published (see LIMITS 3)
0.0001 // additive     -- ASSUMED, not published (see LIMITS 3)


$CMT @annotated
CENT : Central compartment (mg) [ADM, OBS]


$MAIN
// Weight enters linearly because the paper reports CL and V per kilogram.
double CL = TVCL * WT * exp(ETA1 + ETA(1));
double V  = TVV  * WT * exp(ETA2 + ETA(2));


$ODE
dxdt_CENT = - (CL / V) * CENT;


$TABLE
double DV = (CENT/V) * (1 + EPS(1)) + EPS(2);
int i = 0;
while(DV<0 && i <100) {
simeps();
DV = (CENT/V) * (1 + EPS(1)) + EPS(2);
++i;
}


$CAPTURE DV CL V
