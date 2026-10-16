// =========================================================================
// linez_buerger.cpp -- Linezolid, adult ICU severe sepsis / septic shock
//                      (Buerger 2006), unbound-plasma two-compartment model
//
// SOURCE
//   Buerger C, Plock N, Dehghanyar P, Joukhadar C, Kloft C. "Pharmacokinetics
//   of unbound linezolid in plasma and tissue interstitium of critically ill
//   patients after multiple dosing using microdialysis."
//   Antimicrob Agents Chemother 2006 Jul;50(7):2455-2463.
//   doi:10.1128/AAC.01468-05
//   PDF: harmacokinetics of Unbound Linezolid in
//        Plasma and Tissue Buerger.pdf
//   Read here with `pdftotext -layout` (/mingw64/bin, 536-line dump, all 9
//   pages). NOTE FOR THE RECORD: `pdftoppm` is genuinely absent from this
//   machine, so no page could be raster-rendered locally; every number below
//   was taken from the text layer and cross-checked against the printed
//   internal consistency relation 600/9.18 = 65.36 vs printed AUC 65.3.
//   Table 1 reads verbatim: "62.9 (19.2) 9.18 (57.9) 62.1 (87.9) 65.3 (57.9)
//   5.0 (45.5)", footnote "Values in parentheses are cv (%) geometric means
//   (n = 10)."
//
// POPULATION
//   12 adult ICU patients enrolled (9 male, 3 female), Medical University of
//   Vienna. 11/12 (91.7%) SEPTIC SHOCK, 1/12 (8.3%) SEVERE sepsis. All
//   sedated and mechanically ventilated. Median age 62 y (range 51-74).
//   Median weight 81 kg (range 55-133). 600 mg IV over 30 min q12h.
//   These parameters are the STEADY-STATE (visit 2) fit, n = 10.
//   RENAL FUNCTION IS NOT REPORTED AT ALL -- no creatinine, no CrCl/eGFR, no
//   urine output, no dialysis status, no albumin, no APACHE/SOFA. Confirmed
//   by grep over the full text dump (creatinin = 0 hits, renal = 0, albumin
//   = 0, APACHE = 0, SOFA = 0). NOT a bone or prosthetic joint population.
//
// STRUCTURE
//   Open two-compartment model. Zero-order IV infusion into the central
//   compartment, first-order elimination from the central compartment.
//   INDIVIDUAL two-stage curve fitting in WinNonlin v4.0 (Gauss-Newton,
//   1/Cpred^2 weighting). THIS IS NOT A POPULATION PK MODEL -- see LIMITS 1.
//
// PARAMETERS AS PRINTED (Table 1, steady state, n = 10, geometric means)
//   Vss  62.9 L      (cv 19.2%)   = V1 + V2, the split is NOT printed
//   CL    9.18 L/h   (cv 57.9%)
//   CLD2 62.1 L/h    (cv 87.9%)   intercompartmental clearance
//   fAUC0-12  65.3 mg*h/L (cv 57.9%)     t1/2 5.0 h (cv 45.5%)
//   Cmax,ss 16.4 mg/L (cv 30.4%, n=9)    Cmin,ss 1.83 mg/L (cv 188%, n=9)
//   fu   0.866 (cv 7.9%, n=22 determinations, individual range 0.730-0.959)
//   ALL OF THESE ARE ON AN UNBOUND BASIS. Methods: "the compartmental data
//   analysis was based on unbound concentrations for all matrices." So the
//   printed CL is CL_total/fu and the printed Vss is Vss_total/fu. Absolute
//   litres and litres/h -- NOT weight-normalised, NOT per 70 kg. IV data, so
//   this is CL, not CL/F.
//
//   V1/V2 SPLIT -- DERIVED HERE, NOT PUBLISHED (see LIMITS 2).
//   The paper prints only Vss. I solved the two-compartment ODE with the
//   printed CL = 9.18, CLD2 = 62.1 and Vss = 62.9 fixed, sweeping V1:
//       V1 =  5 L -> t1/2 5.35 h, Cmax,ss 24.0 mg/L
//       V1 = 21 L -> t1/2 5.07 h, Cmax,ss 17.8 mg/L
//       V1 = 30 L -> t1/2 4.95 h, Cmax,ss 15.2 mg/L
//       V1 = 60 L -> t1/2 4.75 h, Cmax,ss 11.2 mg/L
//   AUCtau is 65.36 mg*h/L for EVERY V1 (it depends on CL alone) and matches
//   the printed 65.3, which confirms the printed set is internally consistent
//   on exposure. t1/2 is nearly insensitive to V1 because CLD2 >> CL, so
//   rounding to "5.0 h" only brackets V1 to roughly 21-30 L. Matching the
//   printed steady-state Cmax of 16.4 mg/L pins it:
//       V1 = 25.33 L, V2 = 37.57 L  ->  Cmax,ss 16.40, t1/2 5.01 h.
//   That same set predicts Cmin,ss 2.08 mg/L against a printed 1.83 mg/L, a
//   14% mismatch. The mismatch is the PAPER'S, not this file's: Vss, CL,
//   Cmax and Cmin are separately-computed geometric means over different n
//   (n=10 and n=9) and do not form one internally consistent parameter
//   vector. The split is exposed as FRV1 so it can be changed in one place.
//
//   FREE-VS-TOTAL CONVERSION. The rest of this library simulates TOTAL
//   plasma concentrations, so DV here is TOTAL drug. The conversion is done
//   explicitly in $MAIN, never by editing the printed numbers: scaling CL,
//   CLD2, V1 and V2 all by fu leaves every rate constant CL/V unchanged and
//   simply reinterprets the amounts as total drug, giving
//       CL_tot = 9.18 * 0.866 = 7.950 L/h,  Vss_tot = 62.9 * 0.866 = 54.5 L.
//   Wiring 9.18 L/h in as a TOTAL clearance would overstate elimination by
//   9.18/7.95 = 1.155, i.e. ~15%, and understate total AUC by the same
//   factor. That is the CL-vs-CL/F class of error in a free-vs-total guise.
//   Set FU = 1 to recover the paper's native unbound-basis simulation; the
//   unbound concentration is captured as CU regardless.
//
// -------------------------------------------------------------------------
// LIMITS -- every reason this model can mislead in an adult PJI cohort
//
//   1. THIS IS NOT A POPULATION PK MODEL, AND IT DOES NOT CLOSE THE
//      OXAZOLIDINONE GAP. It is an individual (two-stage) WinNonlin analysis:
//      no covariate model, no estimated omega, no residual error model, no
//      NONMEM/nlmixr fit. The authors' own closing sentence concedes it --
//      "Modeling the data by using a population pharmacokinetic approach may
//      help to assess the different types of variability." The library still
//      needs a genuine linezolid PopPK source (Sasaki, already in the
//      Reference folder) before the unseen-drug-class weakness is closed.
//      Treat this file as a DETERMINISTIC TYPICAL-VALUE EXPOSURE BENCHMARK
//      for linezolid in septic adults, not as a population simulator.
//
//   2. THE V1/V2 SPLIT IS NOT PUBLISHED. It is my back-calculation from the
//      printed Vss/CL/CLD2/t1/2/Cmax, documented in full above. AUC24 -- the
//      primary SafeBone output -- is immune to it, because AUC24 = Dose24/CL
//      and CL is printed. Cmax is NOT immune: it is the very quantity used to
//      pin the split, so it is fitted, not predicted, and must never be
//      quoted from this file as independent corroboration of the paper.
//      Anything shape-dependent (fT>MIC, distribution half-life, early-phase
//      concentrations) inherits this assumption.
//
//   3. THE $OMEGA VALUES ARE NOT POPULATION OMEGAS. The bracketed figures in
//      Table 1 are, per the table's own footnote, "cv (%) geometric means" of
//      the INDIVIDUAL WinNonlin estimates. They conflate true between-patient
//      variability with parameter estimation error, they are not
//      shrinkage-corrected, and n is only 10. Worse, the paper never states
//      whether that cv was computed as SD/mean of the raw values or as
//      sqrt(exp(sd_log^2)-1), so even the CV-to-omega arithmetic is
//      undefined. Using them as IIV is an ASSUMPTION, not a transcription.
//      They are almost certainly too WIDE for true IIV (estimation error is
//      folded in) while simultaneously resting on far too few patients.
//
//   4. NO COVARIATES OF ANY KIND -- and the weight range was 55-133 kg. The
//      parameters are absolute litres and litres/h with no allometric term
//      and no centering constant. WT is declared below and DELIBERATELY LEFT
//      UNUSED: adding a (WT/70)^0.75 term would be invention, not
//      transcription. The practical consequence is that this model returns
//      the SAME exposure for a 55 kg and a 133 kg patient. For a PJI cohort
//      spanning a normal adult weight range that is a real and unquantified
//      bias, and it is the single biggest reason not to use this file for
//      individual dose selection.
//
//   5. RENAL FUNCTION IS ENTIRELY UNREPORTED, SO THE MODEL CANNOT RESPOND TO
//      IT. None of CREAT, CREAT2, DIAL or IHD can be linked to these
//      estimates. Linezolid is largely non-renally cleared, which softens
//      this, but its two main metabolites do accumulate in renal impairment
//      and this cohort's renal status is simply unknown -- it is not
//      "normal", it is unmeasured. Any PJI patient with meaningful renal
//      impairment is outside anything this model can be said to cover.
//
//   6. THE DONOR POPULATION'S PHYSIOLOGY INFLATES BOTH VOLUME AND CLEARANCE.
//      The paper states plainly that "the average values we found for volume
//      of distribution and clearance were higher" than published healthy-
//      volunteer values, and that steady-state AUC was significantly LOWER
//      than in healthy volunteers (P = 0.017). Septic shock with fluid
//      resuscitation, capillary leak and vasopressors expands volume and
//      raises clearance. An elective or stable PJI patient is NOT this
//      patient. Transported unchanged, this model will UNDER-predict
//      exposure in a stable PJI cohort -- so it is conservative for a
//      target-attainment question and anti-conservative for a toxicity one.
//
//   7. A SINGLE STATIC PARAMETER SET IS CONTRADICTED BY THE PAPER'S OWN
//      FIGURES. One parameter set could not describe both the first dose and
//      steady state in the same patient (Fig. 2); single-dose and
//      steady-state data had to be fitted separately. CL and Vss both moved
//      within every patient between visits, in opposing directions and by
//      differing amounts (Fig. 3), with shallower terminal slopes in 4/9.
//      Group means did not differ (Vss P=0.878, CL P=0.285, AUC P=0.169), so
//      this is purely INTRA-individual drift that no omega in this file
//      represents. There is no IOV term because none was estimated. A model
//      used across a multi-week PJI course is being asked for exactly the
//      stability the source demonstrates it does not have.
//
//   8. THE SAMPLING WINDOW WAS TOO SHORT TO IDENTIFY THE TERMINAL PHASE.
//      Sampling ran 0-8 h of a 12 h interval (20 points/matrix; every 20 min
//      to 3 h, every 30 min to 8 h). The 8-12 h segment, INCLUDING THE TRUE
//      TROUGH, is extrapolated with the terminal slope
//      (C12 = Cz*exp(-lambda_z*dt)), and fAUC0-24 was obtained simply by
//      DOUBLING fAUC0-12. So t1/2 = 5.0 h, Cmin = 1.83 mg/L and the whole
//      24 h exposure rest on an extrapolation, not on observation, and Fig. 1
//      thins to n = 6 per time point at late times. Trough predictions from
//      this file are the weakest thing it produces.
//
//   9. fu = 0.866 IS COHORT-SPECIFIC AND IS DOING REAL WORK HERE. It was
//      measured on only two plasma samples per patient per visit (22
//      determinations) with a 30 kDa membrane, with no fu covariate model and
//      no albumin measured anywhere. It ranged 0.730-0.959 between patients.
//      It is also materially HIGHER than the ~0.69 free fraction usually
//      quoted for linezolid -- plausibly the hypoalbuminaemia of sepsis. The
//      entire total-drug output of this file scales as 1/fu, so a PJI patient
//      with normal albumin could have total exposure ~25% different from what
//      this file returns for reasons that have nothing to do with clearance.
//
//  10. NO RESIDUAL ERROR MAGNITUDE IS PUBLISHED. WinNonlin 1/Cpred^2
//      weighting implies a constant-CV error, but no value is printed. The
//      house default 0.01 proportional / 0.0001 additive is used and is
//      recorded as ASSUMED. Assay context only: LLOQ 0.2 mg/L plasma and
//      0.8 mg/L ultrafiltrate, HPLC interday CV 6.1%, relative error 3.4%.
//
//  11. NO BONE WAS SAMPLED. THE ISF CAPTURES BELOW ARE NOT BONE. The tissues
//      are SUBCUTANEOUS ADIPOSE and SKELETAL MUSCLE interstitium of the LOWER
//      EXTREMITY, uninfected at the probe site, in patients with SYSTEMIC
//      sepsis. The only three occurrences of the string "bone" in the whole
//      paper are the author surname Bone (ref 3) and two reference titles.
//      FTSC = 0.896 and FTIM = 0.999 MUST NOT be promoted into a linezolid
//      bone penetration ratio -- the Stolle linezolid bone microdialysis
//      paper, already in the Reference folder, is the source for that.
//      Two further cautions on these two numbers even for soft tissue:
//      (a) they are AUC(t0 to LAST SAMPLE, ~8 h) ratios, not AUC0-12 or
//          AUC0-24 ratios;
//      (b) they are MEDIANS of a ~7-fold spread (s.c. 20.2-118%, i.m.
//          24.1-144%, n=10 each). One patient of twelve was below 25% in
//          BOTH matrices simultaneously. A scalar point factor is badly
//          misleading for this drug; values above 1.0 are the measurement
//          noise floor of a free-vs-free comparison, not accumulation.
//      Microdialysis relative recovery was only 53.1% (s.c., cv 31.0%) and
//      59.1% (i.m., cv 17.0%), and 4 of 20 probes had no calibration sample
//      of their own and were assigned the median recovery of the other probes
//      in the same tissue. Recovery error propagates linearly into these
//      ratios. They are captured here ONLY because the paper's own integrated
//      model (Fig. 4) applies FT as a fixed proportionality constant on a
//      disposition compartment -- this file reproduces that pattern, not any
//      estimated tissue kinetics.
//
//  12. THE PK/PD MESSAGE IS THAT STANDARD DOSING IS ALREADY BORDERLINE
//      BEFORE ANY BIOFILM MULTIPLIER. At a MIC of 4 mg/L (the breakpoint,
//      not a measured organism MIC -- no pathogen MICs were determined),
//      Table 3 gives mean fAUC/MIC 29.9 (UF plasma) / 29.1 (s.c.) / 33.0
//      (i.m.) against a Rayner target of >=51 that was reached in only two
//      instances across all patients and matrices, and mean fT>MIC of 59% /
//      54% / 90% with MINIMA of 16% / 0% / 24% -- one patient spent zero
//      time above MIC in subcutaneous interstitium. Even the softer animal
//      fT>MIC >=40% target was missed by 3/9 in plasma, 4/9 s.c. and 2/9
//      i.m. The authors propose considering 600 mg three times daily.
//      Regarding the library's unsourced 8-256x biofilm multiplier: peak
//      unbound here is ~16 mg/L, so even 8x on a 4 mg/L MIC (= 32 mg/L) is
//      never exceeded at any point in the interval. That is CONSISTENT with
//      "systemic therapy alone cannot clear an implant biofilm" but is NOT
//      evidence for the multiplier's magnitude, which this paper does not
//      address at all (grep: biofilm = 0 hits, implant = 0, prosthe = 0).
//
//  13. SMALL, SINGLE-CENTRE, 2006, AND ATTRITED. n = 12 enrolled, 10
//      evaluable at steady state, 9 for the graphical/PD analysis. Two
//      patients (16.7%) had serious adverse events and died; one patient's
//      steady-state sampling stopped at 2 h. Protocol deviations: one patient
//      received the visit-2 infusion over 1.25 h instead of 30 min, another
//      received only three doses at 21 h and 24 h intervals. The paper is
//      also internally inconsistent on its own n=9/n=10 denominators
//      (Results says 3/9 patients had plasma fT>MIC < 40%; the Discussion
//      says 7/10 had effective concentrations), and the Methods PROSE defines
//      FT inverted relative to the DISPLAYED EQUATION -- the equation
//      direction, FT = AUC(ISF)/AUC(plasma), is the one consistent with every
//      reported value and is the one used here.
//
// SMOKE EXPECTATION
//   Linezolid 600 mg IV over 30 min q12h, 80 kg adult, steady state
//   (WT is ignored by this model -- LIMITS 4):
//     TOTAL plasma  : Cmax ~18.9 mg/L, Cmin ~2.4 mg/L, AUC24 ~151 mg*h/L
//     UNBOUND plasma: Cmax ~16.4 mg/L, Cmin ~2.1 mg/L, AUC24 ~131 mg*h/L
//     terminal t1/2 ~5.0 h
//   THIS DOES NOT LAND IN THE HOUSE SMOKE BAND, AND THAT IS THE CORRECT
//   RESULT, NOT A DEFECT. The house band is AUC24 200-400 mg*h/L with a
//   trough of 2-8 mg/L. The trough is in band (~2.4 mg/L, at the bottom).
//   AUC24 of ~151 mg*h/L is roughly 25% BELOW the bottom of the band because
//   this is a septic-shock ICU cohort whose clearance was elevated and whose
//   steady-state AUC the paper itself reports as significantly lower than in
//   healthy volunteers (P = 0.017). No parameter has been adjusted to bring
//   it into range. A compile test should assert ~151 mg*h/L, and should
//   FAIL if it ever returns 200-400 from these numbers.
//   Internal cross-checks a test can also assert:
//     AUC24_unbound = 1200/9.18 = 130.7 (paper prints 2 x 65.3 = 130.6)
//     fAUC24/MIC at MIC 4 mg/L = 32.7 (paper's mean across patients 29.9)
// =========================================================================

[SET] end=72, delta=0.1


$PARAM @annotated
TVCLU  : 9.18  : Unbound-basis clearance, printed Table 1 (L/h)
TVVSSU : 62.9  : Unbound-basis Vss = V1+V2, printed Table 1 (L)
TVQU   : 62.1  : Unbound-basis intercompartmental clearance CLD2 (L/h)
FRV1   : 0.4027 : Fraction of Vss assigned to V1 -- DERIVED, NOT PUBLISHED
FU     : 0.866 : Unbound plasma fraction, measured (set to 1 for unbound sim)
FTSC   : 0.896 : s.c. adipose ISF / plasma AUC ratio -- SOFT TISSUE, NOT BONE
FTIM   : 0.999 : skeletal muscle ISF / plasma AUC ratio -- NOT BONE

ETA1 : 0 : IIV on clearance
ETA2 : 0 : IIV on steady-state volume of distribution
ETA3 : 0 : IIV on intercompartmental clearance


$PARAM @annotated @covariates
// DECLARED FOR HARNESS COMPATIBILITY AND DELIBERATELY UNUSED.
// Buerger 2006 fits no covariate relation of any kind -- no allometry, no
// centering constant, no renal term -- despite a 55-133 kg weight range and
// no reported renal function at all. See LIMITS 4 and 5. Adding a covariate
// term here would be invention, not transcription.
WT    : 81.0 : total body weight (kg) -- cohort median; NOT used by the model
AGE   : 62.0 : age (years) -- cohort median; NOT used by the model
SEX   : 0 : 0 = male, 1 = female; NOT used by the model
CREAT : 88.4 : serum creatinine (umol/L) -- NOT reported in the source at all


// IIV. THESE ARE NOT PUBLISHED OMEGAS -- see LIMITS 3. Table 1 prints
// "cv (%) geometric means" of the individual two-stage estimates (n = 10),
// i.e. COEFFICIENTS OF VARIATION IN PERCENT -- not variances, not SDs in
// natural units, and not log-scale SDs. The correct conversion for a CV% is
// therefore ln(1 + CV^2), applied below:
//   ETA1 CL   : cv 57.9% -> ln(1 + 0.579^2) = 0.289112
//   ETA2 Vss  : cv 19.2% -> ln(1 + 0.192^2) = 0.036201
//   ETA3 CLD2 : cv 87.9% -> ln(1 + 0.879^2) = 0.572471
// The arithmetic is exact; the INPUTS are the assumption. Set all three to 0
// for the deterministic typical-value benchmark this file is really for.
$OMEGA 0.289112 0.036201 0.572471


// ASSUMED -- no residual error magnitude is published (LIMITS 10).
// House default: proportional then additive.
$SIGMA
0.01
0.0001


$CMT @annotated
CENT   : Central compartment (mg) [ADM, OBS]
PERIPH : Peripheral compartment (mg)


$MAIN
// Printed parameters are on an UNBOUND basis (Methods: "the compartmental
// data analysis was based on unbound concentrations for all matrices").
// Scaling CL, Q, V1 and V2 all by FU leaves every CL/V rate constant
// unchanged and reinterprets the compartment amounts as TOTAL drug, so DV
// below is a TOTAL plasma concentration like the rest of this library.
// FU = 1 recovers the paper's native unbound simulation. See LIMITS 9.
double CLU  = TVCLU  * exp(ETA1 + ETA(1));
double VSSU = TVVSSU * exp(ETA2 + ETA(2));
double QU   = TVQU   * exp(ETA3 + ETA(3));

double CL = CLU * FU;                    // 9.18  * 0.866 = 7.950 L/h
double Q  = QU  * FU;                    // 62.1  * 0.866 = 53.78 L/h
double V  = VSSU * FRV1 * FU;            // central,    ~21.9 L total-basis
double V2 = VSSU * (1.0 - FRV1) * FU;    // peripheral, ~32.5 L total-basis


$ODE
dxdt_CENT   = - (CL/V)*CENT - (Q/V)*CENT + (Q/V2)*PERIPH;
dxdt_PERIPH =                 (Q/V)*CENT - (Q/V2)*PERIPH;


$TABLE
double CP = CENT / V;                    // TOTAL plasma concentration (mg/L)
double CU = CP * FU;                     // unbound plasma concentration
// Soft-tissue interstitial concentrations, applied as the paper's own fixed
// proportionality factors on unbound plasma (Fig. 4 pattern). SUBCUTANEOUS
// ADIPOSE AND SKELETAL MUSCLE -- THESE ARE NOT BONE. See LIMITS 11.
double CSC = CU * FTSC;
double CIM = CU * FTIM;

double DV = CP * (1 + EPS(1)) + EPS(2);
int i = 0;
while(DV<0 && i <100) {
simeps();
DV = CP * (1 + EPS(1)) + EPS(2);
++i;
}


$CAPTURE DV CP CU CSC CIM CL Q V V2
