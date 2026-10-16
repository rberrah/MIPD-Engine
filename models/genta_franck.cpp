// ===========================================================================
// genta_franck.cpp -- Gentamicin, adult patients on intermittent haemodialysis
//
// SOURCE
//   Franck B. et al., "Population pharmacokinetics of gentamicin in
//   haemodialysis", Eur J Clin Pharmacol 2020;76:947-955.
//   PDF : Population pharmacokinetics of gentamicin in
//         haemodialysis Benedicte Franck.pdf
//   ESM1: ESM1_Franck Bene.docx  (Pmetrics model file and
//         the 14 non-parametric support points -- without it this model has no
//         reproducible between-subject variability at all)
//
// POPULATION  23 adults on chronic intermittent haemodialysis. Median age
//   79.2 y, no patient below 61.4 y. Median creatinine 631 umol/L. Study was
//   terminated early (50 planned). Internal evaluation only.
//
// STRUCTURE   Two compartments, parameterised in RATE CONSTANTS rather than
//   clearances -- the authors state the choice was arbitrary. Implemented in
//   rate constants here so nothing is a derived reconstruction.
//     Ke = Ke1 * Ke2^DIALYSIS      DIALYSIS is a 0/1 indicator
//     V  = V1 * (WT/73.8)^allo     allo FIXED to 1, so scaling is LINEAR
//   ESM1 confirms creatinine was supplied to the fit (CREATBEF, CREATAFT) and
//   does NOT appear in the #Sec block: it was tested and rejected. This model
//   therefore has NO renal-function mechanism whatsoever.
//
// PARAMETERS  Probability-weighted geometric means over the 14 ESM1 support
//   points, with omega^2 the probability-weighted variance of log(parameter):
//     Ke1  0.0243 /h    omega^2 0.5462   (CV 85%)
//     Ke2  31.44        omega^2 0.3882   (CV 69%)   dialysis multiplier
//     V1   9.77 L       omega^2 0.1366   (CV 38%)   at 73.8 kg
//     KCP  0.895 /h     omega^2 0.6721   (CV 98%)
//     KPC  0.531 /h     omega^2 1.0760   (CV 139%)
//   These differ slightly from the medians printed in Table 2 (Ke1 0.027,
//   V1 11.7) because the paper reports medians of the support distribution and
//   these are probability-weighted geometric means. The support points are the
//   primary source and are used here.
//
//   Implied disposition: off dialysis Ke 0.0243 /h, t1/2 28.5 h, CL 0.24 L/h.
//   On dialysis Ke 0.765 /h, t1/2 0.91 h, CL 7.5 L/h -- a 31-fold swing.
//
// RESIDUAL    Pmetrics polynomial SD = 0.25 + 0.1*C with process noise
//   lambda = 0.1, combined as sqrt(SD^2 + lambda^2). mrgsolve composes its
//   proportional and additive terms in quadrature rather than linearly, so the
//   two cannot be made identical. The $SIGMA below matches the source exactly
//   at C = 0 and runs about 12% low at C = 20 mg/L. Documented, not hidden.
//
// A NOTE ON THE ESM1 OUTPUT EQUATION. ESM1 gives Y(1) = (X(1)/V)*WEIGHT, which
// only balances dimensionally if doses were entered in mg/kg -- then X is
// mg/kg, and (mg/kg)/L * kg = mg/L. Doses here are in mg, so the equivalent
// output is simply CENT/V. Do not transcribe the *WEIGHT term.
//
// ---------------------------------------------------------------------------
// LIMITS -- this model is DISABLED in poppk_library.yaml. Reasons:
//
//   1. THE DIALYSIS EFFECT IS TIME-VARYING WITHIN A SESSION. Ke2 = 31 applies
//      only while the dialyser is running (median session 3.92 h, range
//      1.3-4.5 h). ESM1 declares DIALDUR and DIALDEB precisely to carry the
//      session duration and start time. The SAFEBONE harness supplies DIAL and
//      IHD as STATIC patient-level flags, so driving Ke2 from them would apply
//      a 31-fold clearance increase around the clock. That is the single most
//      likely way to get this model catastrophically wrong, and it is why the
//      model is gated rather than merely annotated.
//
//   2. NO RENAL COVARIATE. Applied to preserved renal function this model gives
//      t1/2 = 28.5 h and would under-predict clearance by roughly an order of
//      magnitude, producing simulated accumulation that does not occur.
//
//   3. IT CANNOT VALIDATE genta_debord AND VICE VERSA. genta_debord's whole
//      elimination mechanism is a non-steady-state Jelliffe creatinine
//      clearance; this model has no creatinine term at all. The two share no
//      covariate space, so holding one out tests nothing about the other. The
//      aminoglycoside external-validation gap is NOT closed by this model.
//
//   4. n = 23, elderly (no patient under 61.4 y), ESRD, internal evaluation
//      only, study terminated early.
//
//   5. The 73.8 kg centering constant is ambiguous in the source: Methods and
//      Limitations say LEAN body weight, Table 1 labels the row "Weight (kg)".
//      Total body weight is used here because that is what the harness supplies
//      as WT; if lean weight was intended, V is overestimated in obese
//      patients. Flagged rather than silently resolved.
// ===========================================================================

[SET] end=168, delta=0.1


$PARAM @annotated
TVKE1 : 0.0243 : Elimination rate constant off dialysis (1/h)
TVKE2 : 31.44  : Multiplicative dialysis effect on Ke (dimensionless)
TVV1  : 9.77   : Central volume at the reference weight (L)
TVKCP : 0.895  : Central to peripheral rate constant (1/h)
TVKPC : 0.531  : Peripheral to central rate constant (1/h)
WTREF : 73.8   : Reference weight (kg) -- see LIMITS 5

ETA1 : 0 : IIV on Ke1
ETA2 : 0 : IIV on Ke2
ETA3 : 0 : IIV on V1
ETA4 : 0 : IIV on KCP
ETA5 : 0 : IIV on KPC


$PARAM @annotated @covariates
WT   : 73.8 : body weight (kg)
IHD  : 0    : haemodialysis session in progress (0/1) -- see LIMITS 1
DIAL : 0    : chronic dialysis patient (0/1), not used by the equations


// Probability-weighted variance of log(parameter) over the 14 ESM1 support
// points: Ke1, Ke2, V1, KCP, KPC.
$OMEGA 0.5462 0.3882 0.1366 0.6721 1.0760


$SIGMA
0.01     // proportional -- Pmetrics C1 = 0.1, so variance 0.1^2
0.0724   // additive     -- sqrt(0.25^2 + 0.1^2) = 0.269, so variance 0.0724


$CMT @annotated
CENT   : Central compartment (mg) [ADM, OBS]
PERIPH : Peripheral compartment (mg)


$MAIN
double KE1 = TVKE1 * exp(ETA1 + ETA(1));
double KE2 = TVKE2 * exp(ETA2 + ETA(2));
double V   = TVV1  * (WT / WTREF) * exp(ETA3 + ETA(3));   // allo fixed to 1
double KCP = TVKCP * exp(ETA4 + ETA(4));
double KPC = TVKPC * exp(ETA5 + ETA(5));

// IHD must be a WITHIN-SESSION indicator, not a patient-level flag. The harness
// does not yet express that, which is why this model is disabled.
double KE = KE1 * pow(KE2, IHD);

double CL = KE * V;   // reported for review only; the ODEs use rate constants


$ODE
dxdt_CENT   = - (KE + KCP) * CENT + KPC * PERIPH;
dxdt_PERIPH =         KCP  * CENT - KPC * PERIPH;


$TABLE
double DV = (CENT/V) * (1 + EPS(1)) + EPS(2);
int i = 0;
while(DV<0 && i <100) {
simeps();
DV = (CENT/V) * (1 + EPS(1)) + EPS(2);
++i;
}


$CAPTURE DV CL V KE KE1 KE2
