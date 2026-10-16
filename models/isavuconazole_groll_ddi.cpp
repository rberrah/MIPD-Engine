$PROB
// Isavuconazole DDI-study mean-profile reconstruction
// Article: Groll AH et al. Clin Pharmacol Drug Dev. 2017;6:76-85.
// DOI: 10.1002/cpdd.284
// Population: Healthy adults in controlled drug-drug interaction studies.
// Implementation: deterministic PK module adapted for standalone TDM.
// The cited study reports noncompartmental DDI profiles, not a population model. CL, V and ka below are rounded mean-profile reconstruction values for simulation only; no claim of article-model reproduction is made.

$PARAM @annotated
TVCL : 2 : Typical clearance (L/h)
TVV  : 300 : Typical volume (L)
TVKA : 1 : First-order absorption rate (1/h)
TVF  : 1 : Oral bioavailability
ETA1 : 0 : MAP random effect on clearance

$OMEGA @annotated
IIV_CL : 0.09 : Engineering prior variance on clearance

$SIGMA @annotated
PROP : 0.04 : Engineering prior proportional residual variance
ADD  : 0 : Additive residual variance fixed to zero

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]

$MAIN
double CL = TVCL * exp(ETA1 + ETA(1));
double V = TVV;
double KA = TVKA;
double F = TVF;

$ODE
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT * F - (CL / V) * CENT;

$TABLE
double IPRED = CENT / V;
double DV = IPRED * (1.0 + EPS(1)) + EPS(2);

$CAPTURE DV CL
