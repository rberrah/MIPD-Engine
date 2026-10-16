$PROB
// Fluconazole review-based mean-profile reconstruction
// Article: Debruyne D, Ryckelynck JP. Clin Pharmacokinet. 1993;24:10-27.
// DOI: 10.2165/00003088-199324010-00002
// Population: Adults from clinical pharmacokinetic studies summarized in the source review.
// Implementation: deterministic PK module adapted for standalone TDM.
// The source is a clinical PK review rather than a PopPK model. Mean clearance (1.17 L/h), bioavailability (0.90) and the volume implied by the reported 31.6 h half-life are used; ka and random effects are simulation assumptions.

$PARAM @annotated
TVCL : 1.17 : Typical clearance (L/h)
TVV  : 53.36 : Typical volume (L)
TVKA : 1.5 : First-order absorption rate (1/h)
TVF  : 0.9 : Oral bioavailability
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
