$PROB
// Cobicistat 50 mg sensitivity module
// Article: Mathias AA et al. Clin Pharmacol Ther. 2010;87:322-329.
// DOI: 10.1038/clpt.2009.228
// Population: Healthy adults in a dose-ranging pharmacokinetic study.
// Implementation: deterministic PK module adapted for standalone TDM.
// The source is a nonlinear dose-ranging NCA study, not a PopPK model. The published 50 mg steady-state apparent clearance and terminal volume are retained for this dose-specific sensitivity profile; ka and random effects are assumptions and extrapolation to 150 mg is unsupported.

$PARAM @annotated
TVCL : 154 : Typical clearance (L/h)
TVV  : 316 : Typical volume (L)
TVKA : 0.8 : First-order absorption rate (1/h)
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
