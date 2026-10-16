$PROB
// Cobicistat 50 mg sensitivity module
// Article: Mathias AA et al. Clin Pharmacol Ther. 2010;87:322-329.
// DOI: 10.1038/clpt.2009.228
// Population: Healthy adults in a dose-ranging pharmacokinetic study.
// Implementation: deterministic PK module adapted for standalone TDM.
// Validated only as a 50 mg once-daily sensitivity module; extrapolation to marketed 150 mg is not supported. OMEGA/SIGMA are engineering priors.

$PARAM @annotated
TVCL : 6.05 : Typical clearance (L/h)
TVV  : 26.2 : Typical volume (L)
TVKA : 1.5 : First-order absorption rate (1/h)
TVF  : 0.1 : Oral bioavailability
ETA1 : 0 : MAP random effect on clearance

$OMEGA @annotated
IIV_CL : 0.09 : Engineering prior variance on clearance

$SIGMA @annotated
PROP : 0.04 : Engineering prior proportional residual variance
ADD  : 0.0  : Additive residual variance fixed to zero

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]

$MAIN
double CL = TVCL * exp(ETA1 + ETA(1));
double V = TVV;

$ODE
dxdt_GUT = -TVKA * GUT;
dxdt_CENT = TVKA * GUT * TVF - (CL / V) * CENT;

$TABLE
double IPRED = CENT / V;
double DV = IPRED * (1.0 + EPS(1));

$CAPTURE DV CL
