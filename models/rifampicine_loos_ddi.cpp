$PROB
// Rifampicin chronic-dose PK module
// Article: Loos U et al. Klin Wochenschr. 1985;63:1205-1211.
// DOI: 10.1007/BF01733779
// Population: Healthy adults receiving oral and intravenous rifampicin during chronic administration.
// Implementation: deterministic PK module adapted for standalone TDM.
// OMEGA and SIGMA below are engineering priors added for MAP compatibility; they were not estimated in the cited source.

$PARAM @annotated
TVCL : 25 : Typical clearance (L/h)
TVV  : 50 : Typical volume (L)
TVKA : 1.5 : First-order absorption rate (1/h)
TVF  : 0.7 : Oral bioavailability
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
