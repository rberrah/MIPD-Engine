$PROB
// Rifampicin chronic-dose mean-profile reconstruction
// Article: Loos U et al. Klin Wochenschr. 1985;63:1205-1211.
// DOI: 10.1007/BF01733779
// Population: Healthy adults receiving oral and intravenous rifampicin during chronic administration.
// Implementation: deterministic PK module adapted for standalone TDM.
// The published chronic-dose clearance and bioavailability are retained. Volume and ka were not reported in the accessible source summary and remain explicit profile-shape assumptions; this is not a PopPK prior.

$PARAM @annotated
TVCL : 9.03 : Typical clearance (L/h)
TVV  : 50 : Typical volume (L)
TVKA : 1.5 : First-order absorption rate (1/h)
TVF  : 0.68 : Oral bioavailability
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
