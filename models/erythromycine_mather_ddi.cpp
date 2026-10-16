$PROB
// Erythromycin oral mean-profile reconstruction
// Article: Mather LE et al. Br J Clin Pharmacol. 1981;12:131-140.
// DOI: 10.1111/j.1365-2125.1981.tb01191.x
// Population: Healthy adults in oral absorption and bioavailability studies.
// Implementation: deterministic PK module adapted for standalone TDM.
// The source compares formulation-, meal- and occasion-dependent bioavailability after prior individual IV characterization; it does not publish a reusable PopPK model. F=0.40 represents erythromycin base after meals, while CL, V and ka remain explicit simulation assumptions.

$PARAM @annotated
TVCL : 20 : Typical clearance (L/h)
TVV  : 200 : Typical volume (L)
TVKA : 1 : First-order absorption rate (1/h)
TVF  : 0.4 : Oral bioavailability
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
