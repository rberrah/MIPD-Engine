$PROB
// Diltiazem oral mean-profile reconstruction
// Article: Hermann P et al. Eur J Clin Pharmacol. 1983;24:349-352.
// DOI: 10.1007/BF00610053
// Population: Healthy adults after intravenous and oral diltiazem.
// Implementation: deterministic PK module adapted for standalone TDM.
// The source reports biphasic IV disposition and summary NCA values, not a reusable PopPK model. The 70 kg mean CL, apparent V and oral F are retained in a one-compartment terminal-profile approximation; ka is assumed.

$PARAM @annotated
TVCL : 89.6 : Typical clearance (L/h)
TVV  : 371 : Typical volume (L)
TVKA : 1 : First-order absorption rate (1/h)
TVF  : 0.42 : Oral bioavailability
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
