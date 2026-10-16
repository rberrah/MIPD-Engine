$PROB
// Long-term oral amiodarone PK module
// Article: Pollak PT et al. Clin Pharmacol Ther. 2000;67:642-652.
// DOI: 10.1067/mcp.2000.107047
// Population: Adults receiving long-term oral amiodarone therapy for cardiac arrhythmia.
// Implementation: deterministic PK module adapted for standalone TDM.
// OMEGA and SIGMA below are engineering priors added for MAP compatibility; they were not estimated in the cited source.

$PARAM @annotated
TVCL : 9 : Typical clearance (L/h)
TVV1 : 66 : Typical central volume (L)
TVQ  : 40 : Typical intercompartmental clearance (L/h)
TVV2 : 3000 : Typical peripheral volume (L)
TVKA : 0.6 : First-order absorption rate (1/h)
TVF  : 0.5 : Oral bioavailability
ETA1 : 0 : MAP random effect on clearance

$OMEGA @annotated
IIV_CL : 0.09 : Engineering prior variance on clearance

$SIGMA @annotated
PROP : 0.04 : Engineering prior proportional residual variance
ADD  : 0.0  : Additive residual variance fixed to zero

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]
PERI : Peripheral compartment

$MAIN
double CL = TVCL * exp(ETA1 + ETA(1));
double V1 = TVV1;
double Q = TVQ;
double V2 = TVV2;

$ODE
dxdt_GUT = -TVKA * GUT;
dxdt_CENT = TVKA * GUT * TVF - ((CL + Q) / V1) * CENT + (Q / V2) * PERI;
dxdt_PERI = (Q / V1) * CENT - (Q / V2) * PERI;

$TABLE
double IPRED = CENT / V1;
double DV = IPRED * (1.0 + EPS(1));

$CAPTURE DV CL
