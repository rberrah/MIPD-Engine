$PROB
// Long-term oral amiodarone population PK model
// Article: Pollak PT et al. Clin Pharmacol Ther. 2000;67:642-652.
// DOI: 10.1067/mcp.2000.107047
// Population: Adults receiving long-term oral amiodarone therapy for cardiac arrhythmia.
// Implementation: deterministic PK module adapted for standalone TDM.
// Published apparent disposition parameters and IIV are retained. The accessible source does not identify an absorption model, so first-order ka is an explicit simulation assumption; this is not a MAP prior.

$PARAM @annotated
TVCL : 9.541667 : Published CL/F, 229 L/day converted to L/h
TVV1 : 882 : Published central V/F (L)
TVQ  : 24.958333 : Published intercompartmental CL/F, 599 L/day converted to L/h
TVV2 : 12700 : Published peripheral V/F (L)
TVKA : 0.6 : Assumed first-order absorption rate (1/h)
ETA1 : 0 : Random effect on CL/F
ETA2 : 0 : Random effect on Q/F
ETA3 : 0 : Random effect on V2/F

$OMEGA @annotated
IIV_CL : 0.091755 : log(1 + 0.31^2)
IIV_Q  : 0.272423 : log(1 + 0.56^2)
IIV_V2 : 0.289783 : log(1 + 0.58^2)

$SIGMA @annotated
PROP : 0.04 : Engineering residual variance; not reported in the accessible source
ADD  : 0.0 : Additive residual variance fixed to zero

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]
PERI : Peripheral compartment

$MAIN
double CL = TVCL * exp(ETA1 + ETA(1));
double Q = TVQ * exp(ETA2 + ETA(2));
double V2 = TVV2 * exp(ETA3 + ETA(3));

$ODE
dxdt_GUT = -TVKA * GUT;
dxdt_CENT = TVKA * GUT - ((CL + Q) / TVV1) * CENT + (Q / V2) * PERI;
dxdt_PERI = (Q / TVV1) * CENT - (Q / V2) * PERI;

$TABLE
double IPRED = CENT / TVV1;
double DV = IPRED * (1.0 + EPS(1)) + EPS(2);

$CAPTURE DV CL Q V2
