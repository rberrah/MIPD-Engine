$PROB
// Ritonavir population PK model
// Article: Kappelhoff BS et al. Br J Clin Pharmacol. 2005;59:174-182.
// DOI: 10.1111/j.1365-2125.2004.02241.x
// Population: Adults living with HIV-1; ritonavir used as booster or antiviral.
// Implementation: deterministic PK module adapted for standalone TDM.
// Published structure, typical values, lag, lopinavir covariate and IIV are retained. IOV on F, the V-ka covariance and the two-population residual mixture are omitted; the dominant residual component is used, so this remains simulation-only.

$PARAM @annotated
TVCL : 10.5 : Apparent clearance without lopinavir (L/h)
TVV  : 96.6 : Apparent volume (L)
TVKA : 0.871 : First-order absorption rate (1/h)
TVLAG : 0.778 : Absorption lag time (h)
LPV_CL : 2.72 : Clearance multiplier with concomitant lopinavir
ETA1 : 0 : Random effect on CL/F
ETA2 : 0 : Random effect on V/F
ETA3 : 0 : Random effect on ka

$PARAM @annotated @covariates
LPV : 0 : Concomitant lopinavir (0 = no, 1 = yes)

$OMEGA @annotated
IIV_CL : 0.136671 : log(1 + 0.383^2)
IIV_V  : 0.494696 : log(1 + 0.800^2)
IIV_KA : 1.350319 : log(1 + 1.69^2)

$SIGMA @annotated
PROP : 0.023716 : Published proportional SD 0.154 squared
ADD  : 0.0036 : Dominant published additive SD 0.060 mg/L squared

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]

$MAIN
double CL = TVCL * (LPV > 0.5 ? LPV_CL : 1.0) * exp(ETA1 + ETA(1));
double V = TVV * exp(ETA2 + ETA(2));
double KA = TVKA * exp(ETA3 + ETA(3));
ALAG_GUT = TVLAG;

$ODE
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - (CL / V) * CENT;

$TABLE
double IPRED = CENT / V;
double DV = IPRED * (1.0 + EPS(1)) + EPS(2);

$CAPTURE DV CL V KA
