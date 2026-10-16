$PROB
// Posaconazole integrated population PK model
// Article: Dolton MJ et al. Antimicrob Agents Chemother. 2014;58:6879-6885.
// DOI: 10.1128/AAC.03777-14
// Population: Healthy adults and adult patients receiving posaconazole for prophylaxis or treatment.
// Implementation: deterministic PK module adapted for standalone TDM.
// Published structure, fixed effects, lag and medication/clinical covariates are retained. Patient-study IIV and residual variance are used as the default; the study-specific error branches are not represented, so this remains simulation-only.

$PARAM @annotated
TVCL : 30.2 : Typical clearance (L/h)
TVV  : 1100 : Typical volume (L)
TVKA : 1.26 : First-order absorption rate (1/h)
TVF  : 1 : Oral bioavailability
TVLAG : 1.79 : Absorption lag time (h)
PHE_CL : 7.21 : CL/F multiplier with phenytoin
RIF_CL : 7.21 : CL/F multiplier with rifampicin
FOS_CL : 1.342 : CL/F multiplier with fosamprenavir
ETA1 : 0 : MAP random effect on clearance

$PARAM @annotated @covariates
PHE : 0 : Phenytoin coadministration
RIF : 0 : Rifampicin coadministration
FOS : 0 : Fosamprenavir coadministration
PPI : 0 : Proton-pump inhibitor
MET : 0 : Metoclopramide
NUT : 0 : Nutritional supplement
MUC : 0 : Mucositis
DIARRHEA : 0 : Diarrhea

$OMEGA @annotated
IIV_CL : 0.285156 : Patient-study IIV on CL/F, 53.4 percent squared

$SIGMA @annotated
PROP : 0.289444 : Patient-study proportional residual SD 53.8 percent squared
ADD  : 0 : Additive residual variance fixed to zero

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]

$MAIN
double CL = TVCL * pow(PHE_CL, PHE) * pow(RIF_CL, RIF) * pow(FOS_CL, FOS) * exp(ETA1 + ETA(1));
double V = TVV;
double KA = TVKA;
double F = TVF * pow(0.549, PPI) * pow(0.655, MET) * pow(2.29, NUT) * pow(0.423, MUC) * pow(0.549, DIARRHEA);
ALAG_GUT = TVLAG;

$ODE
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT * F - (CL / V) * CENT;

$TABLE
double IPRED = CENT / V;
double DV = IPRED * (1.0 + EPS(1)) + EPS(2);

$CAPTURE DV CL
