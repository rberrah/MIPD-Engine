$PROB
// Voriconazole nonlinear PK model (van den Born et al.)
// Article: van den Born DA et al. Int J Antimicrob Agents. 2023;61:106750.
// DOI: 10.1016/j.ijantimicag.2023.106750
// Population: 54 predominantly adult patients treated with voriconazole; mixed hospital setting.
// Implementation: deterministic PK module adapted for standalone TDM.
// Published IIV was retained for Vmax; the residual prior is an engineering adaptation for mapbayr.

$PARAM @annotated
TVKA   : 0.62 : Absorption rate (1/h)
TVVD   : 145 : Volume of distribution (L)
TVKM   : 5.7 : Michaelis-Menten constant (mg/L)
TVVMAX : 86.4 : Maximum elimination rate (mg/h)
TVF    : 0.83 : Oral bioavailability
CRPEK  : -0.0046 : Exponential CRP effect on Vmax
ETA1   : 0 : MAP random effect on Vmax

$PARAM @annotated @covariates
CRP : 50 : C-reactive protein (mg/L)

$OMEGA @annotated
IIV_VMAX : 0.683 : Approximate log-variance derived from published CV 99 percent

$SIGMA @annotated
PROP : 0.04 : Engineering prior proportional residual variance
ADD  : 0.0  : Additive residual variance fixed to zero

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]

$MAIN
double VMAX = TVVMAX * exp(CRPEK * CRP) * exp(ETA1 + ETA(1));

$ODE
double CP = CENT / TVVD;
dxdt_GUT = -TVKA * GUT;
dxdt_CENT = TVKA * GUT * TVF - VMAX * CP / (TVKM + CP);

$TABLE
double IPRED = CENT / TVVD;
double DV = IPRED * (1.0 + EPS(1));

$CAPTURE DV VMAX
