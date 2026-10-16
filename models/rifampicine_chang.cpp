$PROB
# Rifampicin PopPK -- Chang et al. (2015)
#
# Primary source: Chang MJ et al. Tuberculosis (Edinb). 2015;95:54-59.
# DOI: 10.1016/j.tube.2014.10.013
# South Korea, n=54 adults with pulmonary TB, including 21 with diabetes
# mellitus, 450-600 mg administered orally after an overnight fast.
#
# Structure: one compartment, first-order absorption and elimination. The primary
# article explicitly reports ADVAN2/TRANS2; no transit compartment is used.
#
#   CL = 6.10 + (BMI/20.3) * 6.22
#   Vd = 48.0 + DM * 16.2
#   Ka = 1.31 + DM * 1.56
#
# IMPORTANT SOURCE AMBIGUITY: the printed CL equation contains
# `(BMI/20.3)6.22` without an operator, while Table 2 reports CL 6.10 L/h and a
# bootstrap median of 6.22 L/h. The multiplication above is the pre-existing
# repository interpretation; it is retained and documented rather than silently
# replacing the malformed published expression with an invented exponent.
#
# CAVEAT ON TRANSPORTABILITY: pulmonary-tuberculosis population with diabetes, not a
# bone or prosthetic-joint infection population. Using it for PJI is an
# extrapolation.
#
#
# POPULATION CAVEAT -- READ WITH rifampicine_marsot.cpp
# This model comes from a PULMONARY TUBERCULOSIS cohort. The library also contains
# `rifampicine_marsot.cpp`, which is Marsot et al., Br J Clin Pharmacol 2017;83:
# 1039-1047, "Population pharmacokinetics of rifampicin in adult patients with
# OSTEOARTICULAR INFECTIONS: interaction with fusidic acid" -- i.e. the SAFEBONE
# population itself. Marsot is therefore the reference model for bone infection and
# this one is the extrapolation, not the other way round.
#
# Note also that the Ju et al. review used to build this file lists Marsot (2017)
# under "Pulmonary tuberculosis" in its Table 1. That is a mis-annotation: the paper
# is osteoarticular. Provenance metadata taken from that review should be treated
# with corresponding caution.
#
# IIV from Table 2:
#   CL 32.8% exponential -> ln(1 + 0.328^2) = 0.1022
#   Vd 53.7% exponential -> ln(1 + 0.537^2) = 0.2535
#   Ka 49.9% proportional normal -> variance 0.499^2 = 0.2490

$PARAM @annotated
CL_INT  : 6.10  : Clearance intercept (L/h)
CL_BMI  : 6.22  : Clearance slope on BMI/20.3 (L/h)
BMIREF  : 20.3  : Reference BMI (kg/m2)
VD_INT  : 48.0  : Volume intercept (L)
VD_DM   : 16.2  : Additional volume if diabetic (L)
KA_INT  : 1.31  : Absorption rate intercept (1/h)
KA_DM   : 1.56  : Additional absorption rate if diabetic (1/h)

ETA1 : 0 : IIV on clearance (mapbayr)
ETA2 : 0 : IIV on volume (mapbayr)
ETA3 : 0 : IIV on absorption rate (mapbayr)

$PARAM @annotated @covariates
WT  : 54.0  : Body weight (kg)
HT  : 165.0 : Height (cm)
DM  : 0.0   : Diabetes mellitus (0 = no, 1 = yes)

$OMEGA @annotated @diagonal
ETA_CL : 0.1022 : IIV on CL, from 32.8 %CV
ETA_VD : 0.2535 : IIV on Vd, from 53.7 %CV
ETA_KA : 0.2490 : Proportional IIV on Ka, from 49.9 %CV

$SIGMA @annotated @diagonal
PROP : 0.0144 : Proportional residual error ((0.12)^2)
ADD  : 2.0164 : Additive residual error variance ((1.42 mg/L)^2)

$CMT @annotated
DEPOT : Absorption compartment (mg) [ADM]
CENT  : Central compartment (mg) [OBS]

$MAIN
double safe_WT = WT;
if(safe_WT < 30.0)  safe_WT = 30.0;
if(safe_WT > 250.0) safe_WT = 250.0;

double safe_HT = HT;
if(safe_HT < 120.0) safe_HT = 120.0;
if(safe_HT > 220.0) safe_HT = 220.0;

double BMI = safe_WT / pow(safe_HT / 100.0, 2.0);
if(BMI < 12.0) BMI = 12.0;
if(BMI > 60.0) BMI = 60.0;

// LINEAR covariate model, unlike the allometric-power form used elsewhere in this
// library. Kept exactly as published, including the intercept.
double CL = (CL_INT + (BMI / BMIREF) * CL_BMI) * exp(ETA1 + ETA(1));
double V  = (VD_INT + DM * VD_DM)              * exp(ETA2 + ETA(2));
double Ka = (KA_INT + DM * KA_DM)              * (1.0 + ETA3 + ETA(3));

if(CL < 0.001) CL = 0.001;
if(V  < 1.0)   V  = 1.0;
if(Ka < 0.001) Ka = 0.001;

$ODE
dxdt_DEPOT = -Ka * DEPOT;
dxdt_CENT  =  Ka * DEPOT - (CL / V) * CENT;

$TABLE
double IPRED = CENT / V;
double DV = IPRED * (1.0 + EPS(1)) + EPS(2);

int i = 0;
while(DV < 0.0 && i < 1000) {
  simeps();
  DV = IPRED * (1.0 + EPS(1)) + EPS(2);
  ++i;
}

$CAPTURE @annotated
DV    : Simulated concentration with residual error (mg/L)
IPRED : Individual predicted concentration (mg/L)
CL    : Individual clearance (L/h)
V     : Individual volume of distribution (L)
Ka    : Individual absorption rate constant (1/h)
BMI   : Body mass index (kg/m2)
