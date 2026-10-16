$PROB
# Rifampicin PopPK -- Jeremiah et al. (2014)
#
# Primary source: Jeremiah K et al. Antimicrob Agents Chemother.
# 2014;58:3468-3474. DOI: 10.1128/AAC.02307-13
# Tanzania, n=100 adults with pulmonary TB, with or without HIV, 450-600 mg oral.
# The structural parameters were checked against Table 2. This implementation
# does not yet include the published HIV, nutritional-supplementation,
# continuation-phase, or interoccasion effects.
#
# Structure: one compartment, first-order elimination, TRANSIT-compartment
# absorption, and -- the reason this model is worth having -- an explicit
# TIME-DEPENDENT AUTO-INDUCTION of clearance:
#
#   CL(t) = CL7 + (CLss - CL7) * (1 - exp(-ln(2)/t_half_ind * t))
#
# with CL7 = 13.9*(FFM/43)^0.75 (day-7 clearance), CLss = 16.5*(FFM/43)^0.75
# (steady-state clearance) and an induction half-life of 6 DAYS. Rifampicin induces
# its own metabolism, so a model that treats clearance as constant over a 6-12 week
# bone-infection course is wrong in a direction that matters: it overestimates
# exposure late in treatment.
#
# TIME UNITS: `t` in the induction term is in DAYS while the ODE runs in HOURS, so
# the conversion is made explicit below rather than left implicit.
#
# CAVEAT ON TRANSPORTABILITY: this is a pulmonary-tuberculosis population, not a
# bone or prosthetic-joint infection population, and rifampicin PK in PJI patients
# on 600-900 mg twice daily is not established. Using it for PJI is an
# extrapolation and must be reported as one.
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
# IIV is entered as VARIANCE on the log scale, converted from the published %CV by
# omega^2 = ln(1 + CV^2): CL 24.0% -> 0.0560. The paper's %CV footnote convention is
# not fully specified; this is the exact lognormal relation and is stated here so the
# choice is auditable rather than hidden.

$PARAM @annotated
TVCL7   : 13.9 : Day-7 clearance at FFM 43 kg (L/h)
TVCLSS  : 16.5 : Steady-state (fully induced) clearance at FFM 43 kg (L/h)
TVVD    : 55.8 : Volume of distribution at FFM 43 kg (L)
TVKA    : 1.77 : Absorption rate constant (1/h)
TVMTT   : 1.50 : Mean transit time (h)
NN      : 27.6 : Number of transit compartments
THALFI  : 6.0  : Auto-induction half-life (DAYS) - FIXED in the source model
FFMREF  : 43.0 : Reference fat-free mass (kg)
ALLOCL  : 0.75 : Allometric exponent on clearance

ETA1 : 0 : IIV on clearance (mapbayr)

$PARAM @annotated @covariates
WT  : 55.0  : Body weight (kg)
HT  : 165.0 : Height (cm)
SEX : 0.0   : 0 = male, 1 = female

$OMEGA @annotated @diagonal
ETA_CL : 0.0560 : IIV on CL, from 24.0 %CV

$SIGMA @annotated @diagonal
PROP : 0.01876 : Proportional residual error (13.7 %)
ADD  : 0.00174 : Additive residual error variance ((0.0417 mg/L)^2)

$CMT @annotated
TR1  : Transit compartment (mg) [ADM]
DEPOT: Absorption compartment (mg)
CENT : Central compartment (mg) [OBS]

$MAIN
double safe_WT = WT;
if(safe_WT < 30.0)  safe_WT = 30.0;
if(safe_WT > 250.0) safe_WT = 250.0;

double safe_HT = HT;
if(safe_HT < 120.0) safe_HT = 120.0;
if(safe_HT > 220.0) safe_HT = 220.0;

// Fat-free mass, Janmahasatian. The source model is parameterised on FFM, so
// deriving it here keeps the covariate interface the same as the other models in
// this library (WT / HT / SEX) instead of demanding FFM from the caller.
double BMI = safe_WT / pow(safe_HT / 100.0, 2.0);
double FFM = (SEX == 1.0)
    ? (9.27e3 * safe_WT) / (8.78e3 + 244.0 * BMI)
    : (9.27e3 * safe_WT) / (6.68e3 + 216.0 * BMI);
if(FFM < 15.0) FFM = 15.0;

double SIZE = pow(FFM / FFMREF, ALLOCL);

// Auto-induction. TIME is in hours; the induction half-life is in days.
double t_days = TIME / 24.0;
double CL7  = TVCL7  * SIZE;
double CLSS = TVCLSS * SIZE;
double CL   = CL7 + (CLSS - CL7) * (1.0 - exp(-log(2.0) / THALFI * t_days));
CL = CL * exp(ETA1 + ETA(1));

double V  = TVVD * (FFM / FFMREF);
double Ka = TVKA;
double KTR = (NN + 1.0) / TVMTT;

if(CL < 0.001) CL = 0.001;
if(V  < 1.0)   V  = 1.0;

$ODE
dxdt_TR1   = -KTR * TR1;
dxdt_DEPOT =  KTR * TR1 - Ka * DEPOT;
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
CL    : Individual clearance at this time, INDUCTION-DEPENDENT (L/h)
V     : Individual volume of distribution (L)
Ka    : Absorption rate constant (1/h)
FFM   : Fat-free mass (kg)
