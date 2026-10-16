$PROB
# Modèle : Amoxicilline PopPK (Mellon et al. 2020)
# Population : adultes obèses recevant co-amoxiclav.
# Structure : 2 compartiments, élimination linéaire.
# Administration : IV dans CENT, orale dans DEPOT avec bioavailability F.
# Absorption orale : chaîne de transit DEPOT -> TRANSIT -> GUT -> CENT.
# Covariable : effet négatif du poids sur KA.
#
# Utilisation :
# - IV : administrer dans CENT, cmt = 1, par exemple perfusion 1000 mg sur 30 min.
# - Oral : administrer dans DEPOT, cmt = 3.
#
# Adaptation note:
# - The main article does not provide the executable transit-model equations.
#   The fixed two-stage transit chain below is an approximation.

$PARAM @annotated
TVF        : 0.797 : Bioavailability orale typique, fraction
TVMTT      : 1.00  : Mean transit time typique (h)
TVKTR      : 1.60  : Constante de transit typique (1/h)
TVKA       : 1.70  : Constante d'absorption orale typique (1/h)
TVV1       : 9.00  : Volume central typique (L)
TVCL       : 14.6  : Clairance typique (L/h)
TVV2       : 6.40  : Volume périphérique typique (L)
TVQ        : 4.20  : Clairance intercompartimentale typique (L/h)
BETA_WT_KA : -3.10 : Effet du poids sur KA
REF_WT     : 109.3 : Poids de référence, médiane de l'étude (kg)

// Mapbayr ETAs
ETA1 : 0.0 : ETA posterior sur F, échelle logit
ETA2 : 0.0 : ETA posterior sur MTT
ETA3 : 0.0 : ETA posterior sur KTR
ETA4 : 0.0 : ETA posterior sur KA
ETA5 : 0.0 : ETA posterior sur V1
ETA6 : 0.0 : ETA posterior sur CL
ETA7 : 0.0 : ETA posterior sur V2
ETA8 : 0.0 : ETA posterior sur Q

$PARAM @covariates @annotated
WT    : 109.3 : Poids du patient (kg)
AGE   : 51.7  : Âge du patient (années), conservé pour compatibilité
CREAT : 70.0  : Créatinine sérique, conservée pour compatibilité
SEX   : 1.0   : Sexe (0 = Homme, 1 = Femme), conservé pour compatibilité
BSA   : 1.73  : Surface corporelle, conservée pour compatibilité

$OMEGA @block
1.0000
0.0000 0.0625
0.0000 0.0000 0.2401
0.0000 0.0000 0.0000 0.0625
0.0000 0.0000 0.0000 0.0000 0.2500
0.0000 0.0000 0.0000 0.0000 0.0000 0.0729
0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0900
0.0000 0.0000 0.0000 0.0000 0.0000 0.14337 0.0000 0.3481

$SIGMA @annotated @diagonal
PROP : 0.01716 : Erreur résiduelle proportionnelle, variance correspondant à 13.1%
ADD  : 0.14669 : Erreur résiduelle additive, variance correspondant à 0.383 mg/L

$CMT @annotated
CENT    : Compartiment central (mg) [ADM, OBS]
PERIPH  : Compartiment périphérique (mg)
DEPOT   : Dépôt oral, compartiment d'administration orale (mg) [ADM]
TRANSIT : Compartiment de transit oral (mg)
GUT     : Compartiment d'absorption orale (mg)

$MAIN
// -----------------------------------------------------------------------------
// 1. Sécurisation des covariables
// -----------------------------------------------------------------------------
double safe_WT = WT;

if(safe_WT < 30.0)  safe_WT = 30.0;
if(safe_WT > 250.0) safe_WT = 250.0;

// -----------------------------------------------------------------------------
// 2. Bioavailability orale avec distribution logit-normale
// -----------------------------------------------------------------------------
double logit_TVF = log(TVF / (1.0 - TVF));
double F_ORAL = 1.0 / (1.0 + exp(-(logit_TVF + ETA(1) + ETA1)));

if(F_ORAL < 0.001) F_ORAL = 0.001;
if(F_ORAL > 0.999) F_ORAL = 0.999;

// mrgsolve : bioavailability appliquée uniquement au compartiment oral DEPOT
F_DEPOT = F_ORAL;   // mrgsolve pre-declares F_<cmt>; redeclaring it collides

// -----------------------------------------------------------------------------
// 3. Paramètres individuels
// -----------------------------------------------------------------------------
double MTT = TVMTT * exp(ETA(2) + ETA2);
double KTR = TVKTR * exp(ETA(3) + ETA3);

double KA = TVKA * pow(safe_WT / REF_WT, BETA_WT_KA) * exp(ETA(4) + ETA4);

double V1 = TVV1 * exp(ETA(5) + ETA5);
double CL = TVCL * exp(ETA(6) + ETA6);
double V2 = TVV2 * exp(ETA(7) + ETA7);
double Q  = TVQ  * exp(ETA(8) + ETA8);

// Sécurisation numérique
if(MTT < 0.01)  MTT = 0.01;
if(KTR < 0.001) KTR = 0.001;
if(KA  < 0.001) KA  = 0.001;
if(V1  < 0.001) V1  = 0.001;
if(V2  < 0.001) V2  = 0.001;
if(CL  < 0.001) CL  = 0.001;
if(Q   < 0.001) Q   = 0.001;

// Constantes de transfert
double K10 = CL / V1;
double K12 = Q  / V1;
double K21 = Q  / V2;

$ODE
// Absorption orale avec transit
dxdt_DEPOT   = -KTR * DEPOT;
dxdt_TRANSIT =  KTR * DEPOT - KTR * TRANSIT;
dxdt_GUT     =  KTR * TRANSIT - KA * GUT;

// Distribution et élimination
dxdt_CENT   = KA * GUT + K21 * PERIPH - (K10 + K12) * CENT;
dxdt_PERIPH = K12 * CENT - K21 * PERIPH;

$TABLE
double CP = CENT / V1;
double DV = CP * (1.0 + EPS(1)) + EPS(2);

if(DV < 0.0) DV = 0.0;

$CAPTURE @annotated
DV     : Variable dépendante simulée, concentration totale (mg/L)
CP     : Concentration totale prédite sans erreur résiduelle (mg/L)
CL     : Clairance individuelle (L/h)
V1     : Volume central individuel (L)
Q      : Clairance intercompartimentale individuelle (L/h)
V2     : Volume périphérique individuel (L)
K10    : Constante d'élimination individuelle (1/h)
K12    : Constante de transfert central vers périphérique (1/h)
K21    : Constante de transfert périphérique vers central (1/h)
F_ORAL : Bioavailability orale individuelle
MTT    : Mean transit time individuel (h)
KTR    : Constante de transit individuelle (1/h)
KA     : Constante d'absorption orale individuelle (1/h)
