$PROB
# Modèle : Amoxicilline PopPK (Fournier et al. 2018)
# Population : patients brûlés hospitalisés en réanimation.
# Structure : 2 compartiments, administration IV, élimination linéaire.
# Covariables :
#   - CLCR estimée par Cockcroft-Gault à partir de AGE, WT, CREAT, SEX.
#   - Effet de CLCR sur CL.
#   - Effet du poids corporel sur V1.
# Partie PD volontairement retirée : pas de MIC, pas de fT>MIC, pas d'AUC PD.
#
# Utilisation :
# - Administration IV dans CENT, cmt = 1.
# - Exemple : perfusion d'amoxicilline 1000 ou 2000 mg sur 0.5 h à 2 h.

$PARAM @annotated
TVCL       : 13.6  : Clairance typique à CLCR 110 mL/min (L/h)
TVV1       : 9.73  : Volume central typique à 70 kg (L)
TVQ        : 20.1  : Clairance intercompartimentale typique (L/h)
TVV2       : 17.6  : Volume périphérique typique (L)
CLCR_CL    : 0.53  : Effet proportionnel de CLCR sur CL
REF_CLCR   : 110.0 : CLCR de référence du modèle (mL/min)
REF_WT     : 70.0  : Poids de référence du modèle (kg)

// Mapbayr ETA
ETA1       : 0.0   : ETA posterior sur CL

$PARAM @covariates @annotated
AGE   : 50.0 : Âge du patient (années)
WT    : 70.0 : Poids corporel du patient (kg)
CREAT : 70.0 : Créatinine sérique (µmol/L)
SEX   : 0.0  : Sexe (0 = Homme, 1 = Femme)

$OMEGA @annotated @diagonal
ETA_CL : 0.1303 : Variabilité interindividuelle sur CL, variance approximée depuis CV 37.3%

$SIGMA @annotated @diagonal
PROP : 0.1369 : Erreur résiduelle proportionnelle, variance correspondant à 37%
ADD  : 0.0064 : Erreur résiduelle additive, variance correspondant à 0.08 mg/L

$CMT @annotated
CENT   : Compartiment central (mg) [ADM, OBS]
PERIPH : Compartiment périphérique (mg)

$MAIN
// -----------------------------------------------------------------------------
// 1. Sécurisation des covariables
// -----------------------------------------------------------------------------
double safe_AGE   = AGE;
double safe_WT    = WT;
double safe_CREAT = CREAT;

if(safe_AGE < 18.0)  safe_AGE = 18.0;
if(safe_AGE > 120.0) safe_AGE = 120.0;

if(safe_WT < 30.0)   safe_WT = 30.0;
if(safe_WT > 250.0)  safe_WT = 250.0;

if(safe_CREAT <= 0.0) safe_CREAT = 1.0;

// -----------------------------------------------------------------------------
// 2. Cockcroft-Gault en mL/min avec créatinine en µmol/L
//    SEX = 0 : homme
//    SEX = 1 : femme
// -----------------------------------------------------------------------------
double sex_factor = (SEX == 1.0) ? 0.85 : 1.0;

double clcr_calc = sex_factor * (1.23 * (140.0 - safe_AGE) * safe_WT) / safe_CREAT;

if(clcr_calc < 1.0)   clcr_calc = 1.0;
if(clcr_calc > 300.0) clcr_calc = 300.0;

// -----------------------------------------------------------------------------
// 3. Paramètres individuels
// -----------------------------------------------------------------------------
// Relation publiée : CL dépend linéairement de CLCR autour de 110 mL/min.
// Interprétation utilisée : effet proportionnel centré sur REF_CLCR.
double clcr_effect = 1.0 + CLCR_CL * ((clcr_calc - REF_CLCR) / REF_CLCR);

if(clcr_effect < 0.10) clcr_effect = 0.10;

// Effet du poids sur V1
double wt_effect_v1 = safe_WT / REF_WT;

if(wt_effect_v1 < 0.10) wt_effect_v1 = 0.10;

// Paramètres individuels
double CL = TVCL * clcr_effect * exp(ETA(1) + ETA1);
double V1 = TVV1 * wt_effect_v1;
double Q  = TVQ;
double V2 = TVV2;

// Sécurisation numérique
if(CL < 0.001) CL = 0.001;
if(V1 < 0.001) V1 = 0.001;
if(Q  < 0.001) Q  = 0.001;
if(V2 < 0.001) V2 = 0.001;

// Constantes de transfert
double K10 = CL / V1;
double K12 = Q  / V1;
double K21 = Q  / V2;

$ODE
dxdt_CENT   = K21 * PERIPH - (K10 + K12) * CENT;
dxdt_PERIPH = K12 * CENT - K21 * PERIPH;

$TABLE
double CP = CENT / V1;
double DV = CP * (1.0 + EPS(1)) + EPS(2);

if(DV < 0.0) DV = 0.0;

$CAPTURE @annotated
DV        : Variable dépendante simulée, concentration totale (mg/L)
CP        : Concentration totale prédite sans erreur résiduelle (mg/L)
CL        : Clairance individuelle (L/h)
V1        : Volume central individuel (L)
Q         : Clairance intercompartimentale individuelle (L/h)
V2        : Volume périphérique individuel (L)
K10       : Constante d'élimination individuelle (1/h)
K12       : Constante de transfert central vers périphérique (1/h)
K21       : Constante de transfert périphérique vers central (1/h)
clcr_calc : Clairance de la créatinine estimée par Cockcroft-Gault (mL/min)