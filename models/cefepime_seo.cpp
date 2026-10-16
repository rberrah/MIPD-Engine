$PROB
# Modèle : Cefepime PopPK (Seo et al. 2023)
# Patient avec pneumonie nosocomiale en état critique.
# Structure : 2 compartiments, Perfusion IV, Solution Analytique (Optimisée pour VPC)

$PARAM @annotated
TVCL     : 6.60  : Clairance Typique (L/h)
TVV1     : 13.3  : Volume Central Typique (L)
TVQ      : 16.5  : Clairance Inter-compartimentale (L/h)
TVV2     : 13.0  : Volume Périphérique (L)
CL_EXP   : 0.656 : Exposant pour l'effet de la CLCR sur CL
REF_CLCR : 77.21 : CLCR de référence dans l'étude (mL/min)

// Mapbayr ETAs (Required for use_posterior)
ETA1 : 0.0 : ETA on clearance
ETA2 : 0.0 : ETA on V1
ETA3 : 0.0 : ETA on Q
ETA4 : 0.0 : ETA on V2

$PARAM @covariates @annotated
AGE   : 67.0 : Âge du patient (années)
WT    : 60.0 : Poids du patient (kg)
CREAT : 70.0 : Créatinine Sérique (µmol/L)
SEX   : 0.0  : Sexe (0 = Homme, 1 = Femme)

$OMEGA @annotated @diagonal
// Variances approximées par ln(CV^2 + 1)
ETA_CL : 0.10757 : Variabilité Inter-Individuelle CL (CV 33.7%)
ETA_V1 : 0.11000 : Variabilité Inter-Individuelle V1 (CV 34.1%)
ETA_Q  : 0.22957 : Variabilité Inter-Individuelle Q (CV 50.8%)
ETA_V2 : 0.14911 : Variabilité Inter-Individuelle V2 (CV 40.1%)

$SIGMA @annotated @diagonal
PROP : 0.005806 : Erreur résiduelle proportionnelle (CV 7.62%)
ADD  : 0.000000 : Erreur résiduelle additive

$CMT @annotated
CENT   : Compartiment Central (mg) [ADM, OBS]
PERIPH : Compartiment Périphérique (mg)

$MAIN
// 1. Cockcroft-Gault (Sécurisé contre la division par zéro et les valeurs négatives)
double safe_CREAT = (CREAT <= 0.0) ? 1.0 : CREAT;
double sex_factor = (SEX == 1.0) ? 0.85 : 1.0;
double clcr_calc = sex_factor * (1.23 * (140.0 - AGE) * WT) / safe_CREAT;

// Limitation des valeurs extrêmes de covariables pour éviter les instabilités numériques
if(clcr_calc < 1.0) clcr_calc = 1.0;

// Paramètres Individuels
double effet_cov = pow(clcr_calc / REF_CLCR, CL_EXP);

// mrgsolve $PKMODEL requiert spécifiquement la nomenclature CL, V1, Q, et V2.
double CL = TVCL * exp(ETA(1) + ETA1) * effet_cov;
double V1 = TVV1 * exp(ETA(2) + ETA2);
double Q  = TVQ  * exp(ETA(3) + ETA3);
double V2 = TVV2 * exp(ETA(4) + ETA4);

double K10 = CL / V1;
double K12 = Q / V1;
double K21 = Q / V2;

$ODE
// Les dxdt_ utilisent les constantes calculées dans $MAIN
dxdt_CENT   = K21 * PERIPH - (K10 + K12) * CENT;
dxdt_PERIPH = K12 * CENT - K21 * PERIPH;

$TABLE
double CP = CENT / V1;
double DV = CP * (1.0 + EPS(1)) + EPS(2);

// Prévention stricte des concentrations négatives générées par l'erreur résiduelle stochastique
if(DV < 0.0) DV = 0.0;

$CAPTURE @annotated
DV        : Variable Dépendante Simulée (mg/L)
CL        : Clairance Individuelle (L/h)
V1        : Volume Central Individuel (L)
CP        : Concentration Plasmatique (mg/L)
clcr_calc : CLCR calculée par Cockcroft-Gault (mL/min)
