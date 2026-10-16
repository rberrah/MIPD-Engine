$PROB
# Modèle : Cefepime PopPK (An et al. 2023)
# Patient en état critique.
# Structure : 2 compartiments, Perfusion IV

$PARAM @annotated
TVCLR    : 2.00   : Clairance Renale Typique (L/h)
TVCLCRRT : 1.64   : Clairance sous CRRT (L/h)
TVCLNR   : 0.526  : Clairance Non-Renale Typique (L/h)
TVVC     : 13.4   : Volume Central Typique (L) à 90 kg TBW
TVVP     : 7.52   : Volume Périphérique Typique (L)
TVQ      : 12.0   : Clairance Inter-compartimentale Typique (L/h)
REF_CLCR : 54.0   : CLCR_LBW de référence (mL/min)
REF_WT   : 90.0   : Poids Total (TBW) de référence (kg)

ETA1 : 0 : CL
ETA2 : 0 : V1
ETA3 : 0 : V2

$PARAM @covariates @annotated
AGE   : 65   : Âge du patient (années)
WT    : 90   : Poids Total du patient - TBW (kg)
HT    : 170  : Taille du patient (cm)
CREAT : 70   : Créatinine Sérique (µmol/L)
SEX   : 0    : Sexe (0 = Homme, 1 = Femme)
CRRT  : 0    : Epuration continue (0 = Non, 1 = Oui)

$OMEGA @annotated @diagonal
// Variances approximées par ln(CV^2 + 1)
ETA_CLT : 0.08563 : IIV Clairance Totale (CV 29.9%)
ETA_VC  : 0.35608 : IIV Volume Central (CV 65.4%)
ETA_VP  : 0.02949 : IIV Volume Périphérique (CV 17.3%)

$SIGMA @annotated
PROP : 0.050176 : Erreur résiduelle proportionnelle (CV 22.4%)
ADD  : 41.4736 : Variance résiduelle additive (SD 6.44 mg/L)

$CMT @annotated
CENT   : Compartiment Central (mg) [ADM, OBS]
PERIPH : Compartiment Périphérique (mg)

$MAIN
// 1. Poids Maigre (LBW) : formule de Boer
double LBW = 0.0;
if (SEX == 0) { LBW = (0.407 * WT) + (0.267 * HT) - 19.2; }
else { LBW = (0.252 * WT) + (0.473 * HT) - 48.3; }

if (LBW < 30.0) { LBW = 30.0; }

double sex_factor = (SEX == 1) ? 0.85 : 1.0;
double clcr_calc = sex_factor * (1.23 * (140.0 - AGE) * LBW) / CREAT;

double on_crrt = CRRT >= 0.5 ? 1.0 : 0.0;
double TVCL = TVCLR * (clcr_calc / REF_CLCR) * (1.0 - on_crrt) + TVCLNR + TVCLCRRT * on_crrt;
double TVV1 = TVVC  * (WT / REF_WT);

// Individuels
double CL = TVCL * exp(ETA(1) + ETA1);
double VC = TVV1 * exp(ETA(2) + ETA2);
double VP = TVVP * exp(ETA(3) + ETA3);
double Q  = TVQ;

// Constantes de transfert
double k10 = CL / VC;
double k12 = Q  / VC;
double k21 = Q  / VP;

$ODE
dxdt_CENT   = -(k10 + k12) * CENT + k21 * PERIPH;
dxdt_PERIPH =   k12 * CENT - k21 * PERIPH;

$TABLE
double CP = CENT / VC;
double DV = CP * (1 + EPS(1)) + EPS(2);

$CAPTURE @annotated
DV            : DV
CL            : Clairance Individuelle (L/h)
VC            : Volume Central Individuel (L)
CP            : Concentration Plasmatique (mg/L)
LBW           : Poids Maigre Calculé en Interne (kg)
clcr_calc     : CLCR calculée par Cockcroft-Gault sur LBW (mL/min)
on_crrt       : Indicateur CRRT utilisé dans la clairance
