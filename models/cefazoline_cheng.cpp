$PROB
Cefazolin PopPK Model (Cheng et al. 2025)
Reference: Cheng V et al. Anaesth Crit Care Pain Med. 2025; 45:10653.
Population: Pediatric cardiac surgery patients with CPB
Structure: 2 compartments, separate V1 for intra/post-op.
The published inter-occasion variability on CL and V1 is omitted because the
current MAP workflow does not represent peri-operative occasions with
independent random effects.
$PARAM @annotated
TVCL        : 1.84  : Typical Clearance (L/h)
TVV1_INTRA  : 1.76  : Typical Central Volume Intra-op (L)
TVV1_POST   : 2.77  : Typical Central Volume Post-op (L)
TVQ         : 3.07  : Typical Inter-compartmental Clearance (L/h)
TVV2        : 4.47  : Typical Peripheral Volume (L)
CL_WT_EXP   : 0.75  : Allometric exponent of WT on CL
CL_eGFR_EXP : 1.50  : Exponent of eGFR effect on CL
Q_WT_EXP    : 0.78  : Allometric exponent of WT on Q
REF_WT      : 6.8   : Reference WT (kg)
REF_eGFR    : 77.3  : Reference eGFR (mL/min/1.73m2)
ETA1 : 0.0 : Mapbayr ETA BSV on CL
ETA2 : 0.0 : Mapbayr ETA BSV on V1
ETA3 : 0.0 : Mapbayr ETA BSV on Q
ETA4 : 0.0 : Mapbayr ETA BSV on V2
$PARAM @annotated @covariates
WT     : 6.8  : Body Weight (kg)
HT     : 65.0 : Height (cm)
CREAT  : 35.0 : Serum Creatinine (micromol/L)
PERIOD : 0.0  : Surgical Period (0 = Intra-op [0-6h], 1 = Post-op [6-30h])
$CMT @annotated
CENT   : Central Compartment (mg) [ADM, OBS]
PERIPH : Peripheral Compartment (mg)
$OMEGA @annotated @diagonal
ETA_CL : 0.033649 : Variance BSV on CL
ETA_V1 : 0.082319 : Variance BSV on V1
ETA_Q  : 0.436228 : Variance BSV on Q
ETA_V2 : 0.280802 : Variance BSV on V2
$SIGMA @annotated @diagonal
PROP : 0.078400 : Proportional error variance
ADD  : 0.000000 : Additive error variance
$MAIN
double creat_safe = (CREAT <= 0.0) ? 1.0  : CREAT;
double ht_safe    = (HT <= 0.0)    ? 50.0 : HT;
double wt_safe    = (WT <= 0.0)    ? 0.1  : WT;
double creat_mgdl = creat_safe / 88.4;
double eGFR_calc  = 0.413 * ht_safe / creat_mgdl;
eGFR_calc = (eGFR_calc < 1.0)   ? 1.0   : eGFR_calc;
eGFR_calc = (eGFR_calc > 400.0) ? 400.0 : eGFR_calc;
double TVV1 = (PERIOD >= 1.0) ? TVV1_POST : TVV1_INTRA;
double CL = TVCL * pow(wt_safe / REF_WT, CL_WT_EXP) * pow(eGFR_calc / REF_eGFR, CL_eGFR_EXP) * exp(ETA(1) + ETA1);
double V1 = TVV1 * (wt_safe / REF_WT) * exp(ETA(2) + ETA2);
double Q  = TVQ  * pow(wt_safe / REF_WT, Q_WT_EXP) * exp(ETA(3) + ETA3);
double V2 = TVV2 * (wt_safe / REF_WT) * exp(ETA(4) + ETA4);
$ODE
dxdt_CENT   = -(CL/V1)*CENT - (Q/V1)*CENT + (Q/V2)*PERIPH;
dxdt_PERIPH =  (Q/V1)*CENT - (Q/V2)*PERIPH;
$TABLE
double CP = CENT / V1;
double DV = CP * (1.0 + EPS(1)) + EPS(2);
if(DV < 0.0) DV = 0.0;
$CAPTURE @annotated
CP         : Predicted Unbound Concentration (mg/L)
CL         : Individual Clearance (L/h)
V1         : Individual Central Volume (L)
Q          : Individual Inter-compartmental Clearance (L/h)
V2         : Individual Peripheral Volume (L)
eGFR_calc  : Calculated eGFR (mL/min/1.73m2)
DV         : Simulated Unbound Concentration (mg/L)
