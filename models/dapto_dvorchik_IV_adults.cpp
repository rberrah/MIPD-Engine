[SET] end=100, delta=0.1

[PARAM] @annotated
TVCL : 0.807 : Typical value of Clearance (L/h)
TVV1 : 4.80  : Typical value of Central volume (L)
TVQ  : 3.46  : Typical value of Intercompartmental clearance (L/h)
TVV2 : 3.13  : Typical value of Peripheral volume of distribution (L)
TVCL_DIAL : 0.269 : Typical clearance during dialysis (L/h)
FEMALE_CL : 0.801 : Female-to-male clearance ratio
INFECT_V2 : 1.93  : V2 multiplier in subjects with bacterial infection

CRCL : 0.00514 : effect of creatinine clearance on CL
TCL  : 0.14    : effect of body temperature on CL
WTQ  : 0.0593  : effect of body weight on Q
WTV2 : 0.0458  : effect of body weight on V2

ETA1 : 0 : IIV on clearance
ETA2 : 0 : IIV on central volume
ETA3 : 0 : IIV on intercompartmental Clearance
ETA4 : 0 : IIV on peripheral volume


[PARAM] @annotated @covariates
CREAT : 60 : estimated creatinine micromol
TEMP    : 37.2 : Body temperature (°C)
WT      : 75.1 : Body weight (kg)
SEX     : 0 : Sex (0=Male, 1=Female)   // corrected; see [MAIN]
AGE  : 62 : age (year)
DIAL : 0 : Dialysis status (0=No, 1=Yes)
INFECT : 1 : Bacterial infection (0=No, 1=Yes)


[OMEGA]
0.093636
0.3215
0.425104
0.0364810


[SIGMA]
0    // proportional
4.28 // additive variance (LC/MS/MS assay SD 2.07 mg/L)


[CMT] @annotated
CENT   : Central compartment (mg) [ADM, OBS]
PERIPH : Peripheral compartment (mg)


[MAIN]
// SEX CONVENTION. These two lines used to contradict each other. Cockcroft-Gault
// below applies 1.25 when SEX == 0, and 1.25 is the MALE constant, so this file
// treats 0 as male -- as do the 17 other models in this library. But the clearance
// multiplier was written `(0.8 + 0.2*SEX)`, which gives 0.8 when SEX == 0, i.e. it
// treated 0 as FEMALE. One of the two had to be wrong.
//
// Dvorchik settles it: "CL in females was 80% that in males", and the reported
// TVCL is the value "for male subject with median creatinine clearance". So the
// reference is male and the 0.8 belongs to females. Keeping 0 = male (correct for
// the Cockcroft-Gault line and consistent with the rest of the library), the
// multiplier must be `1.0 - 0.2*SEX`.
//
// Left as it was, every simulated patient received the sex effect backwards:
// males cleared daptomycin 20 % too slowly and females 25 % too fast.
double CCL = (SEX == 0 ? 1.25 : 1.04) * WT * (140 - AGE) / CREAT;
double CL_R = DIAL >= 0.5 ? TVCL_DIAL : TVCL + CRCL*(CCL - 91.2);
double sex_multiplier = SEX >= 0.5 ? FEMALE_CL : 1.0;
double CL = ((CL_R + TCL*(TEMP - 37.2)) * sex_multiplier) * exp(ETA1 + ETA(1)) ;
double V1 = TVV1 * exp(ETA2 + ETA(2)) ;
double Q  = (TVQ + WTQ * (WT - 75.1)) * exp(ETA3 + ETA(3)) ;
double infection_multiplier = INFECT >= 0.5 ? INFECT_V2 : 1.0;
double V2 = ((TVV2 + WTV2 * (WT - 75.1)) * infection_multiplier) * exp(ETA4 + ETA(4)) ;


[ODE]
dxdt_CENT   =  Q*PERIPH/V2 -CL*CENT/V1 -Q*CENT/V1 ;
dxdt_PERIPH =  Q*CENT/V1 -Q*PERIPH/V2 ;


[TABLE]
double DV = (CENT/V1) * (1 + EPS(1)) + EPS(2);
int i = 0;
while(DV<0 && i <100) {
simeps();
DV = (CENT/V1) * (1 + EPS(1)) + EPS(2);
++i;
}


[CAPTURE] DV CL V2 V1 Q CL_R sex_multiplier infection_multiplier
