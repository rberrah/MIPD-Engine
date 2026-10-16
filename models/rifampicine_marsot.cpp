[SET] end=24, delta=0.1


[PARAM] @annotated
TVKa       : 1.15 : Absorption rate constant (h-1) - FIXED

TVCL_NOFUS : 13.7 : Clearance without Fusidic Acid (L/h)
TVV_NOFUS  : 61.1 : Volume without Fusidic Acid (L)

TVCL_FUS   : 5.1  : Clearance with Fusidic Acid (L/h)
TVV_FUS    : 23.8 : Volume with Fusidic Acid (L)

ETA1 : 0 : IIV CL (ωCL/F (variance) 0.531 )
ETA2 : 0 : IIV V (ωV/F (variance) 0.349)



[PARAM] @annotated @covariates
FUS : 0 : Fusidic Acid Co-medication (0 = Pas d acide fusidique, 1 = Co-administration)



[OMEGA]
0.531
0.349



[SIGMA]
0.0001
2.256 // σadd (mg/L) σ, variance of εij



[CMT] @annotated
TR    : Transit Compartment [ADM]
DEPOT : Absorption Compartment
CENT  : Central Compartment [OBS]



[MAIN]
double TVCL = (FUS == 1) ? TVCL_FUS : TVCL_NOFUS;
double TVV  = (FUS == 1) ? TVV_FUS  : TVV_NOFUS;

double CL = TVCL * exp(ETA1 + ETA(1));
double V  = TVV  * exp(ETA2 + ETA(2));
double Ka = TVKa;

double Ktr = TVKa;



[ODE]
dxdt_TR    = -Ktr * TR;

dxdt_DEPOT =  Ktr * TR
             - Ka * DEPOT;

dxdt_CENT  =  Ka * DEPOT
              - (CL/V) * CENT;



[TABLE]
double IPRED = CENT / V;
double DV = IPRED * (1 + EPS(1)) + EPS(2);

int i = 0;
while(DV < 0 && i < 1000) {
    simeps();
    DV = IPRED * (1 + EPS(1)) + EPS(2);
    ++i;
}



[CAPTURE] DV CL V Ka