import { promises as fs } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const outputDirectory = path.join(root, 'models');

const victimModels = {
  'tacrolimus_woillard_ddi.cpp': 'DEPOT_TAC',
  'tacrolimus_kamp_envarsus_ddi.cpp': 'DEPOT_TAC',
  'sirolimus_djebli_ddi.cpp': 'DEPOT_SIR',
  'everolimus_moes_ddi.cpp': 'DEPOT_EVE',
  'everolimus_terheine_ddi.cpp': 'DEPOT_EVE',
  'ciclosporine_press_ddi.cpp': 'DEPOT_CSA'
};

const adaptationNote = `// Standalone PK adaptation for TDM.
// The DDI mechanism placeholders were removed: this file describes baseline PK only.
// No patient data or pasted model is persisted by the application.
`;

async function normalizeVictimModel(file, administrationCompartment) {
  const target = path.join(outputDirectory, file);
  let code = await fs.readFile(target, 'utf8');
  code = code
    .replace(/^\/\/ Standalone TDM adaptation from DDI Manager\+\.\r?\n\/\/ The DDI mechanism placeholders were removed: this file describes baseline PK only\.\r?\n\/\/ No patient data or pasted model is persisted by the application\.\r?\n/, '')
    .replace(/^.*\{\{[A-Z_]+\}\}.*(?:\r?\n|$)/gm, '')
    .replace(/@covariate\b/g, '@covariates');
  const compartmentPattern = new RegExp(`^(\\s*${administrationCompartment}\\s*:[^\\r\\n]*?)(?:\\s*\\[ADM\\])?\\s*$`, 'm');
  code = code.replace(compartmentPattern, '$1 [ADM]');
  if (!code.includes('Standalone PK adaptation for TDM.')) code = `${adaptationNote}${code}`;
  await fs.writeFile(target, code, 'utf8');
}

function commonHeader({ title, citation, doi, population, limitation }) {
  return `$PROB
// ${title}
// Article: ${citation}
// DOI: ${doi || 'not available'}
// Population: ${population}
// Implementation: deterministic PK module adapted for standalone TDM.
// ${limitation || 'OMEGA and SIGMA below are engineering priors added for MAP compatibility; they were not estimated in the cited source.'}
`;
}

function oneCompartmentModel(definition) {
  return `${commonHeader(definition)}
$PARAM @annotated
TVCL : ${definition.cl} : Typical clearance (L/h)
TVV  : ${definition.v} : Typical volume (L)
TVKA : ${definition.ka} : First-order absorption rate (1/h)
TVF  : ${definition.f} : Oral bioavailability
ETA1 : 0 : MAP random effect on clearance

$OMEGA @annotated
IIV_CL : 0.09 : Engineering prior variance on clearance

$SIGMA @annotated
PROP : 0.04 : Engineering prior proportional residual variance
ADD  : 0.0  : Additive residual variance fixed to zero

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]

$MAIN
double CL = TVCL * exp(ETA1 + ETA(1));
double V = TVV;

$ODE
dxdt_GUT = -TVKA * GUT;
dxdt_CENT = TVKA * GUT * TVF - (CL / V) * CENT;

$TABLE
double IPRED = CENT / V;
double DV = IPRED * (1.0 + EPS(1));

$CAPTURE DV CL
`;
}

function twoCompartmentModel(definition) {
  return `${commonHeader(definition)}
$PARAM @annotated
TVCL : ${definition.cl} : Typical clearance (L/h)
TVV1 : ${definition.v1} : Typical central volume (L)
TVQ  : ${definition.q} : Typical intercompartmental clearance (L/h)
TVV2 : ${definition.v2} : Typical peripheral volume (L)
TVKA : ${definition.ka} : First-order absorption rate (1/h)
TVF  : ${definition.f} : Oral bioavailability
ETA1 : 0 : MAP random effect on clearance

${definition.covariates || ''}$OMEGA @annotated
IIV_CL : 0.09 : Engineering prior variance on clearance

$SIGMA @annotated
PROP : 0.04 : Engineering prior proportional residual variance
ADD  : 0.0  : Additive residual variance fixed to zero

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]
PERI : Peripheral compartment

$MAIN
double CL = TVCL * exp(ETA1 + ETA(1));
double V1 = ${definition.v1Expression || 'TVV1'};
double Q = TVQ;
double V2 = ${definition.v2Expression || 'TVV2'};

$ODE
dxdt_GUT = -TVKA * GUT;
dxdt_CENT = TVKA * GUT * TVF - ((CL + Q) / V1) * CENT + (Q / V2) * PERI;
dxdt_PERI = (Q / V1) * CENT - (Q / V2) * PERI;

$TABLE
double IPRED = CENT / V1;
double DV = IPRED * (1.0 + EPS(1));

$CAPTURE DV CL
`;
}

function voriconazoleModel() {
  return `${commonHeader({
    title: 'Voriconazole nonlinear PK model (van den Born et al.)',
    citation: 'van den Born DA et al. Int J Antimicrob Agents. 2023;61:106750.',
    doi: '10.1016/j.ijantimicag.2023.106750',
    population: '54 predominantly adult patients treated with voriconazole; mixed hospital setting.',
    limitation: 'Published IIV was retained for Vmax; the residual prior is an engineering adaptation for mapbayr.'
  })}
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
`;
}

function isavuconazoleDesaiModel() {
  return `$PLUGIN tad
${commonHeader({
    title: 'Isavuconazole population PK model (Desai et al.)',
    citation: 'Desai A et al. Antimicrob Agents Chemother. 2016;60:5483-5491.',
    doi: '10.1128/AAC.02819-15',
    population: 'Healthy adults and adults with invasive fungal infections in phase 1 and SECURE phase 3 trials.',
    limitation: 'The published time-dependent Weibull absorption is retained; the residual prior is an engineering adaptation for mapbayr.'
  })}
$PARAM @annotated
TVCL : 2.36 : Typical clearance (L/h)
TVV1 : 49.1 : Typical central volume (L)
TVV2 : 417 : Typical peripheral volume (L)
TVQ : 26.6 : Typical intercompartmental clearance (L/h)
KAMAX : 1.08 : Maximum absorption rate (1/h)
RA : 0.72 : Weibull scale parameter
GAM1 : 4.88 : Weibull shape parameter
EFF_BMI_V2 : 0.060 : Linear BMI effect on peripheral volume
BMI_REF : 24.8 : Reference BMI (kg/m2)
ETA1 : 0 : MAP random effect on clearance

$PARAM @annotated @covariates
BMI : 23.6 : Body mass index (kg/m2)

$OMEGA @annotated
IIV_CL : 0.09 : Engineering prior variance on clearance

$SIGMA @annotated
PROP : 0.04 : Engineering prior proportional residual variance
ADD  : 0.0  : Additive residual variance fixed to zero

$CMT @annotated
GUT : Oral depot [ADM]
CENT : Central compartment [OBS]
PERI : Peripheral compartment

$MAIN
double TDOSE = TIME - self.tad();
double CL = TVCL * exp(ETA1 + ETA(1));
double V2 = TVV2 * (1.0 + EFF_BMI_V2 * (BMI - BMI_REF));

$ODE
double TAD = SOLVERTIME - TDOSE;
double KA = 0.0;
if (TAD > 0.0) KA = KAMAX * (1.0 - exp(-pow(RA * TAD, GAM1)));
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - ((CL + TVQ) / TVV1) * CENT + (TVQ / V2) * PERI;
dxdt_PERI = (TVQ / TVV1) * CENT - (TVQ / V2) * PERI;

$TABLE
double IPRED = CENT / TVV1;
double DV = IPRED * (1.0 + EPS(1));

$CAPTURE DV CL
`;
}

const deterministicModels = {
  'voriconazole_vandenborn_ddi.cpp': voriconazoleModel(),
  'ritonavir_kappelhoff_ddi.cpp': oneCompartmentModel({
    title: 'Ritonavir PK module',
    citation: 'Kappelhoff BS et al. Br J Clin Pharmacol. 2005;59:174-182.',
    doi: '10.1111/j.1365-2125.2004.02241.x',
    population: 'Adults living with HIV-1; ritonavir used as booster or antiviral.',
    cl: 10, v: 80, ka: 1, f: 0.7
  }),
  'rifampicine_loos_ddi.cpp': oneCompartmentModel({
    title: 'Rifampicin chronic-dose PK module',
    citation: 'Loos U et al. Klin Wochenschr. 1985;63:1205-1211.',
    doi: '10.1007/BF01733779',
    population: 'Healthy adults receiving oral and intravenous rifampicin during chronic administration.',
    cl: 25, v: 50, ka: 1.5, f: 0.7
  }),
  'posaconazole_dolton_ddi.cpp': oneCompartmentModel({
    title: 'Posaconazole integrated PK module',
    citation: 'Dolton MJ et al. Antimicrob Agents Chemother. 2014;58:6879-6885.',
    doi: '10.1128/AAC.03777-14',
    population: 'Healthy adults and adult patients receiving posaconazole for prophylaxis or treatment.',
    cl: 5, v: 246, ka: 0.8, f: 0.54
  }),
  'isavuconazole_desai_ddi.cpp': isavuconazoleDesaiModel(),
  'isavuconazole_groll_ddi.cpp': twoCompartmentModel({
    title: 'Isavuconazole DDI study PK module',
    citation: 'Groll AH et al. Clin Pharmacol Drug Dev. 2017;6:76-85.',
    doi: '10.1002/cpdd.284',
    population: 'Healthy adults in controlled drug-drug interaction studies.',
    cl: 2.5, v1: 106, q: 6.5, v2: 250, ka: 1, f: 0.98,
    covariates: '$PARAM @annotated @covariates\nWT : 70 : Body weight (kg)\n\n',
    v1Expression: 'TVV1 * (WT / 70.0)'
  }),
  'fluconazole_debruyne_ddi.cpp': oneCompartmentModel({
    title: 'Fluconazole clinical PK module',
    citation: 'Debruyne D, Ryckelynck JP. Clin Pharmacokinet. 1993;24:10-27.',
    doi: '10.2165/00003088-199324010-00002',
    population: 'Adults from clinical pharmacokinetic studies summarized in the source review.',
    cl: 1, v: 50, ka: 1.5, f: 0.9
  }),
  'erythromycine_mather_ddi.cpp': oneCompartmentModel({
    title: 'Erythromycin oral PK module',
    citation: 'Mather LE et al. Br J Clin Pharmacol. 1981;12:131-140.',
    doi: '10.1111/j.1365-2125.1981.tb01191.x',
    population: 'Healthy adults in oral absorption and bioavailability studies.',
    cl: 20, v: 200, ka: 1, f: 0.6
  }),
  'diltiazem_hermann_ddi.cpp': oneCompartmentModel({
    title: 'Diltiazem oral PK module',
    citation: 'Hermann P et al. Eur J Clin Pharmacol. 1983;24:349-352.',
    doi: '10.1007/BF00610053',
    population: 'Healthy adults after intravenous and oral diltiazem.',
    cl: 40, v: 150, ka: 1, f: 0.8
  }),
  'cobicistat_mathias_ddi.cpp': oneCompartmentModel({
    title: 'Cobicistat 50 mg sensitivity module',
    citation: 'Mathias AA et al. Clin Pharmacol Ther. 2010;87:322-329.',
    doi: '10.1038/clpt.2009.228',
    population: 'Healthy adults in a dose-ranging pharmacokinetic study.',
    limitation: 'Validated only as a 50 mg once-daily sensitivity module; extrapolation to marketed 150 mg is not supported. OMEGA/SIGMA are engineering priors.',
    cl: 6.05, v: 26.2, ka: 1.5, f: 0.1
  }),
  'amiodarone_pollak_ddi.cpp': twoCompartmentModel({
    title: 'Long-term oral amiodarone PK module',
    citation: 'Pollak PT et al. Clin Pharmacol Ther. 2000;67:642-652.',
    doi: '10.1067/mcp.2000.107047',
    population: 'Adults receiving long-term oral amiodarone therapy for cardiac arrhythmia.',
    cl: 9, v1: 66, q: 40, v2: 3000, ka: 0.6, f: 0.5
  })
};

await fs.mkdir(outputDirectory, { recursive: true });
for (const [file, administrationCompartment] of Object.entries(victimModels)) {
  await normalizeVictimModel(file, administrationCompartment);
}
for (const [file, code] of Object.entries(deterministicModels)) {
  await fs.writeFile(path.join(outputDirectory, file), code, 'utf8');
}

console.log(`DDI TDM models ready: ${Object.keys(victimModels).length + Object.keys(deterministicModels).length} standalone files.`);
