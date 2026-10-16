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
${definition.lag == null ? '' : `TVLAG : ${definition.lag} : Absorption lag time (h)\n`}${definition.additionalParameters || ''}ETA1 : 0 : MAP random effect on clearance

${definition.covariates || ''}$OMEGA @annotated
IIV_CL : ${definition.omegaCl ?? 0.09} : ${definition.omegaLabel || 'Engineering prior variance on clearance'}

$SIGMA @annotated
PROP : ${definition.sigmaProp ?? 0.04} : ${definition.sigmaPropLabel || 'Engineering prior proportional residual variance'}
ADD  : ${definition.sigmaAdd ?? 0.0} : ${definition.sigmaAddLabel || 'Additive residual variance fixed to zero'}

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]

$MAIN
double CL = ${definition.clExpression || 'TVCL * exp(ETA1 + ETA(1))'};
double V = ${definition.vExpression || 'TVV'};
double KA = ${definition.kaExpression || 'TVKA'};
double F = ${definition.fExpression || 'TVF'};
${definition.lag == null ? '' : 'ALAG_GUT = TVLAG;\n'}
$ODE
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT * F - (CL / V) * CENT;

$TABLE
double IPRED = CENT / V;
double DV = IPRED * (1.0 + EPS(1)) + EPS(2);

$CAPTURE DV CL
`;
}

function amiodaronePollakModel() {
  return `${commonHeader({
    title: 'Long-term oral amiodarone population PK model',
    citation: 'Pollak PT et al. Clin Pharmacol Ther. 2000;67:642-652.',
    doi: '10.1067/mcp.2000.107047',
    population: 'Adults receiving long-term oral amiodarone therapy for cardiac arrhythmia.',
    limitation: 'Published apparent disposition parameters and IIV are retained. The accessible source does not identify an absorption model, so first-order ka is an explicit simulation assumption; this is not a MAP prior.'
  })}
$PARAM @annotated
TVCL : 9.541667 : Published CL/F, 229 L/day converted to L/h
TVV1 : 882 : Published central V/F (L)
TVQ  : 24.958333 : Published intercompartmental CL/F, 599 L/day converted to L/h
TVV2 : 12700 : Published peripheral V/F (L)
TVKA : 0.6 : Assumed first-order absorption rate (1/h)
ETA1 : 0 : Random effect on CL/F
ETA2 : 0 : Random effect on Q/F
ETA3 : 0 : Random effect on V2/F

$OMEGA @annotated
IIV_CL : 0.091755 : log(1 + 0.31^2)
IIV_Q  : 0.272423 : log(1 + 0.56^2)
IIV_V2 : 0.289783 : log(1 + 0.58^2)

$SIGMA @annotated
PROP : 0.04 : Engineering residual variance; not reported in the accessible source
ADD  : 0.0 : Additive residual variance fixed to zero

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]
PERI : Peripheral compartment

$MAIN
double CL = TVCL * exp(ETA1 + ETA(1));
double Q = TVQ * exp(ETA2 + ETA(2));
double V2 = TVV2 * exp(ETA3 + ETA(3));

$ODE
dxdt_GUT = -TVKA * GUT;
dxdt_CENT = TVKA * GUT - ((CL + Q) / TVV1) * CENT + (Q / V2) * PERI;
dxdt_PERI = (Q / TVV1) * CENT - (Q / V2) * PERI;

$TABLE
double IPRED = CENT / TVV1;
double DV = IPRED * (1.0 + EPS(1)) + EPS(2);

$CAPTURE DV CL Q V2
`;
}

function ritonavirKappelhoffModel() {
  return `${commonHeader({
    title: 'Ritonavir population PK model',
    citation: 'Kappelhoff BS et al. Br J Clin Pharmacol. 2005;59:174-182.',
    doi: '10.1111/j.1365-2125.2004.02241.x',
    population: 'Adults living with HIV-1; ritonavir used as booster or antiviral.',
    limitation: 'Published structure, typical values, lag, lopinavir covariate and IIV are retained. IOV on F, the V-ka covariance and the two-population residual mixture are omitted; the dominant residual component is used, so this remains simulation-only.'
  })}
$PARAM @annotated
TVCL : 10.5 : Apparent clearance without lopinavir (L/h)
TVV  : 96.6 : Apparent volume (L)
TVKA : 0.871 : First-order absorption rate (1/h)
TVLAG : 0.778 : Absorption lag time (h)
LPV_CL : 2.72 : Clearance multiplier with concomitant lopinavir
ETA1 : 0 : Random effect on CL/F
ETA2 : 0 : Random effect on V/F
ETA3 : 0 : Random effect on ka

$PARAM @annotated @covariates
LPV : 0 : Concomitant lopinavir (0 = no, 1 = yes)

$OMEGA @annotated
IIV_CL : 0.136671 : log(1 + 0.383^2)
IIV_V  : 0.494696 : log(1 + 0.800^2)
IIV_KA : 1.350319 : log(1 + 1.69^2)

$SIGMA @annotated
PROP : 0.023716 : Published proportional SD 0.154 squared
ADD  : 0.0036 : Dominant published additive SD 0.060 mg/L squared

$CMT @annotated
GUT  : Oral depot [ADM]
CENT : Central compartment [OBS]

$MAIN
double CL = TVCL * (LPV > 0.5 ? LPV_CL : 1.0) * exp(ETA1 + ETA(1));
double V = TVV * exp(ETA2 + ETA(2));
double KA = TVKA * exp(ETA3 + ETA(3));
ALAG_GUT = TVLAG;

$ODE
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - (CL / V) * CENT;

$TABLE
double IPRED = CENT / V;
double DV = IPRED * (1.0 + EPS(1)) + EPS(2);

$CAPTURE DV CL V KA
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
    limitation: 'Published Weibull absorption, race and patient-status effects, BMI effect and reported IIV are retained. Patient CL variability is the default and the healthy/patient residual branches are not separated, so this remains simulation-only.'
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
ETA2 : 0 : Random effect on Q
ETA3 : 0 : Random effect on peripheral volume
ETA4 : 0 : Random effect on Weibull RA
ETA5 : 0 : Random effect on Weibull GAM1

$PARAM @annotated @covariates
BMI : 23.6 : Body mass index (kg/m2)
ASIAN : 0 : Asian race indicator
PATIENT : 1 : SECURE patient (1) or healthy phase-1 participant (0)

$OMEGA @annotated
IIV_CL : 0.3290 : Patient CL IIV, log(1 + 0.6244^2)
IIV_Q : 0.2158 : Q IIV, log(1 + 0.4909^2)
IIV_V2 : 0.0962 : Vp IIV, log(1 + 0.3178^2)
IIV_RA : 0.1502 : RA IIV, log(1 + 0.4024^2)
IIV_GAM1 : 0.1900 : GAM1 IIV, log(1 + 0.4571^2)

$SIGMA @annotated
PROP : 0.201960 : Published exponential residual SD 0.4494 squared
ADD  : 0.0  : Additive residual variance fixed to zero

$CMT @annotated
GUT : Oral depot [ADM]
CENT : Central compartment [OBS]
PERI : Peripheral compartment

$MAIN
double TDOSE = TIME - self.tad();
double CL = (ASIAN > 0.5 ? 1.51 : TVCL) * exp(ETA1 + ETA(1));
double Q = TVQ * exp(ETA2 + ETA(2));
double V2_BASE = PATIENT > 0.5 ? TVV2 : 260.0;
double V2 = V2_BASE * (1.0 + EFF_BMI_V2 * (BMI - BMI_REF)) * exp(ETA3 + ETA(3));
double RA_IND = RA * exp(ETA4 + ETA(4));
double GAM1_IND = GAM1 * exp(ETA5 + ETA(5));
if (V2 < 1.0) V2 = 1.0;

$ODE
double TAD = SOLVERTIME - TDOSE;
double KA = 0.0;
if (TAD > 0.0) KA = KAMAX * (1.0 - exp(-pow(RA_IND * TAD, GAM1_IND)));
dxdt_GUT = -KA * GUT;
dxdt_CENT = KA * GUT - ((CL + Q) / TVV1) * CENT + (Q / V2) * PERI;
dxdt_PERI = (Q / TVV1) * CENT - (Q / V2) * PERI;

$TABLE
double IPRED = CENT / TVV1;
double DV = IPRED * (1.0 + EPS(1));

$CAPTURE DV CL
`;
}

const deterministicModels = {
  'voriconazole_vandenborn_ddi.cpp': voriconazoleModel(),
  'ritonavir_kappelhoff_ddi.cpp': ritonavirKappelhoffModel(),
  'rifampicine_loos_ddi.cpp': oneCompartmentModel({
    title: 'Rifampicin chronic-dose mean-profile reconstruction',
    citation: 'Loos U et al. Klin Wochenschr. 1985;63:1205-1211.',
    doi: '10.1007/BF01733779',
    population: 'Healthy adults receiving oral and intravenous rifampicin during chronic administration.',
    limitation: 'The published chronic-dose clearance and bioavailability are retained. Volume and ka were not reported in the accessible source summary and remain explicit profile-shape assumptions; this is not a PopPK prior.',
    cl: 9.03, v: 50, ka: 1.5, f: 0.68
  }),
  'posaconazole_dolton_ddi.cpp': oneCompartmentModel({
    title: 'Posaconazole integrated population PK model',
    citation: 'Dolton MJ et al. Antimicrob Agents Chemother. 2014;58:6879-6885.',
    doi: '10.1128/AAC.03777-14',
    population: 'Healthy adults and adult patients receiving posaconazole for prophylaxis or treatment.',
    limitation: 'Published structure, fixed effects, lag and medication/clinical covariates are retained. Patient-study IIV and residual variance are used as the default; the study-specific error branches are not represented, so this remains simulation-only.',
    cl: 30.2, v: 1100, ka: 1.26, f: 1, lag: 1.79,
    additionalParameters: 'PHE_CL : 7.21 : CL/F multiplier with phenytoin\nRIF_CL : 7.21 : CL/F multiplier with rifampicin\nFOS_CL : 1.342 : CL/F multiplier with fosamprenavir\n',
    covariates: '$PARAM @annotated @covariates\nPHE : 0 : Phenytoin coadministration\nRIF : 0 : Rifampicin coadministration\nFOS : 0 : Fosamprenavir coadministration\nPPI : 0 : Proton-pump inhibitor\nMET : 0 : Metoclopramide\nNUT : 0 : Nutritional supplement\nMUC : 0 : Mucositis\nDIARRHEA : 0 : Diarrhea\n\n',
    clExpression: 'TVCL * pow(PHE_CL, PHE) * pow(RIF_CL, RIF) * pow(FOS_CL, FOS) * exp(ETA1 + ETA(1))',
    fExpression: 'TVF * pow(0.549, PPI) * pow(0.655, MET) * pow(2.29, NUT) * pow(0.423, MUC) * pow(0.549, DIARRHEA)',
    omegaCl: 0.285156,
    omegaLabel: 'Patient-study IIV on CL/F, 53.4 percent squared',
    sigmaProp: 0.289444,
    sigmaPropLabel: 'Patient-study proportional residual SD 53.8 percent squared'
  }),
  'isavuconazole_desai_ddi.cpp': isavuconazoleDesaiModel(),
  'isavuconazole_groll_ddi.cpp': oneCompartmentModel({
    title: 'Isavuconazole DDI-study mean-profile reconstruction',
    citation: 'Groll AH et al. Clin Pharmacol Drug Dev. 2017;6:76-85.',
    doi: '10.1002/cpdd.284',
    population: 'Healthy adults in controlled drug-drug interaction studies.',
    limitation: 'The cited study reports noncompartmental DDI profiles, not a population model. CL, V and ka below are rounded mean-profile reconstruction values for simulation only; no claim of article-model reproduction is made.',
    cl: 2.0, v: 300, ka: 1, f: 1
  }),
  'fluconazole_debruyne_ddi.cpp': oneCompartmentModel({
    title: 'Fluconazole review-based mean-profile reconstruction',
    citation: 'Debruyne D, Ryckelynck JP. Clin Pharmacokinet. 1993;24:10-27.',
    doi: '10.2165/00003088-199324010-00002',
    population: 'Adults from clinical pharmacokinetic studies summarized in the source review.',
    limitation: 'The source is a clinical PK review rather than a PopPK model. Mean clearance (1.17 L/h), bioavailability (0.90) and the volume implied by the reported 31.6 h half-life are used; ka and random effects are simulation assumptions.',
    cl: 1.17, v: 53.36, ka: 1.5, f: 0.9
  }),
  'erythromycine_mather_ddi.cpp': oneCompartmentModel({
    title: 'Erythromycin oral mean-profile reconstruction',
    citation: 'Mather LE et al. Br J Clin Pharmacol. 1981;12:131-140.',
    doi: '10.1111/j.1365-2125.1981.tb01191.x',
    population: 'Healthy adults in oral absorption and bioavailability studies.',
    limitation: 'The source compares formulation-, meal- and occasion-dependent bioavailability after prior individual IV characterization; it does not publish a reusable PopPK model. F=0.40 represents erythromycin base after meals, while CL, V and ka remain explicit simulation assumptions.',
    cl: 20, v: 200, ka: 1, f: 0.4
  }),
  'diltiazem_hermann_ddi.cpp': oneCompartmentModel({
    title: 'Diltiazem oral mean-profile reconstruction',
    citation: 'Hermann P et al. Eur J Clin Pharmacol. 1983;24:349-352.',
    doi: '10.1007/BF00610053',
    population: 'Healthy adults after intravenous and oral diltiazem.',
    limitation: 'The source reports biphasic IV disposition and summary NCA values, not a reusable PopPK model. The 70 kg mean CL, apparent V and oral F are retained in a one-compartment terminal-profile approximation; ka is assumed.',
    cl: 89.6, v: 371, ka: 1, f: 0.42
  }),
  'cobicistat_mathias_ddi.cpp': oneCompartmentModel({
    title: 'Cobicistat 50 mg sensitivity module',
    citation: 'Mathias AA et al. Clin Pharmacol Ther. 2010;87:322-329.',
    doi: '10.1038/clpt.2009.228',
    population: 'Healthy adults in a dose-ranging pharmacokinetic study.',
    limitation: 'The source is a nonlinear dose-ranging NCA study, not a PopPK model. The published 50 mg steady-state apparent clearance and terminal volume are retained for this dose-specific sensitivity profile; ka and random effects are assumptions and extrapolation to 150 mg is unsupported.',
    cl: 154, v: 316, ka: 0.8, f: 1
  }),
  'amiodarone_pollak_ddi.cpp': amiodaronePollakModel()
};

await fs.mkdir(outputDirectory, { recursive: true });
for (const [file, administrationCompartment] of Object.entries(victimModels)) {
  await normalizeVictimModel(file, administrationCompartment);
}
for (const [file, code] of Object.entries(deterministicModels)) {
  await fs.writeFile(path.join(outputDirectory, file), code, 'utf8');
}

console.log(`DDI TDM models ready: ${Object.keys(victimModels).length + Object.keys(deterministicModels).length} standalone files.`);
