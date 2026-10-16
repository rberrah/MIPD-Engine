import { promises as fs } from 'node:fs';
import { createHash } from 'node:crypto';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

await import('./generate_ddi_models.mjs');

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const modelsDirectory = path.join(root, 'models');
const outputFile = path.join(root, 'models', 'catalog.json');
const metadataFile = path.join(root, 'models', 'metadata.json');
const englishMetadataFile = path.join(root, 'models', 'i18n.en.json');
const mlRegistryFile = path.join(root, 'ml', 'registry.json');
const mlBenchmarkFile = path.join(root, 'ml', 'catalog-summary.generated.json');

const drugLabels = {
  amik: 'Amikacine',
  amox: 'Amoxicilline',
  cefazoline: 'Cefazoline',
  cefepime: 'Cefepime',
  ciclosporine: 'Ciclosporine',
  cobicistat: 'Cobicistat',
  dapto: 'Daptomycine',
  diltiazem: 'Diltiazem',
  erythromycine: 'Érythromycine',
  everolimus: 'Évérolimus',
  fluconazole: 'Fluconazole',
  genta: 'Gentamicine',
  isavuconazole: 'Isavuconazole',
  levo: 'Levofloxacine',
  linez: 'Linezolide',
  posaconazole: 'Posaconazole',
  rifampicine: 'Rifampicine',
  ritonavir: 'Ritonavir',
  sirolimus: 'Sirolimus',
  tacrolimus: 'Tacrolimus',
  voriconazole: 'Voriconazole',
  amiodarone: 'Amiodarone',
  vanco: 'Vancomycine'
};

const drugLabelsEn = {
  amik: 'Amikacin',
  amox: 'Amoxicillin',
  cefazoline: 'Cefazolin',
  cefepime: 'Cefepime',
  ciclosporine: 'Cyclosporine',
  cobicistat: 'Cobicistat',
  dapto: 'Daptomycin',
  diltiazem: 'Diltiazem',
  erythromycine: 'Erythromycin',
  everolimus: 'Everolimus',
  fluconazole: 'Fluconazole',
  genta: 'Gentamicin',
  isavuconazole: 'Isavuconazole',
  levo: 'Levofloxacin',
  linez: 'Linezolid',
  posaconazole: 'Posaconazole',
  rifampicine: 'Rifampicin',
  ritonavir: 'Ritonavir',
  sirolimus: 'Sirolimus',
  tacrolimus: 'Tacrolimus',
  voriconazole: 'Voriconazole',
  amiodarone: 'Amiodarone',
  vanco: 'Vancomycin'
};

const populationTagEn = {
  'Adulte': 'Adult',
  'Âgé': 'Older adult',
  'Brûlé': 'Burn',
  'Cancer': 'Cancer',
  'Cardiologie': 'Cardiology',
  'CEC': 'CPB',
  'Chirurgie': 'Surgery',
  'Chirurgie orthopédique': 'Orthopedic surgery',
  'CRRT': 'CRRT',
  'Diabète': 'Diabetes',
  'Dialyse': 'Dialysis',
  'Endocardite': 'Endocarditis',
  'Hémodialyse': 'Hemodialysis',
  'Hors dialyse': 'Not on dialysis',
  'Hors ICU': 'Non-ICU',
  'Hors ICU non précisé': 'ICU status not specified',
  'Hospitalisé': 'Hospitalized',
  'ICU': 'ICU',
  'Infection fongique': 'Fungal infection',
  'Infection ostéo-articulaire': 'Bone and joint infection',
  'Insuffisance rénale': 'Renal impairment',
  'Néonatal': 'Neonatal',
  'Obèse': 'Obesity',
  'Pédiatrique': 'Pediatric',
  'Pneumonie': 'Pneumonia',
  'Population à confirmer': 'Population pending verification',
  'Population mixte': 'Mixed population',
  'Sepsis': 'Sepsis',
  'Transplantation rénale': 'Kidney transplantation',
  'Tuberculose': 'Tuberculosis',
  'Ventilé': 'Mechanically ventilated',
  'VIH': 'HIV',
  'Volontaire sain': 'Healthy volunteer'
};

const modelTypeEn = {
  'Module de sensibilité déterministe': 'Deterministic sensitivity module',
  'Module PK adapté': 'Adapted PK module',
  'Module PK déterministe': 'Deterministic PK module',
  'PK publié': 'Published PK',
  'PK publié adapté': 'Adapted published PK',
  'PopPK à vérifier': 'PopPK pending verification',
  'PopPK mécanistique adapté': 'Adapted mechanistic PopPK',
  'PopPK publié': 'Published PopPK',
  'PopPK publié adapté': 'Adapted published PopPK',
  'Prior PK à vérifier': 'PK prior pending verification'
};

const administrationModeLabels = {
  ORAL: { fr: 'Orale', en: 'Oral' },
  IV_INTERMITTENT: { fr: 'IV intermittente', en: 'Intermittent IV' },
  IV_CONTINUOUS: { fr: 'IV continue', en: 'Continuous IV' }
};

const additionallyAdaptedModels = new Set([
  'genta_franck',
  'genta_gomes',
  'genta_hodiamont',
  'genta_rosario_clcr',
  'genta_rosario_clinical',
  'levo_gergs',
  'linez_sasaki'
]);

const drugKeys = Object.keys(drugLabels);

function titleCase(part) {
  if (part.toUpperCase() === 'IV') return 'IV';
  return `${part[0]?.toUpperCase() ?? ''}${part.slice(1)}`;
}

function administrationRoutes(code, stem) {
  const cmtBlock = code.match(
    /^\s*(?:\$CMT|\[CMT\])[^\r\n]*\r?\n([\s\S]*?)(?=^\s*(?:\$[A-Z]|\[[A-Z]+\])|(?![\s\S]))/m
  )?.[1] ?? '';
  const compartments = cmtBlock
    .split(/\r?\n/)
    .map((line) => line.match(/^\s*([A-Za-z][A-Za-z0-9_]*)\s*:\s*(.*?)\s*\[([^\]]+)\]\s*$/))
    .filter(Boolean)
    .filter((match) => match[3].split(',').some((tag) => tag.trim() === 'ADM'))
    .map((match) => ({ name: match[1], description: match[2] }));

  if (!compartments.length) throw new Error(`No annotated [ADM] compartment found for ${stem}.`);

  const routeCompartments = {};
  for (const compartment of compartments) {
    const descriptor = `${compartment.name} ${compartment.description}`;
    const route = /oral|gut|depot|transit|absorp/i.test(descriptor) ? 'Oral' : 'IV';
    routeCompartments[route] ??= compartment.name;
  }

  const routes = ['IV', 'Oral'].filter((route) => routeCompartments[route]);
  return {
    routes,
    ivCmt: routeCompartments.IV ?? null,
    oralCmt: routeCompartments.Oral ?? null
  };
}

function administrationModes(routes, details, stem) {
  const defaults = routes.flatMap((route) => route === 'Oral' ? ['ORAL'] : ['IV_INTERMITTENT']);
  const modes = [...new Set(details.administrationModes ?? defaults)];
  const unknown = modes.filter((mode) => !administrationModeLabels[mode]);
  if (unknown.length) throw new Error(`Unknown administration mode for ${stem}: ${unknown.join(', ')}`);
  const modeRoutes = modes.map((mode) => mode === 'ORAL' ? 'Oral' : 'IV');
  if (modeRoutes.some((route) => !routes.includes(route)) || routes.some((route) => !modeRoutes.includes(route))) {
    throw new Error(`Administration modes do not cover the declared routes for ${stem}.`);
  }
  return {
    administrationModes: modes,
    administrationCategories: modes.map((mode) => administrationModeLabels[mode].fr),
    administrationCategoriesEn: modes.map((mode) => administrationModeLabels[mode].en)
  };
}

function implementationStatus(details, stem) {
  const adapted = additionallyAdaptedModels.has(stem) || /adapté|module/i.test(details.modelType);
  return adapted
    ? { implementationStatus: 'ADAPTED', implementationStatusLabel: 'Adaptation documentée', implementationStatusLabelEn: 'Documented adaptation' }
    : { implementationStatus: 'ARTICLE', implementationStatusLabel: "Implémentation de l'article", implementationStatusLabelEn: 'Article implementation' };
}

function maturityStatus(details, stem) {
  for (const field of ['verificationStatus', 'validationStatus', 'contextOfUse']) {
    if (!Object.hasOwn(details, field)) throw new Error(`Missing explicit ${field} for ${stem}.`);
  }
  const { verificationStatus, validationStatus, contextOfUse } = details;
  if (!['NOT_TESTED', 'SOFTWARE_TESTED', 'NUMERICALLY_REPRODUCED'].includes(verificationStatus)) {
    throw new Error(`Unknown verification status for ${stem}: ${verificationStatus}`);
  }
  if (!['NONE', 'INTERNAL', 'EXTERNAL'].includes(validationStatus)) {
    throw new Error(`Unknown validation status for ${stem}: ${validationStatus}`);
  }
  if (!['TEACHING_RESEARCH', 'EXPLORATORY_MIPD', 'CLINICALLY_VALIDATED'].includes(contextOfUse)) {
    throw new Error(`Unknown context of use for ${stem}: ${contextOfUse}`);
  }
  return { verificationStatus, validationStatus, contextOfUse };
}

function analysisEligibility(details, implementationDetails, stem) {
  if (!Object.hasOwn(details, 'analysisEligible') || typeof details.analysisEligible !== 'boolean') {
    throw new Error(`Missing explicit boolean analysisEligible for ${stem}.`);
  }
  const analysisEligible = details.analysisEligible;
  return analysisEligible
    ? {
        analysisEligible: true,
        analysisEligibilityLabel: "Prior MAP fondé sur l'article",
        analysisEligibilityLabelEn: 'Article-based MAP prior'
      }
    : {
        analysisEligible: false,
        analysisEligibilityLabel: 'Simulation uniquement',
        analysisEligibilityLabelEn: 'Simulation only'
      };
}

function parseModelFile(file, metadata, englishMetadata, code) {
  const stem = file.replace(/\.cpp$/i, '');
  const drugKey = drugKeys.find((key) => stem.startsWith(`${key}_`)) ?? stem.split('_')[0];
  const suffix = stem.slice(drugKey.length + 1);
  const model = suffix.split('_').filter((part) => part && part.toLowerCase() !== 'ddi').map(titleCase).join(' ');
  const drug = drugLabels[drugKey] ?? titleCase(drugKey);

  const details = metadata[stem];
  const english = englishMetadata[stem];
  if (!details) throw new Error(`Missing TDM metadata for ${stem}.`);
  if (!english?.population) throw new Error(`Missing English TDM population for ${stem}.`);
  if ('provenance' in details) throw new Error(`Application provenance is not an article source: ${stem}.`);
  if (!details.citation || !details.doi) {
    throw new Error(`Missing primary article citation or DOI for ${stem}.`);
  }
  if (!/^10\.\d{4,9}\//i.test(details.doi)) throw new Error(`Invalid DOI for ${stem}: ${details.doi}`);

  const displayModel = details.model ?? model;
  if (!details.citation.toLocaleLowerCase('fr').startsWith(displayModel.toLocaleLowerCase('fr'))) {
    throw new Error(`Displayed model name must match the citation first author for ${stem}.`);
  }
  if (details.sourceStatus === 'secondary') {
    throw new Error(`Secondary references cannot be published in the TDM catalog: ${stem}.`);
  }
  const routeDetails = administrationRoutes(code, stem);
  const modeDetails = administrationModes(routeDetails.routes, details, stem);
  const implementationDetails = implementationStatus(details, stem);
  if (implementationDetails.implementationStatus === 'ADAPTED' && !details.note?.trim()) {
    throw new Error(`Adapted model ${stem} must document its implementation difference in note.`);
  }
  const maturityDetails = maturityStatus(details, stem);
  const analysisDetails = analysisEligibility(details, implementationDetails, stem);

  return {
    id: stem,
    file,
    drugKey,
    timeUnit: 'h',
    concentrationUnit: ['ciclosporine', 'everolimus', 'sirolimus', 'tacrolimus'].includes(drugKey) ? 'ng/mL' : 'mg/L',
    drug,
    drugEn: drugLabelsEn[drugKey] ?? drug,
    model: displayModel,
    format: 'mrgsolve/C++',
    href: `/tdm/models/${file}`,
    ...details,
    populationEn: english.population,
    populationTagsEn: (details.populationTags ?? []).map((tag) => populationTagEn[tag] ?? tag),
    modelTypeEn: modelTypeEn[details.modelType] ?? details.modelType,
    noteEn: english.note ?? '',
    ...routeDetails,
    ...modeDetails,
    ...implementationDetails,
    ...maturityDetails,
    ...analysisDetails,
    tags: Array.from(new Set([
      drug,
      model,
      details.modelType,
      details.citation,
      details.doi,
      english.population,
      drugLabelsEn[drugKey],
      ...routeDetails.routes,
      ...modeDetails.administrationModes,
      ...modeDetails.administrationCategories,
      ...modeDetails.administrationCategoriesEn,
      implementationDetails.implementationStatusLabel,
      implementationDetails.implementationStatusLabelEn,
      analysisDetails.analysisEligibilityLabel,
      analysisDetails.analysisEligibilityLabelEn,
      ...(details.populationTags ?? []),
      ...(details.populationTags ?? []).map((tag) => populationTagEn[tag] ?? tag)
    ].filter(Boolean)))
  };
}

async function main() {
  const entries = await fs.readdir(modelsDirectory, { withFileTypes: true });
  const allFiles = entries
    .filter((entry) => entry.isFile() && entry.name.toLowerCase().endsWith('.cpp'))
    .map((entry) => entry.name)
    .sort((a, b) => a.localeCompare(b, 'en'));

  if (!allFiles.length) throw new Error(`No .cpp model found in ${modelsDirectory}`);

  for (const file of allFiles) {
    const filePath = path.join(modelsDirectory, file);
    const code = await fs.readFile(filePath, 'utf8');
    const normalized = code.replace(/[ \t]+$/gm, '');
    if (normalized !== code) await fs.writeFile(filePath, normalized, 'utf8');
  }

  const metadataDocument = JSON.parse(await fs.readFile(metadataFile, 'utf8'));
  const metadata = metadataDocument.models ?? {};
  const englishMetadataDocument = JSON.parse(await fs.readFile(englishMetadataFile, 'utf8'));
  const englishMetadata = englishMetadataDocument.models ?? {};
  const allIds = new Set(allFiles.map((file) => file.replace(/\.cpp$/i, '')));
  const missingMetadata = [...allIds].filter((id) => !metadata[id]);
  if (missingMetadata.length) throw new Error(`Missing TDM metadata for: ${missingMetadata.join(', ')}`);

  const files = allFiles.filter((file) => metadata[file.replace(/\.cpp$/i, '')].listed !== false);
  const modelCode = new Map(await Promise.all(files.map(async (file) => [
    file,
    await fs.readFile(path.join(modelsDirectory, file), 'utf8')
  ])));
  const models = files.map((file) => parseModelFile(file, metadata, englishMetadata, modelCode.get(file)));
  const ids = new Set(models.map((model) => model.id));
  if (ids.size !== models.length) throw new Error('Duplicate TDM model identifiers detected.');
  const extraMetadata = Object.keys(metadata).filter((id) => !allIds.has(id));
  if (extraMetadata.length) throw new Error(`Metadata without a matching model: ${extraMetadata.join(', ')}`);

  const output = `${JSON.stringify({ version: 5, models }, null, 2)}\n`;
  const previous = await fs.readFile(outputFile, 'utf8').catch(() => '');
  if (previous !== output) await fs.writeFile(outputFile, output, 'utf8');


  const mlRegistry = JSON.parse(await fs.readFile(mlRegistryFile, 'utf8'));
  const modelById = new Map(models.map((model) => [model.id, model]));
  const modelHashes = new Map([...modelCode].map(([file, code]) => [
    file.replace(/\.cpp$/i, ''),
    createHash('sha256').update(code.replace(/\r/g, '')).digest('hex')
  ]));
  const registeredArtifacts = mlRegistry.artifacts ?? [];
  const requiredArtifactPaths = [
    'id', 'baseModelId', 'baseModelSha256', 'artifactPath', 'artifactSha256',
    'route', 'administrationMode', 'samplingProtocol.description', 'samplingProtocol.target',
    'featureSchema', 'prediction.metric', 'prediction.unit', 'prediction.horizonHours', 'methodology.sourcePopulation',
    'methodology.regimenConfiguration', 'methodology.populationConfiguration',
    'training.seed', 'training.nVirtualDevelopment', 'training.nTraining', 'training.nUntouchedHoldout',
    'training.split.unit', 'training.split.method', 'training.split.holdoutUntouched',
    'training.software', 'trainingDomain.features', 'validation'
  ];
  const valueAt = (object, dottedPath) => dottedPath.split('.').reduce((value, key) => value?.[key], object);
  for (const artifact of registeredArtifacts) {
    const missing = requiredArtifactPaths.filter((field) => {
      const value = valueAt(artifact, field);
      return value == null || value === '' || (Array.isArray(value) && value.length === 0) ||
        (typeof value === 'object' && !Array.isArray(value) && Object.keys(value).length === 0);
    });
    if (missing.length) throw new Error(`Incomplete ML artifact ${artifact.id ?? '<unknown>'}: ${missing.join(', ')}`);
  }
  const mlArtifacts = registeredArtifacts.map((artifact) => {
    const model = modelById.get(artifact.baseModelId);
    const currentHash = modelHashes.get(artifact.baseModelId);
    const hashMatches = Boolean(currentHash && currentHash === artifact.baseModelSha256);
    return {
      id: artifact.id,
      modelId: artifact.baseModelId,
      model: model?.model ?? artifact.baseModelId,
      drug: artifact.drug,
      route: artifact.route,
      administrationMode: artifact.administrationMode,
      releaseLevel: artifact.releaseLevel,
      analysisEligible: Boolean(model?.analysisEligible),
      hashMatches,
      availableInTdm: Boolean(model?.analysisEligible && hashMatches),
      minimumObservations: artifact.samplingProtocol?.minimumObservations ?? null,
      steadyStateRequired: Boolean(artifact.samplingProtocol?.steadyStateRequired),
      nTraining: artifact.training?.nTraining ?? null,
      nHoldout: artifact.training?.nUntouchedHoldout ?? null,
      nAlternate: artifact.training?.nUnseenPopPk ?? null,
      validation: artifact.validation ?? {}
    };
  });
  const pairedFile = path.join(root, 'ml', 'validation', 'paired-auc-benchmark.json');
  const paired = JSON.parse(await fs.readFile(pairedFile, 'utf8').catch((error) => {
    if (error.code !== 'ENOENT') throw error;
    return '{"results":[]}';
  }));
  if (paired.smoke) throw new Error('A smoke benchmark must not be published.');
  const pairedResults = (paired.results ?? []).map((result) => ({
    ...result,
    model: modelById.get(result.modelId)?.model ?? result.modelId,
    analysisEligible: Boolean(modelById.get(result.modelId)?.analysisEligible),
    hashMatches: result.baseModelSha256 === modelHashes.get(result.modelId)
  }));
  const mlOutput = `${JSON.stringify({
    version: mlRegistry.version,
    benchmarkDate: mlRegistry.benchmarkDate ?? null,
    samplingDesign: mlRegistry.samplingDesign ?? null,
    artifacts: mlArtifacts,
    pairedBenchmark: { ...paired, results: pairedResults }
  }, null, 2)}\n`;
  const previousMlBenchmark = await fs.readFile(mlBenchmarkFile, 'utf8').catch(() => '');
  if (previousMlBenchmark !== mlOutput) await fs.writeFile(mlBenchmarkFile, mlOutput, 'utf8');

  console.log(`TDM catalog OK: ${models.length} models and ${mlArtifacts.length} ML artifacts indexed.`);
}

await main();
