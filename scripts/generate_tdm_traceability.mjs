import { promises as fs } from 'node:fs';
import { createHash } from 'node:crypto';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const catalogFile = path.join(root, 'models', 'catalog.json');
const sourceReviewFile = path.join(root, 'models', 'source-review.json');
const pdfInventoryFile = path.join(root, 'models', 'pdf-inventory.json');
const modelsDirectory = path.join(root, 'models');
const outputFile = path.join(root, 'docs', 'tdm-model-traceability.md');
const manifestFile = path.join(root, 'docs', 'tdm-model-traceability.json');
const checkOnly = process.argv.includes('--check');

const referenceCases = new Map([
  ['amox_carlier', '1 000 mg IV; concentrations at 1, 4, 12 and 24 h'],
  ['vanco_pkjust', '1 000 mg IV; concentrations at 1, 4, 12 and 24 h'],
  ['vanco_goti', '1 000 mg IV; concentrations at 1, 4, 12 and 24 h'],
  ['vanco_roberts', '1 000 mg/24 h by continuous IV infusion; concentrations at 1, 4, 12 and 24 h'],
  ['tacrolimus_woillard_ddi', '5 mg oral; concentrations at 1, 4, 12 and 24 h'],
  ['dapto_dvorchik_IV_adults', '500 mg IV; concentrations at 1, 4, 12 and 24 h'],
  ['cefepime_jonckheere', '2 000 mg by 30-minute IV infusion; concentrations at 1, 4, 12 and 24 h'],
  ['rifampicine_jing', '450 mg oral; concentrations at 1, 4, 12 and 24 h']
]);

function sections(code) {
  const result = [];
  let current = null;
  for (const line of code.replace(/\r/g, '').split('\n')) {
    const header = line.match(/^\s*(?:\$([A-Z]+)\b(.*)|\[([A-Z]+)\](.*))$/);
    if (header) {
      current = { name: header[1] ?? header[3], options: (header[2] ?? header[4]).trim(), lines: [] };
      result.push(current);
    } else if (current) {
      current.lines.push(line);
    }
  }
  return result;
}

function declaredNames(blocks) {
  return [...new Set(blocks.flatMap((block) => block.lines).map((line) => {
    const clean = line.replace(/\/\/.*$/, '').trim();
    return clean.match(/^([A-Za-z][A-Za-z0-9_]*)\s*(?::|=)/)?.[1] ?? null;
  }).filter(Boolean))];
}

function indexedUses(code, token) {
  return [...new Set([...code.matchAll(new RegExp(`\\b${token}\\(\\s*(\\d+)\\s*\\)`, 'g'))]
    .map((match) => Number(match[1])))]
    .sort((a, b) => a - b);
}

function compactList(values) {
  return values.length ? values.map((value) => `\`${value}\``).join(', ') : 'None declared';
}

function escapeCell(value) {
  return String(value ?? '').replace(/\|/g, '\\|').replace(/\s+/g, ' ').trim();
}

function articleStatus(model, review) {
  // The implementation family and its independently reviewed source fidelity
  // are different axes. A historical ARTICLE implementation may be downgraded
  // to unverified without rewriting the frozen C++ or changing its family.
  if (review.result === 'SOURCE_AMBIGUITY_DOCUMENTED' &&
      model.implementationStatus === 'ARTICLE') {
    return 'Article-derived implementation; source fidelity NOT VERIFIED for identified boundary branches.';
  }
  return model.implementationStatus === 'ARTICLE'
    ? 'Article implementation according to the author audit; independent reproduction pending.'
    : 'Documented adaptation; see the implementation note below.';
}

const catalog = JSON.parse(await fs.readFile(catalogFile, 'utf8'));
const sourceReview = JSON.parse(await fs.readFile(sourceReviewFile, 'utf8'));
const pdfInventory = JSON.parse(await fs.readFile(pdfInventoryFile, 'utf8'));
const reviewedSources = new Map(sourceReview.reviewed.map((review) => [review.id, review]));
const reviewIds = [...reviewedSources.keys()];
const catalogIds = catalog.models.map((model) => model.id);
if (new Set(reviewIds).size !== reviewIds.length ||
    reviewIds.some((id) => !catalogIds.includes(id)) ||
    catalogIds.some((id) => !reviewIds.includes(id))) {
  throw new Error('Source-review coverage must contain every catalog model exactly once.');
}
for (const model of catalog.models) {
  const review = reviewedSources.get(model.id);
  const expectedStatus = review.result === 'ARTICLE_ALIGNED' ? 'ARTICLE' : 'ADAPTED';
  const conservativeSourceDowngrade =
    model.implementationStatus === 'ARTICLE' &&
    review.result === 'SOURCE_AMBIGUITY_DOCUMENTED';
  if (model.implementationStatus !== expectedStatus && !conservativeSourceDowngrade) {
    throw new Error(`Catalog/source-review status mismatch for ${model.id}: expected ${expectedStatus}.`);
  }
}
const pdfClassifications = {
  ACTIVE_MODEL_PRIMARY_SOURCE: 'activeModelPrimarySources',
  DUPLICATE_PRIMARY_SOURCE: 'duplicatePrimarySources',
  SECONDARY_REVIEW: 'secondaryReviews',
  APPLICATION_CONTEXT: 'applicationContext',
  BACKGROUND_OR_NONCATALOG_SOURCE: 'backgroundOrNoncatalogSources'
};
const countedPdfs = Object.fromEntries(Object.values(pdfClassifications).map((key) => [key, 0]));
for (const document of pdfInventory.documents) {
  const countKey = pdfClassifications[document.classification];
  if (!countKey) throw new Error(`Unknown PDF classification: ${document.classification}`);
  countedPdfs[countKey] += 1;
  const unknownModels = document.models.filter((id) => !catalogIds.includes(id));
  if (unknownModels.length) throw new Error(`PDF inventory references unknown models: ${unknownModels.join(', ')}`);
}
if (pdfInventory.documents.length !== pdfInventory.counts.documents ||
    Object.entries(countedPdfs).some(([key, value]) => pdfInventory.counts[key] !== value)) {
  throw new Error('PDF inventory counts do not match its documents.');
}
const records = await Promise.all(catalog.models.map(async (model) => {
  const code = await fs.readFile(path.join(modelsDirectory, model.file), 'utf8');
  const parsed = sections(code);
  const parameterBlocks = parsed.filter((block) => block.name === 'PARAM');
  const covariateBlocks = parameterBlocks.filter((block) => /@covariates\b/i.test(block.options));
  const ordinaryBlocks = parameterBlocks.filter((block) => !/@covariates\b/i.test(block.options));
  const covariates = declaredNames(covariateBlocks);
  const parameters = declaredNames(ordinaryBlocks).filter((name) => !/^ETA\d+$/i.test(name));
  const etaUses = indexedUses(code, 'ETA');
  const epsUses = indexedUses(code, 'EPS');
  const blockText = (name) => parsed.filter((block) => block.name === name)
    .map((block) => [block.options, ...block.lines].join('\n').trim()).filter(Boolean);
  const compartments = declaredNames(parsed.filter((block) => block.name === 'CMT'));
  const differentialEquations = parsed.filter((block) => block.name === 'ODE')
    .flatMap((block) => block.lines.map((line) => line.replace(/\/\/.*$/, '').trim()))
    .filter((line) => /^dxdt_[A-Za-z][A-Za-z0-9_]*\s*=/.test(line));
  const outputs = parsed.filter((block) => block.name === 'CAPTURE')
    .flatMap((block) => [block.options, ...block.lines].join(' ').split(/\s+/).filter(Boolean));
  const sourceSha256 = createHash('sha256').update(code.replace(/\r/g, '')).digest('hex');
  return { model, parameters, covariates, etaUses, epsUses, compartments, differentialEquations, outputs,
    omega: blockText('OMEGA'), sigma: blockText('SIGMA'), sourceSha256 };
}));

const articleCount = records.filter(({ model }) => model.implementationStatus === 'ARTICLE').length;
const adaptedCount = records.length - articleCount;
const mapCount = records.filter(({ model }) => model.analysisEligible).length;
const fullTextCount = [...reviewedSources.values()].filter((review) => review.fullTextReviewed).length;
const limitedEvidenceCount = reviewedSources.size - fullTextCount;

const lines = [
  '# TDM model source-to-implementation traceability',
  '',
  '> Generated by `scripts/generate_tdm_traceability.mjs`; do not edit manually.',
  '',
  'This inventory links each shipped model to its declared primary source and to the implementation that is actually executed. It is an implementation audit, **not** proof that the publication has been independently reproduced. Three levels are kept separate:',
  '',
  '1. **Syntax and contract:** the complete library compiles and satisfies the mrgsolve/TDM interface contract.',
  '2. **Deterministic software regression:** selected models have fixed inputs and expected outputs in `tests/model_reference_cases.R`.',
  '3. **Scientific reproduction:** comparison with source equations, source software or published reference outputs; this remains pending unless explicitly documented.',
  '',
  `Current catalog: **${records.length} models** (${articleCount} labelled article implementations, ${adaptedCount} documented adaptations; ${mapCount} exposed as MAP priors). **${referenceCases.size} models** currently have deterministic numeric reference profiles.`,
  `Every model has a recorded source comparison. Primary full text has been checked for **${fullTextCount}/${records.length} models**; **${limitedEvidenceCount}** older or paywalled sources currently have an abstract-level comparison, explicitly identified below. This author-led source audit is not independent reproduction.`,
  '',
  'Population text and administration route define the currently encoded validity context. Numerical age, weight, renal-function and dose ranges are not yet structured for every article; absence of a range below must not be interpreted as permission to extrapolate.',
  '',
  '## Summary',
  '',
  '| Model | Source | Full-text review | Route / mode | Units | Implementation | MAP | Numeric profile |',
  '| --- | --- | --- | --- | --- | --- | --- | --- |'
];

for (const { model } of records) {
  const source = `[${escapeCell(model.citation)}](https://doi.org/${model.doi})`;
  const route = `${model.routes.join(' + ')} / ${model.administrationModes.join(', ')}`;
  const sourceFinding = reviewedSources.get(model.id);
  const sourceStatus = `${sourceFinding.result} / ${sourceFinding.evidence}`;
  lines.push(`| \`${model.id}\` | ${source} | ${sourceStatus} | ${escapeCell(route)} | ${model.timeUnit}; ${model.concentrationUnit} | ${model.implementationStatus} | ${model.analysisEligible ? 'Yes' : 'No'} | ${referenceCases.has(model.id) ? 'Deterministic regression' : 'Compile + contract'} |`);
}

lines.push('', '## Per-model implementation record', '');

for (const { model, parameters, covariates, etaUses, epsUses } of records) {
  const note = model.note?.trim() || 'No adaptation is documented in catalog metadata. Independent source reproduction is still required.';
  const reference = referenceCases.get(model.id);
  const sourceFinding = reviewedSources.get(model.id);
  lines.push(
    `### ${model.drug} - ${model.model} (\`${model.id}\`)`,
    '',
    `- **Primary source:** ${model.citation} DOI: [${model.doi}](https://doi.org/${model.doi}).`,
    `- **Source population:** ${model.population}`,
    `- **Administration and units:** ${model.routes.join(' + ')}; ${model.administrationModes.join(', ')}; time in ${model.timeUnit}; concentration in ${model.concentrationUnit}.`,
    `- **Structural/fixed parameters declared in code:** ${compactList(parameters)}.`,
    `- **Covariates declared in code:** ${compactList(covariates)}.`,
    `- **Random and residual effects used by the equations:** ${etaUses.length ? `ETA indices ${etaUses.join(', ')}` : 'no indexed ETA use'}; ${epsUses.length ? `EPS indices ${epsUses.join(', ')}` : 'no indexed EPS use'}. OMEGA/SIGMA details remain authoritative in [the model source](../models/${model.file}).`,
    `- **MAP prior:** ${model.analysisEligible ? 'exposed in the TDM analysis workflow' : 'not exposed; simulation only'}.`,
    `- **Fidelity status:** ${articleStatus(model, sourceFinding)}`,
    `- **Primary-source review:** ${sourceFinding.result} (${sourceFinding.evidence}, ${sourceReview.reviewDate}). ${sourceFinding.note}`,
    `- **Documented modification or limitation:** ${note}`,
    `- **Reference testing:** ${reference ? `${reference}; zero ETA and residual errors; deterministic software-regression tolerance \`1e-7\`.` : 'compilation and interface-contract test only; no frozen numeric profile yet.'}`,
    ''
  );
}

const output = `${lines.join('\n')}\n`;
const previous = await fs.readFile(outputFile, 'utf8').catch(() => '');
const manifest = {
  version: 1,
  generatedFrom: ['models/catalog.json', 'models/source-review.json', 'models/pdf-inventory.json', 'models/*.cpp'],
  models: records.map(({ model, parameters, covariates, etaUses, epsUses, compartments, differentialEquations, outputs, omega, sigma, sourceSha256 }) => ({
    id: model.id,
    sourceFile: model.file,
    sourceSha256,
    primarySource: {
      citation: model.citation,
      doi: model.doi,
      status: model.sourceStatus ?? 'primary',
      fullTextReview: {
        status: reviewedSources.get(model.id).fullTextReviewed ? 'full-text-reviewed' : 'limited-evidence-reviewed',
        evidence: reviewedSources.get(model.id).evidence,
        date: sourceReview.reviewDate,
        result: reviewedSources.get(model.id).result,
        note: reviewedSources.get(model.id).note
      }
    },
    population: { description: model.population, tags: model.populationTags ?? [] },
    administration: { routes: model.routes, modes: model.administrationModes, timeUnit: model.timeUnit, concentrationUnit: model.concentrationUnit },
    structuralModel: { compartments, differentialEquations, outputs },
    parameters,
    covariates,
    interindividualVariability: { etaIndices: etaUses, omega },
    residualError: { epsIndices: epsUses, sigma },
    pmxAdaptations: model.note?.trim() || 'No PMx adaptation documented.',
    implementationStatus: model.implementationStatus,
    softwareContractStatus: model.verificationStatus,
    numericalCheckStatus: referenceCases.has(model.id) ? 'deterministic-software-regression' : 'not-checked',
    sourceReproductionStatus: model.verificationStatus === 'NUMERICALLY_REPRODUCED' ? 'declared-numerically-reproduced' : 'pending-independent-reproduction',
    mapTdm: { eligible: Boolean(model.analysisEligible), status: model.analysisEligibilityLabel },
    limitations: model.note?.trim() || 'Numerical validity ranges and independent reproduction remain to be documented.'
  }))
};
const requiredPaths = ['primarySource.doi', 'population.description', 'administration.routes', 'structuralModel.compartments',
  'parameters', 'covariates', 'interindividualVariability.omega', 'residualError.sigma', 'pmxAdaptations',
  'implementationStatus', 'softwareContractStatus', 'numericalCheckStatus', 'sourceReproductionStatus', 'mapTdm.status', 'limitations'];
const valueAt = (object, dottedPath) => dottedPath.split('.').reduce((value, key) => value?.[key], object);
for (const record of manifest.models) {
  const missing = requiredPaths.filter((field) => valueAt(record, field) == null);
  if (missing.length) throw new Error(`Incomplete TDM traceability for ${record.id}: ${missing.join(', ')}`);
}
const manifestOutput = `${JSON.stringify(manifest, null, 2)}\n`;
const previousManifest = await fs.readFile(manifestFile, 'utf8').catch(() => '');

if (checkOnly) {
  if (previous !== output || previousManifest !== manifestOutput) {
    console.error('TDM traceability inventory is stale. Run `npm run tdm:traceability`.');
    process.exit(1);
  }
  console.log(`TDM traceability inventory is current: ${records.length} models.`);
} else {
  if (previous !== output) await fs.writeFile(outputFile, output, 'utf8');
  if (previousManifest !== manifestOutput) await fs.writeFile(manifestFile, manifestOutput, 'utf8');
  console.log(`Generated ${path.relative(root, outputFile)} for ${records.length} models.`);
}
