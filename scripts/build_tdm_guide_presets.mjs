// Builds the versioned clinical-preset artifact the Shiny app reads.
// Targets come from the PMx TDM Guide snapshot; dose grids are declared by MIPD,
// because the guide publishes no posology at all. Nothing here adjusts a dose.
import { readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { checkReferenceSnapshot, validateCrosswalk, checkApplicabilityRules,
  assessCandidateApplicability } from './preview_tdm_guide_recommendations.mjs';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
export const PRESET_SCHEMA = 'mipd.tdm.clinical-presets.v1';
const sha256 = (value) => createHash('sha256').update(JSON.stringify(value)).digest('hex');

// A grid may only ever be written into the simulation fields once it is signed off.
export const DOSE_GRID_SIGNOFFS = Object.freeze(['PENDING', 'APPROVED']);
export const doseGridIsApplicable = (grid) => grid?.signoff === 'APPROVED';

export function checkDoseGrids(grids) {
  if (grids?.schema !== 'mipd.tdm.dose-grids.v1' || !grids.grids ||
    typeof grids.grids !== 'object') return false;
  if (!Array.isArray(grids.applicationPolicy?.allowedSignoffValues) ||
    grids.applicationPolicy.appliedWhenSignoffIs !== 'APPROVED') return false;
  return Object.values(grids.grids).every((grid) =>
    Number.isFinite(grid.dose_min) && Number.isFinite(grid.dose_max) &&
    Number.isFinite(grid.dose_step) && grid.dose_min > 0 &&
    grid.dose_max > grid.dose_min && grid.dose_step > 0 &&
    Array.isArray(grid.intervals) && grid.intervals.length &&
    grid.intervals.every((h) => Number.isFinite(h) && h > 0) &&
    Number.isFinite(grid.infusion) && grid.infusion >= 0 &&
    // A grid is a posology statement: it must quote the label it comes from.
    typeof grid.source?.url === 'string' && grid.source.url.startsWith('https://') &&
    typeof grid.source?.statement === 'string' && grid.source.statement.length > 0 &&
    DOSE_GRID_SIGNOFFS.includes(grid.signoff));
}

/**
 * A guide reference may feed a preset only when every condition below holds.
 * Any failure leaves it as documentation in the STP card, never as a preset.
 */
export function presetEligibility({ candidate, rules }) {
  const equivalence = rules.appMetricEquivalence.equivalent;
  const reasons = [];
  const metric = Object.hasOwn(equivalence, candidate.metric)
    ? equivalence[candidate.metric] : null;
  if (!metric) reasons.push('metric_not_equivalent_to_an_app_metric');
  if (candidate.relation !== 'range') reasons.push('not_an_explicit_two_sided_range');
  if (!Number.isFinite(candidate.min) || !Number.isFinite(candidate.max) ||
    candidate.max <= candidate.min) reasons.push('range_bounds_unusable');
  if (candidate.applicability?.mismatches?.length) reasons.push('structural_context_mismatch');
  return { eligible: reasons.length === 0, metric, reasons };
}

export function buildClinicalPresets({ catalogue, snapshot, crosswalk, rules, doseGrids }) {
  if (!Array.isArray(catalogue?.models) || !catalogue.models.length)
    throw new Error('Missing independent model catalog');
  if (!checkReferenceSnapshot(snapshot)) throw new Error('Invalid PMx snapshot or missing immutable SHA256');
  if (!validateCrosswalk(crosswalk, catalogue)) throw new Error('Model-drug crosswalk incomplete or stale');
  if (!checkApplicabilityRules(rules)) throw new Error('Applicability rules missing or stale');
  if (!checkDoseGrids(doseGrids)) throw new Error('MIPD dose grids missing, malformed or already signed off');

  const guide = new Map(snapshot.drugs.map((drug) => [drug.id, drug]));
  const presets = {};
  for (const model of catalogue.models) {
    const tdmDrugId = crosswalk.keys[model.drugKey];
    const entry = guide.get(tdmDrugId);
    if (!entry) continue;
    for (const mode of Array.isArray(model.administrationModes) ? model.administrationModes : []) {
      const key = `${model.drugKey}::${mode}`;
      const declaredGrid = Object.hasOwn(doseGrids.grids, key) ? doseGrids.grids[key] : null;
      if (!presets[key]) presets[key] = {
        drugKey: model.drugKey, mode, tdmDrugId,
        doseGrid: declaredGrid,
        // Stated explicitly so the panel never has to infer it from the signoff string.
        doseGridApplicable: doseGridIsApplicable(declaredGrid),
        candidates: []
      };
      const bucket = presets[key];
      for (const target of entry.targets) {
        // Gate against THIS model, so a context mismatch on any model keeps the
        // reference out of the shared preset list for that molecule and route.
        const applicability = assessCandidateApplicability({ target, model, tdmDrugId, rules });
        const { eligible, metric } = presetEligibility(
          { candidate: { ...target, applicability }, rules });
        if (!eligible) continue;
        if (bucket.candidates.some((existing) => existing.id === target.id)) continue;
        bucket.candidates.push({
          id: target.id,
          target_metric: metric,
          target_low: target.min,
          target_high: target.max,
          target_unit: target.unit,
          indication: target.indication,
          labelFr: target.indicationFr, labelEn: target.indicationEn,
          scopeFr: target.contextFr, scopeEn: target.contextEn,
          citation: target.source?.citation ?? null,
          url: target.source?.url ?? null,
          evidenceStatus: target.status ?? null,
          // Carried through so the panel can state what the guide does not publish.
          undecidable: applicability.undecidable,
          modelApplicabilityVerified: false,
          automaticDoseAdjustmentAllowed: false
        });
      }
    }
  }
  // Only a single qualified reference may ever be written without the reader choosing.
  for (const bucket of Object.values(presets)) {
    bucket.candidates.sort((a, b) => a.id.localeCompare(b.id));
    bucket.autoApply = bucket.candidates.length === 1;
  }
  for (const key of Object.keys(presets))
    if (!presets[key].candidates.length) delete presets[key];

  const body = {
    schema: PRESET_SCHEMA,
    producer: 'MIPD-Engine standalone preset builder',
    sourceSnapshotSha256: snapshot.sha256,
    doseGridSignoff: doseGrids.signoff,
    intendedUse: 'research_education_target_exploration',
    clinicalInterpretationAllowed: false,
    automaticDoseAdjustmentAllowed: false,
    patientDataIncluded: false,
    presets
  };
  return { ...body, sha256: sha256(body) };
}

export async function loadClinicalPresets(snapshotFile) {
  if (!snapshotFile || !snapshotFile.endsWith('.json'))
    throw new Error('Explicit local JSON snapshot path required');
  const [cat, map, gates, grids, snap] = await Promise.all([
    readFile(path.join(root, 'models/catalog.json'), 'utf8'),
    readFile(path.join(root, 'integration/tdm-guide-drug-crosswalk.preview.v1.json'), 'utf8'),
    readFile(path.join(root, 'integration/tdm-guide-applicability-rules.preview.v1.json'), 'utf8'),
    readFile(path.join(root, 'integration/tdm-guide-dose-grids.v1.json'), 'utf8'),
    readFile(snapshotFile, 'utf8')
  ]);
  return buildClinicalPresets({
    catalogue: JSON.parse(cat), crosswalk: JSON.parse(map), rules: JSON.parse(gates),
    doseGrids: JSON.parse(grids), snapshot: JSON.parse(snap)
  });
}

if (process.argv[1]?.endsWith('build_tdm_guide_presets.mjs')) {
  const at = process.argv.indexOf('--snapshot');
  if (at < 0 || !process.argv[at + 1])
    throw new Error('Usage: --snapshot /local/verified-guide.json [--output integration/tdm-guide-clinical-presets.v1.json] [--check]');
  const result = await loadClinicalPresets(process.argv[at + 1]);
  const out = process.argv.indexOf('--output');
  const target = out >= 0 ? process.argv[out + 1] : null;
  if (process.argv.includes('--check')) {
    if (!target) throw new Error('--check needs --output <committed artifact>');
    const committed = JSON.parse(await readFile(target, 'utf8'));
    // The declared digest alone would accept a body edited together with its own
    // sha256 field, so re-hash the committed body and compare the objects in full.
    const { sha256: declared, ...body } = committed;
    if (sha256(body) !== declared)
      throw new Error('Committed preset artifact does not match its own sha256: ' + target +
        ' was edited by hand');
    const { deepStrictEqual } = await import('node:assert');
    deepStrictEqual(committed, result,
      'Committed preset artifact is stale: regenerate ' + target);
    process.stdout.write('Clinical preset artifact is self-consistent and matches the pinned guide snapshot.\n');
  } else if (target) {
    if (!target.endsWith('.json')) throw new Error('Output filename must end with .json');
    await writeFile(target, JSON.stringify(result, null, 2) + '\n', 'utf8');
  } else {
    process.stdout.write(JSON.stringify(result, null, 2) + '\n');
  }
}
