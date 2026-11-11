// The clinical-preset artifact is a posology-adjacent surface: it is tested against
// the committed file and against the exclusions that keep unusable references out.
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { buildClinicalPresets, checkDoseGrids, presetEligibility, PRESET_SCHEMA,
  DOSE_GRID_SIGNOFFS, doseGridIsApplicable }
  from '../scripts/build_tdm_guide_presets.mjs';

const read = async (file) => JSON.parse(await readFile(new URL(file, import.meta.url), 'utf8'));
const cat = await read('../models/catalog.json');
const map = await read('../integration/tdm-guide-drug-crosswalk.preview.v1.json');
const rules = await read('../integration/tdm-guide-applicability-rules.preview.v1.json');
const grids = await read('../integration/tdm-guide-dose-grids.v1.json');
const committed = await read('../integration/tdm-guide-clinical-presets.v1.json');

// --- MIPD-owned dose grids: a posology statement, so it must stay sourced and unsigned ---
assert.equal(checkDoseGrids(grids), true);
assert.equal(grids.signoff, 'PENDING_PHARMACOLOGIST_REVIEW');
assert.deepEqual(Object.keys(grids.grids).sort(),
  ['ciclosporine::ORAL', 'everolimus::ORAL', 'sirolimus::ORAL', 'tacrolimus::ORAL']);
for (const [key, grid] of Object.entries(grids.grids)) {
  assert.equal(grid.signoff, 'PENDING', key + ' must not be pre-signed');
  assert.match(grid.source.url, /^https:\/\//, key);
  assert.ok(grid.source.statement.length > 20, key + ' must quote its label');
  assert.ok(grid.basisFr && grid.basisEn, key + ' must justify its bounds');
}
// A grid that loses its source or its bounds disables the whole layer.
assert.equal(checkDoseGrids({ ...grids,
  grids: { ...grids.grids, 'tacrolimus::ORAL': { ...grids.grids['tacrolimus::ORAL'], source: { url: 'http://x' } } } }), false);
assert.equal(checkDoseGrids({ ...grids,
  grids: { ...grids.grids, 'tacrolimus::ORAL': { ...grids.grids['tacrolimus::ORAL'], dose_max: 0.1 } } }), false);
assert.equal(checkDoseGrids({ ...grids, schema: 'wrong' }), false);
// An unknown sign-off value, or a policy that does not gate on APPROVED, is refused.
assert.equal(checkDoseGrids({ ...grids,
  grids: { ...grids.grids, 'tacrolimus::ORAL': { ...grids.grids['tacrolimus::ORAL'], signoff: 'SORT_OF' } } }), false);
assert.equal(checkDoseGrids({ ...grids, applicationPolicy: undefined }), false);
assert.equal(checkDoseGrids({ ...grids,
  applicationPolicy: { ...grids.applicationPolicy, appliedWhenSignoffIs: 'PENDING' } }), false);
assert.equal(grids.applicationPolicy.appliedWhenSignoffIs, 'APPROVED');
assert.deepEqual(DOSE_GRID_SIGNOFFS, ['PENDING', 'APPROVED']);
// Only a signed-off grid is applicable; an unsigned one is displayed, never written.
assert.equal(doseGridIsApplicable(grids.grids['tacrolimus::ORAL']), false);
assert.equal(doseGridIsApplicable({ ...grids.grids['tacrolimus::ORAL'], signoff: 'APPROVED' }), true);
assert.equal(doseGridIsApplicable(null), false);

// --- The two label constraints that bound these grids, verified against the SmPCs ---
// Rapamune forbids making up a dose from multiples of the 0.5 mg tablet, which is not
// fully bioequivalent to the other strengths on Cmax, so the step must not be 0.5 mg.
assert.equal(grids.grids['sirolimus::ORAL'].dose_step, 1);
assert.ok(/0\.5 mg tablets should not be used as a substitute/
  .test(grids.grids['sirolimus::ORAL'].source.statement));
// Prograf's 0.20-0.30 mg/kg/day is a twice-daily immediate-release regimen: a 24 h
// interval would silently cover Advagraf and Envarsus, which are dosed differently.
assert.deepEqual(grids.grids['tacrolimus::ORAL'].intervals, [12]);
assert.ok(grids.grids['tacrolimus::ORAL'].excludesFormulations.length > 0);
assert.ok(Array.isArray(grids.knownLimitations) && grids.knownLimitations.length >= 2);

// --- Eligibility: a reference the app cannot represent is excluded, never converted ---
const base = { relation: 'range', min: 5, max: 10, applicability: { mismatches: [] } };
assert.equal(presetEligibility({ candidate: { ...base, metric: 'C₀' }, rules }).eligible, true);
assert.equal(presetEligibility({ candidate: { ...base, metric: 'AUC₂₄' }, rules }).metric, 'AUC24');
// AUC over 12 h is not the AUC over 24 h the app optimises.
const auc12 = presetEligibility({ candidate: { ...base, metric: 'AUC₀–₁₂' }, rules });
assert.equal(auc12.eligible, false);
assert.ok(auc12.reasons.includes('metric_not_equivalent_to_an_app_metric'));
assert.equal(presetEligibility({ candidate: { ...base, metric: 'C₂h' }, rules }).eligible, false);
assert.ok(Object.hasOwn(rules.appMetricEquivalence.deliberatelyExcluded, 'AUC₀–₁₂'));
// One-sided and inverted ranges cannot become a two-bound target.
assert.equal(presetEligibility({ candidate: { ...base, metric: 'C₀', relation: 'above', max: null }, rules }).eligible, false);
assert.equal(presetEligibility({ candidate: { ...base, metric: 'C₀', max: 5 }, rules }).eligible, false);
// A reference whose declared context contradicts the model never becomes a preset.
const mismatched = presetEligibility({ candidate: { ...base, metric: 'C₀',
  applicability: { mismatches: [{ code: 'organ_indication_mismatch' }] } }, rules });
assert.equal(mismatched.eligible, false);
assert.ok(mismatched.reasons.includes('structural_context_mismatch'));
// Prototype keys must not reach the equivalence table.
assert.equal(presetEligibility({ candidate: { ...base, metric: 'constructor' }, rules }).eligible, false);

// --- Built against the real pinned snapshot, compared with the committed artifact ---
const snapshotPath = process.argv[2];
if (!snapshotPath) {
  console.log('TDM preset artifact: dose grids and eligibility verified (no snapshot given, committed file not rebuilt).');
  process.exit(0);
}
const snapshot = JSON.parse(await readFile(snapshotPath, 'utf8'));
const built = buildClinicalPresets({ catalogue: cat, crosswalk: map, rules, doseGrids: grids, snapshot });
assert.equal(built.schema, PRESET_SCHEMA);
assert.equal(built.sha256, committed.sha256, 'Committed preset artifact is stale: regenerate it');
assert.deepEqual(built, committed);
assert.equal(built.clinicalInterpretationAllowed, false);
assert.equal(built.automaticDoseAdjustmentAllowed, false);
assert.equal(built.doseGridSignoff, 'PENDING_PHARMACOLOGIST_REVIEW');

const keys = Object.keys(built.presets);
assert.deepEqual(keys.sort(), ['ciclosporine::ORAL', 'everolimus::ORAL', 'linez::IV_INTERMITTENT',
  'linez::ORAL', 'sirolimus::ORAL', 'tacrolimus::ORAL', 'vanco::IV_CONTINUOUS',
  'vanco::IV_INTERMITTENT', 'voriconazole::ORAL'].sort());
const all = Object.values(built.presets).flatMap((entry) => entry.candidates);
// 3 ciclosporine + 3 everolimus + 2 linezolid routes + 2 sirolimus + 3 tacrolimus
// + 2 vancomycin routes + 1 voriconazole
assert.equal(all.length, 16);
assert.ok(all.every((candidate) => candidate.modelApplicabilityVerified === false &&
  candidate.automaticDoseAdjustmentAllowed === false &&
  ['AUC24', 'Cmin', 'Cmax'].includes(candidate.target_metric) &&
  candidate.target_high > candidate.target_low &&
  /^https:\/\//.test(candidate.url)));

// Only a single qualified reference may be auto-applied; a choice never is.
for (const [key, entry] of Object.entries(built.presets))
  assert.equal(entry.autoApply, entry.candidates.length === 1, key);
assert.deepEqual(keys.filter((key) => built.presets[key].autoApply).sort(),
  ['linez::IV_INTERMITTENT', 'linez::ORAL', 'vanco::IV_CONTINUOUS', 'vanco::IV_INTERMITTENT',
    'voriconazole::ORAL']);
// Every auto-applicable key already has a curated MIPD preset, so this artifact
// introduces no new automatic writing. Keep it that way.
const curated = await readFile(new URL('../R/clinical_presets.R', import.meta.url), 'utf8');
for (const key of keys.filter((k) => built.presets[k].autoApply))
  assert.ok(curated.includes('"' + key + '"'),
    'New auto-applied key without a curated MIPD preset: ' + key);

// The immunosuppressants are exactly the multi-target keys, and they carry a grid.
for (const key of ['tacrolimus::ORAL', 'ciclosporine::ORAL', 'everolimus::ORAL', 'sirolimus::ORAL']) {
  assert.ok(built.presets[key].candidates.length > 1, key);
  assert.equal(built.presets[key].autoApply, false, key);
  assert.ok(built.presets[key].doseGrid, key + ' needs a declared grid');
  assert.equal(built.presets[key].doseGrid.signoff, 'PENDING');
}
// The liver-transplant tacrolimus reference is flagged elsewhere and must never be offered here.
assert.ok(!built.presets['tacrolimus::ORAL'].candidates
  .some((candidate) => candidate.id === 'tac-liver-everolimus-fda'),
  'A context-mismatched reference leaked into the preset list');
assert.deepEqual(built.presets['tacrolimus::ORAL'].candidates.map((c) => c.id),
  ['tac-kidney-early', 'tac-kidney-late', 'tac-kidney-mid']);
// AUC0-12 tacrolimus references are excluded, so no preset claims an AUC24 for them.
assert.ok(built.presets['tacrolimus::ORAL'].candidates.every((c) => c.target_metric === 'Cmin'));
// Beta-lactams stay out entirely: fT>MIC is not a two-bound concentration target.
assert.ok(!keys.some((key) => /^(amox|cefepime|cefazoline)::/.test(key)));

// Fail closed on a tampered snapshot or missing declarations.
assert.throws(() => buildClinicalPresets({ catalogue: cat, crosswalk: map, rules,
  doseGrids: grids, snapshot: { ...snapshot, sha256: '0'.repeat(64) } }), /Invalid PMx snapshot/);
assert.throws(() => buildClinicalPresets({ catalogue: cat, crosswalk: map, rules,
  doseGrids: { ...grids, schema: 'wrong' }, snapshot }), /dose grids/);
assert.throws(() => buildClinicalPresets({ catalogue: cat, crosswalk: map,
  rules: { ...rules, schema: 'wrong' }, doseGrids: grids, snapshot }), /Applicability rules/);

console.log('TDM preset artifact: ' + keys.length + ' keys, ' + all.length +
  ' published targets, auto-applied only where a single reference qualifies and a curated preset already exists.');

// --- An unvalidated dose grid must never be applicable ---
import { createHash } from 'node:crypto';
// The committed artifact must match its own digest: a hand edit is caught even if
// the sha256 field was updated to match the tampered body.
const { sha256: declaredSha, ...committedBody } = committed;
assert.equal(createHash('sha256').update(JSON.stringify(committedBody)).digest('hex'),
  declaredSha, 'Committed preset artifact does not match its own sha256');
// Every declared grid is still pending, so none of them may reach the dose inputs.
for (const [key, entry] of Object.entries(built.presets)) {
  assert.equal(entry.doseGridApplicable, entry.doseGrid?.signoff === 'APPROVED', key);
  assert.equal(entry.doseGridApplicable, false, key + ' must not be applicable while pending');
}
assert.equal(Object.values(built.presets).some((entry) => entry.doseGridApplicable), false,
  'No dose grid may be applicable before an explicit pharmacologist sign-off');
// Signing one off is what makes it applicable, and nothing else changes.
const signedOff = buildClinicalPresets({ catalogue: cat, crosswalk: map, rules, snapshot,
  doseGrids: { ...grids, grids: { ...grids.grids,
    'tacrolimus::ORAL': { ...grids.grids['tacrolimus::ORAL'], signoff: 'APPROVED' } } } });
assert.equal(signedOff.presets['tacrolimus::ORAL'].doseGridApplicable, true);
assert.equal(signedOff.presets['ciclosporine::ORAL'].doseGridApplicable, false);
assert.deepEqual(signedOff.presets['tacrolimus::ORAL'].candidates,
  built.presets['tacrolimus::ORAL'].candidates, 'Sign-off must not change the targets');
// The corrected grids are what the artifact carries.
assert.equal(built.presets['sirolimus::ORAL'].doseGrid.dose_step, 1);
assert.deepEqual(built.presets['tacrolimus::ORAL'].doseGrid.intervals, [12]);
assert.equal(built.presets['tacrolimus::ORAL'].doseGrid.dose_max, 14);
console.log('Dose grids: all PENDING and non-applicable, sign-off is the only switch, label constraints enforced.');
