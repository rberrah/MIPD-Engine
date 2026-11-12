import assert from 'node:assert/strict';
import {readFileSync, readdirSync} from 'node:fs';
import {join} from 'node:path';

const read = path => JSON.parse(readFileSync(path, 'utf8'));
const cat = read('models/catalog.json').models;
const cross = read('integration/tdm-guide-drug-crosswalk.preview.v1.json').keys;
const guide = read('integration/tdm-guide-clinical-presets.v1.json');
const listing = readdirSync('integration/tdm-drugs').filter(n => n.endsWith('.json')).sort();
const keys = Object.keys(cross).sort();
assert.deepEqual(listing, keys.map(k=>k+'.json'));
assert.equal(keys.length,23);
assert.equal(cat.length,46);
const profiles = {};
const modeCount = {};
let qualified = 0;
for(const key of keys){
  const d = read(join('integration/tdm-drugs',key+'.json'));
  profiles[key] = d;
  const models = cat.filter(x => x.drugKey === key);
  assert.equal(d.schema,'mipd.tdm.drug-profile.v1');
  assert.equal(d.drugKey,key);
  assert.equal(d.tdmDrugId,cross[key]);
  assert.deepEqual([...d.modelIds].sort(),models.map(m=>m.id).sort());
  assert.equal(new Set(models.map(x=>x.concentrationUnit)).size,1);
  assert.equal(d.concentrationUnit,models[0].concentrationUnit);
  assert.equal(d.doseUnit,'mg');
  assert.equal(d.intendedUse,'TEACHING_RESEARCH_ONLY');
  assert.equal(d.clinicalInterpretationAllowed,false);
  assert.equal(d.automaticDoseAdjustmentAllowed,false);
  assert.equal(d.patientDataIncluded,false);
  assert.ok(d.guideHref.startsWith('https://rberrah.github.io/pharmacometrie/cibles-tdm/?drug='));
  assert.match(d.sourceSnapshotSha256,/^[a-f0-9]{64}$/);
  for(const [mode,entry] of Object.entries(d.guideTargetPresets)){
    modeCount[key+'::'+mode] = true;
    assert.deepEqual(entry,guide.presets[key+'::'+mode],
      'Drug dossier must match the reviewed central evidence snapshot');
    qualified += entry.candidates.length;
    for(const c of entry.candidates){
      assert.equal(c.modelApplicabilityVerified,false);
      assert.equal(c.automaticDoseAdjustmentAllowed,false);
      assert.ok(c.target_high>c.target_low);
      assert.match(c.url,/^https:\/\//);
    }
    if(entry.doseGrid){
      assert.equal(entry.doseGrid.signoff,'PENDING');
      assert.ok(entry.doseGrid.dose_min>=0);
      assert.ok(entry.doseGrid.dose_step>0);
      assert.match(entry.doseGrid.source.url,/^https:\/\//);
    }
  }
  for(const v of d.sourceReferences) {
    assert.equal(v.modelApplicabilityVerified,false);
    assert.equal(v.automaticDoseAdjustmentAllowed,false);
  }
  for(const [mode,preset] of Object.entries(d.simulationPresets)){
    assert.ok(['ORAL','IV_INTERMITTENT','IV_CONTINUOUS','IM'].includes(mode));
    assert.ok(['Cmin','Cmax','AUC24'].includes(preset.target_metric));
    assert.ok(preset.target_low>=0 && preset.target_high>preset.target_low);
    assert.ok(preset.dose_min>=0 && preset.dose_max>preset.dose_min);
    assert.ok(preset.dose_step>0);
    assert.ok([].concat(preset.intervals).every(v=>v>0));
    assert.match('https://doi.org/'+preset.doi,/^https:\/\//);
  }
}
assert.deepEqual(Object.keys(modeCount).sort(),Object.keys(guide.presets).sort());
assert.equal(qualified,guide.presets && Object.values(guide.presets).reduce((a,x)=>a+x.candidates.length,0));
const tac=profiles.tacrolimus;
assert.equal(tac.concentrationUnit,'ng/mL');
assert.equal(tac.modelIds.length,2);
assert.equal(tac.simulationPresets.ORAL,undefined,'Trough is context-specific: no arbitrary default');
assert.equal(tac.guideTargetPresets.ORAL.candidates.length,3);
assert.equal(tac.guideTargetPresets.ORAL.doseGrid.signoff,'PENDING');
assert.deepEqual(tac.guideTargetPresets.ORAL.doseGrid.intervals,[12]);
assert.equal(tac.formulationCovariates.tacrolimus_woillard_ddi.values['1'].product,'Prograf');
assert.equal(tac.formulationCovariates.tacrolimus_woillard_ddi.values['0'].product,'Advagraf');
assert.equal(tac.formulationCovariates.tacrolimus_kamp_envarsus_ddi.product,'Envarsus');
const vanco=profiles.vanco;
assert.equal(vanco.simulationPresets.IV_INTERMITTENT.target_low,400);
assert.equal(vanco.simulationPresets.IV_INTERMITTENT.target_high,600);
assert.equal(vanco.simulationPresets.IV_INTERMITTENT.target_metric,'AUC24');
assert.equal(profiles.linez.simulationPresets.ORAL.target_metric,'Cmin');
assert.equal(profiles.voriconazole.simulationPresets.ORAL.target_high,5.5);

const app=readFileSync('app.R','utf8');
const clinicalR=readFileSync('R/clinical_presets.R','utf8');
const guideR=readFileSync('R/tdm_guide_presets.R','utf8');
assert.ok(!app.includes('No universal clinical preset'));
assert.ok(app.includes('output$tdm_drug_profile_ui <- renderUI'));
assert.ok(app.includes('tdm_drug_profile_grid_eligible'));
assert.ok(app.includes('tdm_drug_profile_metric_unit'));
assert.ok(app.includes('if (length(modes) <= 1L) return(NULL)'));
assert.ok(clinicalR.includes('TDM_DRUG_PROFILES'));
assert.ok(guideR.includes('tdm_drug_profile_guide_entry'));
for(const name of ['tacrolimus_woillard_ddi','tacrolimus_kamp_envarsus_ddi']){
 const cpp=readFileSync('models/'+name+'.cpp','utf8');
 assert.ok(cpp.includes('double CONC = CENT_TAC / (V1/1000.0);'));
 assert.ok(!cpp.includes('CONC = CONC / 1000'), 'No additional conversion on model ng/mL values');
}
console.log('23 drug dossiers cover 46/46 models, all context references and 5 historical presets.');
console.log('Tacrolimus: ng/mL, Prograf ST=1, Advagraf ST=0, Envarsus distinct; source label grids not auto-prescribed.');
