// Covers every model without needing to import, install or deploy PMx Explain.
import assert from 'node:assert/strict';
import {readFile} from 'node:fs/promises';
import {createHash} from 'node:crypto';
import {buildModelTdmPreview,checkReferenceSnapshot,validateCrosswalk,checkApplicabilityRules,assessCandidateApplicability} from '../scripts/preview_tdm_guide_recommendations.mjs';
const cat=JSON.parse(await readFile(new URL('../models/catalog.json',import.meta.url),'utf8'));
const map=JSON.parse(await readFile(new URL('../integration/tdm-guide-drug-crosswalk.preview.v1.json',import.meta.url),'utf8'));
const rules=JSON.parse(await readFile(new URL('../integration/tdm-guide-applicability-rules.preview.v1.json',import.meta.url),'utf8'));
assert.equal(checkApplicabilityRules(rules),true);
assert.equal(cat.models.length,46,'Frozen 46 models must be covered before adding Limoges models');
assert.equal(new Set(cat.models.map(m=>m.id)).size,46);
assert.equal(validateCrosswalk(map,cat),true);
const referenced=new Set(['tac','vancomycin','cefepime','fluconazole']);
const missing=new Set(['amiodarone','cobicistat','diltiazem','erythromycin','ritonavir']);
const drugs=[...new Set(Object.values(map.keys))].filter(d=>!missing.has(d)).map(id=>({
  id,fr:id,en:id,
  index:{key:id==='cefepime'?'fT_MIC':'unknown',fr:'Indice',en:'Index'},
  practicalSampling:{key:'unknown',fr:'À documenter',en:'Not assessed'},
  pkReview:'class_reviewed',
  targets:id==='tac'||id==='vancomycin'?[{
    id:id+'-original-target',metric:'AUC',unit:'mg·h/L',
    indication:'study_only',indicationFr:'Étude',indicationEn:'Study',
    relation:'range',min:10,max:20,
    source:{id:'original',citation:'Synthetic test',url:'https://example.org/ref',year:2026,section:'Synthetic'},
    missingAnalyticalIdentity:['assayMethod'],humanSignoff:'PENDING',
    interpretationAllowed:false,automaticDoseAdjustmentAllowed:false
  }]:[],
  pkpd:id==='cefepime'?[{
    id:'cefepime:icu-beta-continuous',kind:'beta_continuous',metric:'fCss/MIC',
    labelFr:'ICU',labelEn:'ICU',criterionFr:'4 CMI',criterionEn:'4 MIC',
    source:{url:'https://example.org/pd'},isContextualOnly:true,interpretationAllowed:false
  }]:[],
  society:id==='fluconazole'?{recommendation:'neither_recommend_nor_discourage',source:{url:'https://example.org/paper'}}:null
}));
assert.equal(drugs.length,23-missing.size);
const payload={schema:'pmx.tdm.model-references.preview.v1',producer:'synthetic-only',
  intendedUse:'research_education_bibliographic_preview',
  notForClinicalInterpretation:true,patientDataIncluded:false,
  containsDosingAdvice:false,modelNamesIncluded:false,
  referenceApproval:'INDEPENDENT_PHARMACOLOGIST_SIGNOFF_PENDING',drugs};
const snapshot={...payload,sha256:createHash('sha256').update(JSON.stringify(payload)).digest('hex')};
assert.equal(checkReferenceSnapshot(snapshot),true);
const out=buildModelTdmPreview({catalogue:cat,crosswalk:map,rules,snapshot});
assert.equal(out.models.length,46);
assert.equal(new Set(out.models.map(m=>m.modelId)).size,46);
assert.equal(out.enabledInRuntime,false);
assert.equal(out.clinicalInterpretationAllowed,false);
assert.ok(out.models.every(m=>m.recommendation.interpretationAllowed===false));
assert.ok(out.models.every(m=>m.recommendation.automaticDoseAdjustmentAllowed===false));
assert.ok(out.models.every(m=>m.recommendation.evidenceSourceVersion===snapshot.sha256));
assert.equal(out.models.filter(m=>m.recommendation.status==='not_indexed_in_guide').length,5);
assert.equal(out.models.find(m=>m.modelId==='tacrolimus_woillard_ddi')?.recommendation.status,'published_target_candidates');
assert.equal(out.models.find(m=>m.modelId==='vanco_roberts')?.recommendation.status,'published_target_candidates');
assert.equal(out.models.find(m=>m.modelId==='cefepime_an')?.recommendation.status,'contextual_pkpd_only');
assert.equal(out.models.find(m=>m.modelId==='cefepime_an')?.recommendation.pharmacodynamicCandidates[0].status,'different_administration_mode',
  'Continuous-infusion PD guidance must not automatically match intermittent model');
assert.equal(out.models.find(m=>m.modelId==='fluconazole_debruyne_ddi')?.recommendation.status,'society_tdm_context_only');
assert.equal(out.models.find(m=>m.modelId==='ritonavir_kappelhoff_ddi')?.recommendation.status,'not_indexed_in_guide');
assert.equal(out.models.filter(m=>m.simulationOnly===false).length,13);
assert.equal(checkReferenceSnapshot({...snapshot,sha256:'0'.repeat(64)}),false);
assert.throws(()=>buildModelTdmPreview({catalogue:cat,crosswalk:map,rules,snapshot:{...snapshot,sha256:'0'.repeat(64)}}),/Invalid PMx snapshot/);
const tampered={...map,keys:{...map.keys}};
delete tampered.keys.tacrolimus;
assert.equal(validateCrosswalk(tampered,cat),false);
assert.throws(()=>buildModelTdmPreview({catalogue:cat,crosswalk:tampered,rules,snapshot}),/crosswalk/);
const newModel={...cat.models[0],id:'new_ltc_model',drugKey:'new_limoges_molecule'};
assert.throws(()=>buildModelTdmPreview({catalogue:{models:[...cat.models,newModel]},crosswalk:map,rules,snapshot}),/crosswalk/,
  'New Limoges models must explicitly establish their guide drug ID');
console.log('MIPD TDM preview: 46/46 cards, 23 drug keys, 13 MAP flags, no patient-facing interpretation.');

// --- Structural applicability gates: a warning may be added, never an authorisation ---
assert.equal(checkApplicabilityRules({...rules,schema:'wrong'}),false);
assert.equal(checkApplicabilityRules({...rules,organTransplantIndications:undefined}),false);
assert.throws(()=>buildModelTdmPreview({catalogue:cat,crosswalk:map,snapshot,
  rules:{...rules,schema:'wrong'}}),/Applicability rules/);
assert.throws(()=>buildModelTdmPreview({catalogue:cat,crosswalk:map,snapshot,rules:undefined}),
  /Applicability rules/);

const kidneyModel={populationTags:['Adulte','Transplantation rénale'],administrationModes:['ORAL']};
const noOrganModel={populationTags:['Adulte'],administrationModes:['ORAL']};
const continuousModel={populationTags:['Adulte','ICU'],administrationModes:['IV_CONTINUOUS']};
const intermittentModel={populationTags:['Adulte','ICU'],administrationModes:['IV_INTERMITTENT']};
const ref=(over={})=>({metric:'C₀',indication:'kidney',missingAnalyticalIdentity:[],...over});

// A liver-transplant reference on a kidney-transplant model is reported, not silently listed.
const liverOnKidney=assessCandidateApplicability(
  {target:ref({indication:'liver'}),model:kidneyModel,tdmDrugId:'tac',rules});
assert.equal(liverOnKidney.verdict,'structural_context_mismatch');
assert.equal(liverOnKidney.mismatches.length,1);
assert.equal(liverOnKidney.mismatches[0].code,'organ_indication_mismatch');
assert.equal(liverOnKidney.mismatches[0].referenceContext,'liver');
assert.equal(liverOnKidney.mismatches[0].modelContext,'kidney');
// The matching organ raises no warning and is still never marked applicable.
const kidneyOnKidney=assessCandidateApplicability(
  {target:ref(),model:kidneyModel,tdmDrugId:'tac',rules});
assert.deepEqual(kidneyOnKidney.mismatches,[]);
assert.equal(kidneyOnKidney.verdict,'requires_human_context_review');
// An undeclared organ is undecidable, never a mismatch: the gate stays conservative.
const unknownOrgan=assessCandidateApplicability(
  {target:ref({indication:'liver'}),model:noOrganModel,tdmDrugId:'tac',rules});
assert.deepEqual(unknownOrgan.mismatches,[]);
assert.ok(unknownOrgan.undecidable.includes('transplanted_organ_not_declared_by_model'));
// Out-of-scope drugs keep the same indication codes with another meaning (digoxin heart failure).
assert.deepEqual(assessCandidateApplicability(
  {target:ref({indication:'heart'}),model:kidneyModel,tdmDrugId:'digoxin',rules}).mismatches,[]);
// A trough cannot be read on a continuous-infusion-only model.
const troughOnContinuous=assessCandidateApplicability(
  {target:ref({indication:'infection'}),model:continuousModel,tdmDrugId:'vancomycin',rules});
assert.equal(troughOnContinuous.mismatches[0].code,
  'trough_reference_under_continuous_infusion_only_model');
assert.deepEqual(assessCandidateApplicability(
  {target:ref({indication:'infection'}),model:intermittentModel,tdmDrugId:'vancomycin',rules}
).mismatches,[]);
// An AUC reference stays readable under continuous infusion.
assert.deepEqual(assessCandidateApplicability(
  {target:ref({metric:'AUC₂₄',indication:'infection'}),model:continuousModel,
    tdmDrugId:'vancomycin',rules}).mismatches,[]);
// Unstructured identity dimensions are reported as undecidable, never as a match.
assert.deepEqual(assessCandidateApplicability({target:ref(
  {missingAnalyticalIdentity:['assayMethod','formulationId']}),model:kidneyModel,
  tdmDrugId:'tac',rules}).undecidable,['assayMethod','formulationId']);
// No gate may ever conclude that a reference is applicable.
assert.ok(out.models.flatMap(m=>m.recommendation.numericalReferenceCandidates)
  .every(c=>['requires_human_context_review','structural_context_mismatch']
    .includes(c.applicability.verdict)&&c.modelApplicabilityVerified===false));
assert.ok(out.models.every(m=>typeof m.recommendation.structuralMismatchCount==='number'&&
  typeof m.recommendation.administrationModeWarningCount==='number'));
assert.equal(out.models.find(m=>m.modelId==='cefepime_an').recommendation
  .administrationModeWarningCount,1,'Mode warning must be countable for the read-only panel');
console.log('MIPD applicability gates: organ, administration mode and undecidable dimensions reported; none grants applicability.');
