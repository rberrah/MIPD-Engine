// Research / education PREVIEW ONLY. Not loaded by the Shiny runtime.
// The standalone model catalog owns the models; PMx TDM Guide owns STP references.
// This adapter never selects a therapeutic value, interprets a patient, or adjusts a dose.
import assert from 'node:assert/strict';
import { readFile, writeFile } from 'node:fs/promises';
import { createHash } from 'node:crypto';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
export const MODEL_PREVIEW_SCHEMA='mipd.tdm.model-recommendations.preview.v1';
const REFERENCE_SCHEMA='pmx.tdm.model-references.preview.v1';
const sha256=(value)=>createHash('sha256').update(JSON.stringify(value)).digest('hex');

export function checkReferenceSnapshot(snapshot) {
  if (!snapshot || snapshot.schema!==REFERENCE_SCHEMA ||
    snapshot.notForClinicalInterpretation!==true || snapshot.patientDataIncluded!==false ||
    snapshot.containsDosingAdvice!==false || snapshot.modelNamesIncluded!==false ||
    !Array.isArray(snapshot.drugs) || snapshot.referenceApproval!=='INDEPENDENT_PHARMACOLOGIST_SIGNOFF_PENDING') return false;
  const {sha256:declared,...body}=snapshot;
  if (typeof declared!=='string'||!/^[a-f0-9]{64}$/.test(declared)||sha256(body)!==declared) return false;
  const keys=new Set();
  for(const drug of snapshot.drugs) {
    if (typeof drug.id!=='string'||!drug.id||keys.has(drug.id)||
      !Array.isArray(drug.targets)||!Array.isArray(drug.pkpd)||!drug.index||!drug.practicalSampling) return false;
    keys.add(drug.id);
    if (drug.targets.some((target)=>target.interpretationAllowed!==false||
      target.automaticDoseAdjustmentAllowed!==false||target.humanSignoff!=='PENDING'||
      !target.source?.url?.startsWith('https://'))) return false;
    if (drug.pkpd.some((pd)=>pd.isContextualOnly!==true||pd.interpretationAllowed!==false)) return false;
  }
  return true;
}
export function validateCrosswalk(crosswalk,catalogue) {
  if (crosswalk?.schema!=='mipd.tdm.drug-crosswalk.preview.v1'||!crosswalk.keys||
    typeof crosswalk.keys!=='object'||!Array.isArray(catalogue?.models)) return false;
  const used=new Set(catalogue.models.map((m)=>m.drugKey));
  // No silent skip: a future local CHU model MUST receive a reviewed canonical key.
  return used.size===Object.keys(crosswalk.keys).length &&
    [...used].every((key)=>typeof crosswalk.keys[key]==='string'&&crosswalk.keys[key]) &&
    Object.keys(crosswalk.keys).every((key)=>used.has(key));
}
const pdModeWarning=(pd,modes)=>{
  const m=Array.isArray(modes)?modes:[];
  if(pd.kind==='beta_continuous' && !m.includes('IV_CONTINUOUS')) return 'different_administration_mode';
  if(pd.kind==='beta_trough' && !m.includes('IV_INTERMITTENT')) return 'different_administration_mode';
  return 'population_indication_and_assay_not_checked';
};

export function checkApplicabilityRules(rules) {
  const organ=rules?.organTransplantIndications;
  return rules?.schema==='mipd.tdm.applicability-rules.preview.v1' &&
    Array.isArray(organ?.appliesToTdmDrugIds) &&
    !!organ?.organPopulationTags && typeof organ.organPopulationTags==='object' &&
    Object.values(organ.organPopulationTags).every((tags)=>Array.isArray(tags)&&tags.length) &&
    Array.isArray(rules.troughMetricMarkers) && Array.isArray(rules.continuousInfusionModes) &&
    Array.isArray(rules.undecidableIdentityDimensions?.fields);
}

/**
 * Structural, conservative context gates. A gate can only ADD a warning: it never
 * marks a reference as applicable, never picks a value and never touches a dose.
 * Anything the declared metadata cannot settle is reported as undecidable, not as a match.
 */
export function assessCandidateApplicability({target,model,tdmDrugId,rules}) {
  const mismatches=[],undecidable=[];
  const organRule=rules.organTransplantIndications;
  const organTags=organRule.organPopulationTags;
  const tags=Array.isArray(model.populationTags)?model.populationTags:[];
  const modes=Array.isArray(model.administrationModes)?model.administrationModes:[];
  // 1. Grafted organ declared by the reference vs grafted organ declared by the model.
  if(organRule.appliesToTdmDrugIds.includes(tdmDrugId) && organTags[target.indication]) {
    const wanted=organTags[target.indication];
    const declared=Object.entries(organTags)
      .filter(([,list])=>list.some((tag)=>tags.includes(tag))).map(([organ])=>organ);
    if(!wanted.some((tag)=>tags.includes(tag))) {
      if(declared.length) mismatches.push({
        dimension:'transplanted_organ',code:'organ_indication_mismatch',
        referenceContext:target.indication,modelContext:declared.join('+')
      });
      else undecidable.push('transplanted_organ_not_declared_by_model');
    }
  }
  // 2. A trough cannot be read on a model that only supports continuous infusion.
  const isTrough=rules.troughMetricMarkers.some((marker)=>
    typeof target.metric==='string' && target.metric.includes(marker));
  if(isTrough && modes.length && modes.every((mode)=>rules.continuousInfusionModes.includes(mode)))
    mismatches.push({
      dimension:'administration_mode',code:'trough_reference_under_continuous_infusion_only_model',
      referenceContext:target.metric,modelContext:modes.join('+')
    });
  // 3. Dimensions the guide leaves unstructured stay explicitly undecidable.
  for(const field of (target.missingAnalyticalIdentity||[]))
    if(rules.undecidableIdentityDimensions.fields.includes(field)) undecidable.push(field);
  return {
    // Never 'applicable': only a human review of the measurement conditions can conclude.
    verdict:mismatches.length?'structural_context_mismatch':'requires_human_context_review',
    mismatches,undecidable
  };
}
export function buildModelTdmPreview({catalogue,snapshot,crosswalk,rules}) {
  if(!Array.isArray(catalogue?.models)||!catalogue.models.length) throw new Error('Missing independent model catalog');
  if(!checkReferenceSnapshot(snapshot)) throw new Error('Invalid PMx snapshot or missing immutable SHA256');
  if(!validateCrosswalk(crosswalk,catalogue)) throw new Error('Model-drug crosswalk incomplete or stale');
  if(!checkApplicabilityRules(rules)) throw new Error('Applicability rules missing or stale');
  const uniqueModels=new Set();
  const guide=new Map(snapshot.drugs.map((d)=>[d.id,d]));
  const models=catalogue.models.map((model)=>{
    if(uniqueModels.has(model.id)) throw new Error('Duplicate model ID '+model.id);
    uniqueModels.add(model.id);
    const tdmId=crosswalk.keys[model.drugKey],entry=guide.get(tdmId);
    const targetCandidates=(entry?.targets||[]).map((target)=>({
      id:target.id,metric:target.metric,unit:target.unit,
      indication:target.indication,indicationFr:target.indicationFr,indicationEn:target.indicationEn,
      contextFr:target.contextFr,contextEn:target.contextEn,
      analyteId:target.analyteId,analyteBasis:target.analyteBasis,
      matrixId:target.matrixId,matrixFr:target.matrixFr,matrixEn:target.matrixEn,
      formulationId:target.formulationId,samplingTime:target.samplingTime,
      relation:target.relation,min:target.min,max:target.max,
      analyticalIdentityMissing:target.missingAnalyticalIdentity||[],
      source:target.source,status:'bibliographic_candidate_requires_context_review',
      // Structural gates: drug ID alone never establishes applicability.
      applicability:assessCandidateApplicability({target,model,tdmDrugId:tdmId,rules}),
      // Clinical applicability is NEVER asserted by drug ID alone.
      modelApplicabilityVerified:false,automaticDoseAdjustmentAllowed:false
    }));
    const pkpdCandidates=(entry?.pkpd||[]).map((pd)=>({
      id:pd.id,kind:pd.kind,metric:pd.metric,labelFr:pd.labelFr,labelEn:pd.labelEn,
      criterionFr:pd.criterionFr,criterionEn:pd.criterionEn,
      source:pd.source,status:pdModeWarning(pd,model.administrationModes),
      modelApplicabilityVerified:false,automaticDoseAdjustmentAllowed:false
    }));
    const recommendationStatus= !entry ? 'not_indexed_in_guide'
      : targetCandidates.length ? 'published_target_candidates'
      : pkpdCandidates.length ? 'contextual_pkpd_only'
      : entry.society ? 'society_tdm_context_only'
      : 'no_qualified_stp_reference';
    return {
      modelId:model.id,modelDrugKey:model.drugKey,tdmDrugId:tdmId,
      modelName:model.model,drugName:model.drug,
      population:model.population,populationEn:model.populationEn??null,
      populationTags:model.populationTags,
      routes:model.routes,administrationModes:model.administrationModes,
      simulationOnly:!model.analysisEligible,
      verificationStatus:model.verificationStatus,
      clinicalValidationStatus:model.validationStatus,
      recommendation:{
        status:recommendationStatus,
        index:entry?.index??{key:'unknown',fr:'Non référencé',en:'Not indexed'},
        sampling:entry?.practicalSampling??{key:'unknown',fr:'Non référencé',en:'Not indexed'},
        monitoringContext:entry?.society??null,
        numericalReferenceCandidates:targetCandidates,
        pharmacodynamicCandidates:pkpdCandidates,
        // Counts let the read-only panel surface warnings before any number.
        structuralMismatchCount:targetCandidates
          .filter((c)=>c.applicability.mismatches.length).length,
        administrationModeWarningCount:pkpdCandidates
          .filter((c)=>c.status==='different_administration_mode').length,
        evidenceSourceVersion:snapshot.sha256,
        patientContextMatched:false,
        modelFormulationMatched:false,
        analyteAndMatrixMatched:false,
        assayMethodMatched:false,
        interpretationAllowed:false,
        automaticDoseAdjustmentAllowed:false
      }
    };
  });
  if(models.length!==catalogue.models.length) throw new Error('Model coverage failure');
  const body={
    schema:MODEL_PREVIEW_SCHEMA,
    producer:'MIPD-Engine standalone draft adapter',
    sourceSnapshotSha256:snapshot.sha256,
    intendedUse:'research_education_preview_only',
    enabledInRuntime:false,patientDataIncluded:false,
    clinicalInterpretationAllowed:false,automaticDoseAdjustmentAllowed:false,
    models
  };
  return {...body,sha256:sha256(body)};
}
export async function loadModelTdmPreview(snapshotFile) {
  if(!snapshotFile||!snapshotFile.endsWith('.json')) throw new Error('Explicit local JSON snapshot path required');
  const [cat,map,gates,snap]=await Promise.all([
    readFile(path.join(root,'models/catalog.json'),'utf8'),
    readFile(path.join(root,'integration/tdm-guide-drug-crosswalk.preview.v1.json'),'utf8'),
    readFile(path.join(root,'integration/tdm-guide-applicability-rules.preview.v1.json'),'utf8'),
    readFile(snapshotFile,'utf8')
  ]);
  return buildModelTdmPreview({catalogue:JSON.parse(cat),crosswalk:JSON.parse(map),
    rules:JSON.parse(gates),snapshot:JSON.parse(snap)});
}
if(process.argv[1]?.endsWith('preview_tdm_guide_recommendations.mjs')) {
  const a=process.argv.indexOf('--snapshot');
  if(a<0||!process.argv[a+1]) throw new Error('Usage: --snapshot /local/verified-guide.json [--output /tmp/preview.json]');
  const result=await loadModelTdmPreview(process.argv[a+1]);
  const out=process.argv.indexOf('--output');
  if(out>=0) {
    const filename=process.argv[out+1];
    if(!filename?.endsWith('.json')) throw new Error('Output filename must end with .json');
    await writeFile(filename,JSON.stringify(result,null,2)+'\n','utf8');
  } else {
    process.stdout.write(JSON.stringify(result,null,2)+'\n');
  }
}
