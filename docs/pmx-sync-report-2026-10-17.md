# PMx Explain synchronization handoff - 2026-10-17

This document records the MIPD Engine release that PMx Explain may consume. It
does not claim that PMx Explain has already imported or validated the release.

## Qualified public release

| Field | Value |
|---|---|
| Repository | `rberrah/MIPD-Engine` |
| Git SHA | `fec16cbdf8e33ed6c791226ccf236761333f8de8` |
| Git tag | `public-2026-10-17-fec16cb` |
| Public URL | `https://tdmhub.shinyapps.io/MIPD_Engine/` |
| Shinyapps account/application | `tdmhub/MIPD_Engine` |
| Application ID | `17776160` |
| Bundle | `12702795` |
| Deployment task | `1758931558` |
| Image | `15803034` |
| Rollback bundle retained | `12702610` |
| Production flags | `ALLOW_CUSTOM_MODELS=false`; `ENABLE_PMETRICS=false` |

The deployed `/deployment-release.json` returned the exact immutable SHA above.
No new Shinyapps application was created.

## Public acceptance results

All checks below passed against the public URL after bundle `12702795` became
active:

- HTTP application load and release fingerprint;
- WebSocket initialization and two-session isolation;
- French and English initialization, with 46 model cards in each language;
- default synthetic PK simulation and MAP fit;
- MAP outputs including steady-state AUC24 `573.7` and 24-hour AUC `417.2`;
- DDI model generation, compilation and simulation;
- DDI clearance target `TVCL_TAC`, with AUC24 `283.0` without interaction and
  `560.0` with interaction (ratio `1.98`);
- public rejection of arbitrary custom C++ and disabled Pmetrics execution.

The DDI test specifically exercises the production fix in this release:
server-generated DDI C++ compiles, while user-appended C++ after the Lego marker
is discarded before compilation.

These are software acceptance checks. They are not clinical validation or an
independent reproduction of every source publication.

## PMx Explain baseline observed before handoff

At PMx Explain commit `5fdb0c94c91ce101badb0e3c13f2106f2b54847e`, the committed
`integration/mipd-release.json` pins MIPD Engine SHA
`83a3026df7c0881d2725333a90d81ef2ed218f90`.

An uncommitted PMx Explain working tree observed on 2026-10-17 instead pins
`3302020159361edc74dcba4cfd3aa003b0ef90e8` and includes source-review sync work.
That work must be reviewed and rebased or superseded intentionally. It must not
be overwritten blindly. No PMx Explain file was modified while preparing this
handoff.

## Model-source delta from the committed PMx pin

Exactly 18 executable model files differ between the committed PMx pin
`83a3026` and the qualified MIPD release `fec16cb`.

| Model ID | Current source-review result | Material release change or retained limit |
|---|---|---|
| `amiodarone_pollak_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Published disposition values and reported IIV replace mismatched values; absorption and residual terms remain assumptions. |
| `cefazoline_cheng` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Final unbound equations and Q exponent 0.78 are retained; published occasion variability is omitted. |
| `cefazoline_lanoiselee` | `ARTICLE_ALIGNED` | Two compartments, CKD-EPI effect, correlated CL-Vc block, Q/Vp variability and proportional error use the published Monolix scale. |
| `cefepime_an` | `ARTICLE_ALIGNED` | Renal, CRRT and nonrenal clearance, volumes, IIV and combined residual error align with the final model. |
| `cefepime_seo` | `ARTICLE_ALIGNED` | Structure, Cockcroft-Gault effect, fixed effects, IIV and proportional error align with Table 2. |
| `ciclosporine_press_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Published weight, prednisolone, IIV and residual effects retained; F IOV omitted and CYP3A fraction remains interaction-only. |
| `cobicistat_mathias_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Published 50 mg steady-state NCA profile only; ka and random effects remain assumptions, not a PopPK reproduction. |
| `diltiazem_hermann_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Published mean CL, V and F represented by a documented one-compartment approximation; ka is assumed. |
| `erythromycine_mather_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Base post-meal F is 0.40; CL, V and ka remain declared assumptions because the source is not a reusable PopPK model. |
| `fluconazole_debruyne_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Review-derived mean CL, F and half-life-derived V; ka and random effects remain simulation assumptions. |
| `isavuconazole_desai_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Weibull absorption and reported covariates/IIV retained; default patient and residual branches replace study-specific alternatives. |
| `isavuconazole_groll_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Rounded mean-profile reconstruction from an NCA DDI study, not a reusable article PopPK model or MAP prior. |
| `posaconazole_dolton_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Published lag model and medication/clinical covariates replace mismatched values; one study branch supplies default IIV/error. |
| `rifampicine_loos_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Published chronic-dose CL and F retained; V and ka remain explicit assumptions. |
| `rifampicine_marsot` | `SOURCE_AMBIGUITY_DOCUMENTED` | Additive 2.256 mg/L is treated as an SD and squared; source text also mentions log-transformed data. |
| `ritonavir_kappelhoff_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Published values, lag, lopinavir effect and IIV retained; F IOV, V-ka covariance and residual mixture are omitted. |
| `tacrolimus_kamp_envarsus_ddi` | `DOCUMENTED_ADAPTATION_CONFIRMED` | Published fixed effects, allometry, CL-V1 block, KTR variability and lag class retained; CL IOV and residual switching are omitted. |
| `vanco_goti` | `ARTICLE_ALIGNED` | Primary-source age >60 creatinine floor and CrCl cap are corrected; independent synthetic equation and MAP-oracle tests cover both boundaries. |

The authoritative details are in `models/source-review.json`,
`models/metadata.json`, `models/i18n.en.json`, `models/catalog.json` and
`docs/tdm-model-traceability.json` at the release SHA. `ARTICLE_ALIGNED` is not
equivalent to full source reproduction or external validation.

## Scientific qualification at the release SHA

- 46 executable models pass software compilation/contract checks.
- 13 models have experimental MAP priors.
- 11 of those 13 are marked `ARTICLE_ALIGNED` by author-led source review.
- `rifampicine_marsot` is `SOURCE_AMBIGUITY_DOCUMENTED`.
- `tacrolimus_woillard_ddi` is `DOCUMENTED_ADAPTATION_CONFIRMED`.
- Roberts, Revilla and Goti have partial independent source-equation and MAP
  checks on synthetic scenarios.
- 0/46 models have a claimed full independent source reproduction.
- 0 external clinical validations are claimed.

Goti remains source-aligned but not externally validated. Its dependent XGBoost
and DALEX artifacts were retrained after the renal-rule correction, and their
provenance points to the current model hash.

## ML synchronization limits

The release deliberately fail-closes ML artifacts whose model SHA no longer
matches the corrected C++ source. The following 12 artifacts are therefore not
available in TDM until they are retrained and requalified:

`amiodarone_pollak_ddi`, `ciclosporine_press_ddi`,
`cobicistat_mathias_ddi`, `diltiazem_hermann_ddi`,
`erythromycine_mather_ddi`, `fluconazole_debruyne_ddi`,
`isavuconazole_desai_ddi`, `isavuconazole_groll_ddi`,
`posaconazole_dolton_ddi`, `rifampicine_loos_ddi`,
`ritonavir_kappelhoff_ddi`, and `tacrolimus_kamp_envarsus_ddi`.

Do not update those hashes without retraining: that would falsely relabel stale
artifacts. PMx Explain should display current availability from the generated ML
catalog, not infer it from the presence of an `.rds` file.

### ML scientific-readiness distinction (additional audit)

Do **not** interpret `availableInTdm=true` in
`ml/catalog-summary.generated.json` as a passed internal research
validation. The public `fec16cb` manifest has **12 technically available**
ML predictors, but only **3** meet both repeated-CV and untouched-holdout
thresholds: `cefepime_seo`, `dapto_dvorchik_IV_adults` and
`tacrolimus_woillard_ddi`. The remaining **9** can only be described as
*experimental and below the recorded performance thresholds*, not as
fully qualified predictors.

In particular, the retrained Goti XGBoost artifact now has a matching source
SHA, but its untouched-holdout relative RMSE is **19.6028%** (threshold
15%), P20 **73.2%** (threshold 80%), and
`validation.untouchedHoldout.passed=false`. Its alternate-PopPK
transportability test also fails. PMx must display model source alignment,
ML source compatibility, internal-prediction performance and external
clinical validation as **different fields**; never use an
`ARTICLE_ALIGNED` badge as a proxy for predictive ML validation.
Every model currently has `validationStatus=NONE` for clinical validation.

## Required PMx Explain import and checks

1. Resolve the existing uncommitted PMx synchronization work, then pin
   `integration/mipd-release.json` to the full `fec16cb...` SHA.
2. Import the exact 46-file model snapshot and the matching catalog, metadata,
   French details, English localization and source-review records.
3. Import generated traceability and ML summaries from the same SHA; do not
   combine files from different MIPD revisions.
4. Assert byte/hash parity for every imported C++ model and exact ID parity for
   all 46 catalog records.
5. Verify that Goti's source status and renal-boundary note propagate to every
   PMx model view without being promoted to external validation.
6. Verify that the 12 stale ML artifacts above are shown as unavailable rather
   than silently enabled. Separate *available/compatible* from *passed internal
   ML thresholds*: only 3/12 available models pass CV and holdout; Goti remains
   below threshold despite matching model hashes.
7. Run PMx transfer-contract tests for Simulation, Lego, TDM, DDI and PD,
   including explicit unit and unsupported-field reporting.
8. Run browser checks for the public MIPD link, release fingerprint, FR/EN model
   metadata and the public custom-model/Pmetrics restrictions.
9. Keep MIPD Engine executable code out of PMx Explain; PMx should retain only
   the documented source snapshot and integration contract it needs.

## Deployment automation follow-up

The `production` GitHub environment still requires secure manual entry of
`SHINYAPPS_ACCOUNT`, `SHINYAPPS_TOKEN` and `SHINYAPPS_SECRET`. The account value
is `tdmhub`; token and secret values must never be committed or printed. The
workflow must continue to deploy only the existing application and preserve the
safe production flags.

Any deployment after `fec16cb` requires a new authorization and a fresh public
acceptance record. A documentation-only commit after this release does not
change the currently deployed application SHA.
