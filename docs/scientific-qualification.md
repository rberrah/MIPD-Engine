# Scientific qualification: MIPD-Engine

Baseline (2026-10-16): 46 shipped models pass compilation/contract tests;
8 deterministic software regression profiles; 13 experimental MAP priors;
0 full independent source reproductions and 0 external validations claimed.

## Evidence levels (not interchangeable)

| Level | What it demonstrates | What it does not demonstrate |
|---|---|---|
| SOFTWARE_TESTED | Compiles and executes according to its C++ implementation | Agreement with original paper |
| SOFTWARE_REGRESSION | Does not drift from historical predictions | Independent reproduction |
| PARTIAL_SOURCE_EQUATIONS | Named source equations and published constants independently reproduce fixed-effect trajectory/covariate scenarios | Complete stochastic/population model or clinical performance |
| FULL_SOURCE_REPRODUCTION | Complete source model, covariance/error structure, all supported routes and documented features reproduced and independently reviewed | Generalization to a new patient cohort |
| EXTERNAL_VALIDATION | Out-of-development patient-level predictive assessment with preregistered metrics and confidence intervals | Prospective clinical utility |
| CLINICAL_UTILITY | Prospective/clinical impact and safety as appropriate | Regulatory clearance by itself |

**Independent-equation tranche:** Roberts (2011) and Revilla (2010)
vancomycin models, 6 synthetic scenarios in
tests/source_equation_reproduction.R. Independent analytical R closed-form
solutions are compared with mrgsolve predictions. Published CL/V, covariates,
units and C(t) are tested. A third model, Goti (2018), has a separate
independent analytical IV two-compartment test with 8 synthetic scenarios in
tests/goti_source_equations.R. It covers the primary-source age >60
serum-creatinine floor, CrCl cap, fixed parameters, OMEGA/SIGMA, C(t) and AUC.
Thus 3/46 PARTIAL_SOURCE_EQUATIONS; 0/46
FULL_SOURCE_REPRODUCTION; 43 pending. These counts must not be
misrepresented as completely validated models.

### MAP prior and posterior numerical qualification - Roberts, Revilla and Goti

`tests/independent_map_prior_oracles.R` independently calculates concentrations,
Gaussian diagonal prior penalties and mixed/additive residual likelihoods.
It solves two-ETA one-compartment objectives for Roberts and Revilla, and a
three-ETA two-compartment objective for Goti, with base R. It compares each
optimum with `mapbayr` for one synthetic profile per model and checks
`OMEGA`/`SIGMA` against the three primary publications and explicitly prints
the small differences between reported CV-derived log-variances and historical
stored variances where applicable. The Goti case invokes both renal
preprocessing rules. CI acceptance is bounded and reproducible.

This qualifies the **software MAP optimization contract for those synthetic
input scenarios**, not the clinical predictive properties of the original
articles or individual posteriors fitted on the study cohorts. Coverage
remains **3/46 partial**, **0/46 full** and **0 external validations**.

The source-to-model evidence records are in
scientific-validation/source-evidence.json. Every CI run generates the
46-row report scientific-validation/coverage.generated.json with
python3 scripts/audit_scientific_evidence.py (also uploaded as an Actions
artifact). Claims with missing tests, mismatched DOIs or unjustified full
reproduction fail CI.

## Work packages and acceptance criteria

### P0 — Independent source reproduction of 13 MAP priors

Queue:
1. Vancomycin: vanco_roberts, vanco_pkjust, vanco_goti (initial synthetic MAP oracles complete)
2. Beta lactams: amox_carlier, cefazoline_lanoiselee, cefepime_an,
   cefepime_jonckheere, cefepime_seo
3. Remaining: dapto_dvorchik_IV_adults, levo_canoui,
   rifampicine_jing, rifampicine_marsot, tacrolimus_woillard_ddi

For each model, obtain and archive a permitted source reference; record
exact paper page/table/supplement and DOI, provenance, covariate unit and
population bounds, dose/formulation, and all fixed parameters. Separately
implement reference equations (analytical solver, independently written
ODEs, authors' model, or published numeric output). Include several within
source-domain cases per dosing route with representative and boundary
covariates. Pre-specify both absolute and relative tolerances. Compare
CL, V, Q, ka/Ktr as appropriate and concentration-time curves, AUC, Cmax
and tmax. Confirm random effect covariance, residual error and IOV with
separate tests. Reproduce MAP posterior parameters and predictive exposure
using an independent solver/reference code. An independent reviewer must
approve source agreement before the model is promoted to FULL.

### P1 — 46-model numeric regression, without falsely claiming source parity

Extend frozen deterministic software-regression cases from 8 to all 46,
with route, infusion mode, dose, synthetic covariates and two distinct
clinically plausible parameter/covariate scenarios per model. These detect
future drift only. Include independent physical/analytic invariants when
possible: nonnegative mass and concentration, correct units, mass
conservation, AUC = F*dose/CL for valid linear systems, expected changes
in exposure with clearance. Never generate source reproduction oracles
from the software under examination.

### P1 — Carlier 2013 independent structural-invariant tests (not source reproduction)

`tests/amox_carlier_analytical_invariants.R` independently solves the
linear two-compartment intravenous infusion equations for **three**
renal-function/infusion-duration scenarios. For the frozen code it tests
CL/Vc/Vp/Q, concentration-time curves, dose linearity and the physical
identity AUC(0–infinity) = IV dose / CL by dense numerical integration.
These are **mechanistic software/physical regression safeguards**.

The primary DOI is [10.1093/jac/dkt240](https://doi.org/10.1093/jac/dkt240).
Because the published model's complete fixed/random/residual parameterization
has not yet been independently reproduced in this tranche, this test
**does not** promote `amox_carlier` to partial or full primary-source
reproduction. The corresponding source-equation, IIV/SIGMA and MAP tests
remain in G1 (#6); real external patient validation remains in G3 (#8).

### P1 — Review 35 models without confirmed source alignment

These comprise 25 models with documented adaptations and 10 that remain
not verified at the source-fidelity level. The original publication may lack executable code or describe components
not supported by the current engine (IOV, correlations, transit, NPAG
nonparametric priors, DDI simplifications). Correct mismatches where source
permits. Otherwise provide named, dated discrepancy and retain a
DOCUMENTED_ADAPTATION / source-unverifiable outcome. A verified adaptation
is NOT interchangeable with the published original model.

### P2 — Independent external predictive validation

For each intended patient group and MAP sampling design, use a completely
held-out cohort (source centre/study split where possible), not synthetic
patients simulated by the same model. Predefine target and assessment
population, assay uncertainty, withheld observations, uncertainty intervals
and acceptance thresholds. Include bias, RMSE/MAE, coverage/calibration,
patient-bootstrap CIs, subgroup failures, dose-exposure prediction and
comparison with published alternatives. Avoid claiming clinical utility
from a simulated cohort or internal bootstrap.

### P2 — Publication-grade evidence and integration

For every decision archive permanent source citation, independent scripts,
reference input-output fixtures, operating-system/package/environment
versions, frozen commit hash, CI logs, discrepancies, evidence scope and
reviewer signature. Promote statuses only via reviewed PRs. PMx Explain
should consume model metadata and scientific verdicts from MIPD Engine,
never maintain its own duplicate clinical-science claims. Public deployment
remains for research and education; no autonomous clinical dosing.

## Goti model renal-boundary resolution

The primary article was re-reviewed locally. Its Methods specify bounding
Cockcroft-Gault CrCl at 150 mL/min and setting serum creatinine to 1 mg/dL
when it is lower than 1 mg/dL in patients **older than 60 years**. This
primary threshold supersedes the >65 wording in the 2023 secondary
reconstruction. `models/vanco_goti.cpp` now implements both rules, uses the
SI Cockcroft-Gault coefficients and encodes the Table 2 IIV values as exact
log-normal variances. Eight independent source-equation cases and one
independent three-ETA MAP case cover the boundaries. The source-fidelity badge
is `aligned`; this remains partial source reproduction, not external validation.

The dependent Goti XGBoost and DALEX artifacts were retrained from 1,000
synthetic profiles after the source correction. On the untouched internal
holdout, C0 plus one-hour post-infusion rRMSE was 19.60% with 73.2% within
+/-20%. Transport to synthetic profiles generated by Revilla was poor
(rRMSE 203.76%, bias 134.23%), so the artifact remains `experimental` and is
not used for the MAP-BE recommendation. The paired benchmark on 300 untouched
Goti profiles also favored MAP over XGBoost (rRMSE 4.31% versus 4.97% with two
samples). Machine-readable provenance and before/after metrics are stored in
`ml/validation/goti-source-correction-2026-10-17.json` and
`ml/validation/goti-paired-source-correction-2026-10-17.json`.

## Existing independent analytical references

- Roberts, Antimicrob Agents Chemother (2011) 55:2704-2709,
  DOI 10.1128/AAC.01708-10, equations 3/4 and Table 2:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3101407/
- Revilla, Br J Clin Pharmacol (2010) 70:201-212,
  DOI 10.1111/j.1365-2125.2010.03679.x, final equations and Table 4:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2911550/
- Goti, Ther Drug Monit (2018) 40:212-221,
  DOI 10.1097/FTD.0000000000000490; independent parameter cross-check:
  Antibiotics (2023) 12:301, Table 2, DOI 10.3390/antibiotics12020301.

### Exit gate for initial PR

46 audit rows, precisely 3 partial source-equation records, exactly 0 full
source reproductions, and 14 independent source-equation scenarios passing CI.
The Goti model hash and dependent ML provenance must change together; other
model origins and hashes remain unchanged.
