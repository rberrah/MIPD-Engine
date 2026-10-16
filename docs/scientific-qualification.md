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

**First scientific tranche:** Roberts (2011) and Revilla (2010)
vancomycin models, 6 synthetic scenarios in
tests/source_equation_reproduction.R. Independent analytical R closed-form
solutions are compared with mrgsolve predictions. Published CL/V, covariates,
units and C(t) are tested. OMEGA/SIGMA, Bayesian fitting, external patients
and clinical dosing are NOT validated. Thus 2/46 PARTIAL_SOURCE_EQUATIONS;
0/46 FULL_SOURCE_REPRODUCTION; 44 pending. These counts must not be
misrepresented as 2 completely validated models.

The source-to-model evidence records are in
scientific-validation/source-evidence.json. Every CI run generates the
46-row report scientific-validation/coverage.generated.json with
python3 scripts/audit_scientific_evidence.py (also uploaded as an Actions
artifact). Claims with missing tests, mismatched DOIs or unjustified full
reproduction fail CI.

## Work packages and acceptance criteria

### P0 — Independent source reproduction of 13 MAP priors

Queue:
1. Vancomycin: vanco_roberts, vanco_pkjust, vanco_goti
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

### P1 — Review 35 models without confirmed source alignment

These comprise 25 models with documented adaptations and 10 that remain\nnot verified at the source-fidelity level. The original publication may lack executable code or describe components
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

## Existing independent analytical references

- Roberts, Antimicrob Agents Chemother (2011) 55:2704-2709,
  DOI 10.1128/AAC.01708-10, equations 3/4 and Table 2:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3101407/
- Revilla, Br J Clin Pharmacol (2010) 70:201-212,
  DOI 10.1111/j.1365-2125.2010.03679.x, final equations and Table 4:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2911550/

### Exit gate for initial PR

46 audit rows, precisely 2 partial source-equation records, exactly 0 full
source reproductions, and 6 independent numerical scenarios passing CI.
The pre-existing 8 regressions and model origin/hashes must not change.
