# MIPD Engine

Standalone research and educational engine for model-informed precision dosing, extracted from **PMx Explain**.

This repository was initialized from the validated PMx source at:

- source repository: `rberrah/rberrah.github.io`
- source commit: `2f97f109c3f4af257776a595f6546c9cd5f2321d`
- extraction started: 2026-10-16

## Scope

The repository is intended to contain:

- the Shiny MIPD engine;
- the population PK model library and its structured provenance;
- MAP estimation with `mapbayr`;
- optional nonparametric Pmetrics support when the runtime can build it;
- model averaging, exposure and future-dose simulations;
- DDI, PD, oncology and infectiology research workbenches used by the engine;
- software-regression and source-reproduction tests;
- experimental simulation-trained ML artifacts, clearly separated from clinically validated methods.

PMx Explain remains the educational portal and consumer of the deployed engine. The extraction must preserve the current public PMx URLs and contracts until an independently deployable release of this repository has been verified.

## Evidence levels

Model availability is not equivalent to scientific or clinical validation. The migration preserves explicit maturity fields distinguishing software testing, deterministic numerical regression, independent source reproduction and external/clinical validation.

At the extraction baseline:

- 46 models are catalogued;
- 46/46 are software tested;
- 8/46 have deterministic software-regression reference profiles;
- 0/46 are currently marked as independently source-reproduced;
- 13/46 are explicitly eligible for the current experimental MAP workflow.

## Safety and intended use

Research and education only. This software is not a medical device and must not be used as an autonomous clinical dosing system.

Public deployments must keep arbitrary uploaded C++ disabled. Patient-shaped inputs must remain session-scoped and non-persistent.

## Migration status

The frozen scientific engine and model-library baseline has been copied without removing the active PMx consumer. The standalone repository now contains the 46-model library, all 92 frozen experimental ML RDS artifacts, model provenance, Pmetrics fixtures and the engine test suite. Removal from PMx remains gated on an independently verified deployment and cross-repository compatibility checks.

See `MIGRATION.md` and `ORIGIN.md`.


## Standalone qualification status

The extraction is intentionally staged. PMx Explain remains the active consumer until this repository is independently deployable.

Current migration branch qualification:

- origin/hash parity for frozen source files;
- deterministic regeneration of the 46-model catalogue;
- compilation/contract validation of all 46 models;
- 8 deterministic software-regression reference profiles;
- safe Lego import contract and session cleanup;
- Pmetrics 3.2.6 patient-level backend on a compatible runner.

Tests that require fixtures produced by PMx Explain remain cross-repository integration contracts and are not reclassified as standalone tests.

The simulation-trained XGBoost layer is experimental. Its 92 binary RDS artifacts are now migrated bit-for-bit from the frozen source and remain explicitly distinct from source-reproduced or externally validated model evidence.

### Reproducible runtime

Core standalone CI is pinned to R 4.6.1 and the dated CRAN snapshot `2026-10-16`. `runtime/reference-runtime.json` records exact direct package versions (including mrgsolve, mapbayr and XGBoost), and CI fails closed on drift. Pmetrics is separately pinned to 3.2.6 on R 4.6.1. A complete transitive `renv.lock` / offline mirror remains a later reproducibility enhancement rather than an unverified claim.
