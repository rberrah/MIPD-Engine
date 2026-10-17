# MIPD Engine

Standalone research and educational engine for model-informed precision dosing, extracted from **PMx Explain**.

This repository was initialized from the validated PMx source at:

- source repository: `rberrah/rberrah.github.io`
- source commit: `2f97f109c3f4af257776a595f6546c9cd5f2321d`
- extraction started: 2026-10-16

## Public release and consumers

The standalone Shiny app was **deployed and accepted** on 2026-10-17:
[https://tdmhub.shinyapps.io/MIPD_Engine/](https://tdmhub.shinyapps.io/MIPD_Engine/).

- **Deployed source SHA:** `fec16cbdf8e33ed6c791226ccf236761333f8de8` (tag `public-2026-10-17-fec16cb`).
- **Shinyapps:** `tdmhub/MIPD_Engine`, existing app ID `17776160`, bundle `12702795`; rollback bundle `12702610`.
- **Public restrictions:** no arbitrary custom C++ and no Pmetrics on the current host.
- **Status:** software/public smoke acceptance passed; full article-model reproduction and external clinical validation **not claimed**.
- **Deployment method:** authenticated local rsconnect. GitHub `production` environment secrets remain unconfigured, so automated Shinyapps deployment is **not** established.
- **PMx Explain:** separate consumer, still requires a reversible, versioned sync from its historic source pin. Do not infer PMx synchronization from MIPD deployment.

See the [deployment evidence](docs/DEPLOYMENT.md),
[release-to-PMx handoff](docs/pmx-sync-report-2026-10-17.md)
and [scientific qualification](docs/scientific-qualification.md).
A later documentation commit on `main` does **not** advance the deployed source SHA.

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

PMx Explain remains the separately managed educational portal. Its consumer pin and public model snapshot require explicit synchronization and independent CI validation, even after successful standalone Shinyapps publication.

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

The standalone repository contains the 46-model library, 92 experimental ML/DALEX RDS artifacts, model provenance, Pmetrics fixtures and the engine tests. PMx no longer contains an executable engine mirror, but retains an older browser-only snapshot. Its 46-model source pin must be updated in a **separate** qualified PMx PR. Twelve historical ML artifacts are intentionally hash-blocked following model fixes, and the retrained Goti predictor remains below internal research-performance thresholds.

See `MIGRATION.md` and `ORIGIN.md`.

The current article-to-code review coverage is recorded in [`docs/model-source-review.md`](docs/model-source-review.md) and [`models/source-review.json`](models/source-review.json). It distinguishes full-text review from software testing and independent numerical reproduction.


## Standalone qualification status

The standalone public deployment is complete; PMx consumer synchronization and independent scientific qualification remain separate, incomplete workstreams.

Current standalone GitHub qualification:

- origin/hash parity for frozen source files;
- deterministic regeneration of the 46-model catalogue;
- compilation/contract validation of all 46 models;
- 8 deterministic software-regression reference profiles;
- safe Lego import contract and session cleanup;
- Pmetrics 3.2.6 patient-level backend on a compatible runner.

Tests that require fixtures produced by PMx Explain remain cross-repository integration contracts and are not reclassified as standalone tests.

The simulation-trained XGBoost layer remains experimental. Historical artifact snapshots are preserved; the Goti XGBoost/DALEX files were retrained after the primary-source renal-boundary correction. A matching model hash means compatibility, **not** that CV, holdout, cross-model transportability or external-patient performance thresholds have passed. Research-level results must be read individually from the ML registry.

### Reproducible runtime

Core standalone CI is pinned to R 4.6.1 and the dated CRAN snapshot `2026-10-16`. `renv.lock` freezes the complete public dependency graph, while `runtime/reference-runtime.json` independently checks exact direct scientific package versions (including mrgsolve, mapbayr and XGBoost). Pmetrics is excluded from the public lock and separately pinned to 3.2.6 on R 4.6.1. The repository does not claim an offline package mirror or a frozen operating-system toolchain.

### Scientific qualification programme

Scientific source fidelity is tracked separately from regression testing. Independently transcribed equations cover vancomycin models Roberts (2011), Revilla (2010), and Goti (2018), across 14 synthetic analytical scenarios. Goti includes primary-source renal preprocessing, concentrations, AUC, OMEGA/SIGMA and an independent synthetic MAP oracle. **Partial equation reproduction is not full model reproduction, external validation or proof of clinical utility.** See [scientific qualification](docs/scientific-qualification.md) and [evidence records](scientific-validation/source-evidence.json). CI exports a 46-model coverage report and enforces model/ML provenance hashes.
