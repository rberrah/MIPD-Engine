# MIPD Engine

Standalone research and educational engine for model-informed precision dosing,
extracted from **PMx Explain** and now maintained independently.

## Repository authority

This repository was initialized from the validated PMx source at:

- source repository: `rberrah/rberrah.github.io`;
- source commit: `2f97f109c3f4af257776a595f6546c9cd5f2321d`;
- extraction started: 2026-10-16.

Since 2026-10-17, **MIPD-Engine is the sole authoritative source** for the
executable engine, the 46-model library, model metadata, scientific
qualification and MIPD ML registry. PMx Explain is a separate educational
consumer and integration surface.

## Public release and consumers

The standalone Shiny app was deployed and accepted on 2026-10-17:

- **Deployed source SHA:** `fec16cbdf8e33ed6c791226ccf236761333f8de8`
  (tag `public-2026-10-17-fec16cb`);
- **Shinyapps target:** `tdmhub/MIPD_Engine`, app ID `17776160`,
  bundle `12702795`; rollback bundle `12702610`;
- **Public restrictions:** arbitrary custom C++ disabled and Pmetrics disabled
  on the current public host;
- **Status:** public software smoke acceptance passed; full article-model
  reproduction and external clinical validation are not claimed;
- **Deployment method:** authenticated local rsconnect. Automated Shinyapps
  deployment is not established because the GitHub `production` secrets are
  not configured.

PMx Explain currently pins the exact deployed SHA `fec16cb` through its
versioned consumer contract and independently checks the 46-model catalogue,
metadata, source-review information and ML availability. MIPD changes must be
qualified here first; PMx may then update only its consumer contract.

The frozen code/qualification baseline before this documentation-only freeze
update is `b7895d0aaeb54dd847ef39aa171444f6c60cb449`. It is five commits ahead of
the public source SHA; those five commits add governance, documentation and
qualification assets only. The executable runtime and model source files are
unchanged from `fec16cb`, so the public deployment remains correctly identified
by `fec16cb`.

See [deployment evidence](docs/DEPLOYMENT.md),
[PMx consumer handoff](docs/pmx-sync-report-2026-10-17.md) and
[scientific qualification](docs/scientific-qualification.md).

## Scope

The repository contains:

- the Shiny MIPD engine;
- the population PK model library and structured provenance;
- MAP estimation with `mapbayr`;
- optional nonparametric Pmetrics support when the runtime can build it;
- model averaging, exposure and future-dose simulations;
- DDI, PD, oncology and infectiology research workbenches used by the engine;
- software-regression and source-reproduction tests;
- experimental simulation-trained ML artifacts, clearly separated from
  clinically validated methods.

## Evidence levels

Model availability is not equivalent to scientific or clinical validation.
Current standalone qualification distinguishes software testing, frozen
software regression, partial independent source reproduction and external
validation.

At the qualified code baseline `b7895d0`:

- 46/46 models compile and satisfy their software contracts;
- 46/46 models are covered by frozen synthetic numerical software baselines
  across 98 route/profile cases, with 294 frozen AUC/peak/C24 values;
- Roberts, Revilla and Goti have partial independently transcribed
  source-equation checks; 0/46 models are claimed as fully independently
  reproduced from their complete primary publication;
- 13/46 models are eligible for the current experimental MAP workflow and have
  numerical prior/covariance qualification;
- 0/46 models are claimed as externally clinically validated.

These software baselines are regression evidence, not proof that every
published model has been independently reproduced.

## Safety and intended use

Research and education only. This software is not a medical device and must not
be used as an autonomous clinical dosing system.

Public deployments must keep arbitrary uploaded C++ disabled. Patient-shaped
inputs must remain session-scoped and non-persistent.

## Standalone qualification

Current standalone CI covers:

- deterministic regeneration of the 46-model catalogue and traceability;
- compilation and software-contract validation of all 46 models;
- frozen all-model numerical regression across 98 supported profiles;
- 13-model MAP prior numerical matrix;
- DOI and model/ML provenance consistency;
- partial independent primary-source equation oracles;
- safe Lego import and private-session cleanup;
- Pmetrics 3.2.6 patient-level backend on its compatible runner.

Cross-repository compatibility tests that consume PMx-generated fixtures remain
integration contracts; they do not make PMx a runtime dependency of MIPD.

The simulation-trained XGBoost layer remains experimental. Artifact
compatibility, CV/holdout performance, cross-model transportability and
external-patient performance must be interpreted separately from model
software qualification.

### Reproducible runtime

Core standalone CI is pinned to R 4.6.1 and the dated CRAN snapshot
`2026-10-16`. `renv.lock` freezes the public dependency graph, while
`runtime/reference-runtime.json` checks exact direct scientific package
versions. Pmetrics is separately pinned to 3.2.6. The repository does not claim
an offline package mirror or a frozen operating-system/native-compiler
toolchain.

### Scientific qualification programme

Scientific source fidelity is tracked separately from regression testing.
Independently transcribed equations currently cover vancomycin models Roberts
(2011), Revilla (2010) and Goti (2018), across synthetic analytical scenarios.
Goti also includes renal preprocessing, exposure and an independent synthetic
MAP oracle. Partial equation reproduction is not full model reproduction,
external validation or proof of clinical utility.

See [scientific qualification](docs/scientific-qualification.md) and
[scientific evidence records](scientific-validation/source-evidence.json).
