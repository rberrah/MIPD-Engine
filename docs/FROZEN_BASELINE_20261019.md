# MIPD-Engine — frozen research/teaching baseline

**Freeze decision:** 19 October 2026  
**Intended use:** research and education only  
**Repository:** `rberrah/MIPD-Engine`

## Decision

MIPD-Engine is feature-frozen for its current research/teaching scope. The
repository is the authoritative source for the executable engine, 46-model
library, scientific metadata, qualification tests and experimental ML registry.

The merge containing this document is the maintenance baseline once protected
main qualification succeeds.

## Public runtime versus repository baseline

The accepted public Shiny deployment remains:

- source SHA `fec16cbdf8e33ed6c791226ccf236761333f8de8`;
- Shinyapps app ID `17776160`;
- bundle `12702795`;
- rollback bundle `12702610`;
- arbitrary custom C++ disabled;
- public Pmetrics disabled.

The qualified repository baseline immediately before this freeze PR is
`428d38f65244680e8bc6335afd125469e804d805`.

A direct comparison from `fec16cb` to `428d38f` contains **no changes** to
`app.R`, `R/`, `models/`, `runtime/` or `renv.lock`. The intervening
commits are governance, documentation and qualification assets. Consequently,
redeploying Shiny only to make the public commit SHA equal to repository HEAD
would add operational risk without changing the executable/model payload.

## Qualification at freeze

The protected repository currently establishes:

- 46/46 model compilation and software contracts;
- 46/46 frozen numerical software regression across 98 supported
  route/profile cases and 294 AUC/peak/C24 values;
- 13 experimental MAP prior/covariance numerical checks;
- strict DOI/source-provenance qualification;
- independent analytical/physical source checks for selected models;
- Pmetrics 3.2.6 qualification on its compatible runner;
- active main ruleset requiring PRs and `origin-parity`, `core-r`,
  `pmetrics` and `crossref-doi-audit`.

These are software/scientific qualification layers for a research tool. They
do not imply independent reproduction of every primary publication, external
clinical validation, clinical utility or medical-device status.

## Deferred work that does not block the frozen baseline

The following are explicitly deferred research or operational hardening:

- complete independent reproduction of every primary model publication;
- external clinical validation;
- additional review of documented model adaptations and selected covariate
  boundary conventions;
- further qualification/retraining of historical simulation-trained ML
  predictors;
- optional GitHub Actions Shinyapps credentials and automated deployment;
- hosting-provider log, backup and retention review.

These items must not be presented as missing implementation of the frozen
research/teaching engine.

## Maintenance policy

After freeze, changes should be limited to:

1. security or dependency fixes needed to keep the engine buildable;
2. factual/source corrections where the current implementation or metadata
   would otherwise be misleading;
3. reproducibility fixes for a broken qualification test;
4. explicitly requested model/research updates treated as a new maintenance or
   research release;
5. deployment fixes required to restore the already documented public service.

No new major model-library expansion, clinical-product workflow or ML feature
programme is part of this frozen baseline.
