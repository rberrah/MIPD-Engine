## Purpose and evidence

Describe the problem, affected source-model IDs, intended release target, and links to related issues/primary articles. Do **not** include patient data, private PDFs, keys or Shinyapps tokens.

## Risk classification

- [ ] Documentation/governance only; no executable or generated scientific data changes
- [ ] Scientific model / mrgsolve / mapbayr / Pmetrics / DDI / numerical changes
- [ ] ML artifact, training domain or provenance changes
- [ ] Security/privacy/runtime/deployment change

## Qualification checklist

- [ ] Relevant standalone MIPD and DOI workflows are green (or exceptions explained)
- [ ] Changes to model parameters, routes, covariates, OMEGA/SIGMA and residual error are referenced to primary article passages
- [ ] Independent numerical oracle or regression tests added/updated when results can change
- [ ] Corresponding catalogue, source-review, French/English metadata and traceability regenerated
- [ ] Changed source hashes reviewed against all dependent ML artifacts; stale ML stays blocked until retraining
- [ ] Software-regression vs primary-source reproduction vs external-patient validation described honestly
- [ ] Safe production profile maintained: no arbitrary public C++ and no unsupported Pmetrics on Shinyapps
- [ ] If deploying, operator approval and immutable target SHA, public fingerprint/browser tests, rollback identified
- [ ] Any required PMx consumer sync is documented but performed in a separately qualified PMx PR

## Rollback and release impact

State whether this PR should change the **running Shinyapps** application. Merging a Git commit does **not** itself publish it. A docs-only commit does **not** advance the last deployed SHA.
