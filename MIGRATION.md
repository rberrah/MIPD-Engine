# Migration plan

## Phase 1 — frozen baseline
Copy the engine, 46-model library, provenance, tests, Pmetrics fixtures and experimental ML layer from PMx commit `2f97f109c3f4af257776a595f6546c9cd5f2321d`. Do not remove the PMx copy.

## Phase 2 — standalone qualification
A clean installation must compile all 46 models, reproduce the 8 current deterministic reference profiles, pass MAP/privacy/Lego/DDI/PD/oncology/infectiology tests, run Pmetrics independently on a compatible runner, verify ML artifact hashes, and deploy with arbitrary C++ disabled.

## Phase 3 — PMx cut-over
PMx becomes a consumer of the independently deployed engine through versioned contracts. Only then remove duplicated engine/library files from PMx, in a separate reversible PR.

## Evidence baseline
- 46/46 software tested
- 8/46 deterministic software-regression profiles
- 0/46 marked independently source-reproduced
- 13/46 currently eligible for the experimental MAP workflow
