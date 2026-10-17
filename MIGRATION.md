# Migration plan

## Phase 1 — frozen baseline — COMPLETE

The engine, 46-model library, provenance, tests, Pmetrics fixtures and full frozen experimental ML artifact set were copied from PMx commit `2f97f109c3f4af257776a595f6546c9cd5f2321d`. Origin parity is machine-checked by Git blob hashes. PMx Explain has since removed the executable `tdm-engine/` mirror, while retaining the `static/tdm/` browser model snapshot and pinning the external source to immutable SHA `83a3026df7c0881d2725333a90d81ef2ed218f90`; PMx remains a separately maintained consumer.

## Phase 2 — standalone public deployment — COMPLETE; scientific qualification ongoing

The **independent public app is deployed and accepted at immutable MIPD SHA
`fec16cbdf8e33ed6c791226ccf236761333f8de8`**, Shinyapps bundle
`12702795`, existing app ID `17776160`. The original target
`https://tdmhub.shinyapps.io/MIPD_Engine/` was kept, production arbitrary
C++ and Pmetrics were disabled, and public HTTP, FR/EN websocket,
two-session, synthetic PK, MAP and DDI acceptance passed. Prior bundle
`12702610` remains available for rollback. See
[release handoff](docs/pmx-sync-report-2026-10-17.md).

This release was performed through **authenticated local rsconnect**, not
GitHub Actions: the GitHub `production` secrets are still missing. The
earlier Actions attempt failed before deploying. A new code or documentation
commit on `main` **does not** imply the live version changed.

Standalone CI additionally covers 46/46 model compilation/contracts,
8 frozen deterministic software-reference profiles, three *partial*
independent primary-source equation suites, MAP/physical invariants, safe
Lego/DDI, local session cleanup, Pmetrics 3.2.6 on its compatible CI host,
DOI provenance and release manifest. It does **not** prove independent
complete primary-model reproduction, provider-side nonretention, external
patient validation or medical safety.

Outstanding operational/research work, **not a reason to revert the
accepted public release**:
- Configure the three GitHub Shinyapps secrets in `production` to enable
  future workflow-based releases;
- independently audit host logs/backup/retention;
- continue independent primary-source/MAP review and external validation;
- requalify 12 intentionally blocked, stale-hash ML artifacts;
- synchronize the independently versioned PMx consumer after its own CI.

## Phase 3 — PMx source-pin/live release alignment
PMx has already removed the executable R/Shiny `tdm-engine/` mirror and consumes an immutable standalone SHA. Its 46-model static snapshot is a browser UI resource, not a second executable model library. The remaining work is to replace the stale immutable pin with an exactly deployed, independently tested MIPD commit, explicitly refresh source-fidelity badges (Goti is now `ARTICLE_ALIGNED` at the public SHA after primary renal-boundary corrections, still **not** clinically validated), validate interface transfers and privacy, and preserve a reversible rollback to the last known working SHA. Do not delete the browser snapshot unless an equivalent versioned consumer mechanism exists.

## PMx browser snapshot drift observed on 2026-10-17

A fresh Git-blob comparison between standalone `models/*.cpp` and PMx
`static/tdm/models/*.cpp` shows **28 identical** and **18 differing** files
among the 46 active models. PMx additionally carries four retired/unlisted
static files (`genta_debord`, `genta_evers_icu`, `genta_evers_standard`,
`vanco_revilla`). That is expected for a pinned older consumer snapshot but
**blocks advancing** the source pin without coordinated migration.

Current differences: amiodarone_pollak_ddi, cefazoline_cheng,
cefazoline_lanoiselee, cefepime_an, cefepime_seo,
ciclosporine_press_ddi, cobicistat_mathias_ddi, diltiazem_hermann_ddi,
erythromycine_mather_ddi, fluconazole_debruyne_ddi,
isavuconazole_desai_ddi, isavuconazole_groll_ddi,
posaconazole_dolton_ddi, rifampicine_loos_ddi, rifampicine_marsot,
ritonavir_kappelhoff_ddi, tacrolimus_kamp_envarsus_ddi, vanco_goti.

The PMx snapshot derives `aligned/adapted/not verified` badges from older
static metadata. The newly deployed MIPD source re-reviewed and corrected
`vanco_goti` (CrCl cap and age >60 creatinine floor), now recording
`ARTICLE_ALIGNED` on *author-led article review*, with synthetic numerical
source/MAP checks but **zero external validations**. The next separate PMx
integration PR must update all 18 changed browser C++ files, regenerate the
catalogue, import Goti's current source-review/Details, and surface the 12
blocked stale-hash ML artifacts while retaining all 46 active models.
A green PMx CI against old `83a3026` does not qualify the newer MIPD SHA.

## Current evidence counts (source-level maturity)
- 46/46 software tested and 8/46 frozen deterministic reference profiles;
- 3/46 **partial** independent source-equation checks; 0/46 complete primary-model reproductions;
- 11 article-aligned, 33 documented adaptations and 2 source ambiguities in the raw review; source-fidelity UI separately treats missing full text as `not verified`;
- 13/46 experimental MAP priors; 0/46 externally clinically validated;
- 12 older ML artifacts unusable due to changed model-source SHA; the Goti ML was retrained but fails its preregistered internal performance thresholds.


## Test ownership

The following imported tests are **cross-repository PMx → MIPD contract tests**, not standalone engine tests:

- `advanced_lego_test.R`
- `lego_external_exports.R`
- `lego_mlxtran_exports.R`
- `lego_native_exports.R`
- `covariates_native_test.R`
- `input_labs_test.R`
- `teaching_lab_test.R`
- `workshop_pk_test.R`

They require fixtures generated by PMx Explain's Lego/translator/covariate/teaching-lab/workshop tooling. They remain versioned here as compatibility harnesses, but standalone CI does not fabricate their missing upstream fixtures. They must run from the PMx integration pipeline against a released MIPD Engine contract.

`paired_benchmark_test.R` is part of the standalone core suite now that the frozen binary ML artifact set is present.


## Frozen ML artifact migration

Completed on 2026-10-16 from the same frozen PMx commit.

- 92 RDS files copied;
- 5,762,570 bytes total;
- every target Git blob SHA matches the source Git blob SHA recorded in `migration/source-manifest.json`;
- the artifacts remain explicitly experimental and their migration does not alter their validation status.
