# Standalone deployment

The repository can deploy the qualified public MIPD Engine directly to Shinyapps.io without a PMx Explain checkout.

## Public profile

The production profile is deliberately conservative:

- `ALLOW_CUSTOM_MODELS=false`;
- `ENABLE_PMETRICS=false` on the current Shinyapps runtime;
- R 4.6.1;
- dated CRAN snapshot `2026-10-16`;
- complete public dependency graph restored from `renv.lock`;
- exact direct scientific-package versions enforced by `runtime/reference-runtime.json`;
- all 46 model contracts and deterministic reference cases rechecked immediately before deployment.

Pmetrics 3.2.6 remains independently tested in CI. It must not be enabled on Shinyapps until that host provides a compatible runtime; a separate compatible/containerized deployment can expose it later.

## Required GitHub production secrets

Configure these repository or environment secrets before using the manual deployment workflow:

- `SHINYAPPS_ACCOUNT`
- `SHINYAPPS_TOKEN`
- `SHINYAPPS_SECRET`

The workflow uses GitHub environment `production`, so environment protection/review can be added without changing code.

## Manual local deployment

Configure the account locally with `rsconnect::setAccountInfo()`, then inspect the manifest:

```bash
Rscript scripts/deploy_shiny.R
```

Run a credential-free dry-run of the **exact production bundle** (including a
temporary immutable Git SHA fingerprint, then safely removing it):

```bash
Rscript scripts/deploy_shiny.R --verify-bundle
```

Deploy only after the manifest and bundle are reviewed:

```bash
ALLOW_CUSTOM_MODELS=false ENABLE_PMETRICS=false \
Rscript scripts/deploy_shiny.R --deploy
```

The production target is deliberately fixed to the existing application:

- account: `tdmhub`;
- app name: `MIPD_Engine`;
- appId: `17776160`;
- public URL: `https://tdmhub.shinyapps.io/MIPD_Engine/`.

The deployment script refuses a different account and supplies the existing appId explicitly, preventing an accidental parallel application from being created during repository cut-over.

## Versioned public proof and post-deploy acceptance (new)

A successful compilation or CI run is **not** evidence that the current public
Shiny instance is running that commit. The manual
[`Deploy standalone MIPD Engine to Shinyapps`](https://github.com/rberrah/MIPD-Engine/actions/workflows/deploy-shiny.yml)
workflow deploys **only `main`** to the existing application; it must use
`include_pmetrics=false`.

At deployment time `scripts/deploy_shiny.R --deploy` generates
`www/deployment-release.json` from `git rev-parse HEAD`, packages it with the
app and deletes the local temporary marker afterwards. The public URL is
`https://tdmhub.shinyapps.io/MIPD_Engine/deployment-release.json`. The
release contains **no credential or patient data**. It declares the exact
40-character Git SHA, app ID, model count, and the restricted production
profile. Startup code in `R/deployment_profile.R` enforces both restrictions
even if hosting environment variables differ; invalid fingerprints abort
startup.

After upload, deployment must pass BOTH:

1. `python3 scripts/check_live_deployment.py --expected-sha <FULL_MAIN_SHA>`
   verifies HTTP 200, the live Shiny HTML, **the actual deployed commit**,
   model count 46, `ALLOW_CUSTOM_MODELS=false` and Pmetrics disabled.
2. `python3 scripts/check_live_sessions.py` opens two independent FR/EN
   Chromium browser contexts, waits for actual Shiny websocket connections,
   checks distinct session IDs, and verifies one remains connected after
   closing the other. It sends no patient or custom C++ data.

Both checks are mandatory *after* `rsconnect::deployApp` in the production
workflow. An upload that succeeded but fails either test is **not a validated
public release**. A previously available public bundle with no
`deployment-release.json` cannot be identified as the new standalone release.

### Current public release evidence — 2026-10-17

**Current public release:** MIPD Engine commit
`fec16cbdf8e33ed6c791226ccf236761333f8de8`, immutable tag
`public-2026-10-17-fec16cb`, Shinyapps bundle **12702795**, task
`1758931558`, image `15803034`, existing app ID `17776160`.
The authoritative receipt and test evidence are in
[PMx sync/release handoff](pmx-sync-report-2026-10-17.md), which is explicitly
versioned to the **deployed** SHA. A local authenticated deployment was
performed after successful standalone qualification. This was **not** a
successful GitHub Actions deployment run.

The recorded post-release acceptance passed HTTP 200, exact SHA at
`/deployment-release.json`, 46 bilingual model cards, distinct FR/EN
websocket-connected sessions, synthetic PK, MAP and DDI simulation,
prohibition of arbitrary custom C++, and disabled Pmetrics. The generated
DDI comparison reported AUC24 283.0 vs 560.0 (ratio 1.98). These are
software checks with synthetic inputs, **not** external clinical validation.
We do not independently attest to Shinyapps provider logs, backups or data
retention.

**Rollback:** prior public bundle `12702610` (source SHA
`f954b53204335b0b4ee6f30999f4a3e88c968115`) is retained.
Do not remove it without a replacement rollback plan.

**Do not confuse `main` with production.** Subsequent commit
`9645ff4c8d743a350730f31f9aaab0605f4943e5` records documentation
and did **not** change the deployed app. Additional documentation commits
also do not advance the public fingerprint. All read-only live verification
against this release must specify
`--expected-sha fec16cbdf8e33ed6c791226ccf236761333f8de8` or the
equivalent full SHA workflow input; omitting it compares against newer `main`
and can produce a false *stale deployment* alert.

**Automation not yet configured:** the first manual Actions deployment
failed on the step `Configure Shinyapps account` because GitHub
`production` credentials were not available. Subsequent deployment was
completed using an authenticated **local** rsconnect session. This does
not create GitHub secrets. A future **GitHub Actions** deployment requires
the operator to configure `SHINYAPPS_ACCOUNT`, `SHINYAPPS_TOKEN` and
`SHINYAPPS_SECRET` as repository or `production` environment secrets.
No token value should be included in Git or a report.

### Required operator action

The current GitHub app connection does not provide a workflow-dispatch
action or access to GitHub environment secrets. To publish after successful
CI, open the workflow link above, choose **Run workflow**, select **main**,
leave `include_pmetrics` unchecked, and approve the `production` environment
if protection is configured. The following secrets must be available in the
repository or the `production` environment: `SHINYAPPS_ACCOUNT=tdmhub`,
`SHINYAPPS_TOKEN`, `SHINYAPPS_SECRET`. Never paste the token or secret
into an issue, PR, log or chat.

A complete acceptance report must identify the **deployment method**
(GitHub Actions or authenticated local rsconnect), source SHA, Shinyapps
bundle/task, public fingerprint equality and live browser test evidence.
A GitHub Actions run is required to claim **Actions-based** publication,
but is not required for a verified **local** deployment. Green CI alone
never proves a running release.

### Independent read-only re-verification and release artifact

A second GitHub Actions workflow,
[`Reverify public MIPD Shiny release`](https://github.com/rberrah/MIPD-Engine/actions/workflows/verify-public-shiny.yml),
can be run manually on `main` **without Shinyapps credentials**. Leave
`expected_sha` empty to compare the public application to the current
`main` commit, or provide an existing full 40-character released SHA.

It automatically runs again after a **successful** production deployment
workflow. A passing run proves public fingerprint/URL status and two real
websocket-connected, distinct Shiny sessions. It uploads
`public-mipd-live-verification` as a non-sensitive JSON artifact with the
deployed SHA, test date, public URL and proof-run URL.

**Important:** a release-verification run that is skipped, absent or fails
is not evidence of a successful production deployment. The independent
artifact is intended for the PMx consumer integration gate, alongside
source-catalogue parity and browser tests. No workflow creates secret
configuration, deploys without authorization or examines private user data.

The deployment script also rejects locally modified **tracked** Git files
when creating a release fingerprint; a dirty checkout must be committed
and qualified before publishing. Untracked files are excluded from the
`git ls-files` deployment manifest.

### PMx Explain compatibility gate

PMx Explain removed the executable `tdm-engine/` mirror and retains
`static/tdm/models/` as a browser-facing, versioned source snapshot.
At the latest inspected committed PMx revision on 2026-10-17,
`integration/mipd-release.json` still pins old MIPD SHA
`83a3026df7c0881d2725333a90d81ef2ed218f90`. There are **18 C++ source
differences** between that pin and the qualified public MIPD SHA `fec16cb`.
The Goti model is now `ARTICLE_ALIGNED` after independent synthetic renal
boundary tests (not externally validated), and **12 other ML artifacts**
are intentionally blocked because their training source hashes no longer
match. PMx must take C++ models, catalog, localized source-review fields,
traceability and ML availability **from the same published SHA**.
The full handoff is [docs/pmx-sync-report-2026-10-17.md](pmx-sync-report-2026-10-17.md).

Advance the PMx immutable pin only via a separate reversible PR after all
46 model source/hash parity tests and browser/Lego/TDM/DDI/PD consumer
contracts succeed. The PMx owner must reconcile uncommitted work before
replacing snapshots. No PMx source files should be changed from this repo.

### Evidence still outstanding

The real-browser check proves only distinct active sessions and connection
lifecycle. It does **not** establish hosting-provider log retention policy,
absence of sensitive data in memory or `/tmp` after interrupted processes,
or external clinical safety/validity. Follow the separate privacy audit and
human inspection before stating those claims.
