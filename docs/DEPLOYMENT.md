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

### Current public release evidence

The public application was deployed from a clean local `main` checkout on
2026-10-17 after the complete standalone qualification workflow passed:

- source SHA: `f954b53204335b0b4ee6f30999f4a3e88c968115`;
- qualification run: [GitHub Actions 39021387658](https://github.com/rberrah/MIPD-Engine/actions/runs/39021387658);
- Shinyapps bundle `12702610`, task `1758914627`, image `15802838`;
- target: existing app ID `17776160` at
  `https://tdmhub.shinyapps.io/MIPD_Engine/`;
- immutable fingerprint check: HTTP 200, exact SHA, 46 models,
  `ALLOW_CUSTOM_MODELS=false`, `ENABLE_PMETRICS=false`;
- browser check: distinct FR/EN websocket sessions, with the second session
  remaining connected after the first closed.

The first browser probe encountered a transient non-200 response while the
Shinyapps instance was cold. A direct HTTP request then returned 200 and the
repeated browser probe passed. The probe now retries this bounded startup
condition but still fails after three unsuccessful responses.

This was a credentialed local deployment, **not** a GitHub Actions deployment
run. The `production` environment still needs the three documented Shinyapps
secrets before automated deployment can be claimed. The previously recorded
bundle `12698610` remains the rollback reference in the Shinyapps deployment
history. Provider-side logs, backups and retention were not inspected.

### Required operator action

The current GitHub app connection does not provide a workflow-dispatch
action or access to GitHub environment secrets. To publish after successful
CI, open the workflow link above, choose **Run workflow**, select **main**,
leave `include_pmetrics` unchecked, and approve the `production` environment
if protection is configured. The following secrets must be available in the
repository or the `production` environment: `SHINYAPPS_ACCOUNT=tdmhub`,
`SHINYAPPS_TOKEN`, `SHINYAPPS_SECRET`. Never paste the token or secret
into an issue, PR, log or chat.

A complete acceptance report must identify the GitHub Actions deployment
run, deployment commit and fingerprint equality; absence of a run or marker
means **not deployed/verified**, regardless of a green standalone CI.

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

On 2026-10-16, PMx Explain removed its executable `tdm-engine/` mirror and
uses `integration/mipd-release.json` pin
`83a3026df7c0881d2725333a90d81ef2ed218f90`. Its
`static/tdm/models/` remains a 46-model *browser snapshot*, not a second
engine. This pin is older than the current independent scientific-qualification
commits. PMx must **not** advance its immutable consumer pin or remove the
browser snapshot before deployment is proved by exact SHA, its model
snapshot/contract pass against the new SHA, and the new scientific status
of `vanco_goti` is propagated to UI. Changes to PMx must be handled by the
separate PMx owner/agent in a reversible, separately tested PR.

### Evidence still outstanding

The real-browser check proves only distinct active sessions and connection
lifecycle. It does **not** establish hosting-provider log retention policy,
absence of sensitive data in memory or `/tmp` after interrupted processes,
or external clinical safety/validity. Follow the separate privacy audit and
human inspection before stating those claims.
