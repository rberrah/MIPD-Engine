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

Deploy only after the manifest is reviewed:

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
