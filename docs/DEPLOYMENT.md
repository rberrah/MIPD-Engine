# Standalone deployment

The repository can deploy the qualified public MIPD Engine directly to Shinyapps.io without a PMx Explain checkout.

## Public profile

The production profile is deliberately conservative:

- `ALLOW_CUSTOM_MODELS=false`;
- `ENABLE_PMETRICS=false` on the current Shinyapps runtime;
- R 4.6.1;
- dated CRAN snapshot `2026-10-16`;
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
