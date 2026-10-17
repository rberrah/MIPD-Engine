# Security policy

MIPD Engine is research and educational software. It is **not** a medical device and must not be used to make autonomous patient-dosing decisions.

## Supported versions

Security fixes are targeted to the current `main` branch and to explicitly identified immutable releases. Older model and deployment snapshots are preserved for reproducibility, not guaranteed to receive backported fixes. See `docs/DEPLOYMENT.md` for the currently running public SHA.

## Report a vulnerability privately

Please use GitHub's **private vulnerability reporting / security advisory** facility for this repository, if enabled by its maintainers. Otherwise contact the repository owner through an established **private** organizational channel. Do not post a proof of exploit, keys, access tokens, patient data, clinical records, or private scientific articles in a public issue, pull request, or CI log.

A useful report describes the affected commit/version, reproducible **synthetic** input, threat model, impact, and mitigation. Do not test attacks on the live Shinyapps service or other people's sessions without explicit authorization.

If you inadvertently observe personal health information or secret material, do not copy or publish it. Ask the maintainer to arrange a secure response path.

## Security boundaries

- Production **must** set `ALLOW_CUSTOM_MODELS=false`: arbitrary uploaded C++ and host code are not an accepted public feature.
- The public Shinyapps runtime currently **must** set `ENABLE_PMETRICS=false` until its native dependencies are independently qualified on that host.
- Keep uploaded observations in session-scoped workspaces, validate formats and file sizes, and clean up local temporary files. A local cleanup test does **not** prove Shinyapps logs, host backups, or interrupted-process storage are nonpersistent.
- Never put Shinyapps access tokens into the repository, fixtures, GitHub issue comments or build logs. Use GitHub `production` environment/repository secrets and minimal token scope.
- Clinical or patient-level datasets must not be uploaded into this public repository or included in bug reports.
- Pin/verify package versions and source-model fingerprints for reproducibility; updates affecting numerical outputs require numerical requalification.

## Handling and disclosure

The maintainer evaluates reports, triages impact, fixes the issue on a controlled branch, and publishes a non-sensitive summary after remediation. Do not disclose a potential vulnerability publicly before coordinated remediation. There is no promised fixed response time or clinical service level agreement.
