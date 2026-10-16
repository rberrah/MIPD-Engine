# Known-truth parametric versus NPAG benchmark

Status: **synthetic methodological benchmark**, generated 2026-10-11 with
Pmetrics 3.2.6. It is not clinical validation and does not establish universal
superiority of either population representation.

## Design

- One-compartment IV model, dose 100 mg, fixed volume 30 L.
- Gaussian residual SD `0.02 mg/L + 0.10 × prediction`.
- Three clearance populations: log-normal unimodal, asymmetric mixture and
  bimodal mixture.
- 90 training and 180 held-out subjects per scenario, with samples at 1, 4 and
  8 h.
- Population fitting uses all three training samples. Held-out prediction uses
  the 1 h sample to predict the latent concentration and observation at 8 h.
- Parametric fit: one log-normal CL distribution estimated by marginal maximum
  likelihood on a fixed numerical integration grid.
- NPAG: Pmetrics `$fit(algorithm = "NPAG")`, two Sobol initializations (40/80
  points; seeds 23/71), 50 cycles maximum.

Primary metrics are conditional held-out predictive log-likelihood, absolute
bias, RMSE, MAE and 90% predictive coverage. Relative errors are archived but
are not primary because division by low late concentrations makes them unstable.

## Results

NPAG rows below are the first run; the second run gives the same rounded
predictive metrics while the final number of support points can differ.

| True population | Method | Mean predictive log-likelihood | Bias (mg/L) | RMSE (mg/L) | MAE (mg/L) | 90% coverage |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Log-normal | Parametric log-normal | -0.214 | 0.032 | 0.282 | 0.226 | 91.7% |
| Log-normal | NPAG | -0.201 | 0.008 | 0.281 | 0.223 | 90.0% |
| Asymmetric | Parametric log-normal | -0.526 | 0.051 | 0.361 | 0.282 | 83.9% |
| Asymmetric | NPAG | -0.416 | 0.027 | 0.357 | 0.280 | 86.7% |
| Bimodal | Parametric log-normal | -0.515 | 0.042 | 0.454 | 0.378 | 87.2% |
| Bimodal | NPAG | -0.121 | 0.031 | 0.455 | 0.361 | 87.8% |

The negative control matters: under a truly log-normal population, the two
approaches are practically indistinguishable on absolute RMSE and MAE. Under
the bimodal scenario, NPAG improves predictive log-likelihood and MAE but not
RMSE. Method choice therefore remains conditional on purpose and diagnostics.
Support-point count is not used as an estimate of biological subgroup count.

## Reproducibility

- Script: `tdm-engine/benchmarks/pmetrics_known_truth.R`
- Configuration: `tdm-engine/pmetrics/benchmarks/pmetrics-known-truth-config.json`
- Synthetic data: `tdm-engine/pmetrics/benchmarks/pmetrics-known-truth-data.csv`
- Full results: `pmetrics-known-truth-results.csv` and `.json`
- Stability summary: `pmetrics-known-truth-stability.csv`
- Command: `npm run benchmark:pmetrics`

## Official Pmetrics example

`tdm-engine/benchmarks/pmetrics_official_example.R` loads
`Examples/src/bimodal_ke.csv` directly from the installed Pmetrics package. It
does not copy or modify the public tutorial dataset. The parameter bounds,
assay-error polynomial, seed, analytical shortcut and ODE equation reproduce
the two official `bimodal_ke` example scripts shipped with Pmetrics 3.2.6.

The automated smoke run uses fewer initial points and cycles, then checks the
population means from the analytical and ODE formulations within an 8%
tolerance. The full run uses the official 2028 Sobol points, 1000 cycles and a
2% tolerance:

```text
npm run test:pmetrics-official
npm run benchmark:pmetrics-official
```

The generated `pmetrics-official-example.json` records the package version,
source checksum, settings, support-point diagnostics and numerical comparison.
This is a software reproduction check, not an external or clinical validation.

The experimental NPAG/NPOD smoke comparison is separate in
`pmetrics-npod-smoke.json`. It uses two seeds per algorithm and records final
objective, support-point count and runtime. It is deliberately not a
superiority test.
