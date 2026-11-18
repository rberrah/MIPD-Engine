# Two-regimen A/B comparison (MVP-inspired, independently implemented)

Conceptual inspiration: [MVPapp](https://github.com/stevechoy/MVPapp) and
Choy et al. (2025), [doi:10.1002/psp4.70023](https://doi.org/10.1002/psp4.70023).
No MVPapp source code, UI templates, copyrighted media or assets have been reused.
No dependencies or licences were changed.

## Scope

Inside the existing **Posologies** results page, a closed-by-default panel
allows a direct A/B comparison of **two previously ranked dose regimens**.
It uses the existing `fits`, model weights, published model compilation and
`simulate_averaged_regimen()` machinery, with no refitting or reranking.

- The first 20 **distinct**, numerically valid ranked regimens are eligible.
- For predictable interactive costs, this optional display supports intervals
  between 0 and 336 h. Existing recommendation calculations remain unchanged.
- A and B are drawn on the same 24-h-to-longest-interval time horizon.
- The metric table reports AUC0–24, Cmin/Cmax **over each regimen's own
  interval**, and %time above the current concentration threshold when defined.
- Export session-scoped numerical curves and metric summaries as CSV.
- Both selected regimens use the **same** individual model fits and model
  averaging weights. There is no new propensity adjustment, model selection,
  population variability calculation or target attainment inference.

## Safety and interpretation

This display is descriptive and research-only. It must never update inputs,
the candidate ranking, MAP estimates or dose recommendations. Cmin and Cmax
may not have identical observation windows if interval A differs from B.
The 0.25-h drawing grid is deliberately not the same as the main recommendation
grid: small differences in AUC or Cmax due to simulation resolution are possible.
Concentrations use each model's existing units; do not apply an additional
unit conversion to the model output without an independently verified contract.

No data are persisted by this feature. CSV files are explicit, local downloads.
In public use, the existing model source/route/privacy gates continue to apply.

## Verification

- `Rscript tests/regimen_ab_test.R`: candidate filtering, duplication,
  invalid selections, deterministic Roberts population/MAP comparison,
  expected proportionality for dose changes, unchanged recommendation table.
- Existing `.github/workflows/ci.yml` runs the new R test as part of the
  `tests/*.R` regression suite.
- Before merging: perform EN/FR Shiny UI review, verify profiles of distinct
  intervals and long half-life drugs, verify session isolation and real public
  deployment only after a separate release qualification.

## Follow-up backlog (separate PRs)

1. Model parameter perturbation/sensitivity: define a **small source-reviewed
   whitelist** per model (not raw unrestricted mrgsolve `$PARAM`) before OAT
   tornado analyses; test linear and nonlinear models, correlated ETA handling.
2. Securely exported reproducible R runs: require exact source SHA, package
   versions, scenario inputs and an explicit data-exclusion policy.
3. Synthetic/virtual populations: validate reference database applicability,
   IIV/IOV coupling, covariance and sampling limits before any public claim.
4. PMx educational companion: an original OAT panel in the existing four
   teaching laboratories (no new site section).

A/B visualization is a lower-risk, useful first slice. It should not become
an independent application or a new source of clinical model authority.
