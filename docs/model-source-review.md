# Model source review

Review date: 2026-10-16.

This review compares each executable library model with the locally available primary full text. It is an author-led source audit, not an independent reproduction of the original analysis. The private article PDFs are not part of this repository.

## Coverage

- 46 active library models.
- 23 models checked against a locally available primary full text.
- 23 models still require their primary full text before an equation-level conclusion can be recorded.
- 0 models are claimed as independently source-reproduced.

The machine-readable coverage and per-model findings are in [`models/source-review.json`](../models/source-review.json).

## Corrections from this review

- `vanco_goti`: CL between-subject variance is now `log(1 + 0.398^2) = 0.147043196`, consistent with the published 39.8% CV and with the conversion already used for V1 and V2.
- `rifampicine_marsot`: the final residual model is additive only. The published 2.256 mg/L estimate is treated as an SD, so mrgsolve receives the variance `2.256^2 = 5.089536`; proportional variance is zero.
- `cefazoline_cheng`: no numerical change. The final published equation gives a Q body-weight exponent of 0.78, although the generic methods paragraph says 0.75 for clearance terms. The executable model keeps the final equation and now documents the inconsistency and its unbound-concentration output.

## Models awaiting full text

`amiodarone_pollak_ddi`, `amox_fournier`, `cefazoline_lanoiselee`, `cefazoline_schmitz`, `cefepime_an`, `cefepime_seo`, `ciclosporine_press_ddi`, `cobicistat_mathias_ddi`, `diltiazem_hermann_ddi`, `erythromycine_mather_ddi`, `everolimus_moes_ddi`, `everolimus_terheine_ddi`, `fluconazole_debruyne_ddi`, `isavuconazole_desai_ddi`, `isavuconazole_groll_ddi`, `posaconazole_dolton_ddi`, `rifampicine_jeremiah`, `rifampicine_loos_ddi`, `ritonavir_kappelhoff_ddi`, `sirolimus_djebli_ddi`, `tacrolimus_kamp_envarsus_ddi`, `tacrolimus_woillard_ddi`, `voriconazole_vandenborn_ddi`.

Until those texts are available, their DOI, citation and documented implementation status are checked, but their equations must not be described as full-text verified.
