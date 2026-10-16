# Model source review

Review date: 2026-10-16.

This is an author-led comparison of every executable library model with its primary publication at the evidence level recorded in [`models/source-review.json`](../models/source-review.json). The historical `_ddi` suffix only records the model's origin in DDI Manager; it does not remove the model from the scientific audit. This review is not an independent reproduction of the original analyses.

## Coverage

- 46/46 active models have a recorded source comparison.
- 38 models were checked against a primary full text: 23 model entries from private local PDFs and 15 from open full texts.
- 8 older or paywalled sources currently have a limited comparison based on the primary abstract, sometimes cross-checked against a full secondary scientific source.
- 0 private PDFs are committed to Git.
- 0 models are claimed as independently source-reproduced.

The separate [`models/pdf-inventory.json`](../models/pdf-inventory.json) accounts for all 47 local PDFs: 22 unique primary articles linked to active models, 2 duplicate copies, 3 secondary reviews, 2 application-context documents and 18 background or noncatalog studies. One primary article can support more than one executable model, as for the two Rosario variants.

## Corrections in this audit

- Pollak amiodarone: replaced mismatched disposition values with published CL/F, Q/F, V1/F, V2/F and IIV.
- Press ciclosporine: restored the 76 kg reference, prednisolone effects and residual model; removed the unsupported CYP3A5 baseline effect.
- Dolton posaconazole: restored CL/F, V/F, ka, lag and the medication and clinical covariates.
- Kappelhoff ritonavir: restored lag, lopinavir effect, fixed effects and IIV; documented omitted IOV, covariance and mixture residual.
- Kamp Envarsus: restored weight allometry, correlated CL-V1 variability, KTR variability and lag class; removed the copied Woillard variability and unsupported CYP effect.
- Desai isavuconazole: added race, patient-status, BMI and reported variability terms.
- Loos rifampicin, Hermann diltiazem, Mathias cobicistat, Debruyne fluconazole and Mather erythromycin: corrected the available published summary values and labelled the entries as NCA/review reconstructions rather than PopPK article models.
- Lanoiselée cefazolin, An cefepime and Seo cefepime: corrected OMEGA/SIGMA transcription to the source software's squared-SD convention.

## Remaining evidence limit

Equation-level checking is still limited for eight sources whose full text or executable model specification was not available during this pass: Pollak amiodarone, Mathias cobicistat, Hermann diltiazem, Moes everolimus, Debruyne fluconazole, Groll isavuconazole, Loos rifampicin and Djebli sirolimus. These entries are simulation-only except where a separate full-text-verified model is used, and every unverified assumption is stated in code and metadata.

The generated [`docs/tdm-model-traceability.md`](tdm-model-traceability.md) lists the finding, evidence level, implementation parameters, covariates, random effects and engine limitations for each model.
