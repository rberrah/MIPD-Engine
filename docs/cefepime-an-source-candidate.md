# Cefepime An (2023) — independent source-equation candidate (2026-10-18)

**Primary source:** [An et al., *J Antimicrob Chemother* 78:1460–1470](https://pmc.ncbi.nlm.nih.gov/articles/PMC10474939/), DOI [10.1093/jac/dkad106](https://doi.org/10.1093/jac/dkad106). Methods Eq. (2) specifies renal/nonrenal/CRRT clearance and Table 2 provides numerical population and variability estimates.

**New independent candidate test:** `tests/cefepime_an_source_oracle.R`. It transcribes the paper's fixed effects independently, checks published variance scales, and compares six *synthetic* patients' IV concentrations against an analytical two-compartment zero-order-infusion calculation. It includes renal impairment, augmented clearance, women, CRRT and infusions lasting 0.5 h, 3 h or 24 h.

**Unresolved source issue:** The article's main HTML describes creatinine clearance using lean body weight, but it does not independently establish whether the source estimator is *Boer*, nor a minimum LBW clamp. The current C++ uses Boer and a floor at 30 kg; the test stays above this threshold and labels the conversion an *implementation convention*. Before promoting this model's evidence to PARTIAL_SOURCE_EQUATIONS, verify the exact lean-body-weight calculation with the original supplementary code or authors, and independently review how the IIV/residual model was parameterized. The existence of this test is not itself evidence of full population-model reproduction.

The three previously registered models (Roberts, Revilla, Goti) remain the only **adjudicated** partial source reproductions in `scientific-validation/source-evidence.json` until an independent source review is completed. No external-patient validation is claimed.
