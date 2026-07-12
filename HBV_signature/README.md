# HBV prognostic gene signature

Code for building and validating a prognostic gene-expression signature for HBV-related
hepatocellular carcinoma — filtering genes, picking survival-associated ones, and testing how well
the signature separates good vs. poor prognosis patients. Prediction uses Nearest Template
Prediction (NTP) inside a leave-one-out cross-validation loop, and results are tied back to
Kaplan-Meier survival.

### Main flow

1. `1_variation.filter.revised.R` — filter genes by variance/expression before modeling.
2. `2_SurvivalGene.R` — select genes associated with survival (univariate Cox).
3. `3_loocv.nn_cox_GP.R` — LOOCV using NTP, with Cox / nearest-neighbor gene ranking.
4. `4_combine_poor_good_signatures.R` — combine the poor- and good-prognosis gene sets into one signature.

### Helpers

- `Multiple_variation_filtering_to_KM_ver2.R` — sweep filtering thresholds through to Kaplan-Meier.
- `NTP.R` — Nearest Template Prediction. Method credited to **Yujin Hoshida (Broad Institute)**; kept
  with its original attribution.
- `NTPtoKMsummary.R` / `NTPtoKMsummary_batch.R` — turn NTP calls into KM survival summaries (single / batched).
- `strat.split.stef.R` / `strat.split.stef_batch.R` — stratified train/test splitting. Original
  utility credited to **"Stef" (2007)**; kept as-is.
- `summary.R` — collects results across runs.

### Notes

- Inputs are GCT expression matrices + a clinical/survival table. Cluster paths in the scripts are
  from where I ran them; point them at your own data first.
