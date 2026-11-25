# Cognitive control and cortical oscillatory dynamics in athletes with ACL reconstruction

Statistical analysis code for:

> Jiménez-Martínez, J., Gutiérrez-Capote, A., Alarcón, F., & Cárdenas, D. (2026). Altered cognitive control and cortical oscillatory dynamics in athletes with ACL reconstruction during decision-making under uncertainty. *The American Journal of Sports Medicine*. https://doi.org/10.1177/03635465261494772

## What this repository contains

`R/stii_eeg_lmm_analysis.R` reproduces the statistical analysis reported in the article, structured in three levels:

1. **Behavioural performance** on a probabilistic go/no-go task (RTm, PART, IFRT, Accuracy, CER, PIF): linear mixed-effects models `DV ~ Group * Probability + (1 | Subject)`, Type III F-tests, partial eta-squared with 95% CI, Bonferroni-adjusted post hoc contrasts and Cohen's d.
2. **Cortical oscillatory activity** (midfrontal theta, left-central beta, posterior alpha; stimulus- and response-locked): `Power ~ Group * Probability * Congruence + (1 | Subject)`, with the same inference and effect sizes, on correct trials.
3. **Brain–behaviour relationships**: Pearson correlations between subject-level oscillatory power and behaviour, overall and within each group, with 95% CI and Fisher's z test for group differences.

EEG preprocessing and time–frequency decomposition (EEGLAB and FieldTrip in MATLAB, surface Laplacian, Morlet wavelets) were performed as described in the Methods section of the article; this repository covers the statistical analysis of the extracted power values.

## Requirements

R (≥ 4.3) with the packages:

```r
install.packages(c("readxl", "dplyr", "tidyr", "lme4", "lmerTest", "emmeans", "effectsize"))
```

## How to run

1. Place the data file `ST_EEG_ACL.xlsx` in `data/` (see `data/README.md` for the expected structure).
2. Open the repository folder as the working directory (e.g., as an RStudio project) and run:

```r
source("R/stii_eeg_lmm_analysis.R")
```

## Data availability

Participant-level data are not included because they contain health-related information protected under the EU General Data Protection Regulation. The study was approved by the Research Ethics Committee of the University of Granada (3110/CEIH/2022). Data may be requested from the corresponding author for research purposes, subject to ethics approval.

## Citation

If you use this code, please cite the article above (see also `CITATION.cff`).

## Licence

Code released under the MIT Licence (see `LICENSE`).
