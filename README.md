# Metaoverfitpaper — Transportability-Weighted Meta-Analysis (TWMA)

Reproducible R analysis for **"The Overfitting Crisis in Meta-Regression:
Evidence from Optimism-Corrected Transportability."**

## What this is

Meta-regression borrows strength from study-level moderators to explain
between-study heterogeneity. This analysis asks how much of that explained
heterogeneity actually **transports** to new, unseen studies — versus how much
is just **overfitting** (in-sample optimism) that vanishes under honest
out-of-sample evaluation.

The central quantity is a **transportability statistic T²**:

```
T² = 1 - tau²_full / tau²_null
```

i.e. the proportion of residual heterogeneity that the moderators remove. The
paper contrasts the **apparent** T² (fit and evaluated on the same studies)
against an **optimism-corrected** version `T̃²`, where the residual variance is
estimated by **leave-one-study-out cross-validation**:

```
tau²_cv = mean( max(0, (y_i - ŷ_{-i})² - v_i) )
T̃²      = 1 - tau²_cv / tau²_null
```

All resampling (LOO, bootstrap, conformal) is done **at the study level** — each
study is one unit — so there is no row/study leakage in the cross-validation.
Transportability variances and T² values are clamped to be non-negative.

The script applies this to four published meta-analyses from the `metadat`
package (BCG vaccine, teacher-expectancy, St John's Wort, Konstantopoulos
education) and then surveys 16 published meta-regressions ("Wall of Shame") to
score their overfitting risk by k (number of studies) and k/p (studies per
parameter).

## Files

| File | Purpose |
|------|---------|
| `twma_analysis.R` | Main analysis script (clean, parseable R). |
| `test_twma.R` | Smoke test: sources the functions and `stopifnot()`s basic properties. |
| `metaoverfitpaper` | Original source (R wrapped in Markdown code fences — kept for provenance; **not** runnable as-is). |

> **Note on the original file.** `metaoverfitpaper` (no extension) was R code
> wrapped in Markdown ```` ```r ```` / ```` ``` ```` fences plus a prose line. Those
> three non-R lines were stripped to produce the runnable `twma_analysis.R`;
> nothing else in the analysis was changed.

## Requirements

R (tested on 4.6.0) with:

```r
install.packages(c("metafor", "metadat", "knitr", "ggplot2"))
```

- **Required:** `metafor`, `metadat`.
- **Optional:** `knitr` (pretty tables — a base-R fallback is used if absent),
  `ggplot2` (the overfitting-risk scatter figure — skipped if absent).

The script does **not** install packages at runtime; the `install.packages(...)`
call from the original code is commented out.

## How to run

```sh
Rscript twma_analysis.R          # full analysis (a few minutes: bootstrap + nested LOO conformal)
Rscript test_twma.R              # fast smoke test, exits 0 on success
```

## Outputs

Console output includes:

- Standard transportability T² per dataset (with bootstrap 95% CIs).
- Optimism-corrected T̃² and the optimism gap (apparent − corrected).
- Leave-one-study-out cross-validation performance (RMSE, MAE, coverage, correlation).
- Parametric vs **conformal** prediction-interval coverage (target 95%).
- Trim-and-fill publication-bias sensitivity.
- A "Wall of Shame" survey of 16 published meta-regressions by overfitting risk.

It also writes (both git-ignored, regenerated on each run):

- `published_metaregressions.csv` — the Wall-of-Shame table.
- `overfitting_risk.png` — k-vs-p risk scatter (only if `ggplot2` is installed).

### Representative results (seed = 42, R 4.6.0)

| Dataset | k | p | T² (apparent) | T̃² (corrected) | Optimism |
|---------|---|---|---------------|-----------------|----------|
| BCG | 13 | 4 | 0.572 | 0.000 | 0.572 |
| Teacher | 19 | 2 | 0.406 | 0.000 | 0.406 |
| St John's Wort | 22 | 4 | 0.773 | 0.430 | 0.344 |
| Konstantopoulos (simple) | 56 | 2 | 0.000 | 0.000 | 0.000 |
| Konstantopoulos (complex) | 56 | 10 | 0.640 | 0.245 | 0.395 |

Key finding: datasets with k < 20 show **complete overfitting** (T̃² ≈ 0), and
87.5% of the 16 surveyed published meta-regressions fall in the extreme/severe
overfitting-risk band.

## Notes / caveats

- `metafor` emits benign "Fisher scoring … setting tau² = 0" warnings for some
  small-k full models; this is expected when moderators absorb all heterogeneity.
- Bootstrap and conformal sections use their own fixed seeds inside the
  functions, so results are reproducible run-to-run.
- If `knitr`/`ggplot2` are missing, the corresponding tables/figure degrade
  gracefully (base-R table print / figure skipped) — the numeric analysis is
  unaffected.
