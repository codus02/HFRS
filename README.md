# HFRS spatial epidemiology in South Korea (2011–2024)

Reproducibility package for:

> **Spatiotemporal patterns of hemorrhagic fever with renal syndrome and associations with local infrastructure and land use in South Korea, 2011–2024.**
> Chaeyeon Lee, Geonsang Lee, Seongdae Kim. Hongik University. (Prepared for submission to *Health & Place*.)

A Zenodo DOI badge will be added here after the first GitHub release is archived through the GitHub–Zenodo integration.

## What this does

A Bayesian negative-binomial disease-mapping analysis of district-level HFRS incidence across 228 South Korean districts over 3,192 district-years (2011–2024; 5,634 notified cases), fitted with **R-INLA** (23.9.9):

- Age-standardized expected counts as the offset; BYM2 spatial effect (scaled, PC priors), first-order random-walk year effect and an independent district-year effect; fixed-effect prior N(0, precision 4)
- 19 district-level covariates in six domains (land use, water environment, population and households, socioeconomic conditions, healthcare access, education and health; sex ratio and fiscal autonomy retained a priori), coded as original values, log(1+x), binary or ordinal scores and standardized
- Period-specific refits (2011–2015, 2016–2019, 2020–2024) with the full-period scaling; urban-population-share adjustment; model ladder M1–M6 (DIC/WAIC)
- Residual diagnostics: Moran's I of observed SIR, posterior predictive quantile residuals and offset-corrected Pearson residuals; CPO/PIT
- Sensitivity program: imputed-covariate exclusion, NB without district-year effects, Poisson and zero-inflated NB, additive-logratio land composition, eight neighbourhood graphs (queen, rook, k-NN 2–7), alternative PC and fixed-effect priors, exclusion of 2020–2022 and a pandemic indicator, removal of healthcare-access covariates, multicollinearity diagnostics

**Interpretation.** The covariate specification was developed through exploratory model building. The package reproduces the selected models; the associations are ecological and hypothesis-generating.

Repository: https://github.com/codus02/HFRS

## Repository layout

| Path | Content |
|---|---|
| `code/fit_main_and_sensitivity_models.R` | principal model (BASE) and the first sensitivity set (imputation exclusion, NB/Poisson/iid structure, land ALR) |
| `code/residual_diagnostics.R` | posterior predictive quantile residuals, offset-corrected Pearson residuals, random-effect SDs |
| `code/period_and_urban_models.R` | three period-specific refits, urban-share adjustment, optional M1–M6 ladder (`--ladder`) |
| `code/extended_sensitivity_graphs_priors.R` | eight neighbourhood graphs, priors, zero-inflated NB, pandemic years, healthcare-access removal, VIF/condition number, CPO/PIT |
| `code/render_main_figures.py`, `code/render_supplementary_figures.py`, `code/map_panel.py` | manuscript figures from archived results (no fitting) |
| `data/` | compiled district-year panel, codebook, scaling audit, ordinal mappings, land composition, urban covariate, neighbourhood graph, district boundaries |
| `reference_results/` | archived manuscript values used by the figure scripts (`irr_values.csv`, `sir_values.csv`, `table1_values.csv`, `supplement_spatial_values.csv`) |
| `results/` | sensitivity coefficients, fit statistics, hyperparameters and diagnostics; `results/extended/` holds the graph, prior and posterior-predictive outputs; `results/tables/` holds the submitted main and supplementary tables (Table_1, Table_2, Table_S1–S13) as CSV |
| `environment/` | R and Python package versions and session information |

## How to run

1. Clone to a path without non-ASCII characters. Install R (4.3.1 was used), R-INLA, spdep, sf, and Python 3.11 with pandas, numpy, geopandas, matplotlib (`environment/requirements.txt`).
2. From the repository root:
   ```
   Rscript code/fit_main_and_sensitivity_models.R
   Rscript code/residual_diagnostics.R
   Rscript code/period_and_urban_models.R            # add --ladder for M1-M6
   Rscript code/extended_sensitivity_graphs_priors.R
   python code/render_main_figures.py
   python code/render_supplementary_figures.py
   ```
   Fitted INLA objects are cached under `results/model_objects/` (git-ignored). Each script writes CSV outputs into `results/`.
3. Seeds are fixed (1008). INLA is deterministic up to small numerical differences across builds and thread settings.

## Data

See `DATA_PROVENANCE.md` for sources and redistribution notes, and `data/variable_codebook_EN.csv` for the 19 covariates (original Korean variable names, units, selected forms and source groups).

## Citation

See `CITATION.cff`. Please cite the article and the archived software release.

## License

The analysis source code is distributed under the MIT license (see `LICENSE`). Source datasets and boundary layers retain their original source terms; the software license does not grant additional rights over third-party data.
