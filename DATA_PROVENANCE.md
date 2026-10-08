# Data provenance and redistribution note

The analysis uses aggregated district-year counts and administrative indicators; it contains no individual-level records or direct personal identifiers.

## Input domains

- HFRS notifications by district of residence and age group, and resident-registration population denominators: Korea Disease Control and Prevention Agency infectious-disease surveillance (https://dportal.kdca.go.kr) and the Korean Statistical Information Service (https://kosis.kr)
- Land-category areas (miscellaneous, park, warehouse, salt-pan, religious land; elevated-road and riverside green-space areas): national land-use statistics
- Domestic groundwater wells; population using septic tanks outside sewer service areas: groundwater-use and sewerage statistics
- Housing and vacant dwellings, National Pension income reporters, healthcare facilities and workforce, long-term care agencies, primary schools, welfare facilities for older adults, sex ratio, fiscal autonomy: public administrative statistics compiled through KOSIS

Source groups and original Korean variable names for the 19 covariates are listed in `data/variable_codebook_EN.csv`. Missing years within a district were linearly interpolated and endpoints carried; the fill flags are retained in `data/main_covariates_pretransform_with_fill_flags.csv`.

## Compiled analytic files

- `data/model_panel_numeric.csv`: the 3,192-row district-year panel used by every reported model (cases, age-standardized expected counts, standardized covariates V01–V19, graph and time indices).
- `data/neighbors_228.graph`, `data/neighbors_228.csv`: the neighbourhood graph and district order; `data/district_boundaries.gpkg`: district polygons used for maps and for the alternative graphs.
- `data/land_composition_percent.csv`, `data/urban_sensitivity_covariate.csv`: inputs for the compositional and urbanization sensitivity analyses.

Before public release, the depositing authors must verify that redistribution of each derived field and of the boundary layer complies with the terms of its source portal. If any source terms prohibit redistribution, remove the affected file from the public release and retain the codebook and source descriptions so that authorised users can reconstruct it.

## Interpretation

The 19-covariate specification was developed through exploratory model building on the same data. The package reproduces the reported models and diagnostics; it does not turn the analysis into a prospectively specified confirmatory study.
