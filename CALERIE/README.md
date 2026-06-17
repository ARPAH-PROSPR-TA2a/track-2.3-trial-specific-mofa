# CALERIE MOFA analysis workflow

This folder contains R scripts for running a MOFA analysis on CALERIE omics data and then testing associations between the fitted MOFA factors and sample-level metadata variables.

The workflow has two main stages, plus optional model-visualization and factor-interpretation stages:

1. **Fit a MOFA model** from one or more omics data matrices.
2. **Regress a selected metadata outcome** on each fitted MOFA factor, optionally adjusting for metadata covariates, or optionally regress each MOFA factor on a selected metadata variable.
3. **Visualize the fitted MOFA model** in a rendered HTML or PDF report.
4. **Interpret the fitted MOFA factors** using view-specific enrichment analyses of the factor loadings.

## Files in this folder

| File | Purpose |
|---|---|
| `resources.R` | Defines local file paths for the CALERIE input files. Edit this file to specify the file locations for your local compute environment. |
| `MOFA_analysis_calerie.R` | Main MOFA-fitting script. It reads the omics and phenotype files, preprocesses the selected omics views, fits a MOFA model, and writes the fitted model to `model.hdf5`.  |
| `MOFA_plots_report.Rmd` | Parameterized R Markdown report for visual inspection of a fitted `model.hdf5` file. It loads the fitted MOFA object and renders model-level plots, including data overview, variance explained, factor correlations, factor-covariate correlations, selected factor plots, top feature weights, dimensionality reduction, and data heatmaps. |
| `regress_mofa_factors.R` | Loads a fitted `model.hdf5` file, runs one regression per MOFA factor, stores the fitted models in `res`, and generates a PDF report summarizing model performance and factor associations. The outcome/explanatory metadata variable, adjustment covariates, regression direction, and categorical variables are user-selected near the top of the script. |
| `mofa_factor_interpretation.R` | Loads a fitted MOFA object, extracts factor loadings for the available omics views, runs view-specific enrichment analyses, and generates a PDF report with factor-level biological interpretation summaries. For expression, proteomics, and methylation views, the script uses signed ranked `fgsea` analysis; for metabolomics, it uses Metabolon pathway/category overrepresentation analysis. |
| `methods.R` | Helper functions sourced by `MOFA_analysis_calerie.R`, including data readers such as `readPhenoData()`, `readMethylation()`, `readProteomics()`, `readMetabolomics()`, and matrix helpers such as `fillEmptyCols()`.  |

## Required R packages

The scripts explicitly load:

```r
library(MOFA2)
library(dplyr)
library(tidyr)
```

MOFA uses the Python package `mofapy2` under the hood. The MOFA script runs:

```r
run_mofa(..., use_basilisk = TRUE)
```

so MOFA2 will use `basilisk` to manage the Python environment automatically. The first run on a new machine can take longer because the Python/MOFA dependencies may need to be installed.

The optional factor-interpretation script, `mofa_factor_interpretation.R`, uses additional Bioconductor and CRAN packages for annotation, ranked enrichment, and reporting, including `fgsea`, `AnnotationDbi`, `clusterProfiler`, `missMethyl`, `org.Hs.eg.db`, and optionally `future` / `future.apply` for parallel enrichment tasks.

The optional plotting report, `MOFA_plots_report.Rmd`, uses `ggplot2`, `corrplot`, `rmarkdown`, and the plotting functions provided by `MOFA2`. Rendering to HTML requires `rmarkdown`; rendering to PDF also requires a working LaTeX installation.

## Input data

The file `resources.R` defines the input file locations. The visit/phenotype file is required, and the omics views are optional.

### Required input

```r
calerie_clin_fn
```

This should point to the CALERIE visit data file. The MOFA script expects this file to contain subject/sample metadata, including columns used to match samples to omics measurements. In the current script, relevant metadata fields include columns such as `DEID`, `CR`, `fumo`, and omics barcode columns.

### Optional omics inputs

To exclude an omics type from the MOFA analysis, leave the corresponding file-path variable undefined or commented out in `resources.R`. The script automatically uses only the omics views whose file paths are specified.

```r
calerie_methy_fn = ...
calerie_prot_fn = ...
calerie_metab_data_fn = ...
calerie_metab_annot_fn = ...
calerie_clin_fn = ...
```

`calerie_metab_annot_fn` is also defined for metabolomics annotation support and may be used by the metabolomics reader in `methods.R`.

## Stage 1: fit the MOFA model

Open `MOFA_analysis_calerie.R` and edit the configuration block near the top of the file.

```r
code_dir <- "D:/git/misc/prospr-fast/CALERIE"
output_dir <- "D:/R_sessions/PROSPR/MOFA"

timepoint <- 24
use_change_scores <- TRUE
use_groups <- FALSE

use_residuals <- FALSE
covariates <- data.frame(
  var_names = character(0),
  is_factor = logical(0),
  stringsAsFactors = FALSE
)
use_parallel_residualization <- TRUE
n_residualization_cores <- max(1, parallel::detectCores() - 1)
```

### Configuration options

| Variable | Meaning |
|---|---|
| `code_dir` | Folder containing this code, including `resources.R` and `methods.R`. |
| `output_dir` | Folder where MOFA result folders will be created. |
| `timepoint` | Follow-up visit to analyze, `12` or `24` months. Samples are selected using `fumo == timepoint`. |
| `use_change_scores` | If `TRUE`, each omics value is converted to a change score: target visit minus baseline. If `FALSE`, the script uses target-visit values directly. |
| `use_groups` | If `TRUE`, the MOFA model treats the AL and CR arms as separate MOFA groups. If `FALSE`, all samples are analyzed together without MOFA grouping. |
| `use_residuals` | If `TRUE`, residualizes each omics data matrix against selected phenotype covariates before any downstream derivations such as change-score calculation or timepoint extraction. The default is `FALSE`. |
| `covariates` | A two-column data frame specifying phenotype covariates to use for residualization. `covariates$var_names` must contain column names from the phenotype table loaded from `calerie_clin_fn`; `covariates$is_factor` indicates whether each covariate should be coerced to a factor before residualization. |
| `use_parallel_residualization` | If `TRUE`, uses parallel feature-wise residualization when missing omics values prevent fully vectorized residualization. |
| `n_residualization_cores` | Number of worker processes to use for parallel residualization. |

### Optional residualization of omics data

The MOFA script can optionally residualize each omics matrix before the matrices are used in downstream calculations. This is controlled by:

```r
use_residuals <- FALSE

covariates <- data.frame(
  var_names = character(0),
  is_factor = logical(0),
  stringsAsFactors = FALSE
)
```

When `use_residuals` is `FALSE`, no residualization is performed and the analysis proceeds as usual.

When `use_residuals` is `TRUE`, each omics matrix is residualized immediately after it is read from file and before downstream derivations such as change-score calculation, timepoint-specific extraction, variance filtering, feature selection, or MOFA model construction. For each omics view, the script:

1. Identifies the barcode column in the phenotype table that corresponds to the omics matrix.
2. Removes samples from the omics matrix if their barcodes are not present in the phenotype table.
3. Restricts and orders the phenotype table so that its rows match the remaining omics samples.
4. Coerces any covariates listed with `is_factor == TRUE` to factors.
5. Fits a linear model for each omics feature using `covariates$var_names` as regressors.
6. Replaces the original feature values with residuals from the fitted models.

For example:

```r
use_residuals <- TRUE

covariates <- data.frame(
  var_names = c("female", "bmistrat", "deidsite"),
  is_factor = c(TRUE, TRUE, TRUE),
  stringsAsFactors = FALSE
)
```

In this example, each omics feature is residualized against `female`, `bmistrat`, and `deidsite`, with all three variables treated as categorical predictors. The names in `covariates$var_names` must match column names in the phenotype table read from `calerie_clin_fn`.

If a covariate is numeric but should be modeled categorically, set its corresponding `is_factor` value to `TRUE`. If a covariate should be modeled as a continuous predictor, set `is_factor` to `FALSE`.

The script uses a fast vectorized residualization step when possible. If feature-level missing values require separate models per feature, the script can use parallel residualization controlled by:

```r
use_parallel_residualization <- TRUE
n_residualization_cores <- max(1, parallel::detectCores() - 1)
```

To disable parallel residualization, set:

```r
use_parallel_residualization <- FALSE
```

### What the MOFA script does

`MOFA_analysis_calerie.R` performs the following steps:

1. Sources `resources.R` and `methods.R`.
2. Reads the phenotype/clinical data.
3. Reads whichever omics files are defined in `resources.R`.
4. If `use_residuals == TRUE`, residualizes each omics matrix against the phenotype covariates specified in `covariates`.
5. For each available omics view, either:
   - computes change scores as `timepoint - baseline`, or
   - extracts values from the selected follow-up timepoint.
6. Removes zero-variance features.
7. Drops omics views with no samples at the requested timepoint.
8. Optionally infers AL/CR group labels from sample IDs.
9. Filters each omics view to retain highly variable features.
10. Adds missing sample columns filled with `NA`, as required by MOFA when views have different sample coverage.
11. Creates and prepares the MOFA object.
12. Adds sample metadata to the MOFA object using `samples_metadata(mo) <- meta_d`.
13. Fits the model with `run_mofa()`.
14. Writes run outputs to a timestamped result folder.

### Feature filtering

Large omics views are reduced to their most variable features before MOFA fitting. The relevant settings are inside `MOFA_analysis_calerie.R`:

```r
max_factor <- 5
max_features <- 10000
```

The script keeps no more than `max_factor` times the number of features in the smallest available view, and no more than `max_features` features per view.

### MOFA training settings

The current script uses:

```r
train_opts$drop_factor_threshold <- 0.01
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 4000
```

The training settings for `convergence_mode` and `train_opts$maxiter` are set more strictly than the defaults to enable more accurate model training. 

### Running the MOFA script

The script includes code that logs the exact commands executed to ensure complete provenance. For this to work correctly, the script must be run in a fresh R GUI session (i.e., one with no prior commands) by copy-pasting the entire script, rather than sourcing it. 

## MOFA output files

Each MOFA run creates a timestamped result folder under `output_dir`. The folder name records the analysis settings, for example:

```text
2026_May_04_10_41_49_visit=24mo_change_scores_ungrouped_views=me+pr+mb
```

The exact folder name depends on the run time, selected visit, whether change scores were used, whether groups were used, and which omics views were included.

Common output files include:

| Output | Description |
|---|---|
| `model.hdf5` | Fitted MOFA model. This is the main file used by downstream analyses. |
| `mofa_out.txt` | Log file intended to capture output from the MOFA run. |
| `sessionInfo.txt` | R session information for reproducibility. |
| `MOFA_analysis_calerie.R` | Saved R command history. |
| `README.txt` | Short run-specific description generated by the MOFA script. |

The fitted model stored in `model.hdf5` contains the fitted data views, feature weights, predicted sample factor values, and the sample metadata table added before training.

## Optional stage: visualize the fitted MOFA model

After fitting a MOFA model, use `MOFA_plots_report.Rmd` to generate a rendered HTML or PDF report that provides a visual overview of the fitted model. This report is intended as a model-inspection step before, alongside, or after the regression analysis. It does not refit the MOFA model; it reads an existing `model.hdf5` file produced by `MOFA_analysis_calerie.R`.

Edit the report parameters near the top of `MOFA_plots_report.Rmd` or pass them to `rmarkdown::render()`. The most important parameter is the path to the fitted MOFA model:

```r
mofa_res_fname <- "D:/R_sessions/PROSPR/MOFA/.../model.hdf5"
```

The report also includes user-selected plotting settings such as the metadata columns to correlate with factors, the factors to highlight in detailed plots, the metadata column used for sample coloring, the views to use for top-weight plots, and the views to use for heatmaps.

To render the report to HTML, run:

```r
rmarkdown::render(
  input = "MOFA_plots_report.Rmd",
  output_format = "html_document",
  output_file = "MOFA_plots_report.html",
  output_dir = "D:/R_sessions/PROSPR/MOFA/reports"
)
```

To render the report to PDF, use `output_format = "pdf_document"` and set `output_file` to a `.pdf` filename. The output folder is controlled by the `output_dir` argument to `rmarkdown::render()`, not by the fitted MOFA model folder.

The report includes plots for data availability, variance explained by each factor, pairwise factor correlations, factor-covariate correlations, selected factor distributions, pairwise selected-factor scatterplots, top feature weights, t-SNE/UMAP embeddings of samples in factor space, and heatmaps of high-weight features for selected factors and views.

## Stage 2: regress metadata outcomes on MOFA factors

After fitting a model, use `regress_mofa_factors.R` to test whether a selected metadata variable is associated with each MOFA factor. Setting `factors_as_Y <- FALSE` treats the selected metadata variable `Y` as the dependent variable and each factor as the explanatory variable. Alternatively, setting `factors_as_Y <- TRUE` treats each MOFA factor as the dependent variable and tests `Y` as the explanatory variable of interest. In both cases, additional covariates can be specified in the variable `X`.

Open `regress_mofa_factors.R` and edit the user-specified inputs near the top of the file.

```r
mofa_res_fname <- "D:/R_sessions/PROSPR/MOFA/.../model.hdf5"

X <- c("bmistrat")
Y <- "CR"
factors_as_Y <- FALSE

categorical <- c(
  "DEID", "CR", "deidsite", "strata", "female", "ethnic", "race",
  "race3", "bmistrat", "visit", "fu", "fumo", "Visit.ID"
)
```

### Regression inputs

| Variable | Meaning |
|---|---|
| `mofa_res_fname` | Full path to the fitted MOFA model, i.e., the `model.hdf5` file created by `MOFA_analysis_calerie.R`. |
| `X` | Character vector of metadata covariates to include in every model. Use `character(0)` for no covariates. |
| `Y` | Name of the metadata column to use as the dependent variable when `factors_as_Y == FALSE`, or as the explanatory variable of interest when `factors_as_Y == TRUE`. |
| `factors_as_Y` | Boolean variable controlling the regression direction. If `FALSE`, the script fits `Y ~ Factor_i + X_covariates`. If `TRUE`, the script fits `Factor_i ~ Y + X_covariates`. |
| `categorical` | Metadata columns that should be treated as categorical variables. Include `Y` here when `factors_as_Y == FALSE` and the outcome should be modeled using logistic regression. Also include `Y` or any covariates here when they should be treated as categorical explanatory variables in `factors_as_Y == TRUE` mode. |
| `regression_report_fname` | Output path for the PDF regression report. By default, this is written in a sub-folder within the directory where `model.hdf5` is located, and the filename includes the outcome/explanatory variable `Y` and the covariates in `X`. |

### Model fitted by the regression script

The script loads the fitted MOFA model:

```r
mo.tr <- MOFA2::load_model(mofa_res_fname)
```

Then it runs:

```r
res <- runRegression(mo.tr, X, Y, categorical)
```

The regression direction is controlled by `factors_as_Y`.

When:

```r
factors_as_Y <- FALSE
```

the script fits one model per MOFA factor of the form:

```r
Y ~ Factor_i + X_covariates
```

Internally, the current factor is stored in a temporary column named `F`, so the formula is constructed as:

```r
Y ~ F + covariates
```

When:

```r
factors_as_Y <- TRUE
```

the script instead fits one linear model per MOFA factor of the form:

```r
Factor_i ~ Y + X_covariates
```

Internally, the current factor is again stored in a temporary column named `F`, so the formula is constructed as:

```r
F ~ Y + covariates
```

In this mode, the report focuses on the coefficient or term-level p-value for `Y`, because `Y` is the explanatory variable being tested against each factor.

### Logistic versus linear regression

When `factors_as_Y == FALSE`, `runRegression()` chooses the model type from the outcome variable `Y`:

| Outcome type | Model used |
|---|---|
| `Y` listed in `categorical`, or otherwise stored as a factor/character/logical variable | Logistic regression using `glm(..., family = binomial())` |
| Numeric `Y` not listed in `categorical` | Linear regression using `lm()` |

Important: in this  mode, a numeric `0/1` outcome is treated as continuous unless its column name is included in `categorical`. 

For logistic regression, the outcome must have exactly two observed classes after missing values are removed for a given factor model. The first factor level is used as the reference level by R. Use `relevel()` in the metadata before modeling if you need a specific reference category.

When `factors_as_Y == TRUE`, each MOFA factor is the dependent variable, so the script always fits a linear regression model. In this mode, `Y` is an explanatory variable. If `Y` is categorical, include it in `categorical` so it is represented using R's standard factor contrasts. For categorical explanatory variables with more than two levels, the report uses an omnibus term-level p-value for `Y`, rather than a single dummy-coefficient p-value.

### Regression output

The object `res` is a named list with one entry per MOFA factor:

```r
res[["Factor 1"]]
res[["Factor 2"]]
```

When `factors_as_Y == FALSE`, each entry is either a `glm` object for logistic regression or an `lm` object for linear regression. When `factors_as_Y == TRUE`, each entry is an `lm` object because the MOFA factor is the continuous dependent variable.

Useful inspection commands:

```r
attr(res, "model_type")
attr(res, "formula")
summary(res[["Factor 1"]])
lapply(res, summary)
```

For logistic regression in the default `factors_as_Y == FALSE` mode, coefficients are on the log-odds scale. Exponentiate coefficients to obtain odds ratios:

```r
sapply(res, function(model) exp(coef(model)))
```

### Regression PDF report

After running:

```r
res <- runRegression(mo.tr, X, Y, categorical)
```

the script generates a PDF report that quantitatively summarizes the performance of the regression models stored in `res`. The report is intended to provide a quick overview of which MOFA factors are most strongly associated with the selected outcome or explanatory variable, as well as basic model-performance diagnostics.

By default, the report is written to a sub-directory within the same folder as the fitted MOFA model. The report filename includes both the outcome/explanatory variable and the adjustment covariates. For example:

```text
mofa_factor_regression_report_Y=CR_covar=female+race.pdf
```

If no covariates are included in `X`, the filename uses `covar=none`:

```text
mofa_factor_regression_report_Y=CR_covar=none.pdf
```

The report includes explanatory text tailored to the type and direction of regression that was run. When `factors_as_Y == FALSE`, the report summarizes the coefficient for the factor term, its standard error, p-value, and confidence interval. When `factors_as_Y == TRUE`, the report instead summarizes the coefficient or omnibus term-level p-value for `Y`, because `Y` is the explanatory variable of interest. The model-performance summary is sorted by the nominal p-value for the tested term, with the smallest p-values listed first. If the number of factors exceeds the number that can fit on one page, the table automatically continues across additional pages.

For logistic-regression analyses, the report includes classification-oriented summaries such as AUC, accuracy, McFadden pseudo-R², Brier score, AIC, and BIC, along with plots summarizing factor effect estimates, nominal p-values, and model performance metrics. Logistic regression is used only in the default `factors_as_Y == FALSE` mode when `Y` is categorical or binary.

For linear-regression analyses, the report includes continuous-outcome summaries such as R², adjusted R², RMSE, MAE, AIC, and BIC, along with plots summarizing the tested effect estimates, nominal p-values, and model performance metrics. Linear regression is used for continuous outcomes in the default mode and for all models when `factors_as_Y == TRUE`.

The report also includes diagnostic plots for the model with the smallest nominal p-value for the tested term. For linear models, these diagnostics focus on fitted values and residual behavior. For logistic models, they focus on predicted probabilities and observed outcome classes.

If `ggplot2` is installed, the script uses it for report plots. Otherwise, the report is still generated using base R graphics.

## Optional stage: interpret MOFA factors by enrichment analysis

After fitting a MOFA model, use `mofa_factor_interpretation.R` to generate a biological interpretation report for the fitted factors. This script is intended to complement the regression workflow: instead of testing whether a factor is associated with a metadata outcome, it asks which biological pathways or feature categories are represented by the features with strong loadings on each factor.

The script extracts per-factor feature loadings from the fitted MOFA object and analyzes the views present in the model. It expects the view names to be consistent with the MOFA-fitting workflow, such as `proteomics`, `methylation`, `metabolomics`, and `expression`.

For expression, proteomics, and methylation views, the script uses a signed ranked enrichment strategy based on `fgsea`. For each factor and view, it builds one full ranked loading vector and runs enrichment once; the sign of the normalized enrichment score indicates whether a gene set is associated with the positive or negative side of the factor. Proteomics features are mapped from UniProt IDs to Entrez IDs, expression features are interpreted as Entrez IDs or gene symbols, and methylation CpG features are mapped to Entrez IDs using `missMethyl`.

For metabolomics, the script uses the Metabolon annotation file to perform categorical overrepresentation analysis against fields such as `SUPER_PATHWAY`, `SUB_PATHWAY`, `PARENT_CLASSIFICATION`, and `CLASSIFICATION`. This branch remains set-based because the Metabolon pathway/category annotations are not gene sets and are interpreted separately from the ranked gene-based analyses.

A typical call is:

```r
mofa_obj <- MOFA2::load_model(mofa_res_fname)
metabolon_annot <- readr::read_tsv(calerie_metab_annot_fn, show_col_types = FALSE)

res <- run_mofa_integrative_interpretation(
  mofa_obj = mofa_obj,
  metabolon_annotation = metabolon_annot,
  output_dir = "mofa_factor_interpretation_report",
  metabolomics_feature_col = "BIOCHEMICAL",
  expression_id_type = "auto",
  methylation_array_type = "EPIC",
  report_significance_cutoff = 0.05,
  render_report = TRUE,
  include_top_features_in_report = FALSE,
  report_engine = "base_pdf",
  n_cores = 4,
  fgsea_nproc = 1
)
```

The default `base_pdf` report engine writes a PDF directly with base R graphics and does not require Pandoc or LaTeX. The optional `rmarkdown_pdf` engine can produce a more polished PDF, but it requires `rmarkdown`, Pandoc, and a LaTeX installation.

The most important outputs are written to `output_dir` and include:

| Output | Description |
|---|---|
| `mofa_integrative_factor_report.pdf` | PDF report summarizing factor-level enrichment results and biological interpretation text. |
| `enrichment_results.csv` | Machine-readable table of enrichment results across factors and views. |
| `top_factor_features.csv` | Top-loading features per factor and view. These are included in the PDF only when `include_top_features_in_report = TRUE`. |
| `selected_factor_features.csv` | Features selected by the configured loading-selection rule, used primarily for set-based metabolomics enrichment. |
| `all_mofa_loadings.csv` | Full table of extracted MOFA feature loadings. |

## Recommended workflow

1. Edit `resources.R` so all local file paths are correct.
2. Edit the configuration block in `MOFA_analysis_calerie.R`:
   - `code_dir`
   - `output_dir`
   - `timepoint`
   - `use_change_scores`
   - `use_groups`
   - `use_residuals`
   - `covariates`, if residualization is enabled
3. Run the MOFA analysis.
4. Locate the resulting `model.hdf5` file in the timestamped output folder.
5. Optionally edit and render `MOFA_plots_report.Rmd` to inspect the fitted model visually.
6. Edit `mofa_res_fname`, `X`, `Y`, `factors_as_Y`, and `categorical` in `regress_mofa_factors.R`.
7. Run the regression script.
8. Inspect the returned `res` object to review factor-outcome or outcome-factor associations.
9. Review the generated PDF report, whose filename includes the selected outcome/explanatory variable and covariates, for example `mofa_factor_regression_report_Y=CR_covar=female+race.pdf`.
10. Optionally run `mofa_factor_interpretation.R` on the same fitted MOFA model to generate a pathway- and feature-level interpretation report for the MOFA factors.
