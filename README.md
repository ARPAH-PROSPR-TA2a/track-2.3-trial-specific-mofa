# MOFA analysis workflow

This repository contains a reusable workflow for fitting Multi-Omics Factor
Analysis (MOFA) models to study-specific omics data and carrying out
post-MOFA analyses of the learned factors. The analysis code itself is intended to be reusable
across studies; the parameter files specify study-specific input files,
clinical column names, labels, preprocessing choices, and downstream analysis
settings.

The main analysis workflow is:

1. **Run `mofa_analysis.R`** to prepare the selected omics views and fit a MOFA
   model.
2. **Run `regress_mofa_factors.R`** to test associations between MOFA factors
   and selected sample-level metadata variables.
3. **Render `mofa_plots_report.Rmd`** to generate an optional visual report of
   the fitted MOFA model.
4. **Run `mofa_factor_interpretation.R`** to generate an optional biological
   interpretation report for the fitted factors.
5. **Run `mofa_train_glm.R`** to train view-specific LASSO/elastic-net models
   that predict the MOFA factor values learned in step 1.

The files `glm_factor_estimate.R` and `project_mofa_factors.R` are intended for
use by **external studies** after the primary MOFA model has been learned. They
allow a study that did not participate in MOFA fitting to estimate the learned
factor values for its own samples using either trained GLM predictors or the
saved MOFA feature weights. These external-scoring workflows are described
below.

## Files in the repository

| File | Purpose |
|---|---|
| `mofa_parameters.R` | Main configuration file for `mofa_analysis.R`. Defines study-specific settings, analysis parameters, input-file locations, and optional omics-view definitions. |
| `mofa_analysis.R` | Main MOFA workflow. Reads the clinical and omics data, performs the configured preprocessing, fits the MOFA model, and writes the model and derived factor/weight/metadata files. It can also perform a metadata-only run when no omics files are specified. |
| `file_readers.R` | Study-specific file-reader functions used to convert the clinical and omics input files into the formats expected by the workflow. It currently contains `readPhenoData()`, `readMethylation()`, `readProteomics()`, and `readMetabolomics()`. |
| `methods.R` | Shared helper functions used by `mofa_analysis.R`, `mofa_train_glm.R`, `glm_factor_estimate.R`, and `project_mofa_factors.R`, including clinical configuration validation, sample-ID construction, timepoint/change-score handling, omics loading/preparation, and feature scaling. |
| `regression_parameters.R` | Configuration file for `regress_mofa_factors.R`. Selects the regression input source and defines metadata variables, covariate sets, and categorical variables. |
| `regress_mofa_factors.R` | Runs regression between selected metadata variables and every available MOFA factor (in both directions, factor as the predictor and as the response variable) for multiple user-configured covariate sets, generates PDF reports, and writes consolidated Excel and RDS summaries of the fitted models and coefficients. |
| `mofa_plots_report_parameters.R` | Configuration file for `mofa_plots_report.Rmd`. Identifies the fitted MOFA model and controls the metadata variables, factors, views, feature-weight plots, heatmaps, and dimensionality-reduction analyses included in the report. |
| `mofa_plots_report.Rmd` | R Markdown report for visual inspection of a fitted MOFA model. User-configurable settings are read from `mofa_plots_report_parameters.R`. |
| `mofa_factors_interpretation_parameters.R` | Configuration file for `mofa_factor_interpretation.R`. Identifies the fitted MOFA model and optional Metabolon annotation file and controls view/factor selection, feature-selection settings, GO ontologies, significance thresholds, report generation, and parallelization. |
| `mofa_factor_interpretation.R` | Biological interpretation workflow for fitted MOFA factors. Performs ranked GO enrichment for gene-oriented views, Metabolon over-representation analysis, and optional PDF reporting. |
| `mofa_train_glm_parameters.R` | Configuration file for `mofa_train_glm.R`. Identifies the fitted MOFA model and controls factor/view selection, repeated validation, elastic-net settings, parallelization, and the training report. |
| `mofa_train_glm.R` | Trains view-specific LASSO/elastic-net models that predict the learned MOFA factors from individual omics views. The resulting predictor object is used for later external factor estimation. |
| `glm_factor_estimate_parameters.R` | Configuration file for applying trained GLM factor predictors to an external study. Defines the trained predictor, external data files, preprocessing, and prediction options. |
| `glm_factor_estimate.R` | Estimates learned MOFA factor values in an external dataset using the view-specific GLM/LASSO predictors trained by `mofa_train_glm.R`. |
| `project_mofa_factors_parameters.R` | Configuration file for projecting an external dataset directly onto learned MOFA feature weights. |
| `project_mofa_factors.R` | Estimates MOFA factor values in an external dataset by applying saved MOFA feature weights to the corresponding external omics views. |

## Required R packages

The primary MOFA workflow explicitly uses:

```r
library(MOFA2)
library(dplyr)
library(tidyr)
library(xlsx)
```

`mofa_train_glm.R` additionally requires `glmnet` and, when parallelization
is enabled, `future` and `future.apply`. Its training report uses `ggplot2`.

`regress_mofa_factors.R` uses the standard R modeling functions and `xlsx`
for Excel input/output. It also writes an RDS copy of the consolidated
regression-summary tables using base R.

The plotting report uses `ggplot2`, `corrplot`, `rmarkdown`, and MOFA2 plotting
functions. Rendering to PDF requires a working LaTeX installation.

The factor-interpretation workflow additionally uses packages for annotation,
enrichment analysis, Excel output, and report generation, including `fgsea`,
`missMethyl`, `AnnotationDbi`, `org.Hs.eg.db`, `GO.db`, `clusterProfiler`,
`xlsx`, and `rmarkdown`. Parallel execution on Windows uses `future` and
`future.apply`. See the documentation in `mofa_factor_interpretation.R` for
the current requirements.

MOFA2 uses the Python package `mofapy2` under the hood. The main MOFA workflow
runs MOFA with `use_basilisk = TRUE`, allowing `basilisk` to manage the Python
environment. The first run on a new machine may therefore take longer while
the required dependencies are installed.

## Repository organization and study-specific configuration

The repository is designed so that the code can be used, essentially unchanged,
to execute the MOFA workflow across multiple studies. Study-specific details are
specified through settings in the corresponding parameter files rather than
being embedded in the analysis scripts. For example, the following types of 
settings belong in `mofa_parameters.R`:

- study name and clinical column names
- mappings between clinical codes and user-facing labels
- target timepoint and change-score settings
- omics input-file paths
- view-specific barcode columns and reader functions
- feature-selection and MOFA training settings

The downstream workflows follow the same pattern. For example,
`regression_parameters.R`, `mofa_plots_report_parameters.R`,
`mofa_factors_interpretation_parameters.R`, `mofa_train_glm_parameters.R`,
`glm_factor_estimate_parameters.R`, and `project_mofa_factors_parameters.R`
contain the settings specific to those workflows.

The reader functions in `file_readers.R` are the main place where a study's
physical file formats are translated into the standardized matrices expected by
the workflow. The shared processing and validation logic belongs in
`methods.R`.

## Running the scripts

The scripts are designed to be run after configuring the corresponding
parameter file. When a script sources supporting files using relative paths,
run it from the repository directory or use `source(..., chdir = TRUE)` so that
the working directory is temporarily set to the directory containing the
script.

For example:

```r
source("/path_to_mofa_files/mofa_analysis.R", local = TRUE, chdir = TRUE)
```

or from the command line:

```text
Rscript.exe -e "source('/path_to_mofa_files/mofa_analysis.R', local = TRUE, chdir = TRUE)"
```

Running `mofa_analysis.R` from an R session is recommended when MOFA training
output needs to be captured in the log, because `MOFA2::run_mofa()` executes
Python code in a separate environment and its progress messages are not always
routed through R's standard output when launched through `Rscript.exe`.

# 1. Fit the MOFA model with `mofa_analysis.R`

## Configure `mofa_parameters.R`

The first step is to edit `mofa_parameters.R`. The required clinical file and
any desired omics files are specified there. The most important settings are:

| Variable | Meaning |
|---|---|
| `output_dir` | Directory in which timestamped MOFA result folders are created. |
| `study_name` | Study name used in generated documentation and as the worksheet name for the factor-values and sample-metadata Excel files. |
| `subject_id_col` | Column name in the clinical/phenotype file containing subject identifiers. |
| `treatment_col` | Column name in the clinical/phenotype file containing treatment codes. |
| `timepoint_col` | Column name in the clinical/phenotype file containing timepoint codes. |
| `timepoint_labels` | Named mapping between user-facing timepoint labels and the codes present in `timepoint_col`; exactly one entry must be named `baseline`. |
| `treatment_labels` | Named mapping between user-facing treatment labels and the codes present in `treatment_col`. |
| `timepoint` | Target visit to analyze. This must be one among the values found under `timepoint_col` in the clinical/phenotype file |
| `use_change_scores` | If `TRUE`, analyze target-minus-baseline change scores; if `FALSE`, analyze target-visit values directly. |
| `max_features` | Absolute upper limit on the number of retained features in a view. |
| `max_factor` | Factor for capping number of features in large views. A view may have at most `max_factor` times as many features as the smallest view, and never more than `max_features` |
| `scale_views` | Whether MOFA data views should be centered/scaled before model fitting. |
| `mofa_max_iter` | Maximum number of variational-inference iterations used by MOFA training. |
| `resource_vars` | Names of variables that specify file locations in the parameter file. |
| `omics_input_defs` | Definitions of the optional views, including view names, file variables, clinical barcode columns, reader functions, and run-specific documentation. |

For each view in `omics_input_defs`, the `reader` function must return a
numeric **features x samples** matrix. Feature identifiers belong in the row
names and sample barcodes in the column names.

For example:

```r
output_dir <- "D:/R_sessions/PROSPR/MOFA"
study_name <- "CALERIE"

subject_id_col <- "DEID"
treatment_col <- "CR"
timepoint_col <- "fumo"

timepoint_labels <- setNames(
  c(0, 12, 24, 3, 6, 18, 98, 99), 
  c("baseline", "12mo", "24mo", "3mo", "6mo", "18mo", "BLs1", "BLs2")
)

treatment_labels <- setNames(c(0, 1), c("AL", "CR"))

timepoint <- 24
use_change_scores <- TRUE

max_factor <- 5
max_features <- 10000
scale_views <- TRUE
mofa_max_iter <- 4000
```

The corresponding file variables are then defined for the omics and clinical
inputs, for example `methy_fn`, `prot_fn`, `metab_data_fn`,
`metab_annot_fn`, and `clin_fn`.

### Omics views are optional

One or more omics views may be included. To exclude a view, leave its file
variable undefined or comment it out. The script automatically processes only
the configured views.

A special case is a **metadata-only run**. If all omics file variables are
omitted, `mofa_analysis.R` does not fail. Instead, it constructs the target-
visit sample metadata and writes the metadata and provenance files described
below, then completes without running MOFA.

### Sample IDs and metadata

The sample IDs generated by `mofa_analysis.R` have the form:

```text
<subject ID>_<treatment label>_<timepoint label>
```

The clinical columns and code-to-label mappings used to construct these IDs are
controlled by `subject_id_col`, `treatment_col`, `timepoint_col`,
`treatment_labels`, and `timepoint_labels`.

### `file_readers.R`

Edit the reader functions in `file_readers.R` when the physical format of a
study's clinical or omics files differs from the expected format. The main
script calls these readers through the `reader` functions specified in
`omics_input_defs`.

### `methods.R`

The shared helper methods in `methods.R` should normally be changed only when
the common workflow behavior itself needs to change. They are used by multiple
scripts and therefore provide the common implementation of clinical
validation, sample-ID construction, target-visit/change-score processing,
optional view loading, and related operations.

## What `mofa_analysis.R` does

For a normal MOFA run, the script:

1. Loads `mofa_parameters.R`, `file_readers.R`, and `methods.R`.
2. Reads and validates the clinical/phenotype table.
3. Reads every configured omics view.
4. Extracts target-visit values or computes target-minus-baseline change scores,
   according to `use_change_scores`.
5. Removes zero-variance features.
6. Drops views with no samples at the requested visit.
7. Retains the most variable features subject to `max_factor` and
   `max_features`.
8. Adds missing sample columns filled with `NA` when views do not have complete
   sample coverage.
9. Creates and prepares a MOFA object.
10. Adds the target-visit sample metadata to the MOFA object.
11. Trains the model with `MOFA2::run_mofa()`.
12. Exports the learned factor values, feature weights, and sample metadata in
    both Excel and RDS formats. The persisted feature weights can subsequently
    be used by `project_mofa_factors.R` to project the learned MOFA factors in
    external studies.
13. Saves the training log, session information, and source/configuration files
    in the result folder's `logs` subfolder for provenance.

The script fits all samples together as a single MOFA group. Treatment labels
are retained in the sample IDs and metadata but are not used to define separate
MOFA groups.

## Outputs from `mofa_analysis.R`

Each run creates a timestamped folder beneath `output_dir`. The exact
folder name records the run time, target visit, score mode, and included views.

The principal outputs are:

| Output | Description |
|---|---|
| `model.hdf5` | Trained MOFA model. This is the primary input to `mofa_train_glm.R`, `regress_mofa_factors.R` option 1, `mofa_plots_report.Rmd`, and `mofa_factor_interpretation.R`. |
| `logs/mofa_out.txt` | Log containing the training iteration output from `MOFA2::run_mofa()`. |
| `factor_values.xlsx` | Samples x factors factor-value matrix. The worksheet is named using `study_name`. |
| `factor_values.rds` | RDS version of the samples x factors factor-value matrix. |
| `feature_weights.xlsx` | Excel workbook with one worksheet per MOFA view, containing the features x factors weight matrix with feature IDs as row names and factor names as column names. |
| `feature_weights.rds` | Named list of view-specific feature-weight matrices, with names equal to the MOFA view names. |
| `sample_metadata.xlsx` | Sample metadata table. The worksheet is named using `study_name`. |
| `sample_metadata.rds` | RDS version of the sample metadata table. |
| `logs/sessionInfo.txt` | R session information for reproducibility. |
| `logs/mofa_analysis.R` | Exact analysis script used for the run. |
| `logs/mofa_parameters.R` | Exact parameter file used for the run. |
| `logs/file_readers.R` | Reader functions used for the run. |
| `logs/methods.R` | Shared helper functions used for the run. |
| `README.txt` | Short run-specific description of what was analyzed. |

For the feature-weight workbook, each worksheet is named after its MOFA view,
for example `methylation`, `proteomics`, or `metabolomics`. The RDS object has
the equivalent structure:

```r
feature_weights[["proteomics"]]
feature_weights[["methylation"]]
```

### Outputs from a metadata-only run

When no omics file is configured, the script creates a timestamped folder with
`_metadata_only` in its name. It writes:

- `sample_metadata.xlsx`
- `sample_metadata.rds`
- `README.txt`
- `logs/sessionInfo.txt`
- `logs/mofa_analysis.R`
- `logs/mofa_parameters.R`
- `logs/file_readers.R`
- `logs/methods.R`

No MOFA model, factor values, feature weights, or training log are generated
because no model is fitted.

# 2. Regress metadata variables on MOFA factors with `regress_mofa_factors.R`

The second step tests associations between selected metadata variables and the
MOFA factors. The script automatically runs **both regression directions** for
every combination of metadata variable `Y` and covariate set `X`.

The two model directions are:

```text
Y ~ Factor_i + X_covariates
```

and

```text
Factor_i ~ Y + X_covariates
```

The first treats the metadata variable as the outcome. The second treats the
MOFA factor as the outcome and the metadata variable as the explanatory
variable of interest.

## Configure `regression_parameters.R`

Exactly one of these two input modes must be selected.

### Mode 1: use the fitted MOFA model

Set:

```r
mofa_res_fname <- "path/to/model.hdf5"
sample_metadata_fname <- NULL
factor_values_fname <- NULL
```

The script loads the model with `MOFA2::load_model()` and extracts the sample
metadata and factor values directly from the MOFA object.

### Mode 2: use persisted metadata and factor-value files

Set:

```r
mofa_res_fname <- NULL
sample_metadata_fname <- "path/to/sample_metadata.xlsx"
factor_values_fname <- "path/to/factor_values.xlsx"
```

or use the corresponding `.rds` files. The two files may independently be
Excel or RDS files.

The second mode is primarily intended for use by external studies where factor
values are estimated by either applying a trained GLM generated from a MOFA 
model or via direct MOFA factor projection. Accordingly, in this mode,
although the metadata and factor-value tables may be constructed manually 
(provided they follow the required formatting conventions), they will typically
be generated by one of the following scripts:

```r
glm_factor_estimate.R
project_mofa_factors.R
```

The metadata table must contain the column specified by `sample_id_col`, as well
as every column listed in `categorical`. In the factor-values matrix, the first
column must contain the sample IDs.

Before running the regressions, the script verifies that every sample ID in the
factor-values matrix is present in the configured sample-ID column of the
metadata table. If any IDs are missing, execution stops with an informative
error listing them.

For files generated by the external-study workflows, set `factor_values_fname`
to `glm_factor_values.xlsx`/`.rds` or
`mofa_factor_projection.xlsx`/`.rds`, as appropriate. Direct regression input
requires one samples x factors matrix. If `combine_views = "none"` was used,
select or combine a view-specific matrix before running the regression script.

The main regression settings are:

| Variable | Meaning |
|---|---|
| `mofa_res_fname` | Path to `model.hdf5` for input option 1, or `NULL` for option 2. |
| `sample_metadata_fname` | Path to persisted sample metadata (`.rds` or `.xlsx`) for option 2. |
| `factor_values_fname` | Path to persisted factor values (`.rds` or `.xlsx`) for option 2. |
| `sample_id_col` | Metadata column containing sample IDs. The default is `"sample_id"`, but the name can be changed for a different metadata schema. |
| `Y` | Character vector of one or more metadata variables to test. |
| `X` | List of covariate sets. Use `character(0)` for an unadjusted model. |
| `categorical` | Metadata columns that should be modeled as categorical variables. |

The total number of regression batches is:

```r
2 * length(Y) * length(X)
```

### Logistic versus linear regression

In the `Y ~ Factor_i + X` direction, a `Y` variable listed in `categorical`
or stored as a factor/character/logical variable is modeled using logistic
regression. A numeric `Y` that is not listed in `categorical` is modeled using
linear regression.

Therefore, numeric `0/1` outcomes should be included in `categorical` when
logistic regression is desired.

In the `Factor_i ~ Y + X` direction, each MOFA factor is the continuous
outcome, so the model is always linear regression. Categorical explanatory
variables are represented using standard R factor contrasts, with omnibus
term-level tests used when appropriate.

## Outputs from `regress_mofa_factors.R`

The script generates two PDF reports for each `Y`/`X` combination:

```text
regression_analysis/outcome_as_Y/mofa_factor_regression_report_Y=<Y>_covar=<X>.pdf
regression_analysis/factors_as_Y/mofa_factor_regression_report_Y=<Y>_covar=<X>.pdf
```

For example:

```text
regression_analysis/outcome_as_Y/mofa_factor_regression_report_Y=CR_covar=none.pdf
regression_analysis/factors_as_Y/mofa_factor_regression_report_Y=CR_covar=female+race.pdf
```

In addition to the PDF reports, the script writes two consolidated summary
files in the `regression_analysis` directory:

```text
regression_analysis/regression_model_summary.xlsx
regression_analysis/regression_model_summary.rds
```

The Excel workbook contains three worksheets:

- `Models` contains one row per fitted model. The model ID is the actual
  regression formula, with the specific MOFA factor name substituted for the
  temporary internal factor variable. Models are grouped by regression
  direction, metadata variable, and MOFA factor, and are ordered by increasing
  number of adjustment covariates so that nested or progressively adjusted
  models can be inspected together. Covariate sets of the same size retain
  their order in `X`.
- `Coefficients` contains one row per regression coefficient, including
  adjustment-covariate coefficients and categorical contrasts. For logistic
  regression, both log-odds coefficients and exponentiated odds ratios with
  confidence intervals are retained.
- `Dictionary` describes the columns in the `Models` worksheet.

The `Models` worksheet includes the association estimate, standard error,
confidence interval, p-value, sample size, AIC, BIC, log likelihood, and
regression-type-specific fit/performance measures. Linear models include
R-squared, adjusted R-squared, RMSE, and MAE. Logistic models include McFadden
pseudo-R-squared, AUC, accuracy, Brier score, outcome-event counts, and odds
ratios with confidence intervals where a single association coefficient is
defined. Fields that do not apply to a model type are left blank.

The RDS file stores the same three tables as a named list with entries
`Models`, `Coefficients`, and `Dictionary`. Each entry is stored as a data
matrix. Because each table contains both numeric and character fields, the
matrix representation is character-valued.

The script also creates an in-memory object named `res` containing the fitted
models, separated into `res$outcome_as_Y` and `res$factors_as_Y`. The report
summaries and generated PDF paths are similarly retained in memory for
interactive inspection.

The regression-summary RDS file does not contain the fitted model objects in
`res`; those objects remain available in the R session in which the script is
run.

# 3. Visualize the fitted MOFA model with `mofa_plots_report.Rmd` (optional)

The plotting report reads an existing fitted `model.hdf5` file and generates
exploratory plots. It does not refit MOFA.

## Configure `mofa_plots_report_parameters.R`

Before rendering the report, edit the companion
`mofa_plots_report_parameters.R` file. It should be located in the same
directory as `mofa_plots_report.Rmd`; the report sources this file from its
setup chunk.

The most important setting is:

```r
mofa_res_fname <- "path/to/model.hdf5"
```

The remaining settings control:

- the continuous metadata columns to correlate with factors (`col_names_to_use`) in the factors x metadata correlation plot;
- the subset of "interesting" factors to utilize in detailed factor-specific plots (`sel_facs`);
- the metadata variable used in those more focused plots (`metadata_col`);
- the factor and number of features used for top-weight displays (`top_weight_factor` and `top_weight_nfeatures`);
- the views used for top-weight plots (`top_weight_views`);
- the factor and number of features used for sample heatmaps (`heatmap_factor` and `heatmap_features`);
- the views used for heatmaps (`heatmap_views`);
- the metadata columns used to annotate samples in the heatmap (`heatmap_annotation_samples`); and
- whether dimensionality-reduction plots are generated (`run_dimred`).

The parameter file contains user-level comments describing these settings.
Study- or run-specific plotting choices should normally be made there rather
than by editing `mofa_plots_report.Rmd`.

## Render the report

Keep `mofa_plots_report.Rmd` and `mofa_plots_report_parameters.R` in the same
directory. When the report is rendered, its setup code sources the parameter
file and the remainder of the report uses the variables defined there.

For example, to render HTML:

```r
rmarkdown::render(
  input = "mofa_plots_report.Rmd",
  output_format = "html_document",
  output_file = "mofa_plots_report.html",
  output_dir = "path/to/report_directory"
)
```

For PDF, use `output_format = "pdf_document"` and provide a `.pdf` output
filename. PDF rendering requires a working LaTeX installation.

The report contains model-level visualizations such as data availability,
variance explained, factor correlations, factor/metadata correlations,
selected-factor plots, top feature weights, dimensionality reduction, and
feature heatmaps.

The rendered HTML/PDF file is the principal output of this step.

# 4. Interpret MOFA factors biologically with `mofa_factor_interpretation.R`

After the model has been fitted and the relevant factors have been identified,
`mofa_factor_interpretation.R` can be used to generate biological enrichment
results and an optional factor-by-factor interpretation report.

The script extracts feature loadings from the fitted MOFA model and performs
view-specific enrichment analyses. The current implementation includes:

- signed ranked GO enrichment with `fgsea` for expression, proteomics, and
  methylation, using the full signed loading vector for each factor/view;
- mapping of proteomics features from UniProt identifiers to Entrez IDs;
- mapping of methylation CpG identifiers to Entrez IDs using `missMethyl`;
- interpretation of expression identifiers as Entrez IDs or gene symbols;
- separate `fgsea` analyses for each requested GO ontology (`"BP"`, `"MF"`,
  and/or `"CC"`), so multiple-testing adjustment remains ontology-specific; and
- Metabolon pathway/category over-representation analysis (ORA) for the
  positive and negative high-loading metabolite sets.

For the ranked analyses, the report uses the direction label `ranked`, because
positive and negative loadings are analyzed together in one signed ranking.
The sign of the normalized enrichment score (`NES`) indicates the end of the
ranking toward which a GO term is enriched. Metabolomics ORA uses separate 
`positive` and `negative` selected-feature groups.

## Configure `mofa_factors_interpretation_parameters.R`

Before running the script, edit the companion
`mofa_factors_interpretation_parameters.R` file. It must be located in the same
directory as `mofa_factor_interpretation.R`.

The user-configurable settings are:

| Variable | Meaning |
|---|---|
| `mofa_res_fname` | Full path to the `model.hdf5` produced by `mofa_analysis.R`. |
| `metab_ann_fn` | Full path to the Metabolon annotation CSV file. This may be `NULL` when the MOFA model does not contain a `metabolomics` view. |
| `output_dir` | Parent directory in which the `mofa_factor_interpretation` results folder is created. If `NULL`, that folder is created alongside `mofa_res_fname`. |
| `views` | MOFA views to analyze. Use `"all"` for all views in the fitted model, or provide a character vector of view names. |
| `factors` | MOFA factors to analyze. Use `"all"` for all factors or specify the factor selection accepted by `MOFA2::get_weights()`. |
| `selection_mode` | High-loading feature-selection rule used for outputs and metabolomics ORA: `"top_n"`, `"zscore"`, or `"quantile"`. |
| `top_n` | Number of positive and negative features selected per factor/view when `selection_mode = "top_n"`. |
| `z_cutoff` | Absolute loading z-score threshold when `selection_mode = "zscore"`. |
| `quantile_cutoff` | Absolute-loading quantile threshold when `selection_mode = "quantile"`. |
| `go_ontology` | One or more GO namespaces to analyze with `fgsea`: `"BP"`, `"MF"`, and/or `"CC"`. Each ontology is analyzed separately and the results are combined. |
| `qvalue_cutoff` | Maximum adjusted p-value retained in the enrichment results. The cutoff is applied uniformly to `fgsea` and metabolomics ORA results. |
| `report_significance_cutoff` | Adjusted p-value threshold used for enrichment terms shown in the PDF report and included in the factor-level interpretation text. |
| `render_report` | Whether to generate the PDF factor-interpretation report. |
| `n_cores` | Number of task-level parallel workers. Use `1` for serial execution; on Windows, parallel execution requires `future` and `future.apply`. |

Other arguments of `run_mofa_integrative_interpretation()` retain the values
specified by the workflow itself rather than being exposed in the parameter
file. In particular, the standard pipeline uses the `base_pdf` report engine
and `fgsea_nproc = 1`; task-level parallelization is controlled by `n_cores`.

The feature-selection settings (`selection_mode`, `top_n`, `z_cutoff`, and
`quantile_cutoff`) determine the high-loading feature sets used for
metabolomics ORA. Expression, proteomics, and methylation enrichment instead
uses the full signed loading ranking for each factor/view.

## Run `mofa_factor_interpretation.R`

As with the other parameterized scripts in the repository, run the script with
`chdir = TRUE` so that it can source its companion parameter file by name. From
an R session:

```r
source(
  "/path_to_mofa_files/mofa_factor_interpretation.R",
  local = TRUE,
  chdir = TRUE
)
```

or from the command line:

```text
Rscript.exe -e "source('/path_to_mofa_files/mofa_factor_interpretation.R', local = TRUE, chdir = TRUE)"
```

The script loads the fitted MOFA model specified by `mofa_res_fname`, expands
`views = "all"` to the views present in the model, reads the Metabolon
annotation when a metabolomics view is present, and then calls
`run_mofa_integrative_interpretation()` with the configured settings.

The default `base_pdf` report engine writes a landscape PDF directly using base
R graphics and therefore does not require Pandoc or LaTeX. The underlying
analysis function also supports an optional `rmarkdown_pdf` engine, which
requires Pandoc and LaTeX, but this setting is not exposed in the standard
parameter file.

## Enrichment-result filtering and output

`qvalue_cutoff` controls which enrichment results are retained in the analysis
outputs. For GO `fgsea`, adjusted p-values are calculated separately within
each requested ontology. The same cutoff is applied to metabolomics ORA after
the enrichment results are combined, so the persisted enrichment table uses a
consistent adjusted-p-value threshold across methods.

`report_significance_cutoff` is a separate, typically stricter threshold used
only for the PDF report and factor-level interpretation text. Thus, the Excel
and RDS results may contain terms that pass `qvalue_cutoff` but do not pass
`report_significance_cutoff` and therefore do not appear in the report.

The output filenames include the requested GO ontology or ontologies. For
example, with:

```r
go_ontology <- c("BP", "MF", "CC")
```

the suffix is `BP_MF_CC`.

The principal outputs in the `mofa_factor_interpretation` results folder are:

| Output | Description |
|---|---|
| `mofa_integrative_factor_report_<ontologies>.pdf` | Optional landscape PDF report containing factor-level interpretation text, variance-explained information, and significant enrichment results. Generated when `render_report = TRUE`. |
| `enrichment_results_<ontologies>.xlsx` | Enrichment results across the selected factors and views. The workbook includes explanatory comments on recognized column headers. Results are reported using a consolidated schema shared by `fgsea` and metabolomics ORA. |
| `mofa_integrative_results_<ontologies>.rds` | Condensed R results object containing only `enrichment_results`, `interpretations`, and `parameters`. |

The full analysis result list, including loadings, selected features, top
features, and variance explained, remains available in the in-memory object
`res` during an interactive run. The standard scripted workflow does not write
the optional loading/feature Excel files because `save_loadings` retains its
default value of `FALSE`.

# 5. Train external factor predictors with `mofa_train_glm.R`

The fifth step assesses whether the learned MOFA factors can be predicted
from individual omics views. This is useful because the resulting predictor
objects can subsequently be applied to external studies that have comparable
molecular measurements.

## Configure `mofa_train_glm_parameters.R`

The most important settings are:

| Variable | Meaning |
|---|---|
| `mofa_res_fname` | Full path to the `model.hdf5` produced by `mofa_analysis.R`. |
| `factors` | Factor names to train predictors for, or `"all"`. |
| `views` | MOFA views to train predictors for, or `"all"`. |
| `outcome_type` | Type of factor outcome. For MOFA factor values this is normally `"continuous"`. |
| `n_repeats` | Number of repeated internal train/test splits used to evaluate prediction performance. |
| `train_frac` | Fraction of samples used for training in each split. |
| `seed` | Optional seed for reproducibility. |
| `alpha` | Elastic-net mixing parameter; `1` corresponds to LASSO. |
| `lambda_choice` | Penalty selected from each `cv.glmnet()` fit, normally `"lambda.min"` or `"lambda.1se"`. |
| `standardize` | Whether `glmnet` standardizes predictors internally. |
| `nfolds` | Number of folds used within `cv.glmnet()`. |
| `min_samples` | Minimum number of matched samples required for a view/factor predictor. |
| `min_features` | Minimum number of usable features required for a predictor. |
| `use_parallel` | Whether independent view/factor predictor jobs are run in parallel. |
| `n_workers` | Number of parallel worker processes. |
| `verbose` | Whether progress messages are printed. |
| `report_file` | Name of the PDF training report. |
| `include_top_features` | Whether top nonzero features are also exported in the report tables. |
| `top_n_features` | Number of top features to include when `include_top_features` is `TRUE`. |

The fitted MOFA model is loaded directly, and `mofa_train_glm.R` retrieves both
the preprocessed training views and factor values from that model. Thus, the
predictors are trained on the same matrices used to fit MOFA.

## What `mofa_train_glm.R` does

For each selected view/factor combination, the script:

1. Matches the view's samples to the corresponding factor values.
2. Removes all-`NA` sample columns that may have been added to a MOFA view to
   represent missing measurements.
3. Removes unusable features.
4. Runs `n_repeats` internal train/test validations. Each validation uses a `train_frac` 
   fraction of the samples for training and the remaining samples for validation. Within 
   each training subset, `cv.glmnet()` performs `nfolds`-fold cross-validation to select 
   the penalization parameter, after which prediction performance is evaluated on the 
   held-out validation samples.
5. Fits the **final** penalized regression model using all available matched samples. The 
   final model uses `cv.glmnet()` with `nfolds`-fold cross-validation to select the penalization 
   parameter.
6. Records the selected features and coefficients, as well as prediction
   performance metrics.
7. Saves all view/factor models in a single `MofaFactorPredictor` RDS object.

The predictor object is therefore a learned mapping from one molecular view to
one or more MOFA factors. The internal validation results are especially useful
when deciding whether a factor can be estimated reliably from a particular
view.

## Outputs from `mofa_train_glm.R`

Outputs are written to a `glm_model` folder alongside the fitted MOFA model:

| Output | Description |
|---|---|
| `mofa_glm_predictor.rds` | Main `MofaFactorPredictor` object containing settings, validation results, final penalized models, selected features, and coefficients. This is the key file needed by `glm_factor_estimate.R`. |
| `training_metrics_summary.csv` | Summary table of internal validation performance. |
| `training_metrics_long.csv` | Long-format validation metrics used by the report. |
| `nonzero_coefficients.csv` | Nonzero coefficients for the final predictors. |
| `top_nonzero_coefficients.csv` | Top nonzero coefficients, generated when `include_top_features = TRUE`. |
| `mofa_factor_predictor_training_report.pdf` | PDF report summarizing validation performance and predictor behavior; the filename is controlled by `report_file`. |
| `README.txt` | Run-specific description of the training results and logging files. |
| `logs/sessionInfo.txt` | R session information. |
| `logs/mofa_train_glm.R` | Copy of the training script used. |
| `logs/mofa_train_glm_parameters.R` | Copy of the parameter file used. |

# External factor estimation for other studies

The repository contains two workflows specifically intended for external
studies. These studies do **not** refit the MOFA model. Instead, they use the
model artifacts learned in the primary study to estimate the corresponding MOFA
factor values for new samples.

There are two complementary approaches:

1. **GLM/LASSO factor estimation** with `glm_factor_estimate.R`, using the
   `mofa_glm_predictor.rds` object trained by `mofa_train_glm.R`.
2. **Direct factor projection** with `project_mofa_factors.R`, using the
   `feature_weights.rds` object exported by `mofa_analysis.R`.

Both approaches require the external study's data to be made comparable to the
data used to learn the original MOFA factors. In particular, the external
omics matrices must use compatible feature identifiers and preprocessing,
including the relevant normalization, transformations, change-score definition,
and scaling conventions.

## A. Estimate factors using trained GLM/LASSO predictors

### What this approach uses

`mofa_train_glm.R` learns a separate penalized regression predictor for each
view/factor pair. For example, it can learn a predictor of `Factor 1` from the
proteomics view and another predictor of `Factor 1` from the metabolomics view.

The trained object is saved as:

```text
glm_model/mofa_glm_predictor.rds
```

`glm_factor_estimate.R` loads this predictor object and applies the trained
models to the external study's omics data.

This is the preferred path when the training analysis has established that
particular factors can be predicted adequately from particular molecular views,
as assessed by the internal cross-validation metrics generated by
`mofa_train_glm.R`.

### Configure `glm_factor_estimate_parameters.R`

The key settings are (most are analogous to the ones used in `mofa_parameters.R`):

| Variable | Meaning |
|---|---|
| `glm_predictor_fn` | Path to `mofa_glm_predictor.rds` produced by `mofa_train_glm.R`. |
| `output_dir` | Optional output root. If `NULL`, results are written beside the predictor object. |
| `study_name` | Study name used in generated documentation and Excel sheet names. |
| `subject_id_col` | Column name in the clinical/phenotype file containing subject identifiers. |
| `treatment_col` | Column name in the clinical/phenotype file containing treatment codes. |
| `timepoint_col` | Column name in the clinical/phenotype file containing timepoint codes. |
| `timepoint_labels` | Mapping between external clinical timepoint codes and labels; includes the `baseline` definition. |
| `treatment_labels` | Mapping between external clinical treatment codes and labels. |
| `timepoint` | External-study visit to estimate factor values for. |
| `use_change_scores` | Whether external data should be converted to target-minus-baseline change scores before prediction. This must match the definition used in the training workflow. |
| `scale_views` | Whether external features should be z-score scaled before prediction. This must match the scale assumed by the trained predictor. |
| `methy_fn`, `prot_fn`, `metab_data_fn`, `metab_annot_fn`, `clin_fn` | External-study input files. Comment out an omics file to exclude that view. |
| `omics_input_defs` | External view definitions, including barcode columns and reader functions. The view names should match those used in the trained predictor. |
| `views` | Trained views to use, or `"all"`. |
| `combine_views` | How to combine view-specific predictions for the same factor: `"mean"`, `"median"`, `"best_cv"`, or `"none"`. |
| `metric_for_best_cv` | CV metric used when selecting the best view. |
| `lambda_choice` | Optional override of the lambda value stored in each trained model. |

The view-specific external input data are processed using the same general
clinical/timepoint/change-score conventions as the main MOFA workflow. Feature
identifiers must match the identifiers used by the trained predictors.

### Outputs from `glm_factor_estimate.R`

The script creates a descriptively-named folder under the configured output root.
Theprincipal outputs are:

| Output | Description |
|---|---|
| `glm_factor_values.xlsx` | Estimated factor values. With combined views, this is a samples x factors matrix; with `combine_views = "none"`, the workbook contains one worksheet per view. |
| `glm_factor_values.rds` | RDS version of the corresponding factor estimates. With combined views, this is a samples x factors matrix; with `combine_views = "none"`, it is a named list of view-specific samples x factors matrices. |
| `sample_metadata.xlsx` | Sample metadata for the configured target visit, constructed with the same sample-ID convention as `mofa_analysis.R`. The worksheet is named using `study_name`. |
| `sample_metadata.rds` | RDS version of the sample metadata table. |
| `README.txt` | Run-specific description of the external scoring operation. |
| `logs/sessionInfo.txt` | R session information. |
| `logs/glm_factor_estimate.R` | Copy of the estimation script used. |
| `logs/glm_factor_estimate_parameters.R` | Copy of the parameter file used. |
| `logs/file_readers.R` | Reader functions used for the external data. |
| `logs/methods.R` | Shared helper functions used for external-data preparation. |

These estimated factor values can then be analyzed in the external study using
its own clinical metadata and statistical workflow.

## B. Project factors directly using saved MOFA feature weights

`project_mofa_factors.R` provides a second external-scoring strategy that does
not use the GLM predictor object. Instead, it reads the `feature_weights.rds`
file created by `mofa_analysis.R` and projects the external study's molecular
profiles onto the learned view-specific MOFA feature weights.

The default projection is a regularized least-squares calculation. It should be
understood as a practical factor-projection method, not as exact MOFA2
variational posterior inference for new samples.

### Configure `project_mofa_factors_parameters.R`

The key settings are (most are analogous to the ones used in `mofa_parameters.R`):

| Variable | Meaning |
|---|---|
| `feature_weights_fn` | Path to `feature_weights.rds` from the original MOFA run. |
| `glm_predictor_fn` | Optional path to `mofa_glm_predictor.rds`; needed only for `best_cv` or `weighted_cv` view combination. |
| `output_dir` | Optional output root. If `NULL`, results are written beside `feature_weights.rds`. |
| `study_name` | Study name used in generated documentation and Excel sheet names. |
| `subject_id_col` | Column name in the clinical/phenotype file containing subject identifiers. |
| `treatment_col` | EColumn name in the clinical/phenotype file containing treatment codes. |
| `timepoint_col` | Column name in the clinical/phenotype file containing timepoint codes. |
| `timepoint_labels` | Mapping between external timepoint codes and labels. |
| `treatment_labels` | Mapping between external treatment codes and labels. |
| `timepoint` | External-study visit to project. |
| `use_change_scores` | Whether external omics data should be converted to target-minus-baseline change scores. This must match the original MOFA preprocessing. |
| `scale_views` | Whether external features are z-score scaled before projection. This must match the original preprocessing scale. |
| `methy_fn`, `prot_fn`, `metab_data_fn`, `metab_annot_fn`, `clin_fn` | External input files. |
| `omics_input_defs` | Definitions of the external omics views. View names must match the names in `feature_weights.rds`. |
| `combine_views` | How view-specific projections are combined: `"mean"`, `"median"`, `"best_cv"`, `"weighted_cv"`, or `"none"`. |
| `metric` | Internal-CV metric used by `best_cv` or `weighted_cv`. |
| `ridge_lambda` | Ridge regularization strength used in the feature-weight projection. Positive values generally improve numerical stability. |
| `verbose` | Whether projection progress messages are printed. |

The direct-projection workflow is particularly sensitive to preprocessing
compatibility. The external data should be put on the same scale as the data
used to learn the original MOFA model, including the same transformations,
normalization, feature identifiers, change-score definition, and any required
z-scoring.

### Outputs from `project_mofa_factors.R`

The script creates a descriptively-named projection-result folder. The principal
outputs are:

| Output | Description |
|---|---|
| `mofa_factor_projection.xlsx` | Projected factor values. With combined views, one worksheet contains the samples x factors matrix; with `combine_views = "none"`, one worksheet is written per view. |
| `mofa_factor_projection.rds` | RDS version of the corresponding factor projections. With combined views, this is a samples x factors matrix; with `combine_views = "none"`, it is a named list of view-specific samples x factors matrices. |
| `sample_metadata.xlsx` | Sample metadata for the configured target visit, constructed with the same sample-ID convention as `mofa_analysis.R`. The worksheet is named using `study_name`. |
| `sample_metadata.rds` | RDS version of the sample metadata table. |
| `README.txt` | Run-specific description of the projection. |
| `logs/sessionInfo.txt` | R session information. |
| `logs/project_mofa_factors.R` | Copy of the projection script used. |
| `logs/project_mofa_factors_parameters.R` | Copy of the projection parameter file used. |
| `logs/file_readers.R` | Reader functions used for the external data. |
| `logs/methods.R` | Shared helper functions used for external-data preparation. |

# Recommended end-to-end workflow

For a new study, the intended workflow is:

1. **Configure and run `mofa_analysis.R`.** Edit `mofa_parameters.R` to define
   the study-specific clinical columns, input files, target visit, and omics
   views. The result folder contains the fitted `model.hdf5`, factor values,
   feature weights, and sample metadata. Its `logs` subfolder contains the
   model-training log, session information, and source/configuration files used
   for the run.
2. **Configure and run `regress_mofa_factors.R`.** Edit
   `regression_parameters.R` to select either the fitted MOFA model or 
   persisted metadata/factor-value files, set `sample_id_col`, and specify `Y`,
   `X`, and `categorical`. Review both sets of regression PDF reports and the
   consolidated `regression_model_summary.xlsx`/`.rds` model summaries.
3. **Optionally render `mofa_plots_report.Rmd`.** Edit
   `mofa_plots_report_parameters.R` to point to the fitted `model.hdf5` and
   select the factors, metadata variables, views, feature-weight, heatmap, and
   dimensionality-reduction settings of interest.
4. **Configure and run `mofa_factor_interpretation.R`.** Edit
   `mofa_factors_interpretation_parameters.R` to identify the fitted MOFA model,
   select the factors/views and GO ontologies to analyze, set enrichment/report
   thresholds, and configure parallelization. Review the generated enrichment
   workbook and optional PDF interpretation report.
5. **Configure and run `mofa_train_glm.R`.** Edit
   `mofa_train_glm_parameters.R` to point to the `model.hdf5` file and select
   the factors/views to evaluate. Use the training report and internal-CV
   metrics to assess view-specific predictability and retain the resulting
   `mofa_glm_predictor.rds`.

## External-study workflow

Once the primary MOFA study has produced `mofa_glm_predictor.rds` and/or
`feature_weights.rds`, an external study can estimate the same learned factor
values without refitting MOFA:

1. Obtain the model artifact from the primary study (`mofa_glm_predictor.rds`
   for GLM/LASSO estimation, or `feature_weights.rds` for direct projection).
2. Prepare the external omics data using the same feature identifiers and the
   same preprocessing conventions used by the primary study.
3. Configure either `glm_factor_estimate_parameters.R` or
   `project_mofa_factors_parameters.R` with the external clinical and omics
   files.
4. Run the corresponding estimation script.
5. Use the resulting samples x factors matrix in the external study's own
   statistical analyses, while retaining the generated `logs` folder and
   `README.txt` for provenance.
6. **Optionally run `regress_mofa_factors.R`** to test associations between the
   externally estimated factors and the external study's metadata. In
   `regression_parameters.R`, select the persisted-file input mode, set
   `sample_metadata_fname` to the generated `sample_metadata.xlsx` or `.rds`
   file, and set `factor_values_fname` to the generated
   `glm_factor_values.xlsx`/`.rds` or
   `mofa_factor_projection.xlsx`/`.rds` file. This direct input requires a
   single samples x factors matrix; view-specific outputs produced with
   `combine_views = "none"` must first be selected or combined.

The two methods answer slightly different operational needs. GLM/LASSO
estimation applies a supervised predictor learned from the original MOFA
factor values and provides CV-based view-selection options. Direct projection
uses the original MOFA feature weights and a regularized projection formula.
Neither method refits the original MOFA model, and both require careful
harmonization of the external data with the original training scale.
