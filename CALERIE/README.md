# CALERIE MOFA analysis workflow

This folder contains R scripts for running a MOFA analysis on CALERIE omics data and then testing associations between the fitted MOFA factors and sample-level metadata variables.

The workflow has two main stages:

1. **Fit a MOFA model** from one or more omics data matrices.
2. **Regress a selected metadata outcome** on each fitted MOFA factor, optionally adjusting for metadata covariates.

## Files in this folder

| File | Purpose |
|---|---|
| `resources.R` | Defines local file paths for the CALERIE input files. Edit this file to specify the file locations for your local compute environment. |
| `MOFA_analysis_calerie.R` | Main MOFA-fitting script. It reads the omics and phenotype files, preprocesses the selected omics views, fits a MOFA model, and writes the fitted model to `model.hdf5`.  |
| `regress_mofa_factors.R` | Loads a fitted `model.hdf5` file and runs one regression per MOFA factor. The outcome, adjustment covariates, and categorical variables are user-selected near the top of the script. |
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
```

### Configuration options

| Variable | Meaning |
|---|---|
| `code_dir` | Folder containing this code, including `resources.R` and `methods.R`. |
| `output_dir` | Folder where MOFA result folders will be created. |
| `timepoint` | Follow-up visit to analyze, `12` or `24` months. Samples are selected using `fumo == timepoint`. |
| `use_change_scores` | If `TRUE`, each omics value is converted to a change score: target visit minus baseline. If `FALSE`, the script uses target-visit values directly. |
| `use_groups` | If `TRUE`, the MOFA model treats the AL and CR arms as separate MOFA groups. If `FALSE`, all samples are analyzed together without MOFA grouping. |

### What the MOFA script does

`MOFA_analysis_calerie.R` performs the following steps:

1. Sources `resources.R` and `methods.R`.
2. Reads the phenotype/visit data.
3. Reads whichever omics files are defined in `resources.R`.
4. For each available omics view, either:
   - computes change scores as `timepoint - baseline`, or
   - extracts values from the selected follow-up timepoint.
5. Removes zero-variance features.
6. Drops omics views with no samples at the requested timepoint.
7. Optionally infers AL/CR group labels from sample IDs.
8. Filters each omics view to retain highly variable features.
9. Adds missing sample columns filled with `NA`, as required by MOFA when views have different sample coverage.
10. Creates and prepares the MOFA object.
11. Adds sample metadata to the MOFA object using `samples_metadata(mo) <- meta_d`.
12. Fits the model with `run_mofa()`.
13. Writes run outputs to a timestamped result folder.

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

## Stage 2: regress metadata outcomes on MOFA factors

After fitting a model, use `regress_mofa_factors.R` to test whether a selected metadata outcome is associated with each MOFA factor.

Open `regress_mofa_factors.R` and edit the user-specified inputs near the top of the file.

```r
mofa_res_fname <- "D:/R_sessions/PROSPR/MOFA/.../model.hdf5"

X <- c("bmistrat")
Y <- "CR"

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
| `Y` | Name of the metadata column to use as the dependent variable. This will typically be `"CR"`|
| `categorical` | Metadata columns that should be treated as categorical variables. Include `Y` here when the outcome should be modeled using logistic regression. |

### Model fitted by the regression script

The script loads the fitted MOFA model:

```r
mo.tr <- MOFA2::load_model(mofa_res_fname)
```

Then it runs:

```r
res <- runRegression(mo.tr, X, Y, categorical)
```

For each MOFA factor, the script fits one model of the form:

```r
Y ~ Factor_i + X_covariates
```

Internally, the current factor is stored in a temporary column named `F`, so the formula is constructed as:

```r
Y ~ F + covariates
```

### Logistic versus linear regression

`runRegression()` chooses the model type from the outcome variable `Y`:

| Outcome type | Model used |
|---|---|
| `Y` listed in `categorical`, or otherwise stored as a factor/character/logical variable | Logistic regression using `glm(..., family = binomial())` |
| Numeric `Y` not listed in `categorical` | Linear regression using `lm()` |

Important: a numeric `0/1` outcome is treated as continuous unless its column name is included in `categorical`. 

For logistic regression, the outcome must have exactly two observed classes after missing values are removed for a given factor model. The first factor level is used as the reference level by R. Use `relevel()` in the metadata before modeling if you need a specific reference category.

### Regression output

The object `res` is a named list with one entry per MOFA factor:

```r
res[["Factor 1"]]
res[["Factor 2"]]
```

Each entry is either a `glm` object for logistic regression or an `lm` object for linear regression.

Useful inspection commands:

```r
attr(res, "model_type")
attr(res, "formula")
summary(res[["Factor 1"]])
lapply(res, summary)
```

For logistic regression, coefficients are on the log-odds scale. Exponentiate coefficients to obtain odds ratios:

```r
sapply(res, function(model) exp(coef(model)))
```

## Recommended workflow

1. Edit `resources.R` so all local file paths are correct.
2. Edit the configuration block in `MOFA_analysis_calerie.R`:
   - `code_dir`
   - `output_dir`
   - `timepoint`
   - `use_change_scores`
   - `use_groups`
3. Run the MOFA analysis.
4. Locate the resulting `model.hdf5` file in the timestamped output folder.
5. Edit `mofa_res_fname`, `X`, `Y`, and `categorical` in `regress_mofa_factors.R`.
6. Run the regression script.
7. Inspect the returned `res` object to review factor-outcome associations.
