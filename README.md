# MOFA Study Workflows

This repository contains code for running and exploring **Multi-Omics Factor Analysis (MOFA)** models across multiple studies.

The repository is organized with **one folder per study**. Each study folder contains a customized MOFA workflow tailored to that study's input data, metadata, modeling choices, and downstream analyses.

## Repository organization

A typical repository layout is:

```text
.
├── CALERIE/
│   ├── README.md
│   ├── resources.R
│   ├── MOFA_analysis.R
│   ├── regress_mofa_factors.R
│   └── ...
├── AnotherStudy/
│   ├── README.md
│   ├── resources.R
│   ├── MOFA_analysis.R
│   ├── regress_mofa_factors.R
│   └── ...
└── README.md
```

Each study folder is intended to be relatively self-contained. The scripts inside a study folder define the data sources, run the MOFA model, and perform study-specific downstream analyses.

## Per-study customization

Although the overall analysis pattern is expected to be similar across studies, the implementation is intentionally customized within each study folder. This allows each workflow to reflect study-specific details such as:

- available omics data matrices
- sample metadata and outcome variables
- preprocessing decisions
- MOFA model settings
- covariates used in downstream analyses
- study-specific plots, summaries, and result exploration

This organization favors clarity and reproducibility over forcing all studies into a single generalized pipeline.

## Typical contents of a study folder

Most study folders contain some combination of the following files:

| File | Purpose |
|---|---|
| `README.md` | Study-specific documentation describing inputs, outputs, and usage. |
| `resources.R` | Defines file paths and other study-specific resources. |
| `MOFA_analysis.R` | Runs the MOFA analysis for that study. |
| `regress_mofa_factors.R` | Performs regression analyses relating MOFA factors to metadata outcomes. |
| Additional scripts | Optional study-specific scripts for visualization, diagnostics, or downstream analyses. |

The exact filenames may vary slightly by study. Consult the `README.md` inside each study folder for details.

## General workflow

The usual workflow within a study folder is:

1. Edit `resources.R` to point to the relevant input files.
2. Run the study-specific MOFA analysis script.
3. Inspect the fitted MOFA model and generated outputs.
4. Run downstream analysis scripts, such as factor regression or visualization scripts.
5. Review the study-specific `README.md` for additional instructions and assumptions.

## Notes on reproducibility

Because each study folder contains its own customized scripts and documentation, users should treat each folder as a separate analysis unit. Study-specific decisions should be documented in that folder's `README.md` and reflected explicitly in the code.

This structure makes it easier to rerun, audit, and adapt each MOFA workflow without obscuring the assumptions that are unique to a given study.
