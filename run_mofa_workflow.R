# *****************************************************************************
# Run the MOFA analysis, factor-regression, and GLM-training workflows for
# multiple combinations of omics views.
#
# Before running this script:
#   1. Configure mofa_parameters.R, regression_parameters.R, and
#      mofa_train_glm_parameters.R as for their respective standalone scripts.
#      Define the input file paths for every omics view that will be requested.
#   2. Edit configure_runs.txt to list one combination per line, using the
#      case-sensitive view names defined in mofa_parameters.R::omics_input_defs.
#      For example:
#          methylation, proteomics, metabolomics
#          proteomics, metabolomics
#      Blank lines and lines beginning with # are ignored.
#
# For each combination, this script runs mofa_analysis.R, followed by
# regress_mofa_factors.R and mofa_train_glm.R. It passes the newly generated
# model path to both downstream scripts; the parameter files are not edited.
# A CSV summary records the status of each step and any errors. See README.md
# for details, including how failures are handled.
#
# Run from an R session:
#   source("path_to_mofa_code/run_mofa_workflow.R", local = TRUE, chdir = TRUE)
#
# Or from the command line:
#   Rscript.exe -e "source('path_to_mofa_code/run_mofa_workflow.R',
#                         local = TRUE, chdir = TRUE)"
#
# chdir = TRUE temporarily sets the working directory to this script's folder,
# allowing the supporting scripts and the default configuration file to be
# located there. Relative input/output paths use that working directory too.
# *****************************************************************************

#' Run the three-step MOFA workflow for configured omics combinations
#'
#' @param config_file Path to a text file containing one comma-separated list
#'   of omics_input_defs view names per line. Defaults to configure_runs.txt in
#'   the working directory.
#' @return Invisibly returns a data frame with one row per requested run,
#'   including the input line, views, model path, step statuses, and errors.
#'   The same information is written incrementally to a CSV file in output_dir.
#' @details Configuration and input-file checks occur before model training.
#'   If MOFA fails for a combination, its downstream steps are skipped. A
#'   regression failure does not prevent GLM training for that combination.
#'   After a step failure, the script proceeds to the next applicable step or
#'   combination. A status of completed means the script returned without error.
runMofaWorkflow <- function(config_file = "configure_runs.txt") {
  # Read the entire run configuration before starting any model training. The
  # original line numbers are retained for diagnostics and the output summary.
  if (!file.exists(config_file)) stop("Cannot find configuration: ", config_file)
  lines <- readLines(config_file, warn = FALSE, encoding = "UTF-8")
  # Remove a possible UTF-8 byte-order mark, and ignore comments/blank lines.
  lines <- sub("^\ufeff", "", lines)
  indices <- which(nzchar(trimws(lines)) & !grepl("^\\s*#", lines))
  if (!length(indices)) stop("No view combinations found in ", config_file)
  # Parse each remaining line into a character vector of view names. Trim
  # whitespace and reject empty or duplicate entries within a combination.
  runs <- lapply(indices, function(i) {
    line <- trimws(lines[[i]])
    # strsplit drops trailing empty fields, so check for them explicitly.
    if (grepl(",\\s*$", line)) stop("Empty view name on line ", i)
    fields <- trimws(strsplit(line, ",", fixed = TRUE)[[1L]])
    if (any(!nzchar(fields))) stop("Empty view name on line ", i)
    if (anyDuplicated(fields)) stop("Duplicate view name on line ", i)
    fields
  })
  # Sort only for comparison: combinations with the same views in different
  # orders count as duplicates. The original order is retained for execution.
  run_keys <- vapply(runs, function(x) paste(sort(x), collapse = "\r"), "")
  if (anyDuplicated(run_keys)) stop("The same combination appears more than once.")

  # Source the same settings as the analysis script, without reading omics data.
  # Use a fresh environment so values in the interactive workspace do not count
  # as file paths that are missing from the parameters file.
  settings <- new.env(parent = globalenv())
  sys.source("file_readers.R", envir = settings)
  sys.source("mofa_parameters.R", envir = settings)
  sys.source("methods.R", envir = settings)
  # Check definitions and requested names against the parameter file. Validate
  # every requested combination before committing time to any MOFA run.
  settings$validateOmicsInputDefinitions(settings$omics_input_defs)
  view_names <- names(settings$omics_input_defs)
  if (anyNA(view_names) || any(!nzchar(view_names)) || anyDuplicated(view_names))
    stop("omics_input_defs must have unique, non-empty names.")
  unknown <- setdiff(unique(unlist(runs)), view_names)
  if (length(unknown)) stop("Unknown view name(s): ", paste(unknown, collapse = ", "))
  # Each view definition names a variable containing its input path. Require
  # that variable to be defined in the parameter environment itself, and check
  # that it names an existing file. Custom readers may require additional files,
  # which must be configured separately and are not checked by this loop.
  for (view in unique(unlist(runs))) {
    var <- settings$omics_input_defs[[view]]$file_var
    if (!is.character(var) || length(var) != 1L || is.na(var) || !nzchar(var) ||
        !exists(var, envir = settings, inherits = FALSE))
      stop("No file path is defined in mofa_parameters.R for view: ", view)
    path <- get(var, envir = settings, inherits = FALSE)
    if (!is.character(path) || length(path) != 1L || is.na(path) ||
        !nzchar(trimws(path)) || !file.exists(path) || dir.exists(path))
      stop("Invalid or missing input file for view '", view, "': ", paste(path, collapse = ", "))
  }
  # Create the common output directory if necessary. Each individual workflow
  # still creates its usual results and logs under this directory.
  if (!is.character(settings$output_dir) || length(settings$output_dir) != 1L ||
      is.na(settings$output_dir) || !nzchar(trimws(settings$output_dir)))
    stop("output_dir must be a non-empty path.")
  dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(settings$output_dir)) stop("Cannot create output_dir.")
  # tempfile avoids summary collisions when batches are started close together.
  summary_file <- tempfile("workflow_summary_", tmpdir = settings$output_dir, fileext = ".csv")
  # Initialize one row per combination. Persist the summary now and after each
  # step, so completed work and errors remain visible if a later step fails.
  summary <- data.frame(line = indices, views = vapply(runs, paste, "", collapse = ", "),
                        model_file = "", mofa = "pending", regression = "pending",
                        glm = "pending", error = "", stringsAsFactors = FALSE)
  save_summary <- function() write.csv(summary, summary_file, row.names = FALSE, na = "")
  save_summary()
  message("Workflow summary: ", normalizePath(summary_file, winslash = "/"))

  # Process combinations sequentially. Each script is evaluated in a fresh R
  # environment to avoid carrying variables from one step or run into another.
  # These environments still inherit from the global environment; this is not
  # a separate R process and does not isolate process-wide settings or packages.
  for (i in seq_along(runs)) {
    message("Starting run ", i, "/", length(runs), ": ", summary$views[i])
    analysis <- new.env(parent = globalenv())
    # The analysis script recognizes this private variable and selects exactly
    # these views instead of using all views with configured paths. On success,
    # it places the path to its model.hdf5 file in .mofa_run_model.
    analysis$.mofa_run_views <- runs[[i]]
    err <- tryCatch({
      sys.source("mofa_analysis.R", envir = analysis)
      if (!exists(".mofa_run_model", envir = analysis, inherits = FALSE) ||
          !file.exists(analysis$.mofa_run_model)) stop("MOFA model was not produced.")
      NULL
    }, error = function(e) conditionMessage(e))
    # A missing or failed MOFA model prevents both downstream steps. Record
    # the error and move directly to the next configured combination.
    if (!is.null(err)) {
      summary$mofa[i] <- "failed"
      summary$regression[i] <- summary$glm[i] <- "skipped"
      summary$error[i] <- paste("MOFA:", err)
      message(summary$error[i]); save_summary(); next
    }
    summary$mofa[i] <- "completed"
    summary$model_file[i] <- analysis$.mofa_run_model
    save_summary()
    # Pass the exact model path to each downstream script through a private
    # override variable. Each script applies the override after sourcing its
    # parameters, without changing the parameter file on disk. The regression
    # script also selects its fitted-model input mode for this run.
    #
    # Both steps are attempted independently once MOFA has succeeded: a
    # regression error does not stop GLM training or subsequent combinations.
    for (step in c("regression", "glm")) {
      step_env <- new.env(parent = globalenv())
      step_env$.mofa_model_override <- analysis$.mofa_run_model
      script <- if (step == "regression") "regress_mofa_factors.R" else "mofa_train_glm.R"
      err <- tryCatch({ sys.source(script, envir = step_env); NULL },
                      error = function(e) conditionMessage(e))
      summary[[step]][i] <- if (is.null(err)) "completed" else "failed"
      # Preserve errors from both downstream steps if both fail.
      if (!is.null(err)) {
        detail <- paste(step, err, sep = ": ")
        summary$error[i] <- paste(c(summary$error[i][nzchar(summary$error[i])], detail),
                                  collapse = " | ")
        message(detail)
      }
      save_summary()
    }
    # Release references to the completed run's environments before the next
    # combination; garbage collection may reclaim large R objects they contain.
    rm(analysis, step_env)
    gc()
  }
  message("Workflow finished. Summary: ", normalizePath(summary_file, winslash = "/"))
  invisible(summary)
}

# Execute the default configuration when this file is sourced.
runMofaWorkflow()
