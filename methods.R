# *****************************************************************************
# Shared helper methods used by the MOFA, GLM factor-estimation, and MOFA
# factor-projection workflows.
#
# The functions in this file implement common resource-file handling, clinical
# configuration validation, sample-ID construction, timepoint/change-score
# extraction, optional view loading, and zero-variance filtering.
# *****************************************************************************

#' Test whether a configured resource-file variable contains a usable path
#'
#' @param var_name Character name of a variable containing a file path.
#' @param envir Environment in which to look up var_name. Defaults to the caller.
#' @return TRUE when the variable exists and contains one non-empty string.
resourceFileSpecified <- function(var_name, envir = parent.frame()) {
  if (!exists(var_name, envir = envir, inherits = TRUE)) return(FALSE)
  value <- get(var_name, envir = envir, inherits = TRUE)
  is.character(value) && length(value) == 1L && !is.na(value) && nzchar(value)
}

#' Retrieve a configured resource-file path
#'
#' @param var_name Character name of a file-path variable.
#' @param required If TRUE, missing/empty variables trigger an error.
#' @param envir Environment in which to look up var_name. Defaults to the caller.
#' @return The configured path, or NULL for an unspecified optional resource.
getResourceFile <- function(var_name, required = FALSE, envir = parent.frame()) {
  if (!resourceFileSpecified(var_name, envir = envir)) {
    if (required) {
      stop("Required variable '", var_name, "' is not defined or is empty.", call. = FALSE)
    }
    return(NULL)
  }
  get(var_name, envir = envir, inherits = TRUE)
}

#' Validate a named code-to-label mapping
validateLabelMap <- function(labels, var_name, required_label = NULL) {
  label_names <- names(labels)
  label_values <- as.character(unname(labels))

  if (!is.atomic(labels) || length(labels) == 0L || is.null(label_names) ||
      length(label_names) != length(labels) || anyNA(label_names) ||
      any(!nzchar(trimws(label_names)))) {
    stop("'", var_name, "' must be a non-empty named vector.", call. = FALSE)
  }
  if (anyDuplicated(label_names)) {
    stop("Names in '", var_name, "' must be unique.", call. = FALSE)
  }
  if (anyNA(label_values) || any(!nzchar(label_values)) || anyDuplicated(label_values)) {
    stop("Values in '", var_name, "' must be non-missing and unique.", call. = FALSE)
  }
  if (!is.null(required_label) && sum(label_names == required_label) != 1L) {
    stop("'", var_name, "' must contain exactly one entry named '",
         required_label, "'.", call. = FALSE)
  }
}

#' Convert configured clinical codes to their display labels
mapValuesToLabels <- function(values, labels, mapping_name) {
  value_keys <- as.character(values)
  label_index <- match(value_keys, as.character(unname(labels)))
  unmapped <- is.na(label_index) & !is.na(value_keys)

  if (any(unmapped)) {
    stop("Unmapped ", mapping_name, " code(s): ",
         paste(unique(value_keys[unmapped]), collapse = ", "), call. = FALSE)
  }

  mapped_values <- names(labels)[label_index]
  mapped_values[is.na(value_keys)] <- NA_character_
  mapped_values
}

#' Return the configured label for a single clinical code
getLabelForValue <- function(value, labels, mapping_name) {
  mapValuesToLabels(value, labels, mapping_name)[[1L]]
}

#' Construct standardized sample identifiers
#'
#' Sample IDs have the form <subject>_<treatment label>_<timepoint label> and
#' therefore match the convention used by MOFA_analysis.R.
makeSampleIds <- function(subject_ids, treatment_codes, timepoint_value) {
  paste(
    subject_ids,
    mapValuesToLabels(treatment_codes, treatment_labels, "treatment"),
    getLabelForValue(timepoint_value, timepoint_labels, "timepoint"),
    sep = "_"
  )
}

#' Validate the observed values in a clinical-code column
validateClinicalCodes <- function(pheno_df, column_name, labels, mapping_name) {
  observed <- unique(as.character(pheno_df[[column_name]]))
  observed <- observed[!is.na(observed)]
  configured <- as.character(unname(labels))

  if (!setequal(observed, configured)) {
    stop(
      "Values in '", mapping_name, "' must match the unique non-missing values in clinical column '",
      column_name, "'. Observed: ", paste(sort(observed), collapse = ", "),
      "; configured: ", paste(sort(configured), collapse = ", "), ".",
      call. = FALSE
    )
  }
}

#' Validate shared study/clinical settings and derive timepoint context
#'
#' @param pheno_df Clinical data frame returned by readPhenoData().
#' @return A list containing baseline_timepoint and target_timepoint_label.
validateClinicalConfiguration <- function(pheno_df) {
  if (!is.character(study_name) || length(study_name) != 1L ||
      is.na(study_name) || !nzchar(trimws(study_name))) {
    stop("'study_name' must be a non-empty character value.", call. = FALSE)
  }

  clinical_col_settings <- list(
    subject_id_col = subject_id_col,
    treatment_col = treatment_col,
    timepoint_col = timepoint_col
  )
  if (any(!vapply(clinical_col_settings, function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }, logical(1)))) {
    stop("'subject_id_col', 'treatment_col', and 'timepoint_col' must each contain one non-empty column name.", call. = FALSE)
  }

  clinical_cols <- unname(unlist(clinical_col_settings, use.names = FALSE))
  if (anyDuplicated(clinical_cols)) {
    stop("'subject_id_col', 'treatment_col', and 'timepoint_col' must name three different clinical columns.", call. = FALSE)
  }

  missing_clinical_cols <- setdiff(clinical_cols, names(pheno_df))
  if (length(missing_clinical_cols) > 0L) {
    stop("Clinical data are missing configured column(s): ",
         paste(missing_clinical_cols, collapse = ", "), call. = FALSE)
  }

  validateLabelMap(timepoint_labels, "timepoint_labels", "baseline")
  validateLabelMap(treatment_labels, "treatment_labels")

  if (length(timepoint) != 1L || is.na(timepoint)) {
    stop("'timepoint' must contain one non-missing clinical timepoint code.", call. = FALSE)
  }

  validateClinicalCodes(pheno_df, timepoint_col, timepoint_labels, "timepoint_labels")
  validateClinicalCodes(pheno_df, treatment_col, treatment_labels, "treatment_labels")

  baseline <- unname(timepoint_labels[["baseline"]])
  target_label <- getLabelForValue(timepoint, timepoint_labels, "timepoint")

  if (use_change_scores && as.character(timepoint) == as.character(baseline)) {
    stop("'timepoint' must differ from the configured baseline when use_change_scores is TRUE.", call. = FALSE)
  }

  list(
    baseline_timepoint = baseline,
    target_timepoint_label = target_label
  )
}

#' Validate the README-note slots in omics_input_defs
validateOmicsInputDefinitions <- function(omics_input_defs) {
  if (!is.list(omics_input_defs) || is.null(names(omics_input_defs))) {
    stop("'omics_input_defs' must be a named list.", call. = FALSE)
  }

  for (view_name in names(omics_input_defs)) {
    view_def <- omics_input_defs[[view_name]]
    required_slots <- c("moniker", "file_var", "barcode_col", "reader", "readme_notes")
    missing_slots <- setdiff(required_slots, names(view_def))
    if (length(missing_slots) > 0L) {
      stop("Omics view '", view_name, "' is missing required slot(s): ",
           paste(missing_slots, collapse = ", "), call. = FALSE)
    }

    readme_notes <- view_def$readme_notes
    if (!is.list(readme_notes) ||
        !all(c("change_scores", "timepoint_values") %in% names(readme_notes)) ||
        any(!vapply(readme_notes[c("change_scores", "timepoint_values")],
                    function(note) is.character(note) && length(note) == 1L &&
                      !is.na(note) && nzchar(note), logical(1)))) {
      stop("Omics view '", view_name,
           "' must define non-empty README notes for 'change_scores' and 'timepoint_values'.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Read all configured omics views whose file paths are defined
#'
#' @param omics_input_defs Named view-definition list from a parameter file.
#' @param envir Environment containing the file-path variables named by file_var.
#' @return A list with entries loaded_views and barcode_cols.
loadConfiguredOmicsViews <- function(omics_input_defs, envir = parent.frame()) {
  loaded_views <- list()
  barcode_cols <- character(0)

  for (view_name in names(omics_input_defs)) {
    view_def <- omics_input_defs[[view_name]]
    input_file <- getResourceFile(view_def$file_var, envir = envir)
    if (is.null(input_file)) next

    message("Reading ", view_name, " data from: ", input_file)
    loaded_views[[view_name]] <- view_def$reader(input_file)
    barcode_cols[[view_name]] <- view_def$barcode_col
  }

  list(loaded_views = loaded_views, barcode_cols = barcode_cols)
}

#' Remove features with zero variance across selected samples
removeZeroVarianceFeatures <- function(data_mat) {
  if (ncol(data_mat) < 2L) return(data_mat)
  row_vars <- apply(data_mat, 1, stats::var, na.rm = TRUE)
  keep <- !is.na(row_vars) & row_vars != 0
  data_mat[keep, , drop = FALSE]
}

#' Compute target-minus-baseline change scores for one omics view
computeDelta <- function(data_mat, pheno_df, barcode_col, timepoint,
                         remove_zero_var = TRUE, rename_cols = TRUE) {
  tp_label <- getLabelForValue(timepoint, timepoint_labels, "timepoint")

  ph_slice <- pheno_df |>
    dplyr::filter(
      as.character(.data[[timepoint_col]]) %in%
        as.character(c(baseline_timepoint, timepoint)),
      .data[[barcode_col]] %in% colnames(data_mat)
    ) |>
    dplyr::mutate(
      .subject_id = .data[[subject_id_col]],
      .treatment_code = .data[[treatment_col]],
      .timepoint_label = mapValuesToLabels(
        .data[[timepoint_col]], timepoint_labels, "timepoint"
      )
    ) |>
    dplyr::group_by(.subject_id, .treatment_code, .timepoint_label) |>
    dplyr::summarise(barcode = dplyr::first(.data[[barcode_col]]), .groups = "drop") |>
    tidyr::pivot_wider(names_from = .timepoint_label, values_from = barcode)

  if (!("baseline" %in% names(ph_slice))) ph_slice[["baseline"]] <- NA_character_
  if (!(tp_label %in% names(ph_slice))) ph_slice[[tp_label]] <- NA_character_

  ph_slice <- ph_slice |> dplyr::filter(!is.na(.data[["baseline"]]))
  ph_slice_tp <- ph_slice |> dplyr::filter(!is.na(.data[[tp_label]]))

  if (nrow(ph_slice_tp) == 0L) {
    return(list(delta = data_mat[, FALSE, drop = FALSE], slice = ph_slice_tp))
  }

  delta <- data_mat[, ph_slice_tp[[tp_label]], drop = FALSE] -
    data_mat[, ph_slice_tp[["baseline"]], drop = FALSE]

  if (remove_zero_var) delta <- removeZeroVarianceFeatures(delta)

  if (rename_cols) {
    name_map <- setNames(
      makeSampleIds(ph_slice_tp$.subject_id, ph_slice_tp$.treatment_code, timepoint),
      ph_slice_tp[[tp_label]]
    )
    colnames(delta) <- name_map[colnames(delta)]
  }

  list(delta = delta, slice = ph_slice_tp)
}

#' Extract one target timepoint from an omics view without baseline subtraction
extractTimepointData <- function(data_mat, pheno_df, barcode_col, timepoint,
                                 remove_zero_var = TRUE, rename_cols = TRUE) {
  ph_slice_tp <- pheno_df |>
    dplyr::filter(
      as.character(.data[[timepoint_col]]) == as.character(timepoint),
      .data[[barcode_col]] %in% colnames(data_mat)
    ) |>
    dplyr::mutate(
      .subject_id = .data[[subject_id_col]],
      .treatment_code = .data[[treatment_col]]
    ) |>
    dplyr::group_by(.subject_id, .treatment_code) |>
    dplyr::summarise(barcode = dplyr::first(.data[[barcode_col]]), .groups = "drop")

  if (nrow(ph_slice_tp) == 0L) {
    return(list(data = data_mat[, FALSE, drop = FALSE], slice = ph_slice_tp))
  }

  data_tp <- data_mat[, ph_slice_tp$barcode, drop = FALSE]
  if (remove_zero_var) data_tp <- removeZeroVarianceFeatures(data_tp)

  if (rename_cols) {
    name_map <- setNames(
      makeSampleIds(ph_slice_tp$.subject_id, ph_slice_tp$.treatment_code, timepoint),
      ph_slice_tp$barcode
    )
    colnames(data_tp) <- name_map[colnames(data_tp)]
  }

  list(data = data_tp, slice = ph_slice_tp)
}

#' Prepare one omics view according to the configured score mode
prepareOmicsView <- function(data_mat, pheno_df, barcode_col, timepoint,
                             use_change_scores) {
  if (use_change_scores) {
    return(computeDelta(data_mat, pheno_df, barcode_col, timepoint)$delta)
  }
  extractTimepointData(data_mat, pheno_df, barcode_col, timepoint)$data
}

#' Prepare all loaded omics views and drop views with no target samples
prepareConfiguredOmicsViews <- function(loaded_views, pheno_df, barcode_cols,
                                        timepoint, use_change_scores) {
  data_views <- list()
  for (view_name in names(loaded_views)) {
    data_views[[view_name]] <- prepareOmicsView(
      loaded_views[[view_name]], pheno_df, barcode_cols[[view_name]],
      timepoint, use_change_scores
    )
  }

  if (length(data_views) == 0L) return(data_views)

  sample_counts <- vapply(data_views, ncol, integer(1))
  if (any(sample_counts == 0L)) {
    warning(
      "Dropping omics views with no samples at the requested timepoint: ",
      paste(names(data_views)[sample_counts == 0L], collapse = ", ")
    )
    data_views <- data_views[sample_counts > 0L]
  }
  data_views
}

#' Build target-visit sample metadata using the shared sample-ID convention
buildTargetVisitSampleMetadata <- function(pheno_df, timepoint_value = timepoint) {
  out <- pheno_df |>
    dplyr::filter(as.character(.data[[timepoint_col]]) == as.character(timepoint_value)) |>
    dplyr::mutate(
      sample_id = makeSampleIds(
        .data[[subject_id_col]], .data[[treatment_col]], timepoint_value
      )
    ) |>
    dplyr::distinct(sample_id, .keep_all = TRUE)

  rownames(out) <- out$sample_id
  out
}

#' Scale each feature across samples using z-scores
#'
#' This helper is used by the GLM external-scoring workflow when scale_views is
#' TRUE. Constant features are removed before this step by the common
#' preprocessing functions, so finite standard deviations are normally expected.
scaleOmicsViewFeatures <- function(data_mat) {
  if (ncol(data_mat) < 2L) return(data_mat)
  scaled <- t(scale(t(data_mat)))
  if (any(!is.finite(scaled))) {
    stop("View scaling produced non-finite values. Check for zero-variance or missing features.", call. = FALSE)
  }
  scaled
}
