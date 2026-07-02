#' Extract MOFA2 Factor Loadings to an Excel Workbook
#'
#' Extracts factor-loading/weight values for all views from a fitted MOFA2 object 
#' and saves them into an Excel workbook (`mofa_factor_loadings.xlsx`) with one 
#' worksheet per data view. The data is exported in a clean, long format.
#'
#' @param mofa_obj Either a character string giving the full path to a fitted 
#'   `MOFA2` model file (e.g., an `.hdf5` file), or an already-loaded fitted 
#'   `MOFA2` object.
#' @param output_dir A character string specifying the directory where the Excel 
#'   file should be saved. Defaults to \code{NULL}. If \code{NULL}:
#'   \itemize{
#'     \item If \code{mofa_obj} is a file path, saves in the same directory as the file.
#'     \item If \code{mofa_obj} is an object, saves in the current working directory.
#'   }
#'
#' @return Invisibly returns the full normalized path to the generated Excel file.
#' @export
#'
#' @importFrom MOFA2 load_model get_weights views factors
#' @importFrom openxlsx createWorkbook addWorksheet writeData saveWorkbook
#' @importFrom tools file_path_as_absolute
extractMofaFactors <- function(mofa_obj, output_dir = NULL) {
  
  # ---------------------------------------------------------------------------
  # 1. Input Validation & MOFA Object Loading
  # ---------------------------------------------------------------------------
  
  is_path <- is.character(mofa_obj) && length(mofa_obj) == 1
  
  if (is_path) {
    # Ensure the file actually exists before trying to load it
    if (!file.exists(mofa_obj)) {
      stop("The provided mofa_obj file path does not exist: ", mofa_obj)
    }
    
    # Load the MOFA2 object from the HDF5 file
    message("Loading MOFA2 model from file...")
    # Wrap in tryCatch to handle corrupted or non-MOFA HDF5 files
    mofa_model <- tryCatch({
      MOFA2::load_model(mofa_obj)
    }, error = function(e) {
      stop("Failed to load MOFA2 model from file. Ensure it is a valid MOFA2 HDF5 file. Error: ", e$message)
    })
  } else if (inherits(mofa_obj, "MOFA")) {
    mofa_model <- mofa_obj
  } else {
    stop("Argument 'mofa_obj' must be either a character string path to a MOFA2 file or a fitted MOFA object.")
  }
  
  # Verify that the object contains views and factors
  model_views <- MOFA2::views_names(mofa_model)
  model_factors <- MOFA2::factors_names(mofa_model)
  
  if (length(model_views) == 0) {
    stop("The provided MOFA object does not contain any data views.")
  }
  if (length(model_factors) == 0 || (is.list(model_factors) && length(model_factors[[1]]) == 0)) {
    stop("The provided MOFA object does not contain any factors.")
  }
  
  # ---------------------------------------------------------------------------
  # 2. Determine and Prepare Output Directory
  # ---------------------------------------------------------------------------
  
  target_dir <- output_dir
  
  if (is.null(target_dir)) {
    if (is_path) {
      # Use the same directory as the input file
      target_dir <- dirname(mofa_obj)
    } else {
      # Use the current R working directory
      target_dir <- getwd()
    }
  }
  
  # Create directory if it doesn't exist
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
    if (!dir.exists(target_dir)) {
      stop("Failed to create the output directory: ", target_dir)
    }
  }
  
  # Define full path for the output file
  output_file_path <- file.path(target_dir, "mofa_factor_loadings.xlsx")
  
  # ---------------------------------------------------------------------------
  # 3. Extract and Transform Factor Loadings
  # ---------------------------------------------------------------------------
  
  # Extract weights/loadings using official MOFA2 accessor
  # Returns a list of matrices (one per view), dimensions: features x factors
  weights_list <- MOFA2::get_weights(mofa_model, views = "all", factors = "all")
  
  # Initialize openxlsx Workbook
  wb <- openxlsx::createWorkbook()
  
  for (view_name in model_views) {
    # Extract matrix for the specific view
    weight_matrix <- weights_list[[view_name]]
    
    if (is.null(weight_matrix) || nrow(weight_matrix) == 0) {
      warning(paste("No loading data found for view:", view_name, "- Skipping worksheet creation."))
      next
    }
    
    # Standardize Excel Worksheet Names:
    # 1. Excel sheets cannot contain: \ , / , ? , * , : , [ , ]
    # 2. Max length is 31 characters
    clean_sheet_name <- gsub("[\\\\/?*:\\[\\]]", "_", view_name)
    if (nchar(clean_sheet_name) > 31) {
      clean_sheet_name <- substr(clean_sheet_name, 1, 31)
    }
    
    # Convert matrix to a tidy long-format data frame
    # Row names are Features, Column names are Factors (e.g., "Factor1")
    features <- rownames(weight_matrix)
    factors_vec <- colnames(weight_matrix)
    
    # Format factor names strictly as "Factor 1", "Factor 2", etc.
    formatted_factors <- gsub("Factor\\s*([0-9]+)", "Factor \\1", factors_vec)
    
    # Construct the long-format data frame using standard matrix reshaping
    # This avoids a strict dependency on tidyr/melt
    long_list <- lapply(seq_along(factors_vec), function(i) {
      data.frame(
        Factor = formatted_factors[i],
        Feature = features,
        Loading = as.numeric(weight_matrix[, i]),
        stringsAsFactors = FALSE
      )
    })
    
    view_long_df <- do.call(rbind, long_list)
    
    # Ensure correct column ordering and structure
    view_long_df <- view_long_df[, c("Factor", "Feature", "Loading")]
    
    # ---------------------------------------------------------------------------
    # 4. Write to Excel Workbook
    # ---------------------------------------------------------------------------
    
    # Add worksheet and write data
    openxlsx::addWorksheet(wb, sheetName = clean_sheet_name)
    openxlsx::writeData(wb, sheet = clean_sheet_name, x = view_long_df)
  }
  
  # Save workbook to file
  openxlsx::saveWorkbook(wb, file = output_file_path, overwrite = TRUE)
  
  message("Successfully saved MOFA factor loadings to: ", output_file_path)
  
  # Return the absolute path invisibly
  invisible(tools::file_path_as_absolute(output_file_path))
}