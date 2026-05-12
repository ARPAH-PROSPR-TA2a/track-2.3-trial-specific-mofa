# *****************************************************************************
# Open files, prepare data, and run MOFA analysis. Three variables control
# how the analysis is configured:
#   timepoint, use_change_scores, use_groups
#
# The variable timepoint specifies the follow-up visit to analyze
# (identified by fumo == timepoint). This determines which samples are
# used as the target condition.
#
# If use_change_scores is TRUE, omics values are converted to change
# scores by subtracting baseline values (fumo == 0) from the target
# visit (fumo == timepoint).
#
# If use_change_scores is FALSE, omics values from the target visit
# (fumo == timepoint) are used directly, without baseline subtraction.
#
# If use_groups is TRUE, treatment groups (caloric restriction or normal
# diet) are modeled as separate MOFA groups. If FALSE, all samples are
# analyzed together without grouping.
#
# The script supports any combination of input data matrices
# (methylation, proteomics, metabolomics). The file paths for these are
# defined in resources.R. To exclude a data type, simply comment out the
# corresponding file path in resources.R; the script will automatically
# use only the available data matrices.
#
# To ensure that savehistory() captures all executed commands, run this code
# by pasting it into a fresh R GUI session.
# *****************************************************************************

Sys.setenv(R_HISTSIZE = 2000) # to make sure all history is saved at the end.

library(MOFA2)
library(dplyr)
library(tidyr)

########################################
#------------ Configure run by setting the values of the variables below
# Location of source code root directory - modify as needed for your environment
code_dir <- "D:/git/misc/prospr-fast/CALERIE"

# Location of directory where to store results and output - modify as desired
# for your local run environment
output_dir <- "D:/R_sessions/PROSPR/MOFA"

# Specify the target visit timepoint, either 12 or 24
timepoint <- 24

# If TRUE, compute (timepoint - baseline) change scores. If FALSE, use the
# source data directly for samples where fumo == timepoint.
use_change_scores <- TRUE

# If TRUE, run MOFA with AL and CR as separate groups. If FALSE, do not pass
# groups to MOFA, so all samples are analyzed together.
use_groups <- FALSE
#------------------------------------------------------------------------------
########################################

# Source files. 
source(file.path(code_dir, "resources.R"))
source(file.path(code_dir, "methods.R"))
# Names of variables specifying file locations, expected to be defined in 
# resources.R
resource_vars <- c("calerie_methy_fn", "calerie_prot_fn",
		"calerie_metab_data_fn", "calerie_metab_annot_fn",
		"calerie_clin_fn")

# Return TRUE only if the variable is defined in resources.R and contains a
# non-empty file name.
resourceFileSpecified <- function(var_name){
	if (!exists(var_name)){
		return(FALSE)
	}	
	value <- get(var_name)
	is.character(value) && length(value) == 1 && !is.na(value) && nzchar(value)
}

# Get a file name defined in resources.R. Required files trigger an error if
# missing; optional files return NULL if missing.
getResourceFile <- function(var_name, required = FALSE){
	if (!resourceFileSpecified(var_name)){
		if (required){
			stop(paste0("Required variable '", var_name,
							"' is not defined in resources.R."))
		}
		return(NULL)
	}
	get(var_name)
}

# Read phenotypic data. This file is required.
pheno_d <- readPhenoData(getResourceFile("calerie_clin_fn", required = TRUE))

# Definitions for the optional omics data views.
omics_input_defs <- list(
		methylation = list(
				moniker = "me",
				file_var = "calerie_methy_fn",
				barcode_col = "blood_DNAm_barcode",
				reader = function(fn) readMethylation(fn)
		),
		proteomics = list(
				moniker = "pr",
				file_var = "calerie_prot_fn",
				barcode_col = "proteomics_barcode",
				reader = function(fn) readProteomics(fn)
		),
		metabolomics = list(
				moniker = "mb",
				file_var = "calerie_metab_data_fn",
				barcode_col = "metabolomics_barcode",
				reader = function(fn) readMetabolomics(fn)[["data"]]
		)
)

# Read only the omics data views whose file names are defined in resources.R.
loaded_views <- list()
barcode_cols <- character(0)
for (view_name in names(omics_input_defs)){
	view_def <- omics_input_defs[[view_name]]
	input_file <- getResourceFile(view_def$file_var)
	
	if (is.null(input_file)){
		next
	}
	
	message("Reading ", view_name, " data from: ", input_file)
	loaded_views[[view_name]] <- view_def$reader(input_file)
	barcode_cols[view_name] <- view_def$barcode_col
}

if (length(loaded_views) == 0){
	stop("No omics data file names are defined in resources.R.")
}

# Remove features with zero variance across selected samples. If fewer than two
# samples are present, variance cannot be estimated, so the data are returned
# unchanged.
removeZeroVarianceFeatures <- function(data_mat){
	if (ncol(data_mat) < 2){
		return(data_mat)
	}
	
	row_vars <- apply(data_mat, 1, var, na.rm = TRUE)
	keep <- !is.na(row_vars) & row_vars != 0
	data_mat[keep, , drop = FALSE]
}

# Compute delta (timepoint - baseline) for a single data view
#
# Generate change scores for the target visit against baseline. The helper
# extracts samples with both a baseline measurement and a measurement at the
# requested timepoint, then computes per-feature differences.
#
# Columns can be renamed to the format <DEID>_<AL|CR>_<visit>.
computeDelta <- function(data_mat, pheno_df, barcode_col, timepoint,
		remove_zero_var = TRUE, rename_cols = TRUE){
	tp_label <- paste0("mo_", timepoint)
	
	# Build the pheno slice: keep only rows for visits 0, timepoint and where 
	# the barcode for this view is present in the data matrix's columns.
	ph_slice <- pheno_df %>%
			filter(
					as.character(fumo) %in% c("0", as.character(timepoint)),
					.data[[barcode_col]] %in% colnames(data_mat)
			) %>%
			# normalise visit codes to strings like 'mo_0', 'mo_12' etc.
			mutate(fumo = paste0("mo_", fumo)) %>%
			# collapse duplicates by DEID/CR/fumo keeping the first barcode
			group_by(DEID, CR, fumo) %>%
			summarise(barcode = first(.data[[barcode_col]]), .groups = "drop") %>%
			# pivot to wide so each visit becomes its own column (mo_0, mo_12, mo_24)
			pivot_wider(names_from = fumo, values_from = barcode)
	
	if (!("mo_0" %in% names(ph_slice))){
		ph_slice[["mo_0"]] <- NA_character_
	}
	if (!(tp_label %in% names(ph_slice))){
		ph_slice[[tp_label]] <- NA_character_
	}
	
	# Require baseline and timepoint measurements present
	ph_slice <- ph_slice %>% filter(!is.na(mo_0))
	ph_slice_tp <- ph_slice %>% filter(!is.na(.data[[tp_label]]))
	
	if (nrow(ph_slice_tp) == 0){
		return(list(delta = data_mat[, FALSE, drop = FALSE],
						slice = ph_slice_tp))
	}
	
	# Compute delta (timepoint - baseline)
	delta <- data_mat[, ph_slice_tp[[tp_label]], drop = FALSE] -
			data_mat[, ph_slice_tp[["mo_0"]], drop = FALSE]
	
	# Eemove features with zero variance across selected delta
	# samples (primarily useful for proteomics/metabolomics views)
	if (remove_zero_var){
		delta <- removeZeroVarianceFeatures(delta)
	}
	
	# Rename columns to a standardized format: <DEID>_<AL|CR>_<visit>
	if (rename_cols){
		name_map <- setNames(
				paste(ph_slice_tp$DEID,
						ifelse(ph_slice_tp$CR == 0, "AL", "CR"),
						paste0(timepoint, "m"), sep = "_"),
				ph_slice_tp[[tp_label]]
		)
		colnames(delta) <- name_map[colnames(delta)]
	}
	
	list(delta = delta, slice = ph_slice_tp)
}

# Extract target-visit data for a single data view without subtracting baseline.
extractTimepointData <- function(data_mat, pheno_df, barcode_col, timepoint,
		remove_zero_var = TRUE, rename_cols = TRUE){
	ph_slice_tp <- pheno_df %>%
			filter(
					as.character(fumo) == as.character(timepoint),
					.data[[barcode_col]] %in% colnames(data_mat)
			) %>%
			group_by(DEID, CR) %>%
			summarise(barcode = first(.data[[barcode_col]]), .groups = "drop")
	
	if (nrow(ph_slice_tp) == 0){
		return(list(data = data_mat[, FALSE, drop = FALSE],
						slice = ph_slice_tp))
	}
	
	data_tp <- data_mat[, ph_slice_tp$barcode, drop = FALSE]
	
	if (remove_zero_var){
		data_tp <- removeZeroVarianceFeatures(data_tp)
	}
	
	if (rename_cols){
		name_map <- setNames(
				paste(ph_slice_tp$DEID,
						ifelse(ph_slice_tp$CR == 0, "AL", "CR"),
						paste0(timepoint, "m"), sep = "_"),
				ph_slice_tp$barcode
		)
		colnames(data_tp) <- name_map[colnames(data_tp)]
	}
	
	list(data = data_tp, slice = ph_slice_tp)
}

# Prepare one omics view according to the use_change_scores setting.
prepareOmicsView <- function(data_mat, pheno_df, barcode_col, timepoint,
		use_change_scores){
	if (use_change_scores){
		return(computeDelta(data_mat, pheno_df, barcode_col, timepoint)$delta)
	}
	
	extractTimepointData(data_mat, pheno_df, barcode_col, timepoint)$data
}

# Apply the selected preprocessing to all available omics views.
data_views <- list()
for (view_name in names(loaded_views)){
	data_views[[view_name]] <- prepareOmicsView(
			loaded_views[[view_name]], pheno_d, barcode_cols[[view_name]],
			timepoint, use_change_scores
	)
}

# Drop any view that has no samples at the requested timepoint.
sample_counts <- vapply(data_views, ncol, integer(1))
if (any(sample_counts == 0)){
	warning(
			"Dropping omics views with no samples at the requested timepoint: ",
			paste(names(data_views)[sample_counts == 0], collapse = ", ")
	)
	data_views <- data_views[sample_counts > 0]
}

if (length(data_views) == 0){
	stop("No omics views have samples at the requested timepoint.")
}

# Extract the AL/CR label from sample IDs of the form <DEID>_<AL|CR>_<visit>.
# The treatment label is taken from the penultimate field, so this remains
# valid if DEID itself contains underscores.
getSampleGroups <- function(sample_ids){
	groups <- vapply(strsplit(sample_ids, split = "_", fixed = TRUE),
			function(parts){
				if (length(parts) < 3){
					return(NA_character_)
				}
				parts[length(parts) - 1]
			}, character(1))
	
	if (any(is.na(groups) | !nzchar(groups))){
		stop("Could not infer AL/CR labels from all sample IDs.")
	}
	
	names(groups) <- sample_ids
	groups
}

# Return the indices of the most variable features in a matrix. If there are
# fewer than two samples, feature variance cannot be estimated, so the first
# max_features rows are returned.
getTopVariableFeatureIndex <- function(view_data, max_features){
	max_features <- min(max_features, nrow(view_data))
	
	if (max_features >= nrow(view_data)){
		return(seq_len(nrow(view_data)))
	}
	if (ncol(view_data) < 2){
		return(seq_len(max_features))
	}
	
	row_vars <- apply(view_data, 1, var, na.rm = TRUE)
	row_vars[is.na(row_vars)] <- -Inf
	order(row_vars, decreasing = TRUE)[seq_len(max_features)]
}

# Select highly variable features within each group, then take their union.
# If the union exceeds max_features, retain the features with the largest
# within-group variance across all groups. This is how MOFA recommends
# selection of highly variable features be handled when using groups.
getTopVariableFeatureIndexByGroup <- function(view_data, view_groups,
		max_features){
	# Ensure that the group vector matches the number of samples (columns)
	if (length(view_groups) != ncol(view_data)){
		stop("The group vector must contain one label per sample column.")
	}
	
	# Disallow missing or empty group labels
	if (any(is.na(view_groups) | !nzchar(view_groups))){
		stop("Missing group labels are not allowed for grouped HVF selection.")
	}
	
	# Cap max_features to the total number of available features (rows)
	max_features <- min(max_features, nrow(view_data))
	
	# Identify unique group labels (e.g., "AL", "CR")
	group_names <- unique(view_groups)
	
	# For each group:
	#   - subset the data to samples in that group
	#   - compute the top variable features within that subset
	# Combine all selected feature indices across groups
	group_feat_idx <- unlist(lapply(group_names, function(group_name){
				group_cols <- view_groups == group_name
				getTopVariableFeatureIndex(view_data[, group_cols, drop = FALSE],
						max_features)
			}), use.names = FALSE)
	# Take the union of features selected across all groups
	feat_idx <- unique(group_feat_idx)
	
	# If the total number of selected features does not exceed the limit,
	# return them directly
	if (length(feat_idx) <= max_features){
		return(feat_idx)
	}
	
	# For each group, compute the variance of each feature (row-wise)
	group_vars <- lapply(group_names, function(group_name){
				group_data <- view_data[, view_groups == group_name, drop = FALSE]
				if (ncol(group_data) < 2){
					return(rep(NA_real_, nrow(group_data)))
				}
				apply(group_data, 1, var, na.rm = TRUE)
			})
	
	# Combine per-group variance vectors into a matrix:
	# rows = features, columns = groups
	group_var_mat <- do.call(cbind, group_vars)
	group_var_mat[is.na(group_var_mat)] <- -Inf
	
	# For each feature (restricted to those already selected across groups),
	# take the maximum variance observed across all groups
	feature_scores <- apply(group_var_mat[feat_idx, , drop = FALSE], 1, max)
	# Rank features by their maximum group-specific variance and
	# select the top max_features
	feat_idx[order(feature_scores, decreasing = TRUE)[seq_len(max_features)]]
}

# MOFA recommends that data views have a comparable number of features, to avoid
# one view dominating factor fitting. We therefore reduce large views by keeping
# only the most variable features. A view may have at most max_factor times as
# many features as the smallest view, and never more than max_features.
#
# If use_groups is TRUE, highly variable features are selected separately
# within each group and then combined across groups before enforcing the final
# feature cap, per MOFA guidelines.
filterMostVariableFeatures <- function(data_views, max_factor=5, groups = NULL,
		use_groups = FALSE, max_features = 10000){
	# Get number of features (rows) for each omics view
	feature_counts <- vapply(data_views, nrow, integer(1))
	
	if (any(feature_counts == 0)){
		stop("At least one omics view has zero features after preprocessing.")
	}
	if (use_groups && is.null(groups)){
		stop("groups must be provided when use_groups is TRUE.")
	}
	if (use_groups && is.null(names(groups))){
		stop("groups must be a named vector when use_groups is TRUE.")
	}
	
	min_feat_num <- min(feature_counts)
	
	# Iterate over each omics view (data matrix)
	mapply(function(view_data, view_name){
				
				# Compute number of features to retain for this view:
				# limited by:
				#   - max_factor * min features across views
				#   - global max_features cap
				#   - total features available in this view
				vf_num <- min(max_factor * min_feat_num, max_features,
						nrow(view_data))
				vf_num <- floor(vf_num)
				
				if (vf_num >= nrow(view_data)){
					return(view_data)
				}
				
				if (use_groups){
					# Map sample columns to their group labels
					view_groups <- groups[colnames(view_data)]
					# Ensure all samples have valid group assignments
					if (any(is.na(view_groups))){
						stop(paste0("Missing group labels for view '", view_name,
								"'."))
					}
					# Select top variable features using group-aware method
					feat_idx <- getTopVariableFeatureIndexByGroup(view_data,
							view_groups, vf_num)
				} else {
					# Select top variable features across all samples (no grouping)
					feat_idx <- getTopVariableFeatureIndex(view_data, vf_num)
				}
				
				view_data[feat_idx, , drop = FALSE]
			}, data_views, names(data_views), SIMPLIFY = FALSE)
}

# Set scaling factor controlling how many features to retain per view
# (relative to the smallest view) 
max_factor <- 5
# Set an absolute upper bound on the number of features per view - useful
# primarily when running MOFA with only the methylation data view, to
# avoid performance degradation due to hundreds of thousands of features
max_features <- 10000
# Collect the union of all sample IDs across views (column names),
# ensuring that group assignment covers all samples present
sample_ids_for_filter <- Reduce(union, lapply(data_views, colnames))
sample_groups_for_filter <- if (use_groups){
			getSampleGroups(sample_ids_for_filter)
		} else {
			NULL
		}

# Apply feature filtering to each view:
#   - selects most variable features
#   - uses group-aware selection if use_groups == TRUE
#   - enforces both relative (max_factor) and absolute (max_features) limits
data_views <- filterMostVariableFeatures(data_views, max_factor,
		groups = sample_groups_for_filter, use_groups = use_groups,
		max_features = max_features)

# Indicates if the data should be scaled. According to the MOFA paper, scaling
# is not needed. Scaling the data, however, may be helpful for some plots; for
# example, clusters may be more distinguishable in heatmaps when using z-scores.
scale_views <- FALSE

# Prep data matrices by filling in missing columns. This is a requirement of
# the MOFA package: if a data view is missing data for some samples, these need
# to be added with NA values.
sample_ids <- Reduce(union, lapply(data_views, colnames))
if (length(sample_ids) == 0){
	stop("No samples are available after preprocessing the omics views.")
}

# Prep data matrices by filling in missing columns. This is a requirement of
# the MOFA package: if a data view is missing data for some samples, these
# need to be added with NA values.
data_views <- lapply(data_views, function(view_data){
			view_data <- fillEmptyCols(view_data, sample_ids)
			view_data[, sample_ids, drop = FALSE]
		})

# Assign each sample to its AL/CR treatment label. These labels are passed to
# MOFA only when use_groups is TRUE.
grps <- if (use_groups){
			getSampleGroups(sample_ids)
		} else {
			NULL
		}

# Create MOFA object, comprising all available data matrices (aka, views).
if (use_groups){
	mdd <- create_mofa(data_views, groups = grps)
} else {
	mdd <- create_mofa(data_views)
}
views_names(mdd) <- names(data_views)

# Set data options.
data_opts <- get_default_data_options(mdd)
if (scale_views){
	data_opts$center_groups <- TRUE
	data_opts$scale_groups <- TRUE
}
data_opts$use_float32 <- FALSE

# Set model options.
model_opts <- get_default_model_options(mdd)
# model_opts$num_factors <- 10 # If we want to manually set the factor num

# Set training options.
train_opts <- get_default_training_options(mdd)
train_opts$drop_factor_threshold <- 0.01
train_opts$convergence_mode <- "slow"
# increase from default value of 1000, just in case more iters are required 
# for convergence
train_opts$maxiter <- 4000

# Prepare the MOFA object and add the metadata data frame, with one row per
# sample. The sample column is always required; the group column is added only
# when running MOFA in multi-group mode.
mo <- prepare_mofa(object = mdd, data_options = data_opts,
		model_options = model_opts, training_options = train_opts)

meta_d <- pheno_d %>%
		filter(as.character(fumo) == as.character(timepoint)) %>%
		mutate(sample_id = paste(DEID, ifelse(CR == 0, "AL", "CR"),
						paste0(timepoint, "m"), sep = "_")) %>%
		distinct(sample_id, .keep_all = TRUE)
rownames(meta_d) <- meta_d$sample_id

missing_meta <- setdiff(sample_ids, rownames(meta_d))
if (length(missing_meta) > 0){
	stop("Missing phenotypic metadata for samples: ",
			paste(missing_meta, collapse = ", "))
}

meta_d <- meta_d[sample_ids, , drop = FALSE]
meta_d <- data.frame(meta_d, sample = rownames(meta_d))
if (use_groups){
	meta_d$group <- grps[rownames(meta_d)]
}
samples_metadata(mo) <- meta_d

# Train the MOFA model.
#
# Note the parameter setting "use_basilisk = TRUE". This forces a local Python
# environment to be installed in order to download and run the mofapy2 package,
# which performs the MOFA computations.
#
# This setup occurs only once, the first time the code is run on a given
# machine. The initial run may take significantly longer, as it involves
# downloading and configuring several gigabytes of dependencies.
#
# Alternatively, if a suitable Python environment already exists, you can
# install mofapy2 manually and run the model with "use_basilisk = FALSE".

# Tag the output folder with the current date and time, for logging purposes.
# Also save the output of the MOFA method.
score_tag <- ifelse(use_change_scores, "change_scores", "timepoint_values")
group_tag <- ifelse(use_groups, "grouped", "ungrouped")
views_tag <- paste(sapply(views_names(mo), function(vn){
					omics_input_defs[[vn]]$moniker}), 
		collapse="+")
res_folder <- paste(output_dir,
		paste0(format(Sys.time(), "%Y_%b_%d_%H_%M_%S"),
				"_visit=", timepoint, "mo_", score_tag, "_", group_tag,
				"_views=", views_tag),
		sep = "/")
res_file <- paste(res_folder, "model.hdf5", sep = "/")
log_file <- paste(res_folder, "mofa_out.txt", sep = "/")
dir.create(res_folder, recursive = TRUE)
sink(file = log_file)
mo.tr <- run_mofa(mo, res_file, use_basilisk = TRUE)
sink()

# Logging - write out the history of commands executed as well as the session
# info.
session_file <- paste(res_folder, "sessionInfo.txt", sep = "/")
writeLines(capture.output(sessionInfo()), session_file)
hist_file <- paste(res_folder, "MOFA_analysis_calerie.R", sep = "/")
savehistory(hist_file)

# Document what was done.
readme_file <- paste(res_folder, "README.txt", sep = "/")
score_desc <- if (use_change_scores){
			paste0("change scores (", timepoint, "m - baseline)")
		} else {
			paste0("source values from samples where fumo == ", timepoint)
		}

included_views <- paste(names(data_views), collapse = ", ")
data_notes <- character(0)
if ("methylation" %in% names(data_views)){
	data_notes <- c(data_notes, "Methylation data are reported as M-values.")
}
if ("metabolomics" %in% names(data_views)){
	metab_note <- if (use_change_scores){
				paste("Metabolomics data in CALERIE are reported in Metabolon",
						"median-centered units; accordingly, we log2 transform them",
						"prior to change-score calculation.")
			} else {
				paste("Metabolomics data in CALERIE are reported in Metabolon",
						"median-centered units; accordingly, we log2 transform them",
						"prior to analysis.")
			}
	data_notes <- c(data_notes, metab_note)
}
if ("proteomics" %in% names(data_views)){
	data_notes <- c(data_notes,
			"Proteomics data, from Somascan, are already log2 transformed.")
}

group_desc <- if (use_groups){
			paste("Samples were split into two MOFA groups: subjects receiving",
					"the caloric restriction diet (CR group) and subjects on the",
					"alternative diet (AL group).")
		} else {
			paste("Samples were not split into MOFA groups; all samples were",
					"analyzed together. AL/CR labels remain available in the sample",
					"IDs and phenotype metadata.")
		}

readme_text <- paste(c(
				paste("MOFA was run on", score_desc, "."),
				paste("Included omics views:", included_views, "."),
				data_notes,
				"Each included omics data table becomes a separate MOFA view.",
				group_desc
		), collapse = " ")

sink(file = readme_file)
writeLines(strwrap(readme_text, width = 80))
sink()
