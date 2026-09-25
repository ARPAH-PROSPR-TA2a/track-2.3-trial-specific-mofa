# MOFA integrative factor interpretation pipeline
#
# Purpose:
#   Given a MOFA/MOFA2 results object, extract per-factor feature loadings from
#   each available omics view, run view-specific enrichment analyses, and render
#   a factor-by-factor PDF report with top features, enriched terms, and a
#   cautious biological interpretation.
#
# Supported MOFA view names:
#   proteomics     features are UniProt IDs
#   methylation    features are Illumina CpG IDs
#   metabolomics   features are Metabolon metabolite names/descriptions
#   expression     features are Entrez IDs or gene symbols
#
# Core design:
#   1. Extract MOFA weights/loadings.
#   2. Select high-loading features for outputs and for metabolomics ORA.
#   3. Run view-appropriate enrichment:
#        - proteomics: UniProt -> Entrez -> one signed ranked GO fgsea run per factor/view
#        - expression: Entrez or symbol -> Entrez -> one signed ranked GO fgsea run per factor/view
#        - methylation: CpG -> Entrez mapping via missMethyl, followed by one signed ranked GO fgsea run per factor/view
#        - metabolomics: positive/negative Metabolon SUPER_PATHWAY / SUB_PATHWAY / classification ORA
#   4. Render a report.
#
# Notes:
#   - The script is intentionally conservative. It produces plausible biological
#     summaries, but it does not claim causality.
#   - You should review the generated interpretations manually.
#   - For methylation, set methylation_array_type to the appropriate value for
#     your array, usually "EPIC".
#   - The default base_pdf report engine does not require Pandoc or LaTeX.
#   - Enrichment tasks can be parallelized with n_cores > 1.
#   - The optional rmarkdown_pdf report engine requires Pandoc plus a LaTeX
#     installation. The lightweight TinyTeX distribution is usually the easiest
#     option from R:
#       install.packages("tinytex")
#       tinytex::install_tinytex()

suppressPackageStartupMessages({
			library(dplyr)
			library(tidyr)
			library(purrr)
			library(stringr)
			library(readr)
			library(tibble)
			library(rlang)
		})

# -----------------------------------------------------------------------------
# Optional package helpers
# -----------------------------------------------------------------------------

#' Require an R package namespace for a specific analysis step.
#'
#' This helper checks whether a package is installed and stops with a clear
#' message if the package is missing. It is used to make optional dependencies
#' fail only when the corresponding analysis step is requested.
#'
#' @param pkg Character string naming the package namespace to check.
#' @param purpose Optional character string describing why the package is needed.
#'
#' @return Invisibly returns `TRUE` when the package is available; otherwise stops.
require_namespace <- function(pkg, purpose = NULL) {
	# Check optional dependencies at the point where they are needed.
	if (!requireNamespace(pkg, quietly = TRUE)) {
		msg <- paste0("Package '", pkg, "' is required")
		if (!is.null(purpose)) msg <- paste0(msg, " for ", purpose)
		msg <- paste0(msg, ". Please install it before running this step.")
		stop(msg, call. = FALSE)
	}
	invisible(TRUE)
}

#' Test whether an R package namespace is available.
#'
#' This lightweight helper is used for optional analysis branches where a missing
#' package should skip the branch rather than stop the full pipeline.
#'
#' @param pkg Character string naming the package namespace to check.
#'
#' @return Logical value indicating whether the package can be loaded.
is_pkg_available <- function(pkg) {
	# Return a simple TRUE/FALSE value without stopping the pipeline.
	requireNamespace(pkg, quietly = TRUE)
}

#' Normalize a requested CPU count.
#'
#' Converts a user-supplied CPU count to a valid positive integer and caps it at
#' the number of detected logical CPUs when that information is available.
#'
#' @param n_cores Requested number of worker processes.
#'
#' @return Integer number of cores to use.
normalize_n_cores <- function(n_cores = 1L) {
	# Convert invalid or missing requests to a conservative single-core default.
	if (length(n_cores) == 0 || is.na(n_cores)) return(1L)
	n_cores <- max(1L, as.integer(n_cores[1]))
	
	# Cap the request at the number of available logical CPUs when detectable.
	detected <- parallel::detectCores(logical = TRUE)
	if (!is.na(detected)) n_cores <- min(n_cores, detected)
	
	n_cores
}

#' Apply a function in parallel and row-bind tibble results.
#'
#' This helper centralizes task-level parallelization for enrichment jobs. It
#' uses `future.apply` with a multisession backend when available, falls back to
#' `parallel::mclapply()` on Unix-like systems, and otherwise runs serially.
#'
#' @param x Vector or list of task identifiers.
#' @param .f Function applied to each element of `x`.
#' @param n_cores Number of worker processes to use. Use `1` for serial execution.
#'
#' @return A tibble created by row-binding the per-task results.
parallel_map_dfr <- function(x, .f, n_cores = 1L) {
	# Normalize the CPU request so downstream backends receive a valid value.
	n_cores <- normalize_n_cores(n_cores)
	
	# Preserve original serial behavior unless parallel execution is requested.
	if (n_cores <= 1L || length(x) <= 1L) {
		return(purrr::map_dfr(x, .f))
	}
	
	# Prefer future.apply for cross-platform parallel workers when installed.
	if (is_pkg_available("future.apply") && is_pkg_available("future")) {
		old_plan <- future::plan()
		on.exit(future::plan(old_plan), add = TRUE)
		
		# multisession works on Windows, macOS, and Linux.
		future::plan(future::multisession, workers = n_cores)
		
		out <- future.apply::future_lapply(
			x,
			.f,
			future.seed = TRUE,
			future.packages = c("dplyr", "tidyr", "purrr", "stringr", "readr", "tibble", "rlang")
		)
		return(dplyr::bind_rows(out))
	}
	
	# On Unix-like systems, forking is available without extra packages.
	if (.Platform$OS.type != "windows") {
		out <- parallel::mclapply(x, .f, mc.cores = n_cores)
		return(dplyr::bind_rows(out))
	}
	
	# Windows without future.apply cannot fork, so fall back gracefully.
	warning(
		"Parallel execution on Windows requires the 'future.apply' and 'future' packages. ",
		"Falling back to serial execution.",
		call. = FALSE
	)
	purrr::map_dfr(x, .f)
}

# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------

#' Validate and filter MOFA view names.
#'
#' Compares the available MOFA view names with the set supported by this script
#' and returns only the supported views. Unsupported views are ignored with a warning.
#'
#' @param views Character vector of view names extracted from, or supplied for, a MOFA object.
#'
#' @return Character vector containing supported view names.
validate_views <- function(views) {
	# Define the set of omics views handled by the downstream pipeline.
	supported <- c("proteomics", "methylation", "metabolomics", "expression")
	unsupported <- setdiff(unique(views), supported)
	if (length(unsupported) > 0) {
		warning(
				"Unsupported MOFA view(s) found and ignored: ",
				paste(unsupported, collapse = ", "),
				call. = FALSE
		)
	}
	intersect(unique(views), supported)
}

# -----------------------------------------------------------------------------
# MOFA loading extraction
# -----------------------------------------------------------------------------

#' Extract feature loadings from a fitted MOFA object.
#'
#' Uses `MOFA2::get_weights()` to retrieve factor loadings, normalizes column
#' names across possible MOFA2 output variants, and retains only supported omics views.
#'
#' @param mofa_obj A fitted MOFA/MOFA2 object.
#' @param views Character vector of view names to extract.
#' @param factors Factor selection passed to `MOFA2::get_weights()`; use `"all"` for all factors.
#'
#' @return Tibble with columns `view`, `factor`, `feature`, and `value`.
extract_mofa_loadings <- function(mofa_obj,
		views = c("proteomics", "methylation", "metabolomics", "expression"),
		factors = "all") {
	# Loadings are the core inputs used by all downstream feature selection steps.
	require_namespace("MOFA2", "extracting MOFA feature loadings")
	
	# MOFA2::get_weights generally supports as.data.frame = TRUE.
	w <- MOFA2::get_weights(
			object = mofa_obj,
			views = views,
			factors = factors,
			as.data.frame = TRUE
	)
	
	w <- as_tibble(w)
	
	# Make this robust to small naming differences across MOFA2 versions.
	nm <- names(w)
	if (!"view" %in% nm && "View" %in% nm) w <- rename(w, view = View)
	if (!"factor" %in% nm && "Factor" %in% nm) w <- rename(w, factor = Factor)
	if (!"feature" %in% nm && "Feature" %in% nm) w <- rename(w, feature = Feature)
	if (!"value" %in% nm && "weight" %in% nm) w <- rename(w, value = weight)
	if (!"value" %in% nm && "Weight" %in% nm) w <- rename(w, value = Weight)
	
	needed <- c("view", "factor", "feature", "value")
	missing_cols <- setdiff(needed, names(w))
	if (length(missing_cols) > 0) {
		stop(
				"Could not find required loading columns: ",
				paste(missing_cols, collapse = ", "),
				". Inspect MOFA2::get_weights(..., as.data.frame = TRUE) output.",
				call. = FALSE
		)
	}
	
	valid_views <- validate_views(w$view)
	
	w %>%
			filter(view %in% valid_views) %>%
			mutate(
					view = as.character(view),
					factor = as.character(factor),
					feature = as.character(feature),
					value = as.numeric(value)
			) %>%
			filter(!is.na(value))
}

# -----------------------------------------------------------------------------
# Optional MOFA variance explained extraction
# -----------------------------------------------------------------------------

#' Extract variance explained from a fitted MOFA object when available.
#'
#' Attempts to retrieve the per-factor, per-view variance explained using
#' `MOFA2::get_variance_explained()`. If extraction fails or the result format
#' is not recognized, an empty tibble is returned.
#'
#' @param mofa_obj A fitted MOFA/MOFA2 object.
#'
#' @return Tibble with `factor`, `view`, and `variance_explained`, or an empty tibble.
extract_variance_explained <- function(mofa_obj) {
	# Variance-explained extraction is optional because MOFA object formats can vary.
	if (!is_pkg_available("MOFA2")) return(tibble())
	
	out <- tryCatch({
				ve <- MOFA2::get_variance_explained(mofa_obj)
				
				if (is.list(ve) && "r2_per_factor" %in% names(ve)) {
					x <- ve$r2_per_factor
				} else {
					x <- ve
				}
				
				if (is.matrix(x) || is.data.frame(x)) {
					as.data.frame(x) %>%
							rownames_to_column("factor") %>%
							pivot_longer(-factor, names_to = "view", values_to = "variance_explained") %>%
							mutate(
									factor = as.character(factor),
									view = as.character(view),
									variance_explained = as.numeric(variance_explained)
							)
				} else {
					tibble()
				}
			}, error = function(e) tibble())
	
	out
}

# -----------------------------------------------------------------------------
# Feature selection from factor loadings
# -----------------------------------------------------------------------------

#' Select high-loading features for each factor, view, and direction.
#'
#' Identifies the positive and negative features most strongly associated with
#' each MOFA factor. The selection can be based on top absolute loading rank,
#' loading z-score, or absolute-loading quantile.
#'
#' @param loadings Tibble produced by `extract_mofa_loadings()`.
#' @param selection_mode Selection rule: `"top_n"`, `"zscore"`, or `"quantile"`.
#' @param top_n Number of features selected per direction when using `"top_n"`.
#' @param z_cutoff Absolute z-score threshold when using `"zscore"`.
#' @param quantile_cutoff Absolute-loading quantile when using `"quantile"`.
#' @param min_features Minimum number of features required to retain a group.
#'
#' @return Tibble of selected feature loadings.
select_factor_features <- function(loadings,
		selection_mode = c("top_n", "zscore", "quantile"),
		top_n = 100,
		z_cutoff = 2.5,
		quantile_cutoff = 0.95,
		min_features = 5) {
	# Select features separately for positive and negative loading directions.
	selection_mode <- match.arg(selection_mode)
	
	loadings <- loadings %>%
			group_by(view, factor) %>%
			mutate(
					abs_value = abs(value),
					z_value = ifelse(sd(value, na.rm = TRUE) > 0,
							(value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE),
							NA_real_),
					direction = case_when(
							value > 0 ~ "positive",
							value < 0 ~ "negative",
							TRUE ~ "zero"
					)
			) %>%
			ungroup() %>%
			filter(direction %in% c("positive", "negative"))
	
	if (selection_mode == "top_n") {
		selected <- loadings %>%
				group_by(view, factor, direction) %>%
				arrange(desc(abs(value)), .by_group = TRUE) %>%
				slice_head(n = top_n) %>%
				ungroup()
	} else if (selection_mode == "zscore") {
		selected <- loadings %>%
				filter(abs(z_value) >= z_cutoff)
	} else {
		selected <- loadings %>%
				group_by(view, factor) %>%
				mutate(abs_threshold = stats::quantile(abs_value, probs = quantile_cutoff, na.rm = TRUE)) %>%
				ungroup() %>%
				filter(abs_value >= abs_threshold)
	}
	
	selected %>%
			group_by(view, factor, direction) %>%
			filter(n() >= min_features) %>%
			ungroup()
}

# -----------------------------------------------------------------------------
# Identifier mapping
# -----------------------------------------------------------------------------

#' Map UniProt protein identifiers to Entrez gene identifiers.
#'
#' Uses a Bioconductor organism annotation package through `AnnotationDbi` to
#' convert UniProt IDs from the proteomics view into Entrez IDs suitable for
#' gene-set enrichment.
#'
#' @param uniprot_ids Character vector of UniProt identifiers.
#' @param OrgDb Character string naming the organism annotation package.
#'
#' @return Tibble with `feature`, `ENTREZID`, and `SYMBOL`.
map_uniprot_to_entrez <- function(uniprot_ids,
		OrgDb = "org.Hs.eg.db") {
	# Convert proteomics identifiers into a gene namespace used by enrichment tools.
	require_namespace("AnnotationDbi", "UniProt-to-gene mapping")
	require_namespace(OrgDb, "organism annotation mapping")
	
	db <- getExportedValue(OrgDb, OrgDb)
	
	ids <- unique(na.omit(as.character(uniprot_ids)))
	if (length(ids) == 0) return(tibble(feature = character(), ENTREZID = character(), SYMBOL = character()))
	
	suppressMessages({
				mapped <- AnnotationDbi::select(
						db,
						keys = ids,
						keytype = "UNIPROT",
						columns = c("ENTREZID", "SYMBOL")
				)
			})
	
	mapped %>%
			as_tibble() %>%
			rename(feature = UNIPROT) %>%
			distinct(feature, ENTREZID, SYMBOL) %>%
			filter(!is.na(ENTREZID))
}

#' Map expression feature identifiers to Entrez gene identifiers.
#'
#' Converts expression features supplied as Entrez IDs or gene symbols to Entrez
#' identifiers. In `"auto"` mode, mostly numeric IDs are treated as Entrez IDs
#' and other IDs are treated as symbols.
#'
#' @param expression_ids Character vector of expression feature identifiers.
#' @param expression_id_type Identifier type: `"auto"`, `"ENTREZID"`, or `"SYMBOL"`.
#' @param OrgDb Character string naming the organism annotation package.
#'
#' @return Tibble with `feature`, `ENTREZID`, and `SYMBOL`.
map_expression_to_entrez <- function(expression_ids,
		expression_id_type = c("auto", "ENTREZID", "SYMBOL"),
		OrgDb = "org.Hs.eg.db") {
	# Normalize expression identifiers into Entrez IDs for enrichment analysis.
	expression_id_type <- match.arg(expression_id_type)
	
	require_namespace("AnnotationDbi", "expression identifier mapping")
	require_namespace(OrgDb, "organism annotation mapping")
	
	db <- getExportedValue(OrgDb, OrgDb)
	
	ids <- unique(na.omit(as.character(expression_ids)))
	if (length(ids) == 0) return(tibble(feature = character(), ENTREZID = character(), SYMBOL = character()))
	
	if (expression_id_type == "auto") {
		numeric_fraction <- mean(str_detect(ids, "^[0-9]+$"))
		expression_id_type <- ifelse(numeric_fraction > 0.8, "ENTREZID", "SYMBOL")
	}
	
	if (expression_id_type == "ENTREZID") {
		suppressMessages({
					mapped <- AnnotationDbi::select(
							db,
							keys = ids,
							keytype = "ENTREZID",
							columns = c("ENTREZID", "SYMBOL")
					)
				})
		mapped <- mapped %>% rename(feature = ENTREZID)
	} else {
		suppressMessages({
					mapped <- AnnotationDbi::select(
							db,
							keys = ids,
							keytype = "SYMBOL",
							columns = c("ENTREZID", "SYMBOL")
					)
				})
		mapped <- mapped %>% rename(feature = SYMBOL)
	}
	
	mapped %>%
			as_tibble() %>%
			distinct(feature, ENTREZID, SYMBOL) %>%
			filter(!is.na(ENTREZID))
}

# -----------------------------------------------------------------------------
# Gene-set enrichment for expression/proteomics
# -----------------------------------------------------------------------------

#' Run Gene Ontology enrichment for a selected gene set.
#'
#' Performs over-representation analysis with `clusterProfiler::enrichGO()` using
#' a user-supplied background universe from the same MOFA view.
#'
#' @param selected_entrez Entrez IDs selected for a factor/view/direction group.
#' @param background_entrez Entrez IDs defining the tested background universe.
#' @param ontology GO ontology to test, usually `"BP"`, `"MF"`, or `"CC"`.
#' @param OrgDb Character string naming the organism annotation package.
#' @param p_adjust_method Multiple-testing correction method.
#' @param qvalue_cutoff q-value cutoff passed to `enrichGO()`.
#' @param min_gene_set_size Minimum gene-set size to test.
#' @param max_gene_set_size Maximum gene-set size to test.
#'
#' @return Tibble of GO enrichment results, or an empty tibble.
run_go_enrichment <- function(selected_entrez,
		background_entrez,
		ontology = "BP",
		OrgDb = "org.Hs.eg.db",
		p_adjust_method = "BH",
		qvalue_cutoff = 0.25,
		min_gene_set_size = 10,
		max_gene_set_size = 500) {
	# Run enrichment only after intersecting selected genes with the background universe.
	require_namespace("clusterProfiler", "GO enrichment")
	require_namespace(OrgDb, "GO enrichment")
	
	db <- getExportedValue(OrgDb, OrgDb)
	selected_entrez <- unique(na.omit(as.character(selected_entrez)))
	background_entrez <- unique(na.omit(as.character(background_entrez)))
	selected_entrez <- intersect(selected_entrez, background_entrez)
		
	ego <- tryCatch({
				clusterProfiler::enrichGO(
						gene = selected_entrez,
						universe = background_entrez,
						OrgDb = db,
						keyType = "ENTREZID",
						ont = ontology,
						pAdjustMethod = p_adjust_method,
						qvalueCutoff = qvalue_cutoff,
						minGSSize = min_gene_set_size,
						maxGSSize = max_gene_set_size,
						readable = TRUE
				)
			}, error = function(e) NULL)
	
	if (is.null(ego)) return(tibble())
	
	as.data.frame(ego) %>%
			as_tibble() %>%
			mutate(database = paste0("GO_", ontology))
}

#' Run Reactome pathway enrichment for a selected gene set.
#'
#' Performs Reactome over-representation analysis when `ReactomePA` is installed.
#' If the package is unavailable or enrichment fails, an empty tibble is returned.
#'
#' @param selected_entrez Entrez IDs selected for a factor/view/direction group.
#' @param background_entrez Entrez IDs defining the tested background universe.
#' @param organism Organism name passed to `ReactomePA::enrichPathway()`.
#' @param p_adjust_method Multiple-testing correction method.
#' @param qvalue_cutoff q-value cutoff passed to `enrichPathway()`.
#' @param min_gene_set_size Minimum pathway size to test.
#' @param max_gene_set_size Maximum pathway size to test.
#'
#' @return Tibble of Reactome enrichment results, or an empty tibble.
run_reactome_enrichment <- function(selected_entrez,
		background_entrez,
		organism = "human",
		p_adjust_method = "BH",
		qvalue_cutoff = 0.25,
		min_gene_set_size = 10,
		max_gene_set_size = 500) {
	# Reactome enrichment is optional and skipped if ReactomePA is unavailable.
	if (!is_pkg_available("ReactomePA")) return(tibble())
	
	selected_entrez <- unique(na.omit(as.character(selected_entrez)))
	background_entrez <- unique(na.omit(as.character(background_entrez)))
	selected_entrez <- intersect(selected_entrez, background_entrez)
		
	er <- tryCatch({
				ReactomePA::enrichPathway(
						gene = selected_entrez,
						universe = background_entrez,
						organism = organism,
						pAdjustMethod = p_adjust_method,
						qvalueCutoff = qvalue_cutoff,
						minGSSize = min_gene_set_size,
						maxGSSize = max_gene_set_size,
						readable = TRUE
				)
			}, error = function(e) NULL)
	
	if (is.null(er)) return(tibble())
	
	as.data.frame(er) %>%
			as_tibble() %>%
			mutate(database = "Reactome")
}

#' Run bundled gene-set enrichment analyses.
#'
#' Combines GO enrichment and, when available, Reactome enrichment for the same
#' selected Entrez IDs and background universe.
#'
#' @param selected_entrez Entrez IDs selected for testing.
#' @param background_entrez Entrez IDs defining the enrichment background.
#' @param OrgDb Character string naming the organism annotation package.
#' @param organism Organism name for Reactome enrichment.
#' @param ontology GO ontology to test.
#' @param qvalue_cutoff q-value cutoff used by enrichment functions.
#'
#' @return Tibble combining GO and Reactome enrichment results.
run_gene_enrichment_bundle <- function(selected_entrez,
		background_entrez,
		OrgDb = "org.Hs.eg.db",
		organism = "human",
		ontology = "BP",
		qvalue_cutoff = 0.25) {
	# Combine complementary gene-set resources into one results table.
	bind_rows(
					run_go_enrichment(
							selected_entrez = selected_entrez,
							background_entrez = background_entrez,
							ontology = ontology,
							OrgDb = OrgDb,
							qvalue_cutoff = qvalue_cutoff
					),
					run_reactome_enrichment(
							selected_entrez = selected_entrez,
							background_entrez = background_entrez,
							organism = organism,
							qvalue_cutoff = qvalue_cutoff
					)
			) %>%
			mutate(across(any_of(c("pvalue", "p.adjust", "qvalue")), as.numeric)) %>%
			arrange(p.adjust)
}


# -----------------------------------------------------------------------------
# Ranked gene-set enrichment for expression/proteomics with fgsea
# -----------------------------------------------------------------------------

#' Collapse mapped feature loadings to one ranked score per Entrez gene.
#'
#' `fgsea` requires a named numeric vector with one score per gene. This helper
#' resolves many-to-one mappings, such as multiple UniProt accessions mapping to
#' the same Entrez gene, by retaining the feature with the largest absolute MOFA
#' loading for each Entrez ID.
#'
#' @param mapped_loadings Tibble containing at least `ENTREZID` and `value`.
#'
#' @return Named numeric vector of signed MOFA loadings, sorted decreasingly.
collapse_ranked_gene_scores <- function(mapped_loadings) {
	# Keep only valid gene IDs and finite signed loading values.
	collapsed <- mapped_loadings %>%
			filter(!is.na(ENTREZID), !is.na(value), is.finite(value)) %>%
			mutate(ENTREZID = as.character(ENTREZID)) %>%
			group_by(ENTREZID) %>%
			arrange(desc(abs(value)), .by_group = TRUE) %>%
			slice_head(n = 1) %>%
			ungroup()
	
	# Convert the collapsed tibble to the named vector expected by fgsea.
	stats <- collapsed$value
	names(stats) <- collapsed$ENTREZID
	stats <- stats[!duplicated(names(stats))]
	stats <- sort(stats, decreasing = TRUE)
	stats
}

#' Build GO pathway lists for fgsea from an organism annotation package.
#'
#' This function derives GO gene sets directly from the same `OrgDb` namespace
#' used for ID mapping. The resulting list is restricted to the requested GO
#' ontology and can be passed directly to `fgsea::fgsea()`.
#'
#' @param OrgDb Character string naming the organism annotation package.
#' @param ontology GO ontology to include, usually `"BP"`, `"MF"`, or `"CC"`.
#'
#' @return Named list mapping GO IDs to Entrez ID character vectors.
build_go_pathways_for_fgsea <- function(OrgDb = "org.Hs.eg.db",
		ontology = "BP") {
	# GO gene sets are constructed from the organism annotation package.
	require_namespace("AnnotationDbi", "building GO gene sets for fgsea")
	require_namespace(OrgDb, "building GO gene sets for fgsea")
	
	db <- getExportedValue(OrgDb, OrgDb)
	all_entrez <- AnnotationDbi::keys(db, keytype = "ENTREZID")
	
	# Prefer GOALL/ONTOLOGYALL when available so inherited GO annotations are used.
	go_map <- suppressMessages({
		AnnotationDbi::select(
				db,
				keys = all_entrez,
				keytype = "ENTREZID",
				columns = c("GOALL", "ONTOLOGYALL")
		)
	})
	
	go_map %>%
			as_tibble() %>%
			filter(!is.na(GOALL), !is.na(ONTOLOGYALL)) %>%
			filter(ONTOLOGYALL == ontology) %>%
			distinct(GOALL, ENTREZID) %>%
			split(.$GOALL) %>%
			purrr::map(~ unique(as.character(.x$ENTREZID)))
}

#' Retrieve GO term descriptions when GO.db is available.
#'
#' `fgsea` returns pathway identifiers. This helper adds readable GO term names
#' without making GO.db a hard dependency for the entire script.
#'
#' @param go_ids Character vector of GO identifiers.
#'
#' @return Named character vector with GO IDs as names and descriptions as values.
lookup_go_descriptions <- function(go_ids) {
	# Fall back to the GO ID itself if GO.db is not installed.
	go_ids <- unique(na.omit(as.character(go_ids)))
	if (length(go_ids) == 0) return(character())
	if (!is_pkg_available("GO.db")) return(setNames(go_ids, go_ids))
	
	go_terms <- tryCatch({
		AnnotationDbi::select(
				GO.db::GO.db,
				keys = go_ids,
				keytype = "GOID",
				columns = c("TERM")
		) %>%
				as_tibble() %>%
				filter(!is.na(GOID), !is.na(TERM)) %>%
				distinct(GOID, TERM)
	}, error = function(e) tibble(GOID = go_ids, TERM = go_ids))
	
	out <- setNames(go_terms$TERM, go_terms$GOID)
	missing_ids <- setdiff(go_ids, names(out))
	c(out, setNames(missing_ids, missing_ids))
}

#' Run fgsea over GO gene sets using signed MOFA loadings.
#'
#' This ranked enrichment method is used for proteomics, expression, and mapped methylation views.
#' Unlike over-representation analysis, it uses the full signed loading vector
#' rather than an arbitrary top-feature cutoff. Positive NES values indicate GO
#' terms concentrated near the positive-loading end of the ranked vector; negative NES
#' values indicate GO terms concentrated near the negative-loading end. Report
#' tables label these fgsea results as `"ranked"`, because each result is
#' produced from one full signed ranked vector rather than separate positive
#' and negative analyses.
#'
#' @param mapped_loadings Tibble containing mapped Entrez IDs and MOFA loadings.
#' @param OrgDb Character string naming the organism annotation package.
#' @param ontology GO ontology to test.
#' @param p_adjust_method Multiple-testing correction method applied to p-values.
#' @param qvalue_cutoff Adjusted p-value threshold used only for optional filtering.
#' @param min_gene_set_size Minimum gene-set size to test.
#' @param max_gene_set_size Maximum gene-set size to test.
#' @param fgsea_nproc Number of CPUs passed to `fgsea::fgsea()` for within-task
#'   parallelization. When `run_all_enrichments()` is already using multiple
#'   task-level workers, keep this at `1` to avoid oversubscribing CPUs.
#'
#' @return Tibble of fgsea GO results, or an empty tibble.
run_fgsea_go_enrichment <- function(mapped_loadings,
		OrgDb = "org.Hs.eg.db",
		ontology = "BP",
		p_adjust_method = "BH",
		qvalue_cutoff = 0.25,
		min_gene_set_size = 10,
		max_gene_set_size = 500,
		fgsea_nproc = 1L) {
	# fgsea is the primary ranked enrichment method for gene/protein views.
	require_namespace("fgsea", "ranked GO enrichment")
	
	# Validate the requested within-fgsea CPU count.
	fgsea_nproc <- normalize_n_cores(fgsea_nproc)
	
	# Convert feature-level signed loadings to one ranked score per Entrez gene.
	stats <- collapse_ranked_gene_scores(mapped_loadings)
	# fgsea cannot test gene sets larger than the available ranked universe.
	# If the ranked vector is smaller than the requested minimum gene-set size,
	# no pathway can satisfy the size filter.
	if (length(stats) < min_gene_set_size) return(tibble())
	
	# Build GO gene sets and restrict them to genes observed in this MOFA view.
	pathways <- build_go_pathways_for_fgsea(OrgDb = OrgDb, ontology = ontology)
	pathways <- purrr::map(pathways, ~ intersect(.x, names(stats)))
	pathways <- pathways[lengths(pathways) >= min_gene_set_size & lengths(pathways) <= max_gene_set_size]
	if (length(pathways) == 0) return(tibble())
	
	# Run preranked enrichment using signed MOFA loadings.
	fg <- tryCatch({
		fgsea::fgsea(
				pathways = pathways,
				stats = stats,
				minSize = min_gene_set_size,
				maxSize = max_gene_set_size,
				nproc = fgsea_nproc
		)
	}, error = function(e) NULL)
	
	if (is.null(fg) || nrow(fg) == 0) return(tibble())
	
	fg <- as_tibble(fg)
	go_descriptions <- lookup_go_descriptions(fg$pathway)
	
	fg %>%
			mutate(
					database = paste0("GO_", ontology, "_fgsea"),
					ID = as.character(pathway),
					Description = unname(go_descriptions[ID]),
					pvalue = as.numeric(pval),
					p.adjust = p.adjust(pvalue, method = p_adjust_method),
					qvalue = p.adjust,
					Count = as.integer(size),
					GeneRatio = NA_character_,
					BgRatio = NA_character_,
					geneID = purrr::map_chr(leadingEdge, ~ paste(.x, collapse = "/")),
					direction = "ranked"
			) %>%
			filter(p.adjust <= qvalue_cutoff | is.na(qvalue_cutoff)) %>%
			arrange(p.adjust)
}

#' Run ranked gene-set enrichment for a mapped MOFA view.
#'
#' This wrapper standardizes the interface used by `run_all_enrichments()` for
#' expression and proteomics views.
#'
#' @param mapped_loadings Tibble containing feature IDs, Entrez IDs, and loadings.
#' @param OrgDb Character string naming the organism annotation package.
#' @param ontology GO ontology to test.
#' @param qvalue_cutoff Adjusted p-value threshold for reported fgsea results.
#' @param fgsea_nproc Number of CPUs passed to `fgsea::fgsea()`.
#'
#' @return Tibble of ranked GO enrichment results.
run_gene_fgsea_bundle <- function(mapped_loadings,
		OrgDb = "org.Hs.eg.db",
		ontology = "BP",
		qvalue_cutoff = 0.25,
		fgsea_nproc = 1L) {
	# Keep this wrapper separate so additional ranked gene-set collections can be
	# added later without changing the enrichment dispatcher.
	run_fgsea_go_enrichment(
			mapped_loadings = mapped_loadings,
			OrgDb = OrgDb,
			ontology = ontology,
			qvalue_cutoff = qvalue_cutoff,
			fgsea_nproc = fgsea_nproc
	)
}


# -----------------------------------------------------------------------------
# Ranked methylation enrichment with fgsea
# -----------------------------------------------------------------------------

#' Map methylation CpG loadings to Entrez gene-level loadings.
#'
#' This helper prepares methylation factor loadings for ranked enrichment. It
#' uses `missMethyl::getMappedEntrezIDs()` to apply missMethyl's CpG-to-gene
#' annotation logic to the CpG universe being analyzed. Because `fgsea` requires
#' a signed score for each gene, the helper then uses missMethyl's flattened
#' annotation table to retain the explicit CpG-to-Entrez mapping and joins the
#' original MOFA CpG loadings back to the mapped genes.
#'
#' When multiple CpGs map to the same Entrez gene, downstream ranked enrichment
#' collapses those mappings with `collapse_ranked_gene_scores()`, retaining the
#' CpG with the largest absolute MOFA loading for that gene.
#'
#' @param cpg_loadings Tibble containing at least `feature` and `value` columns,
#'   where `feature` contains CpG probe IDs and `value` contains signed MOFA
#'   loadings for one factor.
#' @param all_cpgs Character vector containing the CpG universe for the tested
#'   methylation view. This is passed to `missMethyl::getMappedEntrezIDs()`.
#' @param array_type Illumina methylation array type understood by missMethyl,
#'   for example `"EPIC"` or `"450K"` in this script.
#'
#' @return Tibble with CpG features, Entrez IDs, and signed MOFA loading values.
map_methylation_to_entrez_loadings <- function(cpg_loadings,
		all_cpgs = NULL,
		array_type = c("EPIC", "450K")) {
	# Match the supported missMethyl array type before calling missMethyl helpers.
	array_type <- match.arg(array_type)
	
	# missMethyl supplies the array-aware CpG-to-gene annotation machinery.
	require_namespace("missMethyl", "mapping methylation CpGs to Entrez IDs")
	
	# Keep only valid CpG features and signed loading values for the current factor.
	cpg_loadings <- cpg_loadings %>%
			select(feature, value) %>%
			mutate(feature = as.character(feature), value = as.numeric(value)) %>%
			filter(!is.na(feature), !is.na(value), is.finite(value))
	
	if (nrow(cpg_loadings) == 0) return(tibble())
	
	# Define the tested CpG universe. In this pipeline this should usually be all
	# methylation features present in the MOFA view for the current factor.
	if (is.null(all_cpgs)) all_cpgs <- cpg_loadings$feature
	all_cpgs <- unique(na.omit(as.character(all_cpgs)))
	
	# Use getMappedEntrezIDs() so that CpG filtering and supported annotation
	# behavior follow missMethyl's own enrichment mapping logic.
	mapped_ids <- tryCatch({
		missMethyl::getMappedEntrezIDs(
				sig.cpg = cpg_loadings$feature,
				all.cpg = all_cpgs,
				array.type = array_type
		)
	}, error = function(e) NULL)
	
	# If missMethyl cannot map these CpGs, fgsea cannot be run for this factor.
	if (is.null(mapped_ids) || length(mapped_ids$sig.eg) == 0) return(tibble())
	
	# getMappedEntrezIDs() returns gene sets but not a direct CpG-loading table.
	# The flattened annotation gives the required CpG -> Entrez relationship.
	flat <- tryCatch({
		missMethyl:::.getFlatAnnotation(array_type)
	}, error = function(e) NULL)
	
	if (is.null(flat) || nrow(flat) == 0) return(tibble())
	
	flat_tbl <- flat %>%
			as.data.frame(stringsAsFactors = FALSE) %>%
			as_tibble() %>%
			transmute(
					feature = as.character(cpg),
					ENTREZID = as.character(entrezid)
			) %>%
			filter(!is.na(feature), !is.na(ENTREZID), feature %in% cpg_loadings$feature) %>%
			distinct(feature, ENTREZID)
	
	# Restrict the explicit mapping to Entrez IDs that missMethyl reports as
	# mappable for the CpGs supplied to getMappedEntrezIDs().
	valid_entrez <- unique(as.character(mapped_ids$sig.eg))
	flat_tbl <- flat_tbl %>% filter(ENTREZID %in% valid_entrez)
	
	if (nrow(flat_tbl) == 0) return(tibble())
	
	# Join the original signed CpG loadings to the mapped Entrez IDs. The resulting
	# table has the same structure used by expression/proteomics fgsea helpers.
	flat_tbl %>%
			left_join(cpg_loadings, by = "feature") %>%
			filter(!is.na(value), is.finite(value))
}

#' Prepare a CpG-to-Entrez mapping for methylation fgsea.
#'
#' This helper performs the missMethyl CpG-to-gene annotation step once, before
#' task-level parallelization begins. The resulting plain tibble can then be
#' safely passed to parallel workers, avoiding repeated calls to missMethyl
#' annotation internals inside independent worker R sessions.
#'
#' @param all_cpgs Character vector containing the CpG universe for the
#'   methylation view.
#' @param array_type Illumina methylation array type understood by missMethyl,
#'   for example `"EPIC"` or `"450K"` in this script.
#'
#' @return Tibble with columns `feature` and `ENTREZID`.
prepare_methylation_entrez_mapping <- function(all_cpgs,
		array_type = c("EPIC", "450K")) {
	# Match the supported missMethyl array type in the main R process.
	array_type <- match.arg(array_type)
	
	# missMethyl supplies the array-aware CpG-to-gene annotation machinery.
	require_namespace("missMethyl", "precomputing methylation CpG-to-Entrez mapping")
	
	# Normalize the CpG universe and return a correctly shaped empty tibble if no
	# methylation features are available.
	all_cpgs <- unique(na.omit(as.character(all_cpgs)))
	if (length(all_cpgs) == 0) {
		return(tibble(feature = character(), ENTREZID = character()))
	}
	
	# Call getMappedEntrezIDs() once in the main process. Passing all CpGs as both
	# sig.cpg and all.cpg asks missMethyl to identify all mappable Entrez IDs in
	# the methylation universe used by the MOFA model.
	mapped_ids <- tryCatch({
		missMethyl::getMappedEntrezIDs(
				sig.cpg = all_cpgs,
				all.cpg = all_cpgs,
				array.type = array_type
		)
	}, error = function(e) NULL)
	
	if (is.null(mapped_ids)) {
		return(tibble(feature = character(), ENTREZID = character()))
	}
	
	# Retain the Entrez IDs that missMethyl regards as mappable for this CpG
	# universe. Different missMethyl versions expose both sig.eg and all.eg, so use
	# whichever are present.
	valid_entrez <- unique(na.omit(as.character(unlist(mapped_ids[c("sig.eg", "all.eg")], use.names = FALSE))))
	if (length(valid_entrez) == 0) {
		return(tibble(feature = character(), ENTREZID = character()))
	}
	
	# The flattened annotation provides the explicit CpG -> Entrez relationship
	# needed to attach signed MOFA loadings to mapped genes.
	flat <- tryCatch({
		missMethyl:::.getFlatAnnotation(array_type)
	}, error = function(e) NULL)
	
	if (is.null(flat) || nrow(flat) == 0) {
		return(tibble(feature = character(), ENTREZID = character()))
	}
	
	flat_tbl <- flat %>%
			as.data.frame(stringsAsFactors = FALSE) %>%
			as_tibble()
	
	# Defensively locate the CpG and Entrez columns. missMethyl currently returns
	# lower-case names, but this keeps the helper robust to capitalization changes.
	flat_names <- tolower(names(flat_tbl))
	cpg_col <- names(flat_tbl)[match("cpg", flat_names)]
	entrez_col <- names(flat_tbl)[match("entrezid", flat_names)]
	
	if (is.na(cpg_col) || is.na(entrez_col)) {
		return(tibble(feature = character(), ENTREZID = character()))
	}
	
	flat_tbl %>%
			transmute(
					feature = as.character(.data[[cpg_col]]),
					ENTREZID = as.character(.data[[entrez_col]])
			) %>%
			filter(!is.na(feature), !is.na(ENTREZID)) %>%
			filter(feature %in% all_cpgs, ENTREZID %in% valid_entrez) %>%
			distinct(feature, ENTREZID)
}

#' Join methylation loadings to a precomputed CpG-to-Entrez mapping.
#'
#' This helper is used inside enrichment workers. It does not call missMethyl;
#' it only joins one factor's CpG loadings to the mapping prepared by
#' `prepare_methylation_entrez_mapping()` in the main process.
#'
#' @param cpg_loadings Tibble containing at least `feature` and `value` columns.
#' @param methylation_entrez_mapping Tibble with columns `feature` and `ENTREZID`.
#'
#' @return Tibble with CpG features, Entrez IDs, and signed MOFA loading values.
map_methylation_loadings_with_precomputed_mapping <- function(cpg_loadings,
		methylation_entrez_mapping) {
	# Return a correctly shaped empty tibble if no precomputed mapping is available.
	if (is.null(methylation_entrez_mapping) || nrow(methylation_entrez_mapping) == 0) {
		return(tibble(feature = character(), ENTREZID = character(), value = numeric()))
	}
	
	# Keep only valid CpG features and signed loading values for the current factor.
	cpg_loadings <- cpg_loadings %>%
			select(feature, value) %>%
			mutate(feature = as.character(feature), value = as.numeric(value)) %>%
			filter(!is.na(feature), !is.na(value), is.finite(value))
	
	if (nrow(cpg_loadings) == 0) {
		return(tibble(feature = character(), ENTREZID = character(), value = numeric()))
	}
	
	# Join to the precomputed mapping rather than calling missMethyl inside a
	# parallel worker.
	methylation_entrez_mapping %>%
			inner_join(cpg_loadings, by = "feature") %>%
			filter(!is.na(ENTREZID), !is.na(value), is.finite(value))
}

# -----------------------------------------------------------------------------
# Methylation enrichment via missMethyl
# -----------------------------------------------------------------------------

#' Run methylation-aware enrichment for selected CpG probes.
#'
#' Uses `missMethyl::gometh()` to test whether selected CpGs are enriched for
#' genes in GO or KEGG categories while accounting for probe-number bias.
#'
#' @param selected_cpgs Character vector of selected CpG IDs.
#' @param background_cpgs Character vector of all tested CpG IDs in the MOFA view.
#' @param array_type Illumina array type, `"EPIC"` or `"450K"`.
#' @param collection Annotation collection, usually `"GO"` or `"KEGG"`.
#' @param sig_genes Logical flag passed to `missMethyl::gometh()`.
#' @param prior_prob Logical flag indicating whether to model probe-number bias.
#' @param p_adjust_method Multiple-testing correction method.
#' @param min_selected Minimum required number of selected CpG IDs.
#' @param min_background Minimum required number of all tested CpG IDs.
#'
#' @return Tibble of methylation enrichment results, or an empty tibble.
run_methylation_enrichment <- function(selected_cpgs,
		background_cpgs,
		array_type = c("EPIC", "450K"),
		collection = c("GO", "KEGG"),
		sig_genes = FALSE,
		prior_prob = TRUE,
		p_adjust_method = "BH",
		min_selected = 5,
		min_background = 50) {
	# Use a CpG-aware method to reduce bias from unequal probe coverage across genes.
	array_type <- match.arg(array_type)
	collection <- match.arg(collection)
	
	require_namespace("missMethyl", "methylation-aware enrichment")
	
	selected_cpgs <- unique(na.omit(as.character(selected_cpgs)))
	background_cpgs <- unique(na.omit(as.character(background_cpgs)))
	selected_cpgs <- intersect(selected_cpgs, background_cpgs)
	
	if (length(selected_cpgs) < min_selected || length(background_cpgs) < min_background) 
		return(tibble())
	
	gm <- tryCatch({
				missMethyl::gometh(
						sig.cpg = selected_cpgs,
						all.cpg = background_cpgs,
						collection = collection,
						array.type = array_type,
						sig.genes = sig_genes,
						prior.prob = prior_prob
				)
			}, error = function(e) NULL)
	
	if (is.null(gm) || nrow(gm) == 0) return(tibble())
	
	gm %>%
			rownames_to_column(var = "ROWID") %>%
			mutate(
					database = paste0("missMethyl_", collection),
					ID = coalesce(as.character(ROWID), as.character(TERM)),
					Description = as.character(TERM),
					pvalue = as.numeric(P.DE),
					p.adjust = p.adjust(pvalue, method = p_adjust_method),
					qvalue = p.adjust,
					Count = as.integer(DE),
					GeneRatio = paste0(DE, "/", N)
			) %>% 
			arrange(p.adjust) %>% 
			column_to_rownames("ROWID") %>% as_tibble()
}

# -----------------------------------------------------------------------------
# Metabolon categorical enrichment
# -----------------------------------------------------------------------------

#' Clean Metabolon annotation values.
#'
#' Converts all columns to character and treats literal periods as missing
#' values, matching common Metabolon annotation-file conventions.
#'
#' @param metabolon_annotation Data frame containing Metabolon annotations.
#'
#' @return Tibble with cleaned character columns.
clean_metabolon_annotation <- function(metabolon_annotation) {
	# Standardize Metabolon missing-value encodings before joining or testing terms.
	metabolon_annotation %>%
			as_tibble() %>%
			mutate(across(everything(), ~ na_if(as.character(.x), ".")))
}

#' Run categorical enrichment for Metabolon annotation terms.
#'
#' Tests whether selected metabolites are over-represented in Metabolon pathway
#' or classification categories using a hypergeometric test against the MOFA
#' metabolomics background universe.
#'
#' @param selected_features Character vector of selected metabolite features.
#' @param background_features Character vector of all metabolite features in the MOFA view.
#' @param annotation Cleaned or raw Metabolon annotation table.
#' @param feature_col Column in `annotation` matching MOFA metabolomics feature names.
#' @param term_cols Annotation columns to test as categorical terms.
#' @param min_term_size Minimum background count required for a term to be tested.
#' @param p_adjust_method Multiple-testing correction method.
#' @param min_selected Minimum required number of selected metabolite features.
#' @param min_background Minimum required number of all metabolite features in the MOFA view.
#'
#' @return Tibble of categorical enrichment results.
run_categorical_enrichment <- function(selected_features,
		background_features,
		annotation,
		feature_col,
		term_cols = c("SUPER_PATHWAY", "SUB_PATHWAY", "PARENT_CLASSIFICATION", "CLASSIFICATION"),
		min_term_size = 3,
		p_adjust_method = "BH",
		min_selected = 3,
		min_background = 10) {
	# Treat each Metabolon pathway/class label as a categorical feature set.
	selected_features <- unique(na.omit(as.character(selected_features)))
	background_features <- unique(na.omit(as.character(background_features)))
	selected_features <- intersect(selected_features, background_features)
	
	if (length(selected_features) < min_selected || 
			length(background_features) < min_background) 
		return(tibble())
	
	if (!feature_col %in% names(annotation)) {
		stop("Metabolon feature_col '", feature_col, "' not found in annotation.", call. = FALSE)
	}
	
	term_cols <- intersect(term_cols, names(annotation))
	if (length(term_cols) == 0) return(tibble())
	
	ann_long <- annotation %>%
			mutate(feature = as.character(.data[[feature_col]])) %>%
			filter(feature %in% background_features) %>%
			select(feature, all_of(term_cols)) %>%
			pivot_longer(
					cols = all_of(term_cols),
					names_to = "term_type",
					values_to = "term"
			) %>%
			filter(!is.na(term), term != "") %>%
			distinct(feature, term_type, term)
	
	if (nrow(ann_long) == 0) return(tibble())
	
	bg_n <- length(unique(ann_long$feature))
	sel_n <- length(unique(intersect(selected_features, ann_long$feature)))
	
	ann_long %>%
			group_by(term_type, term) %>%
			summarise(
					term_background_n = n_distinct(feature),
					term_selected_n = n_distinct(feature[feature %in% selected_features]),
					selected_features_in_term = paste(sort(unique(feature[feature %in% selected_features])), collapse = "; "),
					.groups = "drop"
			) %>%
			filter(term_background_n >= min_term_size, term_selected_n > 0) %>%
			mutate(
					pvalue = purrr::pmap_dbl(
							list(term_selected_n, term_background_n),
							function(k, m) {
								stats::phyper(
										q = k - 1,
										m = m,
										n = bg_n - m,
										k = sel_n,
										lower.tail = FALSE
								)
							}
					),
					p.adjust = p.adjust(pvalue, method = p_adjust_method),
					qvalue = p.adjust,
					database = paste0("Metabolon_", term_type),
					ID = paste(term_type, term, sep = ":"),
					Description = term,
					Count = term_selected_n,
					GeneRatio = paste0(term_selected_n, "/", sel_n),
					BgRatio = paste0(term_background_n, "/", bg_n),
					geneID = selected_features_in_term
			) %>%
			arrange(p.adjust)
}

# -----------------------------------------------------------------------------
# Top-feature tables
# -----------------------------------------------------------------------------

#' Create a table of top-loaded features for reporting.
#'
#' Selects the largest absolute loadings by factor, view, and direction and
#' optionally appends Metabolon annotations to metabolomics features.
#'
#' @param loadings Tibble of MOFA loadings.
#' @param n Number of top features to retain per factor/view/direction.
#' @param metabolon_annotation Optional Metabolon annotation table.
#' @param metabolomics_feature_col Column matching metabolomics features.
#'
#' @return Tibble of top-loaded features and optional annotations.
make_top_feature_table <- function(loadings,
		n = 20,
		metabolon_annotation = NULL,
		metabolomics_feature_col = "BIOCHEMICAL") {
	# Build a compact table of the most interpretable features for each factor.
	top <- loadings %>%
			mutate(direction = ifelse(value >= 0, "positive", "negative")) %>%
			group_by(view, factor, direction) %>%
			arrange(desc(abs(value)), .by_group = TRUE) %>%
			slice_head(n = n) %>%
			ungroup()
	
	if (!is.null(metabolon_annotation) && "metabolomics" %in% top$view) {
		ann <- clean_metabolon_annotation(metabolon_annotation)
		if (metabolomics_feature_col %in% names(ann)) {
			keep_cols <- intersect(
					c(metabolomics_feature_col, "SUPER_PATHWAY", "SUB_PATHWAY", "KEGG", "HMDB", "PUBCHEM", "CHEMSPIDER"),
					names(ann)
			)
			ann2 <- ann %>%
					select(all_of(keep_cols)) %>%
					rename(feature = all_of(metabolomics_feature_col)) %>%
					distinct(feature, .keep_all = TRUE)
			
			top <- top %>%
					left_join(ann2, by = "feature")
		}
	}
	
	top
}

# -----------------------------------------------------------------------------
# Main enrichment runner
# -----------------------------------------------------------------------------

#' Run all view-specific enrichment analyses.
#'
#' Iterates over each factor and omics view and dispatches to the appropriate
#' enrichment method. Expression, proteomics, and methylation are analyzed with
#' ranked GO enrichment using `fgsea` and the full signed loading vector after
#' mapping features to Entrez genes. Each of these views is tested once per
#' factor/view combination; the sign of the `fgsea` NES determines whether a
#' term is associated with positive or negative factor loadings. Metabolomics
#' retains top-feature over-representation against Metabolon annotation categories.
#'
#' @param loadings Tibble of all MOFA loadings.
#' @param selected_features Tibble of selected high-loading features, used by
#'   metabolomics enrichment and retained for output/reporting.
#' @param metabolon_annotation Optional Metabolon annotation table.
#' @param metabolomics_feature_col Column matching metabolomics features.
#' @param expression_id_type Expression ID type: `"auto"`, `"ENTREZID"`, or `"SYMBOL"`.
#' @param OrgDb Character string naming the organism annotation package.
#' @param organism Organism name retained for compatibility with previous ORA code.
#' @param go_ontology GO ontology for ranked expression, proteomics, and methylation enrichment.
#' @param methylation_array_type Illumina methylation array type.
#' @param methylation_collection Retained for compatibility with older ORA-based methylation
#'   enrichment code; ranked methylation enrichment uses `go_ontology`.
#' @param methylation_entrez_mapping Optional precomputed CpG-to-Entrez mapping.
#'   If `NULL`, the mapping is prepared once in the main process before
#'   task-level parallelization begins.
#' @param qvalue_cutoff Adjusted p-value cutoff for reported enrichment results.
#' @param metabolon_term_cols Metabolon columns used as categorical terms.
#' @param n_cores Number of task-level parallel workers. Each task is one
#'   factor/view/direction job for metabolomics or one factor/view ranked job
#'   for expression, proteomics, or methylation.
#' @param fgsea_nproc Number of CPUs passed to `fgsea::fgsea()` inside each
#'   gene/protein/methylation fgsea task. Use `1` when `n_cores > 1` to avoid CPU oversubscription.
#'
#' @return Tibble of enrichment results from all available views.
run_all_enrichments <- function(loadings,
		selected_features,
		metabolon_annotation = NULL,
		metabolomics_feature_col = "BIOCHEMICAL",
		expression_id_type = c("auto", "ENTREZID", "SYMBOL"),
		OrgDb = "org.Hs.eg.db",
		organism = "human",
		go_ontology = "BP",
		methylation_array_type = c("EPIC", "450K"),
		methylation_collection = "GO",
		methylation_entrez_mapping = NULL,
		qvalue_cutoff = 0.25,
		metabolon_term_cols = c("SUPER_PATHWAY", "SUB_PATHWAY", "PARENT_CLASSIFICATION", "CLASSIFICATION"),
		n_cores = 1L,
		fgsea_nproc = 1L) {
	# Dispatch enrichment analyses according to the omics view being processed.
	expression_id_type <- match.arg(expression_id_type)
	methylation_array_type <- match.arg(methylation_array_type)
	
	# Normalize parallelization parameters before dispatching tasks.
	n_cores <- normalize_n_cores(n_cores)
	fgsea_nproc <- normalize_n_cores(fgsea_nproc)
	
	# Prepare methylation CpG-to-Entrez mapping before task-level parallelization.
	# This avoids calling missMethyl annotation helpers from independent worker R
	# sessions, which can behave inconsistently across parallel workers.
	if (is.null(methylation_entrez_mapping) && any(loadings$view == "methylation")) {
		methylation_entrez_mapping <- prepare_methylation_entrez_mapping(
				all_cpgs = loadings %>% filter(view == "methylation") %>% pull(feature),
				array_type = methylation_array_type
		)
	}
	
	# Ensure that workers always receive a plain tibble with the expected columns.
	if (is.null(methylation_entrez_mapping)) {
		methylation_entrez_mapping <- tibble(feature = character(), ENTREZID = character())
	} else {
		methylation_entrez_mapping <- methylation_entrez_mapping %>%
				as_tibble() %>%
				transmute(feature = as.character(feature), ENTREZID = as.character(ENTREZID)) %>%
				filter(!is.na(feature), !is.na(ENTREZID)) %>%
				distinct(feature, ENTREZID)
	}
	
	# Ranked gene/protein/methylation analyses use the full signed loading vector
	# exactly once per factor and view. fgsea then assigns positive or negative
	# direction from the sign of the NES, rather than from separate input subsets.
	ranked_tasks <- loadings %>%
			filter(view %in% c("proteomics", "expression", "methylation")) %>%
			distinct(view, factor) %>%
			mutate(direction = "ranked")
	
	# Metabolomics still uses selected top-loading features and therefore retains
	# explicit positive/negative selected-feature groups because it is analyzed by
	# ORA rather than ranked enrichment.
	selected_tasks <- selected_features %>%
			filter(!view %in% c("proteomics", "expression", "methylation")) %>%
			distinct(view, factor, direction)
	
	tasks <- bind_rows(ranked_tasks, selected_tasks) %>%
			arrange(factor, view, direction)
	
	if (nrow(tasks) == 0) return(tibble())
	
	# Process independent enrichment jobs in parallel when n_cores > 1.
	parallel_map_dfr(seq_len(nrow(tasks)), function(i) {
				task <- tasks[i, ]
				view_i <- task$view
				factor_i <- task$factor
				direction_i <- task$direction
				
				selected_i <- selected_features %>%
						filter(view == view_i, factor == factor_i, direction == direction_i)
				
				background_i <- loadings %>%
						filter(view == view_i, factor == factor_i)
				
				result <- tryCatch({
							if (view_i == "proteomics") {
								# Map every protein loading in this factor to Entrez IDs, then run fgsea.
								mapped_loadings <- map_uniprot_to_entrez(background_i$feature, OrgDb = OrgDb) %>%
										left_join(background_i, by = "feature")
								
								run_gene_fgsea_bundle(
										mapped_loadings = mapped_loadings,
										OrgDb = OrgDb,
										ontology = go_ontology,
										qvalue_cutoff = qvalue_cutoff,
										fgsea_nproc = fgsea_nproc
								)
							} else if (view_i == "expression") {
								# Map every expression loading in this factor to Entrez IDs, then run fgsea.
								mapped_loadings <- map_expression_to_entrez(background_i$feature, expression_id_type = expression_id_type, OrgDb = OrgDb) %>%
										left_join(background_i, by = "feature")
								
								run_gene_fgsea_bundle(
										mapped_loadings = mapped_loadings,
										OrgDb = OrgDb,
										ontology = go_ontology,
										qvalue_cutoff = qvalue_cutoff,
										fgsea_nproc = fgsea_nproc
								)
							} else if (view_i == "methylation") {
								# Join this factor's CpG loadings to the CpG-to-Entrez mapping that
								# was prepared once before parallel workers were launched.
								mapped_loadings <- map_methylation_loadings_with_precomputed_mapping(
										cpg_loadings = background_i,
										methylation_entrez_mapping = methylation_entrez_mapping
								)
								
								run_gene_fgsea_bundle(
										mapped_loadings = mapped_loadings,
										OrgDb = OrgDb,
										ontology = go_ontology,
										qvalue_cutoff = qvalue_cutoff,
										fgsea_nproc = fgsea_nproc
								)
							} else if (view_i == "metabolomics") {
								# Metabolomics enrichment remains categorical over-representation against Metabolon terms.
								if (is.null(metabolon_annotation)) return(tibble())
								ann <- clean_metabolon_annotation(metabolon_annotation)
								
								run_categorical_enrichment(
										selected_features = selected_i$feature,
										background_features = background_i$feature,
										annotation = ann,
										feature_col = metabolomics_feature_col,
										term_cols = metabolon_term_cols
								)
							} else {
								tibble()
							}
						}, error = function(e) {
							tibble(
									database = paste0("ERROR_", view_i),
									ID = NA_character_,
									Description = conditionMessage(e),
									pvalue = NA_real_,
									p.adjust = NA_real_,
									qvalue = NA_real_,
									Count = NA_integer_,
									GeneRatio = NA_character_
							)
						})
				
				if (nrow(result) == 0) return(tibble())
				
				# fgsea results are labeled as "ranked" because they come from one full
				# signed loading vector. ORA results inherit the selected-feature direction
				# from the task table.
				if (!"direction" %in% names(result)) {
					result <- result %>% mutate(direction = direction_i)
				}
				
				result %>%
						mutate(
								factor = factor_i,
								view = view_i,
								.before = 1
						) %>%
						relocate(factor, view, direction)
			}, n_cores = n_cores)
}

# -----------------------------------------------------------------------------
# Cautious biological interpretation text
# -----------------------------------------------------------------------------

#' Summarize enrichment results into factor-level interpretation text.
#'
#' Converts significant enrichment terms for a single factor into cautious,
#' human-readable prose for inclusion in the report.
#'
#' @param enrichment_results Tibble produced by `run_all_enrichments()`.
#' @param factor_id Factor identifier to summarize.
#' @param fdr_cutoff FDR threshold for including enrichment terms in prose.
#' @param max_terms_per_view Maximum number of terms per view/direction to include.
#'
#' @return Character string containing the factor-level interpretation.
summarise_factor_interpretation <- function(enrichment_results,
		factor_id,
		fdr_cutoff = 0.10,
		max_terms_per_view = 3) {
	# Convert enriched terms into conservative prose for the report.
	er <- enrichment_results %>%
			filter(factor == factor_id, !is.na(p.adjust), !is.na(qvalue)) %>%
			filter(p.adjust < fdr_cutoff, qvalue < fdr_cutoff) %>%
			group_by(view, direction) %>%
			arrange(p.adjust, .by_group = TRUE) %>%
			slice_head(n = max_terms_per_view) %>%
			ungroup()
	
	if (nrow(er) == 0) {
		return("No enrichment terms passed the selected FDR threshold. Interpret this factor primarily from its top features and variance-explained profile.")
	}
	
	pieces <- er %>%
			mutate(term_label = paste0(Description, " [", database, ", adj.P=", signif(p.adjust, 3), ", q=", signif(qvalue, 3), "]")) %>%
			group_by(view, direction) %>%
			summarise(terms = paste(term_label, collapse = "; "), .groups = "drop") %>%
			mutate(sentence = ifelse(
					direction == "ranked",
					paste0("The ranked loadings in the ", view,
							" view show enrichment for: ", terms, "."),
					paste0("The ", direction, " loadings in the ", view,
							" view are enriched for: ", terms, ".")
			)) %>%
			pull(sentence)
	
	paste(
			c(
					"Plausible interpretation based on enrichment results:",
					pieces,
					"This should be treated as a hypothesis-generating summary rather than a causal conclusion."
			),
			collapse = "\n\n"
	)
}

# -----------------------------------------------------------------------------
# Report rendering
# -----------------------------------------------------------------------------

#' Wrap text for the base-R PDF report.
#'
#' This helper prepares long interpretation strings for drawing into a PDF with
#' base R graphics. It does not require Pandoc, LaTeX, rmarkdown, or knitr.
#'
#' @param x Character string to wrap.
#' @param width Approximate maximum number of characters per line.
#'
#' @return Character vector of wrapped lines.
wrap_pdf_text <- function(x, width = 95) {
	# Prepare long text for direct drawing into a base R PDF device.
	# Replace missing text with an empty string so text wrapping is stable.
	x <- ifelse(is.na(x), "", x)
	
	# Preserve paragraph breaks by wrapping each newline-separated paragraph.
	unlist(lapply(strsplit(x, "
									", fixed = TRUE)[[1]], function(z) {
						if (nchar(z) == 0) return("")
						strwrap(z, width = width)
					}), use.names = FALSE)
}

#' Open a blank page in the base-R PDF report.
#'
#' @param title Optional page title.
#' @param subtitle Optional subtitle printed under the title.
#'
#' @return Invisibly returns NULL.
new_pdf_page <- function(title = NULL, subtitle = NULL) {
	# Create a blank plotting canvas for one PDF page.
	# Initialize a blank page with normalized x/y coordinates from 0 to 1.
	graphics::plot.new()
	graphics::par(mar = c(0, 0, 0, 0))
	graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1))
	
	# Draw title text near the top-left corner.
	if (!is.null(title)) {
		graphics::text(0.04, 0.96, title, adj = c(0, 1), cex = 1.3, font = 2)
	}
	
	# Draw subtitle text, when supplied.
	if (!is.null(subtitle)) {
		graphics::text(0.04, 0.915, subtitle, adj = c(0, 1), cex = 0.85)
	}
	
	invisible(NULL)
}

#' Draw wrapped prose on the current PDF page.
#'
#' @param text Character string to draw.
#' @param x Left x-coordinate in normalized page coordinates.
#' @param y Starting y-coordinate in normalized page coordinates.
#' @param width Approximate wrapping width in characters.
#' @param cex Text size.
#' @param line_height Vertical spacing between lines.
#'
#' @return The y-coordinate immediately below the drawn text.
draw_wrapped_text <- function(text,
		x = 0.04,
		y = 0.88,
		width = 95,
		cex = 0.72,
		line_height = 0.026) {
	# Draw text line-by-line so long summaries fit on the page.
	# Wrap prose before drawing it line by line.
	lines <- wrap_pdf_text(text, width = width)
	
	# Draw each line and decrement the y-coordinate after each line.
	for (line in lines) {
		graphics::text(x, y, line, adj = c(0, 1), cex = cex)
		y <- y - line_height
	}
	
	y
}

#' Convert a data frame to a compact monospaced text table.
#'
#' @param df Data frame to format.
#' @param max_rows Maximum number of rows to print.
#' @param max_col_width Default maximum number of characters per column value.
#' @param col_widths Optional named numeric vector giving per-column character
#'   widths. Values in `col_widths` override `max_col_width` for matching columns.
#'
#' @return Character vector representing the table.
format_pdf_table <- function(df, max_rows = 25, max_col_width = 28, col_widths = NULL) {
	# Convert tabular results into fixed-width text for PDF display.
	# Return a clear placeholder if the input table is empty.
	if (is.null(df) || nrow(df) == 0) return("No rows to display.")
	
	# Keep each rendered table chunk small enough to fit on the page.
	df <- as.data.frame(head(df, max_rows), stringsAsFactors = FALSE)
	
	# Convert values to character and truncate long entries. Named column-specific
	# widths allow important fields, such as gene-set descriptions, to be wider.
	for (col_name in names(df)) {
		this_width <- max_col_width
		if (!is.null(col_widths) && col_name %in% names(col_widths)) {
			this_width <- as.integer(col_widths[[col_name]])
		}
		this_width <- max(4L, this_width)
		col <- as.character(df[[col_name]])
		col[is.na(col)] <- ""
		df[[col_name]] <- ifelse(
				nchar(col) > this_width,
				paste0(substr(col, 1, this_width - 3), "..."),
				col
		)
	}
	
	# Capture the base R data-frame printout for drawing with graphics::text().
	capture.output(print(df, row.names = FALSE, right = FALSE))
}

#' Draw a data frame across one or more PDF pages.
#'
#' This table renderer intentionally uses base R graphics only. The result is
#' less polished than a LaTeX/rmarkdown table, but it is portable and avoids
#' external document-conversion dependencies.
#'
#' @param df Data frame to draw.
#' @param title Page title.
#' @param subtitle Optional page subtitle.
#' @param max_rows_per_page Maximum number of data rows per page.
#' @param max_col_width Default maximum number of characters per table cell.
#' @param col_widths Optional named numeric vector giving per-column character
#'   widths. Use this to widen specific columns in report tables.
#'
#' @return Invisibly returns NULL.
draw_table_pages <- function(df,
		title,
		subtitle = NULL,
		max_rows_per_page = 25,
		max_col_width = 28,
		col_widths = NULL) {
	# Paginate large tables across multiple PDF pages.
	# Draw a one-page placeholder for empty tables.
	if (is.null(df) || nrow(df) == 0) {
		new_pdf_page(title = title, subtitle = subtitle)
		draw_wrapped_text("No rows to display.", y = 0.88)
		return(invisible(NULL))
	}
	
	# Split rows into page-sized chunks.
	page_ids <- ceiling(seq_len(nrow(df)) / max_rows_per_page)
	chunks <- split(df, page_ids)
	
	# Render each chunk on a separate PDF page.
	for (i in seq_along(chunks)) {
		page_suffix <- if (length(chunks) > 1) {
					paste0(" page ", i, " of ", length(chunks))
				} else {
					NULL
				}
		
		new_pdf_page(title = title, subtitle = paste(c(subtitle, page_suffix), collapse = ""))
		
		table_lines <- format_pdf_table(
				chunks[[i]],
				max_rows = max_rows_per_page,
				max_col_width = max_col_width,
				col_widths = col_widths
		)
		
		# Draw the table in a monospaced font so columns remain aligned.
		y <- 0.88
		for (line in table_lines) {
			graphics::text(0.04, y, line, adj = c(0, 1), cex = 0.53, family = "mono")
			y <- y - 0.026
		}
	}
	
	invisible(NULL)
}

#' Render the MOFA interpretation report using only base R graphics.
#'
#' This is the default, no-Pandoc PDF engine. It writes a PDF directly through
#' grDevices::pdf(), so it does not require Pandoc, LaTeX, rmarkdown, or knitr.
#'
#' @param results List returned by run_mofa_integrative_interpretation().
#' @param output_dir Directory where the PDF should be written.
#' @param output_file Name of the PDF file.
#' @param report_significance_cutoff Significance threshold used for report prose and for filtering enrichment rows shown in the PDF report. Rows are shown only when both `p.adjust` and `qvalue` are below this threshold.
#' @param top_enrichment_rows Number of enrichment rows to print per group.
#' @param include_top_features Logical. If `TRUE`, include factor-specific top
#'   loaded feature tables in the PDF report.
#'
#' @return Path to the written PDF file.
render_mofa_pdf_report_base <- function(results,
		output_dir = "mofa_integrative_report",
		output_file = "mofa_integrative_factor_report.pdf",
		report_significance_cutoff = 0.10,
		top_enrichment_rows = 40,
		include_top_features = FALSE) {
	# Render the portable no-Pandoc report by drawing pages directly to PDF.
	# Create output directory and output path.
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
	pdf_path <- file.path(output_dir, output_file)
	
	# Open the PDF device. This is the key step that avoids Pandoc/LaTeX.
	grDevices::pdf(pdf_path, width = 11, height = 8.5, onefile = TRUE)
	
	# Ensure the PDF file is closed properly even if an error occurs.
	on.exit(grDevices::dev.off(), add = TRUE)
	
	# Unpack result components for readability.
	loadings <- results$loadings
	selected_features <- results$selected_features
	enrichment_results <- results$enrichment_results
	top_features <- results$top_features
	variance_explained <- results$variance_explained
	interpretations <- results$interpretations
	factor_ids <- sort(unique(loadings$factor))
	
	# Title page.
	new_pdf_page(
			title = "MOFA Integrative Factor Interpretation Report",
			subtitle = paste0("Generated on ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
	)
	y <- draw_wrapped_text(
			paste(
					"This report summarizes MOFA feature loadings and view-specific enrichment analyses.",
					"Ranked fgsea views are labeled as ranked; ORA-based metabolomics results retain positive/negative labels.",
					"Interpretations are hypothesis-generating and should be manually reviewed.",
					sep = " "
			),
			y = 0.84,
			cex = 0.82
	)
	draw_wrapped_text(
			paste0("Report significance threshold for summary text and enrichment tables: ", report_significance_cutoff),
			y = y - 0.04,
			cex = 0.78
	)
	
	# Overview: number of features per view.
	overview_features <- loadings %>%
			dplyr::count(view, name = "n_features") %>%
			arrange(view)
	draw_table_pages(overview_features, title = "Overview: features by view", max_rows_per_page = 35)
	
	# Omit the selected-feature overview table from the PDF report.
	# The selected features are still written to selected_factor_features.csv.
	
	# Overview: variance explained, if MOFA2 provided it successfully.
	if (!is.null(variance_explained) && nrow(variance_explained) > 0) {
		draw_table_pages(
				variance_explained %>% arrange(factor, view),
				title = "Overview: variance explained",
				max_rows_per_page = 35
		)
	}
	
	# Factor-specific sections.
	for (fac in factor_ids) {
		# Interpretation page.
		new_pdf_page(title = paste("Factor", fac), subtitle = "Summary interpretation")
		interp <- interpretations[[fac]]
		if (is.null(interp) || is.na(interp)) {
			interp <- "No interpretation text was generated for this factor."
		}
		draw_wrapped_text(interp, y = 0.88, cex = 0.72)
		
		# Factor-specific variance explained.
		if (!is.null(variance_explained) && nrow(variance_explained) > 0) {
			ve <- variance_explained %>%
					filter(factor == fac) %>%
					arrange(desc(variance_explained))
			draw_table_pages(ve, title = paste("Factor", fac, "variance explained"), max_rows_per_page = 30)
		}
		
		# Factor-specific top features are optional because they can make the PDF long.
		if (isTRUE(include_top_features)) {
			tf <- top_features %>%
					filter(factor == fac) %>%
					arrange(view, direction, desc(abs(value)))
			draw_table_pages(
					tf,
					title = paste("Factor", fac, "top loaded features"),
					max_rows_per_page = 24,
					max_col_width = 24
			)
		}
		
		# Factor-specific enrichment results. Only terms significant by both adjusted
		# p-value and q-value are included in the PDF report. The full, unfiltered
		# enrichment table remains available in enrichment_results.csv.
		er <- enrichment_results %>%
				filter(factor == fac) %>%
				filter(!is.na(p.adjust), !is.na(qvalue)) %>%
				filter(p.adjust < report_significance_cutoff, qvalue < report_significance_cutoff) %>%
				arrange(view, direction, p.adjust) %>%
				select(any_of(c(
										"view", "direction", "database", "ID", "Description", "GeneRatio",
										"BgRatio", "Count", "pvalue", "p.adjust", "qvalue", "geneID"
								))) %>%
				group_by(view, direction, database) %>%
				slice_head(n = top_enrichment_rows) %>%
				ungroup()
		
		draw_table_pages(
				er,
				title = paste("Factor", fac, "enrichment results"),
				max_rows_per_page = 20,
				max_col_width = 16,
				col_widths = c(Description = 55, database = 22, geneID = 26)
		)
	}
	
	invisible(pdf_path)
}

#' Write an optional R Markdown report template.
#'
#' This function is retained for users who choose report_engine =
#' "rmarkdown_pdf". That engine can produce a more polished PDF but requires
#' rmarkdown, Pandoc, and LaTeX.
#'
#' @param report_rmd_path Path where the .Rmd file should be written.
#'
#' @return Invisibly returns report_rmd_path.
write_report_rmd <- function(report_rmd_path) {
	# Build the optional R Markdown template as literal text.
	# Build the R Markdown report template as a character vector.
	rmd <- c(
			"---",
			"title: \"MOFA Integrative Factor Interpretation Report\"",
			"output:",
			"  pdf_document:",
			"    toc: true",
			"    toc_depth: 3",
			"    number_sections: true",
			"    latex_engine: xelatex",
			"params:",
			"  results_rds: null",
			"  significance_cutoff: 0.10",
			"  include_top_features: false",
			"---",
			"",
			"```{r setup, include=FALSE}",
			"knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)",
			"library(dplyr)",
			"library(tidyr)",
			"library(purrr)",
			"library(stringr)",
			"library(knitr)",
			"res <- readRDS(params$results_rds)",
			"loadings <- res$loadings",
			"selected_features <- res$selected_features",
			"enrichment_results <- res$enrichment_results",
			"top_features <- res$top_features",
			"variance_explained <- res$variance_explained",
			"interpretations <- res$interpretations",
			"factors <- sort(unique(loadings$factor))",
			"```",
			"",
			"# Overview",
			"",
			"This report summarizes feature loadings and enrichment analyses for each MOFA factor. Ranked fgsea views are reported using the label ranked; ORA-based metabolomics results retain positive/negative labels. Enrichment results are hypothesis-generating and should be manually reviewed.",
			"",
			"```{r overview-tables}",
			"loadings %>% count(view, name = 'n_features') %>% kable()",
			"```",
			"",
			"```{r variance-explained}",
			"if (nrow(variance_explained) > 0) {",
			"  cat('## Variance explained by factor and view\n\n')",
			"  variance_explained %>% arrange(factor, view) %>% kable(digits = 3)",
			"}",
			"```",
			"",
			"# Factor-level results",
			"",
			"```{r factor-sections, results='asis'}",
			"for (fac in factors) {",
			"  cat('\n\n## ', fac, '\n\n', sep = '')",
			"  if (fac %in% names(interpretations)) {",
			"    cat('### Summary interpretation\n\n')",
			"    cat(gsub('\n', '\n\n', interpretations[[fac]]), '\n\n')",
			"  }",
			"  if (nrow(variance_explained) > 0) {",
			"    ve <- variance_explained %>% filter(factor == fac) %>% arrange(desc(variance_explained))",
			"    if (nrow(ve) > 0) {",
			"      cat('### Variance explained\n\n')",
			"      print(kable(ve, digits = 3))",
			"      cat('\n\n')",
			"    }",
			"  }",
			"  if (isTRUE(params$include_top_features)) {",
			"    tf <- top_features %>% filter(factor == fac) %>% arrange(view, direction, desc(abs(value)))",
			"    if (nrow(tf) > 0) {",
			"      cat('### Top loaded features\n\n')",
			"      print(kable(tf, digits = 4))",
			"      cat('\n\n')",
			"    }",
			"  }",
			"  er <- enrichment_results %>% filter(factor == fac) %>% filter(!is.na(p.adjust), !is.na(qvalue)) %>% filter(p.adjust < params$significance_cutoff, qvalue < params$significance_cutoff) %>% arrange(view, direction, p.adjust)",
			"  if (nrow(er) > 0) {",
			"    cat('### Enrichment results\n\n')",
			"    er2 <- er %>%",
			"      select(any_of(c('view', 'direction', 'database', 'ID', 'Description', 'GeneRatio', 'BgRatio', 'Count', 'pvalue', 'p.adjust', 'qvalue', 'geneID'))) %>%",
			"      group_by(view, direction, database) %>%",
			"      slice_head(n = 10) %>%",
			"      ungroup()",
			"    print(kable(er2, digits = 4))",
			"    cat('\n\n')",
			"  } else {",
			"    cat('No enrichment results were produced for this factor.\n\n')",
			"  }",
			"}",
			"```"
	)
	
	# Write the report template to disk.
	writeLines(rmd, report_rmd_path)
	invisible(report_rmd_path)
}

#' Render the MOFA interpretation report through rmarkdown.
#'
#' This optional engine produces a more polished PDF, but requires rmarkdown,
#' Pandoc, and LaTeX. Use report_engine = "base_pdf" to avoid those external
#' dependencies.
#'
#' @param results List returned by run_mofa_integrative_interpretation().
#' @param output_dir Directory where output files should be written.
#' @param output_file Name of the PDF file.
#' @param report_significance_cutoff Significance threshold passed into the R Markdown report.
#' @param include_top_features Logical. If `TRUE`, include factor-specific top
#'   loaded feature tables in the PDF report.
#'
#' @return Path to the written PDF file.
render_mofa_report_rmarkdown_pdf <- function(results,
		output_dir = "mofa_integrative_report",
		output_file = "mofa_integrative_factor_report.pdf",
		report_significance_cutoff = 0.10,
		include_top_features = FALSE) {
	# Render the optional polished PDF through rmarkdown.
	# rmarkdown requires Pandoc, and PDF rendering also requires LaTeX.
	require_namespace("rmarkdown", "rmarkdown PDF report rendering")
	
	# Create the output directory.
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
	
	# Save the analysis results for use by the R Markdown template.
	results_rds <- file.path(output_dir, "mofa_integrative_results.rds")
	saveRDS(results, results_rds)
	
	# Write the parameterized R Markdown template.
	report_rmd <- file.path(output_dir, "mofa_integrative_report.Rmd")
	write_report_rmd(report_rmd)
	
	# Render the final PDF through rmarkdown.
	rmarkdown::render(
			input = report_rmd,
			output_file = output_file,
			output_dir = output_dir,
			params = list(
					results_rds = normalizePath(results_rds, winslash = "/", mustWork = TRUE),
					significance_cutoff = report_significance_cutoff,
				include_top_features = include_top_features
			),
			envir = new.env(parent = globalenv()),
			quiet = TRUE
	)
	
	file.path(output_dir, output_file)
}

# -----------------------------------------------------------------------------
# End-to-end wrapper
# -----------------------------------------------------------------------------

#' Run an end-to-end integrative interpretation of MOFA factors.
#'
#' This is the main user-facing function for the script. It extracts feature
#' loadings from a fitted MOFA/MOFA2 object, selects high-loading features for
#' each factor, view, and loading direction, performs view-specific enrichment
#' analyses, uses one full signed-vector fgsea run per factor/view for ranked
#' expression, proteomics, and methylation enrichment, writes machine-readable
#' result tables, and optionally renders a PDF report.
#'
#' The function assumes that MOFA views are named consistently as any subset of
#' `"proteomics"`, `"methylation"`, `"metabolomics"`, and `"expression"`.
#' Proteomics features are assumed to be UniProt IDs, expression features are
#' assumed to be Entrez IDs or gene symbols, methylation features are assumed to
#' be Illumina CpG IDs that can be mapped to Entrez IDs by missMethyl, and
#' metabolomics features are matched to a Metabolon annotation table.
#'
#' @param mofa_obj A fitted MOFA/MOFA2 object.
#' @param metabolon_annotation Optional data frame containing Metabolon
#'   annotations. Required for metabolomics enrichment. The table should include
#'   the metabolite feature column specified by `metabolomics_feature_col` and
#'   pathway/category columns such as `SUPER_PATHWAY` and `SUB_PATHWAY`.
#' @param output_dir Character string giving the directory where result files
#'   and the report should be written.
#' @param views Character vector of MOFA view names to analyze. Unsupported view
#'   names are ignored with a warning.
#' @param factors Factors to extract from the MOFA object. Passed to
#'   `MOFA2::get_weights()`; use `"all"` to analyze all factors.
#' @param selection_mode Character string specifying how high-loading features
#'   are selected. Options are `"top_n"`, `"zscore"`, and `"quantile"`.
#' @param top_n Integer. Number of positive and negative features selected per
#'   factor and view when `selection_mode = "top_n"`.
#' @param z_cutoff Numeric. Absolute z-score threshold used when
#'   `selection_mode = "zscore"`.
#' @param quantile_cutoff Numeric between 0 and 1. Absolute-loading quantile
#'   used when `selection_mode = "quantile"`.
#' @param min_features Integer. Minimum number of selected features required for
#'   a factor/view/direction group to be retained for enrichment.
#' @param top_feature_n Integer. Number of top-loaded features to include in the
#'   top-feature output table. These features are shown in the PDF only when
#'   `include_top_features_in_report = TRUE`.
#' @param metabolomics_feature_col Character string naming the column in
#'   `metabolon_annotation` that matches the metabolomics feature names in the
#'   MOFA object. Common choices are `"BIOCHEMICAL"` and `"chemical_name"`.
#' @param expression_id_type Character string specifying expression feature ID
#'   type. Options are `"auto"`, `"ENTREZID"`, and `"SYMBOL"`.
#' @param OrgDb Character string naming the Bioconductor organism annotation
#'   package used for gene/protein ID mapping, for example `"org.Hs.eg.db"`.
#' @param organism Character string specifying the organism name used by
#'   `ReactomePA::enrichPathway()`, for example `"human"`.
#' @param go_ontology Character string specifying the GO ontology used to build
#'   ranked gene sets for `fgsea`, commonly `"BP"`, `"MF"`, or `"CC"`.
#' @param methylation_array_type Character string specifying the Illumina array
#'   type used by missMethyl to map CpGs to Entrez IDs. Supported values in this
#'   script are `"EPIC"` and `"450K"`.
#' @param methylation_collection Retained for compatibility with older ORA-based
#'   methylation code paths; ranked methylation enrichment uses `go_ontology`.
#' @param qvalue_cutoff Numeric adjusted p-value cutoff used for reported
#'   enrichment results, including fgsea results for expression, proteomics, and
#'   methylation views.
#' @param report_significance_cutoff Numeric significance threshold used when generating the
#'   factor-level interpretation text and filtering enrichment rows in the PDF
#'   report. Rows are shown only when both `p.adjust` and `qvalue` are below
#'   this threshold.
#' @param metabolon_term_cols Character vector naming Metabolon annotation
#'   columns used for metabolomics categorical enrichment.
#' @param render_report Logical. If `TRUE`, render a PDF report after running
#'   the enrichment analyses.
#' @param include_top_features_in_report Logical. If `TRUE`, include factor-specific
#'   top loaded feature tables in the PDF report. Defaults to `FALSE` to keep
#'   the report focused on enrichment summaries.
#' @param save_loadings Logical. If `TRUE`, save information about features and 
#'   feature loadings in a number of CSV files. 
#' @param report_engine Character string specifying the report engine. Use
#'   `"base_pdf"` to generate a PDF directly with base R graphics without
#'   Pandoc or LaTeX. Use `"rmarkdown_pdf"` for a more polished PDF that
#'   requires rmarkdown, Pandoc, and LaTeX.
#' @param n_cores Number of task-level parallel workers for enrichment analyses.
#'   The number of enrichment analyses that will performed equal (1) the number
#'   of factors specified in argument `factors`times (2) the number of data 
#'   views in the MOFA results object passed as the value of `mofa_obj`. Each
#'   of the `n_cores` specified will be dispatched to handle a single enrichment
#'   analysis. Use `n_cores = 1` for serial execution. On Windows, parallel task 
#'   execution requires the optional `future` and `future.apply` packages.
#' @param fgsea_nproc Number of CPUs passed to `fgsea::fgsea()` inside each
#'   expression, proteomics, or methylation enrichment task. Keep this at `1`
#'   when `n_cores > 1` to avoid oversubscribing CPUs. 
#'
#' @return Invisibly returns a named list containing `loadings`,
#'   `variance_explained`, `selected_features`, `top_features`,
#'   `enrichment_results`, `interpretations`, `parameters`, and, when
#'   `render_report = TRUE`, `report_path`. The same core result tables are also
#'   written to CSV files in `output_dir`.
#'
#' @examples
#' # mofa_obj <- readRDS("path/to/mofa_results_object.rds")
#' # metabolon_annot <- readr::read_tsv("path/to/metabolon_annotation.tsv")
#' # res <- run_mofa_integrative_interpretation(
#' #   mofa_obj = mofa_obj,
#' #   metabolon_annotation = metabolon_annot,
#' #   output_dir = "mofa_factor_interpretation_report",
#' #   metabolomics_feature_col = "BIOCHEMICAL",
#' #   expression_id_type = "auto",
#' #   methylation_array_type = "EPIC",
#' #   render_report = TRUE,
#' #   include_top_features_in_report = FALSE,
#' #   report_engine = "base_pdf",
#' #   n_cores = 4,
#' #   fgsea_nproc = 1
#' # )
run_mofa_integrative_interpretation <- function(mofa_obj,
		metabolon_annotation = NULL,
		output_dir = "mofa_integrative_report",
		views = views_names(mofa_obj),
		factors = "all",
		selection_mode = c("top_n", "zscore", "quantile"),
		top_n = 100,
		z_cutoff = 2.5,
		quantile_cutoff = 0.95,
		min_features = 5,
		top_feature_n = 20,
		metabolomics_feature_col = "BIOCHEMICAL",
		expression_id_type = c("auto", "ENTREZID", "SYMBOL"),
		OrgDb = "org.Hs.eg.db",
		organism = "human",
		go_ontology = "BP",
		methylation_array_type = c("EPIC", "450K"),
		methylation_collection = "GO",
		qvalue_cutoff = 0.25,
		report_significance_cutoff = 0.05,
		metabolon_term_cols = c("SUPER_PATHWAY", "SUB_PATHWAY", "PARENT_CLASSIFICATION", "CLASSIFICATION"),
		render_report = TRUE,
		include_top_features_in_report = FALSE,
		save_loadings = FALSE,
		report_engine = c("base_pdf", "rmarkdown_pdf"),
		n_cores = 1L,
		fgsea_nproc = 1L) {
	# Validate user options before running the multi-step pipeline.
	selection_mode <- match.arg(selection_mode)
	expression_id_type <- match.arg(expression_id_type)
	methylation_array_type <- match.arg(methylation_array_type)
	
	# Select the PDF rendering backend. The default base_pdf engine writes the
	# report directly with grDevices::pdf() and therefore does not need Pandoc.
	report_engine <- match.arg(report_engine)
	
	# Normalize parallelization parameters once at the top-level wrapper.
	n_cores <- normalize_n_cores(n_cores)
	fgsea_nproc <- normalize_n_cores(fgsea_nproc)
	
	message("Extracting MOFA loadings...")
	loadings <- extract_mofa_loadings(
			mofa_obj = mofa_obj,
			views = views,
			factors = factors
	)
	
	message("Extracting variance explained, if available...")
	variance_explained <- extract_variance_explained(mofa_obj)
	
	message("Selecting high-loading features...")
	selected_features <- select_factor_features(
			loadings = loadings,
			selection_mode = selection_mode,
			top_n = top_n,
			z_cutoff = z_cutoff,
			quantile_cutoff = quantile_cutoff,
			min_features = min_features
	)
	
	message("Creating top-feature tables...")
	top_features <- make_top_feature_table(
			loadings = loadings,
			n = top_feature_n,
			metabolon_annotation = metabolon_annotation,
			metabolomics_feature_col = metabolomics_feature_col
	)
	
	message("Running view-specific enrichment analyses...")
	enrichment_results <- run_all_enrichments(
			loadings = loadings,
			selected_features = selected_features,
			metabolon_annotation = metabolon_annotation,
			metabolomics_feature_col = metabolomics_feature_col,
			expression_id_type = expression_id_type,
			OrgDb = OrgDb,
			organism = organism,
			go_ontology = go_ontology,
			methylation_array_type = methylation_array_type,
			methylation_collection = methylation_collection,
			qvalue_cutoff = qvalue_cutoff,
			metabolon_term_cols = metabolon_term_cols,
			n_cores = n_cores,
			fgsea_nproc = fgsea_nproc
	)
	
	message("Generating factor-level interpretation text...")
	factor_ids <- sort(unique(loadings$factor))
	interpretations <- setNames(
			map(factor_ids, ~ summarise_factor_interpretation(
							enrichment_results = enrichment_results,
							factor_id = .x,
							fdr_cutoff = report_significance_cutoff
					)),
			factor_ids
	)
	
	results <- list(
			loadings = loadings,
			variance_explained = variance_explained,
			selected_features = selected_features,
			top_features = top_features,
			enrichment_results = enrichment_results,
			interpretations = interpretations,
			parameters = list(
					selection_mode = selection_mode,
					top_n = top_n,
					z_cutoff = z_cutoff,
					quantile_cutoff = quantile_cutoff,
					min_features = min_features,
					metabolomics_feature_col = metabolomics_feature_col,
					expression_id_type = expression_id_type,
					OrgDb = OrgDb,
					organism = organism,
					go_ontology = go_ontology,
					methylation_array_type = methylation_array_type,
					methylation_collection = methylation_collection,
					qvalue_cutoff = qvalue_cutoff,
					report_significance_cutoff = report_significance_cutoff,
					include_top_features_in_report = include_top_features_in_report,
					n_cores = n_cores,
					fgsea_nproc = fgsea_nproc
			)
	)
	
	dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
	saveRDS(results, file.path(output_dir, 
					paste0("mofa_integrative_results_", go_ontology, ".rds")))
	write_csv(enrichment_results, file.path(output_dir, 
					paste0("enrichment_results_", go_ontology, ".csv")))
	if(save_loadings){
		write_csv(loadings, file.path(output_dir, "all_mofa_loadings.csv"))
		write_csv(selected_features, file.path(output_dir, "selected_factor_features.csv"))
		write_csv(top_features, file.path(output_dir, "top_factor_features.csv"))
	}
	
	if (render_report) {
		# Render the final factor-by-factor PDF report using the selected engine.
		if (report_engine == "base_pdf") {
			# Portable option: no Pandoc, no LaTeX, and no rmarkdown dependency.
			message("Rendering PDF report with base R graphics...")
			report_path <- render_mofa_pdf_report_base(
					results = results,
					output_dir = output_dir,
					output_file = paste0("mofa_integrative_factor_report_", go_ontology, ".pdf"),
					report_significance_cutoff = report_significance_cutoff,
					include_top_features = include_top_features_in_report
			)
		} else {
			# Polished option: requires rmarkdown, Pandoc, and LaTeX.
			message("Rendering PDF report with rmarkdown...")
			report_path <- render_mofa_report_rmarkdown_pdf(
					results = results,
					output_dir = output_dir,
					output_file = paste0("mofa_integrative_factor_report_", go_ontology, ".pdf"),
					report_significance_cutoff = report_significance_cutoff,
					include_top_features = include_top_features_in_report
			)
		}
		
		results$report_path <- report_path
		message("Report written to: ", report_path)
	}
	
	invisible(results)
}

# -----------------------------------------------------------------------------
# Run pipeine
# -----------------------------------------------------------------------------

library(MOFA2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(fgsea)
library(missMethyl)
library(rmarkdown)

# Specify file locations - edit according to your local execution environment
mofa_dir <- file.path("D:/R_sessions/PROSPR/MOFA",
		"2026_Jun_28_00_49_09_visit=24mo_change_scores_ungrouped_views=me")
#mofa_dir <- file.path("D:/R_sessions/PROSPR/MOFA",
#		"2026_May_15_00_40_17_visit=24mo_change_scores_ungrouped_views=me+pr+mb")
mofa_res_fn <- file.path(mofa_dir, "model.hdf5")
metab_ann_fn <- file.path("E:/Data/PROSPR-FAST/CALERIE/metabolomics",
		"metabolon_annotations.csv")
output_dir <- mofa_dir # save results at same location of MOFA results

mofa_obj <- MOFA2::load_model(mofa_res_fn)
metabolon_annot <- readr::read_csv(metab_ann_fn, show_col_types = FALSE)

res <- run_mofa_integrative_interpretation(
		mofa_obj = mofa_obj,
		metabolon_annotation = metabolon_annot,
		output_dir = file.path(mofa_dir, "mofa_factor_interpretation"),
		views = views_names(mofa_obj),
		factors = "all",
		selection_mode = "top_n",
		top_n = 100,
		min_features = 5,
		metabolomics_feature_col = "BIOCHEMICAL",
		expression_id_type = "auto",
		OrgDb = "org.Hs.eg.db",
		organism = "human",
		go_ontology = "BP",
		methylation_array_type = "EPIC",
		methylation_collection = "GO",
		qvalue_cutoff = 0.1,
		report_significance_cutoff = 0.05,
		render_report = TRUE,
		report_engine = "base_pdf",
		n_cores = 1, # Change according to your execution environment
		fgsea_nproc = 1 # Best to leave this at 1
)

# Important output files:
#   mofa_factor_interpretation_report/mofa_integrative_factor_report.pdf
#   mofa_factor_interpretation_report/enrichment_results.csv
#   
# To generate the files below, set the argument `save_loading` to TRUE when 
# calling run_mofa_integrative_interpretation().
#   mofa_factor_interpretation_report/top_factor_features.csv
#   mofa_factor_interpretation_report/selected_factor_features.csv
