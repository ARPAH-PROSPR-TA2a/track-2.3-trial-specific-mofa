# =============================================================================
# package_results.R
# =============================================================================
# Archive results using gzip-compressed POSIX pax tar where supported.
# Usage (RGui or Rscript):
#   source("path_to_code/package_results.R", local = TRUE, chdir = TRUE)
# Settings are defined in package_results_parameters.R.
# No additional R packages are required.
# =============================================================================

#' Package a results directory, respecting glob exclusions
#'
#' @param root_dir Directory whose *contents* should be archived.
#' @param exclude_files Character vector of shell-style glob patterns. Patterns
#'   without '/' match basenames anywhere; patterns containing '/' match paths
#'   relative to root_dir. A matching directory excludes its whole subtree.
#'   '*' and '?' are wildcards; matching is case-sensitive.
#' @param gzip_res_file Destination .tar.gz file, or NULL for a .tar.gz file
#'   alongside root_dir named after the directory.
#' @return Invisibly, the absolute destination path.
#' @details Entries are explicitly enumerated before archiving. Only files and
#'   empty directories are passed to tar, preventing tar from recursively
#'   reintroducing excluded files through a parent directory. Symlinks are
#'   skipped. Prefers GNU tar's POSIX pax format for long path support;
#'   otherwise uses R's internal tar writer (which may warn for long paths).
package_results <- function(root_dir, exclude_files = character(),
                            gzip_res_file = NULL) {
  if (!is.character(root_dir) || length(root_dir) != 1L || is.na(root_dir) ||
      !nzchar(root_dir) || !dir.exists(root_dir)) {
    stop("'root_dir' must be an existing directory.")
  }
  if (!is.character(exclude_files) || anyNA(exclude_files) ||
      any(!nzchar(exclude_files))) {
    stop("'exclude_files' must be a character vector of nonempty glob patterns.")
  }
  root_dir <- normalizePath(root_dir, winslash = "/", mustWork = TRUE)
  if (is.null(gzip_res_file)) {
    gzip_res_file <- file.path(dirname(root_dir),
                               paste0(basename(root_dir), ".tar.gz"))
  }
  if (!is.character(gzip_res_file) || length(gzip_res_file) != 1L ||
      is.na(gzip_res_file) || !nzchar(gzip_res_file)) {
    stop("'gzip_res_file' must be NULL or one nonempty filename.")
  }
  if (!dir.exists(dirname(gzip_res_file))) {
    stop("Output directory does not exist: ", dirname(gzip_res_file))
  }
  output <- file.path(normalizePath(dirname(gzip_res_file), winslash = "/"),
                      basename(gzip_res_file))
  if (!grepl("\\.tar\\.gz$", output, ignore.case = TRUE)) {
    stop("Output filename must end in '.tar.gz'.")
  }

  # glob2rx() anchors the pattern. Match basenames for simple patterns and
  # root-relative paths for patterns containing a slash.
  exclude_files <- gsub("\\\\", "/", exclude_files)
  regex <- vapply(exclude_files, utils::glob2rx, character(1), USE.NAMES = FALSE)
  path_pattern <- grepl("/", exclude_files, fixed = TRUE)
  matches_exclusion <- function(path) {
    any(vapply(seq_along(regex), function(i) {
      candidate <- if (path_pattern[i]) path else basename(path)
      grepl(regex[i], candidate)
    }, logical(1)))
  }

  # Include hidden files. List directories separately to identify empty ones.
  paths <- list.files(root_dir, recursive = TRUE, all.files = TRUE,
                      include.dirs = TRUE, no.. = TRUE, full.names = FALSE)
  paths <- gsub("\\\\", "/", paths)
  paths <- unique(paths)
  if (!length(paths)) stop("The results directory is empty.")
  absolute <- file.path(root_dir, paths)
  is_dir <- dir.exists(absolute)
  is_link <- nzchar(Sys.readlink(absolute))

  # Exclude a matching directory's descendants, including when only the
  # directory name (e.g. 'logs') is matched by a glob.
  direct_excluded <- vapply(paths, matches_exclusion, logical(1))
  excluded_dirs <- paths[is_dir & (direct_excluded | is_link)]
  excluded <- direct_excluded | is_link | (absolute == output)
  if (length(excluded_dirs)) {
    for (d in excluded_dirs) {
      excluded <- excluded | startsWith(paths, paste0(d, "/"))
    }
  }
  kept <- paths[!excluded]
  kept_is_dir <- is_dir[!excluded]
  if (!length(kept)) stop("No entries remain after applying exclusions.")

  # IMPORTANT: tar recursively descends into directories by default. Passing
  # nonempty directories would silently add excluded files back to the archive.
  # Therefore archive regular files and only genuinely empty directories.
  # Nonempty directory structure is reconstructed from each file's path.
  keep_dir <- kept_is_dir & !vapply(kept, function(p) {
    any(startsWith(kept, paste0(p, "/")))
  }, logical(1))
  entries <- kept[!kept_is_dir | keep_dir]
  if (!length(entries)) stop("No files or empty directories remain to archive.")

  # Use a temporary archive beside the final output. Preserve the existing
  # destination until the new archive has been successfully verified.
  temporary <- tempfile(".package_results_", tmpdir = dirname(output),
                        fileext = ".tar.gz")
  on.exit(unlink(temporary), add = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(root_dir)
  message("Packaging ", sum(!kept_is_dir), " files and ", sum(keep_dir),
          " empty directories (excluded ", sum(excluded), " entries).")

  # Prefer GNU tar for pax extended headers, which handle paths >100 bytes.
  # -T passes a list of entries; --no-recursion prevents directory expansion.
  tar_cmd <- Sys.which("tar")
  version <- if (nzchar(tar_cmd)) {
    suppressWarnings(tryCatch(system2(tar_cmd, "--version", stdout = TRUE,
                                      stderr = FALSE), error = function(e) ""))
  } else ""
  is_gnu <- length(version) && any(grepl("GNU tar", version, fixed = TRUE))
  if (is_gnu) {
    manifest <- tempfile("package_results_manifest_", fileext = ".txt")
    on.exit(unlink(manifest), add = TRUE)
    # Newline-delimited manifests cannot represent filenames with newlines.
    if (any(grepl("[\r\n]", entries))) {
      stop("Cannot package filenames containing newline characters.")
    }
    writeLines(paste0("./", entries), manifest, useBytes = TRUE)
    args <- c("--format=pax", "--no-recursion", "-czf",
              shQuote(temporary), "-T", shQuote(manifest))
    status <- system2(tar_cmd, args = args)
    if (!identical(as.integer(status), 0L)) stop("GNU tar failed: ", status)
  } else {
    message("GNU tar not available; using R internal tar. Long-path warnings ",
            "may occur on this system.")
    utils::tar(tarfile = temporary, files = entries,
               compression = "gzip", tar = "internal")
  }
  if (!file.exists(temporary) || is.na(file.info(temporary)$size) ||
      file.info(temporary)$size == 0) stop("Archive was not created successfully.")

  # Verify actual archive members, not just the intended input list.
  archived <- utils::untar(temporary, list = TRUE)
  normalize_member <- function(x) {
    x <- gsub("\\\\", "/", x)
    x <- sub("^\\./", "", x)
    sub("/+$", "", x)
  }
  actual <- normalize_member(archived)
  expected <- normalize_member(entries)
  if (!setequal(actual, expected)) {
    stop("Archive verification failed. Missing: ",
         paste(head(setdiff(expected, actual), 5L), collapse = ", "),
         "; unexpected: ",
         paste(head(setdiff(actual, expected), 5L), collapse = ", "))
  }
  if (file.exists(output) && !file.remove(output)) {
    stop("Cannot replace existing archive: ", output)
  }
  if (!file.rename(temporary, output) &&
      !file.copy(temporary, output, overwrite = TRUE)) {
    stop("Cannot move archive to destination: ", output)
  }
  message("Verified archive: ", output, " (",
          format(file.info(output)$size, big.mark = ","), " bytes).")
  invisible(output)
}

# chdir = TRUE makes this path relative to the script, not the caller.
source("package_results_parameters.R", local = TRUE)
package_results(root_dir, exclude_files, gzip_res_file)
