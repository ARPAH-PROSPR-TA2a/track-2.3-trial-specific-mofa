# Check packages used by R scripts and R Markdown documents in this directory.
# Run from RGui with:
#   source("path_to_code/check_environment.R", local = TRUE, chdir = TRUE)
#
# Files are parsed, not executed. .R and .Rmd files in this directory are
# inspected, including R code chunks in .Rmd files; subdirectories, inline R
# expressions, chunk options, and YAML headers are not inspected. Package
# names must appear literally in library(), require(), requireNamespace(),
# loadNamespace(), ::, or ::: calls. Computed names cannot be determined
# reliably without executing the code and are reported for manual review.
# Conditional calls are included even if their condition would be false.

local({
  r_files <- list.files(".", pattern = "\\.(r|rmd)$", full.names = TRUE, ignore.case = TRUE)
  if (length(r_files) == 0L) {
    cat("No R or R Markdown files found in:", normalizePath("."), "\n")
    return(invisible(NULL))
  }

  packages <- character()
  unresolved <- character()
  parse_failures <- character()

  # Extract a literal package name. Unquoted symbols are accepted only where
  # R itself treats them as package names (e.g., library(dplyr)).
  literal_name <- function(x, allow_symbol = FALSE) {
    if (missing(x)) return(NULL)
    if (is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)) {
      return(x)
    }
    if (allow_symbol && is.symbol(x)) return(as.character(x))
    NULL
  }

  inspect <- function(node, file) {
    # Calls such as x[, 1] contain an omitted argument. Passing it through
    # the recursive traversal creates a missing argument, which must be skipped.
    if (missing(node)) return(invisible(NULL))
    if (is.expression(node) || is.pairlist(node)) {
      for (child in node) inspect(child, file)
      return(invisible(NULL))
    }
    if (!is.call(node)) return(invisible(NULL))

    head <- node[[1L]]
    fun <- if (is.symbol(head)) as.character(head) else ""
    # Recognize e.g. base::library(...) as well as library(...).
    if (is.call(head) && length(head) == 3L &&
        is.symbol(head[[1L]]) && as.character(head[[1L]]) %in% c("::", ":::")) {
      fun <- as.character(head[[3L]])
    }

    if (fun %in% c("::", ":::")) {
      name <- literal_name(node[[2L]], allow_symbol = TRUE)
      if (!is.null(name)) packages <<- c(packages, name)
    } else if (fun %in% c("library", "require", "requireNamespace", "loadNamespace")) {
      args <- as.list(node)[-1L]
      arg_names <- names(args)
      if (is.null(arg_names)) arg_names <- rep("", length(args))
      arg_names[is.na(arg_names)] <- ""
      # 'package' is the first formal argument in these functions. If it is
      # explicitly named, use it; otherwise use the first unnamed argument.
      idx <- which(arg_names == "package")
      if (!length(idx)) idx <- which(arg_names == "")[1L]
      char_idx <- which(arg_names == "character.only")
      char_only_false <- !length(char_idx) || identical(args[char_idx[1L]], list(FALSE))
      allow_symbol <- fun %in% c("library", "require") && char_only_false
      name <- if (length(idx) && !is.na(idx)) literal_name(args[[idx[1L]]], allow_symbol) else NULL
      if (is.null(name)) {
        unresolved <<- c(unresolved, paste0(basename(file), ": ", paste(deparse(node, width.cutoff = 120L), collapse = " ")))
      } else {
        packages <<- c(packages, name)
      }
    }
    # Traverse the call head too, so namespace-qualified calls are detected.
    for (child in as.list(node)) inspect(child, file)
    invisible(NULL)
  }

  # Parse each code block independently, so an error in one R Markdown chunk
  # does not prevent checking the others. Report its starting line for review.
  inspect_code <- function(code, label) {
    parsed <- tryCatch(parse(text = code, keep.source = FALSE), error = function(e) {
      parse_failures <<- c(parse_failures, paste0(label, ": ", conditionMessage(e)))
      NULL
    })
    if (!is.null(parsed)) inspect(parsed, label)
  }

  # Identify fenced Markdown blocks and inspect only those with an R engine.
  # Both backtick and tilde fences are supported, including longer fences.
  inspect_rmd <- function(file) {
    lines <- readLines(file, warn = FALSE)
    in_block <- FALSE
    r_block <- FALSE
    fence_char <- ""
    fence_length <- 0L
    code <- character()
    start_line <- 0L

    for (i in seq_along(lines)) {
      line <- lines[i]
      if (!in_block) {
        opening <- regmatches(line, regexec("^ {0,3}(`{3,}|~{3,})(.*)$", line))[[1L]]
        if (length(opening)) {
          in_block <- TRUE
          fence <- opening[2L]
          fence_char <- substr(fence, 1L, 1L)
          fence_length <- nchar(fence)
          r_block <- grepl("^[[:space:]]*\\{[[:space:]]*[rR]([[:space:],}]|$)", opening[3L])
          code <- character()
          start_line <- i
        }
      } else {
        closing_pattern <- paste0("^ {0,3}", fence_char, "{", fence_length,
                                  ",}[[:space:]]*$")
        if (grepl(closing_pattern, line)) {
          if (r_block) inspect_code(code, paste0(file, " (chunk starting at line ", start_line, ")"))
          in_block <- FALSE
        } else if (r_block) {
          code <- c(code, line)
        }
      }
    }
    if (in_block && r_block) {
      parse_failures <<- c(parse_failures, paste0(basename(file),
        ": R chunk starting at line ", start_line, " has no closing fence."))
    }
  }

  for (file in r_files) {
    tryCatch({
      if (grepl("\\.rmd$", file, ignore.case = TRUE)) {
        inspect_rmd(file)
      } else {
        inspect_code(readLines(file, warn = FALSE), file)
      }
    }, error = function(e) {
      parse_failures <<- c(parse_failures, paste0(basename(file), ": ", conditionMessage(e)))
    })
  }

  packages <- sort(unique(packages))
  # installed.packages() searches the library paths available to this R
  # session, including the user's library, rather than a different R install.
  installed <- rownames(utils::installed.packages(lib.loc = .libPaths()))
  missing <- setdiff(packages, installed)

  cat("Checked", length(r_files), "R or R Markdown file(s) in", normalizePath("."), "\n")
  if (length(missing)) {
    cat("\nMissing packages:", paste(missing, collapse = ", "), "\n")
    cat("To install them, run:\n")
    cat("install.packages(c(", paste(encodeString(missing, quote = '"'), collapse = ", "), "))\n", sep = "")
    cat("Note: packages unavailable from your configured repositories may need a different installation method.\n")
  } else if (length(packages)) {
    cat("\nAll", length(packages), "identified packages are installed in the current R session's library paths.\n")
  } else {
    cat("\nNo packages with literal names were identified.\n")
  }

  if (length(unresolved)) {
    cat("\nPackage calls requiring manual review (names could not be determined):\n")
    cat(paste0("  ", unique(unresolved), collapse = "\n"), "\n", sep = "")
  }
  if (length(parse_failures)) {
    cat("\nThe following files or code chunks could not be checked:\n")
    cat(paste0("  ", parse_failures, collapse = "\n"), "\n", sep = "")
  }
  if (length(unresolved) || length(parse_failures)) {
    cat("\nThe check is incomplete; review the items above.\n")
  }
  invisible(NULL)
})
