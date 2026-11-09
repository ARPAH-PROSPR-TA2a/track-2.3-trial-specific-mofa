# =============================================================================
# Parameters for package_results.R
# =============================================================================

# Directory whose contents will be packaged recursively.
root_dir <- "D:/R_sessions/PROSPR/MOFA/CALERIE"

# Output archive. NULL creates <basename(root_dir)>.tar.gz next to root_dir.
gzip_res_file <- NULL


# *****************************************************************************
#
# The values of the following parameters should be left to their default settings
#
# *****************************************************************************

# Glob patterns matched against each file's basename AND its path relative
# to root_dir (forward-slash separators). Matching is case-sensitive.
# Examples:
#   "*.hdf5"       excludes all HDF5 files at any depth
#   "logs"         excludes a directory named logs and its contents
#   "logs/*"       excludes files and folders under the top-level logs folder
#   "plots/*.png"  excludes PNG files directly under plots (and may match
#                  deeper paths because '*' also matches '/' in glob2rx)
# Set character(0) to include everything.
exclude_files <- c("*.hdf5", "sample_metadata*")
