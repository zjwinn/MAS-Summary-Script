#!/usr/bin/env Rscript

# Marker-Assisted Selection (MAS) Calling Script (CLI Runner)
# Author: Zachary Winn
#
# Usage:
#   Rscript MAS-Calling-Script.R --vcf <input.vcf.gz> --key <key.txt> --outdir <output_dir> [options]

# ---------------------------------------------------------
# 1. Environment & Package Setup
# ---------------------------------------------------------

# Add user library path for Linux/HPC/WSL environments if present
user_r_lib <- path.expand("~/R/library")
if (dir.exists(user_r_lib) && !(user_r_lib %in% .libPaths())) {
  .libPaths(c(user_r_lib, .libPaths()))
}

# Define required packages
required_packages <- c("parallel", "foreach", "doParallel", "data.table", "R.utils")

# Check and install missing packages
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# Load packages
suppressPackageStartupMessages({
  library(parallel)
  library(foreach)
  library(doParallel)
  library(data.table)
})

# Determine script directory
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), mustWork = FALSE)))
  }
  sys_frames <- sys.frames()
  for (f in sys_frames) {
    if (!is.null(f$ofile)) {
      return(dirname(normalizePath(f$ofile, mustWork = FALSE)))
    }
  }
  return(getwd())
}

script_dir <- get_script_dir()

# Source modular helper functions from the R/ directory
r_helper_dir <- file.path(script_dir, "R")
if (!dir.exists(r_helper_dir)) {
  r_helper_dir <- file.path(getwd(), "R")
}
if (!dir.exists(r_helper_dir)) {
  stop("Helper function directory not found at: ", r_helper_dir)
}
r_helper_files <- list.files(r_helper_dir, pattern = "\\.R$", full.names = TRUE)
invisible(lapply(r_helper_files, source))

# ---------------------------------------------------------
# 2. Command-Line Argument Parsing
# ---------------------------------------------------------

parse_cli_args <- function(args) {
  parsed <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg %in% c("-h", "--help")) {
      cat("Marker-Assisted Selection (MAS) Summary Script\n")
      cat("Usage: Rscript MAS-Calling-Script.R [options]\n\n")
      cat("Options:\n")
      cat("  -v, --vcf FILE        Input VCF file (.vcf or .vcf.gz) [required for pipeline]\n")
      cat("  -k, --key FILE        Haplotype key file (.txt tab-delimited) [required for pipeline]\n")
      cat("  -o, --outdir DIR      Output directory for results [default: example/]\n")
      cat("  -c, --cores INT       Number of CPU cores for parallel processing [default: detectCores()-1]\n")
      cat("  --verbose [BOOL]      Enable verbose logging (TRUE/FALSE) [default: TRUE]\n")
      cat("  -h, --help            Show this help message and exit\n\n")
      cat("Example:\n")
      cat("  Rscript MAS-Calling-Script.R --vcf sample.vcf.gz --key hap_key.txt --outdir ./results\n\n")
      quit(save = "no", status = 0)
    } else if (arg %in% c("-v", "--vcf") && i < length(args)) {
      parsed$vcf_file <- args[i + 1]
      i <- i + 2
    } else if (grepl("^--vcf=", arg)) {
      parsed$vcf_file <- sub("^--vcf=", "", arg)
      i <- i + 1
    } else if (arg %in% c("-k", "--key", "--mas-key") && i < length(args)) {
      parsed$mas_key_file <- args[i + 1]
      i <- i + 2
    } else if (grepl("^--(key|mas-key)=", arg)) {
      parsed$mas_key_file <- sub("^--(key|mas-key)=", "", arg)
      i <- i + 1
    } else if (arg %in% c("-o", "--outdir", "--output-dir") && i < length(args)) {
      parsed$output_directory <- args[i + 1]
      i <- i + 2
    } else if (grepl("^--(outdir|output-dir)=", arg)) {
      parsed$output_directory <- sub("^--(outdir|output-dir)=", "", arg)
      i <- i + 1
    } else if (arg %in% c("-c", "--cores") && i < length(args)) {
      parsed$cores <- as.integer(args[i + 1])
      i <- i + 2
    } else if (grepl("^--cores=", arg)) {
      parsed$cores <- as.integer(sub("^--cores=", "", arg))
      i <- i + 1
    } else if (arg == "--verbose") {
      parsed$verbose <- TRUE
      i <- i + 1
    } else if (grepl("^--verbose=", arg)) {
      parsed$verbose <- as.logical(sub("^--verbose=", "", arg))
      i <- i + 1
    } else {
      warning("Ignoring unrecognized argument: ", arg)
      i <- i + 1
    }
  }
  return(parsed)
}

cli_args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

# Resolve inputs with fallback to example dataset
vcf_file <- if (!is.null(cli_args$vcf_file) && nzchar(cli_args$vcf_file)) {
  cli_args$vcf_file
} else {
  file.path(script_dir, "example", "Example.vcf.gz")
}

mas_key_file <- if (!is.null(cli_args$mas_key_file) && nzchar(cli_args$mas_key_file)) {
  cli_args$mas_key_file
} else {
  file.path(script_dir, "example", "Example_File.txt")
}

output_directory <- if (!is.null(cli_args$output_directory) && nzchar(cli_args$output_directory)) {
  cli_args$output_directory
} else {
  file.path(script_dir, "example")
}

num_cores <- if (!is.null(cli_args$cores) && as.integer(cli_args$cores) > 0) {
  as.integer(cli_args$cores)
} else {
  max(1, parallel::detectCores() - 1)
}

verbose <- if (!is.null(cli_args$verbose)) as.logical(cli_args$verbose) else TRUE

# Validate input files
if (!file.exists(vcf_file)) {
  stop("Following VCF file does not exist: ", vcf_file)
}
if (!file.exists(mas_key_file)) {
  stop("Following MAS key file does not exist: ", mas_key_file)
}

# Ensure output directory exists
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
}

# Derive base name for outputs
output_name <- sub("\\.[^.]*$", "", basename(mas_key_file))

# ---------------------------------------------------------
# 3. Read VCF and Convert to HapMap
# ---------------------------------------------------------

if (verbose) {
  print_text <- paste("### Initializing script:", Sys.time(), "###")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  cat(print_text, "\n")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  cat("VCF File:         ", vcf_file, "\n")
  cat("MAS Key File:     ", mas_key_file, "\n")
  cat("Output Directory: ", output_directory, "\n")
  cat("Parallel Cores:   ", num_cores, "\n\n")
}

# Read VCF table
vcf <- as.data.frame(data.table::fread(input = vcf_file, skip = count_header_lines(vcf_file)))

# Extract 3-char genotype calls (e.g. 0/0, 0/1, 1/1, ./.)
vcf[, 10:ncol(vcf)] <- lapply(vcf[, 10:ncol(vcf)], function(x) {
  substr(x, 1, 3)
})

# Read MAS key
key <- read.table(mas_key_file, sep = "\t", header = TRUE, na.strings = c("", "NA"))

# Extract unique KIMs from key
key_markers <- unique(unlist(strsplit(key[, "Markers"], split = ":")))

# Subset VCF to key markers
vcf <- vcf[vcf[, "ID"] %in% key_markers, ]

if (nrow(vcf) == 0) {
  stop("None of the markers specified in the MAS key file were found in the VCF!")
}

# Start cluster for HapMap conversion
cl <- makeCluster(num_cores)
clusterCall(cl, function(libs) .libPaths(libs), .libPaths())
registerDoParallel(cl)

if (verbose) {
  print_text <- paste("### HapMap conversion start:", Sys.time(), "###")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  cat(print_text, "\n")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
}

# Convert VCF rows to HapMap format in parallel
hapmap <- foreach(row_idx = 1:nrow(vcf), .combine = "rbind", .export = c("vcf_to_hapmap")) %dopar% {
  vcf_to_hapmap(vcf[row_idx, ])
}

stopCluster(cl)
remove(cl)

# Write HapMap file
hapmap_path <- file.path(output_directory, paste0(output_name, "_hapmap.hmp.txt.gz"))
fwrite(
  as.data.table(hapmap, keep.rownames = FALSE),
  file = hapmap_path,
  sep = "\t"
)

if (verbose) {
  print_text <- paste("### HapMap conversion stop:", Sys.time(), "###")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  cat(print_text, "\n")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
}

# ---------------------------------------------------------
# 4. Generate MAS Reports
# ---------------------------------------------------------

if (verbose) {
  print_text <- paste("### Generating marker report start:", Sys.time(), "###")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  cat(print_text, "\n")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
}

# Check for duplicated loci
if (any(duplicated(key[, "Locus"]))) {
  stop("There appears to be duplicated loci in the key file. Please check your inputs!")
}

# Start cluster for locus haplotype calling
cl <- makeCluster(num_cores)
clusterCall(cl, function(libs) .libPaths(libs), .libPaths())
registerDoParallel(cl)

# Call haplotypes for each locus in parallel
reports <- foreach(
  locus_name = unique(key[, "Locus"]),
  .export = c("call_locus_haplotypes", "auto_call_haplotypes")
) %dopar% {
  call_locus_haplotypes(locus_name = locus_name, key = key, hapmap = hapmap)
}

stopCluster(cl)
remove(cl)

# Merge parallel locus reports
merged_outputs <- merge_locus_reports(reports)
marker_report <- merged_outputs$marker_report
auto_suggestions <- merged_outputs$auto_suggestions
remove(merged_outputs)

if (verbose) {
  print_text <- paste("### Generating marker report stop:", Sys.time(), "###")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  cat(print_text, "\n")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  print_text <- paste("### Writing outputs start:", Sys.time(), "###")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  cat(print_text, "\n")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
}

# 1. Full report
write.csv(
  marker_report,
  file.path(output_directory, paste0(output_name, "_Marker_Report_Full.csv")),
  row.names = FALSE,
  quote = TRUE
)

# 2. Summary report (annotated calls only)
columns_to_pull <- c("Line", colnames(marker_report)[grep("Summary_Annotated", colnames(marker_report))])
marker_summary <- marker_report[, columns_to_pull]
colnames(marker_summary) <- gsub("_Summary_Annotated", "", colnames(marker_summary))

write.csv(
  marker_summary,
  file.path(output_directory, paste0(output_name, "_Marker_Report_Summary.csv")),
  row.names = FALSE,
  quote = TRUE
)

# 3. Frequency statistics
marker_freq_report <- calculate_marker_statistics(marker_report)
write.csv(
  marker_freq_report,
  file.path(output_directory, paste0(output_name, "_Marker_Report_Statistics.csv")),
  row.names = FALSE,
  quote = TRUE
)

# 4. Unidentified haplotypes (if any)
missing_haplotypes <- find_unidentified_haplotypes(
  marker_report = marker_report,
  marker_summary = marker_summary,
  key = key,
  key_markers = key_markers
)

if (!is.null(missing_haplotypes)) {
  write.csv(
    missing_haplotypes,
    file.path(output_directory, paste0(output_name, "_Unidentified_Haplotypes.csv")),
    row.names = FALSE,
    quote = TRUE
  )
}

# 5. Auto-suggested haplotypes (if any)
if (!is.null(auto_suggestions)) {
  write.csv(
    auto_suggestions,
    file.path(output_directory, paste0(output_name, "_Auto_Suggested_Haplotypes.csv")),
    row.names = FALSE,
    quote = TRUE
  )
}

if (verbose) {
  print_text <- paste("### Writing outputs stop:", Sys.time(), "###")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  cat(print_text, "\n")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  print_text <- paste("### Script finished time:", Sys.time(), "###")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
  cat(print_text, "\n")
  cat(paste(rep("#", nchar(print_text)), collapse = ""), "\n")
}
