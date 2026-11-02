# Merge parallel locus reports, calculate summary statistics, and identify missing haplotypes

merge_locus_reports <- function(reports) {
  marker_report <- NULL
  auto_suggestions <- NULL

  for (report_idx in seq_along(reports)) {
    locus_result <- reports[[report_idx]]

    # If the report is a warning message (character)
    if (is.character(locus_result) == TRUE) {
      cat("############################\n")
      cat(locus_result)
      cat("############################\n")
      next
    } else if (is.data.frame(locus_result) == TRUE) {
      if (is.null(marker_report)) {
        marker_report <- locus_result
      } else {
        marker_report <- suppressWarnings(merge(marker_report, locus_result, by = "Line"))
      }
    } else if (is.list(locus_result) == TRUE) {
      if (!is.null(locus_result$Warnings) && length(locus_result$Warnings) > 0) {
        for (warn_msg in locus_result$Warnings) {
          cat("############################\n")
          cat(warn_msg)
          cat("############################\n")
        }
      }
      if (is.null(auto_suggestions)) {
        auto_suggestions <- locus_result[["AutoSuggestions"]]
      } else {
        auto_suggestions <- rbind(auto_suggestions, locus_result[["AutoSuggestions"]])
      }
      if (is.null(marker_report)) {
        marker_report <- locus_result[["Report"]]
      } else {
        marker_report <- suppressWarnings(merge(marker_report, locus_result[["Report"]], by = "Line"))
      }
    }
  }

  if (is.null(marker_report)) {
    stop("No locus reports could be generated! Check that the markers in your key file match your VCF.")
  }

  # Convert string "NA" into NA
  marker_report[marker_report == "NA"] <- NA

  return(list(
    marker_report = marker_report,
    auto_suggestions = auto_suggestions
  ))
}

calculate_marker_statistics <- function(marker_report) {
  # Pull summary column names without annotation
  columns_to_pull <- c("Line", colnames(marker_report)[grep("Summary", colnames(marker_report))])
  columns_to_pull <- columns_to_pull[-grep("Annotated", columns_to_pull)]

  marker_calls <- marker_report[, columns_to_pull, drop = FALSE]
  marker_freq_report <- NULL

  for (summary_col in colnames(marker_calls)[2:ncol(marker_calls)]) {
    locus_call_vec <- marker_calls[, c(summary_col)]
    call_counts_table <- data.table::as.data.table(table(locus_call_vec))

    locus_stats_row <- data.frame(
      gsub("_Summary", "", summary_col),
      ifelse(nrow(call_counts_table[call_counts_table$locus_call_vec == "POS", "N"]) == 0, 0, call_counts_table[call_counts_table$locus_call_vec == "POS", "N"]),
      ifelse(nrow(call_counts_table[call_counts_table$locus_call_vec == "HET", "N"]) == 0, 0, call_counts_table[call_counts_table$locus_call_vec == "HET", "N"]),
      ifelse(nrow(call_counts_table[call_counts_table$locus_call_vec == "NEG", "N"]) == 0, 0, call_counts_table[call_counts_table$locus_call_vec == "NEG", "N"]),
      sum(is.na(locus_call_vec)),
      ifelse(
        nrow(call_counts_table[call_counts_table$locus_call_vec == "UNIDENTIFIED", "N"]) == 0,
        0,
        call_counts_table[call_counts_table$locus_call_vec == "UNIDENTIFIED", "N"]
      ),
      length(locus_call_vec)
    )
    colnames(locus_stats_row) <- c(
      "Locus",
      "Positive_Count",
      "Heterozygous_Count",
      "Negative_Count",
      "Missing_Count",
      "Unidentified_Count",
      "Total_Sampled"
    )

    # Calculate frequencies
    locus_stats_row$Positive_Freq <- round(locus_stats_row$Positive_Count / locus_stats_row$Total_Sampled, 2)
    locus_stats_row$Heterozygous_Freq <- round(locus_stats_row$Heterozygous_Count / locus_stats_row$Total_Sampled, 2)
    locus_stats_row$Negative_Freq <- round(locus_stats_row$Negative_Count / locus_stats_row$Total_Sampled, 2)
    locus_stats_row$Missing_Freq <- round(locus_stats_row$Missing_Count / locus_stats_row$Total_Sampled, 2)
    locus_stats_row$Unidentified_Freq <- round(locus_stats_row$Unidentified_Count / locus_stats_row$Total_Sampled, 2)

    if (is.null(marker_freq_report)) {
      marker_freq_report <- locus_stats_row
    } else {
      marker_freq_report <- rbind(marker_freq_report, locus_stats_row)
    }
  }

  return(marker_freq_report)
}

find_unidentified_haplotypes <- function(marker_report, marker_summary, key, key_markers) {
  unidentified_report <- c()

  for (locus_name in colnames(marker_summary)[2:ncol(marker_summary)]) {
    locus_markers <- unlist(strsplit(key[key[, "Locus"] == locus_name, "Markers"], ":"))
    expected_cols <- c(
      "Line",
      paste(locus_name, locus_markers, sep = "_"),
      paste(locus_name, "_Haplotype", sep = ""),
      paste(locus_name, "_Summary", sep = ""),
      paste(locus_name, "_Summary_Annotated", sep = "")
    )
    present_cols <- intersect(expected_cols, colnames(marker_report))

    if (length(present_cols) <= 1) {
      next
    }

    locus_df <- as.data.frame(marker_report[, present_cols, drop = FALSE])

    # Strip locus prefix from column names
    escaped_locus <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", locus_name)
    colnames(locus_df) <- sub(paste("^", escaped_locus, "_", sep = ""), "", colnames(locus_df))

    # Filter to UNIDENTIFIED rows
    unidentified_rows <- na.omit(locus_df[locus_df[, "Summary"] == "UNIDENTIFIED", , drop = FALSE])

    if (nrow(unidentified_rows) != 0) {
      locus_key_row <- key[key[, "Locus"] == locus_name, ]
      unique_unidentified_haps <- unique(unidentified_rows[, "Haplotype"])

      locus_unidentified_df <- data.frame(
        locus_key_row,
        Unidentified_Haplotype = unique_unidentified_haps,
        row.names = NULL
      )

      unidentified_report <- rbind(unidentified_report, locus_unidentified_df)
    }
  }

  return(unidentified_report)
}
