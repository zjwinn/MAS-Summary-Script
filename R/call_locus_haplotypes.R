# Call haplotypes across all lines in the HapMap for a single locus
call_locus_haplotypes <- function(locus_name, key, hapmap) {
  # Pull key information for this locus
  locus_key_row <- key[key[, "Locus"] == locus_name, ]

  # Parse defining marker IDs
  locus_markers <- unlist(strsplit(locus_key_row[, "Markers"], ":"))

  # Pull markers from HapMap
  locus_hapmap <- as.data.frame(hapmap[hapmap[, "rs#"] %in% locus_markers, ])
  rownames(locus_hapmap) <- locus_hapmap[, 1]

  # Extract marker allele string from the second column ('alleles')
  marker_alleles <- c(locus_hapmap[, 2])
  names(marker_alleles) <- locus_hapmap[, 1]

  # Subset to sample genotype columns in marker order
  locus_genotypes <- locus_hapmap[, 12:ncol(locus_hapmap), drop = FALSE]
  locus_genotypes <- locus_genotypes[locus_markers, , drop = FALSE]

  # Check if all markers for this locus are present in the HapMap
  if (all(locus_markers %in% rownames(locus_genotypes)) == FALSE) {
    warning_message <- paste(
      "Warning: Markers identified for ",
      locus_name,
      " are not found HapMap.\n",
      " Not producing call for locus ",
      locus_name,
      " and moving to the next locus.\n",
      " If this is in error, check your VCF input marker names and your haplotyping key!\n",
      sep = ""
    )
    return(warning_message)
  }

  # Parse user-specified haplotype strings
  pos_haplotypes <- if (is.na(locus_key_row[, "Pos"])) NA else unlist(strsplit(locus_key_row[, "Pos"], "_"))
  het_haplotypes <- if (is.na(locus_key_row[, "Het"])) NA else unlist(strsplit(locus_key_row[, "Het"], "_"))
  neg_haplotypes <- if (is.na(locus_key_row[, "Neg"])) NA else unlist(strsplit(locus_key_row[, "Neg"], "_"))
  missing_haplotypes <- if (is.na(locus_key_row[, "Missing"])) NA else unlist(strsplit(locus_key_row[, "Missing"], "_"))

  auto_suggestions <- c()
  auto_warnings <- character(0)

  # Check if automatic haplotype calling is requested
  if (locus_key_row[, "Auto"] == TRUE) {
    auto_result <- auto_call_haplotypes(
      locus_key_row = locus_key_row,
      locus_markers = locus_markers,
      locus_genotypes = locus_genotypes,
      marker_alleles = marker_alleles,
      pos_haplotypes = pos_haplotypes,
      het_haplotypes = het_haplotypes,
      neg_haplotypes = neg_haplotypes,
      missing_haplotypes = missing_haplotypes,
      locus_name = locus_name
    )
    locus_key_row <- auto_result$updated_key_row
    auto_warnings <- auto_result$warnings

    # Re-parse expanded haplotype strings
    pos_haplotypes <- if (is.na(locus_key_row[, "Pos"])) NA else unlist(strsplit(locus_key_row[, "Pos"], "_"))
    het_haplotypes <- if (is.na(locus_key_row[, "Het"])) NA else unlist(strsplit(locus_key_row[, "Het"], "_"))
    neg_haplotypes <- if (is.na(locus_key_row[, "Neg"])) NA else unlist(strsplit(locus_key_row[, "Neg"], "_"))
    missing_haplotypes <- if (is.na(locus_key_row[, "Missing"])) NA else unlist(strsplit(locus_key_row[, "Missing"], "_"))

    auto_suggestions <- locus_key_row
  }

  # Pull optional annotations
  pos_annotation <- locus_key_row[, "Pos_Annot"]
  het_annotation <- locus_key_row[, "Het_Annot"]
  neg_annotation <- locus_key_row[, "Neg_Annot"]
  missing_annotation <- locus_key_row[, "Missing_Annot"]

  locus_report <- NULL

  # Classify every individual in the HapMap
  for (line_id in colnames(locus_genotypes)) {
    line_marker_calls <- locus_genotypes[, line_id]
    line_haplotype <- paste(line_marker_calls, collapse = ":")

    line_call <- ifelse(line_haplotype %in% pos_haplotypes, "POS",
      ifelse(line_haplotype %in% het_haplotypes, "HET",
        ifelse(line_haplotype %in% neg_haplotypes, "NEG",
          ifelse(line_haplotype %in% missing_haplotypes, "NA", "UNIDENTIFIED")
        )
      )
    )

    line_report_row <- data.frame(
      Line = line_id,
      t(line_marker_calls),
      line_haplotype,
      line_call,
      line_call
    )

    # Apply custom annotations
    if (is.na(line_report_row[, ncol(line_report_row)]) == TRUE & is.na(missing_annotation) == FALSE) {
      line_report_row[, ncol(line_report_row)] <- missing_annotation
    } else if (as.character(line_report_row[, ncol(line_report_row)]) == "POS" & is.na(pos_annotation) == FALSE) {
      line_report_row[, ncol(line_report_row)] <- pos_annotation
    } else if (as.character(line_report_row[, ncol(line_report_row)]) == "HET" & is.na(het_annotation) == FALSE) {
      line_report_row[, ncol(line_report_row)] <- het_annotation
    } else if (as.character(line_report_row[, ncol(line_report_row)]) == "NEG" & is.na(neg_annotation) == FALSE) {
      line_report_row[, ncol(line_report_row)] <- neg_annotation
    }

    colnames(line_report_row) <- c(
      "Line",
      paste(locus_name, rownames(locus_genotypes), sep = "_"),
      paste(locus_name, "_Haplotype", sep = ""),
      paste(locus_name, "_Summary", sep = ""),
      paste(locus_name, "_Summary_Annotated", sep = "")
    )

    if (is.null(locus_report)) {
      locus_report <- line_report_row
    } else {
      locus_report <- rbind(locus_report, line_report_row)
    }
  }

  if (is.data.frame(auto_suggestions) == TRUE) {
    return(list(
      Report = locus_report,
      AutoSuggestions = locus_key_row,
      Warnings = auto_warnings
    ))
  } else {
    return(locus_report)
  }
}
