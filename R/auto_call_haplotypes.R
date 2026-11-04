# Perform automatic haplotype expansion and calling for a multi-marker locus
auto_call_haplotypes <- function(
    locus_key_row,
    locus_markers,
    locus_genotypes,
    marker_alleles,
    pos_haplotypes,
    het_haplotypes,
    neg_haplotypes,
    missing_haplotypes,
    locus_name
) {
  auto_warnings <- character(0)

  # Check if positive calls are given
  if (is.na(pos_haplotypes[1])) {
    stop(
      "It appears there are no positive haplotypes given for ",
      locus_name,
      ". There must be at least one positive haplotype defined for the 'Auto' feature.",
      " Please check your inputs and resubmit!"
    )
  }

  # Pull number of tolerated missing markers
  max_missing_allowed <- locus_key_row[, "N_Missing"]

  # Check if number of missing is defined
  if (is.na(max_missing_allowed) | !is.numeric(max_missing_allowed)) {
    stop(
      "Number of tolerated missing calls must be defined and must be a number for ",
      locus_name,
      "! Check the Column 'N_Missing' and resubmit!"
    )
  }

  # Check if the locus only has one marker
  if (length(locus_markers) == 1) {
    stop(
      "Only one defining marker detected for ",
      locus_name,
      ". The 'Auto' function is only available for loci with multiple markers! ",
      "Please define your haplotypes, set auto to 'FALSE', and resubmit."
    )
  }

  # Split user-supplied haplotypes by marker position to check if '-' or '0' is explicitly used
  all_user_haps <- na.omit(c(pos_haplotypes, het_haplotypes, neg_haplotypes))
  split_user_haps <- strsplit(all_user_haps, split = ":")
  valid_len_haps <- split_user_haps[vapply(split_user_haps, length, integer(1)) == length(locus_markers)]
  user_hap_matrix <- if (length(valid_len_haps) > 0) {
    matrix(unlist(valid_len_haps), ncol = length(locus_markers), byrow = TRUE)
  } else {
    matrix(character(0), ncol = length(locus_markers))
  }
  colnames(user_hap_matrix) <- locus_markers

  # Extract REF/ALT alleles of each marker in IUPAC format
  ref_alt <- c()

  for (marker_id in rownames(locus_genotypes)) {
    allele_tokens <- unlist(strsplit(marker_alleles[marker_id], "/"))
    allele_g <- allele_tokens[grep("G", allele_tokens)]
    allele_t <- allele_tokens[grep("T", allele_tokens)]
    allele_c <- allele_tokens[grep("C", allele_tokens)]
    allele_a <- allele_tokens[grep("A", allele_tokens)]
    allele_gap <- allele_tokens[grep("-", allele_tokens)]

    # If a marker has '-' in VCF alleles (e.g., 'A/-' or '-/A'), check if '-' or '0' is ever
    # observed in the HapMap or specified in the key file. If neither, '-' is a placeholder for
    # a dominant presence/absence amplicon where absence is called 'N'.
    if (length(allele_gap) > 0) {
      obs_vals <- unique(as.character(locus_genotypes[marker_id, ]))
      user_vals <- if (nrow(user_hap_matrix) > 0) unique(user_hap_matrix[, marker_id]) else character(0)
      if (!any(c("-", "0") %in% c(obs_vals, user_vals))) {
        allele_gap <- character(0)
      }
    }

    marker_allele_df <- data.frame(
      Marker = marker_id,
      G = ifelse(length(allele_g) == 0, NA, allele_g),
      `T` = ifelse(length(allele_t) == 0, NA, allele_t),
      C = ifelse(length(allele_c) == 0, NA, allele_c),
      A = ifelse(length(allele_a) == 0, NA, allele_a),
      gap = ifelse(length(allele_gap) == 0, NA, allele_gap)
    )
    marker_allele_df <- marker_allele_df[, colSums(is.na(marker_allele_df)) != nrow(marker_allele_df), drop = FALSE]

    # Check if the marker is biallelic
    if (ncol(marker_allele_df) != 3 & ncol(marker_allele_df) != 2) {
      stop(
        "Marker ",
        marker_id,
        " associated with locus ",
        locus_name,
        " appears to be multiallelic rather than biallelic.",
        " Calls can only be made in Auto mode with biallelic markers!",
        " Please check your inputs and resubmit!"
      )
    } else if (ncol(marker_allele_df) == 2) {
      message("#################################")
      message(
        "Marker ",
        marker_id,
        " associated with locus ",
        locus_name,
        " appears to be monomorphic rather than biallelic.",
        " Moving forward with assumption that the marker is dominant..."
      )
      marker_allele_df <- cbind(marker_allele_df, NA)
      colnames(marker_allele_df) <- c("Marker", "Allele_1", "Allele_2")
      ref_alt <- rbind(ref_alt, marker_allele_df)
    } else {
      colnames(marker_allele_df) <- c("Marker", "Allele_1", "Allele_2")
      ref_alt <- rbind(ref_alt, marker_allele_df)
    }
  }

  # Assign heterozygous call
  ref_alt$Allele_Het <- paste(
    ref_alt[, "Allele_1"],
    "/",
    ref_alt[, "Allele_2"],
    sep = ""
  )
  ref_alt[ref_alt == "A/G"] <- "R"
  ref_alt[ref_alt == "G/A"] <- "R"
  ref_alt[ref_alt == "C/T"] <- "Y"
  ref_alt[ref_alt == "T/C"] <- "Y"
  ref_alt[ref_alt == "G/C"] <- "S"
  ref_alt[ref_alt == "C/G"] <- "S"
  ref_alt[ref_alt == "A/T"] <- "W"
  ref_alt[ref_alt == "T/A"] <- "W"
  ref_alt[ref_alt == "G/T"] <- "K"
  ref_alt[ref_alt == "T/G"] <- "K"
  ref_alt[ref_alt == "A/C"] <- "M"
  ref_alt[ref_alt == "C/A"] <- "M"

  # Assign indel heterozygous calls
  ref_alt[ref_alt == "A/-"] <- "0"
  ref_alt[ref_alt == "-/A"] <- "0"
  ref_alt[ref_alt == "T/-"] <- "0"
  ref_alt[ref_alt == "-/T"] <- "0"
  ref_alt[ref_alt == "G/-"] <- "0"
  ref_alt[ref_alt == "-/G"] <- "0"
  ref_alt[ref_alt == "C/-"] <- "0"
  ref_alt[ref_alt == "-/C"] <- "0"

  # Assign dominant calls
  ref_alt[ref_alt == "A/NA"] <- NA
  ref_alt[ref_alt == "G/NA"] <- NA
  ref_alt[ref_alt == "C/NA"] <- NA
  ref_alt[ref_alt == "T/NA"] <- NA
  ref_alt[ref_alt == "-/NA"] <- NA

  # Define positive allele candidates
  pos_candidates <- data.frame(Haplotype = pos_haplotypes, stringsAsFactors = FALSE)
  pos_candidates <- cbind(pos_candidates, matrix(
    data = unlist(strsplit(pos_haplotypes, split = ":")),
    nrow = nrow(pos_candidates),
    byrow = TRUE
  ), stringsAsFactors = FALSE)
  colnames(pos_candidates) <- c("Haplotype", ref_alt[, "Marker"])

  # First look for strictly complete homozygous positive haplotypes (no N, no 0, and no IUPAC het codes)
  positive <- pos_candidates
  positive[positive == "N"] <- NA
  positive[positive == "0"] <- NA
  positive[positive == "R"] <- NA
  positive[positive == "Y"] <- NA
  positive[positive == "S"] <- NA
  positive[positive == "W"] <- NA
  positive[positive == "K"] <- NA
  positive[positive == "M"] <- NA
  positive <- na.omit(positive)

  # If no strictly non-N positive case exists, check if a dominant/null-positive ('N') haplotype was supplied
  if (nrow(positive) == 0) {
    dom_positive <- pos_candidates
    dom_positive[dom_positive == "0"] <- NA
    dom_positive[dom_positive == "R"] <- NA
    dom_positive[dom_positive == "Y"] <- NA
    dom_positive[dom_positive == "S"] <- NA
    dom_positive[dom_positive == "W"] <- NA
    dom_positive[dom_positive == "K"] <- NA
    dom_positive[dom_positive == "M"] <- NA
    dom_positive <- na.omit(dom_positive)

    if (nrow(dom_positive) > 0) {
      not_all_n <- apply(dom_positive[, 2:ncol(dom_positive), drop = FALSE], 1, function(r) any(r != "N"))
      if (any(not_all_n)) {
        dom_positive <- dom_positive[not_all_n, , drop = FALSE]
      }
    }

    if (nrow(dom_positive) == 0) {
      stop(
        "It appears that no positive cases for ",
        locus_name,
        " are complete case!",
        " Positive cases must contain at least one complete, homozygous sequence (e.g., A:G:T).",
        " Check your inputs and resubmit!"
      )
    }

    positive <- dom_positive[1, , drop = FALSE]
    auto_warnings <- c(
      auto_warnings,
      paste0(
        "Warning: Locus '", locus_name, "' uses a null/dominant ('N') allele in its positive haplotype (",
        positive[1, "Haplotype"], "). Treating 'N' at that marker position as a dominant null-positive state.\n"
      )
    )
  } else {
    # If multiple positive haplotypes were provided in Pos, warn the user and use the first complete case
    if (length(pos_haplotypes) > 1) {
      auto_warnings <- c(
        auto_warnings,
        paste0(
          "Warning: Multiple positive haplotypes were provided for locus '", locus_name, "' (",
          paste(pos_haplotypes, collapse = ", "), "). Using the first complete positive case (",
          positive[1, "Haplotype"], ") for automatic haplotype calling; all other provided positive cases will override automatic calls.\n"
        )
      )
    }
    positive <- positive[1, , drop = FALSE]
  }

  # Define negative allele from the single seed positive haplotype
  negative <- positive
  null_pos_markers <- character(0)

  for (marker_id in colnames(positive)[2:ncol(positive)]) {
    pos_allele <- negative[1, marker_id]
    marker_info <- ref_alt[ref_alt[, "Marker"] == marker_id, ]

    if (pos_allele == "N") {
      # Null-positive dominant marker: negative state is the amplified allele
      cand_alleles <- na.omit(c(marker_info[, "Allele_1"], marker_info[, "Allele_2"]))
      cand_alleles <- cand_alleles[cand_alleles != "-"]
      if (length(cand_alleles) > 1) {
        obs_non_n <- setdiff(unique(as.character(locus_genotypes[marker_id, ])), c("N", "-"))
        if (length(obs_non_n) == 1 && obs_non_n %in% cand_alleles) {
          neg_allele <- obs_non_n
        } else if (nrow(user_hap_matrix) > 0) {
          user_non_n <- intersect(unique(user_hap_matrix[, marker_id]), cand_alleles)
          neg_allele <- if (length(user_non_n) >= 1) user_non_n[1] else cand_alleles[1]
        } else {
          neg_allele <- cand_alleles[1]
        }
      } else {
        neg_allele <- cand_alleles[1]
      }
      null_pos_markers <- c(null_pos_markers, marker_id)
      # For a null-positive dominant marker, heterozygotes also amplify the present allele
      ref_alt[ref_alt[, "Marker"] == marker_id, "Allele_Het"] <- neg_allele
    } else {
      neg_allele <- ifelse(
        marker_info[, "Allele_1"] == pos_allele,
        marker_info[, "Allele_2"],
        marker_info[, "Allele_1"]
      )
      neg_allele <- ifelse(is.na(neg_allele), "N", neg_allele)
    }

    negative[1, marker_id] <- neg_allele
  }
  negative[1, "Haplotype"] <- paste(as.character(negative[1, 2:ncol(negative)]), collapse = ":")

  # Define heterozygous cases
  het_allele_vec <- as.character(ref_alt[, "Allele_Het"])
  heterozygous <- data.frame(
    Haplotype = paste(het_allele_vec, collapse = ":"),
    t(ref_alt[, "Allele_Het"]),
    stringsAsFactors = FALSE
  )
  colnames(heterozygous)[2:ncol(heterozygous)] <- ref_alt[, "Marker"]

  # If there are NAs which indicate that there is a dominant marker
  if (length(unique(is.na(heterozygous))) != 1) {
    het_vec <- c()
    for (row_idx in 1:nrow(ref_alt)) {
      if (is.na(ref_alt[row_idx, "Allele_Het"])) {
        het_vec <- c(het_vec, ref_alt[row_idx, "Allele_1"])
      } else {
        het_vec <- c(het_vec, ref_alt[row_idx, "Allele_Het"])
      }
    }
    heterozygous <- data.frame(
      Haplotype = paste(het_vec, collapse = ":"),
      t(het_vec),
      stringsAsFactors = FALSE
    )
    colnames(heterozygous)[2:ncol(heterozygous)] <- ref_alt[, "Marker"]
  }

  # Define missing case
  missing <- data.frame(
    Haplotype = paste(rep("N", length(2:ncol(heterozygous))), collapse = ":"),
    t(rep("N", length(2:ncol(heterozygous)))),
    stringsAsFactors = FALSE
  )
  colnames(missing)[2:ncol(missing)] <- ref_alt[, "Marker"]

  # Check if negative cases have N (amplified-positive dominant markers)
  dominant_markers <- character(0)
  if ("N" %in% negative[, 2:ncol(negative)]) {
    neg_marker_cols <- negative[, 2:ncol(negative), drop = FALSE]
    dominant_markers <- colnames(neg_marker_cols)[as.character(neg_marker_cols[1, ]) == "N"]
  }
  all_dominant_markers <- unique(c(dominant_markers, null_pos_markers))

  # Combine the 4 canonical cases (Positive, Heterozygous, Negative, Missing)
  cases <- data.frame(
    Type = c("Positive", "Heterozygous", "Negative", "Missing"),
    rbind(positive, heterozygous, negative, missing),
    stringsAsFactors = FALSE
  )

  # Expand every possible haplotype combination across markers
  haplotypes <- expand.grid(cases[, 3:ncol(cases), drop = FALSE], stringsAsFactors = FALSE)
  haplotypes[, 1:ncol(haplotypes)] <- lapply(haplotypes[, 1:ncol(haplotypes), drop = FALSE], as.character)
  haplotypes <- unique(haplotypes)
  haplotypes$Haplotype <- apply(haplotypes[, 1:ncol(haplotypes), drop = FALSE], 1, paste, collapse = ":")

  # Ensure any extra user-specified positive haplotypes with valid marker length are included
  extra_pos_haps <- setdiff(na.omit(unique(pos_haplotypes)), haplotypes$Haplotype)
  if (length(extra_pos_haps) > 0) {
    split_extra <- strsplit(extra_pos_haps, split = ":")
    keep_extra <- vapply(split_extra, length, integer(1)) == length(locus_markers)
    extra_pos_haps <- extra_pos_haps[keep_extra]
    split_extra <- split_extra[keep_extra]
    if (length(extra_pos_haps) > 0) {
      extra_mat <- as.data.frame(
        matrix(unlist(split_extra), nrow = length(extra_pos_haps), byrow = TRUE),
        stringsAsFactors = FALSE
      )
      colnames(extra_mat) <- colnames(haplotypes)[1:(ncol(haplotypes) - 1)]
      extra_mat$Haplotype <- extra_pos_haps
      haplotypes <- rbind(haplotypes, extra_mat)
    }
  }

  # Give an initial call to haplotypes explicitly defined by the user
  haplotypes[, "Call"] <- ifelse(haplotypes[, "Haplotype"] %in% pos_haplotypes, "POS",
    ifelse(haplotypes[, "Haplotype"] %in% het_haplotypes, "HET",
      ifelse(haplotypes[, "Haplotype"] %in% neg_haplotypes, "NEG",
        ifelse(haplotypes[, "Haplotype"] %in% missing_haplotypes, "MISSING", NA)
      )
    )
  )

  # Count number of missing ('N') markers in each haplotype
  marker_cols <- ref_alt[, "Marker"]
  codom_cols <- setdiff(marker_cols, all_dominant_markers)
  all_pos_is_n <- all(as.character(positive[1, marker_cols]) == "N")

  if (length(all_dominant_markers) == 0) {
    N_Missing <- vapply(haplotypes[, "Haplotype"], function(h) {
      sum(unlist(gregexpr("N", h, fixed = TRUE)) > 0)
    }, integer(1))
  } else {
    N_Missing <- apply(haplotypes[, marker_cols, drop = FALSE], 1, function(row_vals) {
      if (all(row_vals == "N")) {
        # If the positive haplotype itself is all 'N' (full deletion), do not count all-'N' as missing
        return(if (all_pos_is_n) 0L else length(marker_cols))
      } else if (length(codom_cols) == 0) {
        return(0L)
      } else {
        return(sum(row_vals[codom_cols] == "N"))
      }
    })
  }

  haplotypes[, "N_Missing"] <- N_Missing

  # Assign MISSING if N_Missing exceeds max_missing_allowed
  haplotypes[, "Call"] <- ifelse(
    haplotypes[, "N_Missing"] > max_missing_allowed,
    "MISSING",
    haplotypes[, "Call"]
  )

  # Classify positive + missing combinations as POS
  pos_candidates_grid <- expand.grid(cases[c(1, 4), c(3:ncol(cases)), drop = FALSE], stringsAsFactors = FALSE)
  pos_candidates_grid[, 1:ncol(pos_candidates_grid)] <- lapply(
    pos_candidates_grid[, 1:ncol(pos_candidates_grid), drop = FALSE],
    as.character
  )
  pos_candidates_grid[, "Haplotype"] <- apply(
    pos_candidates_grid[, 1:ncol(pos_candidates_grid), drop = FALSE],
    1,
    paste,
    collapse = ":"
  )

  if (length(dominant_markers) > 0) {
    for (dm in dominant_markers) {
      pos_candidates_grid <- pos_candidates_grid[pos_candidates_grid[, dm] == positive[1, dm], , drop = FALSE]
    }
  }

  already_called_haps <- haplotypes[!is.na(haplotypes[, "Call"]), "Haplotype"]
  pos_candidates_grid <- pos_candidates_grid[!pos_candidates_grid[, "Haplotype"] %in% already_called_haps, , drop = FALSE]
  haplotypes[haplotypes[, "Haplotype"] %in% pos_candidates_grid[, "Haplotype"], "Call"] <- "POS"

  # Classify positive + heterozygous + missing combinations as HET
  if (!all_pos_is_n) {
    het_candidates_grid <- expand.grid(cases[c(1, 2, 4), c(3:ncol(cases)), drop = FALSE], stringsAsFactors = FALSE)
    het_candidates_grid[, 1:ncol(het_candidates_grid)] <- lapply(
      het_candidates_grid[, 1:ncol(het_candidates_grid), drop = FALSE],
      as.character
    )
    het_candidates_grid[, "Haplotype"] <- apply(
      het_candidates_grid[, 1:ncol(het_candidates_grid), drop = FALSE],
      1,
      paste,
      collapse = ":"
    )
    het_candidates_grid <- het_candidates_grid[het_candidates_grid[, "Haplotype"] != negative[1, "Haplotype"], , drop = FALSE]

    already_called_haps <- haplotypes[!is.na(haplotypes[, "Call"]), "Haplotype"]
    het_candidates_grid <- het_candidates_grid[!het_candidates_grid[, "Haplotype"] %in% already_called_haps, , drop = FALSE]
    haplotypes[haplotypes[, "Haplotype"] %in% het_candidates_grid[, "Haplotype"], "Call"] <- "HET"
  }

  # Classify remaining combinations containing negative alleles as NEG
  neg_candidates_grid <- expand.grid(cases[c(1, 2, 3, 4), c(3:ncol(cases)), drop = FALSE], stringsAsFactors = FALSE)
  neg_candidates_grid[, 1:ncol(neg_candidates_grid)] <- lapply(
    neg_candidates_grid[, 1:ncol(neg_candidates_grid), drop = FALSE],
    as.character
  )
  neg_candidates_grid[, "Haplotype"] <- apply(
    neg_candidates_grid[, 1:ncol(neg_candidates_grid), drop = FALSE],
    1,
    paste,
    collapse = ":"
  )

  already_called_haps <- haplotypes[!is.na(haplotypes[, "Call"]), "Haplotype"]
  neg_candidates_grid <- neg_candidates_grid[!neg_candidates_grid[, "Haplotype"] %in% already_called_haps, , drop = FALSE]
  if (nrow(neg_candidates_grid) > 0) {
    haplotypes[haplotypes[, "Haplotype"] %in% neg_candidates_grid[, "Haplotype"], "Call"] <- "NEG"
  }

  # Ensure all user-provided positive cases override automatic calls
  if (length(pos_haplotypes) > 1) {
    haplotypes[haplotypes[, "Haplotype"] %in% pos_haplotypes, "Call"] <- "POS"
  }

  # Reconcatenate classified haplotypes back into the locus key row
  for (call_category in unique(haplotypes[, "Call"])) {
    category_haps <- as.character(haplotypes[haplotypes[, "Call"] == call_category, "Haplotype"])
    collapsed_haps <- paste(category_haps, collapse = "_")
    if (call_category == "POS") {
      locus_key_row[, "Pos"] <- collapsed_haps
    } else if (call_category == "HET") {
      locus_key_row[, "Het"] <- collapsed_haps
    } else if (call_category == "NEG") {
      locus_key_row[, "Neg"] <- collapsed_haps
    } else if (call_category == "MISSING") {
      locus_key_row[, "Missing"] <- collapsed_haps
    }
  }

  return(list(
    updated_key_row = locus_key_row,
    warnings = auto_warnings
  ))
}
