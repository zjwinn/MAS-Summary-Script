# Convert a single VCF row into HapMap format
vcf_to_hapmap <- function(vcf_row, verbose = FALSE) {
  ref_allele <- vcf_row[, "REF"]
  alt_allele <- vcf_row[, "ALT"]
  marker_id <- vcf_row[, "ID"]
  chrom_name <- vcf_row[, "#CHROM"]
  marker_pos <- vcf_row[, "POS"]
  genotype_calls <- vcf_row[, 10:ncol(vcf_row), drop = FALSE]

  # Normalize phased genotypes (|) to unphased (/)
  genotype_calls[, 1:ncol(genotype_calls)] <- lapply(
    genotype_calls[, 1:ncol(genotype_calls), drop = FALSE],
    function(col) gsub("\\|", "/", col)
  )

  # If the marker is a SNP
  if (ref_allele %in% c("A", "T", "G", "C") & alt_allele %in% c("A", "T", "G", "C")) {
    if (verbose == TRUE) {
      print("SNP")
    }
    genotype_calls[, 1:ncol(genotype_calls)] <- lapply(
      genotype_calls[, 1:ncol(genotype_calls), drop = FALSE],
      function(col) gsub(0, ref_allele, col)
    )
    genotype_calls[, 1:ncol(genotype_calls)] <- lapply(
      genotype_calls[, 1:ncol(genotype_calls), drop = FALSE],
      function(col) gsub(1, alt_allele, col)
    )
    genotype_calls[genotype_calls == "A/A"] <- "A"
    genotype_calls[genotype_calls == "T/T"] <- "T"
    genotype_calls[genotype_calls == "C/C"] <- "C"
    genotype_calls[genotype_calls == "G/G"] <- "G"
    genotype_calls[genotype_calls == "A/G"] <- "R"
    genotype_calls[genotype_calls == "G/A"] <- "R"
    genotype_calls[genotype_calls == "C/T"] <- "Y"
    genotype_calls[genotype_calls == "T/C"] <- "Y"
    genotype_calls[genotype_calls == "G/C"] <- "S"
    genotype_calls[genotype_calls == "C/G"] <- "S"
    genotype_calls[genotype_calls == "A/T"] <- "W"
    genotype_calls[genotype_calls == "T/A"] <- "W"
    genotype_calls[genotype_calls == "G/T"] <- "K"
    genotype_calls[genotype_calls == "T/G"] <- "K"
    genotype_calls[genotype_calls == "A/C"] <- "M"
    genotype_calls[genotype_calls == "C/A"] <- "M"
    genotype_calls[genotype_calls == "./."] <- "N"

    hapmap_row <- data.frame(
      `rs#` = marker_id,
      alleles = paste(ref_allele, "/", alt_allele, sep = ""),
      chrom = chrom_name,
      pos = marker_pos,
      strand = "+",
      `assembly#` = "NA",
      center = "NA",
      protLSID = "NA",
      assayLSID = "NA",
      panelLSID = "NA",
      QCcode = "NA",
      genotype_calls,
      check.names = FALSE,
      row.names = NULL
    )
    return(hapmap_row)

    # If the marker is an INDEL
  } else if (nchar(ref_allele) < nchar(alt_allele) | nchar(ref_allele) > nchar(alt_allele)) {
    if (verbose == TRUE) {
      print("INDEL")
    }
    genotype_calls[, 1:ncol(genotype_calls)] <- lapply(
      genotype_calls[, 1:ncol(genotype_calls), drop = FALSE],
      function(col) gsub(0, ref_allele, col)
    )
    genotype_calls[, 1:ncol(genotype_calls)] <- lapply(
      genotype_calls[, 1:ncol(genotype_calls), drop = FALSE],
      function(col) gsub(1, alt_allele, col)
    )

    if (nchar(ref_allele) > nchar(alt_allele)) {
      inserted_base <- substr(ref_allele, start = nchar(alt_allele) + 1, stop = nchar(alt_allele) + 1)

      hom_ref_str <- paste(ref_allele, ref_allele, sep = "/")
      genotype_calls[genotype_calls == hom_ref_str] <- inserted_base

      het_str_1 <- paste(ref_allele, alt_allele, sep = "/")
      genotype_calls[genotype_calls == het_str_1] <- 0

      het_str_2 <- paste(alt_allele, ref_allele, sep = "/")
      genotype_calls[genotype_calls == het_str_2] <- 0

      hom_alt_str <- paste(alt_allele, alt_allele, sep = "/")
      genotype_calls[genotype_calls == hom_alt_str] <- "-"

      genotype_calls[genotype_calls == "./."] <- "N"

      hapmap_row <- data.frame(
        `rs#` = marker_id,
        alleles = paste(inserted_base, "/-", sep = ""),
        chrom = chrom_name,
        pos = marker_pos,
        strand = "+",
        `assembly#` = "NA",
        center = "NA",
        protLSID = "NA",
        assayLSID = "NA",
        panelLSID = "NA",
        QCcode = "NA",
        genotype_calls,
        check.names = FALSE,
        row.names = NULL
      )
      return(hapmap_row)
    } else if (nchar(ref_allele) < nchar(alt_allele)) {
      inserted_base <- substr(alt_allele, start = nchar(ref_allele) + 1, stop = nchar(ref_allele) + 1)

      hom_alt_str <- paste(alt_allele, alt_allele, sep = "/")
      genotype_calls[genotype_calls == hom_alt_str] <- inserted_base

      het_str_1 <- paste(ref_allele, alt_allele, sep = "/")
      genotype_calls[genotype_calls == het_str_1] <- 0

      het_str_2 <- paste(alt_allele, ref_allele, sep = "/")
      genotype_calls[genotype_calls == het_str_2] <- 0

      hom_ref_str <- paste(ref_allele, ref_allele, sep = "/")
      genotype_calls[genotype_calls == hom_ref_str] <- "-"

      genotype_calls[genotype_calls == "./."] <- "N"

      hapmap_row <- data.frame(
        `rs#` = marker_id,
        alleles = paste("-/", inserted_base, sep = ""),
        chrom = chrom_name,
        pos = marker_pos,
        strand = "+",
        `assembly#` = "NA",
        center = "NA",
        protLSID = "NA",
        assayLSID = "NA",
        panelLSID = "NA",
        QCcode = "NA",
        genotype_calls,
        check.names = FALSE,
        row.names = NULL
      )
      return(hapmap_row)
    }

    # If the reference is a null allele
  } else if (ref_allele == ".") {
    if (verbose == TRUE) {
      print("DOMINANT")
    }
    genotype_calls[, 1:ncol(genotype_calls)] <- lapply(
      genotype_calls[, 1:ncol(genotype_calls), drop = FALSE],
      function(col) gsub(0, ref_allele, col)
    )
    genotype_calls[, 1:ncol(genotype_calls)] <- lapply(
      genotype_calls[, 1:ncol(genotype_calls), drop = FALSE],
      function(col) gsub(1, alt_allele, col)
    )
    genotype_calls[genotype_calls == "A/A"] <- "A"
    genotype_calls[genotype_calls == "T/T"] <- "T"
    genotype_calls[genotype_calls == "C/C"] <- "C"
    genotype_calls[genotype_calls == "G/G"] <- "G"
    genotype_calls[genotype_calls == "A/G"] <- "R"
    genotype_calls[genotype_calls == "G/A"] <- "R"
    genotype_calls[genotype_calls == "C/T"] <- "Y"
    genotype_calls[genotype_calls == "T/C"] <- "Y"
    genotype_calls[genotype_calls == "G/C"] <- "S"
    genotype_calls[genotype_calls == "C/G"] <- "S"
    genotype_calls[genotype_calls == "A/T"] <- "W"
    genotype_calls[genotype_calls == "T/A"] <- "W"
    genotype_calls[genotype_calls == "G/T"] <- "K"
    genotype_calls[genotype_calls == "T/G"] <- "K"
    genotype_calls[genotype_calls == "A/C"] <- "M"
    genotype_calls[genotype_calls == "C/A"] <- "M"
    genotype_calls[genotype_calls == "./."] <- "N"

    hapmap_row <- data.frame(
      `rs#` = marker_id,
      alleles = alt_allele,
      chrom = chrom_name,
      pos = marker_pos,
      strand = "+",
      `assembly#` = "NA",
      center = "NA",
      protLSID = "NA",
      assayLSID = "NA",
      panelLSID = "NA",
      QCcode = "NA",
      genotype_calls,
      check.names = FALSE,
      row.names = NULL
    )
    return(hapmap_row)

    # If the alternate allele is a null allele
  } else if (alt_allele == ".") {
    if (verbose == TRUE) {
      print("DOMINANT")
    }
    genotype_calls[, 1:ncol(genotype_calls)] <- lapply(
      genotype_calls[, 1:ncol(genotype_calls), drop = FALSE],
      function(col) gsub(0, ref_allele, col)
    )
    genotype_calls[, 1:ncol(genotype_calls)] <- lapply(
      genotype_calls[, 1:ncol(genotype_calls), drop = FALSE],
      function(col) gsub(1, alt_allele, col)
    )
    genotype_calls[genotype_calls == "A/A"] <- "A"
    genotype_calls[genotype_calls == "T/T"] <- "T"
    genotype_calls[genotype_calls == "C/C"] <- "C"
    genotype_calls[genotype_calls == "G/G"] <- "G"
    genotype_calls[genotype_calls == "A/G"] <- "R"
    genotype_calls[genotype_calls == "G/A"] <- "R"
    genotype_calls[genotype_calls == "C/T"] <- "Y"
    genotype_calls[genotype_calls == "T/C"] <- "Y"
    genotype_calls[genotype_calls == "G/C"] <- "S"
    genotype_calls[genotype_calls == "C/G"] <- "S"
    genotype_calls[genotype_calls == "A/T"] <- "W"
    genotype_calls[genotype_calls == "T/A"] <- "W"
    genotype_calls[genotype_calls == "G/T"] <- "K"
    genotype_calls[genotype_calls == "T/G"] <- "K"
    genotype_calls[genotype_calls == "A/C"] <- "M"
    genotype_calls[genotype_calls == "C/A"] <- "M"
    genotype_calls[genotype_calls == "./."] <- "N"

    hapmap_row <- data.frame(
      `rs#` = marker_id,
      alleles = ref_allele,
      chrom = chrom_name,
      pos = marker_pos,
      strand = "+",
      `assembly#` = "NA",
      center = "NA",
      protLSID = "NA",
      assayLSID = "NA",
      panelLSID = "NA",
      QCcode = "NA",
      genotype_calls,
      check.names = FALSE,
      row.names = NULL
    )
    return(hapmap_row)
  }
}
