# Count lines that start with ## in a VCF file
count_header_lines <- function(x) {
  # Open connection to file
  con <- file(x, "r")
  # Set the number for first line
  header_lines <- 0

  # While
  while (TRUE) {
    # Read in the line = 1
    line <- readLines(con, n = 1)

    # If the line equal to zero
    if (length(line) == 0) break # Exit loop when end of file is reached

    # If the line contains ##
    if (substr(line, 1, 2) == "##") {
      # Add to length
      header_lines <- header_lines + 1

      # Else if the line does not contain ##
    } else {
      # Exit loop when non-header lines start
      break
    }
  }

  # Close the file you are reading
  close(con)

  # Return the number of header lines
  return(header_lines)
}
