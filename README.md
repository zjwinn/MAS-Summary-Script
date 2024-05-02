# Introduction
This script was written for the express purpose of taking a variant calling format (VCF) file and a haplotyping file (tab delimited file) to produce a marker assisted selection (MAS) report. Below is a summary of the code and its functions.

# Workflow
The script in this repository follows the proposed workflow:

1. Read in a VCF with known informative markers (KIMs).
2. Use a haplotyping file (format explained in later secitons) to identify haplotypes of interest for specific major loci.
3. Report a detailed summary of major loci, associated KIMs, and produce a debugging dataframe to identify non-specified haplotypes for major loci.

# Sections of the script
There are three distinct sections of the script:

1. 
