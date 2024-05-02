# Introduction
This script was written for the express purpose of taking a variant calling format (VCF) file and a haplotyping file (tab delimited file) to produce a marker assisted selection (MAS) report. Below is a summary of the code and its functions.

# Workflow
The script in this repository follows the proposed workflow:

1. Read in a VCF with known informative markers (KIMs)
2. Use a haplotyping file (format explained in later secitons) to identify haplotypes of interest for specific major loci
3. Report a detailed summary of major loci, associated KIMs, and produce a debugging dataframe to identify non-specified haplotypes for major loci

# Sections of the script
There are three distinct sections of the script:

1. Setting input files
2. Converting VCF to hap-map like file
3. Producing MAS report for major loci

Both steps 2 and 3 of this pipeline are automated and require no user input. The most important inputs from users are the input files. Below we will discuss input formats.

# Input files
There are two input files for this script:

1. A variant calling format (VCF) genotyping file
2. A uniquly formatted haplotyping file (tab delmimited file)

Variant calling format files have specific and rigorous standards. To review proper VCF format [see this link](https://samtools.github.io/hts-specs/VCFv4.2.pdf). The custom haplotype formatted file has a unique structure that allows for a single input file to make haplotype "calls". Below is a displayed example of a haplotyping file:



