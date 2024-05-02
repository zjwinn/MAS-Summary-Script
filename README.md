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

![An image of the haplotyping file header](https://github.com/zjwinn/MAS-Summary-Script/blob/master/Example_File_Head.png)

This file contains seven columns and some are required to run the code properly:

1. Locus (Manditory)
    - This is a unique identifyer that will follow your major locus designation. This must be unique from all other loci in the haplotyping file. There can be as many loci as desired, meaning that the haplotying file will have n rows based on n number of loci reported.
2. N_Expected (Optional)
    - This is a boolean column with the option of TRUE and FALSE. This column can be present or absent in the file. 
3. Markers (Manditory)
    - This is a vector of all markers (in order of haplotype) associated with the major locus of interest. All marker names are **case sensitive and must be present in the VCF**. To make this vector, users must take marker names and concatenate together using the operator ":". This allows the script to parse the marker names for use in haplotyping.
4. Pos (Manditory)
5. Het (Manditory)
6. Neg (Manditory)
7. Missing (Manditory)
    - Pos, Het, Neg, and Missing are vectors of all nucleotide sequences associated which indicate a "positive", "heterozygous", "negative", or "missing" call for a locus. For simplicity, we will only be looking at the Pos column in this following explination. In the case of the example provided with the script, "positive" in this case referes to a line possessing the *Fusarium* head blight resistance locus *Fhb1*. Users will take the [IUPAC nucleotide designation](https://en.wikipedia.org/wiki/Nucleic_acid_notation#cite_note-iupac1-1) to deliniate a haplotype associated with a "positive" call. To indicate a haplotype, users should delimit nucleotide designations with a ":" operator (e.g., 'C:T'). If there are multiple haplotypes associated with a "positive" call, the user can separate the haplotype calls with "_" to indicate a separate haplotype which indicates the same catigory (e.g., 'C:T_G:T'). If there is no known heterozygous call, then users may leave the "Het" column either blank or with a "NA". If a haplotype is identified by the script which does not appear in the provided haplotyping infromation, the call will come out as "UNIDENTIFIED".
   
All columns are case sensitive and columns labled as "manditory" must be found in the haplotyping file to get the code to run properly.

# Output files
There are four potential outputs from the script:

1. vcf_to_hapmap.hmp.gz
    - A hapmap file of selected markers for haplotyping. 
2. Marker_Report_Full.csv (comma seperated value file)
    - This file contains the markers used to make major locus haplotype calls, the resultant haplotype string associated with a major locus call, and the summary of the call (POS, HET, NEG, NA, UNIDENTIFIED).
3. Marker_Report_Summary.csv (comma seperated value file)
    - This file contains only the summary of haplotype calls for each major locus (POS, HET, NEG, NA, UNIDENTIFIED).
4. Unidentified_Haplotypes.csv (comma seperated value file)
    - This file will only be produced if there are "UNIDENTIFIED" haplotypes associated with a major locus. This file is meant to help individuals debugging their own haplotype formatted file.
  
# Identifying haplotypes associated with your locus
Suppose that there are some haplotypes that come up as "UNIDENTIFIED" but you are unsure of how to 
