# Marker Assisted Selection (MAS) Summary Script 

## Table of Contents

- [Introduction](#Introduction)
- [Overall Workflow](#Workflow)
- [Sections of the script](#Sections-of-the-script)
    - [1. Setting input files](#1.-Setting-input-files)
    - [2. Converting VCF to HapMap file](#2.-Converting-VCF-to-HapMap-file)
    - [3. Producing MAS report for major loci](#3.-Producing-MAS-report-for-major-loci)
- [Input files](#Input-files)
    - [1. VCF input file](#1.-VCF-input-file)
    - [2. Haplotyping key file](#2.-Haplotyping-key-file)
- [Output files](#Output-files)
- [Automatic haplotype calling: how does it work?](#Automatic-haplotype-calling:-how-does-it-work)

## Introduction
This script was written for the express purpose of taking a variant calling format (VCF) file and a haplotyping file (tab delimited file) to produce a marker assisted selection (MAS) report. Below is a summary of the code and its functions.

## Overall Workflow
The script in this repository follows the proposed workflow:

1. Read in a VCF with known informative markers (KIMs)
2. Use a haplotyping file (format explained in later secitons) to identify haplotypes of interest for specific major loci
3. Report a detailed summary of major loci, associated KIMs, and produce a debugging dataframe to identify non-specified haplotypes for major loci

## Sections of the script
There are three distinct sections of the script:

1. Setting input files
2. Converting VCF to HapMap file
3. Producing MAS report for major loci

### 1. Setting input files

In the R script there are three lines that need to be changed by the user. These lines are found at line 44, 47, and 84 (as visualized using [R studio](https://posit.co/)). Esentially these lines must be changed to match the user inputs. Other than these three changes, no other part of the script should be edited. Both steps 2 and 3 of this pipeline are automated and require no user input. The most important inputs from users are the input files. We discuss the input files in a [later section](#Input-files).

### 2. Converting VCF to HapMap file

In this section of the R script from lines 108-477, the provided VCF is converted to a HapMap and the output is then written to the user specificed out directory. 

### 3. Producing MAS report for major loci

In this section of the R script from lines 485-1335, the provided HapMap created in step two is then used in tandem with the haplotyping key file to create summaries of the allelic states of each line in the provided VCF. This summary then creates several outputs for users that they may use as they see fit.  

## Input files
There are two input files for this script:

1. A variant calling format (VCF) genotyping file
2. A uniquly formatted haplotyping key file (tab delmimited file)

### 1. VCF input file

Variant calling format files have specific and rigorous standards. To review proper VCF format [see this link](https://samtools.github.io/hts-specs/VCFv4.2.pdf). A haplotype "call" will be provided for every individual in this VCF. 

### 2. Haplotyping key file

The custom haplotyping key file has a unique structure that allows for a single input file to make haplotype "calls". This file contains twelve columns and all are required to run the code properly:

1. Locus
    - This is a unique identifyer that will follow your major locus designation. This must be unique from all other loci in the haplotyping file. There can be as many loci as desired, meaning that the haplotying file will have n rows based on n number of loci reported.
2. Auto
    - This is a boolean column with the option of TRUE and FALSE. This option indicates if the user wishes to use the auto function for locus calling ([auto is explained in a later section](#Automatic-haplotype-calling:-how-does-it-work?)). This column can must be present in the file and the options must be either TRUE or FALSE. 
3. N_Missing
    - This is a numeric argument which may be interpreted as the following: "If I am using the auto function to call my loci, how many missing loci am I willing to tollerate". For example, if the user is calling a major locus using four markers, they may wish to still provide a call individuals missing data for one marker. Thus N_Missing would equal one in this case. If this collumn is left blank, then it will be considered zero.  
4. Markers
    - This is a vector of all markers (in order of haplotype) associated with the major locus of interest. All marker names are **case sensitive and must be present in the VCF**. To make this vector, users must take marker names and concatenate together using the operator ":". This allows the script to parse the marker names for use in haplotyping.
5. Pos 
6. Het 
7. Neg 
8. Missing
    - Pos, Het, Neg, and Missing are vectors of all nucleotide sequences associated which indicate a "positive", "heterozygous", "negative", or "missing" call for a locus. For simplicity, we will only be looking at the Pos column in this following explination. In the case of the example provided with the script, "positive" in this case referes to a line possessing the *Fusarium* head blight resistance locus *Fhb1*. Users will take the [IUPAC nucleotide designation](https://en.wikipedia.org/wiki/Nucleic_acid_notation#cite_note-iupac1-1) to deliniate a haplotype associated with a "positive" call. To indicate a haplotype, users should delimit nucleotide designations with a ":" operator (e.g., 'C:T'). If there are multiple haplotypes associated with a "positive" call, the user can separate the haplotype calls with "_" to indicate a separate haplotype which indicates the same catigory (e.g., 'C:T_G:T'). If there is no known heterozygous call, then users may leave the "Het" column either blank or with a "NA". If a haplotype is identified by the script which does not appear in the provided haplotyping infromation, the call will come out as "UNIDENTIFIED", unless auto is indicated as TRUE (see section below for description of auto function).
9. Pos_Annot
10. Het_Annot
11. Neg_Annot
12. Missing_Annot
    - Pos, Het, Neg, and Missing "_Annot" are annotation the user wishes to assign to the positive, heterozygous, negative, and missing cases. These annotations may be left blank if the user wants the default assignments (e.g., POS, HET, NEG, and NA) however this may not be appropriate for all cases. For resistance loci, like *Fhb1*, perhaps the default assignments make sense. However, a locus where the deliniation between states is not straight forward (e.g., The 3A color locus for wheat which deliniates white vs. red seed color) may require further explanation. These columns allow for POS, HET, NEG, and NA to be reported in a more intuitive way (e.g, 'Red', 'Pink', 'White', and 'Did not amplify').   

All columns are case sensitive and all columns must be present in the haplotype key file to get this code to run properly. [An example file of a 'haplotype key file' is provided along with this script.](https://github.com/zjwinn/MAS-Summary-Script/blob/master/Example_File.txt)

## Output files

The R script functions by taking **the name of your haplotyping key file and assigning that as the string prior to a file name**. In the below script where you see **i**, this will stand in place of the name of your haplotyping key file:

1. **i**_hapmap.hmp.gz (compressed hapmap file)
    - A HapMap file ([here is an example of a HapMap format](https://statgen-esalq.github.io/Hapmap-and-VCF-formats-and-its-integration-with-onemap/#:~:text=The%20Hapmap%20file%20format%20is%20a%20table%20which,Hapmap%20file%20by%20chromosome%20or%20a%20general%20file.)) of the selected markers for haplotyping indicated by the user's haplotyping key input file. 
2. **i**_Marker_Report_Full.csv (comma seperated value file)
    - This file contains the markers used to make major locus haplotype calls, the resultant haplotype string associated with a major locus call, the summary of the call (POS, HET, NEG, NA, UNIDENTIFIED), and the annotated summary of the call.
3. **i**_Marker_Report_Summary.csv (comma seperated value file)
    - This file contains only the annotated summary of haplotype calls for each major locus.
4. **i**_Unidentified_Haplotypes.csv (comma seperated value file)
    - This file will only be produced if there are "UNIDENTIFIED" haplotypes associated with a major locus. This file is meant to help individuals debugging their own haplotype formatted file.
5. **i**_Auto_Suggested_Haplotypes.csv (comma seperated value file)
    - If the user indicates that any locus must be automatically called ([see automatic haplotype calling section for further explanation](#Automatic-haplotype-calling:-how-does-it-work?)), then this file will be produced. This file contains all the automatic calls made and all suggested haplotype classifications. 
6. **i**_Marker_Report_Statistics.csv (comma seperated value file)
    - This output shows the counts and frequencies of haplotypes assigned for a major locus.

## Automatic haplotype calling: how does it work?

