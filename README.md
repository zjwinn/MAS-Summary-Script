# Marker Assisted Selection (MAS) Summary Script 

## Table of Contents

- [Introduction](#Introduction)
- [Overall Workflow](#Workflow)
- [Sections of the script](#Sections-of-the-script)
    - [1. Setting input files](#1-Setting-input-files)
    - [2. Converting VCF to HapMap file](#2-Converting-VCF-to-HapMap-file)
    - [3. Producing MAS report for major loci](#3-Producing-MAS-report-for-major-loci)
- [Input files](#Input-files)
    - [1. VCF input file](#1-VCF-input-file)
    - [2. Haplotype key file](#2-Haplotype-key-file)
- [Output files](#Output-files)
- [Automatic haplotype calling: how does it work?](#automatic-haplotype-calling-how-does-it-work)

## Introduction
This script was written for the express purpose of taking a variant calling format (VCF) file and a haplotyping file (tab delimited file) to produce a marker assisted selection (MAS) report. Below is a summary of the code and its functions.

## Overall Workflow
The script in this repository follows the proposed workflow:

1. Read in a VCF with known informative markers (KIMs)
2. Use a haplotyping file ([format explained in later section](#2-haplotype-key-file)) to identify haplotype of interest for specific major loci
3. Report a detailed summary of major loci, associated KIMs, and produce a debugging dataframe to identify non-specified haplotypes for major loci

## Sections of the script
There are three distinct sections of the script:

1. Setting input files
2. Converting VCF to HapMap file
3. Producing MAS report for major loci

### 1. Setting input files

In the R script there are three lines that need to be changed by the user. These lines are found at line 44, 47, and 84 (as visualized using [R studio](https://posit.co/)). Essentially these lines must be changed to match the user inputs. Other than these three changes, no other part of the script should be edited. Both steps 2 and 3 of this pipeline are automated and require no user input. The most important inputs from users are the input files. We discuss the input files in a [later section](#Input-files).

### 2. Converting VCF to HapMap file

In this section of the R script from lines 108-477, the provided VCF is converted to a HapMap and the output is then written to the user specified out directory. 

### 3. Producing MAS report for major loci

In this section of the R script from lines 485-1335, the provided HapMap created in step two is then used in tandem with the haplotyping key file to create summaries of the allelic states of each line in the provided VCF. This summary then creates several outputs for users that they may use as they see fit.  

## Input files
There are two input files for this script:

1. A variant calling format (VCF) genotyping file
2. A unique formatted haplotyping key file (tab delimited file)

### 1. VCF input file

Variant calling format files have specific and rigorous standards. To review proper VCF format [see this link](https://samtools.github.io/hts-specs/VCFv4.2.pdf). A haplotype "call" will be provided for every individual in this VCF. 

### 2. Haplotype key file

The custom haplotyping key file has a unique structure that allows for a single input file to make haplotype "calls". This file contains twelve columns and all are required to run the code properly:

1. Locus
    - This is a unique identifier that will follow your major locus designation. This must be unique from all other loci in the haplotyping file. There can be as many loci as desired, meaning that the haplotyping file will have n rows based on n number of loci reported.
2. Auto
    - This is a boolean column with the option of TRUE and FALSE. This option indicates if the user wishes to use the auto function for locus calling ([auto is explained in a later section](#automatic-haplotype-calling-how-does-it-work)). This column can must be present in the file and the options must be either TRUE or FALSE. 
3. N_Missing
    - This is a numeric argument which may be interpreted as the following: "If I am using the auto function to call my loci, how many missing loci am I willing to tolerate". For example, if the user is calling a major locus using four markers, they may wish to still provide a call individuals missing data for one marker. Thus N_Missing would equal one in this case. If this column is left blank, then it will be considered zero.  
4. Markers
    - This is a vector of all markers (in order of haplotype) associated with the major locus of interest. All marker names are **case sensitive and must be present in the VCF**. To make this vector, users must take marker names and concatenate together using the operator ":". This allows the script to parse the marker names for use in haplotyping.
5. Pos 
6. Het 
7. Neg 
8. Missing
    - Pos, Het, Neg, and Missing are vectors of all nucleotide sequences associated which indicate a "positive", "heterozygous", "negative", or "missing" call for a locus. For simplicity, we will only be looking at the Pos column in this following explanation. In the case of the example provided with the script, "positive" in this case refers to a line possessing the *Fusarium* head blight resistance locus *Fhb1*. Users will take the [IUPAC nucleotide designation](https://en.wikipedia.org/wiki/Nucleic_acid_notation#cite_note-iupac1-1) to delineate a haplotype associated with a "positive" call. To indicate a haplotype, users should delimit nucleotide designations with a ":" operator (e.g., 'C:T'). If there are multiple haplotypes associated with a "positive" call, the user can separate the haplotype calls with "_" to indicate a separate haplotype which indicates the same category (e.g., 'C:T_G:T'). If there is no known heterozygous call, then users may leave the "Het" column either blank or with a "NA". If a haplotype is identified by the script which does not appear in the provided haplotyping information, the call will come out as "UNIDENTIFIED", unless auto is indicated as TRUE (see section below for description of auto function).
9. Pos_Annot
10. Het_Annot
11. Neg_Annot
12. Missing_Annot
    - Pos, Het, Neg, and Missing "_Annot" are annotation the user wishes to assign to the positive, heterozygous, negative, and missing cases. These annotations may be left blank if the user wants the default assignments (e.g., POS, HET, NEG, and NA) however this may not be appropriate for all cases. For resistance loci, like *Fhb1*, perhaps the default assignments make sense. However, a locus where the delineation between states is not straight forward (e.g., The 3A color locus for wheat which delineates white vs. red seed color) may require further explanation. These columns allow for POS, HET, NEG, and NA to be reported in a more intuitive way (e.g, 'Red', 'Pink', 'White', and 'Did not amplify').   

All columns are case sensitive and all columns must be present in the haplotype key file to get this code to run properly. [An example of a 'haplotype key file' is provided along with this script.](https://github.com/zjwinn/MAS-Summary-Script/blob/master/Example_File.txt)

## Output files

The R script functions by taking **the name of your haplotyping key file and assigning that as the string prior to a file name**. In the below script where you see **i**, this will stand in place of the name of your haplotyping key file:

1. **i**_hapmap.hmp.gz (compressed hapmap file)
    - A HapMap file ([here is an explanation of what HapMap format is](https://statgen-esalq.github.io/Hapmap-and-VCF-formats-and-its-integration-with-onemap/#:~:text=The%20Hapmap%20file%20format%20is%20a%20table%20which,Hapmap%20file%20by%20chromosome%20or%20a%20general%20file.)) of the selected markers for haplotyping indicated by the user's haplotyping key input file. 
2. **i**_Marker_Report_Full.csv (comma separated value file)
    - This file contains the markers used to make major locus haplotype calls, the resultant haplotype string associated with a major locus call, the summary of the call (POS, HET, NEG, NA, UNIDENTIFIED), and the annotated summary of the call.
3. **i**_Marker_Report_Summary.csv (comma separated value file)
    - This file contains only the annotated summary of haplotype calls for each major locus.
4. **i**_Unidentified_Haplotypes.csv (comma separated value file)
    - This file will only be produced if there are "UNIDENTIFIED" haplotypes associated with a major locus. This file is meant to help individuals debugging their own haplotype formatted file.
5. **i**_Auto_Suggested_Haplotypes.csv (comma separated value file)
    - If the user indicates that any locus must be automatically called ([see automatic haplotype calling section for further explanation](#automatic-haplotype-calling-how-does-it-work)), then this file will be produced. This file contains all the automatic calls made and all suggested haplotype classifications. 
6. **i**_Marker_Report_Statistics.csv (comma separated value file)
    - This output shows the counts and frequencies of haplotypes assigned for a major locus.

## Automatic haplotype calling: how does it work?

A custom algorithm for calling biallelic loci was derived for this exact purpose in the provided R code. As N number of loci goes to infinity so too does the number of potential combinations. In fact the number of potential haplotypes at a biallelic locus (including missing data) can become gigantic quickly!    

Here is an example of some R code that can show you how many combinations you get with defining a haplotype with only three markers:  

```R
temp1 <- data.frame(m1=c("A","R","G","N"),                    
                    m2=c("A","R","G","N"),
                    m3=c("A","R","G","N"))

temp1 <- expand.grid(temp1)

temp1$haplotype <- paste(temp1$m1, temp1$m2, temp1$m3, sep = ":")

print(length(temp1$haplotype))
```

As you can see, if you only have three markers in a region to define a gene or quantitative trait locus (QTL) the number of potential haplotypes inflates to 64. This becomes worse with four, five, or more markers to define a locus. Clearly, no one has the time to notate every possible haplotype for a QTL containing three or more markers, so there must be some way to handle that. 

In the case of any larger region that a user is characterizing, there are several assumptions we have made with our algorithm:

1. The only true positive case is a complete case haplotype with no missing data.
    - E.G., A:G:T:G:A:C:C
2. There is some missing level of data that the user is willing to tolerate to make a positive call.
    - E.G., A:G:N:G:A:C:C is still called positive even with a single missing marker at position three.
3. A single heterozygous marker in a region may cause recombination, which means that the major locus must be called as heterozygous if it has some combination of positive and heterozygous markers. 
    - E.G., A:G:Y:W:A:C:C is called heterozygous even though all other loci are in the true positive conformation because positions three and four are in the heterozygous state.
4. A single marker associated with the alternative allele of the positive case within a region could indicate that recombination has occurred and the individual may not possess the causal polymorphism which confers the trait of interest. Therefore, to take a conservative approach to classification, we call that individual negative.
    - E.G., A:G:C:G:A:C:C is called negative even though most markers are in the positive conformation because the third position is C instead of T and C is associated with the negative case.
5. If there is some level of missing data greater that the threshold defined by the user, then that call is considered missing.
    - E.G., A:G:N:N:N:N:C has four missing markers which may be higher than a user is willing to tolerate.

In the case of the algorithm we have proposed, the user only needs to define **a single complete-case positive example in the 'Pos' column**. This means to use the algorithm we have proposed above that the user must provide a sequence free of heterozygous loci and missing data. Here are some unacceptable examples:

- N:G:G
- Y:G:G
- G:G:G:T:G:N:C:C

In the case of setting a locus to Auto==TRUE, the user must provide this complete-case positive example and set the tolerable number of missing data in the N_Missing column. In simplest terms, the auto function does the following with this information provided:

1. Identify the true positive case
2. Identify the true negative case
3. Impute the true heterozygous case
4. Intuit the number of missing markers possible
5. Create every possible haplotype at the locus
6. Identify which haplotypes have too many missing markers and call them missing
7. Identify which haplotype are positive cases that have a tolerable number of missing markers and call them positive
8. Identify which heterozygous cases contain only missing, positive, and heterozygous marker and call them heterozygous
9. Identify all cases which contain at least one negative case marker and call them negative

Once this process has been completed, the output is then used to delineate what haplotype an individual has and assign it to them. 
