# Marker Assisted Selection (MAS) Summary Script 

## Table of Contents

- [Introduction](#introduction)
- [Repository Structure](#repository-structure)
- [Running the Pipeline](#running-the-pipeline)
    - [Option A: Command-Line / Pipeline Execution (Recommended)](#option-a-command-line--pipeline-execution-recommended)
    - [Option B: Interactive RStudio / R Markdown Execution](#option-b-interactive-rstudio--r-markdown-execution)
- [Overall Workflow](#overall-workflow)
- [Sections of the script](#sections-of-the-script)
    - [1. Setting input files](#1-setting-input-files)
    - [2. Converting VCF to HapMap file](#2-converting-vcf-to-hapmap-file)
    - [3. Producing MAS report for major loci](#3-producing-mas-report-for-major-loci)
- [Input files](#input-files)
    - [1. VCF input file](#1-vcf-input-file)
    - [2. Haplotype key file](#2-haplotype-key-file)
- [Output files](#output-files)
- [Automatic haplotype calling: how does it work?](#automatic-haplotype-calling-how-does-it-work)
    - [1. Standard Codominant Biallelic Loci (SNPs and Indels)](#1-standard-codominant-biallelic-loci-snps-and-indels)
    - [2. Multiple Positive (`Pos`) Haplotypes in Auto Mode](#2-multiple-positive-pos-haplotypes-in-auto-mode)
    - [3. Dominant Markers, Deletions, and Translocations (`N` Alleles)](#3-dominant-markers-deletions-and-translocations-n-alleles)

## Introduction
This script was written for the express purpose of taking a variant calling format (VCF) file and a haplotyping file (tab delimited file) to produce a marker assisted selection (MAS) report. Below is a summary of the repository layout, the code, and its functions.

## Repository Structure

The repository is organized into a modular R project:

```text
MAS-Summary-Script/
├── LICENSE
├── MAS-Calling-Script.Rmd        # Main R Markdown driver script (interactive & parameterized)
├── MAS-Calling-Script.R          # Standalone CLI driver script for bash/SLURM pipelines
├── MAS-Summary-Script.Rproj      # RStudio Project file
├── README.md
├── R/                            # Modular helper functions sourced by the driver scripts
│   ├── count_header_lines.R      # Counts '##' metadata header lines in a VCF file
│   ├── vcf_to_hapmap.R           # Converts SNP, INDEL, and DOMINANT VCF rows to HapMap format
│   ├── auto_call_haplotypes.R    # Expands and classifies haplotypes when Auto == TRUE
│   ├── call_locus_haplotypes.R   # Calls locus haplotypes and annotations across all lines
│   └── compile_reports.R         # Merges parallel locus outputs, calculates stats, finds unidentified calls
└── example/                      # Example inputs and reference outputs
    ├── Example.vcf.gz
    ├── Example_File.txt
    ├── Example_File_Auto_Suggested_Haplotypes.csv
    ├── Example_File_Marker_Report_Full.csv
    ├── Example_File_Marker_Report_Statistics.csv
    ├── Example_File_Marker_Report_Summary.csv
    └── Example_File_hapmap.hmp.txt.gz
```

## Running the Pipeline

You can run the pipeline either from the **command line / bash scripts** (ideal for high-throughput pipelines) or interactively via **RStudio / R Markdown**.

### Option A: Command-Line / Pipeline Execution (Recommended)

Run `MAS-Calling-Script.R` directly using `Rscript`:

```bash
Rscript MAS-Calling-Script.R \
  --vcf path/to/genotypes.vcf.gz \
  --key path/to/haplotype_key.txt \
  --outdir path/to/output_dir \
  --cores 8
```

#### CLI Options:
- `-v, --vcf FILE`: Input VCF file (`.vcf` or `.vcf.gz`) [defaults to `example/Example.vcf.gz`]
- `-k, --key FILE`: Tab-delimited MAS haplotype key file [defaults to `example/Example_File.txt`]
- `-o, --outdir DIR`: Directory where report outputs will be written [defaults to `example/`]
- `-c, --cores INT`: Number of parallel CPU cores to use [default: `detectCores() - 1`]
- `--verbose [BOOL]`: Enable verbose logging (TRUE/FALSE) [default: `TRUE`]
- `-h, --help`: Display usage summary and exit

### Option B: Interactive RStudio / R Markdown Execution

1. Open `MAS-Summary-Script.Rproj` in RStudio.
2. Open `MAS-Calling-Script.Rmd`.
3. Set your input file paths in lines 60–85, or pass them dynamically when calling `rmarkdown::render()`:
   ```R
   rmarkdown::render(
     "MAS-Calling-Script.Rmd",
     params = list(
       vcf_file = "path/to/genotypes.vcf.gz",
       mas_key_file = "path/to/haplotype_key.txt",
       output_directory = "path/to/output_dir",
       cores = 8
     )
   )
   ```
4. Run all chunks or click **Knit** to execute the pipeline and generate an HTML summary.


## Overall Workflow
The script in this repository follows the proposed workflow:

1. Read in a VCF with known informative markers (KIMs)
2. Use a haplotyping file ([format explained in later section](#2-haplotype-key-file)) to identify haplotypes of interest for specific major loci
3. Report a detailed summary of major loci, associated KIMs, and produce a debugging dataframe to identify non-specified haplotypes for major loci

## Sections of the script
When opened through `MAS-Summary-Script.Rproj`, `MAS-Calling-Script.Rmd` automatically loads required packages and sources all helper functions in `./R/`. There are three main sections of the script:

1. Setting input files
2. Converting VCF to HapMap file
3. Producing MAS report for major loci

### 1. Setting input files

In `MAS-Calling-Script.Rmd` there are three lines that need to be changed by the user. These lines are found at lines 52, 55, and 80 (as visualized using [RStudio](https://posit.co/)). Essentially these lines must be changed to match the user inputs (`vcf_file`, `mas_key_file`, and `output_directory`). Other than these three changes, no other part of the script needs to be edited. Both steps 2 and 3 of this pipeline are automated and require no user input. The most important inputs from users are the input files, which are discussed in the [Input files](#input-files) section.

### 2. Converting VCF to HapMap file

In this section of `MAS-Calling-Script.Rmd` (lines 99–171), the provided VCF is subset to the markers listed in the key file, converted in parallel to HapMap format via `vcf_to_hapmap()` (supporting both unphased `/` and phased `|` genotypes across SNPs, indels, and dominant markers), and written to the user-specified output directory.

### 3. Producing MAS report for major loci

In this section of `MAS-Calling-Script.Rmd` (lines 177–293), the HapMap table created in step two is used in tandem with the haplotyping key file to call and summarize the allelic states of each line in parallel via `call_locus_haplotypes()`, `auto_call_haplotypes()`, and the reporting helpers in `R/compile_reports.R`.

## Input files
There are two input files for this script:

1. A variant calling format (VCF) genotyping file
2. A uniquely formatted haplotyping key file (tab delimited file)

### 1. VCF input file

Variant calling format files have specific and rigorous standards. To review proper VCF format [see this link](https://samtools.github.io/hts-specs/VCFv4.2.pdf). A haplotype "call" will be provided for every individual in this VCF. 

### 2. Haplotype key file

The custom haplotyping key file has a unique structure that allows for a single input file to make haplotype "calls". This file contains twelve columns and all are required to run the code properly:

1. Locus
    - This is a unique identifier that will follow your major locus designation. This must be unique from all other loci in the haplotyping file. There can be as many loci as desired, meaning that the haplotyping file will have n rows based on n number of loci reported.
2. Auto
    - This is a boolean column with the option of `TRUE` and `FALSE`. This option indicates if the user wishes to use the auto function for locus calling ([auto is explained in a later section](#automatic-haplotype-calling-how-does-it-work)). This column must be present in the file and the options must be either `TRUE` or `FALSE`. 
3. N_Missing
    - This is a numeric argument which may be interpreted as the following: "If I am using the auto function to call my loci, how many missing loci am I willing to tolerate". For example, if the user is calling a major locus using four markers, they may wish to still provide a call for individuals missing data for one marker. Thus `N_Missing` would equal `1` in this case.
4. Markers
    - This is a vector of all markers (in order of haplotype) associated with the major locus of interest. All marker names are **case sensitive and must be present in the VCF**. To make this vector, users must take marker names and concatenate together using the operator `":"`. This allows the script to parse the marker names for use in haplotyping.
5. Pos 
6. Het 
7. Neg 
8. Missing
    - `Pos`, `Het`, `Neg`, and `Missing` are vectors of all nucleotide sequences which indicate a "positive", "heterozygous", "negative", or "missing" call for a locus. For simplicity, we will only be looking at the `Pos` column in this following explanation. In the case of the example provided with the script, "positive" refers to a line possessing the *Fusarium* head blight resistance locus *Fhb1*. Users provide the [IUPAC nucleotide designation](https://en.wikipedia.org/wiki/Nucleic_acid_notation#cite_note-iupac1-1) (or `-` / `0` for indels, and `N` for null/missing calls) to delineate a haplotype associated with a "positive" call. To indicate a haplotype, users delimit nucleotide designations with a `":"` operator (e.g., `'C:T'`). If there are multiple haplotypes associated with a "positive" call, the user can separate the haplotype calls with `"_"` to indicate a separate haplotype which indicates the same category (e.g., `'C:T_G:T'`). If there is no known heterozygous call, then users may leave the `"Het"` column either blank or with `"NA"`. If a haplotype is identified by the script which does not appear in the provided haplotyping information, the call will come out as `"UNIDENTIFIED"`, unless `Auto` is indicated as `TRUE` (see section below for description of the auto function).
9. Pos_Annot
10. Het_Annot
11. Neg_Annot
12. Missing_Annot
    - `Pos`, `Het`, `Neg`, and `Missing` `"_Annot"` are annotations the user wishes to assign to the positive, heterozygous, negative, and missing cases. These annotations may be left blank if the user wants the default assignments (e.g., `POS`, `HET`, `NEG`, and `NA`); however, this may not be appropriate for all cases. For resistance loci, like *Fhb1*, perhaps the default assignments make sense. However, a locus where the delineation between states is not straightforward (e.g., the 3A color locus for wheat which delineates white vs. red seed color) may require further explanation. These columns allow for `POS`, `HET`, `NEG`, and `NA` to be reported in a more intuitive way (e.g., `'Red'`, `'Pink'`, `'White'`, and `'Did not amplify'`).   

All columns are case sensitive and all columns must be present in the haplotype key file to get this code to run properly. [An example of a 'haplotype key file' is provided along with this script.](https://github.com/zjwinn/MAS-Summary-Script/blob/master/example/Example_File.txt)

## Output files

The R script functions by taking **the name of your haplotyping key file and assigning that as the string prior to a file name**. In the list below where you see **i**, this stands in place of the name of your haplotyping key file:

1. **i**_hapmap.hmp.txt.gz (compressed hapmap file)
    - A HapMap file ([here is an explanation of what HapMap format is](https://statgen-esalq.github.io/Hapmap-and-VCF-formats-and-its-integration-with-onemap/#:~:text=The%20Hapmap%20file%20format%20is%20a%20table%20which,Hapmap%20file%20by%20chromosome%20or%20a%20general%20file.)) of the selected markers for haplotyping indicated by the user's haplotyping key input file. 
2. **i**_Marker_Report_Full.csv (comma separated value file)
    - This file contains the markers used to make major locus haplotype calls, the resultant haplotype string associated with a major locus call, the summary of the call (`POS`, `HET`, `NEG`, `NA`, `UNIDENTIFIED`), and the annotated summary of the call.
3. **i**_Marker_Report_Summary.csv (comma separated value file)
    - This file contains only the annotated summary of haplotype calls for each major locus.
4. **i**_Unidentified_Haplotypes.csv (comma separated value file)
    - This file will only be produced if there are `"UNIDENTIFIED"` haplotypes associated with a major locus. This file is meant to help individuals debug their own haplotype formatted file.
5. **i**_Auto_Suggested_Haplotypes.csv (comma separated value file)
    - If the user indicates that any locus must be automatically called ([see automatic haplotype calling section for further explanation](#automatic-haplotype-calling-how-does-it-work)), then this file will be produced. This file contains all the automatic calls made and all suggested haplotype classifications. 
6. **i**_Marker_Report_Statistics.csv (comma separated value file)
    - This output shows the counts and frequencies of haplotypes assigned for a major locus.

## Automatic haplotype calling: how does it work?

As the number of markers $M$ defining a locus increases, the number of potential haplotype combinations ($4^M$, including homozygotes, heterozygotes, and missing calls) grows exponentially. Here is an example in R showing how a 3-marker locus produces 64 possible haplotypes:

```R
marker_grid <- data.frame(
  m1 = c("A", "R", "G", "N"),
  m2 = c("A", "R", "G", "N"),
  m3 = c("A", "R", "G", "N")
)

marker_grid <- expand.grid(marker_grid)
marker_grid$haplotype <- paste(marker_grid$m1, marker_grid$m2, marker_grid$m3, sep = ":")

print(length(marker_grid$haplotype)) # 64
```

With 4 or 5 markers, enumerating every possible haplotype by hand quickly becomes impractical. The `Auto = TRUE` algorithm (`R/auto_call_haplotypes.R`) automates this expansion across three common genotyping scenarios:

### 1. Standard Codominant Biallelic Loci (SNPs and Indels)

For a region characterized by codominant biallelic markers (SNPs or indels where `-` is the gap allele and `0` is the heterozygous indel state), the algorithm makes the following core assumptions:

1. **Seed Positive Case:** A complete-case homozygous haplotype defines the baseline positive state (e.g., `A:G:T:G:A:C:C`).
2. **Tolerated Missingness (`N_Missing`):** Haplotypes matching the positive allele at all non-missing positions with at most `N_Missing` `"N"` calls are still called `POS` (e.g., `A:G:N:G:A:C:C` when `N_Missing >= 1`).
3. **Heterozygous Recombination Risk:** Any haplotype composed of positive, heterozygous (IUPAC `R, Y, S, W, K, M` or indel `0`), and tolerated missing calls—with at least one heterozygous marker and no negative alleles—is called `HET` (e.g., `A:G:Y:W:A:C:C`).
4. **Conservative Negative Classification:** A single marker carrying the alternative (negative) allele could indicate recombination separating the marker from the causal polymorphism; therefore, any haplotype with $\le \text{N\_Missing}$ missing markers that carries at least one negative allele is called `NEG` (e.g., `A:G:C:G:A:C:C`).
5. **Excess Missing Data:** Any haplotype with more than `N_Missing` missing calls at codominant markers is called `MISSING` (e.g., `A:G:N:N:N:N:C`).

### 2. Multiple Positive (`Pos`) Haplotypes in Auto Mode

In practice, a user may set `Auto = TRUE` on a locus while providing multiple positive haplotypes separated by `"_"` in the `Pos` column (for example, `C:T_G:T_N:T` when more than one background haplotype is linked to the favorable allele). When multiple `Pos` haplotypes are supplied:

1. The script emits a warning noting that multiple positive haplotypes were provided for that locus.
2. It selects the **first complete homozygous positive haplotype** (e.g., `C:T`) as the seed to generate the automatic `POS`, `HET`, `NEG`, and `MISSING` haplotype combinations.
3. **User Override:** All additional haplotypes explicitly provided by the user in `Pos`, `Het`, `Neg`, or `Missing` override the automatically generated calls so that no manually curated positive haplotype is lost.

### 3. Dominant Markers, Deletions, and Translocations (`N` Alleles)

In targeted amplicon sequencing (e.g., AgriSeq / GBS) and KASP assays, markers frequently behave as **dominant (presence/absence)** markers or span **deletions and alien translocations** where `"N"` (lack of amplification) is an informative biological state rather than random missing data. `auto_call_haplotypes()` handles three dominant/null cases:

1. **Amplified-Positive Dominant Markers (e.g., `Pos = A`, `Neg = N`):**
   - When a marker in the VCF is monomorphic/single-allele (`REF = "."` or `ALT = "."`) or is an amplicon presence/absence contig (`A/-` or `-/A` where `-` is never observed in the population), the amplified base is the positive allele and `"N"` is the negative (null) allele.
2. **Partial Null-Positive / Deletion Haplotypes (e.g., `Pos = A:N`, `N:A`, or `A:N:N:G`):**
   - Many loci combine an amplicon on the insertion/translocation with an amplicon on the wild-type sequence (for example, `Tamyb10-A1` with `Pos = A:N`, `Vrn-B1a/b` with `Pos = N:A`, or `Vrn-D1` with `Pos = G:N`).
   - When the positive haplotype contains `"N"` at marker $m$, `auto_call_haplotypes()` treats `"N"` at marker $m$ as a **null-positive (deletion) state**, assigns the amplified nucleotide from the VCF as both the `Negative` and `Heterozygous` state at marker $m$ (since a heterozygote carries one copy of the present sequence and therefore amplifies), and does not penalize `"N"` at dominant markers toward `N_Missing` (unless every marker in the locus is `"N"`).
   - For example, given only `Pos = A:N` for a two-marker dominant locus where Marker 1 amplifies `A` and Marker 2 amplifies `T`, `Auto = TRUE` automatically derives:
     - `Pos = A:N`
     - `Het = A:T` (both amplicons present in the heterozygote)
     - `Neg = N:T`
     - `Missing = N:N`
3. **Full-Locus Deletions / Translocation Absences (`Pos = N:N:N:N`):**
   - If the positive haplotype is defined entirely by null calls across all markers in the locus (`N:N:...:N`), `auto_call_haplotypes()` preserves `N:N:...:N` as `POS` (rather than overwriting it as `MISSING`) and classifies haplotypes that amplify the wild-type alleles as `NEG`.
