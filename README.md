# RNA-Seq analysis of _Fusobacterium nucleatum_ upon PNA treatment

- Project name: An antisense antibiotic candidate with unexpected bactericidal activity against several different 
                fusobacteria species 

- Experiments: Valentina Cosi, Linda Popella, Chandradhish Ghosh

- Supervision: Jörg Vogel, Lars Barquist

- Data analysis/RNA-Seq analysis:  Jakob J. Jung

- Start: 2022 

  

## Introduction

RNA-Seq analysis of transcriptomes of Fusobacterium nucleatum nucleatum and Fusobacterium nucleatum vincentii treated 
with the peptide nucleic acid (PNA) PNA-79, scrambed PNA, and untreated control. 2 time points (30 min and 2h) were 
analyzed, and experiments were done in triplicate, so in total 36 samples were analyzed in total.

## Directory structure

The project is divided in 3 main directories:

-   [data](data) : contains all raw, intermediate and final data of the project.  
-   [analysis](analysis): contains analysis files, such as figures, plots, tables, etc. 
-   [scripts](scripts): contains scripts to process and analyze data from the data directory.

Some directories have their own README.md file with information on the respective files. 



## Workflow

Here I describe the workflow, which can be followed to fully reproduce the results.



### 1. Prerequisites

For running the whole analysis, one needs following packages/tools/software:

- **BBMap** (v38.84) & **BBDuk**

- **R** (v.4.1.1) along with packages from Bioconductor/CRAN 

- **Linux** shell (we used Ubuntu 20.04) for commands & bash scripts

- **featureCounts** (v2.0.1) from Subread package

- **bedtools** (v2.27.1) 

- **samtools** (v1.12)

- **mafft**  for multiple sequence alignment 
- **fasttree**  for phylogenetic tree construction 
- **barrnap** for 16S rRNA annotation

  



### 2. Mapping

All raw FastQ files should be located in the folder [./data/fastq](data/fastq) . Details on samples and setup of the 
experiment can be found in the methods section of the manuscript. Fastq files etc can be found in the GEO repository
under accession number [GSE284320](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE284320) and downloaded.
Navigate to 
[quality_stats_by_lib](./data/fastq/2023-12-11_Valentina_Cosi_PR23158_raw_FASTQ/quality_stats_by_lib) to 
find fastQC quality statistics of the raw reads. 
To run the mapping, run the bash script [./scripts/trimm_map_BB.sh](./scripts/trimm_map_BB.sh) . 
The script loops through the fastq-files, trims off adapters using BBDuk, maps against the reference _Fusobacterium_ genomes 
(reference fasta and gff files can be found in [./data/reference_sequences/](./data/reference_sequences/)) and counts 
the reads mapped to coding sequences and sRNAs using featureCounts.
Trimming, mapping and counting statistics are stored in the log file [./scripts/stdout_slurm_mapping](./scripts/stdout_slurm_mapping) . 
The directory [./data/rna_align](./data/rna_align) includes all bam-alignment files as well as the count 
tables [counttable_fnn23.txt](./data/rna_align/counttable_fnn23.txt) and [counttable_fnv.txt](./data/rna_align/counttable_fnv.txt) . 
These count table are imported into R for Differential expression analysis, as described below.



### 3. Differential expression analysis

To run the differential expression analysis, run the R  
script [./scripts/fuso_rnaseq.R](./scripts/fuso_rnaseq.R) . 
This outputs all figures of the manuscript, which are saved as PDF and/or SVG files to 
the [./analysis](./analysis) directory. It might take up to 10 minutes to run this script on a low-memory laptop. 


### 4. Phylogenetic tree construction
To construct a phylogenetic tree of the _Fusobacterium_ species, run the bash script 
[./data/phylog_tree/script_to_extract_16s.sh](./data/phylog_tree/script_to_extract_16s.sh).
Final visualization was done using iTOL.

