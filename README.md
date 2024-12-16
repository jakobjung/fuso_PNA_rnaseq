# Title

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

- BBMap (v38.84) & BBDuk

- R (v.4.1.1) along with packages from Bioconductor/CRAN 

- Linux shell (we used Ubuntu 20.04) for commands & bash scripts

- featureCounts (v2.0.1) from Subread package

- bedtools (v2.27.1) 

- samtools (v1.12)

  



### 2. Steps to reproduce the analysis

#### 2.1. Data 


