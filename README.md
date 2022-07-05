## Summary
This repository contains all scripts used for the paper "DNA methylation differences during development distinguish sympatric morphs of Arctic charr (*Salvelinus alpinus*)", published in the journal Molecular Ecology.

Analyses involving the MethylKit R package were done under Version 1.1.456 – © 2009-2018 RStudio, Inc.   
The rest of the analyses were done under R version 4.1.2 (2021-11-01) -- "Bird Hippie"

RRBS reads are available on ENA with the accession number PRJEB45551.

## I. To get the methylation data from raw RRBS reads:
## I.1 Assessment of C-T SNPs from whole genome data:
We used whole genome sequencing data from 8 individuals (1 male and 1 female for each morph) to call for C-T SNPs with **gettingSNPs.sh**     
       
This script maps the WGS reads to the genome, changes some formats, filters the reads and selects for C-T SNPs. 

The output is: **filteredCT.vcf**

## I.2 Mask genome with these SNPs
Download Salvelinus sp. genome in fasta format from https://www.ncbi.nlm.nih.gov/assembly/GCF_002910315.2/

"mask genome for SNPs with bedtools maskfasta function" -> add this to script
"index genome with bismark_genome_preparation"  -> ditto

## I.3 Trim and merge reads from RRBS Illumina sequencing and align them to the masked genome to get methylation coverage files:
See **getmethylcoverage.sh**, use **PhiX.fasta** and **adapters.fasta**.

## I.4 Use MethylKit to further filter the coverage files in order to keep only CpGs that have >10X coverage in every one of the 48 samples.
See **xxxxxxx.R**
The output is: 

This file was used for most methylation analyses: whether they be PCA based, generalized linear models or even DSS.

## II. Analyse methylation data at single CpG sites:
A number of analyses were done on this methylation data ( **xxxxxxxx**) in order to identify methylation differences.


## II.1 PCA

## II.2 GLM

# II.3 DSS
At a reviewer's request, we accounted for coverage biases in the glm analyses with the DSS package, see **XXXXX.R**    
The results of DSS with basic models were very similar to the basic glm models, suggesting that coverage was not a major concern in this study.        
Because of the impossibility to use missing data with the DSS package (i.e. Sex information had some NA values), we kept the glm analysis in the paper, and this DSS analyses is thus only provided here for information purposes.

## III. Analyse methylation data at the region level
## III.1 


## II.3 Perform analysis in Rstudio 
download the annotation in .gff3 from https://www.ncbi.nlm.nih.gov/assembly/GCF_002910315.2/
gff2gtl.pl
gtf2bed.pl
"kind of modifycovfiles.sh"

R scripts

R scripts to make the graphs

## III.1 Extract coverage information for specific loci. 
extractgeneinfobetter.sh 
Needs to be run in the directory where the coverage files are (generated in II.2).
Use DMRsbetweenmorphsforextraction.tsv / DMRsbetweenstageforextraction.tsv and DMRsqPCRgenesforextraction.tsv to extract the respective loci.

# III.2 Use those coverage files to create graphs at both residue and DMR level.

R scripts

R scripts
