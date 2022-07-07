## Summary
This repository contains all scripts used for the paper "DNA methylation differences during development distinguish sympatric morphs of Arctic charr (*Salvelinus alpinus*)", published in the journal Molecular Ecology.

Analyses with the MethylKit R package were done under Version 1.1.456 – © 2009-2018 RStudio, Inc.   
The rest of the analyses were done under R version 4.1.2 (2021-11-01) -- "Bird Hippie"

RRBS reads are available on ENA with the accession number **PRJEB45551**.

## I. Get methylation data from raw RRBS reads:
#### Assessment of C-T SNPs from whole genome data:
We used whole genome sequencing data from 8 individuals (1 male and 1 female for each morph) to call for C-T SNPs with **gettingSNPs.sh**       
This script maps the WGS reads to the genome, changes some formats, filters the reads and selects for C-T SNPs.        
The output is: **filteredCT.vcf**

#### Mask genome with these SNPs
Download Salvelinus sp. genome in fasta format from https://www.ncbi.nlm.nih.gov/assembly/GCF_002910315.2/      
"mask genome for SNPs with bedtools maskfasta function" -> add this to script       
"index genome with bismark_genome_preparation"  -> ditto       
see **maskgenomeandindex.sh**      

#### Trim and merge reads from RRBS Illumina sequencing and align them to the masked genome to get methylation coverage files:
See **getmethylcoverage.sh**, use **PhiX.fasta** and **adapters.fasta**.

### Use MethylKit to further filter the coverage files in order to keep only CpGs that have >10X coverage in every one of the 48 samples.
See **xxxxxxx.R**
The output is: 

This file was used for most methylation analyses: whether they be PCA based, generalized linear models or even DSS.

## II. Analyse methylation data at single CpG sites:
A number of analyses were done on this methylation data ( **xxxxxxxx**) in order to identify methylation differences.

#### PCA

### GLM


## III. Analyse methylation data at the region level
###  


## IV. qPCR Analysis
The raw data files from the qPCR experiment are available in the "/raw_qPCR_data/" folder.        
The script **qPCR_graphmaking.R** uses these files to create recaps of the Relative Expression Ratios (RErs) for each gene and make qqplots/bar graphs.   
REr tables can be found in the /tables/ folder in RDS format.  
The script **qPCRanalysis.R** uses the REr files to run linear models on the gene expression data.       
Outputs can be found in the /results/ folder.    

## V. GO analysis




















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
