## Summary
This repository contains all scripts used for the paper "DNA methylation differences during development distinguish sympatric morphs of Arctic charr (*Salvelinus alpinus*)", published in the journal Molecular Ecology.

Analyses with the MethylKit R package were done under Version 1.1.456 – © 2009-2018 RStudio, Inc.   
The rest of the analyses were done under R version 4.1.2 (2021-11-01) -- "Bird Hippie"

RRBS reads are available on ENA with the accession number **PRJEB45551**.

## I. Genome masking - Read alignment - Filtering:
#### Genome masking:
We used whole genome sequencing data from 8 individuals (1 male and 1 female for each morph) to call for C-T SNPs.  
The list of 4 659 191 SNPs is **mergedCT.vcf** (currently, this file is a bit too big for github, I need to find a place where to make it accessible).   
Download Salvelinus sp. genome in fasta format from https://www.ncbi.nlm.nih.gov/assembly/GCF_002910315.2/          
Use the **maskgenomeandindex.sh**  to mask the genome with the SNPs and index it with Bismark.    

#### Trim, merge and align reads:
Use **1cleanandmerge.sh** with **PhiX.fasta** and **adapters.fasta** to remove PhiX and adapter sequences, as well as trim for quality and merge the paired reads.  
Use **2aligntogenome.sh** to align the merged reads to the masked Salvelinsu sp. genome.  
Use **3getcoverage.sh** to extract methylation information (confusingly called "coverage" in Bismark).  
Use **4change-string-to-chr.sh** to change the "NC_" or "NW_" string to a "chr" string. This is because of MethylKit requirements.

### Filter the coverage files:
We used the MethylKit package to further filter the coverage files in order to keep only CpGs that have >10X coverage in every one of the 48 samples. 
See **BentvsLimn.R**  
The output is: **methmin1_48samples.csv**   
This csv file is used for most methylation analyses: whether they be PCA or glm based.  

## II. Analyse methylation data at single CpG sites:
A number of analyses were done on this methylation data (**methmin1_48samples.csv**) in order to identify methylation differences.

#### PCA
**linearModelsPCA.R** uses the methmin1_48samples.csv to run ANOVA analysis on each principal component.  
**ResiduesBehindPCs.R** uses the methmin1_48samples.csv to extract the weight of each CpG for each PC.

#### GLM
Use the script glm_methylation.R with the file **methmin1_48samples.csv** to create tables of the statistical significance of each CpG's methylation for each variable (Morph, Time, Sex, Morph by Time). 
These tables can be found in /results/

## III. Analyse methylation data at the region level
###  


## IV. qPCR Analysis
The raw data files from the qPCR experiment are available in the "/raw_qPCR_data/" folder.        
The script **qPCR_graphmaking.R** uses these files to create recaps of the Relative Expression Ratios (RErs) for each gene and make qqplots/bar graphs.   
REr tables can be found in the /tables/ folder in RDS format.  
The script **qPCRanalysis.R** uses the REr files to run linear models on the gene expression data.       
Outputs can be found in the /results/ folder.    

## V. GO analysis
GO analysis can be run on genes located close to CpGs of interest.  
Tables of CpGs of interest can be found in the outputs of GLM or PCA analysis.
Using the script **GOanalysis.R**, alongside these CpGs of interest, the annotation from the Salvelinus sp. assembly in .bed format, **canada_genome_overview_protein.tsv** (Salvelinus sp. protein data) and **geneidgo.db** (GO IDs for genes in the assembly), run GO analysis.  
The outputs are available in the /results/ folder.


## VI. Repeat Element Analysis
















## II.3 Perform analysis in Rstudio 
download the annotation in .gff3 from https://www.ncbi.nlm.nih.gov/assembly/GCF_002910315.2/
gff2gtl.pl
gtf2bed.pl
"kind of modifycovfiles.sh"



## III.1 Extract coverage information for specific loci. 
extractgeneinfobetter.sh 
Needs to be run in the directory where the coverage files are (generated in II.2).
Use DMRsbetweenmorphsforextraction.tsv / DMRsbetweenstageforextraction.tsv and DMRsqPCRgenesforextraction.tsv to extract the respective loci.

# III.2 Use those coverage files to create graphs at both residue and DMR level.

R scripts

R scripts
