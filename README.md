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
Use **maskgenome.sh** to mask the genome with the SNPs.   
Use **indexgenome.sh** to index the genome with Bismark.  

#### Trim, merge and align reads:
Use **1cleanandmerge.sh** with **PhiX.fasta** and **adapters.fasta** to remove PhiX and adapter sequences, as well as trim for quality and merge the paired reads.  
Use **2aligntogenome.sh** to align the merged reads to the masked Salvelinsu sp. genome.  
Use **3getcoverage.sh** to extract methylation information (confusingly called "coverage" in Bismark).  
Use **4change-string-to-chr.sh** to change the "NC_" or "NW_" string to a "chr" string. This is because of MethylKit requirements.

#### Filter the "coverage" files: (again, this is a weird Bismark syntax. "coverage" files refer to methylation information files)
Use **Filter10Xcov.R** with **Matrices_for_analysis_in_R.xlsx**.     
The output is: **methmin1_48samples.csv**   
This csv file is used for most methylation analyses: whether they be PCA or glm based.  

## II. Analyse methylation data at single CpG sites:
A number of analyses were done on this methylation data (**methmin1_48samples.csv**) in order to identify methylation differences.

#### PCA
**linearModelsPCA.R** uses the methmin1_48samples.csv to run ANOVA analysis on each principal component.  
**ResiduesBehindPCs.R** uses the methmin1_48samples.csv to extract the weight of each CpG for each PC.

#### GLM
Use the script **glm_methylation.R** with the file **methmin1_48samples.csv** to create tables of the statistical significance of each CpG's methylation for each variable (Morph, Time, Sex, Morph by Time). 
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
Use the script **GOanalysis.R**, alongside these CpGs of interest, the annotation from the Salvelinus sp. assembly in .bed format, **canada_genome_overview_protein.tsv** (Salvelinus sp. protein data) and **geneidgo.db** (GO IDs for genes in the assembly), to perform GO analysis.  
The outputs are available in the /results/ folder.

NOTE: to get the annotation in .bed format, download the annotation in .gff3 from https://www.ncbi.nlm.nih.gov/assembly/GCF_002910315.2/
and use two custom scripts: **gff2gtf.pl** and **gtf2bed.pl**


## VI. Repeat Element Analysis

use the text files: **DMRsbetweenmorphsforFASTA.txt**, **DMRsbetweenstageforFASTA.txt**, **nonDMRsbetweenmorphsforFASTA.tsv** and **nonDMRsbetweenstagesforFASTA.tsv**
to have a list of DMRs or non DMRs of interest, with their scaffold and start position. 
  
Then use the script **extractFASTA_from_list.pl** to extract the corresponding FASTA sequences from the unmasked and unindexed Salvelinus sp. genome. 
These FASTA sequences can be found under /outputFASTA/  
  
Repeat enrichment analysis was done using the online tool RepeatMasker with ”rmblast” as a search engine and the database Dfam3.0.  
DNA sequences in FASTA format were tested against the Dfam3.0 Danio rerio transposable elements database. 







## III.1 Extract coverage information for specific loci. 
extractgeneinfobetter.sh 
Needs to be run in the directory where the coverage files are (generated in II.2).
Use DMRsbetweenmorphsforextraction.tsv / DMRsbetweenstageforextraction.tsv and DMRsqPCRgenesforextraction.tsv to extract the respective loci.

# III.2 Use those coverage files to create graphs at both residue and DMR level.

R scripts

R scripts
