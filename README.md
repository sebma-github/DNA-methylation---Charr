## DNA-methylation_Charr
Contains all scripts used for the paper "DNA methylation differences during development distinguish sympatric morphs of Arctic charr (*Salvelinus alpinus*)".

All R analyses were done under
Version 1.1.456 – © 2009-2018 RStudio, Inc.

## The workflow for methylation assesment can be resumed as follows:
## I.1 Assessment of C-T SNPs from fastq files:
Used whole genome sequencing data from 8 individuals: 1 male and 1 female for each morph.   
     
See **gettingSNPs.sh**     
       
This script maps the WGS sequencing reads to the genome, changes some formats, filters the reads and selects for C-T SNPs. 

The output is: filteredCT.vcf

## II.1 Mask genome with those SNPs
Download Canadian charr genome in fasta format from https://www.ncbi.nlm.nih.gov/assembly/GCF_002910315.2/

"mask genome for SNPs with bedtools maskfasta function" -> add this to script
"index genome with bismark_genome_preparation"  -> ditto

## II.2 Trim and merge reads from Illumina sequencing and align them to masked genome to get methylation coverage files:
See **getmethylcoverage.sh**, use **PhiX.fasta** and **adapters.fasta**.

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
