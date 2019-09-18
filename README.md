# DNA-methylation---Charr
Contains all scripts used for papers "DNA methylation in Salvelinus alpinus: Epigenetics to understand polymorphism" 
and "Methylation reprogramming and early development epigenetic setup in Salmonids" as well as the raw sequencing files.

All R analyses were done under
Version 1.1.456 – © 2009-2018 RStudio, Inc.

The workflow for paper 1 can be resumed as follows:
I.1 Assessment of C-T SNPs from fastq files:
1bwaCharr_new.sh
2mpileupCharr.sh
3bcftools_merge.sh
4vcftools_filter.sh
"take out SNPs of interest with grep" -> need a script here

II.1 Mask genome with those SNPs
download Canadian charr genome in fasta format from https://www.ncbi.nlm.nih.gov/assembly/GCF_002910315.2/
"mask genome for SNPs with bedtools maskfasta function" -> add this to script
"index genome with bismark_genome_preparation"  -> ditto

II.2 Trim and merge reads from Illumina sequencing and align them to masked genome to get methylation coverage files:
1clean_and_mergeBETTER.sh
2aligntogenomeBETTER.sh
3getcoverage.sh
modifycovfiles.sh

II.3 Perform analysis in Rstudio 
download the annotation in .gff3 from https://www.ncbi.nlm.nih.gov/assembly/GCF_002910315.2/
gff2gtl.pl
gtf2bed.pl
"kind of modifycovfiles.sh"

R scripts

R scripts to make the graphs


The workflow for paper 2 can be resumed as follows:
Using files from II.3:
Different R analyses.
