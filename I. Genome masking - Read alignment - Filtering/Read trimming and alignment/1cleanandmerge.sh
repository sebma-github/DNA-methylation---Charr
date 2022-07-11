#!/bin/bash

##Clean and merge reads from RRBS
#Trim for base quality
/programs/bbmap/bbduk.sh in=~/read1 in2=~/read2 out=R1clean.fq out2=R2clean.fq qtrim=r trimq=10 &&

#Trim for adapter and PhiX sequences
adapters = ~/adapters.fasta
PhiX = ~/PhiX.fasta

/programs/bbmap/bbduk.sh in=R1clean.fq in2=R2clean.fq out=R1_1clean.fq out2=R2_1clean.fq ref=$adapters,$PhiX ktrim=r
k=15 tbo=f tpe=f  &&

#Merge the reads
/programs/bbmap/bbmerge.sh in1=R1_1clean.fq in2=R2_1clean.fq out=S01_merged.fq outu1=S01-1_unmerged.fq outu2=S01-2_unmerged.fq  &&

#Remove the temporary output files
rm R1clean.fq R2clean.fq R1_1clean.fq R2_1clean.fq &&
