#!/bin/bash

YW='\033[1;33m'
NC='\033[0m'

echo

echo -e "${YW}Enter the path to the first read of the paired end data${NC}"

read read1

##example: read1="/data1/RRBS_results/Samples_01-06/MiSeq_V3_RRBS-70291267/FASTQ_Generation_2018-03-29_01_00_33Z-87951917/17LB1200_L001-ds.128ae403dfb34ff8acfd6870b5f82cb8/17LB1200_S4_L001_R1_001.fasta"

echo

echo -e "${YW}Enter the path to the second read of the paired end data${NC}"

read read2

echo

echo -e "${YW}Analysis starting.${NC}"

##sets directories in which to find sequences to trim from the reads
adapters=/data1/RRBS_results/Samples_01-06/Trimmed_reads/adapters.fasta

PhiX=/data1/RRBS_results/PhiX_genome/PhiX.fasta

##Removes ".fasta" so there is no dots 
##y1="/data1/RRBS_results/Samples_01-06/MiSeq_V3_RRBS-70291267/FASTQ_Generation_2018-03-29_01_00_33Z-87951917/17LB1200_L001-ds.128ae403dfb34ff8acfd6870b5f82cb8/17LB1200_S4_L001_R1_001"
y1=${read1%.*}

##Takes y1 and removes all the beginning of the path. Subdirectories 5 and 6 are not always th same, so *
##y11="17LB1200_S4_L001_R1_001"
y11=${y1#/data1/RRBS_results/Samples_01-06/MiSeq_V3_RRBS-70291267/*/*/}

##same for second read. All this helps having better names for files and directories at the end
y2=${read2%.*}

y21=${y2#/data1/RRBS_results/Samples_01-06/MiSeq_V3_RRBS-70291267/*/*/}

##takes the path from read 1 and removes the last 18 characters
## z="/data1/RRBS_results/Samples_01-06/MiSeq_V3_RRBS-70291267/FASTQ_Generation_2018-03-29_01_00_33Z-87951917/17LB1200_L001-ds.128ae403dfb34ff8acfd6870b5f82cb8/17LB1200_S4"
z=${read1%??????????????????}


##takes z and removes the first 6 directories in path.
##w="17LB1200_S4"
w=${z#/data1/RRBS_results/*/*/*/*/}


#Trim for base quality
/programs/bbmap/bbduk.sh in=$read1 in2=$read2 out=R1clean.fq out2=R2clean.fq qtrim=r trimq=10 &&


##Trim for adapter and PhiX sequences
/programs/bbmap/bbduk.sh in=R1clean.fq in2=R2clean.fq out=R1_1clean.fq out2=R2_1clean.fq ref=$adapters,$PhiX ktrim=r 
k=15 tbo=f tpe=f  &&


##Merge the reads
/programs/bbmap/bbmerge.sh in1=R1_1clean.fq in2=R2_1clean.fq out="$w"_merged.fq outu1="$y11"_unmerged.fq outu2="$y21"_unmerged.fq  &&


##Remove the temporary output files  
rm R1clean.fq R2clean.fq R1_1clean.fq R2_1clean.fq &&


mkdir /data1/RRBS_results/Samples_01-06/Trimmed_reads/100ts/"$w" &&

mv -t /data1/RRBS_results/Samples_01-06/Trimmed_reads/100ts/"$w" "$w"_merged.fq "$y11"_unmerged.fq "$y21"_unmerged.fq &&

echo

echo

echo -e "${YW}You are the best !  \e[5my(n,n)y${NC}"

echo

echo















