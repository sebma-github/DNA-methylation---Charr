#!/bin/bash

#Index with Bismark using 'bismark_genome_preparation'. I think it automatically reds in the genome that is in the folder.
#This makes a genome file where all C->T and another file where all G->A.

/programs/Bismark_v0.19.0/bismark_genome_preparation --path_to_bowtie /programs/bowtie2/ --verbose /data1/canada_genome/Bismark_prepared_genome/
#NOTE: The genome file must be in .fa or .fasta. So you have to rename it if it is .fna

