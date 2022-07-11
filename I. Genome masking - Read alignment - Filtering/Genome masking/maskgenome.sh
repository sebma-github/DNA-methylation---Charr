#!/bin/bash

#mask genome with the SNP file
#bedtools maskfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF> -fo <output FASTA>
/programs/bedtools2/bin/bedtools maskfasta -fi ~/canada_genome.fna -bed ~/filteredCT.vcf  -fo ~/testmaskedwithmergedCT.fa





