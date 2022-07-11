#!/usr/bin/perl

#Script to extract FASTA sequences from an .fna file.
#Needs a text file with seqName, begin and end.
#For human, use genome at /data1/WGBS_human/genome/GRCh38_latest_genomic.fna
#For Salvelinus, use genome at /data1/RRBS_results/SNP_stuff/canada_genome.fna
#Scaffolds need to be written as "NC_" and "NW_" and not "chr"
#Written by Lea.


open my $POSITIONS, '<', '/home/sebastien/tablesforextraction/nonDMRsbetweenstagesforFASTA.tsv' or die $!;
open my $BLAST, '>', '/data1/RRBS_results/FASTAfromDMRs/FASTAfromnonDMRstage.fa' or die $!;


while(<$POSITIONS>){
    	chomp;
	my ($seqName,$begin,$end) = split(/\t/);
    	open(SAMTOOLS,"/usr/local/bin/samtools faidx /data1/RRBS_results/SNP_stuff/canada_genome.fna $seqName:$begin-$end |");
    	while(my $line = <SAMTOOLS>){
		chomp;
    		print $BLAST "$line";
    	}
    	close(SAMTOOLS);
}


close($POSITIONS);


