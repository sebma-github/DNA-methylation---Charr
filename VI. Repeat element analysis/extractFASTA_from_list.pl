#!/usr/bin/perl

#Script to extract FASTA sequences from an .fna file.
#Needs a text file with seqName, begin and end.
#Use the unmasked unindexed genome ~/canada_genome.fna
#Scaffolds need to be written as "NC_" and "NW_" and not "chr"
#Written by Lea.

open my $POSITIONS, '<', '~/nonDMRsbetweenstagesforFASTA.tsv' or die $!;
open my $BLAST, '>', '~/FASTAfromDMRs/FASTAfromnonDMRstage.fa' or die $!;


while(<$POSITIONS>){
    	chomp;
	my ($seqName,$begin,$end) = split(/\t/);
    	open(SAMTOOLS,"/usr/local/bin/samtools faidx ~/canada_genome.fna $seqName:$begin-$end |");
    	while(my $line = <SAMTOOLS>){
		chomp;
    		print $BLAST "$line";
    	}
    	close(SAMTOOLS);
}


close($POSITIONS);


