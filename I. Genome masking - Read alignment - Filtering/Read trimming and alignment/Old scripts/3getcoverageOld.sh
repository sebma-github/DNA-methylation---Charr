#!/bin/bash

##the script must be run in the directory where the alignments in .bam are (created from the 2aligntogenome.sh script).

##the files need to be processed in sam format and not bam
##Transform all .bam into .sam
ls *bam | xargs -I {} -n 1 samtools view -h {} -o {}.sam

wait

##run bismark methylation_extractor on all files ending in .sam in the directory. Produces a lot of things, including the coverage reports.

for f in ./*.sam; do
	/programs/Bismark_v0.19.0/bismark_methylation_extractor --comprehensive --merge_non_CpG --samtools_path /programs/samtools-1.7/ --bedgraph --scaffolds --cytosine_report --genome_folder ~/genomemaskedwithCTSNPs "$f"
done

wait

## order all the "useless" output files so it is easier to see where we're at.
mv *M-bias* -t ./Mbiastrash 
wait
mv *.txt -t ./getcoveragereportfiles/100ts
wait
mv *bedGraph* -t ./bedgraphfiles/
wait

##uncompress all the coverage report files for further analysis
gzip -k -d *cov.gz
wait 

##finish ordering everything
mkdir ./coveragefiles/100ts       
wait

mv *.cov* -t ./coveragefiles/100ts
wait

mv *.bam* -t ./bamsamfiles

## The .cov files are the one we will use for further analysis.


