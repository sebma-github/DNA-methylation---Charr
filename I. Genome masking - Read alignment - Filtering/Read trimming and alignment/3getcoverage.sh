#!/bin/bash

##the script must be run in the directory: "/data1/RRBS_results/Alignements/withallCT"

##the files need to be processed in sam format and not bam
##looks in a directory and transforms all .bam into .sam
ls *bam | xargs -I {} -n 1 samtools view -h {} -o {}.sam

wait

##run bismark methylation_extractor on all files ending in .sam in the directory. Produces a lot of things, including the coverage reports.

for f in ./*.sam; do
	/programs/Bismark_v0.19.0/bismark_methylation_extractor --comprehensive --merge_non_CpG --samtools_path /programs/samtools-1.7/ --bedgraph --scaffolds --cytosine_report --genome_folder /data1/RRBS_results/SNP_stuff/genomemaskedwithCT "$f"
done

wait

## order all the "useless" output files so it is easier to see where we're at.

mv *M-bias* -t ./Mbiastrash 
wait

mv *.txt -t ./getcoveragereportfiles/50ts
wait

mv *bedGraph* -t ./bedgraphfiles/
wait

##uncompress all the coverage report files for further analysis
gzip -k -d *cov.gz
wait 

##finish ordering everything
mkdir ./coveragefiles/50ts       
wait

mv *.cov* -t ./coveragefiles/50ts
wait


mv *.bam* -t ./bamsamfiles


