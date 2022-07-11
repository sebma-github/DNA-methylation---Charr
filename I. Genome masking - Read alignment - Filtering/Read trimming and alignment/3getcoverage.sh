#!/bin/bash

##the script must be run in the directory where the bismark transformed genome is.
##the files need to be processed in sam format and not bam
##looks in a directory and transforms all .bam into .sam
ls *bam | xargs -I {} -n 1 samtools view -h {} -o {}.sam
wait

##run bismark methylation_extractor on all files ending in .sam in the directory. Produces a lot of things, including the coverage reports.
for f in ./*.sam; do
        /programs/Bismark_v0.19.0/bismark_methylation_extractor --comprehensive --merge_non_CpG --samtools_path /programs/samtools-1.7/ --bedgraph --scaffolds --cytosine_report --genome_folder ~/genomemaskedwithCT "$f"
done
wait

## order all the "useless" output files so it is easier to see where we're at.
mv *M-bias* -t ./Mbiastrash
mv *.txt -t ./getcoveragereportfiles/
mv *bedGraph* -t ./bedgraphfiles/
mv *.bam* -t ./bamsamfiles
wait

##uncompress all the coverage report files for further analysis
gzip -k -d *cov.gz
wait

##finish ordering everything
mkdir ./coveragefiles/
mv *.cov* -t ./coveragefiles/
wait
