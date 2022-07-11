#!/bin/bash



#find in the folder all files that end in *_merged.fq, then use the function bismark on all of them one after the other.

find /data1/RRBS_results/Samples_01-06/Trimmed_reads/50ts/ -type f -name \*_merged.fq -exec /programs/Bismark_v0.19.0/bismark --genome /data1/RRBS_results/SNP_stuff/genomemaskedwithCT {} \; &&


#once it's done, create a new folder. Put all the generated report files in it.

#mkdir reportfiles & 

mv  *_report.txt -t /data1/RRBS_results/Alignements/withallCT/alignmentsreportfile/50ts



