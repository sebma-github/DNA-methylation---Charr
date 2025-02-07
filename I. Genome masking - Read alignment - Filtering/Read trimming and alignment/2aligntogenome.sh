#!/bin/bash

#find in the folder all files that end in *_merged.fq, then use the function bismark on all of them one after the other
#align them to the masked genome
find ~/Trimmed_reads/ -type f -name \*_merged.fq -exec /programs/Bismark_v0.19.0/bismark --genome ~/genomemaskedwithCT {} \; &&

#once it's done, create a new folder. Put all the generated report files in it. Those files aren't necessary for the rest.
mkdir reportfiles &
mv  *_report.txt -t ~/reportfile/
