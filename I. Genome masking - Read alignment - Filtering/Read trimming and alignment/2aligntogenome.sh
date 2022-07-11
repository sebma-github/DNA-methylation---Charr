#!/bin/bash

#find in the folder all files that end in *_merged.fq, then use the function bismark on all of them one after the other.
find ~/100ts/ -type f -name \*_merged.fq -exec /programs/Bismark_v0.19.0/bismark --genome ~/genomemaskedwithCTSNPs {} \; &&

#once it's done, create a new folder. Put all the generated report files that we do not want in it.
mkdir reportfiles & 

mv  *_report.txt -t ~/alignmentsreportfile/100ts



