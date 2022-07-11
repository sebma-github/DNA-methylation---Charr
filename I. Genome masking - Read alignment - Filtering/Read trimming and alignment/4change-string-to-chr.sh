#!/bin/bash

##To use the MethylKit R package, the coverage files need to have a "chr" string at the beginning of the scaffold name. 
##needs to be run in the directory where the coverage files are.
##for each file ending in bismark.cov remove the first 3 characters of each line. In this case, the "NC_" string
for f in ./*bismark.cov; do
        cut -c 1,2,3 --complement "$f" > "$f"temp
done
wait

##for each file ending in bismark.covtemp, add the "chr" string to the beginning of each line
for f in ./*bismark.covtemp; do
        awk 'BEGIN{ OFS="\t"}$1="chr"$1' "$f" > "$f"final
done
wait

##Remove the temporary output files
rm *.covtemp
wait

##move all the "*.cov" files to the goodcoveragefile directory
mv *.cov -t ~/goodcoveragefiles/
