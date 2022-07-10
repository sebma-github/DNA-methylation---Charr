#!/bin/bash

##Clean and merge reads from RRBS
#Trim for base quality
/programs/bbmap/bbduk.sh in=~/read1 in2=~/read2 out=R1clean.fq out2=R2clean.fq qtrim=r trimq=10 &&

#Trim for adapter and PhiX sequences
adapters = ~/adapters.fasta
PhiX = ~/PhiX.fasta

/programs/bbmap/bbduk.sh in=R1clean.fq in2=R2clean.fq out=R1_1clean.fq out2=R2_1clean.fq ref=$adapters,$PhiX ktrim=r
k=15 tbo=f tpe=f  &&

#Merge the reads
/programs/bbmap/bbmerge.sh in1=R1_1clean.fq in2=R2_1clean.fq out=S01_merged.fq outu1=S01-1_unmerged.fq outu2=S01-2_unmerged.fq  &&

#Remove the temporary output files
rm R1clean.fq R2clean.fq R1_1clean.fq R2_1clean.fq &&

##########################################################

#find in the folder all files that end in *_merged.fq, then use the function bismark on all of them one after the other
#align them to the masked genome
find ~/Trimmed_reads/ -type f -name \*_merged.fq -exec /programs/Bismark_v0.19.0/bismark --genome ~/genomemaskedwithCT {} \; &&

#once it's done, create a new folder. Put all the generated report files in it. Those files aren't necessary for the rest.
mkdir reportfiles &
mv  *_report.txt -t ~/reportfile/

##########################################################

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

###########################################################

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


