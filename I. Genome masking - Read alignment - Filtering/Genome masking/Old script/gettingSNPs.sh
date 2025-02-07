#!/bin/sh

index=~/ASM291031v2.fa
input=~/fastqreadsfolder/
output=~/bamfiles/

###Map whole genome sequencing reads to assembly, merge both reads, do that for all 8 samples
bwa mem -t 12 ${index} ${input}S08qual_R1.fastq.gz ${input}S08qual_R2.fastq.gz  > ${output}S08.sam

###Remove unmapped and not properly aligned. Sort the sam. Convert to bam, do that for all 8 samples
samtools view -bS -F 2308 ${output}S08.sam | samtools sort > ${output}S08.bam

###Index the alignments, for all 8 samples
samtools index ${output}S08.bam

#####################################################################

genome=~/ASM291031v2.fa
input=~/bamfiles/
output=~/vcffiles/

# the -t gives information on variance.

###convert .bam to .vcf, for all 8 samples
samtools mpileup -g -t DP -f ${genome} ${input}S01.bam > ${output}S01.bcf &
samtools mpileup -g -t DP -f ${genome} ${input}S02.bam > ${output}S02.bcf &
samtools mpileup -g -t DP -f ${genome} ${input}S03.bam > ${output}S03.bcf &
samtools mpileup -g -t DP -f ${genome} ${input}S04.bam > ${output}S04.bcf &
wait
bcftools call -c -v ${output}S01.bcf > ${output}S01.vcf &
bcftools call -c -v ${output}S02.bcf > ${output}S02.vcf &
bcftools call -c -v ${output}S03.bcf > ${output}S03.vcf &
bcftools call -c -v ${output}S04.bcf > ${output}S04.vcf &
wait

#####################################################################


input=~/vcffiles/
output=~/mergedvcf/

###bgzip input files
bgzip -i ${input}S01.vcf &
bgzip -i ${input}S02.vcf &
bgzip -i ${input}S03.vcf &
bgzip -i ${input}S04.vcf &
bgzip -i ${input}S05.vcf &
bgzip -i ${input}S06.vcf &
bgzip -i ${input}S07.vcf &
bgzip -i ${input}S08.vcf &
wait


###create tabix index
tabix ${input}S01.vcf.gz &
tabix ${input}S02.vcf.gz &
tabix ${input}S03.vcf.gz &
tabix ${input}S04.vcf.gz &
tabix ${input}S05.vcf.gz &
tabix ${input}S06.vcf.gz &
tabix ${input}S07.vcf.gz &
tabix ${input}S08.vcf.gz &
wait

###merge all files 
bcftools merge -o ${output}merged.vcf -O v ${input}S01.vcf.gz ${input}S02.vcf.gz ${input}S03.vcf.gz ${input}S04.vcf.gz ${input}S05.vcf.gz ${input}S06.vcf.gz ${input}S07.vcf.gz ${input}S08.vcf.gz


###Filter VCF file. 
##raw filtering (this basically doesn't do anything, it keeps every SNP. It just removes one column with unneeded data in the file.)
vcftools --vcf merged.vcf --recode --out filtered --remove-filtered-all

##advanced filtering (this is a real filtering, but in the end I do believe that I didn't do any filtering; the logic being that I would rather mask too many false negative rather than missing false positives)
# vcftools --vcf ${input}merged.vcf --remove-filtered-all --max-missing 0.95 --hwe 0.0001 --min-alleles 2 --max-alleles 2 --maf 0.01 --out filtered --recode

####################################################################

###Identify C-T SNPs
grep -P "\t[CT]\t[CT]\t" merged.vcf > mergedCT.vcf 


