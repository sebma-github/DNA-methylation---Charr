library(methylKit)


##Import the matrix of choice (.xlsx) and be sure that treatment and gender are in num and not chr
##Matrix must be 5 columns in this order: RealName  FileName  SimplifiedName  Treatment  Sex
library(readxl)
samplerecapdf <- read_excel("/Users/sebmatlosz/Documents/important_excel_files/Matrices_for_analysis_in_R.xlsx", sheet="BentvsLimn",
                            col_types = c("text", "text", "text", 
                                          "numeric", "numeric", "text"))
samplerecapmatrix <- as.matrix(samplerecapdf)
#To remove PCs that have interaction
niche <- c("Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent")
newmatrix <- cbind(samplerecapmatrix,niche)
samplerecapmatrix <- newmatrix

##set working directory where coverage files are
setwd("/Users/sebmatlosz/Desktop/coveragefiles/goodcoveragefiles")

##Need to do this to have the location list in a proper format
file.list=list(as.list(samplerecapmatrix[,2]))
bestfile.list <- file.list[[1]]

##Need to do this to have the name list in a proper format
sample.id=list(as.list(samplerecapmatrix[,3]))
bestsample.id <- sample.id[[1]]

##Change treatment values from chr to num (for some reason the top of the script doesn't change anything)
treatment=c(samplerecapmatrix[,4])
numtreatment <- as.numeric(treatment)

##Run the normal methRead function  control=0 treatment=1
myobj=methRead(bestfile.list,
sample.id=bestsample.id,
assembly="charr",
pipeline="bismarkCoverage",
treatment=numtreatment,
context="CpG",
header=FALSE,
mincov=10)

##Filters the bases. Only bases covered 10X will be kept for each sample. Also gets rid of PCR bias (99,9th percentile)
filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL,
hi.count=NULL,hi.perc=99.9)

##Then normally I should normalize coverage !     TEST      TEST      TEST      TEST
filtered.myobj2=normalizeCoverage(filtered.myobj)

##In order to use the batch effect part, cannot use the relaxed merging version..... AAAAAAAAAAHHHH  also not good to use relaxed when 3samples
meth.min1=unite(filtered.myobj2, destrand=FALSE)
#saveRDS(meth.min1, file = "/Users/sebmatlosz/Desktop/methmin1_48samples.rds")

#PCASamples(meth.min1, adj.lim=c(0.01,0.1), comp=c(1,3))
batch=c(samplerecapmatrix[,6])
batchdf <- data.frame()


as=assocComp(mBase=meth.min1,batchdf)
View(as[["association"]])

##This shows which PCs are induced by Batch effect (here the devstage). 
## Lets remove all of them (below 0.05) in the analysis ==> NOT POSSIBLE with Relax Merging
newObj=removeComp(meth.min1,comp=c(1,2,3))


##Window tile analysis. 1000bp regions having at least 3CpG
tiles.minwbatch=tileMethylCounts(newObj,win.size=1000,step.size=1000, cov.bases=3)

##Get the methylation differences while correcting for overdispersion (there is a difference when we do this)
myDiffoverdispwbatch=calculateDiffMeth(tiles.minwbatch, overdispersion="MN",test="Chisq",mc.cores=1)
##Get regions that have more than 5%methylation difference and qvalue<0,01
myDiff5poverdispwbatch=getMethylDiff(myDiffoverdispwbatch,difference=5,qvalue=0.01)

library(genomation)

gene.obj=readTranscriptFeatures("/Users/sebmatlosz/Documents/annotated_charr_genome/goodcharrgenome.bed")

diffAnn=annotateWithGeneParts(as(myDiff5poverdispwbatch,"GRanges"),gene.obj)



##At this point, all the data we need is here, but we need to create a clean matrix with everything in it in order
##All the important data from myDiff5p
scaffold <- myDiff5poverdispwbatch[["chr"]]
windowstart <- myDiff5poverdispwbatch[["start"]]
windowend <- myDiff5poverdispwbatch[["end"]]
pvalue <- myDiff5poverdispwbatch[["pvalue"]]
qvalue <- myDiff5poverdispwbatch[["qvalue"]]
methdiff <- myDiff5poverdispwbatch[["meth.diff"]]

#scaffold is a factor of a lot of elements. need to put it back to chr.
chrscaffold <- as.character(scaffold)

resumemydiff5p <- cbind(chrscaffold,windowstart,windowend,pvalue,qvalue,methdiff)

##And all the important data from diffAnn (problem, most of it is not the same number of rows because of deleted scaffolds during gff to bed conversion)
targetrow <- diffAnn@dist.to.TSS[["target.row"]]
disttofeature <- diffAnn@dist.to.TSS[["dist.to.feature"]]
featurename <- diffAnn@dist.to.TSS[["feature.name"]]
resumediffAnn <- cbind(targetrow,disttofeature,featurename)

regiontype <- diffAnn@members
resume2 <- cbind(resumemydiff5p, regiontype)

#Now need to merge the 4rows resumediffAnn with the 7rows resume2
x=resume2
y=resumediffAnn
x=as.data.frame(x)
y<- as.data.frame(y)
x<-cbind(x,row.names(x))

finaltable <- merge(x,y,by.y="targetrow",by.x="row.names(x)",all=TRUE)

##Save the result as an R object but also as an excel file. 
###                !!! CHANGE NAME !!!
saveRDS(finaltable, "/Users/sebmatlosz/Desktop/RRBSresults/RRBS_results/2morphanalysis/Withoutsexascov/BentvsLimn/BentvsLimn.Rda")
write.table(finaltable, "/Users/sebmatlosz/Desktop/RRBSresults/RRBS_results/2morphanalysis/Withoutsexascov/BentvsLimn/BentvsLimn.tsv", sep="\t", row.names=F, quote = F)

################################# IMPORTANT ################ IMPORTANT ################## IMPORTANT ###############################

#It is of interest to know how many of the window tiles were NOT differentially methylated.
#And how many of those were located close to genes etc..
#In order to have those numbers, need to reuse the object "myDiffoverdispwbatch", map it to the annotated genome and reanalyze
library(genomation)
gene.obj=readTranscriptFeatures("/Users/sebmatlosz/Documents/annotated_charr_genome/goodcharrgenome.bed")
diffAnn=annotateWithGeneParts(as(myDiffoverdispwbatch,"GRanges"),gene.obj)

scaffold <- myDiffoverdispwbatch[["chr"]]
windowstart <- myDiffoverdispwbatch[["start"]]
windowend <- myDiffoverdispwbatch[["end"]]
chrscaffold <- as.character(scaffold)

resumeallregions <- cbind(chrscaffold,windowstart,windowend)

targetrow <- diffAnn@dist.to.TSS[["target.row"]]
disttofeature <- diffAnn@dist.to.TSS[["dist.to.feature"]]
featurename <- diffAnn@dist.to.TSS[["feature.name"]]
resumediffAnn <- cbind(targetrow,disttofeature,featurename)

regiontype <- diffAnn@members
resume2 <- cbind(resumeallregions, regiontype)

x=resume2
y=resumediffAnn
x=as.data.frame(x)
y<- as.data.frame(y)
x<-cbind(x,row.names(x))

alltiles <- merge(x,y,by.y="targetrow",by.x="row.names(x)",all=TRUE)

write.table(alltiles, "/Users/sebmatlosz/Desktop/RRBSresults/RRBS_results/2morphanalysis/Withoutsexascov/BentvsLimn/BentvsLimnALLTILES.tsv", sep="\t", row.names=F, quote = F)

