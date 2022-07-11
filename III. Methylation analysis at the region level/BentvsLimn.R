library(methylKit)
library(genomation)


##Import the matrix of choice (.xlsx) and be sure that treatment and gender are in num and not chr
#There are different sheets in the Matrices_for_analysis_in_R.xlsx file. Make sure you select the "BentvsLimn".
##Matrix must be 6 columns in this order: RealName  FileName  SimplifiedName  Treatment  Sex Stage
#You will have to rename the files based on how you saved them.
library(readxl)
samplerecapdf <- read_excel("~/Matrices_for_analysis_in_R.xlsx", sheet="BentvsLimn",
                            col_types = c("text", "text", "text", 
                                          "numeric", "numeric", "text"))
samplerecapmatrix <- as.matrix(samplerecapdf)

##set working directory where coverage files are
setwd("~/goodcoveragefiles")

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

##Normalize coverage
filtered.myobj2=normalizeCoverage(filtered.myobj)

##In order to use the batch effect part
meth.min1=unite(filtered.myobj2, destrand=FALSE)

##Check which PCs are induced by Batch effect (here the developmental stage is in column 6). 
batch=c(samplerecapmatrix[,6])
batchdf <- data.frame()

as=assocComp(mBase=meth.min1,batchdf)
View(as[["association"]])

## Lets remove all the significant ones (below 0.05) in the analysis
newObj=removeComp(meth.min1,comp=c(1,2,3))

##Window tile analysis. 1000bp regions having at least 3CpG
tiles.minwbatch=tileMethylCounts(newObj,win.size=1000,step.size=1000, cov.bases=3)

##Get the methylation differences while correcting for overdispersion
myDiffoverdispwbatch=calculateDiffMeth(tiles.minwbatch, overdispersion="MN",test="Chisq",mc.cores=1)

##Get regions that have more than 5%methylation difference and qvalue<0,01
myDiff5poverdispwbatch=getMethylDiff(myDiffoverdispwbatch,difference=5,qvalue=0.01)

gene.obj=readTranscriptFeatures("~/goodcharrgenome.bed")
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

##And all the important data from diffAnn 
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

write.table(finaltable, "~/BentvsLimn_DMRs.tsv", sep="\t", row.names=F, quote = F)

################################# IMPORTANT ################ IMPORTANT ################## IMPORTANT ###############################
#It is also of interest to know how many of the window tiles were NOT differentially methylated.
#And how many of those were located close to genes etc..
#In order to have those numbers, need to reuse the object "myDiffoverdispwbatch", map it to the annotated genome and reanalyze
library(genomation)
gene.obj=readTranscriptFeatures("~/goodcharrgenome.bed")
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

write.table(alltiles, "~/BentvsLimn/BentvsLimnALLTILES.tsv", sep="\t", row.names=F, quote = F)

