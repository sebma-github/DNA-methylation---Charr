library(methylKit)

##Import the matrix of choice (.xlsx) and be sure that treatment and gender are in num and not chr
#There are different sheets in the Matrices_for_analysis_in_R.xlsx file. Make sure you select the "BentvsLimn".
##Matrix must be 5 columns in this order: RealName  FileName  SimplifiedName  Treatment  Sex
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

##Unite all samples so that only the CpGs with >10X coverage from all samples are retained. 
meth.min1=unite(filtered.myobj2, destrand=FALSE)

saveRDS(meth.min1, file = "~/methmin1_48samples.rds")

