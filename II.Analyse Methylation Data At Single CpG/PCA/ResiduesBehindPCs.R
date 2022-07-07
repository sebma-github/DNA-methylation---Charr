### This script is for looking into which residues have the most weight in specific Principal Components
#meth.min1 objects generated using the BentvsLimn.R script
#library(methylKit)
library(ggplot2)
library(stringr)
library(dplyr)

                                        #######For PCA on 48 samples
#Read in csv file directly. Cannot install methylkit on this version of R
meth.min1df <- read.table("/Users/sebmatlosz/Desktop/methylation/methmin1_48samples.csv", sep=",", header=TRUE)

##Use the annotations to create an ID
annotation <- meth.min1df[,c(1,2,3,4)]
annotation$ID <- paste(annotation$chr, annotation$start, sep="_")

#remove the first columns that I don't want.
meth.min1df <- meth.min1df[,-c(1,2,3,4)]
#Remove the number of Ts columns (every 3 rows)
meth.min1df <- meth.min1df[,-seq(0,144,3)]

#Separate the data in two
numCdf <- meth.min1df[,seq(0,96,2)]
coveragedf <- meth.min1df[,seq(1,96,2)]

total <- ((numCdf*100)/coveragedf)
#Put back the annotation in terms of row name. 
rownames(total) <- annotation$ID     

#Should invert the df to plot based on samples. Otherwise I plot the entire 10 000 residues.
totalt <- t(total)

df_PCA <- prcomp(totalt)

#  df_out <- as.data.frame(df_PCA$x)   #This gives back the value of each sample for each PCA
df_residue <- as.data.frame(df_PCA$rotation) #This should give back the value of each residue for each PCA
# Now, the good thing is that the rowname for each residue is still the ID

# df_residue$ID <- row.names(df_residue)
# write.csv(df_residue,"/Users/sebmatlosz/Desktop/methylation/weightsonPC_48samplesPCA.csv", row.names = FALSE)



#Keep only the PCAs I am interested in: PC1, PC4 and PC23
df_residueCrop <- df_residue[,c(1:23)]
#rownames(df_residueCrop)
#colnames(df_residueCrop)

#Subset the data. 
#First, need to pass the row names as a new column.
names <- as.data.frame(str_split_fixed(rownames(df_residueCrop),"_",2))
df_residueCrop$chrom <- names$V1
df_residueCrop$start <- names$V2

#Now take the top 5%: 517 residues for each PC that matters
#Make a function for that
top5percent <- function(x) {
y <- data.frame(PC = abs(df_residueCrop[,x]), chrom = as.character(df_residueCrop$chrom), start = as.numeric(as.character(df_residueCrop$start)))
#Change chr into NC_ or NW_
y$temp <- gsub("chr","",y$chrom)
y$NW <- ifelse(nchar(y$temp)==11, "NW_", "NC_")
y$NCBI <- paste(y$NW,y$temp, sep="")
v <- data.frame(PC = y$PC,y$NCBI,y$start)
z <- top_n(v, 517, PC)
colnames(z) <- c(print(paste0("PC",x)), "chrom","start")
return(z)
}

PC1top5p <- top5percent(1)
PC2top5p <- top5percent(2)
PC3top5p <- top5percent(3)
PC4top5p <- top5percent(4)
PC7top5p <- top5percent(7)
PC10top5p <- top5percent(10)
PC23top5p <- top5percent(23)

write.table(PC2top5p, "/Users/sebmatlosz/Desktop/methylation/mostweightPC2.csv", sep=",")
write.table(PC3top5p, "/Users/sebmatlosz/Desktop/methylation/mostweightPC3.csv", sep=",")
write.table(PC10top5p, "/Users/sebmatlosz/Desktop/methylation/mostweightPC10.csv", sep=",")
write.table(PC7top5p, "/Users/sebmatlosz/Desktop/methylation/mostweightPC7.csv", sep=",")



#############################################################################################################
                                #### For PCA on 27 samples (no PI, no 100)
#meth.min1 <- readRDS(file = "/Users/sebmatlosz/Desktop/ArnarMeeting050321/methmin1_27_noPIno100.rds")
#meth.min1df <- getData(meth.min1)

#Read in csv file directly. Cannot install methylkit on this version of R
meth.min1df <- read.table("/Users/sebmatlosz/Desktop/methylation/methmin1_27_noPIno100.csv", sep=",", header=TRUE)
#OR
meth.min1df <- read.table("/Users/sebmatlosz/Desktop/methylation/methmin1_27_noPIno100_noNW.csv", sep=",", header=TRUE)

##Use the annotations to create an ID
annotation <- meth.min1df[,c(1,2,3,4)]
annotation$ID <- paste(annotation$chr, annotation$start, sep="_")

#remove the first columns that I don't want.
meth.min1df <- meth.min1df[,-c(1,2,3,4)]
#Remove the number of Ts columns (every 3 rows)
meth.min1df <- meth.min1df[,-seq(0,81,3)]

#Separate the data in two
numCdf <- meth.min1df[,seq(0,54,2)]
coveragedf <- meth.min1df[,seq(1,54,2)]

total <- ((numCdf*100)/coveragedf)
#Put back the annotation in terms of row name. 
rownames(total) <- annotation$ID     

#Should invert the df to plot based on samples. Otherwise I plot the entire 14 000 residues.
totalt <- t(total)

df_PCA <- prcomp(totalt)

#  df_out <- as.data.frame(df_PCA$x)   #This gives back the value of each sample for each PCA
df_residue <- as.data.frame(df_PCA$rotation) #This should give back the value of each residue for each PCA
# Now, the good thing is that the rowname for each residue is still the ID
#saveRDS(df_residue, "/Users/sebmatlosz/Desktop/Figuresforpapers/weightofresiduesforeachPC_27samples.rds")

#df_residue <- readRDS("/Users/sebmatlosz/Desktop/Figuresforpapers/weightofresiduesforeachPC_27samples.rds")
#Keep only the PCAs I am interested in: PC2, PC3 and PC5 have a significant anova p.value with morph.
# PC3 and PC18 have morph*time interaction. PC4, PC16 and PC19 have morph*sex interaction. Tbh, none of these look like anything in graphs.
#To know this, use the script linearModelsPCA.R or check directly the table "/Users/sebmatlosz/Desktop/Figuresforpapers/linearmodelsPCA_27samples.csv"
#So the most interesting is PC3.
df_residueCrop <- df_residue[,c(1:5)]

#Subset the data. 
#First, need to pass the row names as a new column.
names <- as.data.frame(str_split_fixed(rownames(df_residueCrop),"_",2))
df_residueCrop$chrom <- names$V1
df_residueCrop$start <- names$V2

#Now take the top 5%: 722/14427 or 486/9719  residues for each PC that matters
#Make a function for that
top5percent <- function(x) {
  y <- data.frame(PC = abs(df_residueCrop[,x]), chrom = as.character(df_residueCrop$chrom), start = as.numeric(as.character(df_residueCrop$start)))
  #Change chr into NC_ or NW_
  y$temp <- gsub("chr","",y$chrom)
  y$NW <- ifelse(nchar(y$temp)==11, "NW_", "NC_")
  y$NCBI <- paste(y$NW,y$temp, sep="")
  v <- data.frame(PC = y$PC,y$NCBI,y$start)
  z <- top_n(v, 486, PC)
  colnames(z) <- c(print(paste0("PC",x)), "chrom","start")
  return(z)
}

#PCs that separate morphs.
PC2top5p <- top5percent(2)
PC3top5p <- top5percent(3)
PC5top5p <- top5percent(5)

write.table(PC2top5p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC2_5p.csv", sep=",")
write.table(PC3top5p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC3_5p.csv", sep=",")
write.table(PC5top5p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC5_5p.csv", sep=",")

#Merge these 3 df and remove duplicates, to have a list of all residues that are in top 5% of some PCs influencing morphs
PC2 <- read.table("/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC2_5p.csv", sep=",")
PC3 <- read.table("/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC3_5p.csv", sep=",")
PC5 <- read.table("/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC5_5p.csv", sep=",")

PCmorphs <- rbind(PC2[,-1],PC3[,-1],PC5[,-1])
PCmorphsuniq <- PCmorphs[!duplicated(PCmorphs[ , c("chrom", "start")]), ]

write.table(PCmorphsuniq, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC-morphs_5p.csv", sep=",")

#PCs that separate time.
PC1top5p <- top5percent(1)
#PC2top5p <- top5percent(2)
write.table(PC1top5p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC1_5p.csv", sep=",")

PC1 <- read.table("/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC1_5p.csv", sep=",")
PCtime <- rbind(PC1[,-1],PC2[,-1])
PCtimeuniq <- PCtime[!duplicated(PCtime[ , c("chrom", "start")]), ]

write.table(PCtimeuniq, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC-time_5p.csv", sep=",")




#Now take the top 1%: 98/9719  residues for each PC that matters
#Make a function for that
top1percent <- function(x) {
  y <- data.frame(PC = abs(df_residueCrop[,x]), chrom = as.character(df_residueCrop$chrom), start = as.numeric(as.character(df_residueCrop$start)))
  #Change chr into NC_ or NW_
  y$temp <- gsub("chr","",y$chrom)
  y$NW <- ifelse(nchar(y$temp)==11, "NW_", "NC_")
  y$NCBI <- paste(y$NW,y$temp, sep="")
  v <- data.frame(PC = y$PC,y$NCBI,y$start)
  z <- top_n(v, 98, PC)
  colnames(z) <- c(print(paste0("PC",x)), "chrom","start")
  return(z)
}

#PCs that separate morphs.
PC2top1p <- top1percent(2)
PC3top1p <- top1percent(3)
PC5top1p <- top1percent(5)

write.table(PC2top1p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC2_1p.csv", sep=",")
write.table(PC3top1p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC3_1p.csv", sep=",")
write.table(PC5top1p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC5_1p.csv", sep=",")

#Merge these 3 df and remove duplicates, to have a list of all residues that are in top 5% of some PCs influencing morphs
PC2 <- read.table("/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC2_1p.csv", sep=",")
PC3 <- read.table("/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC3_1p.csv", sep=",")
PC5 <- read.table("/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC5_1p.csv", sep=",")

PCmorphs <- rbind(PC2[,-1],PC3[,-1],PC5[,-1])
PCmorphsuniq <- PCmorphs[!duplicated(PCmorphs[ , c("chrom", "start")]), ]

write.table(PCmorphsuniq, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC-morphs_1p.csv", sep=",")

#PCs that separate time.
PC1top1p <- top1percent(1)
#PC2top1p <- top1percent(2)
write.table(PC1top1p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC1_1p.csv", sep=",")

PC1 <- read.table("/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC1_1p.csv", sep=",")
PCtime <- rbind(PC1[,-1],PC2[,-1])
PCtimeuniq <- PCtime[!duplicated(PCtime[ , c("chrom", "start")]), ]

write.table(PCtimeuniq, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC-time_1p.csv", sep=",")


#Now take the top 1%: 1458/9719  residues for each PC that matters
#Make a function for that
top15percent <- function(x) {
  y <- data.frame(PC = abs(df_residueCrop[,x]), chrom = as.character(df_residueCrop$chrom), start = as.numeric(as.character(df_residueCrop$start)))
  #Change chr into NC_ or NW_
  y$temp <- gsub("chr","",y$chrom)
  y$NW <- ifelse(nchar(y$temp)==11, "NW_", "NC_")
  y$NCBI <- paste(y$NW,y$temp, sep="")
  v <- data.frame(PC = y$PC,y$NCBI,y$start)
  z <- top_n(v, 1458, PC)
  colnames(z) <- c(print(paste0("PC",x)), "chrom","start")
  return(z)
}

PC1top15p <- top15percent(1)
PC2top15p <- top15percent(2)
PC3top15p <- top15percent(3)
PC5top15p <- top15percent(5)

write.table(PC1top15p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC1_15p.csv", sep=",")
write.table(PC2top15p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC2_15p.csv", sep=",")
write.table(PC3top15p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC3_15p.csv", sep=",")
write.table(PC5top15p, "/Users/sebmatlosz/Desktop/methylation/27samples_noNW_mostweightPC5_15p.csv", sep=",")




############################# OLD STUFF ###################################

#PC1
PC1matters <- df_residueCrop[,c(1,4,5)]
PC1matters$PC1 <- abs(PC1matters$PC1)
PC1topweight <- top_n(PC1matters, 517, PC1)
PC1botweight <- top_n(PC1matters, -9823, PC1)
PC1topweight$FST <- 0.5   ###ATH this isn't a real FST value. But as I want to plot it in the manhattan plot I will use this as a yvalue
PC1botweight$FST <- 0
PC1resume <- rbind(PC1topweight,PC1botweight)
PC1resume$chrom <- as.character(PC1resume$chrom)

#check the overlaps between PC1 resume and PC4 resume: there are some but not much
PC1topweight$name <- paste(PC1topweight$chrom,PC1topweight$start)
PC4topweight$name <- paste(PC4topweight$chrom,PC4topweight$start)
PC23topweight$name <- paste(PC23topweight$chrom,PC23topweight$start)
PC1botweight$name <- paste(PC1botweight$chrom,PC1botweight$start)

z <- PC1topweight$name %in% PC4topweight$name
sum(z)
z2 <- PC1botweight$name %in% PC1topweight$name
sum(z2)
z3 <- PC1topweight$name %in% PC23topweight$name
sum(z3)
#PC4
PC4matters <- df_residueCrop[,c(2,4,5)]
PC4matters$PC4 <- abs(PC4matters$PC4)
PC4topweight <- top_n(PC4matters, 517, PC4)
PC4botweight <- top_n(PC4matters, -9823, PC4)
PC4topweight$FST <- 0.5   ###ATH this isn't a real FST value. But as I want to plot it in the manhattan plot I will use this as a yvalue
PC4botweight$FST <- 0
PC4resume <- rbind(PC4topweight,PC4botweight)
PC4resume$chrom <- as.character(PC4resume$chrom)
#PC23
PC23matters <- df_residueCrop[,c(3,4,5)]
PC23matters$PC23 <- abs(PC23matters$PC23)
PC23topweight <- top_n(PC23matters, 517, PC23)
PC23botweight <- top_n(PC23matters, -9823, PC23)
PC23topweight$FST <- 0.5    ###ATH this isn't a real FST value. But as I want to plot it in the manhattan plot I will use this as a yvalue
PC23botweight$FST <- 0
PC23resume <- rbind(PC23topweight,PC23botweight)
PC23resume$chrom <- as.character(PC23resume$chrom)
#PC3

#PC7

#PC10

#Replace the chr by "NW_" or "NC_" to match NCBI annotation
chr2NCBI <- function(dataframe) {
  for (i in 1:nrow(dataframe)) {
    x <- nchar(dataframe[i,2], type = "chars", allowNA = FALSE, keepNA = NA)
    if (x == 11) {
      dataframe[i,2] <- str_replace(dataframe[i,2], "chr", "NC_")
    } else {
      dataframe[i,2] <- str_replace(dataframe[i,2], "chr", "NW_")
    }
  }
return(dataframe)
}

PC1resumeNCBI <- chr2NCBI(PC1resume)
PC4resumeNCBI <- chr2NCBI(PC4resume)
PC23resumeNCBI <- chr2NCBI(PC23resume)
#Save those as .rds objects for further work
saveRDS(PC1resumeNCBI, "/Users/sebmatlosz/Desktop/Report090321/PC1resumeNCBI.rds")
saveRDS(PC4resumeNCBI, "/Users/sebmatlosz/Desktop/Report090321/PC4resumeNCBI.rds")
saveRDS(PC23resumeNCBI, "/Users/sebmatlosz/Desktop/Report090321/PC23resumeNCBI.rds")

str(df_residueCrop)
str(test)



### Those .RDS object were then used in the script makeManhattan.R








#PC1matters <- subset(df_residueCrop, PC1 >= 0.02 | PC1 <= -0.02, select=c(PC1))  # 486 out of 10340   NOT USING THIS ANYMORE
#PC4matters <- subset(df_residueCrop, PC4 >= 0.02 | PC4 <= -0.02, select=c(PC4))  # 453 out of 10340
#PC23matters <- subset(df_residueCrop, PC23 >= 0.02 | PC23 <= -0.02, select=c(PC23)) #537 out of 10340

#Now I should have all the residues that weigh the most for each of the two components that are interesting
#I can try looking into where they are located in terms of intron, exon, prom, etc..
library(stringr)
#chrstart <- as.data.frame(str_split_fixed(rownames(df_residueCrop),"_",2))
#chrstart <- as.data.frame(str_split_fixed(rownames(PC1matters),"_",2))
chrstart <- as.data.frame(str_split_fixed(rownames(PC1matters),"_",2))
chrstart$chrom <- chrstart$V1
chrstart$start <- as.numeric(as.character(chrstart$V2))
chrstart$end <- chrstart$start+1

#From there I can create .rds objects that only contain NCBI and START
PC1unplaced <- chrstart[c(1:167),c(3,4)]
PC1placed <-  chrstart[-c(1:167),c(3,4)]
PC1unplaced$chrom <- str_replace(PC1unplaced$chrom, "chr", "NW_")
PC1placed$chrom <- str_replace(PC1placed$chrom, "chr", "NC_")

library(GenomicRanges)
library(genomation)
gene.obj=readTranscriptFeatures("/Users/sebmatlosz/Documents/annotated_charr_genome/goodcharrgenome.bed")

dfGR <- makeGRangesFromDataFrame(chrstart,
                                 keep.extra.columns=TRUE,
                                 ignore.strand=FALSE,
                                 seqinfo=NULL,
                                 seqnames.field=c("seqnames", "seqname",
                                                  "chromosome", "chrom",
                                                  "chr", "chromosome_name",
                                                  "seqid"),
                                 start.field="start",
                                 end.field=c("end", "stop"),
                                 strand.field="strand",
                                 starts.in.df.are.0based=FALSE)

diffAnndf=annotateWithGeneParts(as(dfGR,"GRanges"),gene.obj, intersect.chr = FALSE)

#Repartition: 
#All 971126 CpGs covered 10X once: 20,26% prom / 13,59% exon / 6,55% intron / 71,47% intergenic
#All 10340 CpGs covered 10X in all (PCA): 7,85% prom / 3,81% exon / 3,74% intron / 87,79% intergenic
#All 486 CpGs that matter for PC1: 14,40% prom / 10,29% exon / 2,47% intron / 81,69% intergenic
#All 453 CpGs that matter for PC4: 18,76% prom / 9,49% exon / 2,21% intron / 78,81% intergenic
#Do a table for that:


scaffold <- chrstart[["chrom"]]
windowstart <- chrstart[["start"]]
windowend <- chrstart[["end"]]
chrscaffold <- as.character(scaffold)
regiontype <- diffAnndf@members
resdf <- cbind(chrscaffold,windowstart,windowend,regiontype)

targetrow <- diffAnndf@dist.to.TSS[["target.row"]]  
disttofeature <- diffAnndf@dist.to.TSS[["dist.to.feature"]]
featurename <- diffAnndf@dist.to.TSS[["feature.name"]]
resumediffAnndf <- cbind(targetrow,disttofeature,featurename)

#At this point, the two df won't have the same number of rows because of loss of data in the annotation if on unplaced scaff.
#Need to put them back together correctly
x=resdf
y=resumediffAnndf
x=as.data.frame(x)
y<- as.data.frame(y)
x<-cbind(x,row.names(x))

betterdf <- merge(x,y,by.y="targetrow",by.x="row.names(x)",all=TRUE)

#Let's save the tables as .Rds files for later.
#saveRDS(betterdf,"/Users/sebmatlosz/Desktop/Figuresforpapers/ResidueRepartPC4.rds")
#saveRDS(betterdf,"/Users/sebmatlosz/Desktop/Figuresforpapers/ResidueRepartPC1.rds")

#Plot distance to TSS ?
betterdf$distToTSSabs <- abs(as.numeric(as.character(betterdf$disttofeature)))
hist(betterdf$distToTSSabs, breaks=20)


#I can look into how much of them are part of the DMRs I identified
PC1Resdata <- readRDS(file="/Users/sebmatlosz/Desktop/Figuresforpapers/ResidueRepartPC1.rds")
PC4Resdata <- readRDS(file="/Users/sebmatlosz/Desktop/Figuresforpapers/ResidueRepartPC4.rds")
DMRmorphs <- read.table(file = '/Users/sebmatlosz/Documents/DMRsbetweenmorphsRecap.csv', sep = '\t', header = FALSE)
DMRstages <- read.table(file = '/Users/sebmatlosz/Documents/DMRsbetweenstageRecap.csv', sep = '\t', header = FALSE)

PC1Resdata$windowstart <- as.numeric(as.character(PC1Resdata$windowstart))
DMRstages$V3 <- as.numeric(as.character(DMRstages$V3))
DMRstages <- DMRstages[-1,] 

PC4Resdata$windowstart <- as.numeric(as.character(PC4Resdata$windowstart))
DMRmorphs$V3 <- as.numeric(as.character(DMRmorphs$V3))
DMRmorphs <- DMRmorphs[-1,] 

#Need to write a function that checks for every line in the PCA res files: 
      #Is the chrscaffold the same as a line in DMRrecaps
      #Is the residue position in PCAresfiles between windowstart and (windowstart +1000) in DMRrecap

for (row in 1:nrow(PC1Resdata)) {
  chrscaff<- PC1Resdata[row,2]
  position <- PC1Resdata[row,3]
        for (i in 1:nrow(DMRstages)) {
          chrscaff2 <- DMRstages[i,1]
          windowstart <- DMRstages[i,3]
          windowend = windowstart + 1000
          
            if (chrscaff == chrscaff2 && windowstart < position && position < windowend) {
              DMRstages[i,12] <- "TRUE"
              PC1Resdata[row,10] <- "TRUE"}
        }
}

matchbtwPC1resandDMRstage <- dplyr::filter(DMRstages, grepl("TRUE", V12)) #84 DMRs out of 388 have some of the PC1 res
matchbtwPC1resandDMRstage2 <- dplyr::filter(PC1Resdata, grepl("TRUE", V10)) #325 residues out of 486 are in those DMRs

#Do the same with PC4 and morphs DMR

for (row in 1:nrow(PC4Resdata)) {
  chrscaff<- PC4Resdata[row,2]
  position <- PC4Resdata[row,3]
  for (i in 1:nrow(DMRmorphs)) {
    chrscaff2 <- DMRmorphs[i,1]
    windowstart <- DMRmorphs[i,3]
    windowend = windowstart + 1000
    
    if (chrscaff == chrscaff2 && windowstart < position && position < windowend) {
      DMRmorphs[i,13] <- "TRUE"
      PC4Resdata[row,10] <- "TRUE"}
  }
}

matchbtwPC4resandDMRmorph <- dplyr::filter(DMRmorphs, grepl("TRUE", V13)) #32 DMRs out of 140 have some of the PC4 res
matchbtwPC4resandDMRmorph2 <- dplyr::filter(PC4Resdata, grepl("TRUE", V10)) #145 residues out of 453 are in those DMRs

#I can look at how much of them are part of the DMRs I looked at with qPCR
# Some of them are :D
#Mainly histones and RASSF4. Otherwise, there is a lot of functional RNAs there.



##################### PCA on Benthic ####################################
#As Han suggested, I can redo PCA by morph to see what happens. Here, I will just check if sex comes up in early PCs
# Using the same 10340 residues

meth.min1 <- readRDS(file = "/Users/sebmatlosz/Desktop/ArnarMeeting050321/methmin1_48samples.rds")
meth.min1df <- getData(meth.min1)

##Use the annotations to create an ID
annotation <- meth.min1df[,c(1,2,3,4)]
annotation$ID <- paste(annotation$chr, annotation$start, sep="_")

#remove the first columns that I don't want.
meth.min1df <- meth.min1df[,-c(1,2,3,4)]
#Remove the number of Ts columns (every 3 rows)
meth.min1df <- meth.min1df[,-seq(0,144,3)]

#Separate the data in two
numCdf <- meth.min1df[,seq(0,96,2)]
coveragedf <- meth.min1df[,seq(1,96,2)]

total <- ((numCdf*100)/coveragedf)
#Put back the annotation in terms of row name. 
rownames(total) <- annotation$ID     

Benthic <- total[,c(4,5,6,10,11,12,16,17,18,22,23,24,28,29,30,34,35,36,40,41,42,46,47,48)]
#Should invert the df to plot based on samples. Otherwise I plot the entire 10 000 residues.
Benthict <- t(Benthic)

df_PCA <- prcomp(Benthict)

df_out <- as.data.frame(df_PCA$x)
#This was mainly for curiosity
df_out$Morph <- c("LB","LB","LB","SB","SB","SB",
                  "LB","LB","SB","SB","SB","LB",
                  "LB","LB","SB","SB","SB","LB",
                  "LB","LB","SB","SB","SB","LB")
df_out$Stage <- c("200ts","200ts","200ts","200ts","200ts","200ts",
                  "150ts","150ts","150ts","150ts","150ts","150ts",
                  "100ts","100ts","100ts","100ts","100ts","100ts",
                  "50ts","50ts","50ts","50ts","50ts","50ts")
df_out$Libraries <- c("Chip1","Chip1","Chip1","Chip2","Chip2","Chip2",
                      "Chip3","Chip3","Chip3","Chip4","Chip4","Chip4",
                      "Chip5","Chip5","Chip5","Chip6","Chip6","Chip6",
                      "Chip7","Chip7","Chip7","Chip8","Chip8","Chip8")
df_out$Sex <- c("M","F","M","M","F","F",
                "M","M","F","F","M","F",
                "F","M","M","M","NA","F",
                "NA","NA","NA","NA","NA","NA")

#releveling for legend
levels(df_out$Stage)
df_out$Stage <- factor(df_out$Stage)
df_out$Stage <- relevel(df_out$Stage, "50ts")

#calculating percentage and putting it in form for axis title
percentage <- round(df_PCA$sdev / sum(df_PCA$sdev) * 100, digits=2)
percentage <- paste(colnames(df_out),"(", paste(as.character(percentage),"%", ")", sep="") )

#Sex ? 12,13,15
p2<-ggplot(df_out,aes(x=PC1,y=PC15,color=Sex,shape=Stage)) + geom_point(size=3) + 
  labs(x=percentage[1],y=percentage[15],title = "PCA analysis on 24 Benthic arctic charr methylomes") +
  theme(legend.title = element_text(size=20), legend.text = element_text(size=15))
p
p2

p3 <- ggarrange(p,p2,labels=c("A","B"),nrow = 2, ncol = 1, common.legend = TRUE, legend="right")
p3






############################### Keep only the CpGs that are outside the two sigmas for each distribution ################
library(dplyr)

PCAfinal <- read.csv("/Users/sebmatlosz/Desktop/methylation/weightsonPC_48samplesPCA.csv")

#PC1
min <- mean(PCAfinal$PC1) - 2*sd(PCAfinal$PC1)
max <- mean(PCAfinal$PC1) + 2*sd(PCAfinal$PC1)
PC1outliers <- filter(PCAfinal, PC1 < min | PC1 > max)

#PC4. I guess the distribution here isn't symmetrical so it makes sense that I don't end up with 5% of the data. 
min <- mean(PCAfinal$PC4) - 2*sd(PCAfinal$PC4)
max <- mean(PCAfinal$PC4) + 2*sd(PCAfinal$PC4)
PC4outliers <- filter(PCAfinal, PC4 < min | PC4 > max)

#PC2. For redoing PCA without cytosines impacting library
min <- mean(PCAfinal$PC2) - 2*sd(PCAfinal$PC2)
max <- mean(PCAfinal$PC2) + 2*sd(PCAfinal$PC2)
PC2outliers <- filter(PCAfinal, PC2 < min | PC2 > max)


#Save those
write.csv(PC1outliers, "/Users/sebmatlosz/Desktop/methylation/PC1outliers_2sigmas_100222.csv", row.names = FALSE)
write.csv(PC4outliers, "/Users/sebmatlosz/Desktop/methylation/PC4outliers_2sigmas_100222.csv", row.names = FALSE)
write.csv(PC2outliers, "/Users/sebmatlosz/Desktop/methylation/PC2outliers_2sigmas_240322.csv", row.names = FALSE)




