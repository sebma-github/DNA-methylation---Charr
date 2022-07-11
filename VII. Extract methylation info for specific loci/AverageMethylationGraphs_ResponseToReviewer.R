#Script where I do analyses to answer reviewer comments or questions

####################################### Look at interindividual variation #####################################

#Make average methylation graphs
#This script was made to calculate the inter-individual variation, compare it to the intermorph and time variation.
#To prove that the inter individual variation is not as big, which would validate our use of different fishes for qPCR.
#I will only do it on the genes where we have enough coverage.
#Obviously, ultimately, the reason behind us not using the same samples was because of a lack of tissue.
#So we could very much validate it this way. But without this proof, it could be considered stupid or not as interesting.

#Despite these genes having the most coverage, only 3 of them (H2A, H3l and RASSF4) have residues in all samples.
#So for the other 4, the statistical power is not as strong.

library(tidyverse) #to select columns
library(plyr) #to merge df ?
library(ggplot2) #for graph stuff
library(plotly) #for graph stuff
library(ggpubr) #for ggarrange


#Genes for which I could calculate averages based on residues that were everywhere.
tableHiH3l <- read.table("/Users/sebmatlosz/Desktop/DMRbetweenmorphsinfo/LOC111971647.covdf", sep="\t", header=TRUE) #Close to HistH3-like 
tableHiH2A <- read.table("/Users/sebmatlosz/Desktop/DMRbetweenmorphsinfo/LOC111971648.covdf", sep="\t", header=TRUE) #Close to HistH2A-like
tableRASSF4 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/RASSF4.covdf", sep="\t", header=TRUE)
tableH1 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/H1.covdf", sep="\t", header=TRUE)
tableH2B <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/H2B.covdf", sep="\t", header=TRUE)

#
tableARL16 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/ARL16.covdf", sep="\t", header=TRUE)
tableMEGF9 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/MEGF9.covdf", sep="\t", header=TRUE)
tableMeis1 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/Meis1.covdf", sep="\t", header=TRUE)
tableNFIX <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/NFIX.covdf", sep="\t", header=TRUE)

tableMAGUKp55 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/MAGUKp55.covdf", sep="\t", header=TRUE)
tableARMC1 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/ARMC12000.covdf", sep="\t", header=TRUE)
tableGli3 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/Gli3.covdf", sep="\t", header=TRUE)
tableLmtk2 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/Lmtk2.covdf", sep="\t", header=TRUE)
tableNkx23 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/Nkx23.covdf", sep="\t", header=TRUE)
tableRhoguanin37 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/Rhoguanin37.covdf", sep="\t", header=TRUE)
tableSLC9A3R2 <- read.table("/Users/sebmatlosz/Desktop/DMRqPCRinfo/SLC9A3R2.covdf", sep="\t", header=TRUE)

#Should I remove all rows where there are NAs?
#If I do that, I will lose a lot of positions. I could also remove the categories where there isn't enough data.
#But keep others? Would that even make sense?

#But let's try first by removing all rows where there are NAs.
getmethylationaverage <- function(x) {
  unmeth <- x %>% select(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97)
  meth <- x %>% select(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96)
  methmat <- as.matrix(meth)
  unmethmat <- as.matrix(unmeth)
  total <- methmat + unmethmat
  perc <- (methmat*100)/total
  perc <- as.data.frame(perc)
  perc <- perc[,-1]
  percnoNA <- perc[complete.cases(perc), ] #Remove rows where NAs.
  return(percnoNA)
}

getmethylationaverage2 <- function(x) {
  unmeth <- x %>% select(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97)
  meth <- x %>% select(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96)
  methmat <- as.matrix(meth)
  unmethmat <- as.matrix(unmeth)
  total <- methmat + unmethmat
  perc <- (methmat*100)/total
  perc <- as.data.frame(perc)
  perc <- perc[,-1]
  #percnoNA <- perc[complete.cases(perc), ] #Remove rows where NAs.
  return(perc)
}

#I could also decide to keep the rows with NA, then do a colMeans anyway. It could still give me decent results.


HiH3l <- getmethylationaverage(tableHiH3l)
HiH2A <- getmethylationaverage(tableHiH2A)
RASSF4 <- getmethylationaverage(tableRASSF4)
H1 <- getmethylationaverage(tableH1)
H2B <- getmethylationaverage(tableH2B)

ARL16 <- getmethylationaverage2(tableARL16)    # 0 res left if I remove NA lines
MEGF9 <-getmethylationaverage2(tableMEGF9)     # 0 res left if I remove NA lines
Meis1 <- getmethylationaverage2(tableMeis1)    # 0 res left if I remove NA lines
NFIX <-getmethylationaverage2(tableNFIX)       # 0 res left if I remove NA lines
MPP3 <- getmethylationaverage2(tableMAGUKp55)
ARMC1 <- getmethylationaverage2(tableARMC1)
GLI3 <- getmethylationaverage2(tableGli3)
LMTK2 <- getmethylationaverage2(tableLmtk2)
RHOGUA <-getmethylationaverage2(tableRhoguanin37)
SLC9A3R2 <- getmethylationaverage2(tableSLC9A3R2)
NKX23 <- getmethylationaverage2(tableNkx23)




#Now this is for all positions.
#But I can calculate the average over the region:   # ATH, might want to put back na.rm = FALSE
resumeintable <- function(x){
    df <- data.frame("Meth"=colMeans(x), 
                     "Morph"=c("LB","LB","LB","LB","LB","LB","LB","LB","LB","LB","LB","LB",
                               "PI","PI","PI","PI","PI","PI","PI","PI","PI","PI","PI","PI",
                               "PL","PL","PL","PL","PL","PL","PL","PL","PL","PL","PL","PL",
                               "SB","SB","SB","SB","SB","SB","SB","SB","SB","SB","SB","SB"), 
                     "Time"=c("100","150","200","50","100","150","200","50","100","150","200","50",
                              "100","150","200","50","100","150","200","50","100","150","200","50",
                              "50","50","50","100","150","200","100","150","200","100","150","200",
                              "150","200","50","100","100","150","200","50","100","150","200","50"))
    df$MT <- paste(df$Morph,df$Time,sep="")
    return(df)
}

resumeintable2 <- function(x){
  df <- data.frame("Meth"=colMeans(x, na.rm=TRUE), 
                   "Morph"=c("LB","LB","LB","LB","LB","LB","LB","LB","LB","LB","LB","LB",
                             "PI","PI","PI","PI","PI","PI","PI","PI","PI","PI","PI","PI",
                             "PL","PL","PL","PL","PL","PL","PL","PL","PL","PL","PL","PL",
                             "SB","SB","SB","SB","SB","SB","SB","SB","SB","SB","SB","SB"), 
                   "Time"=c("100","150","200","50","100","150","200","50","100","150","200","50",
                            "100","150","200","50","100","150","200","50","100","150","200","50",
                            "50","50","50","100","150","200","100","150","200","100","150","200",
                            "150","200","50","100","100","150","200","50","100","150","200","50"))
  df$MT <- paste(df$Morph,df$Time,sep="")
  return(df)
}

HiH2A_clean <- resumeintable(HiH2A)
HiH3l_clean <- resumeintable(HiH3l)
H1_clean <- resumeintable(H1)
H2B_clean <- resumeintable(H2B)

RASSF4_clean <- resumeintable(RASSF4)

Meis1_clean <- resumeintable2(Meis1)
NFIX_clean <- resumeintable2(NFIX)
MEGF9_clean <- resumeintable2(MEGF9)
ARL16_clean <- resumeintable2(ARL16)
MPP3_clean <- resumeintable2(MPP3)
ARMC1_clean <- resumeintable2(ARMC1)
GLI3_clean <- resumeintable2(GLI3)
LMTK2_clean <- resumeintable2(LMTK2)
Rhogua_clean <- resumeintable2(RHOGUA)
SLC9A3R2_clean <- resumeintable2(SLC9A3R2)
NKX23_clean <- resumeintable2(NKX23)

recaptable <- function(df) {
  y <- tapply(df$Meth, df$MT, mean, na.rm=TRUE)
  return(y)
}

LMTK2 <- recaptable(LMTK2_clean)
H2A <- recaptable(HiH2A_clean)
SLC9A3R2 <- recaptable(SLC9A3R2_clean)
NKX23 <- recaptable(NKX23_clean)
NFIX <- recaptable(NFIX_clean)
RASSF4 <- recaptable(RASSF4_clean)
ARL16 <- recaptable(ARL16_clean)
MEIS1 <- recaptable(Meis1_clean)
H3 <- recaptable(HiH3l_clean)
ARMC1 <- recaptable(ARMC1_clean)
GLI3 <- recaptable(GLI3_clean)
MEGF9 <- recaptable(MEGF9_clean)
MPP3 <- recaptable(MPP3_clean)
Rhogua <- recaptable(Rhogua_clean)

methrecap <- as.data.frame(rbind(LMTK2,H2A,SLC9A3R2,NKX23,NFIX,RASSF4,ARL16,MEIS1,H3,ARMC1,GLI3,MEGF9,MPP3,Rhogua))
methrecapt <- as.data.frame(t(methrecap))

write.csv(methrecapt, "/Users/sebmatlosz/Desktop/AllFilesForSubmission/averagemethrecap_270222.csv", row.names = TRUE)









#Now I should plot this as well.
colors <- c("LB" = "#00BA38FF", "SB" = "#619CFFFF", "PL" = "#F8766DFF", "PI" = "darkmagenta" )

#Actually, the best way might be to plot the mean and standard deviation.
#Calculate Mean and sd
calculateStats <- function(df) {
    newdf <- data.frame("Mean"= tapply(df$Meth, df$MT, mean, na.rm=TRUE), "Sd"=tapply(df$Meth, df$MT, sd)) 
    newdf$Morph <-  substr(row.names(newdf),1,2)
    newdf$Time <- substr(row.names(newdf),3,5)

    newdf$Time <- factor(newdf$Time)
    newdf$Time <- relevel(newdf$Time, "200")
    newdf$Time <- relevel(newdf$Time, "150")
    newdf$Time <- relevel(newdf$Time, "100")
    newdf$Time <- relevel(newdf$Time, "50")
    return(newdf)
}

HiH2A_stats <- calculateStats(HiH2A_clean)
HiH3l_stats <- calculateStats(HiH3l_clean)
H1_stats <- calculateStats(H1_clean)
H2B_stats <- calculateStats(H2B_clean)

RASSF4_stats <- calculateStats(RASSF4_clean)
Meis1_stats <- calculateStats(Meis1_clean)
NFIX_stats <- calculateStats(NFIX_clean)
MEGF9_stats <- calculateStats(MEGF9_clean)
ARL16_stats <- calculateStats(ARL16_clean)
MPP3_stats <- calculateStats(MPP3_clean)
ARMC1_stats <- calculateStats(ARMC1_clean)
GLI3_stats <- calculateStats(GLI3_clean)
LMTK2_stats <- calculateStats(LMTK2_clean)
Rhogua_stats <- calculateStats(Rhogua_clean)
SLC9A3R2_stats <- calculateStats(SLC9A3R2_clean)
NKX23_stats <- calculateStats(NKX23_clean)

#Plot mean and sd
plotting <- function(df, name) {
  
  plot <- ggplot(df, aes(x=Time, y=Mean, color=Morph)) +geom_point(position=position_dodge(width=0.5)) + 
                geom_errorbar(aes(ymin=Mean-Sd, ymax=Mean+Sd), width=.2,position=position_dodge(.5)) +
                labs(x = "Developmental stage (ts)", y = "Methylation %age", color = "Morph", title = name) + 
                scale_color_manual(values = colors) + theme_bw() + theme(plot.title = element_text(hjust = 0.5))
return(plot)
}

HiH2A_plot <- plotting(HiH2A_stats, "H2A-like")
HiH3l_plot <- plotting(HiH3l_stats, "H3-like")
RASSF4_plot <- plotting(RASSF4_stats, "RASSF4-like")
H1_plot <- plotting(H1_stats, "H1")
H2B_plot <- plotting(H2B_stats, "H2B")

#Not sure whether I should do anything with these four.
#The average methylation over the regions was calculated no matter the position of the residues.
#And then the Mean between each biological replicate could not always be done. #That last part shouldn't be too much of an issue.
Meis1_plot <- plotting(Meis1_stats, "MEIS1-like")
NFIX_plot <- plotting(NFIX_stats, "NFIX")
MEGF9_plot <- plotting(MEGF9_stats, "MEGF9")
ARL16_plot <- plotting(ARL16_stats, "ARL16")
MPP3_plot <- plotting(MPP3_stats, "MPP3")
ARMC1_plot <- plotting(ARMC1_stats, "ARMC1")
GLI3_plot <- plotting(GLI3_stats, "GLI3-like")
LMTK2_plot <- plotting(LMTK2_stats, "LMTK2")
Rhogua_plot <- plotting(Rhogua_stats, "ARHGEF37-like")
SLC9A3R2_plot <- plotting(SLC9A3R2_stats, "SLC9A3R2-like")
NKX23_plot <- plotting(NKX23_stats, "NKX23-like")

# figure <- ggarrange(LMTK2_plot,ARMC1_plot,SLC9A3R2_plot,NKX23_plot,NFIX_plot,RASSF4_plot,ARL16_plot,Meis1_plot,HiH3l_plot,
#                     HiH2A_plot,GLI3_plot,MEGF9_plot,MPP3_plot,Rhogua_plot, ncol=3, nrow=5, common.legend = TRUE, legend="right")

figure <- ggarrange(H1_plot, HiH2A_plot, H2B_plot, HiH3l_plot, ncol=2, nrow=2, common.legend = TRUE, legend="right")


# pdf("/Users/sebmatlosz/Desktop/AllFilesForSubmission/Figure5_Modified_230222.pdf", 10,15)
pdf("/Users/sebmatlosz/Desktop/AllFilesForSubmission/Figure8_Histoneaverage_230222.pdf", 7,6)
figure
dev.off()

fourOtherGenes <- ggarrange(Meis1_plot,NFIX_plot, MEGF9_plot,ARL16_plot, ncol=2, nrow=2)
fourOtherGenes

pdf("/Users/sebmatlosz/Desktop/Figuresforpapers/fourOtherGenes_methylationSd_060122.pdf")
fourOtherGenes
dev.off()



              ################### Compare 5% top residues of PC2 to the other PCs ######################
library(UpSetR)
library(tidyr)

PC1 <- read.csv("/Users/sebmatlosz/Desktop/methylation/mostweightPC1.csv")
PC2 <- read.csv("/Users/sebmatlosz/Desktop/methylation/mostweightPC2.csv")
PC3 <- read.csv("/Users/sebmatlosz/Desktop/methylation/mostweightPC3.csv")
PC4 <- read.csv("/Users/sebmatlosz/Desktop/methylation/mostweightPC4.csv")

colnames(PC1) <- c("Weight","chrom","start")
colnames(PC2) <- c("Weight","chrom","start")
colnames(PC3) <- c("Weight","chrom","start")
colnames(PC4) <- c("Weight","chrom","start")

PC1$ID <- paste(PC1$chrom,PC1$start,sep="_")
PC2$ID <- paste(PC2$chrom,PC2$start,sep="_")
PC3$ID <- paste(PC3$chrom,PC3$start,sep="_")
PC4$ID <- paste(PC4$chrom,PC4$start,sep="_")

PC1$PC <- "PC1"
PC2$PC <- "PC2"
PC3$PC <- "PC3"
PC4$PC <- "PC4"

test <- rbind(PC1,PC2,PC3,PC4)

#The following is adapted from a code made by Lea.
#I do not understand everything yet.

smalldatatable <- data.frame(ID=test$ID,PC=test$PC,VALUE=test$Weight)

smalldatatable[is.na(smalldatatable)] <- 0
smalldatatable$VALUE[smalldatatable$VALUE>0] <- 1

wide_smalldatatable <- smalldatatable %>%
  pivot_wider(names_from = PC, values_from = VALUE)
wide_smalldatatable <- as.data.frame(wide_smalldatatable)

wide_smalldatatable$PC1 <- as.integer(wide_smalldatatable$PC1)
wide_smalldatatable$PC2 <- as.integer(wide_smalldatatable$PC2)
wide_smalldatatable$PC3 <- as.integer(wide_smalldatatable$PC3)
wide_smalldatatable$PC4 <- as.integer(wide_smalldatatable$PC4)

wide_smalldatatable$ID <- factor(wide_smalldatatable$ID)

table_significant <- wide_smalldatatable
table_significant[is.na(table_significant)] <- 0


color <- c("PC1" = "#3CA6D0", "PC2" = "#4CCFFA", "PC3" = "#D06D21","PC4"="#7B3514")

png('/Users/sebmatlosz/Desktop/Figuresforpapers/UpsetGraph_48samples_5p.png')
upset(table_significant,
      nsets = 4,
      sets =c("PC4","PC3","PC2","PC1"),
      sets.bar.color = color,
      order.by = "freq",
      # c("Intersection size", "y-axis numbers", "Set Size","Numbers of set size", "Groups", "numbers on bars")
      text.scale = c(2.5,1.5,1,1.1,1.5,2.2),
      keep.order = TRUE)
dev.off()


################################# Plot PCA without 100 to see if there is still a library effect ############################
library(ggplot2)
library(ggpubr)

#Read in csv file directly. Cannot install methylkit on this version of R.
#Decide whether to use the full 14 427 residue or only the 9719 ones that are on placed scaffold (or 4708 on NW_)
meth.min1df <- read.table("/Users/sebmatlosz/Desktop/methylation/methmin1_36samples_no100.csv", sep=",", header=TRUE)

#remove the first columns that I don't want
meth.min1df <- meth.min1df[,-c(1,2,3,4)]
#Remove the number of Ts columns (every 3 rows)
meth.min1df <- meth.min1df[,-seq(0,108,3)]

#Separate the data in two
numCdf <- meth.min1df[,seq(0,72,2)]
coveragedf <- meth.min1df[,seq(1,72,2)]

total <- ((numCdf*100)/coveragedf)

#Should invert the df to plot based on samples. Otherwise I plot the entire 10 000 residues.
totalt <- t(total)

df_PCA <- prcomp(totalt)
df_out <- as.data.frame(df_PCA$x)
#This was mainly for curiosity
df_out$morph <- c("PL","PL","PL","LB","LB","LB","PI","PI","PI","SB","SB","SB",
                  "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB",
                  "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB")


df_out$stage <- c("200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts",
                  "150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts",
                  "50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts")


df_out$stagenum <- c(200,200,200,200,200,200,200,200,200,200,200,200,
                     150,150,150,150,150,150,150,150,150,150,150,150,
                     50,50,50,50,50,50,50,50,50,50,50,50)


levels(df_out$stage)
df_out$stage <- factor(df_out$stage)
df_out$stage <- relevel(df_out$stage, "200ts")
df_out$stage <- relevel(df_out$stage, "150ts")
df_out$stage <- relevel(df_out$stage, "50ts")

levels(df_out$morph)
df_out$morph <- factor(df_out$morph)
df_out$morph <- relevel(df_out$morph,"PL")
df_out$morph <- relevel(df_out$morph,"SB")
df_out$morph <- relevel(df_out$morph,"LB")

percentage <- round(df_PCA$sdev / sum(df_PCA$sdev) * 100, digits=2)
percentage <- paste(colnames(df_out),"(", paste(as.character(percentage),"%", ")", sep="") )
colors <- c("LB" = "#00BA38FF", "SB" = "#619CFFFF", "PL" = "#F8766DFF", "PI" = "darkmagenta" )
#colors <- c("LB" = "#00BA38FF", "SB" = "#619CFFFF", "PL" = "#F8766DFF")
#colors <- c("F" = "red", "M"="blue", "NA"="grey")

PC1og2 <- ggplot(df_out,aes(x=PC1,y=PC2,color=morph, shape=stage)) + geom_point(size=3) + 
  labs(x=percentage[1],y=percentage[2],title = "36samples_no100") + theme_bw() + scale_color_manual(values = colors)

PC1og3 <- ggplot(df_out,aes(x=PC1,y=PC3,color=morph, shape=stage)) + geom_point(size=3) + 
  labs(x=percentage[1],y=percentage[3],title = "36samples_no100") + theme_bw() + scale_color_manual(values = colors)

PC1og4 <- ggplot(df_out,aes(x=PC1,y=PC4,color=morph, shape=stage)) + geom_point(size=3) + 
  labs(x=percentage[1],y=percentage[4],title = "36samples_no100") + theme_bw() + scale_color_manual(values = colors)


PCA_36samples <- ggarrange(PC1og2,PC1og3, PC1og4, ncol=2, nrow=2,common.legend = TRUE, legend="right")
PCA_36samples

pdf("/Users/sebmatlosz/Desktop/Figuresforpapers/PCA_36samples_090122.pdf")
PCA_36samples
dev.off()



############################### Do PCA on 48 samples without the top 5% residues in PC2 ###################

#Load all residues
meth.min1df <- read.table("/Users/sebmatlosz/Desktop/methylation/methmin1_48samples.csv", sep=",", header=TRUE)

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
total$ID <- row.names(total)

#Load top residues outside 2sigmas
PC2_5p <- read.table("/Users/sebmatlosz/Desktop/methylation/PC2outliers_2sigmas_240322.csv", sep=",", header=TRUE)
              # PC2_5p$ID <- paste(PC2_5p$chrom, PC2_5p$start, sep="_")
              # 
              # #replace NW_ and NC_ by chr                             #those lines were for looking at the old file, which was 5% and not 2 sigmas
              # PC2_5p$ID <- gsub("NC_","chr",PC2_5p$ID)
              # PC2_5p$ID <- gsub("NW_","chr",PC2_5p$ID)

#Remove lines in big table where ID is the same as the top 5% residues
top5pvector <- PC2_5p$ID

#Remove the rows where the string in column number x matches a string in the DMRmorphvector
noPC2_5p <- total[!grepl(paste(top5pvector, collapse="|"), total$ID),]
#write.table(noPC2_5p, "/Users/sebmatlosz/Desktop/methylation/methmin1_48samples_no5pofPC2.csv", quote=F,row.names = F, sep=",")

newtotal <- noPC2_5p[,-49]
totalt <- t(newtotal)

df_PCA <- prcomp(totalt)
df_out <- as.data.frame(df_PCA$x)
#This was mainly for curiosity
df_out$morph <- c("PL","PL","PL","LB","LB","LB","PI","PI","PI","SB","SB","SB",
                  "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB",
                  "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB",
                  "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB")
df_out$stage <- c("200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts",
                  "150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts",
                  "100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts",
                  "50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts")
df_out$stagenum <- c(200,200,200,200,200,200,200,200,200,200,200,200,
                     150,150,150,150,150,150,150,150,150,150,150,150,
                     100,100,100,100,100,100,100,100,100,100,100,100,
                     50,50,50,50,50,50,50,50,50,50,50,50)


levels(df_out$stage)
df_out$stage <- factor(df_out$stage)
df_out$stage <- relevel(df_out$stage, "200ts")
df_out$stage <- relevel(df_out$stage, "150ts")
df_out$stage <- relevel(df_out$stage, "100ts")
df_out$stage <- relevel(df_out$stage, "50ts")

levels(df_out$morph)
df_out$morph <- factor(df_out$morph)
df_out$morph <- relevel(df_out$morph,"PI")
df_out$morph <- relevel(df_out$morph,"PL")
df_out$morph <- relevel(df_out$morph,"SB")
df_out$morph <- relevel(df_out$morph,"LB")

percentage <- round(df_PCA$sdev / sum(df_PCA$sdev) * 100, digits=2)
percentage <- paste(colnames(df_out),"(", paste(as.character(percentage),"%", ")", sep="") )
colors <- c("LB" = "#00BA38FF", "SB" = "#619CFFFF", "PL" = "#F8766DFF", "PI" = "darkmagenta" )
#colors <- c("LB" = "#00BA38FF", "SB" = "#619CFFFF", "PL" = "#F8766DFF")
#colors <- c("F" = "red", "M"="blue", "NA"="grey")

PC1og2 <- ggplot(df_out,aes(x=PC1,y=PC2,color=morph, shape=stage)) + geom_point(size=3) + 
  labs(x=percentage[1],y=percentage[2],title = "No cytosines lying outside 2sigmas of PC2") + theme_bw() + scale_color_manual(values = colors)

PC1og3 <- ggplot(df_out,aes(x=PC1,y=PC3,color=morph, shape=stage)) + geom_point(size=3) + 
  labs(x=percentage[1],y=percentage[3],title = "No cytosines lying outside 2sigmas of PC2") + theme_bw() + scale_color_manual(values = colors)

PC1og4 <- ggplot(df_out,aes(x=PC1,y=PC4,color=morph, shape=stage)) + geom_point(size=3) + 
  labs(x=percentage[1],y=percentage[4],title = "No cytosines lying outside 2sigmas of PC2") + theme_bw() + scale_color_manual(values = colors)

PC1og5 <- ggplot(df_out,aes(x=PC1,y=PC5,color=morph, shape=stage)) + geom_point(size=3) + 
  labs(x=percentage[1],y=percentage[5],title = "No cytosines lying outside 2sigmas of PC2") + theme_bw() + scale_color_manual(values = colors)


PCA_48samples_no5pofPC2 <- ggarrange(PC1og2,PC1og3, PC1og4, PC1og5, ncol=2, nrow=2,common.legend = TRUE, legend="right")
PCA_48samples_no5pofPC2

pdf("/Users/sebmatlosz/Desktop/Figuresforpapers/PCA_48samples_nooutside2sigmasofPC2_240322.pdf")
PCA_48samples_no5pofPC2
dev.off()

#Rerun analysis on the first PCs
PCrecap <- function(i) {
  x <- data.frame(cbind(df_out[,i],df_out$morph,df_out$stage,df_out$stagenum))
  colnames(x) <- c("perc","morph","stage","stagenum")
  x$perc <- as.numeric(as.character(x$perc))
  return(x)
}

PC1recap <- PCrecap(1)
PC2recap <- PCrecap(2)
PC3recap <- PCrecap(3)
PC4recap <- PCrecap(4)

linearfunctions <- function(datatable,name) {
  
  #transform the data
  #lambdaValue <- BoxCox.lambda(datatable$perc, method = "guerrero")
  #datatable$transformed <- BoxCox(datatable$perc, lambdaValue)
  datatable$transformed <- datatable$perc
  
  #create the two linear models
  #fixed effect
  fixedeffect <- lm(transformed ~ morph + stagenum, data=datatable)
  #interaction
  interaction <- lm(transformed ~ morph * stagenum, data=datatable)
  
  # check if there is difference between the fixed effect and the interaction
  anovatest <- anova(fixedeffect,interaction)
  
  # if the difference is significant
  if(anovatest$`Pr(>F)`[2] < 0.05) {
    # do an AICc test on both the linear regressions
    aicc.fixed <- AICc(fixedeffect,return.K = FALSE)
    aicc.inter <- AICc(interaction,return.K = FALSE) 
    
    #if fixed effect is better, use it
    if(aicc.fixed < aicc.inter) {
      a <- anova(fixedeffect) #we will use the fixed effect model
      i <- anova(interaction) #Still use the interaction model for interactions
      
      #create dataframe with all the values we want
      test <- data.frame("PC" = name, "Morph(pvalue)" = a$`Pr(>F)`[1], "Time(pvalue)" = a$`Pr(>F)`[2], "Sex(pvalue)" = a$`Pr(>F)`[3],
                         "Morph*Time"=i$`Pr(>F)`[4],"Morph*Sex"=i$`Pr(>F)`[5],"Time*Sex"=i$`Pr(>F)`[6],"Morph*Time*Sex"=i$`Pr(>F)`[7])
    }
    
    
    #if interaction model is better, use it
    else if(aicc.fixed > aicc.inter) {
      c <- anova(interaction) #we will use the interaction model
      
      
      #create dataframe with all the values we want
      test <- data.frame("PC" = name, "Morph(pvalue)" = c$`Pr(>F)`[1], "Time(pvalue)" = c$`Pr(>F)`[2], "Sex(pvalue)" = c$`Pr(>F)`[3] )
    }
  }
  
  #If no difference between the models, we choose fixedeffect (basically same thing than on line #65)
  else if(anovatest$`Pr(>F)`[2] > 0.05) {
    a <- anova(fixedeffect) 
    i <- anova(interaction) #Still use the interaction model for interactions
    
    #create dataframe with all the values we want
    test <- data.frame("PC" = name, "Morph(pvalue)" = a$`Pr(>F)`[1], "Time(pvalue)" = a$`Pr(>F)`[2], "Sex(pvalue)" = a$`Pr(>F)`[3],
                       "Morph*Time"=i$`Pr(>F)`[4],"Morph*Sex"=i$`Pr(>F)`[5],"Time*Sex"=i$`Pr(>F)`[6],"Morph*Time*Sex"=i$`Pr(>F)`[7])
  }
  return(test)
}


#Run the function on the data for all the PCs
PC1table <- linearfunctions(PC1recap,"PC1")
PC2table <- linearfunctions(PC2recap,"PC2")
PC3table <- linearfunctions(PC3recap,"PC3")
PC4table <- linearfunctions(PC4recap,"PC4")


################################# Check positions of DMRs that are close to genes #####################################

df <- read.table("/Users/sebmatlosz/Desktop/coveragefiles/goodcoveragefiles/summaryuniqRRBScharr.cov", sep="\t", header=FALSE)

#Now need to get the region for each one of those positions
library(GenomicRanges)
library(genomation)
gene.obj=readTranscriptFeatures("/Users/sebmatlosz/Documents/annotated_charr_genome/goodcharrgenome.bed")

names(df)[1] <- "chrom"
names(df)[2] <- "start"
df$end <- df$start+1

dfGR <- makeGRangesFromDataFrame(df,
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

##Could just read the percentage of each feature here. 971126 CpG that are covered 10times
##Maybe this is the easiest way. However it doesn't give me the numbers.
#20,26% prom / 196749 instances
#13,59% exon / 55322 instances
#6,55% intron / 25031 instances
#71,47% intergenic / 694024 instances



#Otherwise, see below. But this might take also a lot of time for R to do.
scaffold <- df[["chrom"]]
windowstart <- df[["start"]]
windowend <- df[["end"]]
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

#ATH: Columns might be one of because I didn't include the meth.diff :ATH 
for (row in 1:nrow(betterdf)) {
  prom <- betterdf[row,5]
  exon <- betterdf[row,6]
  intron <- betterdf[row,7]
  if(prom == 1){betterdf[row,10] <- "prom"}
  if(prom == 0 && exon == 1){betterdf[row,10] <- "exon"}
  if(prom == 0 && exon == 0 && intron == 1){betterdf[row,10] <- "intron"}
  if(prom == 0 && exon == 0 && intron == 0){betterdf[row,10] <- "intergenic"
  }
}




