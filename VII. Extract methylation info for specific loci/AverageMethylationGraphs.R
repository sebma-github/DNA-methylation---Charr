library(tidyverse) 
library(plyr) 
library(ggplot2) 
library(plotly) 
library(ggpubr) 

#Make average methylation graphs

#Genes for which I could calculate averages based on residues that were everywhere.
tableHiH3l <- read.table("~/LOC111971647.covdf", sep="\t", header=TRUE) #Close to HistH3-like 
tableHiH2A <- read.table("~/LOC111971648.covdf", sep="\t", header=TRUE) #Close to HistH2A-like
tableRASSF4 <- read.table("~/RASSF4.covdf", sep="\t", header=TRUE)
tableH1 <- read.table("~/H1.covdf", sep="\t", header=TRUE)
tableH2B <- read.table("~/H2B.covdf", sep="\t", header=TRUE)

#
tableARL16 <- read.table("~/ARL16.covdf", sep="\t", header=TRUE)
tableMEGF9 <- read.table("~/MEGF9.covdf", sep="\t", header=TRUE)
tableMeis1 <- read.table("~/Meis1.covdf", sep="\t", header=TRUE)
tableNFIX <- read.table("~/NFIX.covdf", sep="\t", header=TRUE)

tableMAGUKp55 <- read.table("~/MAGUKp55.covdf", sep="\t", header=TRUE)
tableARMC1 <- read.table("~/ARMC12000.covdf", sep="\t", header=TRUE)
tableGli3 <- read.table("~/Gli3.covdf", sep="\t", header=TRUE)
tableLmtk2 <- read.table("~/Lmtk2.covdf", sep="\t", header=TRUE)
tableNkx23 <- read.table("~/Nkx23.covdf", sep="\t", header=TRUE)
tableRhoguanin37 <- read.table("~/Rhoguanin37.covdf", sep="\t", header=TRUE)
tableSLC9A3R2 <- read.table("/~/SLC9A3R2.covdf", sep="\t", header=TRUE)


#Remove all rows where there are NAs.
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

#For genes where there is not enough coverage, keep rows with NAs
getmethylationaverage2 <- function(x) {
  unmeth <- x %>% select(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97)
  meth <- x %>% select(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96)
  methmat <- as.matrix(meth)
  unmethmat <- as.matrix(unmeth)
  total <- methmat + unmethmat
  perc <- (methmat*100)/total
  perc <- as.data.frame(perc)
  perc <- perc[,-1]
  return(perc)
}

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


#To calculate the average over the region:
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

write.csv(methrecapt, "~/averagemethrecap_270222.csv", row.names = TRUE)



#Now to plot things
colors <- c("LB" = "#00BA38FF", "SB" = "#619CFFFF", "PL" = "#F8766DFF", "PI" = "darkmagenta" )

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

figure <- ggarrange(LMTK2_plot,ARMC1_plot,SLC9A3R2_plot,NKX23_plot,NFIX_plot,RASSF4_plot,ARL16_plot,Meis1_plot,HiH3l_plot,
                     HiH2A_plot,GLI3_plot,MEGF9_plot,MPP3_plot,Rhogua_plot, ncol=3, nrow=5, common.legend = TRUE, legend="right")

# figure <- ggarrange(H1_plot, HiH2A_plot, H2B_plot, HiH3l_plot, ncol=2, nrow=2, common.legend = TRUE, legend="right")


pdf("~/Figure5_Modified_230222.pdf", 10,15)
# pdf("~/Figure8_Histoneaverage_230222.pdf", 7,6)
figure
dev.off()

