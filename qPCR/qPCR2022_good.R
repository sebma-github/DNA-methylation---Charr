library(dplyr)
library(reshape2)
library(ggplot2)
          
####################### Load the qPCR data ##########################
#Sometimes, one ofthe technical replicates was bad, so I removed it here and replaced the mean by the value of the good technical replicate.
#I could have done it directly in Excel, but I just did it here.

# t100 - Lmtk2 & HiH2A
b1.100 <- read.csv("~/Plate8_100_Lmtk2HiH2a.csv",header=TRUE)
# t150 - Lmtk2 & HiH2a
b1.150 <- read.csv("~/Plate5_150_Lmtk2HiH2a.csv",header=TRUE)
# t200 - Lmtk2 & HiH2a
b1.200 <- read.csv("~/Plate1_200_Lmtk2HiH2a.csv",header=TRUE)
b1.200 <- b1.200[b1.200$Well!="G5",]
b1.200$Ct.Mean[b1.200$Well=="H5"] <- b1.200$Ct[b1.200$Well=="H5"] 

# t100 - MyelTF1 & SLC9A3R2
b2.100 <- read.csv("~/Plate9_100_MyelTF1SLC9A3R2.csv",header=TRUE)
b2.100 <- b2.100[b2.100$Well!="C12",]
b2.100$Ct.Mean[b2.100$Well=="D12"] <- b2.100$Ct[b2.100$Well=="D12"] 
b2.100 <- b2.100[b2.100$Well!="H12",]
b2.100$Ct.Mean[b2.100$Well=="G12"] <- b2.100$Ct[b2.100$Well=="G12"]
# t150 - MyelTF1 & SLC9A3R2
b2.150 <- read.csv("~/Plate6_150_MyelTF1SLC9A3R2.csv",header=TRUE)
b2.150 <- b2.150[b2.150$Well!="D10",]
b2.150$Ct.Mean[b2.150$Well=="C10"] <- b2.150$Ct[b2.150$Well=="C10"]
# t200 - MyelTF1
b2.200Myel <- read.csv("~/Plate2_200_MyelTF1LRP11.csv",header=TRUE)
# t200 - SLC9A3R2
b2.200SLC <- read.csv("~/Plate3_200_SLC9A3R2ZFP2.csv",header=TRUE)
b2.200SLC <- b2.200SLC[b2.200SLC$Well!="G3",]
b2.200SLC$Ct.Mean[b2.200SLC$Well=="H3"] <- b2.200SLC$Ct[b2.200SLC$Well=="H3"] 
b2.200SLC <- b2.200SLC[b2.200SLC$Well!="B4",]
b2.200SLC$Ct.Mean[b2.200SLC$Well=="A4"] <- b2.200SLC$Ct[b2.200SLC$Well=="A4"] 

# t100 - Nkx23 and Ets2
b3.100 <- read.csv("~/Plate10_100_Nkx23Ets2.csv",header=TRUE)
b3.100 <- b3.100[b3.100$Well!="F3",]
b3.100$Ct.Mean[b3.100$Well=="E3"] <- b3.100$Ct[b3.100$Well=="E3"] 
b3.100 <- b3.100[b3.100$Well!="F4",]
b3.100$Ct.Mean[b3.100$Well=="E4"] <- b3.100$Ct[b3.100$Well=="E4"]
# t150 - Nkx23 and Ets2
b3.150 <- read.csv("~/Plate7_150_Nkx23Ets2.csv",header=TRUE)
b3.150 <- b3.150[b3.150$Well!="E3",]
b3.150$Ct.Mean[b3.150$Well=="F3"] <- b3.150$Ct[b3.150$Well=="F3"] 
b3.150 <- b3.150[b3.150$Well!="C7",]
b3.150$Ct.Mean[b3.150$Well=="D7"] <- b3.150$Ct[b3.150$Well=="D7"] 
b3.150 <- b3.150[b3.150$Well!="B9",]
b3.150$Ct.Mean[b3.150$Well=="A9"] <- b3.150$Ct[b3.150$Well=="A9"]
# t200 - Nkx23 and Ets2
b3.200 <- read.csv("~/Plate4_200_Nkx23Ets2.csv",header=TRUE)

# t100 - NFiX-1 and RASSF4-2               H3 and H7 aren't good
b4.100 <- read.csv("~/Plate16_100_NFiX-1RASSF4-2.csv",header=TRUE)
b4.100 <- b4.100[b4.100$Well!="H3",]
b4.100$Ct.Mean[b4.100$Well=="G3"] <- b4.100$Ct[b4.100$Well=="G3"]
b4.100 <- b4.100[b4.100$Well!="H7",]
b4.100$Ct.Mean[b4.100$Well=="G7"] <- b4.100$Ct[b4.100$Well=="G7"]

# t150 - NFiX-1 and RASSF4-2               F12 ain't good
b4.150 <- read.csv("~/Plate15_150_NFiX-1RASSF4-2.csv",header=TRUE)
b4.150 <- b4.150[b4.150$Well!="F12",]
b4.150$Ct.Mean[b4.150$Well=="E12"] <- b4.150$Ct[b4.150$Well=="E12"]

# t200 - NFiX-1 and RASSF4-2              H12 ain't good
b4.200 <- read.csv("~/Plate14_200_NFiX-1RASSF4-2.csv",header=TRUE)
b4.200 <- b4.200[b4.200$Well!="H12",]
b4.200$Ct.Mean[b4.200$Well=="G12"] <- b4.200$Ct[b4.200$Well=="G12"]

# t100 - ARF16-1 and MEiS1l-2           Not good: D9, D10, D11, D12, F10, H12
b5.100 <- read.csv("~/Plate19_100_ARF16-1MEiS1l-2.csv", header=TRUE)
#because there were some "Undetermined" Ct values, here Ct isn't numeric. Need to convert it
b5.100$Ct <- as.numeric(as.character(b5.100$Ct))

b5.100 <- b5.100[b5.100$Well!="D9",]
b5.100$Ct.Mean[b5.100$Well=="C9"] <- b5.100$Ct[b5.100$Well=="C9"]
b5.100 <- b5.100[b5.100$Well!="D10",]
b5.100$Ct.Mean[b5.100$Well=="C10"] <- b5.100$Ct[b5.100$Well=="C10"]
b5.100 <- b5.100[b5.100$Well!="D11",]
b5.100$Ct.Mean[b5.100$Well=="C11"] <- b5.100$Ct[b5.100$Well=="C11"]
b5.100 <- b5.100[b5.100$Well!="D12",]
b5.100$Ct.Mean[b5.100$Well=="C12"] <- b5.100$Ct[b5.100$Well=="C12"]
b5.100 <- b5.100[b5.100$Well!="F10",]
b5.100$Ct.Mean[b5.100$Well=="E10"] <- b5.100$Ct[b5.100$Well=="E10"]
b5.100 <- b5.100[b5.100$Well!="H12",]
b5.100$Ct.Mean[b5.100$Well=="G12"] <- b5.100$Ct[b5.100$Well=="G12"]

# t150 - ARF16-1 and MEiS1l-2           D12, F12 and H12 aren't good
b5.150 <- read.csv("~/Plate18_150_ARF16-1MEiS1l-2.csv", header=TRUE)
b5.150 <- b5.150[b5.150$Well!="D12",]
b5.150$Ct.Mean[b5.150$Well=="C12"] <- b5.150$Ct[b5.150$Well=="C12"]
b5.150 <- b5.150[b5.150$Well!="F12",]
b5.150$Ct.Mean[b5.150$Well=="E12"] <- b5.150$Ct[b5.150$Well=="E12"]
b5.150 <- b5.150[b5.150$Well!="H12",]
b5.150$Ct.Mean[b5.150$Well=="G12"] <- b5.150$Ct[b5.150$Well=="G12"]

# t200 - ARF16-1 and MEiS1l-2          D3 and H8 ain't good
b5.200 <- read.csv("~/Plate17_200_ARF16-1MEiS1l-2.csv", header=TRUE)
b5.200 <- b5.200[b5.200$Well!="D3",]
b5.200$Ct.Mean[b5.200$Well=="C3"] <- b5.200$Ct[b5.200$Well=="C3"]
b5.200 <- b5.200[b5.200$Well!="H8",]
b5.200$Ct.Mean[b5.200$Well=="G8"] <- b5.200$Ct[b5.200$Well=="G8"]

# t100 - HiH3l-2 and ARMC1-1          E4 and F4 are really different, but I can't do anything about it. 
b6.100 <- read.csv("~/Plate22_100_HiH3l2ARMC1-1.csv", header=TRUE)

# t150 - HiH3l-2 and ARMC1-1          H12 ain't good
b6.150 <- read.csv("~/Plate21_150_HiH3l2ARMC1-1.csv", header=TRUE)
b6.150 <- b6.150[b6.150$Well!="H12",]
b6.150$Ct.Mean[b6.150$Well=="G12"] <- b6.150$Ct[b6.150$Well=="G12"]

# t200 - HiH3l-2 and ARMC1-1          G2 not good. H2 slightly better
b6.200 <- read.csv("~/Plate20_200_HiH3l2ARMC1-1.csv", header=TRUE)
b6.200 <- b6.200[b6.200$Well!="G2",]
b6.200$Ct.Mean[b6.200$Well=="H2"] <- b6.200$Ct[b6.200$Well=="H2"]

# t100 - GLi3-1 and MEGF9-1        Not good: A2, A3, A4, A5, A6, A7, A8, F12, D11, D12              
#Hmm. Lots of weird things there. Fortunately always one technical good. I suspect that the new plate seals aren't the best
b7.100 <- read.csv("~/Plate25_100_GLi3-1MEGF9-1.csv", header=TRUE)
#because there were some "Undetermined" Ct values, here Ct isn't numeric. Need to convert it
b7.100$Ct <- as.numeric(as.character(b7.100$Ct))

b7.100 <- b7.100[b7.100$Well!="A2",]
b7.100$Ct.Mean[b7.100$Well=="B2"] <- b7.100$Ct[b7.100$Well=="B2"]
b7.100 <- b7.100[b7.100$Well!="A3",]
b7.100$Ct.Mean[b7.100$Well=="B3"] <- b7.100$Ct[b7.100$Well=="B3"]
b7.100 <- b7.100[b7.100$Well!="A4",]
b7.100$Ct.Mean[b7.100$Well=="B4"] <- b7.100$Ct[b7.100$Well=="B4"]
b7.100 <- b7.100[b7.100$Well!="A5",]
b7.100$Ct.Mean[b7.100$Well=="B5"] <- b7.100$Ct[b7.100$Well=="B5"]
b7.100 <- b7.100[b7.100$Well!="A6",]
b7.100$Ct.Mean[b7.100$Well=="B6"] <- b7.100$Ct[b7.100$Well=="B6"]
b7.100 <- b7.100[b7.100$Well!="A7",]
b7.100$Ct.Mean[b7.100$Well=="B7"] <- b7.100$Ct[b7.100$Well=="B7"]
b7.100 <- b7.100[b7.100$Well!="A8",]
b7.100$Ct.Mean[b7.100$Well=="B8"] <- b7.100$Ct[b7.100$Well=="B8"]
b7.100 <- b7.100[b7.100$Well!="F12",]
b7.100$Ct.Mean[b7.100$Well=="E12"] <- b7.100$Ct[b7.100$Well=="E12"]
b7.100 <- b7.100[b7.100$Well!="D11",]
b7.100$Ct.Mean[b7.100$Well=="C11"] <- b7.100$Ct[b7.100$Well=="C11"]
b7.100 <- b7.100[b7.100$Well!="D12",]
b7.100$Ct.Mean[b7.100$Well=="C12"] <- b7.100$Ct[b7.100$Well=="C12"]

# t150 - GLi3-1 and MEGF9-1          H4 not good
b7.150 <- read.csv("~/Plate24_150_GLi3-1MEGF9-1.csv", header=TRUE)
b7.150 <- b7.150[b7.150$Well!="H4",]
b7.150$Ct.Mean[b7.150$Well=="G4"] <- b7.150$Ct[b7.150$Well=="G4"]

# t200 - GLi3-1 and MEGF9-1          
b7.200 <- read.csv("~/Plate23_200_GLi3-1MEGF9-1.csv", header=TRUE)

# t100 - MAGUK-2 and Rhoguanin-1     A4 and A5 not good
b8.100 <- read.csv("~/Plate28_100_MAGUK-2Rhoguanin-1.csv", header=TRUE)
b8.100 <- b8.100[b8.100$Well!="A4",]
b8.100$Ct.Mean[b8.100$Well=="B4"] <- b8.100$Ct[b8.100$Well=="B4"]
b8.100 <- b8.100[b8.100$Well!="A5",]
b8.100$Ct.Mean[b8.100$Well=="B5"] <- b8.100$Ct[b8.100$Well=="B5"]

# t150 - MAGUK-2 and Rhoguanin-1      A4 not good
b8.150 <- read.csv("~/Plate27_150_MAGUK-2Rhoguanin-1.csv", header=TRUE) 
b8.150 <- b8.150[b8.150$Well!="A4",]
b8.150$Ct.Mean[b8.150$Well=="B4"] <- b8.150$Ct[b8.150$Well=="B4"]

# t200 - MAGUK-2 and Rhoguanin-1      A4 and A5 not good
b8.200 <- read.csv("~/Plate26_200_MAGUK-2Rhoguanin-1.csv", header=TRUE)
b8.200 <- b8.200[b8.200$Well!="A4",]
b8.200$Ct.Mean[b8.200$Well=="B4"] <- b8.200$Ct[b8.200$Well=="B4"]
b8.200 <- b8.200[b8.200$Well!="A5",]
b8.200$Ct.Mean[b8.200$Well=="B5"] <- b8.200$Ct[b8.200$Well=="B5"]

# Lists of Morph=Sample 
# order of morphs in the samples for the 100 timepoint (s25,s26,s27,s28,s29,s30,s31,s32,s33,s34,s35,s36)
s100 <- c("PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB")
# order of morphs in the samples for the 150 timepoint - second batch (13,14,15,16,17,18,19,20,21,22,23,24)
s150 <- c("PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB")
# order of morphs in the samples for the 200 timepoint (s1,s10,s11,s12,s2,s3,s4,s5,s6,s7,s8,s9)
s200 <- c("PL","SB","SB","SB","PL","PL","LB","LB","LB","PI","PI","PI")

# timepoints
timepoint100 <- "t100"
timepoint150 <- "t150"
timepoint200 <- "t200"

## Function to extract certain columns and values from the tables
# Compared to the previous qPCRs I did in 2019, they seemed to have changed "Cq.Mean" into "Ct.Mean"

tafla <- function(t,s,c) {
  t <- data.frame(t$Sample,t$Target,t$Ct.Mean)
  colnames(t) <- c("Sample","Target","Ct.Mean")
  t <- t[!is.na(t$Ct.Mean),]
  t <- t %>% distinct()
  t1 <- dcast(t, t$Sample ~ t$Target, value.var="Ct.Mean")
  colnames(t1)[1]<-"Sample"
  t1$Reference <- sqrt(t1$actb*t1$Ub2l3)
  t1$Morph <- s
  t1$Timepoint <- c
  return(t1)
}

b1.100a <- tafla(b1.100,s100,timepoint100)
b1.150a <- tafla(b1.150,s150,timepoint150)
b1.200a <- tafla(b1.200,s200,timepoint200)

b2.100a <- tafla(b2.100,s100,timepoint100)
b2.150a <- tafla(b2.150,s150,timepoint150)
b2.200Myela <- tafla(b2.200Myel,s200,timepoint200)
b2.200SLCa <- tafla(b2.200SLC,s200,timepoint200)

b3.100a <- tafla(b3.100,s100,timepoint100)
b3.150a <- tafla(b3.150,s150,timepoint150)
b3.200a <- tafla(b3.200,s200,timepoint200)

b4.100a <- tafla(b4.100,s100,timepoint100)
b4.150a <- tafla(b4.150,s150,timepoint150)
b4.200a <- tafla(b4.200,s200,timepoint200)

b5.100a <- tafla(b5.100,s100,timepoint100)
b5.150a <- tafla(b5.150,s150,timepoint150)
b5.200a <- tafla(b5.200,s200,timepoint200)

b6.100a <- tafla(b6.100,s100,timepoint100)
b6.150a <- tafla(b6.150,s150,timepoint150)
b6.200a <- tafla(b6.200,s200,timepoint200)

b7.100a <- tafla(b7.100,s100,timepoint100)
b7.150a <- tafla(b7.150,s150,timepoint150)
b7.200a <- tafla(b7.200,s200,timepoint200)

b8.100a <- tafla(b8.100,s100,timepoint100)
b8.150a <- tafla(b8.150,s150,timepoint150)
b8.200a <- tafla(b8.200,s200,timepoint200)

#################### Specifically make recaps so that it is easier to use #############################
#do not bother with the $Reference column. It will not be used later in the script.
b1 <- rbind(b1.100a, b1.150a, b1.200a)

#Need to tweak things a bit for b2
# First, we don't need b2.200Myel. 
# Then we need b2.200SLC, but the column names are different, so I need to change that.
# This is scuffed but I will change the title of the ZFP2 column into "MyelTF1" so that rbind works. I am not using this ZFP2/MyelTF1 data anyway.
colnames(b2.200SLCa) <- c("Sample","actb","SLC9A3R2","Ub2l3","MyelTF1","Reference", "Morph","Timepoint")
b2 <- rbind(b2.100a, b2.150a, b2.200SLCa)

b3 <- rbind(b3.100a, b3.150a, b3.200a)
b4 <- rbind(b4.100a, b4.150a, b4.200a)
b5 <- rbind(b5.100a, b5.150a, b5.200a)
b6 <- rbind(b6.100a, b6.150a, b6.200a)
b7 <- rbind(b7.100a, b7.150a, b7.200a)
b8 <- rbind(b8.100a, b8.150a, b8.200a)


###################################################################################################
#Here, I use the model from Hendelmans et al., which is derived from the Pfaffl equation to use multiple ref genes:

        # b1.100a <- readRDS("/Users/sebmatlosz/Desktop/b1_100a.rds")
        # b1.150a <- readRDS("/Users/sebmatlosz/Desktop/b1_150a.rds")
        # b1.200a <- readRDS("/Users/sebmatlosz/Desktop/b1_200a.rds")
        # 
        # b1 <- rbind(b1.100a, b1.150a, b1.200a)

#First, calculate deltaCt (control - samples). For both of the reference genes and the goi.
#Let's start by choosing one of the PI100 as a control. We can decide later whether we want it to be a specific one or not.
#Nevermind, even better, take the average mean of all values. Wait but then I would need to take the averages for each gene so it is somehow unrelated?


#Here, x is the data frame, Eff1 is the Efficiency of the first gene, Eff2 the efficiency of the second gene
CalculateREr <- function(x,Eff1,Eff2) {

    #Keep gene names in variable.
    gene1 <- colnames(x)[3]
    gene2 <- colnames(x)[4]

    #Chose a sample of Reference (here S32)
    actb_ctrl <- x[,2][x$Sample=="S32"]
    Ub2l3_ctrl <- x[,5][x$Sample=="S32"]
    gene1_ctrl <- x[,3][x$Sample=="S32"]
    gene2_ctrl <- x[,4][x$Sample=="S32"]

    #Calculate the first delta by substracting each sample from the sample of reference for each gene.
    x$delta_actb <-  actb_ctrl - x[,2]   
    x$delta_Ub2l3 <-  Ub2l3_ctrl - x[,5]
    x[,paste("delta",gene1,sep="")] <-  gene1_ctrl - x[,3]
    x[,paste("delta",gene2,sep="")] <-  gene2_ctrl - x[,4]


    #Calculate relative quantity values for each gene. RQ=E^deltaCt
    Effactb <- 1.95
    EffUb2l3 <- 1.93
    
    x$RQ_actb <-  Effactb ^ x$delta_actb
    x$RQ_Ub2l3 <-  EffUb2l3 ^ x$delta_Ub2l3
    x[,paste("RQ_",gene1,sep="")] <- Eff1 ^ x[,paste("delta",gene1,sep="")]
    x[,paste("RQ_",gene2,sep="")] <- Eff2 ^ x[,paste("delta",gene2,sep="")]
    
    #Calculate the geometric mean of the reference genes RQ values
    x$GEOMEAN_REF <- sqrt(x$RQ_actb*x$RQ_Ub2l3)
    
    #Calculate relative gene expression values
    x[,paste("RE_",gene1,sep="")] <- x[,paste("RQ_",gene1,sep="")] / x$GEOMEAN_REF
    x[,paste("RE_",gene2,sep="")] <- x[,paste("RQ_",gene2,sep="")] / x$GEOMEAN_REF
    
return(x)
    }


b1_bis <- CalculateREr(b1,1.99,1.9) #HiH2A, Lmtk2
b2_bis <- CalculateREr(b2,2.0,1.95) #The first one doesn't matter, it is MyelTF1, then SLC9A3R2
b3_bis <- CalculateREr(b3,2.0,1.93)   #The first one doesn't matter, it is Ets2, then Nkx23
b4_bis <- CalculateREr(b4,2.17,2.21) #NFIX, RASSF4
b5_bis <- CalculateREr(b5,2.49,2.2) #ARL16, Meis
b6_bis <- CalculateREr(b6,1.92,2.17) #ARMC1, HiH3l
b7_bis <- CalculateREr(b7,2.06,1.95) #Gli3, MEGF9
b8_bis <- CalculateREr(b8,2.02,1.92) #MAGUK, Rhoguanin

#Make recaps out of these tables for each gene.
REr.Lmtk2 <- b1_bis[,-c(3,11,15,18)]
REr.HiH2A <- b1_bis[,-c(4,12,16,19)]
REr.SLC9A3R2 <- b2_bis[,-c(3,11,15,18)]
REr.Nkx232 <- b3_bis[,-c(3,11,15,18)]
REr.NFiX1 <- b4_bis[,-c(4,12,16,19)]
REr.RASSF42 <- b4_bis[,-c(3,11,15,18)]
REr.ARF161 <- b5_bis[,-c(4,12,16,19)]
REr.MEiS1l2 <- b5_bis[,-c(3,11,15,18)]
REr.HiH3l2 <- b6_bis[,-c(3,11,15,18)]
REr.ARMC11 <- b6_bis[,-c(4,12,16,19)]
REr.GLi31 <- b7_bis[,-c(4,12,16,19)]
REr.MEGF91 <- b7_bis[,-c(3,11,15,18)]
REr.MAGUK2 <- b8_bis[,-c(4,12,16,19)]
REr.Rhoguanin1 <- b8_bis[,-c(3,11,15,18)]

### Just load the data if I want to skip everything ###

REr.Lmtk2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_Lmtk2.rds")
REr.HiH2A <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_HiH2A.rds")
REr.SLC9A3R2 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_SLC9A3R2.rds")
REr.Nkx232 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_Nkx23.rds")
REr.NFiX1 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_NFIX.rds")
REr.RASSF42 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_RASSF4.rds")
REr.ARF161 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_ARL16.rds")
REr.MEiS1l2 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_MEIS1.rds")
REr.HiH3l2 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_H3-like.rds")
REr.ARMC11 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_ARMC1.rds")
REr.GLi31 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_Gli3.rds")
REr.MEGF91 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_MEGF9.rds")
REr.MAGUK2 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_MAGUK2.rds")
REr.Rhoguanin1 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_Rhoguanin.rds")

#### Check normality with Shapiro
library(dplyr)
library(ggpubr)

shapiro.test(REr.Rhoguanin1$`RE_Rhoguanin-1`)$p.value

#Transform the ones that are not normal.
lambdaValue <- BoxCox.lambda(REr.ARF161$`RE_ARF16-1`, method = "guerrero")
REr.ARF161$transformed <- BoxCox(REr.ARF161$`RE_ARF16-1`, lambdaValue)

lambdaValue <- BoxCox.lambda(REr.HiH2A$RE_HiH2A, method = "guerrero")
REr.HiH2A$transformed <- BoxCox(REr.HiH2A$RE_HiH2A, lambdaValue)

lambdaValue <- BoxCox.lambda(REr.Lmtk2$RE_Lmtk2, method = "guerrero")
REr.Lmtk2$transformed <- BoxCox(REr.Lmtk2$RE_Lmtk2, lambdaValue)

lambdaValue <- BoxCox.lambda(REr.MAGUK2$`RE_MAGUK-2`, method = "guerrero")
REr.MAGUK2$transformed <- BoxCox(REr.MAGUK2$`RE_MAGUK-2`, lambdaValue)

lambdaValue <- BoxCox.lambda(REr.NFiX1$`RE_NFiX-1`, method = "guerrero")
REr.NFiX1$transformed <- BoxCox(REr.NFiX1$`RE_NFiX-1`, lambdaValue)

lambdaValue <- BoxCox.lambda(REr.SLC9A3R2$RE_SLC9A3R2, method = "guerrero")
REr.SLC9A3R2$transformed <- BoxCox(REr.SLC9A3R2$RE_SLC9A3R2, lambdaValue)

shapiro.test(REr.ARF161$`RE_ARF16-1`)$p.value





#Make grahs of the normal distribution of these RErs
     # HIH2Agraph <- ggplot(REr.HiH2A, aes(sample = transformed)) + stat_qq() + stat_qq_line() + theme_bw() + 
     #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("H2A-like REr")
     # LMTK2graph <- ggplot(REr.Lmtk2, aes(sample = transformed)) + stat_qq() + stat_qq_line() + theme_bw() + 
     #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("LMTK2 REr")
     # SLC9A3R2graph <- ggplot(REr.SLC9A3R2, aes(sample = transformed)) + stat_qq() + stat_qq_line() + theme_bw() + 
     #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("SLC9A3R2-like REr")
    # NKX23graph <- ggplot(b3_bis, aes(sample = `RE_Nkx23-2`)) + stat_qq() + stat_qq_line() + theme_bw() + 
    #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("NKX23-like REr")
     # NFIXgraph <- ggplot(REr.NFiX1, aes(sample = transformed)) + stat_qq() + stat_qq_line() + theme_bw() + 
     #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("NFIX REr")
    # RASSF4graph <- ggplot(b4_bis, aes(sample = `RE_RASSF4-2`)) + stat_qq() + stat_qq_line() + theme_bw() + 
    #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("RASSF4-like REr")
    # Meisgraph <- ggplot(b5_bis, aes(sample = `RE_MEiS1l-2`)) + stat_qq() + stat_qq_line() + theme_bw() + 
    #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("MEIS1-like REr")
     # ARL16graph <- ggplot(REr.ARF161, aes(sample = transformed)) + stat_qq() + stat_qq_line() + theme_bw() + 
     #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("ARL16 REr")
    # ARMC1graph <- ggplot(b6_bis, aes(sample = `RE_ARMC1-1`)) + stat_qq() + stat_qq_line() + theme_bw() + 
    #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("ARMC1 REr")
    # H3graph <- ggplot(b6_bis, aes(sample = `RE_HiH3l-2`)) + stat_qq() + stat_qq_line() + theme_bw() + 
    #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("H3-like REr")
    # GLI3graph <- ggplot(b7_bis, aes(sample = `RE_GLi3-1`)) + stat_qq() + stat_qq_line() + theme_bw() + 
    #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("GLI3-like REr")
    # MEGF9graph <- ggplot(b7_bis, aes(sample = `RE_MEGF9-1`)) + stat_qq() + stat_qq_line() + theme_bw() + 
    #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("MEGF9-like REr")
    # Rhoguagraph <- ggplot(b8_bis, aes(sample = `RE_Rhoguanin-1`)) + stat_qq() + stat_qq_line() + theme_bw() + 
    #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("ARHGEF37-like REr")
     # MPP3graph <- ggplot(REr.MAGUK2, aes(sample = transformed)) + stat_qq() + stat_qq_line() + theme_bw() + 
     #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("MPP3 REr")
    # 
    # 
     # figure <- ggarrange(LMTK2graph,SLC9A3R2graph,NFIXgraph,ARL16graph,HIH2Agraph,MPP3graph, ncol=3, nrow=2)
    # figure <- ggarrange(LMTK2graph,ARMC1graph,SLC9A3R2graph,NKX23graph,NFIXgraph,RASSF4graph,ARL16graph,Meisgraph,H3graph,
    #                     HIH2Agraph,GLI3graph,MEGF9graph,MPP3graph,Rhoguagraph, ncol=3, nrow=5)
    # 
     # pdf("/Users/sebmatlosz/Desktop/AllFilesForSubmission/qqNorm_RER_posttransformation_210222.pdf", 10,6)
     # figure
     # dev.off()

      #Since I am at it, check the normality of dCt and ddCt data as well?
        # ddCt.Lmtk2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-Lmtk2.rds")
        # ddCt.HiH2A <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-HiH2A.rds")
        # ddCt.SLC9A3R2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-SLC9A3R2.rds")
        # ddCt.Nkx232 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-Nkx232.rds")
        # ddCt.NFiX1 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-NFiX1.rds")
        # ddCt.RASSF42 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-RASSF42.rds")
        # ddCt.ARF161 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-ARF161.rds")
        # ddCt.MEiS1l2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-MEiS1l2.rds")
        # ddCt.HiH3l2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-HiH3l2.rds")
        # ddCt.ARMC11 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-ARMC11.rds")
        # ddCt.GLi31 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-GLi31.rds")
        # ddCt.MEGF91 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-MEGF91.rds")
        # ddCt.MAGUK2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-MAGUK2.rds")
        # ddCt.Rhoguanin1 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/noNorm/ddCt-Rhoguanin1.rds")
        # 
        # #Check for normality with shapiro:
        # shapiro.test(ddCt.Rhoguanin1$ddCt)$p.value

        # 
        # #Change aes(sample = dCt) to aes(sample = ddCt)  
        # HIH2Agraph <- ggplot(ddCt.HiH2A, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() + 
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("H2A-like ddCt")
        # LMTK2graph <- ggplot(ddCt.Lmtk2, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("LMTK2 ddCt")
        # SLC9A3R2graph <- ggplot(ddCt.SLC9A3R2, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("SLC9A3R2-like ddCt")
        # NKX23graph <- ggplot(ddCt.Nkx232, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("NKX23-like ddCt")
        # NFIXgraph <- ggplot(ddCt.NFiX1, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("NFIX ddCt")
        # RASSF4graph <- ggplot(ddCt.RASSF42, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("RASSF4-like ddCt")
        # Meisgraph <- ggplot(ddCt.MEiS1l2, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("MEIS1-like ddCt")
        # ARL16graph <- ggplot(ddCt.ARF161, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("ARL16 ddCt")
        # ARMC1graph <- ggplot(ddCt.ARMC11, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("ARMC1 ddCt")
        # H3graph <- ggplot(ddCt.HiH3l2, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("H3-like ddCt")
        # GLI3graph <- ggplot(ddCt.GLi31, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("GLI3-like ddCt")
        # MEGF9graph <- ggplot(ddCt.MEGF91, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("MEGF9-like ddCt")
        # Rhoguagraph <- ggplot(ddCt.Rhoguanin1, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("ARHGEF37-like ddCt")
        # MPP3graph <- ggplot(ddCt.MAGUK2, aes(sample = ddCt)) + stat_qq() + stat_qq_line() + theme_bw() +
        #   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("MPP3 ddCt")
        # 
        # 
        # figure <- ggarrange(LMTK2graph,ARMC1graph,SLC9A3R2graph,NKX23graph,NFIXgraph,RASSF4graph,ARL16graph,Meisgraph,H3graph,
        #                     HIH2Agraph,GLI3graph,MEGF9graph,MPP3graph,Rhoguagraph, ncol=3, nrow=5)
        # 
        # pdf("/Users/sebmatlosz/Desktop/AllFilesForSubmission/qqNorm_ddCt_160222.pdf", 10,15)
        # figure
        # dev.off()








################################################################
#Plot these relative expression ratios.
color <- c("#00BA38FF","#619CFFFF","#F8766DFF","darkmagenta")
names(color) <- c("LB","SB","PL","PI")


#Use labs(x= "", y = "Relative expression") for ARL16
#And labs(x = "Developmental stage (ts)", y = "") for Rhoguanin
onebar <- function(df,s){
  df$morph <- factor(df$morph)
  df$morph <- relevel(df$morph,"PI")
  df$morph <- relevel(df$morph,"PL")
  df$morph <- relevel(df$morph,"SB")
  df$morph <- relevel(df$morph,"LB")
  
  
  graph <- ggplot(df, aes(x=df$time, y=df$mean, fill=df$morph)) + 
    geom_bar(stat="identity", position = position_dodge()) +
    labs(title=s,x="Developmental stage (ts)",y="", fill="Morph") +
    theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
    scale_fill_manual(values=color)+
    geom_errorbar(aes(ymin=df$mean-df$sd, ymax=df$mean+df$sd), width=.5,position=position_dodge(.9))
}

bargraphddCtbetter <- function(t,s,column) {

  t100 <- t[t$Timepoint=="t100",]
  t150 <- t[t$Timepoint=="t150",]
  t200 <- t[t$Timepoint=="t200",]
  
  # 100
  meangene100 <- tapply(t100[,column], t100$Morph, mean)
  sdgene100 <- tapply(t100[,column], t100$Morph, sd)
  df100<- data.frame(meangene100,sdgene100)
  df100$morph <- row.names(df100)
  colnames(df100)<- c("mean","sd","morph")
  df100$time <- "100"
  
  # 150
  meangene150 <- tapply(t150[,column], t150$Morph, mean)
  sdgene150 <- tapply(t150[,column], t150$Morph, sd)
  df150 <- data.frame(meangene150,sdgene150)
  df150$morph <- row.names(df150)
  colnames(df150)<- c("mean","sd","morph")
  df150$time <- "150"
  
  # 200
  meangene200 <- tapply(t200[,column], t200$Morph, mean)
  sdgene200 <- tapply(t200[,column], t200$Morph, sd)
  df200 <- data.frame(meangene200,sdgene200)
  df200$morph <- row.names(df200)
  colnames(df200)<- c("mean","sd","morph")
  df200$time <- "200"
  
  dftotal <- rbind(df100,df150,df200)
  graphtotal <- onebar(dftotal,s)
  
}

bar.Lmtk2 <- bargraphddCtbetter(b1_bis,"LMTK2","RE_Lmtk2")
bar.HiH2A <- bargraphddCtbetter(b1_bis,"H2A-like","RE_HiH2A")
bar.SLC9A3R2 <- bargraphddCtbetter(b2_bis,"SLC9A3R2-like","RE_SLC9A3R2")
bar.Nkx23 <- bargraphddCtbetter(b3_bis,"NKX23-like","RE_Nkx23-2")
bar.NFIX <- bargraphddCtbetter(b4_bis,"NFIX","RE_NFiX-1")
bar.RASSF4 <- bargraphddCtbetter(b4_bis,"RASSF4-like","RE_RASSF4-2")
bar.ARL16 <- bargraphddCtbetter(b5_bis,"ARL16","RE_ARF16-1")
bar.MEIS <- bargraphddCtbetter(b5_bis,"MEIS1-like","RE_MEiS1l-2")
bar.ARMC1 <- bargraphddCtbetter(b6_bis,"ARMC1","RE_ARMC1-1")
bar.HiH3l <- bargraphddCtbetter(b6_bis,"H3-like","RE_HiH3l-2")
bar.Gli3 <- bargraphddCtbetter(b7_bis,"GLI3-like","RE_GLi3-1")
bar.MEGF9 <- bargraphddCtbetter(b7_bis,"MEGF9","RE_MEGF9-1")
bar.MAGUK <- bargraphddCtbetter(b8_bis,"MPP3","RE_MAGUK-2")
bar.Rhoguanin <- bargraphddCtbetter(b8_bis,"ARHGEF37-like","RE_Rhoguanin-1")

#Extract legend from one graph to use it as an object.
leg <- get_legend(bar.Lmtk2)

#Remove the legend from all plots
removelegend <- function(x) {
  y <- x + theme(legend.position = "none")
  return(y)
}

bar.Lmtk2 <- removelegend(bar.Lmtk2)
bar.HiH2A <- removelegend(bar.HiH2A)
bar.SLC9A3R2 <- removelegend(bar.SLC9A3R2)
bar.Nkx232 <- removelegend(bar.Nkx23)
bar.NFiX1 <- removelegend(bar.NFIX)
bar.RASSF42 <- removelegend(bar.RASSF4)
bar.ARL161 <- removelegend(bar.ARL16)
bar.MEiS1l2 <- removelegend(bar.MEIS)
bar.HiH3l2 <- removelegend(bar.HiH3l)
bar.ARMC11 <- removelegend(bar.ARMC1)
bar.GLi31 <- removelegend(bar.Gli3)
bar.MEGF91 <- removelegend(bar.MEGF9)
bar.MAGUK2 <- removelegend(bar.MAGUK)
bar.Rhoguanin1 <- removelegend(bar.Rhoguanin)


figure <- ggarrange(bar.Lmtk2,bar.ARMC11,bar.SLC9A3R2,bar.Nkx232,bar.NFiX1,bar.RASSF42,bar.ARL161,bar.MEiS1l2,bar.HiH3l2,bar.HiH2A,
                    bar.GLi31,bar.MEGF91,bar.MAGUK2,bar.Rhoguanin1, leg, ncol=3,nrow=5)

figure


pdf("/Users/sebmatlosz/Desktop/AllFilesForSubmission/qPCR_EffNorm_190122.pdf", 10,15)
figure
dev.off()


############################################################################################
#Pfaffl equation  (OLD)
# b1[,paste("RE_",gene1,sep="")] <- ((Eff1 ^ b1$delta_HiH2A)/(EffRef ^ b1$delta_ref))
# b1[,paste("RE_",gene2,sep="")] <- ((Eff2 ^ b1$delta_Lmtk2)/(EffRef ^ b1$delta_ref))
