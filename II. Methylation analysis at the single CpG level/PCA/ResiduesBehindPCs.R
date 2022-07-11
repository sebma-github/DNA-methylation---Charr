### This script is for looking into which residues have the most weight in specific Principal Components
#meth.min1 objects generated using the BentvsLimn.R script
library(ggplot2)
library(stringr)
library(dplyr)

#Read in the table
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

df_residue <- as.data.frame(df_PCA$rotation) #This should give back the value of each residue for each PCA

# df_residue$ID <- row.names(df_residue)
# write.csv(df_residue,"~/results/weightsonPC_48samplesPCA.csv", row.names = FALSE)


############################### Keep only the CpGs that are outside the two sigmas for each distribution ################
PCAfinal <- read.csv("~/results/weightsonPC_48samplesPCA.csv")

#PC1
    min <- mean(PCAfinal$PC1) - 2*sd(PCAfinal$PC1)
    max <- mean(PCAfinal$PC1) + 2*sd(PCAfinal$PC1)
    PC1outliers <- filter(PCAfinal, PC1 < min | PC1 > max)

#PC4 
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


