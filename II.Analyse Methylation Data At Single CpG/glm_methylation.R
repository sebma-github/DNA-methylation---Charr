####### glm for my methylation data#######
#### So I need to have 10340 columns + 1 column morph (or can I use just a vector of the morphs?)####

meth.min1df <- read.table("/Users/sebmatlosz/Documents/Desktop_methwork/methylation/methmin1_48samples.csv", sep=",", header=TRUE)

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

#Transpose the data so that it is in the right way.
totalt <- as.data.frame(t(total))

#Add a column for morph, timepoint and sex.
#In the Morph column, choose whether you want to specify the Niche or the Morphs specifically
totalt$Morph <- c("Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent",
                  "Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent",
                  "Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent",
                  "Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent")
# totalt$Morph <- c("PL","PL","PL","LB","LB","LB","PI","PI","PI","SB","SB","SB",
#                   "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB",
#                   "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB",
#                   "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB")
totalt$Timepoint <- c("200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts",
                  "150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts",
                  "100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts",
                  "50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts")
totalt$Sex <- c("F","M","M","M","F","M","M","M","M","M","F","F","F",NA,"F","M","M","F",
                "M","F","M","F","M","F","F","F","M","F","M","M","M","F","F","M",NA,"F",
                "F","F",NA,NA,NA,NA,NA,NA,"F",NA,NA,NA)




#It should be something such as lm(methylationres1 ~ morph + timepoint + sex)
#So the easiest is probably to have a df with 10340 columns + 3 columns for these variables.

#To write the residue string in the formula, use get()
#Takes about 5 minutes of running when only taking the pvalues out.
######################################################
#First function, with glms with Sex
glmallresidues <- function() {
  df <- data.frame(matrix(ncol = 7, nrow = 0))
  colnames(df) <- c("Residue","Morph_pvalue","Timepoint_pvalue","Sex_pvalue","MorphXTimepoint_pvalue","MorphXSex_pvalue","TimepointXSex_pvalue")
  
  for (i in 1:10340) {
    residue <- colnames(totalt)[i]
    glmmodel <- glm(get(residue) ~ Morph + Timepoint + Sex + Morph * Timepoint + Morph * Sex + Timepoint * Sex, data=totalt, family = gaussian())
    A <- anova(glmmodel, test="Chisq")
    
    df2 <- data.frame("Residue"=residue,"Morph_pvalue"=A$`Pr(>Chi)`[2],"Timepoint_pvalue"=A$`Pr(>Chi)`[3],"Sex_pvalue"=A$`Pr(>Chi)`[4],
                      "MorphXTimepoint_pvalue"=A$`Pr(>Chi)`[5],"MorphXSex_pvalue"=A$`Pr(>Chi)`[6],"TimepointXSex_pvalue"=A$`Pr(>Chi)`[7])
    df <- rbind(df,df2)
  }
  return(df)
}

dfFinal <- glmallresidues()
#write.csv(dfFinal,"/Users/sebmatlosz/Desktop/methylation/glm_10340res_pvalues.csv", row.names = FALSE)


########################################################
#Second function, taking the Sex out of the equation
glmallresidues2 <- function() {
  df <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(df) <- c("Residue","Morph_pvalue","Timepoint_pvalue","MorphXTimepoint_pvalue")
  
  for (i in 1:10340) {
    residue <- colnames(totalt)[i]
    glmmodel <- glm(get(residue) ~ Morph + Timepoint + Morph:Timepoint, data=totalt, family = gaussian())
    A <- anova(glmmodel, test="Chisq")
    
    df2 <- data.frame("Residue"=residue,"Morph_pvalue"=A$`Pr(>Chi)`[2],"Timepoint_pvalue"=A$`Pr(>Chi)`[3],"MorphXTimepoint_pvalue"=A$`Pr(>Chi)`[4])
    df <- rbind(df,df2)
  }
  return(df)
}

dfFinal <- glmallresidues2()
write.csv(dfFinal,"/Users/sebmatlosz/Documents/Desktop_methwork/methylation/glm_10340res_pvalues_NoSex_gaussian_colon.csv", row.names = FALSE)

#Create a list of residues that are top 5% in each category
#dfFinal <- read.csv("/Users/sebmatlosz/Desktop/methylation/glm_10340res_pvalues.csv")
dfFinal <- read.csv("/Users/sebmatlosz/Documents/Desktop_methwork/methylation/glm_10340res_pvalues_NoSex_gaussian.csv")

library(stringr)
library(dplyr)

#Reformat the annotations
names <- as.data.frame(str_split_fixed(dfFinal$Residue,"_",2))
dfFinal$chrom <- names$V1
dfFinal$start <- names$V2

dfFinal$temp <- gsub("chr","",dfFinal$chrom)
dfFinal$NW <- ifelse(nchar(dfFinal$temp)==11, "NW_", "NC_")
dfFinal$NCBI <- paste(dfFinal$NW,dfFinal$temp, sep="")

#v <- dfFinal[,-c(8,10,11)]
v <- dfFinal[,-c(5,7,8)]

#Then, two different ways
        # FIRST POSSIBILITY
        #If I want to take the 5percent with the most p.value. 
        top5pMorph <- top_n(v, -517, Morph_pvalue)
        top5pTimepoint <- top_n(v, -517, Timepoint_pvalue)
        top5pSex <- top_n(v, -517, Sex_pvalue)
        
        # write.csv(top5pMorph,"/Users/sebmatlosz/Desktop/methylation/glm_top5p_Morphs.csv", row.names = FALSE)
        # write.csv(top5pTimepoint,"/Users/sebmatlosz/Desktop/methylation/glm_top5p_Timepoints.csv", row.names = FALSE)
        # write.csv(top5pSex,"/Users/sebmatlosz/Desktop/methylation/glm_top5p_Sex.csv", row.names = FALSE)
        
        
        #SECOND POSSIBILITY
        #If I want to take all residues that have a significant p.value after Bonferroni correction
        #First correct with Bonferroni
        v$Morph_pvalue_corrected1 <- p.adjust(v$Morph_pvalue, method = "bonferroni", n = length(v$Morph_pvalue))
        v$Timepoint_pvalue_corrected1 <- p.adjust(v$Timepoint_pvalue, method = "bonferroni", n = length(v$Timepoint_pvalue))
        v$Sex_pvalue_corrected1 <- p.adjust(v$Sex_pvalue, method = "bonferroni", n = length(v$Sex_pvalue))
        v$MorphXTimepoint_pvalue_corrected1 <- p.adjust(v$MorphXTimepoint_pvalue, method = "bonferroni", n = length(v$MorphXTimepoint_pvalue))

        
        #OR WITH Hochberg
        v$Morph_pvalue_corrected2 <- p.adjust(v$Morph_pvalue, method ="hochberg", n = length(v$Morph_pvalue))
        v$Timepoint_pvalue_corrected2 <- p.adjust(v$Timepoint_pvalue, method = "hochberg", n = length(v$Timepoint_pvalue))
        v$Sex_pvalue_corrected2 <- p.adjust(v$Sex_pvalue, method = "hochberg", n = length(v$Sex_pvalue))
        v$MorphXTimepoint_pvalue_corrected2 <- p.adjust(v$MorphXTimepoint_pvalue, method = "hochberg", n = length(v$MorphXTimepoint_pvalue))
        
        
        #Hochberg only adds adds 8 CpGs to the 848 signif between time. And nothing changes for morph or sex.
        
        
        
        
        # colnames(v) <- c("Residue","Niche_pvalue","Timepoint_pvalue","Sex_pvalue","NicheXTimepoint_pvalue","NicheXSex_pvalue","TimepointXSex_pvalue",
         #                "start","NCBI","Niche_pvalue_corrected","Timepoint_pvalue_corrected","Sex_pvalue_corrected","NicheXTimepoint_pvalue_corrected")
        #Save this df for further use
        # write.csv(v,"/Users/sebmatlosz/Desktop/methylation/glm_10340res_pvalues_corrected.csv", row.names = FALSE)
        # write.csv(v,"/Users/sebmatlosz/Desktop/methylation/glm_10340res_Niche_pvalues_corrected.csv", row.names = FALSE)
        
        #Only keep the CpGs with corrected p.value < 0.05
        signifMorph <- filter(v,Morph_pvalue_corrected1 <= 0.05)
        signifMorph2 <- filter(v,Morph_pvalue_corrected2 <= 0.05)
        
        signifTime <- filter(v,Timepoint_pvalue_corrected1 <= 0.05)
        signifSex <- filter(v,Sex_pvalue_corrected1 <= 0.05)
        signifMorphXTime <- filter(v,MorphXTimepoint_pvalue_corrected1 <= 0.05)
        
        signifTime2 <- filter(v,Timepoint_pvalue_corrected2 <= 0.05)
        signifSex2 <- filter(v,Sex_pvalue_corrected2 <= 0.05)
        signifMorphXTime2 <- filter(v,MorphXTimepoint_pvalue_corrected2 <= 0.05)
        
         # write.csv(signifMorph,"/Users/sebmatlosz/Desktop/methylation/glm_signif_Morphs_090222.csv", row.names = FALSE)
         # write.csv(signifTime,"/Users/sebmatlosz/Desktop/methylation/glm_signif_Timepoints_090222.csv", row.names = FALSE)
         # write.csv(signifSex,"/Users/sebmatlosz/Desktop/methylation/glm_signif_Sex_090222.csv", row.names = FALSE)
         # write.csv(signifMorphXTime,"/Users/sebmatlosz/Desktop/methylation/glm_signif_MorphXTime_090222.csv", row.names = FALSE)
        

         

        
        
        
        
        
        
        
        
         
         



##########################################################################################################################

      # For testing purposes
      # 
      # df <- data.frame(matrix(ncol = 7, nrow = 0))
      # colnames(df) <- c("Residue","Morph_pvalue","Timepoint_pvalue","Sex_pvalue","MorphXTimepoint_pvalue","MorphXSex_pvalue","TimepointXSex_pvalue")
      # 
      # firstresidue <- glm(chr000861.1_522 ~ Morph + Timepoint + Sex + Morph * Timepoint + Morph * Sex + Timepoint * Sex, data=totalt, family = gaussian())
      # A <- anova(firstresidue, test="Chisq")
      # 
      # residue <- colnames(totalt)[1]
      # 
      # 
      # glmmodel <- glm(get(residue) ~ Morph + Timepoint + Sex + Morph * Timepoint + Morph * Sex + Timepoint * Sex, data=totalt, family = gaussian())
      # 
      # potato <- data.frame("Residue"=residue,"Morph_pvalue"=A$`Pr(>Chi)`[2],"Timepoint_pvalue"=A$`Pr(>Chi)`[3],"Sex_pvalue"=A$`Pr(>Chi)`[4],
      #            "MorphXTimepoint_pvalue"=A$`Pr(>Chi)`[5],"MorphXSex_pvalue"=A$`Pr(>Chi)`[6],"TimepointXSex_pvalue"=A$`Pr(>Chi)`[7])
      # 
      # df <- rbind(df,potato)



