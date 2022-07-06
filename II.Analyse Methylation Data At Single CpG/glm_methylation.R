####### glm for my methylation data#######
library(stringr)
library(dplyr)

#Read in the file with methylation info for the 10340 CpGs that have >10X coverage in all samples
meth.min1df <- read.table("~/methmin1_48samples.csv", sep=",", header=TRUE)

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
      # totalt$Morph <- c("Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent",
      #                  "Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent",
      #                  "Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent",
      #                  "Limn","Limn","Limn","Bent","Bent","Bent","Limn","Limn","Limn","Bent","Bent","Bent")
      totalt$Morph <- c("PL","PL","PL","LB","LB","LB","PI","PI","PI","SB","SB","SB",
                         "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB",
                         "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB",
                         "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB")
      totalt$Timepoint <- c("200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts",
                        "150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts",
                        "100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts",
                        "50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts")
      totalt$Sex <- c("F","M","M","M","F","M","M","M","M","M","F","F","F",NA,"F","M","M","F",
                      "M","F","M","F","M","F","F","F","M","F","M","M","M","F","F","M",NA,"F",
                      "F","F",NA,NA,NA,NA,NA,NA,"F",NA,NA,NA)


#Make glms 
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
#write.csv(dfFinal,"~/glm_10340res_pvalues.csv", row.names = FALSE)

#Reformat the annotations
        names <- as.data.frame(str_split_fixed(dfFinal$Residue,"_",2))
        dfFinal$chrom <- names$V1
        dfFinal$start <- names$V2

        dfFinal$temp <- gsub("chr","",dfFinal$chrom)
        dfFinal$NW <- ifelse(nchar(dfFinal$temp)==11, "NW_", "NC_")
        dfFinal$NCBI <- paste(dfFinal$NW,dfFinal$temp, sep="")

#Only keep residues that have a significant p.value after Bonferroni correction
#First correct with Bonferroni
        v$Morph_pvalue_corrected <- p.adjust(v$Morph_pvalue, method = "bonferroni", n = length(v$Morph_pvalue))
        v$Timepoint_pvalue_corrected <- p.adjust(v$Timepoint_pvalue, method = "bonferroni", n = length(v$Timepoint_pvalue))
        v$Sex_pvalue_corrected <- p.adjust(v$Sex_pvalue, method = "bonferroni", n = length(v$Sex_pvalue))
        v$MorphXTimepoint_pvalue_corrected <- p.adjust(v$MorphXTimepoint_pvalue, method = "bonferroni", n = length(v$MorphXTimepoint_pvalue))

#Only keep the CpGs with corrected p.value < 0.05
        signifMorph <- filter(v,Morph_pvalue_corrected <= 0.05)
        signifTime <- filter(v,Timepoint_pvalue_corrected <= 0.05)
        signifSex <- filter(v,Sex_pvalue_corrected <= 0.05)
        signifMorphXTime <- filter(v,MorphXTimepoint_pvalue_corrected <= 0.05)

        
# write.csv(signifMorph,"/Users/sebmatlosz/Desktop/methylation/glm_signif_Morphs_090222.csv", row.names = FALSE)
# write.csv(signifTime,"/Users/sebmatlosz/Desktop/methylation/glm_signif_Timepoints_090222.csv", row.names = FALSE)
# write.csv(signifSex,"/Users/sebmatlosz/Desktop/methylation/glm_signif_Sex_090222.csv", row.names = FALSE)
# write.csv(signifMorphXTime,"/Users/sebmatlosz/Desktop/methylation/glm_signif_MorphXTime_090222.csv", row.names = FALSE)
        


