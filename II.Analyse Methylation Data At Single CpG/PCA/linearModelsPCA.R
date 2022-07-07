
###################################### Make Chi square tests on Repeat Elements #####################################
#Lets just make a table that looks like what is shown (Need number of elements for Xsquare to work)
M <- as.table(rbind(c(170,5519,2588,4419,185,70,5268,1281), c(4,16,18,10,9,0,28,8),c(111,3649,1798,3017,88,30,3553,899),c(3,56,79,56,33,2,134,23)))
dimnames(M) <- list(regiontype = c("Non-DMRMorph", "DMRMorph","Non-DMRStage","DMRStage"),
                    elementtype = c("RT-SINES","RT-LINES", "RT-LTR","Transp","sRNA","Sat","SR","LC"))
M
N <- M[1:2,1:3]
Xsq <- chisq.test(N)
fisher.test(N)
Xsq
Xsq$observed
Xsq$expected



#This is not going to be used
                  df <- as.data.frame(rbind(NonDMRMorph,DMRMorph,NonDMRStage,DMRStage))
                  colnames(df) <- c("RT-SINES","RT-LINES", "RT-LTR","Transp","sRNA","Sat","SR","LC")
                  df$total <- rowSums(df)
                  
                  NonDMRMorphPerc <- NonDMRMorph *100 /19500
                  DMRMorphPerc <- DMRMorph * 100 / 93
                  NonDMRStagePerc <- NonDMRStage * 100 / 13145
                  DMRStagePerc <- DMRStage * 100 / 386
                  df2 <- as.data.frame(rbind(NonDMRMorphPerc,DMRMorphPerc,NonDMRStagePerc,DMRStagePerc))
                  colnames(df2) <- c("RT-SINES","RT-LINES", "RT-LTR","Transp","sRNA","Sat","SR","LC")

#I don't understand how chi squares work.

###################################################### PCA methylome charr #############################################################
#meth.min1 objects generated using the BentvsLimn.R script
#library(methylKit)
library(ggplot2)

                                              #######For PCA on 48 samples
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

#Should invert the df to plot based on samples. Otherwise I plot the entire 10 000 residues.
totalt <- t(total)

df_PCA <- prcomp(totalt)
df_out <- as.data.frame(df_PCA$x)
#This was mainly for curiosity
df_out$morph <- c("PL","PL","PL","LB","LB","LB","PI","PI","PI","SB","SB","SB",
                  "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB",
                  "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB",
                  "PL","PL","PI","LB","LB","SB","PI","PI","PL","SB","SB","LB")
df_out$niche <- c("Limnetic","Limnetic","Limnetic","Benthic","Benthic","Benthic","Limnetic","Limnetic","Limnetic","Benthic","Benthic","Benthic",
                  "Limnetic","Limnetic","Limnetic","Benthic","Benthic","Benthic","Limnetic","Limnetic","Limnetic","Benthic","Benthic","Benthic",
                  "Limnetic","Limnetic","Limnetic","Benthic","Benthic","Benthic","Limnetic","Limnetic","Limnetic","Benthic","Benthic","Benthic",
                  "Limnetic","Limnetic","Limnetic","Benthic","Benthic","Benthic","Limnetic","Limnetic","Limnetic","Benthic","Benthic","Benthic")
df_out$stage <- c("200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts",
                  "150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts",
                  "100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts","100ts",
                  "50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts")
df_out$Libraries <- c("Chip1","Chip1","Chip1","Chip1","Chip1","Chip1","Chip2","Chip2","Chip2","Chip2","Chip2","Chip2",
                      "Chip3","Chip3","Chip3","Chip3","Chip3","Chip3","Chip4","Chip4","Chip4","Chip4","Chip4","Chip4",
                      "Chip5","Chip5","Chip5","Chip5","Chip5","Chip5","Chip6","Chip6","Chip6","Chip6","Chip6","Chip6",
                      "Chip7","Chip7","Chip7","Chip7","Chip7","Chip7","Chip8","Chip8","Chip8","Chip8","Chip8","Chip8")
df_out$Sex <- c("F","M","M","M","F","M","M","M","M","M","F","F","F","NA","F","M","M","F",
                "M","F","M","F","M","F","F","F","M","F","M","M","M","F","F","M","NA","F",
                "F","F","NA","NA","NA","NA","NA","NA","F","NA","NA","NA")
df_out$stagenum <- c(200,200,200,200,200,200,200,200,200,200,200,200,
                  150,150,150,150,150,150,150,150,150,150,150,150,
                  100,100,100,100,100,100,100,100,100,100,100,100,
                  50,50,50,50,50,50,50,50,50,50,50,50)


#Make grahs of the normal distribution of the data on PCs
PC1 <- ggplot(df_out, aes(sample = PC1)) + stat_qq() + stat_qq_line() + theme_bw() + 
   labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC1")
PC2 <- ggplot(df_out, aes(sample = PC2)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC2")
PC3 <- ggplot(df_out, aes(sample = PC3)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC3")
PC4 <- ggplot(df_out, aes(sample = PC4)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC4")
PC5 <- ggplot(df_out, aes(sample = PC5)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC5")
PC6 <- ggplot(df_out, aes(sample = PC6)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC6")
PC7 <- ggplot(df_out, aes(sample = PC7)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC7")
PC8 <- ggplot(df_out, aes(sample = PC8)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC8")
PC9 <- ggplot(df_out, aes(sample = PC9)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC9")
PC10 <- ggplot(df_out, aes(sample = PC10)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC10")
PC11 <- ggplot(df_out, aes(sample = PC11)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC11")
PC12 <- ggplot(df_out, aes(sample = PC12)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC12")
PC13 <- ggplot(df_out, aes(sample = PC13)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC13")
PC14 <- ggplot(df_out, aes(sample = PC14)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC14")
PC15 <- ggplot(df_out, aes(sample = PC15)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC15")
PC16 <- ggplot(df_out, aes(sample = PC16)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC16")
PC17 <- ggplot(df_out, aes(sample = PC17)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC17")
PC18 <- ggplot(df_out, aes(sample = PC18)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC18")
PC19 <- ggplot(df_out, aes(sample = PC19)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC19")
PC20 <- ggplot(df_out, aes(sample = PC20)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC20")
PC21 <- ggplot(df_out, aes(sample = PC21)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC21")
PC22 <- ggplot(df_out, aes(sample = PC22)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC22")
PC23 <- ggplot(df_out, aes(sample = PC23)) + stat_qq() + stat_qq_line() + theme_bw() + 
  labs(x="Theoretical Quantiles", y= "Sample Quantiles") + ggtitle("PC23")

figure <- ggarrange(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,
                    ncol=5, nrow=5)


pdf("/Users/sebmatlosz/Desktop/AllFilesForSubmission/qqNorm_PCA_210222.pdf", 16,15)
figure
dev.off()




                                      #######  For PCA on 27 samples
#library(methylKit)
library(ggplot2)

#meth.min1 <- readRDS(file = "/Users/sebmatlosz/Desktop/ArnarMeeting050321/methmin1_27_noPIno100.rds")
#meth.min1df <- getData(meth.min1)

#Read in csv file directly. Cannot install methylkit on this version of R.
#Decide whether to use the full 14 427 residue or only the 9719 ones that are on placed scaffold (or 4708 on NW_)
meth.min1df <- read.table("/Users/sebmatlosz/Desktop/methylation/methmin1_27_noPIno100.csv", sep=",", header=TRUE)
#OR
meth.min1df <- read.table("/Users/sebmatlosz/Desktop/methylation/methmin1_27_noPIno100_noNW.csv", sep=",", header=TRUE)
#OR
meth.min1df <- read.table("/Users/sebmatlosz/Desktop/methylation/methmin1_27_noPIno100_noNC.csv", header=TRUE, sep=",")

#remove the first columns that I don't want
meth.min1df <- meth.min1df[,-c(1,2,3,4)]
#Remove the number of Ts columns (every 3 rows)
meth.min1df <- meth.min1df[,-seq(0,81,3)]

#Separate the data in two
numCdf <- meth.min1df[,seq(0,54,2)]
coveragedf <- meth.min1df[,seq(1,54,2)]

total <- ((numCdf*100)/coveragedf)

#Should invert the df to plot based on samples. Otherwise I plot the entire 10 000 residues.
totalt <- t(total)

df_PCA <- prcomp(totalt)
df_out <- as.data.frame(df_PCA$x)
#This was mainly for curiosity
df_out$morph <- c("PL","PL","PL","LB","LB","LB","SB","SB","SB",
                  "PL","PL","LB","LB","SB","PL","SB","SB","LB",
                  "PL","PL","LB","LB","SB","PL","SB","SB","LB")

df_out$niche <- c("Limnetic","Limnetic","Limnetic","Benthic","Benthic","Benthic","Benthic","Benthic","Benthic",
                  "Limnetic","Limnetic","Benthic","Benthic","Benthic","Limnetic","Benthic","Benthic","Benthic",
                  "Limnetic","Limnetic","Benthic","Benthic","Benthic","Limnetic","Benthic","Benthic","Benthic")

df_out$stage <- c("200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts","200ts",
                  "150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts","150ts",
                  "50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts","50ts")

df_out$Sex <- c("F","M","M","M","F","M","M","F","F",
                "F","M","M","M","F","M","F","M","F",
                "NA","NA","NA","NA","NA","NA","NA","NA","NA")

df_out$stagenum <- c(200,200,200,200,200,200,200,200,200,
                     150,150,150,150,150,150,150,150,150,
                     50,50,50,50,50,50,50,50,50)

############################## Make ANOVA test on the first 8 PCAs. Meaning on each of the first 8. #########################
library(dplyr)
library(reshape2)
library(DataCombine)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(EnvStats)
library(tidyverse)
library(agricolae)
library(forecast)
library(AICcmodavg)

#Make a recap for each PC
PCrecap <- function(i) {
  x <- data.frame(cbind(df_out[,i],df_out$morph,df_out$niche,df_out$stage,df_out$Sex,df_out$stagenum))
  colnames(x) <- c("perc","morph","niche","stage","sex","stagenum")
  x$perc <- as.numeric(as.character(x$perc))
  return(x)
}

PC1recap <- PCrecap(1)
PC2recap <- PCrecap(2)
PC3recap <- PCrecap(3)
PC4recap <- PCrecap(4)
PC5recap <- PCrecap(5)
PC6recap <- PCrecap(6)
PC7recap <- PCrecap(7)
PC8recap <- PCrecap(8)
PC9recap <- PCrecap(9)
PC10recap <- PCrecap(10)
PC11recap <- PCrecap(11)
PC12recap <- PCrecap(12)
PC13recap <- PCrecap(13)
PC14recap <- PCrecap(14)
PC15recap <- PCrecap(15)
PC16recap <- PCrecap(16)
PC17recap <- PCrecap(17)
PC18recap <- PCrecap(18)
PC19recap <- PCrecap(19)
PC20recap <- PCrecap(20)
PC21recap <- PCrecap(21)
PC22recap <- PCrecap(22)
PC23recap <- PCrecap(23)
PC24recap <- PCrecap(24)
PC25recap <- PCrecap(25)
PC26recap <- PCrecap(26)
PC27recap <- PCrecap(27)
PC28recap <- PCrecap(28)
PC29recap <- PCrecap(29)
PC30recap <- PCrecap(30)
PC31recap <- PCrecap(31)
PC32recap <- PCrecap(32)
PC33recap <- PCrecap(33)
PC34recap <- PCrecap(34)
PC35recap <- PCrecap(35)
PC36recap <- PCrecap(36)
PC37recap <- PCrecap(37)
PC38recap <- PCrecap(38)
PC39recap <- PCrecap(39)
PC40recap <- PCrecap(40)
PC41recap <- PCrecap(41)
PC42recap <- PCrecap(42)
PC43recap <- PCrecap(43)
PC44recap <- PCrecap(44)
PC45recap <- PCrecap(45)
PC46recap <- PCrecap(46)
PC47recap <- PCrecap(47)











#Check if we want to transform the data
              #qqnorm(PC1recap$perc)
              #qqline(PC1recap$perc)
              
              #qqnorm(PC2recap$perc)
              #qqline(PC2recap$perc)
              
              #qqnorm(PC3recap$perc)
              #qqline(PC3recap$perc)
              #I guess transforming is better. But somehow it gives error message ? Still seem to work though.


#test just to try things out
# datatable <- PC1recap
# lambdaValue <- BoxCox.lambda(datatable$perc, method = "guerrero")
# datatable$transformed <- BoxCox(datatable$perc, lambdaValue)

#interaction <- lm(transformed ~ morph * stagenum * sex, data=datatable)
#i <- anova(interaction)
#i


#Now make a function that looks at every one of those PC and performs linear regression.
linearfunctions <- function(datatable,name) {
  
  #transform the data
  #lambdaValue <- BoxCox.lambda(datatable$perc, method = "guerrero")
  #datatable$transformed <- BoxCox(datatable$perc, lambdaValue)
  datatable$transformed <- datatable$perc
  
  #create the two linear models
  #fixed effect
  fixedeffect <- lm(transformed ~ morph + stagenum + sex, data=datatable)
  #interaction
  interaction <- lm(transformed ~ morph * stagenum * sex, data=datatable)
  
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
PC5table <- linearfunctions(PC5recap,"PC5")
PC6table <- linearfunctions(PC6recap,"PC6")
PC7table <- linearfunctions(PC7recap,"PC7")
PC8table <- linearfunctions(PC8recap,"PC8")
PC9table <- linearfunctions(PC9recap,"PC9")
PC10table <- linearfunctions(PC10recap,"PC10")
PC11table <- linearfunctions(PC11recap,"PC11")
PC12table <- linearfunctions(PC12recap,"PC12")
PC13table <- linearfunctions(PC13recap,"PC13")
PC14table <- linearfunctions(PC14recap,"PC14")
PC15table <- linearfunctions(PC15recap,"PC15")
PC16table <- linearfunctions(PC16recap,"PC16")
PC17table <- linearfunctions(PC17recap,"PC17")
PC18table <- linearfunctions(PC18recap,"PC18")
PC19table <- linearfunctions(PC19recap,"PC19")
PC20table <- linearfunctions(PC20recap,"PC20")
PC21table <- linearfunctions(PC21recap,"PC21")
PC22table <- linearfunctions(PC22recap,"PC22")
PC23table <- linearfunctions(PC23recap,"PC23")
PC24table <- linearfunctions(PC24recap,"PC24")
PC25table <- linearfunctions(PC25recap,"PC25")
PC26table <- linearfunctions(PC26recap,"PC26")
PC27table <- linearfunctions(PC27recap,"PC27")
PC28table <- linearfunctions(PC28recap,"PC28")
PC29table <- linearfunctions(PC29recap,"PC29")
PC30table <- linearfunctions(PC30recap,"PC30")
PC31table <- linearfunctions(PC31recap,"PC31")
PC32table <- linearfunctions(PC32recap,"PC32")
PC33table <- linearfunctions(PC33recap,"PC33")
PC34table <- linearfunctions(PC34recap,"PC34")
PC35table <- linearfunctions(PC35recap,"PC35")
PC36table <- linearfunctions(PC36recap,"PC36")
PC37table <- linearfunctions(PC37recap,"PC37")
PC38table <- linearfunctions(PC38recap,"PC38")
PC39table <- linearfunctions(PC39recap,"PC39")
PC40table <- linearfunctions(PC40recap,"PC40")
PC41table <- linearfunctions(PC41recap,"PC41")
PC42table <- linearfunctions(PC42recap,"PC42")
PC43table <- linearfunctions(PC43recap,"PC43")
PC44table <- linearfunctions(PC44recap,"PC44")
PC45table <- linearfunctions(PC45recap,"PC45")
PC46table <- linearfunctions(PC46recap,"PC46")
PC47table <- linearfunctions(PC47recap,"PC47")

recap <- rbind(PC1table,PC2table,PC3table,PC4table,PC5table,PC6table,PC7table,PC8table,PC9table,PC10table,PC11table,PC12table,
               PC13table,PC14table,PC15table,PC16table,PC17table,PC18table,PC19table,PC20table,PC21table,PC22table,PC23table,
               PC24table,PC25table,PC26table,PC27table,PC28table,PC29table,PC30table,PC31table,PC32table,PC33table,PC34table,
               PC35table,PC36table,PC37table,PC38table,PC39table,PC40table,PC41table,PC42table,PC43table,PC44table,PC45table,
               PC46table,PC47table)


#Correct for multiple testing with Bonferroni
recap$Morph_pvalue_corrected <- p.adjust(recap$Morph.pvalue., method = "bonferroni", n = length(recap$Morph.pvalue.))
recap$Time_pvalue_corrected <- p.adjust(recap$Time.pvalue., method = "bonferroni", n = length(recap$Time.pvalue.))
recap$Sex_pvalue_corrected <- p.adjust(recap$Sex.pvalue., method = "bonferroni", n = length(recap$Sex.pvalue.))
recap$MorphxTime_pvalue_corrected <- p.adjust(recap$Morph.Time, method = "bonferroni", n = length(recap$Morph.Time))



write.table(recap,"/Users/sebmatlosz/Desktop/AllFilesForSubmission/linearmodelsPCA_48samples_notrans_bonfcorrection_230322.csv", row.names = FALSE, col.names = TRUE, sep= ",")

                      #Just some plotting
                      #Need to find a way to see whether niches are involved more than morphs, or which morph is separated. Maybe just plotting it would work.
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
                      #colors <- c("LB" = "#00BA38FF", "SB" = "#619CFFFF", "PL" = "#F8766DFF", "PI" = "darkmagenta" )
                      colors <- c("LB" = "#00BA38FF", "SB" = "#619CFFFF", "PL" = "#F8766DFF")
                      #colors <- c("F" = "red", "M"="blue", "NA"="grey")
                      p<-ggplot(df_out,aes(x=PC1,y=PC4,color=morph, shape=stage)) + geom_point(size=3) + 
                        labs(x=percentage[1],y=percentage[4],title = "") + theme_bw() + scale_color_manual(values = colors)
                      p
                      
                      pdf("/Users/sebmatlosz/Desktop/Figuresforpapers/27samplesnoNC_PC1og4_241121.pdf")
                      print(p)
                      dev.off() 

library(DT)
datatable(recap, rownames = FALSE) %>%
    formatStyle(columns = c(1:7), 
              color = styleInterval(c(0,0.05,1), c("black", "red","black","black")))

#Extralong recap
PC10recap <- PCrecap(10)
PC11recap <- PCrecap(11)
PC12recap <- PCrecap(12)
PC13recap <- PCrecap(13)
PC14recap <- PCrecap(14)
PC15recap <- PCrecap(15)
PC16recap <- PCrecap(16)
PC17recap <- PCrecap(17)
PC18recap <- PCrecap(18)
PC19recap <- PCrecap(19)
PC20recap <- PCrecap(20)
PC21recap <- PCrecap(21)
PC22recap <- PCrecap(22)
PC23recap <- PCrecap(23)
PC24recap <- PCrecap(24)
PC25recap <- PCrecap(25)
PC26recap <- PCrecap(26)
PC27recap <- PCrecap(27)
PC28recap <- PCrecap(28)
PC29recap <- PCrecap(29)
PC30recap <- PCrecap(30)
PC31recap <- PCrecap(31)
PC32recap <- PCrecap(32)
PC33recap <- PCrecap(33)
PC34recap <- PCrecap(34)
PC35recap <- PCrecap(35)
PC36recap <- PCrecap(36)
PC37recap <- PCrecap(37)
PC38recap <- PCrecap(38)
PC39recap <- PCrecap(39)
PC40recap <- PCrecap(40)
PC41recap <- PCrecap(41)
PC42recap <- PCrecap(42)
PC43recap <- PCrecap(43)
PC44recap <- PCrecap(44)
PC45recap <- PCrecap(45)
PC46recap <- PCrecap(46)
PC47recap <- PCrecap(47)

PC10table <- linearfunctions(PC10recap,"PC10")
PC11table <- linearfunctions(PC11recap,"PC11")
PC12table <- linearfunctions(PC12recap,"PC12")
PC13table <- linearfunctions(PC13recap,"PC13")
PC14table <- linearfunctions(PC14recap,"PC14")
PC15table <- linearfunctions(PC15recap,"PC15")
PC16table <- linearfunctions(PC16recap,"PC16")
PC17table <- linearfunctions(PC17recap,"PC17")
PC18table <- linearfunctions(PC18recap,"PC18")
PC19table <- linearfunctions(PC19recap,"PC19")
PC20table <- linearfunctions(PC20recap,"PC20")
PC21table <- linearfunctions(PC21recap,"PC21")
PC22table <- linearfunctions(PC22recap,"PC22")
PC24table <- linearfunctions(PC24recap,"PC24")
PC25table <- linearfunctions(PC25recap,"PC25")
PC26table <- linearfunctions(PC26recap,"PC26")
PC27table <- linearfunctions(PC27recap,"PC27")
PC28table <- linearfunctions(PC28recap,"PC28")
PC29table <- linearfunctions(PC29recap,"PC29")
PC30table <- linearfunctions(PC30recap,"PC30")
PC31table <- linearfunctions(PC31recap,"PC31")
PC32table <- linearfunctions(PC32recap,"PC32")
PC33table <- linearfunctions(PC33recap,"PC33")
PC34table <- linearfunctions(PC34recap,"PC34")
PC35table <- linearfunctions(PC35recap,"PC35")
PC36table <- linearfunctions(PC36recap,"PC36")
PC37table <- linearfunctions(PC37recap,"PC37")
PC38table <- linearfunctions(PC38recap,"PC38")
PC39table <- linearfunctions(PC39recap,"PC39")
PC40table <- linearfunctions(PC40recap,"PC40")
PC41table <- linearfunctions(PC41recap,"PC41")
PC42table <- linearfunctions(PC42recap,"PC42")
PC43table <- linearfunctions(PC43recap,"PC43")
PC44table <- linearfunctions(PC44recap,"PC44")
PC45table <- linearfunctions(PC45recap,"PC45")
PC46table <- linearfunctions(PC46recap,"PC46")
PC47table <- linearfunctions(PC47recap,"PC47")



BIGrecap <- rbind(PC1table,PC2table,PC3table,PC4table,PC5table,PC6table,PC7table,PC8table,PC9table,PC10table,
             PC11table,PC12table,PC13table,PC14table,PC15table,PC16table,PC17table,PC18table,PC19table,PC20table,
             PC21table,PC22table,PC23table,PC24table,PC25table,PC26table,PC27table,PC28table,PC29table,PC30table,
             PC31table,PC32table,PC33table,PC34table,PC35table,PC36table,PC37table,PC38table,PC39table,PC40table,
             PC41table,PC42table,PC43table,PC44table,PC45table,PC46table,PC47table)

colnames(BIGrecap) <- c("PC","Morph","Time","Sex","Morph X Time","Morph X Sex", "Time X Sex", "Morph X Time X Sex")

library(DT)
datatable(BIGrecap[,1:5], rownames = FALSE) %>%
  formatStyle(columns = c(1:5), 
              color = styleInterval(c(0,0.05,1), c("black", "red","black","black")))

#Save table as .csv
write.csv(BIGrecap[,1:5],"/Users/sebmatlosz/Desktop/Figuresforpapers/FinalFiguresForSubmissions/LinearModelsPCA.csv", row.names = FALSE)
write.csv(BIGrecap,"/Users/sebmatlosz/Desktop/Figuresforpapers/FinalFiguresForSubmissions/LinearModelsPCA_allint.csv", row.names = FALSE)




## We could do the same thing using aov() instead of lm(). DOESNT CHANGE ANYTHING

linearfunctionsaov <- function(datatable,name) {
  
  #transform the data
  lambdaValue <- BoxCox.lambda(datatable$perc, method = "guerrero")
  datatable$transformed <- BoxCox(datatable$perc, lambdaValue)
  
  #create the two linear models
  #fixed effect
  fixedeffect <- aov(transformed ~ morph + stage + sex, data=datatable)
  #interaction
  interaction <-aov(transformed ~ morph * stage * sex, data=datatable)
  
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
      
      #create dataframe with all the values we want
      test <- data.frame("PC" = name, "Morph(pvalue)" = a$`Pr(>F)`[1], "Time(pvalue)" = a$`Pr(>F)`[2], "Sex(pvalue)" = a$`Pr(>F)`[3] )
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
    
    
    #create dataframe with all the values we want
    test <- data.frame("PC" = name, "Morph(pvalue)" = a$`Pr(>F)`[1], "Time(pvalue)" = a$`Pr(>F)`[2], "Sex(pvalue)" = a$`Pr(>F)`[3] )
  }
  return(test)
}

#Run the function on the data for all the genes
PC1table <- linearfunctionsaov(PC1recap,"PC1")
PC2table <- linearfunctionsaov(PC2recap,"PC2")
PC3table <- linearfunctionsaov(PC3recap,"PC3")
PC4table <- linearfunctionsaov(PC4recap,"PC4")
PC5table <- linearfunctionsaov(PC5recap,"PC5")
PC6table <- linearfunctionsaov(PC6recap,"PC6")
PC7table <- linearfunctionsaov(PC7recap,"PC7")
PC8table <- linearfunctionsaov(PC8recap,"PC8")
PC23table <- linearfunctionsaov(PC23recap,"PC23")

recap <- rbind(PC1table,PC2table,PC3table,PC4table,PC5table,PC6table,PC7table,PC8table,PC23table)

datatable(recap, rownames = FALSE) %>%
  formatStyle(columns = c("Morph.pvalue.","Time.pvalue.","Sex.pvalue."), 
              color = styleInterval(c(0,0.05,1), c("black", "red","black","black")))

#Doesn't change anything




#OLD STUFF I KEEP IN CASE I NEED THE CODE
#Change the chr into NC_ and NW_
meth.min1df$temp <- gsub("chr","",meth.min1df$chr)
meth.min1df$NW <- ifelse(nchar(meth.min1df$temp)==11, "NW_", "NC_")

meth.mindfnoNWx <- dplyr::filter(meth.min1df, !grepl("NW_",NW))
meth.mindfnoNW <- meth.mindfnoNWx[,1:85]
#write.table(meth.mindfnoNW, "/Users/sebmatlosz/Desktop/methylation/methmin1_27_noPIno100_noNW.csv", row.names = F, quote = F, sep=",")