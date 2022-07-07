# Final script to interpret qPCR data
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

#Load the data made by taking into account the primer efficiencies:
#Only use the REr columns.
ddCt.Lmtk2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_Lmtk2.rds")
ddCt.HiH2A <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_HiH2A.rds")
ddCt.SLC9A3R2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_SLC9A3R2.rds")
ddCt.Nkx232 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_Nkx23.rds")
ddCt.NFiX1 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_NFIX.rds")
ddCt.RASSF42 <- readRDS( "/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_RASSF4.rds")
ddCt.ARF161 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_ARL16.rds")
ddCt.MEiS1l2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_MEIS1.rds")
ddCt.HiH3l2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_H3-like.rds")
ddCt.ARMC11 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_ARMC1.rds")
ddCt.GLi31 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_Gli3.rds")
ddCt.MEGF91 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_MEGF9.rds")
ddCt.MAGUK2 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_MAGUK2.rds")
ddCt.Rhoguanin1 <- readRDS("/Users/sebmatlosz/Desktop/qPCR_dCtobjects/primerEffNorm/REr_Rhoguanin.rds")

##Check normality with KS to see if we need to transform the data
qqnorm(ddCt.NFiX1$dCt)
qqline(ddCt.NFiX1$dCt)

##ATH ARL16, H2A, LMTK2, MPP3, NFIX1, SLC9A3R2 need to be transformed
lambdaValue <- BoxCox.lambda(ddCt.ARF161$`RE_ARF16-1`, method = "guerrero")
ddCt.ARF161$`RE_ARF16-1` <- BoxCox(ddCt.ARF161$`RE_ARF16-1`, lambdaValue)

lambdaValue <- BoxCox.lambda(ddCt.HiH2A$RE_HiH2A, method = "guerrero")
ddCt.HiH2A$RE_HiH2A <- BoxCox(ddCt.HiH2A$RE_HiH2A, lambdaValue)

lambdaValue <- BoxCox.lambda(ddCt.Lmtk2$RE_Lmtk2, method = "guerrero")
ddCt.Lmtk2$RE_Lmtk2 <- BoxCox(ddCt.Lmtk2$RE_Lmtk2, lambdaValue)

lambdaValue <- BoxCox.lambda(ddCt.MAGUK2$`RE_MAGUK-2`, method = "guerrero")
ddCt.MAGUK2$`RE_MAGUK-2` <- BoxCox(ddCt.MAGUK2$`RE_MAGUK-2`, lambdaValue)

lambdaValue <- BoxCox.lambda(ddCt.NFiX1$`RE_NFiX-1`, method = "guerrero")
ddCt.NFiX1$`RE_NFiX-1` <- BoxCox(ddCt.NFiX1$`RE_NFiX-1`, lambdaValue)

lambdaValue <- BoxCox.lambda(ddCt.SLC9A3R2$RE_SLC9A3R2, method = "guerrero")
ddCt.SLC9A3R2$RE_SLC9A3R2 <- BoxCox(ddCt.SLC9A3R2$RE_SLC9A3R2, lambdaValue)


#### ATH, change datatable$dCt into datatable$ddCt or datatable$REr based on what you are looking at
#### ATH Check whether you want to transform the data.
linearfunctions <- function(datatable,name) {
  
  #transform the data
  #lambdaValue <- BoxCox.lambda(datatable$ddCt, method = "guerrero")
  #datatable$transformed <- BoxCox(datatable$ddCt, lambdaValue)
  
  colnames(datatable)[15] <- "REr"   #Only run this line when looking at REr
  colnames(datatable)[6] <- "morphs" #Only run this line when looking at REr
  colnames(datatable)[7] <- "times"  #Only run this line when looking at REr
  datatable$transformed <- datatable$REr
  
  #create the two linear models
  #fixed effect
  fixedeffect <- lm(transformed ~ morphs + times, data=datatable)
  #interaction
  interaction <- lm(transformed ~ morphs * times, data=datatable)
  
  # check if there is difference between the fixed effect and the interaction
  anovatest <- anova(fixedeffect,interaction)
  
  # if the difference is significant
  if(anovatest$`Pr(>F)`[2] < 0.05) {
    # do an AICc test on both the linear regressions
    aicc.fixed <- AICc(fixedeffect,return.K = FALSE)
    aicc.inter <- AICc(interaction,return.K = FALSE) 
    
    #if fixed effect is better, use it
    if(aicc.fixed < aicc.inter) {
      a <- anova(fixedeffect) #we will use the fixed effect model to get the morphs and times pvalue
      print("fixed effect")
      morph <- HSD.test(fixedeffect, "morphs", console=TRUE)
      morph1 <- morph$groups #makes it easier to select the proper rows. Morphs might not be ordered thesame otherwise
      time <- HSD.test(fixedeffect, "times", console=TRUE)
      time1 <- time$groups
      
      b <- anova(interaction) #we will use the interaction model to get the interaction pvalue
      
      #create dataframe with all the values we want
      test <- data.frame("GeneName" = name, "Morph(pvalue)" = a$`Pr(>F)`[1], "LB" = morph1["LB","groups"], 
                         "SB" = morph1["SB","groups"],"PI" = morph1["PI","groups"],"PL" = morph1["PL","groups"],
                         "Time(pvalue)" = a$`Pr(>F)`[2], "100" = time1["t100","groups"],"150" = time1["t150","groups"],
                         "200" = time1["t200","groups"],"Interaction(pvalue)" = b$`Pr(>F)`[3] )
    }
    
    
    #if interaction model is better, use it
    else if(aicc.fixed > aicc.inter) {
      c <- anova(interaction) #we will use the interaction model to get the morphs, times AND interaction pvalue
      print("interaction")
      morph <- HSD.test(interaction, "morphs", console=TRUE)
      morph1 <- morph$groups
      time <- HSD.test(interaction, "times", console=TRUE)
      time1 <- time$groups
      
      #create dataframe with all the values we want
      test <- data.frame("GeneName" = name, "Morph(pvalue)" = c$`Pr(>F)`[1], "LB" = morph1["LB","groups"], 
                         "SB" = morph1["SB","groups"],"PI" = morph1["PI","groups"],"PL" = morph1["PL","groups"],
                         "Time(pvalue)" = c$`Pr(>F)`[2], "100" = time1["t100","groups"],"150" = time1["t150","groups"],
                         "200" = time1["t200","groups"],"Interaction(pvalue)" = c$`Pr(>F)`[3] )
    }
  }
  
  #If no difference between the models, we choose fixedeffect (basically same thing than on line #65)
  else if(anovatest$`Pr(>F)`[2] > 0.05) {
    a <- anova(fixedeffect) 
    print("fixed effect")
    morph <- HSD.test(fixedeffect, "morphs", console=TRUE)
    morph1 <- morph$groups #makes it easier to select the proper rows. Morphs might not be ordered thesame otherwise
    time <- HSD.test(fixedeffect, "times", console=TRUE)
    time1 <- time$groups
    
    b <- anova(interaction) #we will STILL use the interaction model to get the interaction pvalue
    
    #create dataframe with all the values we want
    test <- data.frame("GeneName" = name, "Morph(pvalue)" = a$`Pr(>F)`[1], "LB" = morph1["LB","groups"], 
                       "SB" = morph1["SB","groups"],"PI" = morph1["PI","groups"],"PL" = morph1["PL","groups"],
                       "Time(pvalue)" = a$`Pr(>F)`[2], "100" = time1["t100","groups"],"150" = time1["t150","groups"],
                       "200" = time1["t200","groups"],"Interaction(pvalue)" = b$`Pr(>F)`[3] )
  }
  return(test)
}


#Run the function on the data for all the genes
RASSF42Recap <- linearfunctions(ddCt.RASSF42,"RASSF42")
MEGF91Recap <- linearfunctions(ddCt.MEGF91,"MEGF91")
Nkx232Recap <- linearfunctions(ddCt.Nkx232,"Nkx232")
ARMC11Recap <- linearfunctions(ddCt.ARMC11,"ARMC11")
MEiS1l2Recap <- linearfunctions(ddCt.MEiS1l2,"MEiS1l2")
HiH3l2Recap <- linearfunctions(ddCt.HiH3l2 ,"HiH3l2")
Gli31Recap <- linearfunctions(ddCt.GLi31,"Gli31")
Rhoguanin1Recap <- linearfunctions(ddCt.Rhoguanin1,"Rhoguanin1")

#These next genes need to be transformed so change the function at the beginning
#Unless I have already transformed those who need transforming.
Lmtk2Recap <- linearfunctions(ddCt.Lmtk2,"Lmtk2")
HiH2ARecap <- linearfunctions(ddCt.HiH2A,"HiH2A")
SLC9A3R2Recap <- linearfunctions(ddCt.SLC9A3R2,"SLC9A3R2")
ARF161Recap <- linearfunctions(ddCt.ARF161,"ARF161")
MAGUK2Recap <- linearfunctions(ddCt.MAGUK2,"MAGUK2")
NFiX1Recap <- linearfunctions(ddCt.NFiX1,"NFiX1")


Recap <- rbind(Lmtk2Recap,HiH2ARecap,SLC9A3R2Recap,Nkx232Recap,NFiX1Recap,RASSF42Recap,ARF161Recap,
               MEiS1l2Recap,HiH3l2Recap,ARMC11Recap,Gli31Recap,MEGF91Recap,MAGUK2Recap,Rhoguanin1Recap)

colnames(Recap) <- c("Name","Morph p.value", "LB","SB","PI","PL","Time p.value","100ts","150ts","200ts","Interaction p.value")

write.csv(Recap,"/Users/sebmatlosz/Desktop/AllFilesForSubmission/qPCRAnalysisRecap_REr_withtransformation_210222.csv", row.names = FALSE)

#At this point, we have a recap of everything.
#However, in some cases, differences between morphs could be invisible in this data but visible in the other data, 
# with Lea's way of calculating the ddCt. So we should add this as well.
#NEVERMIND, It gives pretty much the same results, if not less interesting

#Those were made with LEA's way
#ddCt.Lmtk2 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-Lmtk2.rds")
#ddCt.HiH2A <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-HiH2A.rds")
#ddCt.SLC9A3R2 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-SLC9A3R2.rds")
#ddCt.Nkx232 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-Nkx232.rds")
#ddCt.NFiX1 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-NFiX1.rds")
#ddCt.RASSF42 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-RASSF42.rds")
#ddCt.ARF161 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-ARF161.rds")
#ddCt.MEiS1l2 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-MEiS1l2.rds")
#ddCt.HiH3l2 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-HiH3l2.rds")
#ddCt.ARMC11 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-ARMC11.rds")
#ddCt.GLi31 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-GLi31.rds")
#ddCt.MEGF91 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-MEGF91.rds")
#ddCt.MAGUK2 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-MAGUK2.rds")
#ddCt.Rhoguanin1 <- readRDS("/Users/sebmatlosz/Desktop/ddCtobjectsforLN/LEAsWAY/ddCt-Rhoguanin1.rds")


# Need to reformat the Recap table. Good column names and round numbers
round(Recap$Morph.pvalue., digits=2)
options("scipen"=-100, "digits"=3)
#Recap$Morph.pvalue. <- as.character(as.numeric(Recap$Morph.pvalue.)) #Doesn't work. Might need to change it in Excel or manually..

#Somehow it works to convert to character once I reimport it from Excel.
Recap2 <- read.csv("/Users/sebmatlosz/Desktop/qPCRanalysisRecapforDT.csv", sep = ',')
#Recap2 <- Recap2[,-1]
Recap2$Morph.p.value <- as.character(Recap2$Morph.p.value)
Recap2$Time.p.value <- as.character(Recap2$Time.p.value)
Recap2$Interaction.p.value <- as.character(Recap2$Interaction.p.value)
#Change colnames
colnames(Recap2) <- c("Name","Morph p.value", "LB","SB","PI","PL","Time p.value","100ts","150ts","200ts","Interaction p.value")


# Once it is changed to character, can change the colours of specific cells in datatable.
#Actually I can't or I forgot.. So I am going to do this in a scuffed way by specifying every string that needs to be changed..
library(DT)
datatable(Recap2, rownames = FALSE) %>%
  formatStyle(columns = c("LB","SB","PI","PL","100ts","150ts","200ts"), 
              backgroundColor = styleEqual(c("a", "ab","b","bc","c"), c('#ca0020', '#f4a582', '#f7f7f7', '#92c5de', '#0571b0'))) %>%
  formatStyle(columns = c("Time p.value"), 
              color = styleEqual(c("1.51e-02", "2e-06","3.27e-10","3.87e-05","1.63e-30","6.56e-04","3.16e-07","1.83e-04","8.24e-10","3.68e-02","2.85e-07","5.39e-14","7.26e-08"), c("red","red","red","red","red","red","red","red","red","red","red","red","red"))) %>%
  formatStyle(columns = c("Interaction p.value"), 
              color = styleEqual(c("8.32e-03", "8.13e-03","3.53e-02","9.88e-03"), c("red","red","red","red"))) %>%
  formatStyle(columns = c("Morph p.value"), 
              color = styleEqual(c("1.99e-05", "9.27e-03","4.12e-07","2.27e-04","2.59e-05","3.88e-02","2.29e-06","1.48e-05"), c("red","red","red","red","red","red","red","red")))




##################################################################################################################

#I am keeping the following code for use if I ever have to color code using numbers
#Use the diverging color scale from the colour palette paper
library(DT)
datatable(Recap, rownames = FALSE) %>%
  #formatStyle(columns = c("Morph.pvalue.","Time.pvalue.","Interaction.pvalue."), 
              # background = styleInterval(c(0,0.05,1), c("white", "salmon","white","white"))) %>%
  formatStyle(columns = c("LB","SB","PI","PL","X100","X150","X200"), 
               backgroundColor = styleEqual(c("a", "ab","b","bc","c"), c('#ca0020', '#f4a582', '#f7f7f7', '#92c5de', '#0571b0'))) %>%
  formatStyle(columns = c("Morph.pvalue.","Time.pvalue.","Interaction.pvalue."), 
              color = styleInterval(c(0,0.05,1), c("black", "red","black","black"))) %>%
  formatSignif(columns = c("Morph.pvalue.","Time.pvalue.","Interaction.pvalue."), digits=3)






 





