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
    ddCt.Lmtk2 <- readRDS("~/REr_Lmtk2.rds")
    ddCt.HiH2A <- readRDS("~/REr_HiH2A.rds")
    ddCt.SLC9A3R2 <- readRDS("~/REr_SLC9A3R2.rds")
    ddCt.Nkx232 <- readRDS("~/REr_Nkx23.rds")
    ddCt.NFiX1 <- readRDS("~/REr_NFIX.rds")
    ddCt.RASSF42 <- readRDS( "~/REr_RASSF4.rds")
    ddCt.ARF161 <- readRDS("~/REr_ARL16.rds")
    ddCt.MEiS1l2 <- readRDS("~/REr_MEIS1.rds")
    ddCt.HiH3l2 <- readRDS("~/REr_H3-like.rds")
    ddCt.ARMC11 <- readRDS("~/REr_ARMC1.rds")
    ddCt.GLi31 <- readRDS("~/REr_Gli3.rds")
    ddCt.MEGF91 <- readRDS("~/REr_MEGF9.rds")
    ddCt.MAGUK2 <- readRDS("~/REr_MAGUK2.rds")
    ddCt.Rhoguanin1 <- readRDS("~/REr_Rhoguanin.rds")

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
  
  colnames(datatable)[15] <- "REr"   
  colnames(datatable)[6] <- "morphs" 
  colnames(datatable)[7] <- "times"  
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
    Lmtk2Recap <- linearfunctions(ddCt.Lmtk2,"Lmtk2")
    HiH2ARecap <- linearfunctions(ddCt.HiH2A,"HiH2A")
    SLC9A3R2Recap <- linearfunctions(ddCt.SLC9A3R2,"SLC9A3R2")
    ARF161Recap <- linearfunctions(ddCt.ARF161,"ARF161")
    MAGUK2Recap <- linearfunctions(ddCt.MAGUK2,"MAGUK2")
    NFiX1Recap <- linearfunctions(ddCt.NFiX1,"NFiX1")


Recap <- rbind(Lmtk2Recap,HiH2ARecap,SLC9A3R2Recap,Nkx232Recap,NFiX1Recap,RASSF42Recap,ARF161Recap,
               MEiS1l2Recap,HiH3l2Recap,ARMC11Recap,Gli31Recap,MEGF91Recap,MAGUK2Recap,Rhoguanin1Recap)

colnames(Recap) <- c("Name","Morph p.value", "LB","SB","PI","PL","Time p.value","100ts","150ts","200ts","Interaction p.value")

# write.csv(Recap,"/Users/sebmatlosz/Desktop/AllFilesForSubmission/qPCRAnalysisRecap_REr_withtransformation.csv", row.names = FALSE)





