# Make graphs of pvalues by glm ~ weight on PC
library(ggplot2)
library(ggpubr)
library(sqldf)

#Load glm data
    glmfinal <- read.csv("~/GLM/results/glm_10340res_pvalues_withandwithout_Bonfcorrection.csv")

#Load PCA data
    PCAfinal <- read.csv("~/PCA/results/weightsonPC_48samplesPCA.csv")
  
#Merge data frames based on residues
#Only keep the few PCs of interest
    df <- sqldf('SELECT Residue, Morph_pvalue, Morph_pvalue_corrected, Timepoint_pvalue, Timepoint_pvalue_corrected, Sex_pvalue, 
                  Sex_pvalue_corrected, PC1, PC4, PC23,ID
                 FROM glmfinal, PCAfinal
                 WHERE Residue=ID')

#Change the PC weights into absolute values
    df$PC1 <- abs(df$PC1)
    df$PC4 <- abs(df$PC4)
    df$PC23 <- abs(df$PC23)

#Change the pvalues into -log
    df$Morph_pvalue_corrected_log <- -log(df$Morph_pvalue_corrected)
    df$Timepoint_pvalue_corrected_log <- -log(df$Timepoint_pvalue_corrected)
    df$Sex_pvalue_corrected_log <- -log(df$Sex_pvalue_corrected)

    df$Morph_pvalue_log <- -log(df$Morph_pvalue)
    df$Timepoint_pvalue_log <- -log(df$Timepoint_pvalue)
    df$Sex_pvalue_log <- -log(df$Sex_pvalue)


###Check the correlations

correlation <- function(x,y,t) {
  title=paste0(t)
  z <- ggscatter(df, x = x, y = y, 
                 add = "reg.line", conf.int = TRUE, 
                 cor.coef = FALSE, cor.method = "kendall",
                 xlab = "", ylab = "") + ggtitle(title) + theme_bw() + theme(plot.title = element_text(hjust=0.5))
  return(z)
}

    Morphs <- correlation(x = 'PC4',y = 'Morph_pvalue_corrected_log', "Correlation PCA vs glm / Morphs")
    Morphs
    Morphcorrnumber <- cor.test(df$PC4, df$Morph_pvalue_corrected_log, method = "kendall")
    Morphcorrnumber$p.value
    Morphcorrnumber$estimate

    Times <- correlation(x = 'PC1',y = 'Timepoint_pvalue_corrected_log', "Correlation PCA vs glm / Timepoints")
    Times
    Timecorrnumber <- cor.test(df$PC1, df$Timepoint_pvalue_corrected_log, method = "kendall")
    Timecorrnumber$p.value
    Timecorrnumber$estimate


#Plot things
    timepointgraph <- ggplot(df, aes(x=PC1, y=Timepoint_pvalue_corrected_log)) + geom_point() + theme_bw() + 
      ggtitle("Residues affecting timepoints - Model comparison") + labs(x="Weight on PC1",y="Timepoint -log(pvalue) in glm")
    timepointgraph

    morphgraph <- ggplot(df, aes(x=PC4, y=Morph_pvalue_corrected_log)) + geom_point() + theme_bw() + 
      ggtitle("Residues affecting Niches - Model comparison") + labs(x="Weight on PC4",y="Niche -log(pvalue) in glm")
    morphgraph

    sexgraph <- ggplot(df, aes(x=PC23, y=Sex_pvalue_corrected_log)) + geom_point() + theme_bw() + 
      ggtitle("Residues affecting sex - Model comparison") + labs(x="Weight on PC23",y="Sex -log(pvalue) in glm")
    sexgraph

    figure <- ggarrange(timepointgraph,morphgraph,sexgraph, ncol=2, nrow=2)
    figure

pdf("~/results/ModelComparisons_minuslog_corrected_glmPCA_090222.pdf")
figure
dev.off()
