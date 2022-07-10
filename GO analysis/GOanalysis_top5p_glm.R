#### Look at the residues with the lowest pvalue in glm, that explain a variable  ####
library(stringr)
library(dplyr)

#Load the data from the glm (the 517 residues / 10340)
#top5pMorph <- read.csv("/Users/sebmatlosz/Desktop/filesforGOanalysis/glm_top5p_Morphs.csv")
top5pMorph <- read.csv("/Users/sebmatlosz/Desktop/filesforGOanalysis/glm_signif_Morphs_090222.csv")
    #OR
#top5pTimepoint <- read.csv("/Users/sebmatlosz/Desktop/filesforGOanalysis/glm_top5p_Timepoints.csv")
#top5pTimepoint <- read.csv("/Users/sebmatlosz/Desktop/filesforGOanalysis/glm_signif_Timepoints_090222.csv")
    #OR
#top5pSex <- read.csv("/Users/sebmatlosz/Desktop/filesforGOanalysis/glm_top5p_Sex.csv")

#Actually for once I need to get them in chr format
top5pMorph$NCBI <- str_replace(top5pMorph$NCBI, "NW_", "chr")
top5pMorph$NCBI <- str_replace(top5pMorph$NCBI, "NC_", "chr")

          #Or load the data of the overlaps
          #overlap5p <- read.csv("/Users/sebmatlosz/Desktop/filesforGOanalysis/glm_pca_5poverlap_morphs.csv")
          #OR
          #overlap5p <- read.csv("/Users/sebmatlosz/Desktop/filesforGOanalysis/glm_pca_5poverlap_timepoints.csv")
          #OR
          #overlap5p <- read.csv("/Users/sebmatlosz/Desktop/filesforGOanalysis/glm_pca_5poverlap_sex.csv")
          
          overlap5p$wide_smalldatatable.ID <- str_replace(overlap5p$wide_smalldatatable.ID, "NW_", "chr")
          overlap5p$wide_smalldatatable.ID <- str_replace(overlap5p$wide_smalldatatable.ID, "NC_", "chr")
          
          #OR the PCA data
          overlap5p <- read.csv("/Users/sebmatlosz/Desktop/filesforGOanalysis/PC1outliers_2sigmas_100222.csv")
          #overlap5p <- read.csv("/Users/sebmatlosz/Desktop/filesforGOanalysis/PC4outliers_2sigmas_100222.csv")
          
          
          #Separate scaffold from start
          library(stringr)
          names <- as.data.frame(str_split_fixed(overlap5p$ID, "_",2))
          overlap5p$chrom <- names$V1
          overlap5p$start <- names$V2
          overlap5p$start <- as.numeric(as.character(overlap5p$start))

#Now get distance to tss
df <- overlap5p

library(GenomicRanges)
library(genomation)

gene.obj=readTranscriptFeatures("/Users/sebmatlosz/Documents/annotated_charr_genome/goodcharrgenome.bed")

#ATH make sure you specify the correct columns, it can change depending on the df
names(df)[50] <- "chrom"
names(df)[51] <- "start"
df$start <- as.numeric(df$start)
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

#How many of them are on unplaced ? i.e. with NA in dist to feature.
noNA <- betterdf[complete.cases(betterdf$disttofeature),] # top5pMorphs: 77 on unplaced / top5pTimepoint: 67 / top5pSex: 39
#Then out of those, how many of them are less than 10kb away
noNA$disttofeature <- as.numeric(as.character(noNA$disttofeature))
noNA$absDisttoTSS <- abs(noNA$disttofeature)
lessthan10kb <- filter(noNA, absDisttoTSS <10000) # top5pMorphs: 180 / top5pTimepoint: 214 / top5pSex: 192

#Now, let's fine the ID of all unique genes that are <10kb away from those residues.
library(sqldf)
toppar <- data.frame("chrscaffold"=lessthan10kb$chrscaffold, "Tilestart"=lessthan10kb$windowstart)
#need to change the "chr" string into a "NW" or "NC" string
toppar$temporary <- gsub("chr","",toppar$chrscaffold)
toppar$NW <- ifelse(nchar(toppar$temporary)==11, "NW_", "NC_")
toppar$Scaffold<-paste(toppar$NW,toppar$temporary, sep="")
#Load some big data
stortData <- read.table(file='/Users/sebmatlosz/Desktop/canada_genome_overview_protein.tsv', sep = '\t', header=TRUE)

toppar$Tilestart <- as.numeric(as.character(toppar$Tilestart))
toppar$nedra <- toppar$Tilestart-10000
toppar$efra <- toppar$Tilestart+10000
toppar$nedra[toppar$nedra<0] <- 1

nytafla <- sqldf('SELECT chromosome_RefSeq, transcript_id, gene_id, start, stop, Tilestart,nedra, efra 
                 FROM stortData, toppar 
                 WHERE chromosome_RefSeq=Scaffold AND ((start<nedra AND stop>nedra) OR (start<efra AND stop>efra) OR (start>nedra AND stop<efra) OR (start<nedra AND stop>efra))')

listi <- nytafla$gene_id
listi <- sort(listi)
listi <- unique(listi)
length(listi)  #uniq genes: top5pMorph: 102 / top5pTimepoint: 83 / top5pSex: 143

#Now, let's do some GO analysis on those genes.
library(topGO)
library(dplyr)
library(genefilter)
library(Rgraphviz)
library(sqldf)

canada <- readMappings(file = "/Users/sebmatlosz/Desktop/geneidgo.db")

# oll gene_id sem eru ?? listanum
genanofn <- names(canada)
head(genanofn)

#listi is generated above
genalisti <- factor(as.integer(genanofn %in% listi)) 
names(genalisti) <- genanofn

# gerum GO-greiningu a gognunum             #But shouldn't I do this only on the rows where x=1 ? Actually no I think it uses this to calculate the significance
# gerum rad fyrir ad gera BP (getum notad CC eda MF lika)
GOdata <- new("topGOdata", description = "Genes located close to the residues top5pMorph in glm", 
                       ontology = "BP", allGenes = genalisti, nodeSize = 5,
                       annot = annFUN.gene2GO, gene2GO = canada)

#### GO-test ####
# Fisher - classic or weight test. Weight is better because if groups the nodes that relate to each other
resultFisher <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
resultFisher

allRes <- GenTable(GOdata, weightedFisher = resultFisher, topNodes=13)
write.table(allRes, "/Users/sebmatlosz/Desktop/filesforGOanalysis/GOresults_PC1_2sigma_270222.tsv", quote = FALSE, sep="\t", row.names = FALSE)





########################################################################
# allresPC1 <- read.table("/Users/sebmatlosz/Documents/toppar_GO_PC1_020621.tsv", sep="\t")
# allresPC1Trunc <- allresPC1[1:9,]
# write.table(allresPC1Trunc, "/Users/sebmatlosz/Documents/toppar_GO_PC1trunc_020621.tsv", quote = FALSE, sep="\t", row.names = FALSE)


#Cannot do weight test on KS
# Kolmogorov-Smirnov (KS) - classic test
resultKSPC1 <- runTest(GOdataPC1, algorithm = "classic", statistic = "ks")
resultKSPC1

##Order the results
allResPC1 <- GenTable(GOdataPC1, classicFisher = resultFisherPC1, classisKS = resultKSPC1, 
                            orderBy = "classicFisher", ranksOf = "classisKS", topNodes = 29)

#Save this as a .tsv for further look
write.table(allResPC1, "/Users/sebmatlosz/Documents/toppar_GO_PC1.tsv", quote = FALSE, sep="\t", row.names = FALSE)

########################################################################
#Check which genes are the ones that come up as Nucleosome Assembly:
genenames <- as.character(unique(nytafla$gene_id))
#Sort the canada list to only keep the elements that match the gene names
interestinggenes <- canada[genenames]
#Filter the elements where I find "GO:0006334"
#GO <- c("GO:0006334","GO:0006335","GO:0006336","GO:0034080","GO:1903098","GO:1903099")
test <- interestinggenes[lengths(Map(intersect, "GO:0006334", interestinggenes)) > 0]

#Result:
#GO:0006334:
    # 111971647 Histone H3-like
    # 111971885 Histone H1
    # 111971649 Histone H2B
    # 112080403 Histone H3-like
#GO:0006335 and GO:0006335:
    # 111971341 Histone H4

vectorgenes <- c("111971647","111971885","111971649","112080403","111971341")
#Find which gene corresponds to this
goodgenes <- nytafla[nytafla$gene_id %in% vectorgenes,]


#All five of them are histone genes.. But not exactly the same loci ?