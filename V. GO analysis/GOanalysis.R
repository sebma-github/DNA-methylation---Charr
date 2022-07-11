#### Look at the residues with the lowest pvalue in glm, that explain a variable  ####
library(stringr) #in case you need to separate the start and chromosome strings from ID numbers
library(dplyr)
library(GenomicRanges)
library(genomation)
library(sqldf)

#Load the data of interest, for example here, the CpGs statistically significant btw morphs in the glm analysis.
    top5pMorph <- read.csv("~/glm_signif_Morphs.csv")

#Replace the "NW_" and "NC_" strings into "chr"
    top5pMorph$NCBI <- str_replace(top5pMorph$NCBI, "NW_", "chr")
    top5pMorph$NCBI <- str_replace(top5pMorph$NCBI, "NC_", "chr")

#Get distance to transcription start site
    df <- top5pMorph
    gene.obj=readTranscriptFeatures("~/goodcharrgenome.bed") #get the annotation from the Salvelinus sp. assembly, make sure you change it into .bed format
    names(df)[9] <- "chrom"   #ATH make sure you specify the correct columns, it can change depending on the df!
    names(df)[8] <- "start"   #ATH make sure you specify the correct columns, it can change depending on the df!
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

# Mine the important data
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

# At this point, the two df won't have the same number of rows because of loss of data in the annotation if on unplaced scaffold.
# Need to put them back together correctly
    x=resdf
    y=resumediffAnndf
    x=as.data.frame(x)
    y<- as.data.frame(y)
    x<-cbind(x,row.names(x))

betterdf <- merge(x,y,by.y="targetrow",by.x="row.names(x)",all=TRUE)

# Remove the CpGs that are on unplaced scaffolds (i.e. with NA in dist to feature).
    noNA <- betterdf[complete.cases(betterdf$disttofeature),]
# Out of the remaining ones look for genes that are less than 10kb away
    noNA$disttofeature <- as.numeric(as.character(noNA$disttofeature))
    noNA$absDisttoTSS <- abs(noNA$disttofeature)
    lessthan10kb <- filter(noNA, absDisttoTSS <10000)

# Get the ID of all unique genes that are <10kb away from those CpGs.
    toppar <- data.frame("chrscaffold"=lessthan10kb$chrscaffold, "Tilestart"=lessthan10kb$windowstart)
#need to change the "chr" string into a "NW" or "NC" string (fortunately, NC and NW scaffolds always have strings of different lengths)
    toppar$temporary <- gsub("chr","",toppar$chrscaffold)
    toppar$NW <- ifelse(nchar(toppar$temporary)==11, "NW_", "NC_")
    toppar$Scaffold<-paste(toppar$NW,toppar$temporary, sep="")

#Load Salvelinus sp. protein data
    stortData <- read.table(file='~/canada_genome_overview_protein.tsv', sep = '\t', header=TRUE)

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

#At this point, we should have the list of all unique genes <10kb away from the CpGs of interest.

#################################################
#Now, let's do some GO analysis on those genes.
library(topGO)
library(dplyr)
library(genefilter)
library(Rgraphviz)
library(sqldf)

#Load the reference of gene ontology for genes in this assembly
    canada <- readMappings(file = "~/geneidgo.db")
    genenames <- names(canada)
    genalisti <- factor(as.integer(genenames %in% listi)) #listi is generated above
    names(genalisti) <- genenames

    GOdata <- new("topGOdata", description = "Genes located close to the residues top5pMorph in glm", 
                       ontology = "BP", allGenes = genalisti, nodeSize = 5,
                       annot = annFUN.gene2GO, gene2GO = canada)

#### GO-test ####
# Fisher - classic or weight test. Weight is better because it groups the nodes that relate to each other
    resultFisher <- runTest(GOdata, algorithm = "weight", statistic = "fisher")
    resultFisher

    allRes <- GenTable(GOdata, weightedFisher = resultFisher, topNodes=25)
    write.table(allRes, "~/GOresults_glmsignif_morphs_270222.tsv", quote = FALSE, sep="\t", row.names = FALSE)



