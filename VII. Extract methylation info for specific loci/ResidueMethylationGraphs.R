library(tidyverse) 
library(plyr) 
library(ggplot2) 
library(plotly) 
library(ggpubr) 


#You have to work one gene at a time
table <- read.table("~/HistoneCluster.covdf", sep="\t", header=TRUE)

#table <- read.table("~/Lmtk2.covdf", sep="\t", header=TRUE) #LMTK2
#table <- read.table("~/SLC9A3R2.covdf", sep="\t", header=TRUE) #SLC9A3R2
#table <- read.table("~/Nkx23.covdf", sep="\t", header=TRUE) #NKX23
#table <- read.table("~/NFIX.covdf", sep="\t", header=TRUE) #NFIX
#table <- read.table("~/RASSF4.covdf", sep="\t", header=TRUE) #RASSF4
#table <- read.table("~/ARL16.covdf", sep="\t", header=TRUE) #ARL16
#table <- read.table("~/Meis1.covdf", sep="\t", header=TRUE) #Meis1
#table <- read.table("~/LOC111971648.covdf", sep="\t", header=TRUE) #Histone H2Al
#table <- read.table("~/LOC111971647.covdf", sep="\t", header=TRUE) #Histone H3l
#table <- read.table("~/ARMC1.covdf", sep="\t", header=TRUE) #ARMC1
#table <- read.table("~/Gli3.covdf", sep="\t", header=TRUE) #GLI3
#table <- read.table("~/MEGF9.covdf", sep="\t", header=TRUE) #MEGF9
#table <- read.table("~/MAGUKp55.covdf", sep="\t", header=TRUE) #MPP3
#table <- read.table("~/Rhoguanin37.covdf", sep="\t", header=TRUE) #ARHGEF37


unmeth <- table %>% select(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59,61,63,65,67,69,71,73,75,77,79,81,83,85,87,89,91,93,95,97)
meth <- table %>% select(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96)
methmat <- as.matrix(meth)
unmethmat <- as.matrix(unmeth)
total <- methmat + unmethmat

#create matrix percentage. Here I create methylation percentage for each individual sample. And then I will take the mean between them later.
#This is done instead of adding the raw numbers between replicates first and doing the percentage on those totals later.
#In my head both methods work but I chose this one.
perc <- (methmat*100)/total
#Reset the good positions
perc <- as.data.frame(perc)
perc$POS <- table$POS

#Specify the positions, i.e. start and end of the region of interest
template <- data.frame("POS"=50963000:50970000) #For Histone Cluster

#template <- data.frame("POS"=7669001:7670001) #LMTK2 on NC_036841.1
#template <- data.frame("POS"=5001:6001) #SLC9A3R2 on NW_019946553.1
#template <- data.frame("POS"=49001:50001) #NKX23 on NW_019944284.1
#template <- data.frame("POS"=163001:164001) #NFIX on NW_019942922.1
#template <- data.frame("POS"=17000:20000) #RASSF4 on NW_019957652.1
#template <- data.frame("POS"=4001:5001) #ARL16 on NW_019948814.1
#template <- data.frame("POS"=74001:75001) #Meis1 on NW_019944660.1
#template <- data.frame("POS"=50965001:50966001) #Histone H2Al on NC_036853.1
#template <- data.frame("POS"=50964001:50965001)  #Histone H3l on NC_036853.1
#template <- data.frame("POS"=2055001:2056001) #ARMC1 on NC_036867.1 
#template <- data.frame("POS"=13665001:13666001) #GLI3 on NC_036867.1
#template <- data.frame("POS"=1104001:1105001) #MEGF9 on NW_019957565.1
#template <- data.frame("POS"=7897001:7898001) #MPP3 on NC_036860.1
#template <- data.frame("POS"=26482001:26483001) #ARHGEF37 on NC_036845.1


total <- merge(template,perc,by="POS",all.x=T)
average <- data.frame("POS"=total$POS)

LB200 <- data.frame(total$X17LB1200meth,total$X17LB2200meth, total$X17LB3200meth)
average$LB200= rowMeans(LB200, na.rm = TRUE) 
LB150 <- data.frame(total$X17LB1150meth,total$X17LB2150meth, total$X17LB3150meth)
average$LB150 = rowMeans(LB150, na.rm = TRUE)
LB100 <- data.frame(total$X17LB1100meth,total$X17LB2100meth, total$X17LB3100meth)
average$LB100 = rowMeans(LB100, na.rm = TRUE)
LB50 <- data.frame(total$X17LB150meth,total$X17LB250meth, total$X17LB350meth)
average$LB50 = rowMeans(LB50, na.rm = TRUE)

SB200 <- data.frame(total$X17SB1200meth,total$X17SB3200meth, total$X17SB35200meth)
average$SB200 = rowMeans(SB200, na.rm = TRUE)
SB150 <- data.frame(total$X17SB1150meth,total$X17SB3150meth, total$X17SB35150meth)
average$SB150 = rowMeans(SB150, na.rm = TRUE)
SB100 <- data.frame(total$X17SB20100meth,total$X17SB3100meth, total$X17SB35100meth)
average$SB100 = rowMeans(SB100, na.rm = TRUE)
SB50 <- data.frame(total$X17SB150meth,total$X17SB350meth, total$X17SB3550meth)
average$SB50 = rowMeans(SB50, na.rm = TRUE)

PI200 <- data.frame(total$X17PI1200meth,total$X17PI3200meth, total$X17PI5200meth)
average$PI200 = rowMeans(PI200, na.rm = TRUE)
PI150 <- data.frame(total$X17PI1150meth,total$X17PI3150meth, total$X17PI5150meth)
average$PI150 = rowMeans(PI150, na.rm = TRUE)
PI100 <- data.frame(total$X17PI1100meth,total$X17PI3100meth, total$X17PI5100meth)
average$PI100 = rowMeans(PI100, na.rm = TRUE)
PI50 <- data.frame(total$X17PI150meth,total$X17PI350meth, total$X17PI550meth)
average$PI50 = rowMeans(PI50, na.rm = TRUE)

PL200 <- data.frame(total$X17PL30200meth,total$X17PL53200meth, total$X17PL54200meth)
average$PL200 = rowMeans(PL200, na.rm = TRUE)
PL150 <- data.frame(total$X17PL54150meth,total$X17PL30150meth, total$X17PL53150meth)
average$PL150 = rowMeans(PL150, na.rm = TRUE)
PL100 <- data.frame(total$X17PL30100meth,total$X17PL34100meth, total$X17PL54100meth)
average$PL100 = rowMeans(PL100, na.rm = TRUE)
PL50 <- data.frame(total$X17PL2350meth,total$X17PL2550meth, total$X17PL1150meth)
average$PL50 = rowMeans(PL50, na.rm = TRUE)


#Now the other thing is that I would like to know when I have averages that are only from one or two replicates.
#And then identify those differently on the graph.
#Always calculate the average before putting the shape
# ATH: this is outdated in the sense that I do not use the shape in graphs. But you still have to run it!

setting_shape <- function(x) {
  x$shape <- rowSums(is.na(x))
  x$shape <- gsub("1","partial",x$shape)
  x$shape <- gsub("2","partial",x$shape)
  x$shape <- gsub("3","partial",x$shape)
  x$shape <- gsub("0","full",x$shape)
  return(x)
}

average <- data.frame("POS"=total$POS)

LB200 <- data.frame(total$X17LB1200meth,total$X17LB2200meth, total$X17LB3200meth)
average$LB200= rowMeans(LB200, na.rm = TRUE) 
LB200 <- setting_shape(LB200)
LB150 <- data.frame(total$X17LB1150meth,total$X17LB2150meth, total$X17LB3150meth)
average$LB150 = rowMeans(LB150, na.rm = TRUE)
LB150 <- setting_shape(LB150)
LB100 <- data.frame(total$X17LB1100meth,total$X17LB2100meth, total$X17LB3100meth)
average$LB100 = rowMeans(LB100, na.rm = TRUE)
LB100 <- setting_shape(LB100)
LB50 <- data.frame(total$X17LB150meth,total$X17LB250meth, total$X17LB350meth)
average$LB50 = rowMeans(LB50, na.rm = TRUE)
LB50 <- setting_shape(LB50)

SB200 <- data.frame(total$X17SB1200meth,total$X17SB3200meth, total$X17SB35200meth)
average$SB200 = rowMeans(SB200, na.rm = TRUE)
SB200 <- setting_shape(SB200)
SB150 <- data.frame(total$X17SB1150meth,total$X17SB3150meth, total$X17SB35150meth)
average$SB150 = rowMeans(SB150, na.rm = TRUE)
SB150 <- setting_shape(SB150)
SB100 <- data.frame(total$X17SB20100meth,total$X17SB3100meth, total$X17SB35100meth)
average$SB100 = rowMeans(SB100, na.rm = TRUE)
SB100 <- setting_shape(SB100)
SB50 <- data.frame(total$X17SB150meth,total$X17SB350meth, total$X17SB3550meth)
average$SB50 = rowMeans(SB50, na.rm = TRUE)
SB50 <- setting_shape(SB50)

PI200 <- data.frame(total$X17PI1200meth,total$X17PI3200meth, total$X17PI5200meth)
average$PI200 = rowMeans(PI200, na.rm = TRUE)
PI200 <- setting_shape(PI200)
PI150 <- data.frame(total$X17PI1150meth,total$X17PI3150meth, total$X17PI5150meth)
average$PI150 = rowMeans(PI150, na.rm = TRUE)
PI150 <- setting_shape(PI150)
PI100 <- data.frame(total$X17PI1100meth,total$X17PI3100meth, total$X17PI5100meth)
average$PI100 = rowMeans(PI100, na.rm = TRUE)
PI100 <- setting_shape(PI100)
PI50 <- data.frame(total$X17PI150meth,total$X17PI350meth, total$X17PI550meth)
average$PI50 = rowMeans(PI50, na.rm = TRUE)
PI50 <- setting_shape(PI50)

PL200 <- data.frame(total$X17PL30200meth,total$X17PL53200meth, total$X17PL54200meth)
average$PL200 = rowMeans(PL200, na.rm = TRUE)
PL200 <- setting_shape(PL200)
PL150 <- data.frame(total$X17PL54150meth,total$X17PL30150meth, total$X17PL53150meth)
average$PL150 = rowMeans(PL150, na.rm = TRUE)
PL150 <- setting_shape(PL150)
PL100 <- data.frame(total$X17PL30100meth,total$X17PL34100meth, total$X17PL54100meth)
average$PL100 = rowMeans(PL100, na.rm = TRUE)
PL100 <- setting_shape(PL100)
PL50 <- data.frame(total$X17PL2350meth,total$X17PL2550meth, total$X17PL1150meth)
average$PL50 = rowMeans(PL50, na.rm = TRUE)
PL50 <- setting_shape(PL50)



 #If I use every single data point, then I might end up with some dots on the graph which are not present in all samples. 
      # This is not necessarily a bad thing, but it might cause problems on the graphs where there is a lot of coverage (eg Histones)
      #So if I want to keep only the ones that are covered in all positions: 
      LBtemporary <- data.frame("POS"=total$POS,"LB50"=average$LB50,"LB100"=average$LB100,"LB150"=average$LB150,"LB200"=average$LB200,
                                "LB50_shape"=LB50$shape,"LB100_shape"=LB100$shape,"LB150_shape"=LB150$shape,"LB200_shape"=LB200$shape)
     LBallcov <- LBtemporary[complete.cases(LBtemporary), ] #Use for regular graphs
     # LBallcov <- LBtemporary  #Use if you want to have every residue, no matter coverage
      
      LB50bis <- data.frame("POS"=LBallcov$POS,"data"=LBallcov$LB50,"coverage"=LBallcov$LB50_shape)
      LB50bis$stage <- '50ts'
      LB100bis <- data.frame("POS"=LBallcov$POS,"data"=LBallcov$LB100,"coverage"=LBallcov$LB100_shape)
      LB100bis$stage <- '100ts'
      LB150bis <- data.frame("POS"=LBallcov$POS,"data"=LBallcov$LB150,"coverage"=LBallcov$LB150_shape)
      LB150bis$stage <- '150ts'
      LB200bis <- data.frame("POS"=LBallcov$POS,"data"=LBallcov$LB200,"coverage"=LBallcov$LB200_shape)
      LB200bis$stage <- '200ts'
      #Add data.frames on top of each other  
      LB <- rbind (LB50bis,LB100bis,LB150bis,LB200bis)
      
      #SB
      SBtemporary <- data.frame("POS"=total$POS,"SB50"=average$SB50,"SB100"=average$SB100,"SB150"=average$SB150,"SB200"=average$SB200,
                                "SB50_shape"=SB50$shape,"SB100_shape"=SB100$shape,"SB150_shape"=SB150$shape,"SB200_shape"=SB200$shape)
      SBallcov <- SBtemporary[complete.cases(SBtemporary), ]
      #SBallcov <- SBtemporary #Use if you want to have every residue, no matter coverage
      
      SB50bis <- data.frame("POS"=SBallcov$POS,"data"=SBallcov$SB50,"coverage"=SBallcov$SB50_shape)
      SB50bis$stage <- '50ts'
      SB100bis <- data.frame("POS"=SBallcov$POS,"data"=SBallcov$SB100,"coverage"=SBallcov$SB100_shape)
      SB100bis$stage <- '100ts'
      SB150bis <- data.frame("POS"=SBallcov$POS,"data"=SBallcov$SB150,"coverage"=SBallcov$SB150_shape)
      SB150bis$stage <- '150ts'
      SB200bis <- data.frame("POS"=SBallcov$POS,"data"=SBallcov$SB200,"coverage"=SBallcov$SB200_shape)
      SB200bis$stage <- '200ts'
      #Add data.frames on top of each other  
      SB <- rbind (SB50bis,SB100bis,SB150bis,SB200bis)
      
      #PI
      PItemporary <- data.frame("POS"=total$POS,"PI50"=average$PI50,"PI100"=average$PI100,"PI150"=average$PI150,"PI200"=average$PI200,
                                "PI50_shape"=PI50$shape,"PI100_shape"=PI100$shape,"PI150_shape"=PI150$shape,"PI200_shape"=PI200$shape)
      PIallcov <- PItemporary[complete.cases(PItemporary), ]
      #PIallcov <- PItemporary #Use if you want to have every residue, no matter coverage
      
      PI50bis <- data.frame("POS"=PIallcov$POS,"data"=PIallcov$PI50,"coverage"=PIallcov$PI50_shape)
      PI50bis$stage <- '50ts'
      PI100bis <- data.frame("POS"=PIallcov$POS,"data"=PIallcov$PI100,"coverage"=PIallcov$PI100_shape)
      PI100bis$stage <- '100ts'
      PI150bis <- data.frame("POS"=PIallcov$POS,"data"=PIallcov$PI150,"coverage"=PIallcov$PI150_shape)
      PI150bis$stage <- '150ts'
      PI200bis <- data.frame("POS"=PIallcov$POS,"data"=PIallcov$PI200,"coverage"=PIallcov$PI200_shape)
      PI200bis$stage <- '200ts'
      #Add data.frames on top of each other  
      PI <- rbind (PI50bis,PI100bis,PI150bis,PI200bis)
      
      #PL
      PLtemporary <- data.frame("POS"=total$POS,"PL50"=average$PL50,"PL100"=average$PL100,"PL150"=average$PL150,"PL200"=average$PL200,
                                "PL50_shape"=PL50$shape,"PL100_shape"=PL100$shape,"PL150_shape"=PL150$shape,"PL200_shape"=PL200$shape)
      PLallcov <- PLtemporary[complete.cases(PLtemporary), ]
      #PLallcov <- PLtemporary #Use if you want to have every residue, no matter coverage
      
      PL50bis <- data.frame("POS"=PLallcov$POS,"data"=PLallcov$PL50,"coverage"=PLallcov$PL50_shape)
      PL50bis$stage <- '50ts'
      PL100bis <- data.frame("POS"=PLallcov$POS,"data"=PLallcov$PL100,"coverage"=PLallcov$PL100_shape)
      PL100bis$stage <- '100ts'
      PL150bis <- data.frame("POS"=PLallcov$POS,"data"=PLallcov$PL150,"coverage"=PLallcov$PL150_shape)
      PL150bis$stage <- '150ts'
      PL200bis <- data.frame("POS"=PLallcov$POS,"data"=PLallcov$PL200,"coverage"=PLallcov$PL200_shape)
      PL200bis$stage <- '200ts'
      #Add data.frames on top of each other  
      PL <- rbind (PL50bis,PL100bis,PL150bis,PL200bis)
      
      
#Once those four objects are made, it is time to build the graphs.
#Need to reshape data from long to wide format
Reshape <- function(data) {
  x <- reshape(data, idvar = "POS", timevar = "stage", direction = "wide")
  return(x)
}

LBwide <- Reshape(LB)
SBwide <- Reshape(SB)
PIwide <- Reshape(PI)
PLwide <- Reshape(PL)

#Let's make graphs with colourblind palette
colors <- c("050ts" = "#3CA6D0", "100ts" = "#4CCFFA", "150ts" = "#D06D21", "200ts" = "#7B3514" )            
#Create function for plotting.  add shape=coverage in the aes if I want the partial/full
plotting <- function(data,titles) {
  xaxis= paste("Position on NC_036853.1")
  
  x <- ggplot(data, aes(x = POS)) + 
  geom_point(aes(y=data.50ts, color="050ts")) + geom_point(aes(y=data.100ts, color="100ts")) +
  geom_point(aes(y=data.150ts, color="150ts")) + geom_point(aes(y=data.200ts, color="200ts")) +
  
  labs(x=xaxis, y= "Methylation %age",title = titles, color = "Stage")+
  ylim(0,100)+xlim(50963000,50970000) + scale_color_manual(values = colors) + theme_bw() + theme(legend.title=element_text(size=13), 
                                                                                                 legend.text=element_text(size=10))

  return(x)
}


LBgraph <- plotting(LBwide,"Histone Cluster / LB")       
SBgraph <- plotting(SBwide,"Histone Cluster / SB")
PLgraph <- plotting(PLwide,"Histone Cluster / PL")          
PIgraph <- plotting(LBwide,"Histone Cluster / PI")           

figure <- ggarrange(LBgraph,SBgraph,PIgraph,PLgraph,labels = c("","","",""),
                    ncol=2,nrow=2, common.legend = TRUE, legend = "right")

figure
pdf("~/HistoneCluster_colorblind_270222.pdf",15,7)
print(figure)
dev.off()

#Rince and repeat for each gene

