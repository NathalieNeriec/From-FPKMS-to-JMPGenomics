##PACKAGES : gtable
# install.packages("ggplot2")
# install.packages("gtable")
# install.packages("grid")
# install.packages("gridExtra")
# install.packages("reshape2")
# install.packages("scales")
# install.packages("dplyr")
# install.packages("Rcpp")

#Load libraries:
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(reshape2)
library(scales)
library(dplyr)
library(plyr)

#Get the data
mergedTPM<-read.csv("data-2018-02-06.csv")

# Underconstruction : make it with an automatic output folder to define here at the beginning
output<-"C:/Users/nn503/Dropbox/Desktop/FromTPMToJMP"
# #png(paste0(output,"/lowfiltergraph/",plotname,".png"))


#########################################################################################
### PLOT AVG VALUES OF EACH CONDITION TO DEFINE "LOW" EXPRESSION FOR EACH CONDITION######
#########################################################################################


###### 1 ___Removing TPM and FPKMS from the colnames ###################################
colnames(mergedTPM)<-gsub("TPM","",colnames(mergedTPM))
colnames(mergedTPM)<-gsub("FPKMS","",colnames(mergedTPM))
colnames(mergedTPM)
ncolPrimary<-ncol(mergedTPM)

####### 2 ___OUR CONDITIONS  _ UNDER CONSTRUCTION ######################################
#Find a way to Rep them from "everything before _"
#In the mean time, for now:
Conditions <- c("MEFFTOKno", "MEFFTOwt",
                "MEFshnRNPABHet","MEFshnRNPABKno",
                "OlineoFTOCtrol", "OlineoFTOsiRNA",
                "OlineohnRNPABCtrol","OlineohnRNPABKno")

####### 3 ___Adding average condition per column: ###################################### 

#FUNCTION TO ADD COLUMNS WITH AVERAGE PER CONDITION
Add_Avg_Col_percondi <- function(condition){
  #Store colnames
  columnames <- colnames(mergedTPM)
  #Find column containing control CTRL
  condition_columnname <- columnames[grep(condition,columnames)]
  #Get only CTRL subsets of counts 
  condition_columns <- mergedTPM[,condition_columnname]
  #calculaterow averages and insert into a new column
  mergedTPM$condition_Avg=rowMeans(condition_columns)
  #Change the name  of the newly added column to match the pattern
  colnames(mergedTPM)[colnames(mergedTPM) =="condition_Avg"]<-paste0(condition,"_avg")
  mergedTPM<-return(mergedTPM)
}

#Looping to add the column into mergedTPM 
for (condition in Conditions){
  mergedTPM<-Add_Avg_Col_percondi(condition)
}

######## 4___Plotting our data #########################################################
#New dataframe for ploting: Only Average
columnames <- colnames(mergedTPM)
avg_columnname <- columnames[grep("_avg",columnames)]
mergedTPM_avg <- mergedTPM[,avg_columnname]
mergedTPM_avg[mergedTPM_avg==0]<-NA

#Ranking ever ycolumn :
for (i in 1:ncol(mergedTPM_avg)){
  mergedTPM_avg[,i]<-sort(mergedTPM_avg[,i],na.last = TRUE)
}
# Plotting every column, into a new file with y valies of "yvalues"
mergedavgcolnameslist<-colnames(mergedTPM_avg)
yvalues<-c(10,50,100,1000)

for(namecondition in 1:length(mergedTPM_avg)){
    namey<-mergedTPM_avg[,namecondition]
    condition<-mergedavgcolnameslist[[namecondition]]
   for (maxofy in yvalues){
      plotname<-paste0("plot_",condition,"_",maxofy)
      png(paste0(output,"/lowfiltergraph/",plotname,".png"))
      print(
        ggplot(
          mergedTPM_avg,
          aes(x=rownames(mergedTPM_avg),
              y=namey),
              na.rm=TRUE)+
        ggtitle(plotname)+
        geom_point()+
        ylim(0,maxofy)+
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()))
      dev.off()
    }
  }

######### 4_ MAKING FOUR CATEGORIES #########################################################
######### Keep only if at least one condition has ALL over 5###################
interestingtreshold<-c(1,10)
for(treshold in interestingtreshold){
mergedTPM_backup<-mergedTPM
#recorind #columns before we add any
ncolWithAvg<-ncol(mergedTPM)

#Loop to add a FALSE/TRUE column per condition meeting or not "all above 5"
for (condition in Conditions) {
  #Get the columname matchig the condition
  columnames <- colnames(mergedTPM[,1:ncolPrimary])
  condition_columnname <- columnames[grep(condition,columnames)]
  #Add a newcolumn with the percentage in each condition above your treshold here 5
  mergedTPM$newcolumn<-
    apply(mergedTPM[, condition_columnname], 
          MARGIN = 1, 
          function(x) (sum(x>treshold)/length(condition_columnname)))
  #Renaming the new column on a condition based name so that you just add
  colnames(mergedTPM)[colnames(mergedTPM) =="newcolumn"]<-paste0("testabove",treshold,"_",condition)
}

#Gathering Filter  
mergedTPM_backup2<-mergedTPM
columnames2 <- colnames(mergedTPM)
testabove_columnname<-columnames2[grep("testabove",columnames2)]
mergedTPM<-
  mergedTPM[apply(mergedTPM[,testabove_columnname],MARGIN=1,function(x)any(x>0.50)),]

mergedTPM$test<-apply(mergedTPM[,testabove_columnname],MARGIN=1,function(x)any(x>0.50))

Ngeneswithatleast1<-nrow(mergedTPM)
write.csv(mergedTPM,file=paste0("mergedTPM_above",treshold,".csv"),row.names = FALSE)


######### 5_SUBSETTING #########################################################
######### Keep wether # conditions ABOVE OR BELOW 5###################
#gives me the number of TRUE!!

mergedTPM$NTotalConditionsabove0<-apply(mergedTPM[, testabove_columnname], MARGIN = 1, function(x) sum(x>0))
l<-list(1,length(Conditions))
z=1
NcolwithAvg_withSubset<-ncol(mergedTPM)

GenePerConditions<-matrix(nrow=length(Conditions),ncol=2)
colnames(GenePerConditions)<-c("NconditionAbove0","Ngenes")

for (i in 1:length(Conditions)){
  mergedTPM_new<-mergedTPM[mergedTPM$NTotalConditionsabove0 == i,c(1:ncolWithAvg,NcolwithAvg_withSubset)]
  GenePerConditions[i,]<-c(i,nrow(mergedTPM_new))
  write.csv(mergedTPM_new[,1:ncolPrimary],paste0("Filterwith_",i,"_conditionabove",treshold,".csv"),row.names = FALSE)
  assign(paste0("Filterwith_",i,"_conditionabove",treshold),mergedTPM_new)
  l[[z]]<-paste0("Filterwith_",i,"_conditionabove",treshold)
  z=z+1
}


# save the number of genes
write.csv(GenePerConditions,file=paste0("GenePerconditions_withtreshold",treshold,".csv"),row.names = FALSE)
}



# ######### 5_SUBSETTING #########################################################
# #SO we have 8 cases : 
# 
#  "Filterwith_1_conditionabove5" "Filterwith_2_conditionabove5"
#  "Filterwith_3_conditionabove5" "Filterwith_4_conditionabove5"
#  "Filterwith_5_conditionabove5" "Filterwith_6_conditionabove5"
#  "Filterwith_7_conditionabove5" "Filterwith_8_conditionabove5"
# 
#Express in ALL
 write.csv(Filterwith_8_conditionabove5[,1:ncolPrimary],paste0("ExpressinAll_above",treshold,".csv"),row.names = FALSE)


  # GenePerConditions<-matrix(nrow=length(Conditions),ncol=2)
colnames(GenePerConditions)<-c("NconditionAbove0","Ngenes")
i=1
for(i in length(l)){
  nrowCon<-paste0("Filterwith_",i,"_conditionabove5")
  GenePerConditions[1,]<-c(i,nrow(nrowCon))
  i=i+1
}


#Extract Olineo
columnames2 <- colnames(mergedTPM[1:ncolPrimary])
Olineo_columnname<-columnames2[grep("Olin",columnames2)]

Olineo<-mergedTPM[,c("gene.ids","gene.names",Olineo_columnname)]
ncolPrimaryOlineo<-ncol(Olineo)
Conditions=c("OlineoFTOCtrol","OlineoFTOsiRNA","OlineohnRNPABCtrol","OlineohnRNPABKno")
for (condition in Conditions) {
  #Get the columname matchig the condition
  columnames <- colnames(Olineo[,1:ncolPrimaryOlineo])
  condition_columnname <- columnames[grep(condition,columnames)]
  #Add a newcolumn with the percentage in each condition above your treshold here 5
  Olineo$newcolumn<-
    apply(Olineo[, condition_columnname], 
          MARGIN = 1, 
          function(x) (sum(x>treshold)/length(condition_columnname)))
  #Renaming the new column on a condition based name so that you just add
  colnames(Olineo)[colnames(Olineo) =="newcolumn"]<-paste0("testabove",treshold,"_",condition)
}

#Gathering Filter  
Olineo_backup2<-Olineo
columnames2 <- colnames(Olineo)
testabove_columnname<-columnames2[grep("testabove",columnames2)]
Olineo<-  Olineo[apply(Olineo[,testabove_columnname],MARGIN=1,function(x)any(x>0.50)),]

Olineo$test<-apply(Olineo[,testabove_columnname],MARGIN=1,function(x)any(x>0.50))

Ngeneswithatleast1<-nrow(Olineo)
write.csv(Olineo,file=paste0("Olineo_above",treshold,".csv"),row.names = FALSE)
library(xlsx)
write.xlsx(Olineo[,1:ncolPrimaryOlineo],file=paste0("Olineo_above",treshold,".csv"),row.names = FALSE)

#HOW TO REMOVE THE 
 
# #2. ADD +1 TO EVERYTHING
# mergedTPM<-cbind(mergedTPM[,1:2],(mergedTPM[,3:ncol(mergedTPM)]+1))
# 
# 
# library(ggplot2)
 
olineoAnovaResults<-read.xlsx("olineoabove10_onewayanovaresults.xlsx",sheetName = "Sheet1")
