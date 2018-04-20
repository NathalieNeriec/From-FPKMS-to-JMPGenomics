
########## libraries###########
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(reshape2)
library(scales)
library(dplyr)
library(plyr)
library(xlsx)

########### 1_ Prepare your input data ###############

#Read your input
MEFFTO<-read.csv("MEFFTO.csv")
MEFhnRNPAB<-read.csv("MEFhnRNPAB.csv")
OlineoFTO<-read.csv("OlineoFTO.csv")
OlineohnRNPAB<-read.csv("OlineohnRNPAB.csv")


#MEFOnly
MEfs<-cbind(MEFFTO,MEFhnRNPAB)

#make Xp1 and Xp2 of fto
OlineoFTORep2and4XP1<-OlineoFTO
OlineoFTORep1and3XP2<-OlineoFTO
OlineoFTORep1and3XP2$X<-NULL
OlineoFTORep2and4XP1$X<-NULL
OlineoFTORep1and3XP2$OlineoFTOCtrl_2<-NULL
OlineoFTORep1and3XP2$OlineoFTOCtrl_4<-NULL
OlineoFTORep1and3XP2$OlineoFTOsiRNA_2<-NULL
OlineoFTORep1and3XP2$OlineoFTOsiRNA_4<-NULL
OlineoFTORep2and4XP1$OlineoFTOCtrl_1<-NULL
OlineoFTORep2and4XP1$OlineoFTOCtrl_3<-NULL
OlineoFTORep2and4XP1$OlineoFTOsiRNA_1<-NULL
OlineoFTORep2and4XP1$OlineoFTOsiRNA_3<-NULL

dataframeofinterest<-OlineosWithoutRplicate3
names(dataframeofinterest)[names(dataframeofinterest) == 'OlineoFTOCtrl_1'] <- 'OlineoFTOCtrlXp1_1'
names(dataframeofinterest)[names(dataframeofinterest) == 'OlineoFTOCtrl_2'] <- 'OlineoFTOCtrlXp2_2'
names(dataframeofinterest)[names(dataframeofinterest) == 'OlineoFTOCtrl_3'] <- 'OlineoFTOCtrlXp2_3'
names(dataframeofinterest)[names(dataframeofinterest) == 'OlineoFTOCtrl_4'] <- 'OlineoFTOCtrlXp1_4'
names(dataframeofinterest)[names(dataframeofinterest) == 'OlineoFTOsiRNA_1'] <- 'OlineoFTOsiRNAXp1_1'
names(dataframeofinterest)[names(dataframeofinterest) == 'OlineoFTOsiRNA_2'] <- 'OlineoFTOsiRNAXp2_2'
names(dataframeofinterest)[names(dataframeofinterest) == 'OlineoFTOsiRNA_3'] <- 'OlineoFTOsiRNAXp1_3'
names(dataframeofinterest)[names(dataframeofinterest) == 'OlineoFTOsiRNA_4'] <- 'OlineoFTOsiRNAXp2_4'

Conditions <-c("MEFsFTOMinus","MEFsFTOPlus","MEFshnRNPABKno","MEFshnRNPABHet")
Conditionsofinterest<-Conditions
dataframeofinterest<-MEfs 
####### 2_Find your treshold : Plotting the data##########

#Get colnames
MEFFTO_colnames<-colnames(MEFFTO)
MEFhnRNPAB_colnames<-colnames(MEFhnRNPAB)
OlineoFTO_colnames<-colnames(OlineoFTO)
OlineohnRNPAB_colnames<-colnames(OlineohnRNPAB)

#Make them into one file and not divided by "comparisons"
cbindall<-cbind(MEFFTO,MEFhnRNPAB,OlineoFTO,OlineohnRNPAB)
colgeneidsnames<-cbindall[,1:2]
cbindall$gene.ids<-NULL
cbindall$gene.ids<-NULL
cbindall$gene.ids<-NULL
cbindall$gene.ids<-NULL
cbindall$gene.names<-NULL
cbindall$gene.names<-NULL
cbindall$gene.names<-NULL
cbindall$gene.names<-NULL
cbindall<-cbind(colgeneidsnames,cbindall)


#Save your mega dataframe under "MergedTPM and prepare other necessary input to run "TReshold"
mergedTPM<-cbindall
getwd()
output="C:/Users/nn503/Dropbox/Desktop/FromCountsToDEsEQ2"
dir.create(paste0(output,"/lowfiltergraph/"))
#Get your conditions
#UNder construction : How to get them "automatically"
colnames(cbindall[,3:ncol(cbindall)])
Conditions <- c("MEFsFTOMinus", "MEFsFTOPlus","MEFshnRNPABHet","MEFshnRNPABKno",
                "OlineoFTOCtrl", "OlineoFTOsiRNA","OlineohnRNPABCtrl","OlineohnRNPABsiRNA")

Conditions <- c("OlineoFTOCtrl", "OlineoFTOsiRNA","OlineohnRNPABCtrl","OlineohnRNPABsiRNA")


#Adding Avg condition columns 
    #define function
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

#New dataframe with only the average
    columnames <- colnames(mergedTPM)
    avg_columnname <- columnames[grep("_avg",columnames)]
    mergedTPM_avg <- mergedTPM[,avg_columnname]
    mergedTPM_avg[mergedTPM_avg==0]<-NA

#Ranking ever column :
    for (i in 1:ncol(mergedTPM_avg)){
      mergedTPM_avg[,i]<-sort(mergedTPM_avg[,i],na.last = TRUE)
    }
# Plotting every column, into a new file with y valies of "yvalues"
    mergedavgcolnameslist<-colnames(mergedTPM_avg)
    yvalues<-c(10,50,100,1000)

#Looping the plot and saing in /lowerfiltergraph/ plotname.png for the yvalues predefined
    #Takes a bit of time
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

####### 3_Isolate four categories#############
    
    #Necessary input
    
    # Isolate_four_categories("OlineoFTORep1and3XP2","OlineoFTOCtrl","OlineoFTOsiRNA",10)
    # Isolate_four_categories("OlineoFTORep2and4XP1","OlineoFTOCtrl","OlineoFTOsiRNA",10)
    # 
    #Function to automatically isolate the four categories and save as csv that you can open in DESEQ2 tsar
    
    dataframeofinterest<-mergedTPM
    Conditionsofinterest<-Conditions
      
  ncolofinterest<-ncol(dataframeofinterest)
  #MAke the new columns
  #recoring #columns before we add any
  ncolWithAvg<-ncol(dataframeofinterest)
  #Loop to add a FALSE/TRUE column per condition meeting or not "all above 5" and one with the avg
  for (condition in Conditionsofinterest) {
    #Get the columname matchig the condition
    columnames <- colnames(dataframeofinterest[,1:ncolofinterest])
    condition_columnname <- columnames[grep(condition,columnames)]
    #Add a newcolumn with the percentage in each condition above your treshold here 5
    dataframeofinterest$newcolumn<-
      apply(dataframeofinterest[, condition_columnname], 
            MARGIN = 1, 
            function(x) (sum(x>treshold)/length(condition_columnname)))
    colnames(dataframeofinterest)[colnames(dataframeofinterest) =="newcolumn"]<-paste0("testabove",treshold,"_",condition)
  }
  
  #Make a new data frame with the Filter  
  columnames2 <- colnames(dataframeofinterest)
  testabove_columnname<-columnames2[grep("testabove",columnames2)]
  #Only keep the rows where ANY have at least 50% above the treshold (remove categorya)
  dataframeofinterest<-
    dataframeofinterest[apply(dataframeofinterest[,testabove_columnname],MARGIN=1,function(x)any(x>0.50)),]
  Ngeneswithatleast1<-nrow(dataframeofinterest)
  
  # Add a column per condition witht he average
  for (condition in Conditionsofinterest) {
    columnames <- colnames(dataframeofinterest[,1:ncolofinterest])
    condition_columnname <- columnames[grep(condition,columnames)]
    dataframeofinterest$newcolumn2<-
      apply(dataframeofinterest[, condition_columnname], 
            MARGIN = 1, 
            function(x) (mean(x)))
    colnames(dataframeofinterest)[colnames(dataframeofinterest) =="newcolumn2"]<-paste0("avg_",condition)
  }
  
  #Getting names of avg columns
  columnames2 <- colnames(dataframeofinterest)
  avg_columnname<-columnames2[grep("avg",columnames2)]
  
  #Isolating the Category d ON/OFF GENES  
  #Only keep the where ANY of the avg is 0 into a NEWdataFrame
  dataframeofinterest_Categoryd_ONOFF<-
    dataframeofinterest[apply(
      dataframeofinterest[,avg_columnname],
      MARGIN=1,
      function(x)any(x==0)),
      1:ncolofinterest]
  NgenesONOFF<-nrow(dataframeofinterest_Categoryd_ONOFF)
  
  #Resave data frame wihtout the on off : only the b and c
  dataframeofinterest<-
    dataframeofinterest[apply(
      dataframeofinterest[,avg_columnname],
      MARGIN=1,
      function(x)!(any(x==0))),]

  
  #Only 
  dataframeofinterestonlyvalues<-dataframeofinterest[,1:ncolofinterest]
  write.csv(dataframeofinterestonlyvalues,"MEfstogether.csv",row.names = FALSE)
  library(xlsx)
  write.xlsx(dataframeofinterestonlyvalues,"Olineostogether_FTOSplitinXps_ABNoreplicates3.xlsx",row.names = FALSE)
  
  # Separate category b  TESTABLE
  #When they dont all>0.5
  dataframeofinterest_categoryb_testable<-
    dataframeofinterest[apply(
      dataframeofinterest[,testabove_columnname],
      MARGIN=1,
      function(x)!(all(x>0.5))),
      1:ncolofinterest]
  NgenesTestable<-nrow(dataframeofinterest_categoryb_testable)
  # Separate category c DE GENES
  #When they  all>0.5
  dataframeofinterest_categoryc_DEGenes<-
    dataframeofinterest[apply(
      dataframeofinterest[,testabove_columnname],
      MARGIN=1,
      function(x)(all(x>0.5))),
      1:ncolofinterest]
  NgenesDEGenes<-nrow(dataframeofinterest_categoryc_DEGenes)
  
  #Resavewith the new isolated data frame to a new name
  write.csv(dataframeofinterest_categoryb_testable,paste0(comparison,"_Categoryb_testable.csv"),row.names = FALSE)
  write.csv(dataframeofinterest_categoryc_DEGenes,paste0(comparison,"_Categoryc_DEGENES.csv"),row.names = FALSE)
  write.csv(dataframeofinterest_Categoryd_ONOFF,paste0(comparison,"_Categoryd_ONOFF.csv"),row.names = FALSE)
  
  test<-rbind(dataframeofinterest_categoryb_testable,dataframeofinterest_categoryc_DEGenes)
  write.csv(test,paste0(comparison,"_categorybc.csv"),row.names = FALSE)
}



Isolate_four_categories("MEFFTO","MEFsFTOPlus","MEFsFTOMinus",10)
Isolate_four_categories("MEFhnRNPAB","MEFshnRNPABHet","MEFshnRNPABKno",10)
Isolate_four_categories("OlineoFTO","OlineoFTOCtrl","OlineoFTOsiRNA",10)
Isolate_four_categories("OlineohnRNPAB","OlineohnRNPABCtrl","OlineohnRNPABsiRNA",10)


Comparisons=c("MEFFTO","MEFhnRNPAB","OlineoFTO","OlineohnRNPAB")
for (comparison in Comparisons){
  write.csv(get(paste0(comparison,"_categoryb_testable")),paste0(comparison,"_Categoryb_testable.csv"),row.names = FALSE)
  write.csv(get(paste0(comparison,"_categoryc_DEGenes")),paste0(comparison,"_Categoryc_DEGENES.csv"),row.names = FALSE)
  write.csv(get(paste0(comparison,"_categoryd_ONOFF")),paste0(comparison,"_Categoryd_ONOFF.csv"),row.names = FALSE)
  
  test<-rbind(get(paste0(comparison,"_categoryb_testable")),get(paste0(comparison,"_categoryc_DEGenes")))
  write.csv(test,paste0(comparison,"_categorybc.csv"),row.names = FALSE)
}