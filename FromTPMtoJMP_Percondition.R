library(xlsx)
library(haven)
library(UpSetR)

library("RColorBrewer")


setwd("C:/Users/nn503/Dropbox/Desktop/ComparisonMagnus/FromTPMtoJMP_Percondition/Combined")

MEFFTO<-read.csv("./OriginalMergedTPM/MEF_FTO.csv")
MEFhnRNPAB<-read.csv("./OriginalMergedTPM/MEF_AB.csv")
OlineoFTO<-read.csv("./OriginalMergedTPM/Olineo_FTO.csv")
OlineohnRNPAB<-read.csv("./OriginalMergedTPM/Olineo_AB.csv")

###############1 Find the treshold IS SOMEWHERE ELSE
###############2 Remove low expressed

comparison="MEFFTO"
condition1="MEFFTOwt"
condition2="MEFFTOKno"

Removelowexpressed<-function(comparison,condition1,condition2,treshold){
  dataframeofinterest<-get(comparison)
  Conditionsofinterest<-c(condition1,condition2)
  ncolofinterest<-ncol(dataframeofinterest)
  treshold=1
  
  #Make the new columns with x>treshold 
        for (condition in Conditionsofinterest) {
          #Get the columname matchig the condition
            columnames <- colnames(dataframeofinterest[,1:ncolofinterest])
            condition_columnname <- columnames[grep(condition,columnames)]
            ncolofinterest<-ncol(dataframeofinterest)
            
          #Add a newcolumn with the percentage in each condition above your treshold here 5
            dataframeofinterest$newcolumn<-
                apply(dataframeofinterest[,condition_columnname], 
                      MARGIN = 1, 
                      function(x) (sum(x>treshold)/length(condition_columnname)))
            colnames(dataframeofinterest)[colnames(dataframeofinterest) =="newcolumn"]<-paste0("testabove",treshold,"_",condition)
        }
  
  #Remove the not expressed
        columnames2 <- colnames(dataframeofinterest)
        testabove_columnname<-columnames2[grep("testabove",columnames2)]
        
  #Only keep the rows where ANY have at least 50% above the treshold (remove categorya)
        dataframeofinterest2<-
            dataframeofinterest[
                apply(dataframeofinterest[,testabove_columnname],
                      MARGIN=1,
                      function(x)any(x>0.50)),]
        Ngeneswithatleast1<-nrow(dataframeofinterest)

  write.table(dataframeofinterest2,paste0(comparison,"_tresholdof_",treshold,".csv"),
            sep=",",
            row.names = FALSE,
            quote=FALSE)
  dataframeofinterest2<-dataframeofinterest2[,1:ncolofinterest]
  #dataframeofinterest > all without the low expressed > Theone you want
  log2dataframeofinterest <-cbind(dataframeofinterest2[,1:2],
                                   (dataframeofinterest2[,3:ncolofinterest]+1))
  log2dataframeofinterest2<- cbind(log2dataframeofinterest[,1:2],
                                     log(log2dataframeofinterest[3:ncolofinterest],2) )        
  assign(paste0(comparison,"_Plus1_log2transformed"))     
  ass
  write.table(log2dataframeofinterest2,paste0(comparison,"_Plus1_log2transformed.csv"),
              quote=FALSE,
              row.names = FALSE)
   
}
  
Removelowexpressed("MEFFTO","MEFFTOwt","MEFFTOKno",1)
Removelowexpressed("MEFhnRNPAB","MEFshnRNPABHet","MEFshnRNPABKno",1)
Removelowexpressed("OlineoFTO","OlineoFTOCtrol","OlineoFTOsiRNA",1)
Removelowexpressed("OlineohnRNPAB","OlineohnRNPABCtrol","OlineohnRNPABKno",1)

  
######### REMOVING LOW EXPRESS : TRESHOLD AND ISOLATE DIFFERENT CATEGORIES
#Function Isolate_four_Categories: 
# Define the "comparison" = the dataframe to work with
# condition1 and two : like "mutant" / "wt . Correspond to the common string between the columns
# treshold : the low value treshold to work with. 


Isolate_four_categories<-function(comparison,condition1,condition2,treshold){
  dataframeofinterest<-get(comparison)
  Conditionsofinterest<-c(condition1,condition2)
  ncolofinterest<-ncol(dataframeofinterest)
  
  #Make the new columns
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
  
  #Add a column per condition with he average
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
  
  ##### ISOLATING THE ON OFF GENES  
       #Only keep the where ANY of the avg is 0 into a new data frame
        dataframeofinterest_Categoryd_ONOFF<-
             dataframeofinterest[apply(
             dataframeofinterest[,avg_columnname],
                 MARGIN=1,
                 function(x)any(x==0)),
                 1:ncolofinterest]
        NgenesONOFF<-nrow(dataframeofinterest_Categoryd_ONOFF)
  
      #Remove  the onoff from the expressed genes (DE+ testable)
        dataframeofinterest<-
             dataframeofinterest[apply(
             dataframeofinterest[,avg_columnname],
             MARGIN=1,
             function(x)!(any(x==0))),]
        dataframeofinterest_ForJMP<-dataframeofinterest[,1:ncolofinterest]
        NgenesForJMP<-nrow(dataframeofinterest_ForJMP)
        
    #### LOOKING AT LOG2FC
        #Add columns to make the log2 fold change
            dataframeofinterest$log2FC<-
              log2(dataframeofinterest[,paste0("avg_",condition1)])-
              log2(dataframeofinterest[,paste0("avg_",condition2)])
       
        #remove theone with abs(log2)<1
            dataframeofinterest_ForJMP_withLog2FCAbove1<-
               dataframeofinterest[
                 which(
                    abs(dataframeofinterest[,ncol(dataframeofinterest)])>1),
                 ]
            dataframeofinterest_ForJMP_withLog2FCAbove1<-dataframeofinterest_ForJMP_withLog2FCAbove1[1:ncolofinterest]
               
            NgenesForJMP_log2Above1<-nrow(dataframeofinterest_ForJMP_withLog2FCAbove1)
               


  #### EXTRACT THE TESTABLE "b"category
        #When they dont all>0.5
        dataframeofinterest_categoryb_testable<-
             dataframeofinterest[apply(
             dataframeofinterest[,testabove_columnname],
             MARGIN=1,
             function(x)!(all(x>0.5))),
             1:ncolofinterest]
        NgenesTestable<-nrow(dataframeofinterest_categoryb_testable)
        
        
  #### EXTRACT THE DIFF EXPRESSED "c"category
       # When they  all>0.5
        dataframeofinterest_categoryc_DEGenes<-
        dataframeofinterest[apply(
        dataframeofinterest[,testabove_columnname],
          MARGIN=1,
          function(x)(all(x>0.5))),
          1:ncolofinterest]
         NgenesDEGenes<-nrow(dataframeofinterest_categoryc_DEGenes)
         
         

  
  #Make a table with the different number of genes for each conditions
         Ngenes<-c(NgenesONOFF,NgenesDEGenes, NgenesTestable,NgenesForJMP,NgenesForJMP_log2Above1)
         Ngenes<-as.data.frame(Ngenes)
        row.names(Ngenes)<-c("NgenesONOFF","NgenesDEGenes","NgenesTestable","NgenesForJMP","genesForJMP_log2Above1")
  
  #MAKE FILES OUT OF THE FUNCTION
  write.xlsx(dataframeofinterest_categoryb_testable,paste0(comparison,"_Categoryb_testable_tresholdof",treshold,".xlsx"),row.names = FALSE)
  write.xlsx(dataframeofinterest_categoryc_DEGenes,paste0(comparison,"_Categoryc_DEGENES_tresholdof",treshold,".xlsx"),row.names = FALSE)
  write.xlsx(dataframeofinterest_Categoryd_ONOFF,paste0(comparison,"_Categoryd_ONOFF_tresholdof",treshold,".xlsx"),row.names = FALSE)
  write.xlsx(dataframeofinterest_ForJMP,paste0(comparison,"_ForJMP_tresholdof",treshold,".xlsx"),row.names = FALSE)
  write.xlsx(Ngenes,paste0(comparison,"_NofGenes_tresholdof",treshold,".xlsx"),row.names=TRUE)
  write.xlsx(dataframeofinterest_ForJMP_withLog2FCAbove1,paste0(comparison," _ForJMP_withLog2FCAbove1_",treshold,".xlsx"),row.names=FALSE)

  }

for (Tresholdloop in c(1,5)){
Isolate_four_categories("MEFFTO","MEFFTOwt","MEFFTOKno",Tresholdloop)
Isolate_four_categories("MEFhnRNPAB","MEFshnRNPABHet","MEFshnRNPABKno",Tresholdloop)
Isolate_four_categories("OlineoFTO","OlineoFTOCtrol","OlineoFTOsiRNA",Tresholdloop)
Isolate_four_categories("OlineohnRNPAB","OlineohnRNPABCtrol","OlineohnRNPABKno",Tresholdloop)
}


