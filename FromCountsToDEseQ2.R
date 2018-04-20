MEFFTO<-read.csv("MEFFTO.csv")
MEFhnRNPAB<-read.csv("MEFhnRNPAB.csv")
OlineoFTO<-read.csv("OlineoFTO.csv")
OlineohnRNPAB<-read.csv("OlineohnRNPAB.csv")

###############1 Find the treshold
#Get colnames
MEFFTO_colnames<-colnames(MEFFTO)
MEFhnRNPAB_colnames<-colnames(MEFhnRNPAB)
OlineoFTO_colnames<-colnames(OlineoFTO)
OlineohnRNPAB_colnames<-colnames(OlineohnRNPAB)

#Make conditions and mergedTPM
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

colnames(cbindall)
Conditions <- c("MEFsFTOMinus", "MEFsFTOPlus","MEFshnRNPABHet","MEFshnRNPABKno",
                "OlineoFTOCtrl", "OlineoFTOsiRNA","OlineohnRNPABCtrl","OlineohnRNPABsiRNA")

mergedTPM<-cbindall
getwd()
output="C:/Users/nn503/Dropbox/Desktop/FromCountsToDEsEQ2"
dir.create(paste0(output,"/lowfiltergraph/"))



### RUN LIKE REMOVING LOW VALUES now > RUNNING #3 and #4 (as is!)

#Get the "Interesting treshold" 
#Now we want them separated per conditions to mAke four categories:

interestingtreshold<-c(1,10)

Comparisons=c("MEFFTO","MEFhnRNPAB","OlineoFTO","OlineohnRNPAB")
# Conditions = c
#  "MEFsFTOMinus"       "MEFsFTOPlus"        "MEFshnRNPABHet"     "MEFshnRNPABKno"    
# "OlineoFTOCtrl"      "OlineoFTOsiRNA"     "OlineohnRNPABCtrl"  "OlineohnRNPABsiRNA"
condition1="MEFsFTOPlus"
condition2="MEFsFTOMinus"
comparison="MEFFTO"
treshold=10
Isolate_four_categories<-function(comparison,condition1,condition2,treshold){
    dataframeofinterest<-get(comparison)
    Conditionsofinterest<-c(condition1,condition2)
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