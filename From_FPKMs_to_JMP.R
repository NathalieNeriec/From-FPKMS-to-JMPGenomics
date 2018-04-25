#####NOTE IN THIS ONE YOU NEED TO HAVE MERGED DATA >> MAKE A FUNCTION MERGEDATA FROM AYMAN S CODE
setwd("X:/LAB_MEMBERSstorageplace/NATHALIE_Analyses/CBFA_FTO_InOlineo_ForMagnus/PART1_TOJMP/")
#Load libraries:
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(reshape2)
library(scales)
library(dplyr)
library(plyr)
library(parallel)
library(xlsx)

############################################################################################
######### USER DEFINED #####################################################################
############################################################################################
#Your files need to have a name that contains your comparisons "Olineo" "OlineoFTO" etc
# Directory where the FPKMs ARE
MyFPKMDirectory<-"./1_FPKMS/"

Comparisons<-c("OlineoFTO","OlineohnRNPAB")

Conditions <- c("OlineoFTOCtrol", "OlineoFTOsiRNA",
                "OlineohnRNPABsiRNA","OlineohnRNPABCtrol")

############################################################################################
### 1. FROM FPKM TO TPM#####################################################################
############################################################################################

MyTPMDirectory<-"./1_TPMS/"
dir.create(MyTPMDirectory)

#FUNCTION FPKM_to_TPM
FPKM_to_TPM <- function(TheFilename){
  # TPM = FPKM / (Sum(FPKM)/10^6)
  #Getfrom the FPKM file   
  OrigDat<-read.table(paste0(MyFPKMDirectory, TheFilename),
                      header=TRUE)
  #EXTRACT FPKM info
  FPKMExtract <- data.frame(OrigDat$gene_id,OrigDat$FPKM)
  # CALCULATE SUM(FPKM)
  Total_FPKMS<- sum(FPKMExtract$OrigDat.FPKM)
  #Divide Total FPKM by 10^6 (so dont have to multiple everythin by 10^6)
  Total_FPKMS_divided<-Total_FPKMS/1000000
  # Add a new TPM column and process for each Row from FPKM to TPM
  FPKMExtract$TPMs<-FPKMExtract$OrigDat.FPKM/Total_FPKMS_divided
  #REMOVE THE RPKMS COLUMN
  FPKMExtract$OrigDat.FPKM<-NULL
  
  TheFilename<-gsub("FPKMS","",TheFilename)
  
  # SaveAsAFILE
  write.table(FPKMExtract,file=paste0(MyTPMDirectory,TheFilename),
              row.names = FALSE,
              col.names =FALSE,
              sep="\t")
}
#ENF OF FUNCTION

#Look for the files in your FPKM Directory and remove the .txt (Change .txt if you have under another extension)
MyFPKMsNames<-list.files(path = MyFPKMDirectory)
#Runs the FPKMtoTPM on all the files present in your FPKM Directory
lapply(MyFPKMsNames,FPKM_to_TPM)

print("The outputs are TPM txt files in the ./TPMs folder")
print("1. FROM FPKM TO TPM: All Done")


############################################################################################
### 2. MERGE THE DATA PER COMPARISON ########################################################
############################################################################################

MyTPMNames<-list.files(path = MyTPMDirectory)
#FUNCTION
MakeaMerge<-function(TheComparison){ 
     #Get the names of the files to merge (All the files that contain "TheCOmparison" in it)
     MyComparisonFilenames<-MyTPMNames[grep(TheComparison,MyTPMNames)]
    
     #create the final dataframe and already adding the sorted genes names
     MergedTPM <- read.csv(paste0(MyTPMDirectory,MyComparisonFilenames[1]),header=FALSE,sep='\t')
     MergedTPM <- as.data.frame(MergedTPM[order(MergedTPM[,1]),1])

     #Create the names of the final Merged File  
     listofnames<-as.vector(1:(length(MyComparisonFilenames)+1))
     listofnames[1]<-"gene_id"
  
     #loop through the comparison files to extract the values and append them to MergedTPM
     for (i in 1:length(MyComparisonFilenames)){
          #save the name of the condition
          listofnames[i+1]<-gsub(".txt","",MyComparisonFilenames[i])
          #Read the file 
          fileContent<-read.csv(paste0(MyTPMDirectory,MyComparisonFilenames[i]),header=FALSE,sep='\t')
          #Sort by gene_id incase they are not sorted
          fileContent = fileContent[order(fileContent[,1]),]
          MergedTPM<-cbind(as.data.frame(MergedTPM),as.data.frame(fileContent[,2]))
          }
      #Replace col names by sample name
      colnames(MergedTPM)<-listofnames  
  
      ###Add the gene symbol/names from the ENSEMBLE ID
          load("./Rdas/Mus_musculus.GRCm38.82.Rda")
          # Calculate the number of cores
          no_cores <- detectCores() - 1
          # Initiate cluster
          cl <- makeCluster(no_cores)
          #levelsList = character(length(ensNames))
          levelsList = parallel::parLapply(cl,MergedTPM$gene.ids, function(x){
                return(geneid2name[geneid2name$gene_id == as.character(x),]$gene_name)
                })
          stopCluster(cl)    
          #Add the newly formed gene2idname key to the mergedTPM  
          MergedTPM<-  merge(geneid2name[, c("gene_id", "gene_name")],MergedTPM , by="gene_id")
      ###Returns:
          return(MergedTPM)
          
}
#Loop throughout Our comparisons
for (i in Comparisons){
  LoopingMerge<-MakeaMerge(i)
  assign(paste0(i,"_Merged"),LoopingMerge)
}

print("2. MERGE THE DATA: All Done")
print("The outputs are dataframe Named TheComparison_Merged")


############################################################################################
### 3. Add+1 and Log2 trasnformation #######################################################
############################################################################################

#FUNCTION
Add1AndLog2Transform<- function(TheComparison){
  #Get the file to add 1 and log2transform
  TempMergedFile<-get(paste0(TheComparison,"_Merged"))
  #Add+1
  TempMergedFile[,!(names(TempMergedFile) %in% c("gene_id","gene_name"))] = TempMergedFile[,!(names(TempMergedFile) %in% c("gene_id","gene_name"))] + 1
  
  #Log2
  MergedFile = cbind(TempMergedFile[,1:2],log(TempMergedFile[,!(names(TempMergedFile) %in% c("gene_id","gene_name"))],2))
  return(MergedFile)}
#Loop throughout Our comparisons
for (i in Comparisons){
  LoopTransfo<-Add1AndLog2Transform(i)
  assign(paste0(i,"_MergedTransformed"),LoopTransfo)
}


print("3. Add+1 and Log2 trasnformation : ALL DONE")
print("The outputs are dataframe Named TheComparison_MergedTransformed")


############################################################################################
### 4. LOWEXPRESSION PLOTING  ################################################
############################################################################################

#THREE PARTS
    #1. Add average columns
    #2. Plots
Conditions <- c("OlineoFTOCtrol", "OlineoFTOsiRNA",
                "OlineohnRNPABsiRNA","OlineohnRNPABsiRNA")

#1. Add average columns
     #FUNCTION
     Add_Avg_Col_percondi <- function(mergedTPM,condition){
              #Store colnames
              columnames <- colnames(mergedTPM)
              #Find column containing control CTRL
              condition_columnname <- columnames[grep(condition,columnames)]
             
              if(length(condition_columnname)==0) {
              } else {
                     #print(paste0("condition",condition,"found in the dataset",TheComparison))
                     #Get only CTRL subsets of counts
                     condition_columns <- mergedTPM[,condition_columnname]
                     #calculaterow averages and insert into a new column
                     mergedTPM$condition_Avg=rowMeans(condition_columns)
                     #Change the name  of the newly added column to match the pattern
                     colnames(mergedTPM)[colnames(mergedTPM) =="condition_Avg"]<-paste0(condition,"_avg")
                     } 
              return(mergedTPM) 
              
            }

     #Looping through the comparisona nd through the conditions to add a column of average at the end      
     for (TheComparison in Comparisons){
              mergedTPM<-get(paste0(TheComparison,"_MergedTransformed"))
              for (condition in Conditions){
                  mergedTPM<-Add_Avg_Col_percondi(mergedTPM,condition)
                  } 
              assign(paste0(TheComparison,"_Avg"),mergedTPM)
     }
     
     print("_Avg dataframe have been created")
     
#2_ Generate the plots
      #yvalues can bring that one outside if we want
      yvalues<-c(10,50,100,1000)                

      #FUNCTION 
      lowexpressionplot<- function(mergedTPM){
            #Create an directory to Store the Plots
                dir.create(file.path("./2_LowExpressionPlots"), showWarnings = FALSE)
                out_dir_Plotdata<-"./2_LowExpressionPlots/"
            
            #Generating a new average df only df with NA instead of zero
                #get only the columns that contain _avg
                columnames<-colnames(mergedTPM)
                avg_columnname <- columnames[grep("_avg",colnames(mergedTPM))]
                #Generate a df with only the avg
                mergedTPM_avg <- mergedTPM[,avg_columnname]
                #Replacing 0 by NA (for ignoring during plotting)
                mergedTPM_avg[mergedTPM_avg==0]<-NA
            
            
            # Plotting every column, into a new file with y valies of "yvalues"
                for(i in 1:ncol(mergedTPM_avg)){
                    #ranking the column
                    namey<-sort(mergedTPM_avg[,i],na.last = TRUE)
                    #saving the name for the plot
                    plottitle<-avg_columnname[[i]]

                    for (maxofy in yvalues){
                         plotname<-paste0("plot_",plottitle,"_",maxofy)
                         png(paste0(out_dir_Plotdata,plotname,".png"))
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
}
     
       #Looping to generate the plots
      for (TheComparison in Comparisons){
            mergedTPM<-get(paste0(TheComparison,"_Avg"))
            lowexpressionplot(mergedTPM)
      }
      
      
      
      print("Ignore the warnings.")
      print("4. LOWEXPRESSION PLOTING: All Done")
      print("The graphs are available in ./2_LowExpressionPlots")
      
      print("Do you want to keep a treshold of 1? Reminder: log2(1+1)=1")
       
############################################################################################
### 5. FILTERING   ###########################################################################
############################################################################################      
      # fun <- function(){
      #   treshold <- readline("What value do you want to use for treshold?Reminder: log2(1+1)=1")  
      #   return(treshold)
      # }
      # if(interactive()) treshold<-fun()
      # print(paste0("treshold has been set to ",treshold))
      treshold =1


#FUNCTION : This function takes a comparison dataframe, look for what conditions amond the "Conditions"
      #then figure out which genes dont have AT LEAST one condition with over treshold and remove them
      #Then it figure out which genes have AT LEAST one condition that is ALL zero(log2FC(0+1)=0) and extract them
      #Then it returnsA LIST THAT CONTAINS TWO DATAFRAME
          # List[1] the ON OFF genes 
          # List[2] the Filtered expressed in one condition genes 
FilteringLowExpressAndONOFF<-function(TheComparison,Conditions,treshold){
        MergedTPM<-get(paste0(TheComparison,"_MergedTransformed"))
        NcolumnOriginal<-ncol(MergedTPM)
        columnames<-colnames(MergedTPM)

        #Find the list of SelectedConditions existing within the file that are withing our columnames
        test<-as.data.frame(sapply(Conditions,function (y) sapply (columnames, function (x) grepl(y,x))))        
        test<-rbind(test,colSums(test))
        test<-test[nrow(test),]
        SelectedConditions<-colnames(test[,test[1,]!=0])
        
        #Create a column per condition that sum the percentage of condition being above treshold  
        for (condition in SelectedConditions) {
          #Get the columname matchig the condition
          condition_columnname <- columnames[grep(condition,columnames)]
          #Add a newcolumn with the percentage in each condition above your treshold here 5
          MergedTPM$newcolumn<-
            apply(as.data.frame(MergedTPM[, condition_columnname]), 
                  MARGIN = 1, 
                  function(x) (sum(x>treshold)/length(condition_columnname)))
          colnames(MergedTPM)[colnames(MergedTPM) =="newcolumn"]<-paste0("testabove",treshold,"_",condition)
        }
        
        
        
        
        #Make a new data frame with the Filter  
        columnames2 <- colnames(MergedTPM)
        testabove_columnname<-columnames2[grep("testabove",columnames2)]
        
        #Only keep the rows where ANY have at least 50% above the treshold (remove categorya)
        MergedTPM<-
          MergedTPM[apply(MergedTPM[,testabove_columnname],MARGIN=1,function(x)any(x>0.50)),]
        Ngeneswithatleast1<-nrow(MergedTPM)
        
        #Add a column per condition with the average
        for(condition in SelectedConditions){
        MergedTPM<-Add_Avg_Col_percondi(MergedTPM,condition)
        }
        
        #Getting names of avg columns
        columnames2 <- colnames(MergedTPM)
        avg_columnname<-columnames2[grep("avg",columnames2)]
        
        #ISOLATING THE ON OFF GENES  
        #Only keep the where ANY of the avg is 0 into a new data frame
        MergedTPM_ONOFF<-
          MergedTPM[apply(
            MergedTPM[,avg_columnname],
            MARGIN=1,
            function(x)any(x==0)),
            1:NcolumnOriginal]
        NgenesONOFF<-nrow(MergedTPM_ONOFF)
        
        #Remove  the onoff from the expressed genes (DE+ testable)
        MergedTPM_Expressed<-
          MergedTPM[apply(
            MergedTPM[,avg_columnname],
            MARGIN=1,
            function(x)!(any(x==0))),1:NcolumnOriginal]
        NgenesExpressed<-nrow(MergedTPM_Expressed)
        
        #RETURN
        #Make a list of the two dataframe of interest
        my.list<-list("MergedTPM_ONOFF"=MergedTPM_ONOFF,"MergedTPM_Expressed"=MergedTPM_Expressed)
        return(my.list)
}


#Looping through the Comparisons we have 
dir.create("./3_ONOFF/")
dir.create("./3_Expressed/")

for(TheComparison in Comparisons){
  #Run the Filtering
  Templist<-FilteringLowExpressAndONOFF(TheComparison,Conditions,treshold)
  
  #ONOFF Save as a dataframe and an excel file
  MergedTPM_ONOFF<-as.data.frame(Templist[1])
  colnames(MergedTPM_ONOFF)<-gsub("MergedTPM_ONOFF."," ",colnames(MergedTPM_ONOFF))
  write.xlsx(MergedTPM_ONOFF,paste0("./3_ONOFF/",TheComparison,"_ONOFF.xlsx"),row.names = FALSE)
  assign(paste0(TheComparison,"_ONOFF"),MergedTPM_ONOFF)

  #Expressed Save as a dataframe and an excel file
  MergedTPM_Expressed<-as.data.frame(Templist[1])
  colnames(MergedTPM_Expressed)<-gsub("MergedTPM_Expressed."," ",colnames(MergedTPM_Expressed))
  write.xlsx(MergedTPM_Expressed,paste0("./3_Expressed/",TheComparison,"_Expressed.xlsx"),row.names = FALSE)
  assign(paste0(TheComparison,"_Expressed"),MergedTPM_ONOFF)
  
}


print("5. FILTERING: All Done")
print("Outputs are two dataframes per comparison : TheComparison_ONOFF and TheComparison_Expressed")

############################################################################################
### 6. ExperimentalDesignFile   ############################################################
############################################################################################   


So we have TheComparison_Expressed 
We want a file that has one column label
1	TPMOlineoFTOCtrolFPKMS_1	control
2	TPMOlineoFTOCtrolFPKMS_2	control
3	TPMOlineoFTOCtrolFPKMS_3	control
4	TPMOlineoFTOCtrolFPKMS_4	control
5	TPMOlineoFTOsiRNAFPKMS_1	FTOshRNA
6	TPMOlineoFTOsiRNAFPKMS_2	FTOshRNA
7	TPMOlineoFTOsiRNAFPKMS_3	FTOshRNA
8	TPMOlineoFTOsiRNAFPKMS_4	FTOshRNA