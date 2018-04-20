#TO RUN IT : 
#1. CHANGE THE INPUT AND OUTPUT DIRECTORY
#
#2A. copy it all and press RUN
#OR
#2B. Copy it all and paste it into console


##THINGS TO CHANGE : 
#INPUT : ENTER WHERE YOUR FPKMS ARE
MyFPKMDirectory<-"C:/Users/nn503/Dropbox/Desktop/FPKMtoTPM/FPKMs/"
#OUTPUT : ENTER WHERE YOU WANT YOUR TPM TO BE 
MyTPMDirectory<-"C:/Users/nn503/Dropbox/Desktop/FPKMtoTPM/TPMs/"
dir.create(MyTPMDirectory)


#JUST COPY PASTE THE FOLLOWING


####THIS JUST HAS TO BE FIRST
####DEFINE THE FPKM TO TPM FUNCTION

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
  # SaveAsAFILE
  write.table(FPKMExtract,file=paste0(MyTPMDirectory,"TPM",TheFilename),
              row.names = FALSE,
              col.names =FALSE,
              sep="\t")
}
#ENF OF FUNCTION


#Look for the files in your FPKM Directory and remove the .txt (Change .txt if you have under another extension)
MyFPKMsNames<-list.files(path = MyFPKMDirectory)
#Runs the FPKMtoTPM on all the files present in your FPKM Directory
lapply(MyFPKMsNames,FPKM_to_TPM)

print("All Done")