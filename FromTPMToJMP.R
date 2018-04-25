##PACKAGES : gtable
#install.packages("gtable")
#install.packages("ggplot2")
#install.packages("grid")
#install.packages("gridExtra")

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
mergedTPM<-read.csv("MergedTPM.csv")


#########################################################################################
#1. PLOT AVG VALUES OF EACH CONDITION TO DEFINE "LOW" EXPRESSION FOR EACH CONDITION######
#########################################################################################


###### 1A ___Removing TPM and FPKMS from the colnames ###################################
colnames(mergedTPM)<-gsub("TPM","",colnames(mergedTPM))
colnames(mergedTPM)<-gsub("FPKMS","",colnames(mergedTPM))
colnames(mergedTPM)

####### 1B ___OUR CONDITIONS  _ UNDER CONSTRUCTION ######################################
#Find a way to Rep them from "everything before _"
#In the mean time, for now:
Conditions <- c("MEFFTOKno", "MEFFTOwt","MEFshnRNPABHet","MEFshnRNPABKno",
                "OlineoFTOCtrol", "OlineoFTOsiRNA","OlineohnRNPABCtrol","OlineohnRNPABCtrol")

####### 1C ___Adding average condition per column: ###################################### 

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

  
#colnames should have extra with thenameofourconditions_avg at the end 
colnames(mergedTPM)

######### 1D___Plotting: #########################################################

# Need to plot the histogram or density per condition(avged) with the y axis as a percentage

#We need to "reshape" and "melt" inorder to be able to plot all densities together

#Replacing 0 by NA > To melt it and address 0
mergedTPM_0NA <- mergedTPM
mergedTPM_0NA[mergedTPM_0NA ==0]<- NA
columnames <- colnames(mergedTPM_0NA)
avg_columnname <- columnames[grep("_avg",columnames)]
#New dataframe for ploting
mergedTPM_melted=melt(mergedTPM_0NA[,avg_columnname])
mergedTPM_melted=mergedTPM_melted[complete.cases(mergedTPM_melted),]


#plot in fonction of maximum y : 
maxofy= 0.05
  #Make First plot: Show the general shape of the data Violin. 
  p1<-ggplot(data=mergedTPM_melted,aes(x=variable,y=value)) +
              geom_violin()+
              ylim(0,maxofy)
  
  #Make Second plot: Show the histogram repartition
  p2<-ggplot(data=mergedTPM_melted,aes(x=value,color=variable))+
             geom_density()+
              xlim(0,maxofy) 
  # Display them together
  paste0("plot",maxofy)<-plot0.25<-grid.arrange(p1, p2, nrow = 2)
  
#From the plots it looks like let say 0.03 could be a good number
# need to be able to 
  #rank them and define what is the number which above of you have lets say 80% of the genes
  # BUT ONLY THE GENES THAT ARE NOT ZEROS!!! SO LIKE 80% OF THE NON ZERO GENES

  #define the quantile function
  f<-numcolwise(quantile,0.30,na.rm=TRUE)
  quantile2<- f(mergedTPM_0NA)
  #transpose
  quantile2.t<-t(quantile2)
  quantile2.df<-as.data.frame(quantile2.t)
  View(quantile2)
  #Plot the quantiles
  ggplot(quantile2.df,
         aes(x=rownames(quantile2.df),
             y=V1)) + 
        geom_point(color="blue")
  
##########1D Now you need to decide how you want to define low expression value 
      ###### and how you want to remove then
      ###### Should remove per condition? How do you define

#REMOVE ALL ZEROS
#1. REMOVE LINES WITH ALL ZEROS
#Adds a sum column (Starting at 3 because 1 and 2 = gene names=non numeric)
mergedTPM$Sum <-rowSums(mergedTPM[,3:ncol(mergedTPM)])
#Removing rows(genes) where ALL the values are 0
mergedTPM<-mergedTPM[mergedTPM$Sum!=0,]
mergedTPM$Sum <- NULL

#2. ADD +1 TO EVERYTHING
mergedTPM<-cbind(mergedTPM[,1:2],(mergedTPM[,3:ncol(mergedTPM)]+1))


library(ggplot2)
