# source("http://bioconductor.org/biocLite.R") 
# biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
# install.packages("WGCNA")

library(impute)
library(dynamicTreeCut)
library(qvalue)
library(flashClust)
library(Hmisc)
library(WGCNA)


####GGETTING THE DATA READY###
#The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Take a quick look at what is in the data set:
FilteredMergedTPM<-read.csv("FilteredMergedTPM.csv")
FilteredMergedTPM<-FilteredMergedTPM[,1:33]

#This not only only extract the value but TRANSPOSE 
  #The column ar ehte genes
datExpr0 = as.data.frame(t(FilteredMergedTPM[, -c(1:2)]))
names(datExpr0) = FilteredMergedTPM$gene.ids;
rownames(datExpr0) = names(FilteredMergedTPM)[-c(1:2)];

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


####CLUSTERING OUR DATA AND SPOT ANYOUTLIERS OR SO THAT WE WANT TO REMOVE######
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)
#IF OUTLIERS
# # Plot a line to show the cut
# abline(h = 10, col = "red");
# # Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 200000, minSize = 10)
# table(clust)
# # clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#####ADDING CLINICAL OR WHATEVER DESIGN INFO AND PLOT IT UNDER THE DENDOGRAM OF OUR SAMPLEs BASE DON COLORS #######
##NEEDS TO BE NUMERICAL##

traitData  <-read.csv("EXPDesignstuff.csv")
allTraits<-traitData
cellSamples=rownames(datExpr0)
traitRows<-match(cellSamples,allTraits[,1])
datTraits<-allTraits[traitRows,-1]
rownames(datTraits) = allTraits[traitRows, 1]
collectGarbage()

  # Re-cluster samples
  sampleTree2 = hclust(dist(datExpr), method = "average")
  # Convert traits to a color representation: white means low, red means high, grey means missing entry
  traitColors = numbers2colors(datTraits, signed = FALSE);
  # Plot the sample dendrogram and the colors underneath.
  plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

######CHOSING SOFT TRESHOLD POWER #####
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  #A little time
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


  # Plot the results: Choose the one number the highest as softpower
      sizeGrWindow(9, 5)
      par(mfrow = c(1,2));
      cex1 = 0.9;
      # Scale-free topology fit index as a function of the soft-thresholding power
      plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
           ylim=c(-1,1),
           main = paste("Scale independence"));
      text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           labels=powers,cex=cex1,col="red");
      # this line corresponds to using an R^2 cut-off of h
      abline(h=0.90,col="blue")
      # Mean connectivity as a function of the soft-thresholding power
      plot(sft$fitIndices[,1], sft$fitIndices[,5],
           xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
           main = paste("Mean connectivity"))
      text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 6;
adjacency = adjacency(datExpr, power = softPower);


##########Turn adjacency into topological overlap######
#TAKESFOREVER!!!!!!!!!!!!1
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

###########PlotTREE and MAkes modules and MErgeModules#######

#####PlotTree
    # Call the hierarchical clustering function
      geneTree = hclust(as.dist(dissTOM), method = "average")
    # Plot the resulting clustering tree (dendrogram)
      sizeGrWindow(12,9)
      plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
             labels = FALSE, hang = 0.04);

#####MakeMOdules
    # We like large modules, so we set the minimum module size relatively high:
      minModuleSize = 50;
    # Module identification using dynamic tree cut:
      dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);

      table(dynamicMods)
    # Convert numeric lables into colors
      dynamicColors = labels2colors(dynamicMods)
      table(dynamicColors)
    # Plot the dendrogram and colors underneath
      sizeGrWindow(8,6)
      plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                         dendroLabels = FALSE, hang = 0.03,
                         addGuide = TRUE, guideHang = 0.05,
                         main = "Gene dendrogram and module colors")
      
#####Merge Similar Modules
      # Calculate eigengenes
          #I ADDED EXCLUDEGREY=TRUE
          MEList = moduleEigengenes(datExpr,excludeGrey = TRUE, colors = dynamicColors)
          MEs = MEList$eigengenes
      # Calculate dissimilarity of module eigengenes
          MEDiss = 1-cor(MEs);
      # Cluster module eigengenes
          METree = hclust(as.dist(MEDiss), method = "average");
      # Plot the result
          sizeGrWindow(7, 6)
          plot(METree, main = "Clustering of module eigengenes",
               xlab = "", sub = "")
  
              #DEFINE based on how you want to merge modules
        abline(h = 0.42, col = "red");
        MEDissThres = 0.42
      
    # Plot the cut line into the dendrogram
        abline(h=MEDissThres, col = "blue")
    # Call an automatic merging function
    merge = mergeCloseModules(datExpr,
                              dynamicColors,
                              cutHeight = MEDissThres, 
                              verbose = 3)

    # The merged module colors
    mergedColors = merge$colors;
    # Eigengenes of the new merged modules:
    mergedMEs = merge$newMEs;
    sizeGrWindow(12, 9)
    pdf(file = "geneDendro-3.pdf", wi = 9, he = 6)
    plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                        c("Dynamic Tree Cut", "Merged dynamic"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05)
    # Rename to moduleColors
    moduleColors = mergedColors
    # Construct numerical labels corresponding to the colors
    colorOrder = c("grey", standardColors(50));
    moduleLabels = match(moduleColors, colorOrder)-1;
    MEs = mergedMEs
    
    
    
    ###################PART 3############
    
    # Define numbers of genes and samples
    nGenes = ncol(datExpr);
    nSamples = nrow(datExpr);
    # Recalculate MEs with color labels
    MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)
    moduleTraitCor = cor(MEs, datTraits, use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    
    sizeGrWindow(10,6)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3));
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = names(datTraits),
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = greenWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    
 ###Gene relationship to trait and important modules: Gene Significance and Module Membership
    
    # Define variable variabletest containing the variabletest column of datTrait
    variabletest = as.data.frame(datTraits$mutantsonly);
    names(variabletest) = "mutantsonly"
    # names (colors) of the modules
    modNames = substring(names(MEs), 3)
    geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
    MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
    
    names(geneModuleMembership) = paste("MM", modNames, sep="");
    names(MMPvalue) = paste("p.MM", modNames, sep="");
    geneTraitSignificance = as.data.frame(cor(datExpr, variabletest, use = "p"));
    GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
    names(geneTraitSignificance) = paste("GS.", names(variabletest), sep="");
    names(GSPvalue) = paste("p.GS.", names(variabletest), sep="");
    
    

    ##Intramodular analysis: identifying genes with high GS and MM
    for (module in ColorsME){
    column = match(module, modNames);
    moduleGenes = moduleColors==module;
    sizeGrWindow(7, 7);
    par(mfrow = c(1,1));
    png(file=paste0(module,"_toFtoKnoonly.png"))
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for both mutants",
                       main = paste("Module membership vs. gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2)
    dev.off()
    }
# modulecolor = "blue"
#     for (modulecolor in ColorsME){
#       tempdataframe<-names(datExpr)[moduleColors==modulecolor]
#       write.csv(tempdataframe,paste0("Genesinmodule",modulecolor,".csv"),row.names = FALSE)
#       assign(paste0("Genesinmodule",modulecolor),tempdataframe)
#       
#     }
    names(datExpr)[moduleColors=="floralwhite"]
    
# annot = read.csv(file = "GeneAnnotation.csv");
# probes2annot = match(probes, annot$substanceBXH)
# We now create a data frame holding the following information for all probes: 
#   probe ID, gene symbol, Locus Link ID (Entrez code), module color, 
# gene significance for variabletest, 
colnames(MEs)

"MEorangered4"     "MEblue"           "MEfloralwhite"    "MEdarkorange2"    "MEbrown4"         "MEskyblue3"       "MEyellowgreen"   
[8] "MEbisque4"        "MEdarkolivegreen" "MEsienna3"        "MEmidnightblue"   "MEorange"         "MEivory"          "MElightgreen"    
[15] "MEdarkturquoise"  "MEdarkgreen
