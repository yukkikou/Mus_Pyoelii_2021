#setwd("D:/_Plasmodium_yoelii/")
#after yoelii.R

library(WGCNA)
library(reshape2)
library(stringr)
library(data.table)
library(RColorBrewer)
display.brewer.all()
cols <- brewer.pal(6, "Pastel2")
cols2 <- rev(brewer.pal(11, "RdYlBu"))
options(stringsAsFactors = FALSE)

# gene count matrix and filtering
yoeExp = expGene # filtering
musExp = readCount %>% 
  filter_all(all_vars(.>50))

datExpr = rbind(musExp, yoeExp)

# PCA: no sample outlier
expr0 <- t(as.matrix(datExpr0))
expr.pca0 <- PCA(expr0, ncp = 2, scale.unit = TRUE, graph = FALSE)
plot(expr.pca0)


# normalization
rpk <- datExpr0 / (nrow(datExpr0) / 1000)
datExpr0 <- t(t(rpk)/colSums(rpk) * 1000000)
head(datExpr0)

# abstruct exp matrix
datExpr0 <- t(datExpr0)

# check missing and outlier
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
#If the last statement returns TRUE, all genes have passed the cuts.
#If not, remove the offending genes and samples from the data as following
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
dim(datExpr0)

#cluster the samples to see if there are any obvious outliers
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
#show F2-221 is one outlier
#Choose a height cut that will remove the offending sample
# Plot a line to show the cut
abline(h = 30000, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 30000, minSize = 10)
#save picture!
table(clust)
# clust 1 contains the samples we want to keep
keepSamples = (clust==0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#The variable datExpr now contains the expression data ready for network analysis  

########################################################################step 2: phonetype data
#Loading clinical trait data
allTraits <- sampleGroup
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples = rownames(datExpr)
traitRows = match(femaleSamples, allTraits$Sample)
datTraits = allTraits[traitRows, ]
rownames(datTraits) = allTraits[traitRows, 1]
collectGarbage()
#expression data in the variable datExpr
#corresponding clinical traits in the variable datTraits

########################################################################step 3: break and rest a moment :)
##visualize how the clinical traits relate to the sample dendrogram
#Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low value, red means high, grey means missing entry

traitColors = numbers2colors(as.numeric(datTraits$Group), signed = FALSE)

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = "figure/WGCNA/FemaleLiver-01-dataInput.RData")


#2.a Automatic network construction and module detection
########################################################################step 1: Choosing the soft-thresholding power
lnames = load(file = "figure/WGCNA/FemaleLiver-01-dataInput.RData")
lnames
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.7420,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



#2.b Step-by-step network construction and module detection
########################################################################step 1: Choosing the soft-thresholding power
########################################################################step 2: Co-expression similarity and adjacency
softPower = 7
adjacency = adjacency(datExpr, power = softPower)
#Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
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

########################################################################step 3: Merging of modules similar
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

#as figure shows
MEDissThres = 0.05
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
sizeGrWindow(12, 9)
#pdf(file = "geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "figure/WGCNA/FemaleLiver-02-networkConstruction-stepByStep.RData")

#3. Relating modules to external information and identifying important
########################################################################step 1: Quantifying module–trait associations
# Load the expression and trait data saved in the first part
lnames = load(file = "figure/WGCNA/FemaleLiver-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "figure/WGCNA/FemaleLiver-02-networkConstruction-stepByStep.RData")
lnames
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
datTraits = rep(c(0,1),c(3,5))
#datTraits <- datTraits[,2]

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(15,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = "condition",
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#shows several significant module–trait associations

highCorModule = moduleTraitCor %>%
  as.data.frame() %>%
  rownames_to_column(var = "module") %>%
  setNames(c("module", "cor")) %>% 
  filter(abs(cor) > 0.8)

labeledHeatmap(Matrix = as.data.frame(highCorModule$cor),
               xLabels = "condition",
               yLabels = highCorModule$module,
               ySymbols = highCorModule$module,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = round(highCorModule$cor,3),
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationshis (> 0.8)"))

########################################################################step 2:  Gene Significance and Module Membership
# Define variable weight containing the weight column of datTrait
condition = datTraits
names(condition) = "condition"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
geneTraitSignificance = as.data.frame(cor(datExpr, condition, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", "condition", sep="")
names(GSPvalue) = paste("p.GS.", "condition", sep="")

#take module brown as example, identifying genes with high GS and MM
highCorModule %>% filter(abs(cor) > 0.9) %>% arrange(-abs(cor))

module = "brown"
column = match(module, modNames)
moduleGenes = moduleColors==module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))

#floralwhite too thin: col = "navy"
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
#GS and MM are highly correlated
#illustrating that genes highly significantly associated with a trait are often 
#also the most important (central) elements of modules associated with the trait.

#Network visualization using WGCNA functions
lnames = load(file = "figure/WGCNA/FemaleLiver-02-networkConstruction-stepByStep.RData")
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 20)
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function
sizeGrWindow(9,9)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

nSelect = 400
# For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
#TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

#Visualizing the network of eigengenes
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, condition))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)


###********** Special Genesets **********###
# Exporting to Cytoscape
probes = colnames(datExpr)
hiCorYoeliiModule = highCorModule %>% 
  filter(abs(cor) > 0.9) %>% 
  arrange(-abs(cor)) %>%
  separate(module, c("prefix", "color"), "E") %>%
  select(color, cor) %>%
  filter(color %in% unique(moduleColors[probes %in% rownames(yoeExp)]))

hiCorYoeliiModule$color
#[1] "turquoise"    "midnightblue" "lightyellow"  "magenta"      "yellow" 

# Select module
modules = "midnightblue" 
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
modGenes = probes[inModule]

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("figure/WGCNA/CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("figure/WGCNA/CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.3,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])

#Top edge and node
nTop = 100
IMConn = softConnectivity(datExpr[, modProbes])
top = (rank(-IMConn) <= nTop)
topmod <-modTOM[top, top]
hiYoeGene = rownames(topmod)[rownames(topmod) %in% row.names(yoeExp)]
write.table(hiYoeGene, paste("figure/WGCNA/HighYoeliiGene-", paste(modules, collapse="-"), ".txt", sep=""), quote = F, col.names = F, row.names = F)

