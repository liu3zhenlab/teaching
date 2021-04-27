### INSTALL WGCNA ###
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
## For R version < 3.5.0 ##
source("http://bioconductor.org/biocLite.R") 
biocLite(c("GO.db", "preprocessCore", "impute"))
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg"); 
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6)); 
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep=""); 
biocLite(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
install.packages("/Users/chenghe/Downloads/WGCNA_1.66.tgz", repos = NULL, lib=.Library, type = "source",dependencies = TRUE)
install.packages("robust")

## For R version > 3.5.0 ##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("GO.db", "preprocessCore", "impute"))
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg"); 
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6)); 
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep=""); 
BiocManager::install(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
install.packages("/Users/chenghe/Downloads/WGCNA_1.66.tgz", repos = NULL, lib=.Library, type = "source",dependencies = TRUE)
install.packages("robust")


### LOAD WGCNA ###
library(WGCNA)
allowWGCNAThreads(nThreads = 12)

### READ & FORMAT EXP INPUT ###
exp = read.table("GSE48213_gene_expression_fpkm.txt", stringsAsFactors = F,header = FALSE)
datexp = as.data.frame(apply(exp[-1,-1],1,as.numeric))
dim(datexp)
names(datexp) = exp$V1[-1]
rownames(datexp) = as.character(exp[1,-1])
datexp[1:10,1:10]

### READ TRAIT INPUT ###
trait = read.table("GSE48213_trait.txt",sep="\t",header=TRUE)
line = rownames(datexp)
traitrows = match(line,trait$gsm)
dattrait = as.data.frame(trait[traitrows,-1])
rownames(dattrait) = trait[traitrows,1]
colnames(dattrait)[1] <- "trait"
collectGarbage()

### SAMPLE CLUSTERING TO IDENTIFY OUTLIERS ###
sampleTree = hclust(dist(datexp), method = "average")
traitColors = numbers2colors(as.numeric(dattrait$trait), signed = FALSE)
sizeGrWindow(12,9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
pdf("sample_dendrogram.pdf",height = 6,width = 10)
plotDendroAndColors(sampleTree, traitColors,groupLabels = names(dattrait),main = "Sample dendrogram and trait heatmap",sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)  
dev.off()

### DEFINE THE SOFT THRESHOLDING POWER (β) ###
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sftdatexp = pickSoftThreshold(datexp, powerVector = powers, verbose = 3)
sizeGrWindow(9, 5)
pdf("define_soft_thresholding_power.pdf",height = 6,width = 10)
par(mfrow = c(1,2))
cex1 = 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sftdatexp$fitIndices[,1], -sign(sftdatexp$fitIndices[,3])*sftdatexp$fitIndices[,2], 
     xlab = "soft threshold (power)", ylab = "scale free topology model fit, signed R^2", type = "n",
     main = paste("scale independence"))
text(sftdatexp$fitIndices[,1], -sign(sftdatexp$fitIndices[,3])*sftdatexp$fitIndices[,2],
     labels=powers, cex=cex1, col = "red")

# Mean connectivity as a function of the soft-thresholding power
plot(sftdatexp$fitIndices[,1], sftdatexp$fitIndices[,5],
     xlab = "soft threshold (power)", ylab = "mean connectivity", type = "n",
     main = paste("Mean Connectivity"))
text(sftdatexp$fitIndices[,1], sftdatexp$fitIndices[,5], labels=powers, cex=cex1, col = "red") ## select the first stable power
dev.off()

### CONSTRUCT COEXP MODULE ###
net_module = blockwiseModules(datexp, power = 4,
                              TOMType = "unsigned", minModuleSize = 100, maxBlockSize = 6000,
                              reassignThreshold = 0, mergeCutHeight = 0.25,
                              numericLabels = TRUE, pamRespectsDendro = FALSE,
                              saveTOMs = FALSE, verbose = 5)
mergedcolors = labels2colors(net_module$colors)
pdf("WGCNA_module_P4M200.pdf",height = 6,width = 8)
plotDendroAndColors(net_module$dendrograms[[1]], mergedcolors[net_module$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

### SAVE MODULE RESULTS ###
moduleLabels = net_module$colors
moduleColors = labels2colors(net_module$colors)
MEs = net_module$MEs
geneTree = net_module$dendrograms[[1]]
t <- table(moduleColors)
out_module = cbind(exp[-1,1],mergedcolors,moduleLabels)
colnames(out_module) = c("id","module","ME")
write.table(out_module,file="Coexp_modules_for_each_gene_P4M100.txt",sep="\t",quote=FALSE,row.names = FALSE)
write.table(t,file="Summary_of_coexp_modules_P4M100.txt",sep="\t",quote=FALSE,row.names = FALSE)
write.table(MEs,file="Eigengene_for_each_module_P4M100.txt",sep="\t",quote=FALSE)

### Quantifying module–trait associations ###
# Define numbers of genes and samples
nGenes = ncol(datexp)
nSamples = nrow(datexp)
design=model.matrix(~0+ dattrait$trait)
colnames(design)=levels(dattrait$trait)
moduleColors <- labels2colors(net_module$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datexp, moduleColors)$eigengenes
MEs = orderMEs(MEs0);
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
pdf("module_trait_associations_P4M100.pdf",height = 9,width = 8)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               cex.lab.y = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait associations"))
dev.off()

### Exporting Files for Cytoscape Visualization ###
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datexp, power = 4)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("Visulaxation_P4M100_Cyt-edges_p0.2.txt", sep=""),
                               nodeFile = paste("Visulaxation_P4M100_Cyt-nodes_p0.2.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = colnames(datexp),
                               nodeAttr = moduleColors)
