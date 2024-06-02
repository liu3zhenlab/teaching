# DESeq2 R analysis

# data link
data_url="https://raw.githubusercontent.com/liu3zhenlab/teaching/master/RNA-Seq-Workshop/2022/"

# install packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages("matrixStats", repos="http://cran.us.r-project.org")
if (!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2") # DESeq2
if (!require("goseq", quietly = TRUE)) BiocManager::install("goseq") # GOSeq
if (!require("GO.db", quietly = TRUE)) BiocManager::install("GO.db", force=T) # GO.db

# preload modules
# panel.cor2; rnaseq.pca; normalization modules
pls=paste0(data_url, "/utils/load.R")
source(pls)

# load count data
rc <- paste0(data_url, "/data/rc.txt")
grc <- read.delim(rc)
nrow(grc)  # the number of rows/lines
head(grc, 1)

# read normalization (RPKM)
datacols <- c("ck_rep1", "ck_rep2", "ck_rep3",
              "trt_rep1", "trt_rep2", "trt_rep3")
lib.sizes <- colSums(grc[, datacols])
exon.sizes <- grc$ExonSize
grcn <- normalization(counts = grc,
                      normcolname = datacols,
                      libsize = lib.sizes,
                      exonsize = exon.sizes,
                      methods = "RPKM")
head(grcn)

# create a matrix containing count data
geneid <- grc$Gene
in.data <- as.matrix(grc[, 3:8])
head(in.data, 1)

# sample names and grouping information (treatment)
sample.ids <- colnames(in.data)
treatment <- c("ck", "ck", "ck", "trt", "trt", "trt")
sample.info <- data.frame(row.names=sample.ids, trt=treatment)
sample.info

# organize data and design
dds <- DESeqDataSetFromMatrix(countData=in.data,
                              colData=sample.info,
                              formula(~trt))

# DE test
dds <- DESeq(dds)
# Wald test: use the estimated se of a LFC to test if it is equal to zero

# re-organize DE result
res <- results(object = dds)
res <- data.frame(res)
res$Gene <- geneid
res <- res[,c("Gene","baseMean","log2FoldChange","pvalue","padj")]
nrow(res)

# Merge the normalized result with the DE result
out <- merge(grcn, res, by = "Gene")
out <- data.frame(out)
head(out, 2)

# count numbers of significant genes at different FDRs
sum(!is.na(out$padj) & out$padj < 0.05)

# extract data of significant DEG
sig <- out[!is.na(out$padj) & out$padj < 0.05, ]

# plot p-value histogram
pvals <- out$pvalue
hist(pvals,main="Histogram",xlab="p-values",ylab="Number of genes")

# plot scatter plot
x <- log(out$ck_rep1.RPKM)
y <- log(out$ck_rep2.RPKM)
drange <- range(x, y, finite=T)
plot(x, y, main = "log of RPKM",
     xlab = "control rep1", ylab = "control rep2",
     xlim = drange, ylim = drange)
abline(a = 0, b = 1, col = "orange")

# pair-wise scatter plot
logrpkm <- log(out[, 9:14])
pairs(logrpkm, lower.panel=panel.smooth, upper.panel=panel.cor2)

# plot PCA
rnaseq.pca(out, norm.feature="RPKM", group.feature=c("ck", "trt"),
           cex=3, shape.code=c(5, 16), mean.cutoff=0.1, scaling = T,
           title="PCA-all genes", colors=c("dark green","orange"))

# MA plot
aval <- log(out$baseMean)
mval <- out$log2FoldChange
sig.set <- which(out$padj < 0.05)
plot(aval, mval, main="MA plot",
     xlab="log(mean of exp)", ylab="log2FC",
     ylim=c(-4, 6), pch=19, cex=0.3, col="blue")
# add points and a horizontal line
points(aval[sig.set], mval[sig.set], pch=19,cex=0.3,col="red")
abline(h = 0, cex = 1, col = "dark green")

# Volcano plot
log2val <- out$log2FoldChange  # log2 fold change
mlogP <- -log10(out$pvalue) # transformed p-values
### plot
plot(log2val, mlogP, main="Volcano plot",
     xlab="log2FC", ylab="-log10(p-value)",
     pch=19, cex=0.4,col="grey", xlim=c(-6,6), ylim=c(0,16))
### highlight
up <- which(out$padj<0.05 & log2val > 0) # up
down <- which(out$padj<0.05 & log2val < 0)  # down
points(log2val[up], mlogP[up], pch=19,cex=0.4,col="brown")
points(log2val[down], mlogP[down],pch=19,cex=0.4,col="blue")

# GOSeq
gdbf=paste0(data_url, "/data/go.txt")
godb <- read.delim(gdbf)
geneid <- as.character(out$Gene)  # gene vector
# a vector to indicate if the gene is DE (0 or 1)
de.vector <- as.integer(!is.na(out$padj) & out$padj < 0.05)
names(de.vector) <- geneid
countbias <- out$baseMean # total raw reads per gene
# bias fitting
pwf.counts <- nullp(DEgenes=de.vector, bias.data=countbias)

go <- goseq(pwf=pwf.counts, gene2cat=godb, method="Sampling",
            repcnt = 1000, use_genes_without_cat = F)  # GO enrichment

# check the description of a GO
example.go <- GOTERM[["GO:0004175"]]  # GO information
Definition(example.go)  # return GO definition
