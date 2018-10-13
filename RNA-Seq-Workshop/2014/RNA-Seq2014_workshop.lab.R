################################################################################
### RNA-Seq Workshop
### Sanzhen Liu
### 5/29/2014
################################################################################

### setup the working path
setwd(".")

### load function
source(file="http://129.130.89.83/RNA-Seq_Workshop_2014/RNA-Seq2014_workshop.preload.R")

### data set
d <- read.delim("http://129.130.89.83/RNA-Seq_Workshop_2014/data/example.txt")
ncol(d); head(d,2)

################################################################################
### Part 1: Histogram
################################################################################
### p-value histogram
hist(d$pval, main="p-value histogram",
     xlab="p-values", ylab="Number of genes")

### Since we are working on the p-value distribution, let us examine how it look
### if the p-values are derived from the null hypothesis.
### ttpvals was generated using the preloaded code
hist(ttpvals, main="Null Hypothesis",
     xlab="p-values", ylab="Number of genes")

### p.adjust is an R function to convert p-values to q-values to account for
### multiple test correction.
### convert p-value to q-value
qvals <- p.adjust(d$pval, method="BH")

### extract the significant gene subset
sig <- d[d$qval<0.05, ]
sig <- d[d$qval<0.05 & abs(d$log2FC)>1, ]
nrow(sig)

###
### scatter plot 1: raw count data
###
xdata.raw <- d$ck_rep1.RPKM
ydata.raw <- d$ck_rep2.RPKM
drange <- range(xdata.raw, ydata.raw, finite=T)
### plot a scatter plot between two samples
plot(xdata.raw, ydata.raw, xlab="control rep1", ylab="control rep2",
     xlim=drange, ylim=drange, main="Raw counts scatter plot")
abline(a=0, b=1, col="orange")

###
### scatter plot 2: normalized data
###
xdata <- log2(d$ck_rep1.RPKM)
ydata <- log2(d$ck_rep2.RPKM)
drange <- range(xdata, ydata, finite=T)
### plot a scatter plot between two samples
plot(xdata, ydata, xlab="control rep1", ylab="control rep2",
     xlim=drange, ylim=drange, main="RPKM scatter plot")
abline(a=0, b=1, col="orange")

###
### multiple scatter plots
###
pairs(log2(d[, 8:13]), lower.panel = panel.smooth, upper.panel = panel.cor2)

###
### PCA plot (principal component analysis)
###
### 1. full gene set, standardized 
rnaseq.pca(d, norm.feature="RPKM", group.feature=c("ck", "trt"), cex=3,
           shape.code=c(5, 16), mean.cutoff=0.1, title="PCA - full gene set",
           colors=c("dark green", "orange"), scaling=T)

### 2. full gene set, unstandardized
rnaseq.pca(d, norm.feature="RPKM", group.feature=c("ck", "trt"), cex=3,
           shape.code=c(5, 16), mean.cutoff=0.1,
           title="PCA - full gene set, no standarization",
           colors=c("dark green", "orange"), scaling=F)

### 3. significant gene set, standardized
rnaseq.pca(sig, norm.feature="RPKM", group.feature=c("ck", "trt"), cex=3,
           shape.code=c(5, 16), mean.cutoff=0.1, title="PCA - sig gene set",
           colors=c("dark green", "orange"), scaling=T)

### 4. significant gene set, unstandardized
rnaseq.pca(sig, norm.feature="RPKM", group.feature=c("ck", "trt"), cex=3,
           shape.code=c(5, 16), mean.cutoff=0.1,
           title="PCA - sig gene set, , no standarization",
           colors=c("dark green", "orange"), scaling=F)

###
### volcano plot
###
log2val <- d$log2FC  # log2 fold change
minusLogPvals <- -log10(d$pval) # transformed p-values
plot(log2val, minusLogPvals, xlab="log2FC", ylab="-log10(p-value)",
    pch=19, cex=0.4, col="grey",xlim=c(-6,6),ylim=c(0,16),
    main="Volcano plot") ### plot
sig.up <- which(d[,"qval"]<0.05 & log2val > 0) # significant UP set
sig.down <- which(d[,"qval"]<0.05 & log2val < 0)  # significnat DOWN set
points(log2val[sig.up], minusLogPvals[sig.up],pch=19,cex=0.4,col="brown") # highlight UP
points(log2val[sig.down], minusLogPvals[sig.down],pch=19,cex=0.4,col="blue") # highlight DOWN

###
### MA plot
###
aval <- log((d$ck.mean + d$trt.mean)/2)
mval <- d$log2FC
sig.set <- which(d[,"qval"]<0.05)
plot(aval, mval, xlab="log(mean of expression)", ylab="log2FC",
     ylim=c(-4, 6), pch=19, cex=0.3, col="blue", main="MA plot")
points(aval[sig.set], mval[sig.set], pch=19,cex=0.3,col="red")
abline(h=0,cex=1,col="dark green")

###
### GO term enrichment test
###
# "Sampling" uses random sampling to approximate the true distribution
# and uses it to calculate the p-values for over (and under) representation of categories.
library("goseq")
geneid <- as.character(d$Gene)  # gene vector
de.vector <- as.integer(d$qval<0.05) # a vector to indicate if the gene is DE (0 or 1)
names(de.vector) <- geneid  # the corresponding gene names of the vector
countbias <- rowSums(d[, 2:7])  # total raw reads per gene
godb <- read.delim("./data/maize.gene.go.txt")  # GO database
pwf.counts <- nullp(DEgenes=de.vector, bias.data=countbias) # bias fitting
go <- goseq(pwf=pwf.counts, gene2cat=godb, method="Sampling",
            repcnt = 1000, use_genes_without_cat=F)  # GO enrichment test
head(go)
example.go <- GOTERM[["GO:0004175"]]  # extract detailed information of a GO
Term(example.go)  # return GO term
Definition(example.go)  # return GO definition
