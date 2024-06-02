setwd("/homes/liu3zhen/public_html/PLPTH813/labs/lab10_DE/4_DE")
library("DESeq2")

### Parameters - Subject to change
datapath <- "../3_aln"
suffix <- "ReadsPerGene.out.tab"
count.files <- dir(path = datapath, pattern = suffix)

### merge all counts
allcounts <- NULL
for (cf in count.files) {
	counts <- read.delim(paste0(datapath, "/", cf), header = F, stringsAsFactors = F, skip = 4)
	base <- gsub(suffix, "", cf)
	counts <- counts[, 1:2]
	colnames(counts) <- c("Gene", base)

	### merge data
	if (is.null(allcounts)) {
		allcounts <- counts
	} else {
		allcounts <- merge(allcounts, counts, by = "Gene")
	}
}

#################################################################
# count information
geneid <- allcounts$Gene
in.data <- allcounts[, 2:7]
in.data <- as.matrix(in.data)
rownames(in.data) <- geneid

# sample names and grouping information (treatment)
sample.ids <- colnames(in.data)
treatment <- c("cold", "cold", "cold", "norm", "norm", "norm")
sample.info <- data.frame(row.names=sample.ids, trt=treatment)
dds <- DESeqDataSetFromMatrix(countData=in.data,
       colData=sample.info, formula(~trt))
dds <- DESeq(dds, "Wald")
#################################################################

### DE
# load modules
source("../utils/scripts/DESeq2.single.trt.R")
source("../utils/scripts/DE.summary.R")
# DE parameters
fdr.cutoff <- 0.05
# data reformat:
input <- allcounts[, 2:7]
rownames(input) <- allcounts[, 1]
# DE statistical analysis:
DE.out <- DESeq2.single.trt(input.matrix = input,
		min.mean.reads = 5,
		group1.col = 1:3,
		group2.col = 4:6,
		comparison = c("norm", "cold"),
		geneID = rownames(input),
		fdr = fdr.cutoff,
		logpath = ".",
		logfile = "cold-norm.log.md")

# merge DE with counts and output DE result:
DE.out <- data.frame(DE.out)
final.out <- merge(allcounts, DE.out, by.x = "Gene", by.y = "GeneID")
head(final.out)
write.table(final.out, "cold-norm.DESeq2.txt", sep="\t", quote=F, row.names=F)

# check p-values
hist(DE.out$cold_norm.pval, xlab = "p-values",
    ylab = "Number of genes", main = "cold vs. norm")

# summary of the DE result
de.summary <- DE.summary(DE.path=".",
  DE.files="cold-norm.DESeq2.txt",
  qval.feature=".qval",
  log2FC.feature=".log2FC",
  fdr=fdr.cutoff,
  out.path=".",
  out.file="cold-norm.DESeq2.summary.txt")

