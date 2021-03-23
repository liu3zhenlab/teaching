#source("../scripts/phylogenetic.tree.R")
source("../scripts/phylogenetic.bootstrap.tree.R")

library("ape")

stdfile <- "../3_GATK/4o_gatk.filt.txt"

### build a tree
gd <- read.delim(stdfile, stringsAsFactors=F, comment.char="#",
                 as.is=T, check.names=F)
head(gd)
colnames(gd)

#############################################################################################
### All strains
#############################################################################################
selected.cols <- 5:ncol(gd)
pdf("1o-genome.phylo.tree.pdf", width = 7, height=7)
ape.tree(genodata=gd, colranges=selected.cols, cex=0.75, bootstrap=T,
         newickoutfile="1o_tree.phylo.newick")
dev.off()

#############################################################################################


