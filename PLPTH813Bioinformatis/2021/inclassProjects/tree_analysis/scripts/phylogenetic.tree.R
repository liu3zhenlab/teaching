######################################################################################
### Phylogenetic tree
######################################################################################
ape.tree <- function(genodata, colranges, newickoutfile=NULL, name.replace=NULL, ...) {
  ## convert data format to ape format:
  phylotree.data.file <- ".phylo.tree.input.data.tmp.txt"
  cat(length(colranges), nrow(genodata), "\n", file=phylotree.data.file)
  for (i in colranges) {
    geno.name <- colnames(genodata)[i]
    if (!is.null(name.replace)) {
      geno.name <- gsub(name.replace, "", geno.name)
    }
    cat(geno.name, " ", file=phylotree.data.file, append=T)
    genos <- paste(genodata[, i], collapse="")
    genos.recode <- gsub("0", "N", genos)
    genos.recode <- gsub("1", "A", genos.recode)
    genos.recode <- gsub("2", "G", genos.recode)
    genos.recode <- gsub("3", "R", genos.recode)
    cat(genos.recode, "\n", file=phylotree.data.file, append=T)
  }
  
  ### read "ape" formatted data
  dna <- read.dna(phylotree.data.file, format = "sequential")
  dist.result <- dist.dna(x=dna, model="RAW", pairwise.deletion=T)
  
  # Neighbor-Joining Tree Estimation
  phylo <- njs(dist.result)
  
  ### plotting
  par(mfrow=c(1,1), mar=c(0.1,0.1,0.1,0.1))
  plot.phylo(phylo, type = "phylogram", use.edge.length=T, show.tip.label=T, no.margin=T, ...)
  
  ### newick out
  if (!is.null(newickoutfile)) {
  	### newick output
	write.tree(phylo, file=newickoutfile)
  }

  ### cleaning:
  system("rm phylo.tree.input.data.tmp.txt")
  phylo
}
