######################################################################################
### Phylogenetic tree
### Sanzhen Liu
### March 2018
######################################################################################
ape.tree <- function(genodata, colranges, newickoutfile=NULL,
                     name.replace=NULL, bootstrap=F, ...) {
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
  
  # Neighbor-Joining Tree Estimation
  tree_method <- function(x) njs(dist.dna(x, model="RAW", pairwise.deletion=T))
  
  phylo <- tree_method(dna)
 
  if (bootstrap) {
	## get 100 bootstrap trees:
  	boot.tree <- boot.phylo(phylo, dna, FUN=tree_method, trees=T)$trees
	## get proportions of each clad:
	boot_values <- prop.clades(phylo, boot.tree)
	phylo$node.label <- boot_values
  }
  
  # plotting
  par(mfrow=c(1,1), mar=c(0.1,0.1,0.1,0.1))
  plot.phylo(phylo, type = "phylogram", use.edge.length=T,
             show.tip.label=T, no.margin=T, show.node.label=T, ...)
  
  ### newick out
  if (!is.null(newickoutfile)) {
  	### newick output
	write.tree(phylo, file=newickoutfile)
  }
  
  ### cleaning:
  unlink(phylotree.data.file)
  invisible(phylo)
}
