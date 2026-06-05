################## summary.table ####################
# summary.table
#----------------------------------------------
# prepare contrast info:
DE.summary <- function(DE.path=".", DE.files, qval.feature=".qval",
                       log2FC.feature=".log2FC", fdr=0.05,
                       output=T, out.path=".", out.file="project.DE",
					   col.names = T, output.append = F) {
  
  out.file <- paste(out.path, "/", out.file, sep="")
  out <- NULL
  for (each.DE.file in DE.files) {
    DE <- read.delim(paste(DE.path, each.DE.file, sep="/"))
    # significant DEs
  	qval.colnames <- colnames(DE)[grep(qval.feature, colnames(DE))]
	  table.header <- c("comparison", "informative genes",
                      "significant genes", "UPall", "UP1-2fc",
                      "UP2-4fc", "UPgte4fc", "DOWNall", "DOWN1-2fc",
                      "DOWN2-4fc", "DOWNgte4fc");
  	for (i in qval.colnames) {
	  	test.contrast <- gsub(qval.feature, "", i)
	  	print(test.contrast)
	  	log2fc.colname <- paste(test.contrast, log2FC.feature, sep="")
  		informative.set <- DE[!is.na(DE[, i]), ]
		  num.info <- nrow(informative.set)

	  	sig.set <- informative.set[informative.set[, i] <= fdr, ]
  		num.sig <- nrow(sig.set)
      
		  if (log2fc.colname %in% colnames(DE)) {
	  		up.all <- sum(sig.set[, log2fc.colname] > 0)
  			up.0_1 <- sum(sig.set[, log2fc.colname] >= 0 & sig.set[, log2fc.colname] < 1)
    		up.1_2 <- sum(sig.set[, log2fc.colname] >= 1 & sig.set[, log2fc.colname] < 2)
  			up.gte2 <- sum(sig.set[, log2fc.colname] >=2)

			down.all <- sum(sig.set[, log2fc.colname] <= 0)
		  	down.0_1 <- sum(sig.set[, log2fc.colname] <= 0 & sig.set[, log2fc.colname] > (-1))
	  		down.1_2 <- sum(sig.set[, log2fc.colname] <= (-1) & sig.set[, log2fc.colname] > (-2))
			down.gte2 <- sum(sig.set[, log2fc.colname] <= (-2))
		  } else {
	  		up.all <- NA
  			up.0_1 <- NA
	  		up.1_2 <- NA
  			up.gte2 <- NA
			down.all <- NA
			down.0_1 <- NA
			down.1_2 <- NA
		  	down.gte2 <- NA
	  	} # end of if
		
  		single.out <- c(test.contrast, num.info, num.sig, 
		  				up.all, up.0_1, up.1_2, up.gte2,
	  					down.all, down.0_1, down.1_2, down.gte2)
  		names(single.out) <- table.header
		out <- rbind(out, single.out)
	  } # end of second for
  } # end of first for
    # summary.table DE.output:
  # out
  rownames(out) <- 1:nrow(out)
  write.table(out, out.file, row.names=F, quote=F, sep="\t",
              col.names = col.names, append=output.append)
  print(out)
  out
}
