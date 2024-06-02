################################################################################
### RNA-Seq Workshop preloading R codes
### Sanzhen Liu
################################################################################
###
### t.test on 10,000 data sets under the null hypothesis
###
ttpvals <- NULL
nullnum <- 10000
for (i in 1:nullnum) {
  tt <- t.test(x=rnorm(n=5,mean=1,sd=0.2), y=rnorm(n=5,mean=1,sd=0.2))
  ttpval <- tt[[3]]
  ttpvals <- c(ttpvals, ttpval)
}

###
### a module to be used in the "pairs" function to plot pair-wise scatter plots
###
panel.cor2 <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  #usr <- par("usr"); on.exit(par(usr))
  x[is.infinite(x)] <- NA
  y[is.infinite(y)] <- NA
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

######################################################################################
# principal component analysis and ploting:
#-------------------------------------------------------------------------------------
rnaseq.pca <- function(norm.data, norm.feature="RPKM", group.feature, title="",
                       shape.code=NULL, mean.cutoff=0.1, colors=NULL, scaling=T, ...) {
### Sanzhen Liu
### Kansas State University
### 5/25/2014
### the function for a PCA plot
  # example:
  # group.feature <- c("mut","wt")
  # shape.code <- c(5,11)
  # colors <- ("dark green","red")
  
  if (is.null(shape.code) | is.null(colors) |
      length(shape.code) != length(group.feature) |
      length(colors) != length(group.feature)) {
    stop("Warnings: elements in shape.code and colors should correspond 
        elements in group.features\n")
  }
  
  norm.data <- norm.data[, grep(norm.feature, colnames(norm.data))]
  new.colnames <- colnames(norm.data)
  data.columns <- NULL
  for (eachsample in group.feature) {
    data.columns <- c(data.columns, grep(eachsample, new.colnames))
  }
  norm.data <- norm.data[, data.columns]
  rmean <- rowMeans(norm.data)
  rrange <- apply(norm.data, 1, max) - apply(norm.data, 1, min)
  
  ### filer data:
  filter.criteria <- (rmean >= mean.cutoff & rrange > 0)
  norm.data <- norm.data[filter.criteria , ]
  rmean <- rmean[filter.criteria]
  rrange <- rrange[filter.criteria]
  
  ### scaling:
  if (scaling) {
    for (i in 1:ncol(norm.data)) {
      norm.data[, i] <- (norm.data[, i] - rmean)/rrange
    }
  }
  
  ### principal component analysis
  pr <- prcomp(norm.data)
  xval <- pr$rotation[,1]; xoffset <- (max(xval)-min(xval))/10
  yval <- pr$rotation[,2]; yoffset <- (max(yval)-min(yval))/10
  range.x <- c(min(xval)-xoffset*5,max(xval)+xoffset)
  range.y <- c(min(yval)-yoffset,max(yval)+yoffset*2)
  
  ## PCA plot:
  plot(NULL,NULL,xlim=range.x,ylim=range.y,xlab="PC1",ylab="PC2", main=title)
  rn <- rownames(pr$rotation)
  for (i in 1:length(group.feature)) {
    points(pr$rotation[grep(group.feature[i],rn),1],
           pr$rotation[grep(group.feature[i],rn),2],
           col=colors[i], pch=shape.code[i], ...)
  }
  # rep information:
  for (i in 1:3) {
    biorep <- paste("rep",i,sep="")
    text(pr$rotation[grep(biorep,rn),1],
         pr$rotation[grep(biorep,rn),2],
         labels=i,
         col="black",cex=1)
  }
  ### legend:
  legend("topleft",group.feature,pch=shape.code,cex=1.2,col=colors,bty="n")
}

###
### load DESeq2
###
library("DESeq2")

###########################################################################################
# normalization methods 
###########################################################################################
normalization <- function (counts, normcolname, libsize, exonsize, methods=c("RPKM","RPM","TMM","Q75")) {
# This R section is for normalization of different datasets
# author: Sanzhen Liu
# 7/29/2011
# "counts" is a data frame or matrix of gene counts
# "normcol" is column name of counting data
# "libsize" is an vector of the libary size with vector names, (bp)
# "exonsize" is an vector of the exon size (bp)
	
	allmethods <- c("RPKM","RPM","TMM","Q75")
	raw.data <- counts[,normcolname]
	print(head(raw.data))
	libsize <- as.numeric(as.character(libsize[normcolname]))
	
	if (length(methods) == 4) {
		methods = "RPM"
		print("The method is not specified. Therefore, the default method is used: RPM.")
	} else {
		if (length(methods) == 1 & sum(methods %in% allmethods)>0) {
			print("The method selected for the normalization is:")
			print(methods)
		} else {
			print("one of methods below must be selected")
			print(allmethods)
			stop("ERROR: please re-setup the Methods parameter")
		}
	}

	label <- paste(".", methods, sep="")
### method 1:
	if (sum(methods == "RPKM")>0) {
		norm.factor = libsize/1000000
		otheradjust = 1000/exonsize
		label <- ".RPKM"
	}
### method 2:
	if (sum(methods == "RPM")>0) {
		norm.factor = libsize/1000000
		otheradjust = 1
		label <- ".RPM"
	}
### method 3:
	if (sum(methods == "TMM")>0) {
		library(edgeR)
		options(warn=-1) # turn off warnings
	# normalization
		f <- calcNormFactors(as.matrix(raw.data))
		f <- f/exp(mean(log(f))) # compute the normalizatioin factor
		f <- round(f,3)
		print(paste(normcolname,"normalization coefficient",f))
	# uniquely mapped reads:
		libsize <- libsize*f
		norm.factor <- libsize/exp(mean(log(libsize)))
		otheradjust = 1
	}

### method 4:
	if (sum(methods == "Q75")>0) {
		q.75 <- apply(raw.data,2,quantile,0.75)
		norm.factor <- q.75/exp(mean(log(q.75)))
		otheradjust = 1
	}
# output
	norm.data <- raw.data
	for (i in 1:ncol(raw.data)) {i
		norm.data[,i] <- raw.data[,i]/norm.factor[i]*otheradjust
		norm.data[,i] <- round(norm.data[,i],3)
	}
	# add colnames:
	colnames(norm.data) <- paste(normcolname, label, sep="")
	# normalized combined data:	
	out <- cbind(counts, norm.data)
	return(out)
}

