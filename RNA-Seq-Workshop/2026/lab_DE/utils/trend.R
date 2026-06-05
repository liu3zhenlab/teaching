datasource <- "http://129.130.89.83/tmp/public/RNASeq/RNASeq2017/data/seqcost2015_4.csv"
ngs.cost <- read.csv(datasource)

mb.cost <- gsub("\\$", "", as.character(ngs.cost$Cost.per.Mb))
mb.cost <- gsub(",", "", mb.cost)
mb.cost <- as.numeric(mb.cost)
mb.cost.log10 <- log10(mb.cost)

#png(filename = "cost.png", width = 3, height = 2.5,
#    unit="in", res=1200, pointsize=4)

par(mar=c(6, 4, 3, 0))
plot(NULL, NULL, xaxt="n", yaxt="n",
     xlim = c(1, nrow(ngs.cost)), ylim = range(mb.cost.log10),
     xlab="", ylab="", bty="n",
     main = "cost per megabase")

axis(side = 1, at = 1:nrow(ngs.cost), labels = ngs.cost$Date, las=2,  cex.axis=0.9)
axis(side = 2, at = c(-2, log10(0.05), -1, 0, 1, 2, 3, log10(5000)),tick = F, 
labels = c("$0.01", "$0.05", "$0.1", "$1", "$10", "$100", "$1000", "$5000"), las=2,  cex.axis=1.2, xpd = T)

abline(h = c(-2, log10(0.05), -1, 0, 1, 2, 3, log10(5000)), col="light grey")

points(1:nrow(ngs.cost), mb.cost.log10, col="blue", pch=19)

#dev.off()
