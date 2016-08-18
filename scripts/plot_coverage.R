#!/usr/bin/Rscript --vanilla

# script: plot_coverage.R
# required to add the table to the plot
# install once with 'install.packages("plotrix")'
library("plotrix")

setwd("/data1/SRR833244")

folder <- "BDGP5.78_bedtools_genomeCvg"
file <- paste(folder, "BDGP5.78_bwa-mapping-histo.txt", sep="/")
name <- "BDGP5.78_bwa"

outfolder <- paste(folder, "/", name, "-coverage_plots", sep="")
dir.create(outfolder)

data <- read.delim(file, header=F)
colnames(data) <- c("chr", "depth", "count", "chrlen", "cov%")

chr.list <- as.vector(unique(data$chr))

for(chr in chr.list){
  subset <- data[data$chr==chr,]

  # find X-coordinates with 0-coverage
  zero <- subset[subset$depth==0,]
  zero

  # find X-coordinates for 5%, Q1, Q3, 95%
  sub <- subset[subset$depth>0,]
  
  # print summary
  summary(sub)
  
  vec <- rep(sub$depth, sub$count)
  Q001 <- quantile(vec, 0.01)
  Q005 <- quantile(vec, 0.05)
  Q1 <- summary(vec)[2]
  Q2 <- summary(vec)[3]
  Q3 <- summary(vec)[5]
  Q95 <- quantile(vec, 0.95)

  # create table for plot
  leg.table <- data.frame(depth=c(Q001, Q005, Q1, Q2, Q3, Q95))
  leg.table

  # create pictures with coverage distribution
  picname <- paste(outfolder, "/", name, "-", chr, "_coverage-plot.png", sep="")
  title <- paste("coverage depth (", name, "-", chr, ")", sep="")

  png(filename = picname,
  width = 450, height = 450, units = "px", pointsize = 12,
  bg = "white")

  # limits for all and sample datasets
  xmax=2*Q2
  max.count <- max(sub[sub$depth<xmax,]$count)
  ymax <- max.count
  plot(sub$count ~ sub$depth, 
  xlim=c(0, xmax),
  ylim=c(0, ymax),
  xlab=title,
  ylab="bp-count", 
  type="l",
  col="blue",
  lwd=3,
  cex=0.5)

  # add vertical line at pic
  abline(v=Q001, col="grey", lwd=1)
  abline(v=Q005, col="grey", lwd=2)
  abline(v=Q1, col="grey", lwd=2)
  abline(v=Q2, col="red", lwd=2)
  abline(v=Q3, col="grey", lwd=2)
  abline(v=Q95, col="grey", lwd=2)


  # show most of the options
  addtable2plot((70/100*xmax), 
    (60/100*ymax), 
    leg.table,
    bg="white",
    bty="o", 
    display.rownames=TRUE,
    hlines=TRUE, 
    vlines=TRUE,
    cex=1)

  #dev.off()

}
