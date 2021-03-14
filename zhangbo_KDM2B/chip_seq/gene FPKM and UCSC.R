rm(list = ls())
setwd("~/project/zhangbo/rna-seq/")

file <- read.delim("all_compare.txt",sep = "\t",header = TRUE)
fpkm <- file[file$gene_name == "Crabp2",c(5,10,14,18,22,26,30)]
fpkm <- file[file$gene_name == "Crabp2",c(8:13,19)]


library(reshape2)
library(ggplot2)
library(ggpubr)
plotFPKM <- function(fpkm,gene,limits,breaks,labely){
  data <- melt(fpkm,id.vars = "gene_name")
  data <- cbind(data,type=rep(c("ctrl","cKO"),each=3))
  
  p <- ggplot(data,aes(x=type,y=value,fill=type)) +
    stat_boxplot(geom = "errorbar",linetype=1,width=0.3,position = "identity") +
    geom_boxplot(outlier.fill = NA,outlier.shape=NA,width=0.6) +
    stat_compare_means(method = "t.test",paired = FALSE,comparisons = list(c("cKO","ctrl")),label = "p.signif",label.y = labely) +
    theme_bw(base_size = 18,base_family = "sans",base_line_size = 1.1) +
    scale_y_continuous(expand = c(0.03,0.03),limits = limits,breaks = breaks) + 
    #scale_x_discrete(expand = c(0.3,0.3),breaks = c("Normal","Tumor"),labels=c("Normal\nn=42","Tumor\nn=599")) +
    #scale_fill_manual(values = c("Normal" = "#D9D9D9","Tumor" = "#A50F15")) +
    labs(title = gene,x=NULL,y="Expression (FPKM)",fill=NULL) +
    theme(aspect.ratio = 1/0.75,
          plot.title = element_text(size = 18,hjust = 0.5,vjust = 0.5,family = "sans"),
          axis.text = element_text(colour = "black",size = 18,family = "sans"),
          axis.text.y = element_text(size = 18,colour = "black",margin = margin(l = 0.1,r = 0.1,unit = "cm")),
          axis.ticks.length = unit(2,'mm'),
          panel.border = element_rect(size = 1.1),
          panel.grid = element_blank()
    ) +
    guides(fill=FALSE)
  ggsave(paste(gene,"fpkm.png"),p,units = "in",height = 3.51,width = 2.62)
  p
}

fpkm <- file[file$gene_name == "Crabp2",c(8:13,19)]
plotFPKM(fpkm = fpkm,gene = "Crabp2",limits = c(0.8,3),breaks = c(1,2,3),labely = 2.5)
  
fpkm <- file[file$gene_name == "Bcor",c(8:13,19)]
plotFPKM(fpkm = fpkm,gene = "Bcor",limits = c(5,25),breaks = c(5,10,15,20,25),labely = 23)

fpkm <- file[file$gene_name == "Col2a1",c(8:13,19)]
plotFPKM(fpkm = fpkm,gene = "Col2a1",limits = c(3,25),breaks = c(5,10,15,20,25),labely = 23)

fpkm <- file[file$gene_name == "Eomes",c(8:13,19)]
plotFPKM(fpkm = fpkm,gene = "Eomes",limits = c(10,40),breaks = c(10,20,30,40),labely = 38)

fpkm <- file[file$gene_name == "Neurog2",c(8:13,19)]
plotFPKM(fpkm = fpkm,gene = "Neurog2",limits = c(8,35),breaks = c(10,20,30),labely = 32)

fpkm <- file[file$gene_name == "Pax6",c(8:13,19)]
plotFPKM(fpkm = fpkm,gene = "Pax6",limits = c(2,6),breaks = c(2,4,6),labely = 5.8)

fpkm <- file[file$gene_name == "Wnt8b",c(8:13,19)]
plotFPKM(fpkm = fpkm,gene = "Wnt8b",limits = c(1,9),breaks = c(3,6,9),labely = 8.5)

fpkm <- file[file$gene_name == "Lef1",c(8:13,19)]
plotFPKM(fpkm = fpkm,gene = "Lef1",limits = c(1,5),breaks = c(1,3,5),labely = 4)

fpkm <- file[file$gene_name == "Camk2a",c(8:13,19)]
plotFPKM(fpkm = fpkm,gene = "Camk2a",limits = c(30,80),breaks = c(30,50,70),labely = 75)

fpkm <- file[file$gene_name == "Camk2b",c(8:13,19)]
plotFPKM(fpkm = fpkm,gene = "Camk2b",limits = c(80,120),breaks = c(80,100,120),labely = 117)

# UCSC plot ---------------------------------------------------------------
rm(list = ls())
library(karyoploteR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BiocFileCache)

ub <- c(ko87_ub <- "~/project/zhangbo/chip/bigwig/87-cKO-H2AK119Ub.bw",
        ko88_ub <- "~/project/zhangbo/chip/bigwig/88-cKO-H2AK119Ub.bw",
        wt89_ub <- "~/project/zhangbo/chip/bigwig/89-Kff-H2AK119Ub.bw",
        wt90_ub <- "~/project/zhangbo/chip/bigwig/90-Kf+H2AK119Ub.bw")

methy <- c(ko87_me3 <- "~/project/zhangbo/chip/bigwig/87-cKO-H3K27me3.bw",
           ko88_me3 <- "~/project/zhangbo/chip/bigwig/88-cKO-H3K27me3.bw",
           wt89_me3 <- "~/project/zhangbo/chip/bigwig/89-Kff-H3K27me3.bw",
           wt90_me3 <- "~/project/zhangbo/chip/bigwig/90-Kf+H3K27me3.bw")

plotUCSC <- function(region) {
  crabp2.region <- toGRanges(region)
  kp <- plotKaryotype(zoom = crabp2.region,genome = "mm10")
  
  genes.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                      karyoplot = kp,
                                      plot.transcripts = TRUE,
                                      plot.transcripts.structure = TRUE)
  
  library(org.Mm.eg.db)
  genes.data <- addGeneNames(genes.data)
  genes.data <- mergeTranscripts(genes.data)
  
  kp <- plotKaryotype(zoom = crabp2.region,genome = "mm10",cex=1)
  kpPlotGenes(kp,data = genes.data,r0 = 0,r1 = 0.1,gene.name.cex = 1)
  
  # plot Crabp2 
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$leftmargin <- 0.35
  pp$topmargin <- 5
  pp$bottommargin <- 10
  pp$ideogramheight <- 0
  pp$data1inmargin <- 5
  
  kp <- plotKaryotype(zoom = crabp2.region,genome = "mm10",cex=1,cytobands = NULL,plot.params = pp)
  kpPlotGenes(kp,data = genes.data,r0 = 0,r1 = 0.05,gene.name.cex = 1)
  kpAddBaseNumbers(kp, tick.dist = 10000,minor.tick.dist = 2000,cex = 1,tick.len = 3,digits = 6,add.units = TRUE)
  
  total.tracks <- length(ub) + length(methy)
  out.at <- autotrack(1:length(ub), total.tracks,margin = 0.3,r0 = 0.1)
  kpAddLabels(kp,labels = "H2AK119Ub",r0 = out.at$r0,r1 = out.at$r1,cex=1,srt=0,pos = 1,label.margin = 0.14)
  kpAddLabels(kp,labels = "cKO",r0 = out.at$r0,r1 = out.at$r1,cex=1,srt=0,pos = 1,label.margin = 0.14)
  kpAddLabels(kp,labels = "WT",r0 = out.at$r0,r1 = out.at$r1,cex=1,srt=0,pos = 1,label.margin = 0.14)
  
  for (i in seq_len(length(ub))) {
    at <- autotrack(i,length(ub),r0 = out.at$r0,r1 = out.at$r1,margin = 0.1)
    if (i < 3) 
      kp <- kpPlotBigWig(kp,data = ub[i],ymax = 500,r0 = at$r0,r1 = at$r1,col = c("#92C5DE"),border = NA) 
    else
      kp <- kpPlotBigWig(kp,data = ub[i],ymax = 500,r0 = at$r0,r1 = at$r1,col = c("#4393C3"),border = NA)
    #kp <- kpPlotBigWig(kp,data = ub[i],ymax = 500,r0 = at$r0,r1 = at$r1,col = c(color))
    computed.max <- ceiling(kp$latest.plot$computed.values$ymax)
    kpAxis(kp,ymin = 0,ymax = 500,numticks = 2,r0 = at$r0,r1 = at$r1)
    kpAddLabels(kp,labels = "rep1",r0 = at$r0,r1 = at$r1,cex=1,label.margin = 0.035)
  }
  
  out.at <- autotrack((length(ub)+1):total.tracks,total.tracks,margin = 0.1,r0 = 0.05)
  kpAddLabels(kp,labels = "H3K27me3",r0 = out.at$r0,r1 = out.at$r1,cex=1,srt=0,pos = 1,label.margin = 0.14)
  kpAddLabels(kp,labels = "cKO",r0 = out.at$r0,r1 = out.at$r1,cex=1,srt=0,pos = 1,label.margin = 0.14)
  kpAddLabels(kp,labels = "WT",r0 = out.at$r0,r1 = out.at$r1,cex=1,srt=0,pos = 1,label.margin = 0.14)
  
  for (i in seq_len(length(methy))) {
    at <- autotrack(i,length(methy),r0 = out.at$r0,r1 = out.at$r1,margin = 0.1)
    if (i < 3) 
      kp <- kpPlotBigWig(kp,data = methy[i],ymax = 1600,r0 = at$r0,r1 = at$r1,col = c("#F4A582"),border = NA)
    else
      kp <- kpPlotBigWig(kp,data = methy[i],ymax = 1600,r0 = at$r0,r1 = at$r1,col = c("#D6604D"),border = NA)
    computed.max <- ceiling(kp$latest.plot$computed.values$ymax)
    kpAxis(kp,ymin = 0,ymax = 1600,numticks = 2,r0 = at$r0,r1 = at$r1)
    kpAddLabels(kp,labels = "rep2",r0 = at$r0,r1 = at$r1,cex=1,label.margin = 0.035,data.panel = 1)
  }
}

library(RColorBrewer)
brewer.pal(name = "RdBu",n = 9)
plotUCSC(region = "chr3:87946666-87955376") # Crabp2
plotUCSC(region = "chr9:118475212-118489132") # Eomes
plotUCSC(region = "chr3:127631135-127637631") # Neurog2

# plot one by one
kpPlotBigWig(kp,data = ko87_ub,ymax = 400,r0 = 0.2,r1 = 0.3)
kpAddLabels(kp,labels = "rep1",r0 = 0.2,r1 = 0.3,cex=1,label.margin = 0.035)
kpAxis(kp,ymin = 0,ymax = 400,r0 = 0.2,r1 = 0.3)

kpPlotBigWig(kp,data = ko88_ub,ymax = 400,r0 = 0.35,r1 = 0.45)
kpAddLabels(kp,labels = "rep2",r0 = 0.35,r1 = 0.45,cex=1,label.margin = 0.035)
kpAxis(kp,ymin = 0,ymax = 400,r0 = 0.35,r1 = 0.45)

for (i in seq_len(length(ub))) {
  at <- autotrack(i,length(ub),r0 = 0.2,r1 = 1)
  kpPlotBigWig(kp,data = ub[i],ymax = 400,r0 = at$r0,r1 = at$r1)
  kpAxis(kp,ymin = 0,ymax = 400,r0 = at$r0,r1 = at$r1)
  kpAddLabels(kp,labels = basename(ub)[i],r0 = at$r0,r1 = at$r1,cex=1,label.margin = 0.035)
}

plotDefaultPlotParams(plot.type = 1)

