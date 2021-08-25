rm(list = ls())
library(GSVA)
library(CancerSubtypes)
library(ggplot2)
library(RTCGA)
library(RTCGA.mRNA)
library(biomaRt)

setwd("~/project/wangfangyu/")
file <- read.table("GBMcells_featureCounts.txt",sep = "\t",header = TRUE,check.names = FALSE)
data <- file[rowSums(file > 3) >2,]
names(file) <- c("geneid","chr","start","end","strand","length","GBM1-1","GBM1-2","GBM1-3","GBM1-4","GBM2-1","GBM2-2","GBM2-3","GBM2-4")
rownames(file) <- file$geneid

df <- file[,c(1,7:14)]
df <- scale(file[,c(7:14)])
#file$geneid <- rownames(file)

# correlation -------------------------------------------------------------
genelist1 <- read.csv("data/class1_Mesenchymal.csv")
genelist2 <- read.csv("data/class2_Proneural.csv")
genelist3 <- read.csv("data/class3_Classical.csv")


humanTomouse <- read.table("~/project/colon cancer/TCGA/human_to_mouse_symbol_ensembl.txt",sep = "\t",header = TRUE)
library(dplyr)
data1 <- inner_join(genelist1,humanTomouse,by=c("GeneSymbol"="Gene.name"))
data2 <- inner_join(genelist2,humanTomouse,by=c("GeneSymbol"="Gene.name"))
data3 <- inner_join(genelist3,humanTomouse,by=c("GeneSymbol"="Gene.name"))

# library(stringr)
# list1 <- str_to_title(genelist1$GeneSymbol)
# list2 <- str_to_title(genelist2$GeneSymbol)
# list3 <- str_to_title(genelist3$GeneSymbol)

df <- file[,7:14]
d <- df[rowSums(df>3) > 2,]
geneset <- list("Mesenchymal" = data1$Gene.name.1,"Proneural"=data2$Gene.name.1,"Classical"=data3$Gene.name.1)
a <- c(data1$Gene.name.1,data2$Gene.name.1,data3$Gene.name.1)
markergene2 <- df[which(rownames(df) %in% data2$Gene.name.1),]
markgene1 <- df[data1$Gene.name.1,]
markgene2 <- df[data2$Gene.name.1,]
markgene3 <- df[data3$Gene.name.1,]

gs <- lapply(geneset, function(x) x[!is.na(x)])
ssgsea_score <- gsva(mark_rna,gs, method="ssgsea",ssgsea.norm = TRUE,verbose = TRUE)
write.table(ssgsea_score,file = "~/project/wangfangyu/ssGSEA_score.txt",sep = "\t",quote = FALSE,row.names = TRUE)

ht <- Heatmap(as.matrix(ssgsea_score),cluster_columns = FALSE)
pdf(file = "~/project/wangfangyu/ssGSEA.pdf",width = 6,height = 3)
draw(ht)
dev.off()

exp <- as.matrix(df[data1$Gene.name.1,])
score <- gsva(exp, list(data1$Gene.name.1),
              method= "ssgsea",
              # rnaseq=TRUE,
              # abs.ranking=FALSE,
              # min.sz=1,
              # max.sz=Inf,
              # no.bootstraps=0,
              # bootstrap.percent = .632,
              # parallel.sz=0,
              # parallel.type="SOCK",
              # mx.diff=TRUE,
              # tau=switch(method, gsva=1, ssgsea=0.25, NA),
              # kernel=TRUE,
              ssgsea.norm=TRUE,
              verbose=TRUE)


library(circlize)
library(ComplexHeatmap)
Heatmap(as.matrix(allgene[,c(2:4)]))
allgene <- rbind(markgene1,markergene2,markgene3)
Heatmap(as.matrix(allgene),cluster_rows = FALSE,col = colorRamp2(c(-2,0,2),c("green","black","red")))
        #row_km = 3)



# method1:correlation -------------------------------------------------------------
# correlation with human 3 cluster gene
# human_cor <- data.frame(Mesenchymal = c(0.06786704,0.07657803,0.01776889,0.03866529, 0.06523356,0.04417379, -0.01300331,-0.01402527),
#                         Proneural = c(0.06516839,0.06531471,0.11662193,0.08412602, 0.03187258,0.009551572,-0.06908419,-0.05798164),
#                         Classical = c(-0.2095700,-0.2014459,-0.10864033,-0.1620750,-0.23086128,-0.1755755,-0.1739921,-0.144041),
#                         row.names = c("GBM1-1","GBM1-2","GBM1-3","GBM1-4","GBM2-1","GBM2-2","GBM2-3","GBM2-4"))
# Heatmap(as.matrix(t(human_cor)),column_names_side = "top",column_names_rot = 90)

library(DESeq2)
file <- read.table("GBMcells_featureCounts.txt",sep = "\t",header = TRUE,check.names = FALSE,row.names = 1)
file <- file[,c(6:13)]
data <- file[rowSums(file > 3) >2,]
names(data) <- c("GBM1-1","GBM1-2","GBM1-3","GBM1-4","GBM2-1","GBM2-2","GBM2-3","GBM2-4")
condition <- factor(c(rep("GBM1",4),rep("GBM2",4)),levels = c("GBM1","GBM2"))
cData <- data.frame(row.names = colnames(data),condition)
dds <- DESeqDataSetFromMatrix(countData = data,colData = cData,design = ~ condition)
dds <- DESeq(dds)
dds

res <- results(dds,contrast = c("condition","GBM2","GBM1"))
res <- res[order(res$pvalue),]
summary(res)
sig <- res[!is.na(res$pvalue) & res$pvalue < 0.05,]
sig.deseq <- rownames(sig)
#gene_list = paste(file_path,time,"/res_control_",time,".txt",sep = "")
#write.table(x = sig,file = gene_list,sep = '\t',row.names = TRUE,quote = FALSE)
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds)),by='row.names',sort=FALSE)

table(res$padj < 0.05)
resdata$change <- as.factor(
  ifelse(
    resdata$pvalue < 0.05 & abs(resdata$log2FoldChange) > 1,
    ifelse(resdata$log2FoldChange > 1,'Up','Down'),
    'Nd'
  )
)

table <- as.data.frame(table(resdata$change))
table

vsd <- getVarianceStabilizedData(dds)
heatmap(cor(vsd),cexCol=0.75,cexRow=0.75)

library(RColorBrewer)
pr <- prcomp(t(vsd))
# first
plot(pr$x,col="white",main="PC plot",xlim=c(-150,100),tck=-0.01,las=2,mgp=c(1.5,0.5,0))
points(pr$x[,1],pr$x[,2],pch=19,type = "p",col=rep(brewer.pal(5,"Dark2"),each=3),
       cex=1.5)
text(pr$x[,1],pr$x[,2],labels=colnames(vsd),
     cex=0.7)

#### 2010_TCGA_core_sample.txt
tcga_core = read.table("TCGA_unified_CORE_ClaNC840.txt", sep = "\t", header = TRUE, 
                     check.names = FALSE, stringsAsFactors = FALSE)
tcga_core <- tcga_core[-1,]
pn <- tcga_core[which(tcga_core$Genes.highly.expressed.in.each.subtype.used.for.Gene.Ontology %in% "PN"),]
cl <- tcga_core[which(tcga_core$Genes.highly.expressed.in.each.subtype.used.for.Gene.Ontology %in% "NL"),]
mes <- tcga_core[which(tcga_core$Genes.highly.expressed.in.each.subtype.used.for.Gene.Ontology %in% "MES"),]

library(dplyr)
data1 <- inner_join(humanTomouse,mes,by=c("Gene.name"="CLID"))
data2 <- inner_join(humanTomouse,pn,by=c("Gene.name"="CLID"))
data3 <- inner_join(humanTomouse,cl,by=c("Gene.name"="CLID"))

frame1 <- df[which(rownames(df) %in% data1$Gene.name.1),]
frame2 <- df[which(rownames(df) %in% data2$Gene.name.1),]
frame3 <- df[which(rownames(df) %in% data3$Gene.name.1),]

frame <- rbind(frame1,frame2,frame3)
frame <- frame[,-5]
Heatmap(as.matrix(frame),cluster_rows = FALSE,col = colorRamp2(c(-2,0,2),c("green","black","red")))

geneset <- list("Mesenchymal" = data1$Gene.name.1,"Proneural"=data2$Gene.name.1,"Classical"=data3$Gene.name.1)
gs <- lapply(geneset, function(x) x[!is.na(x)])
score <- gsva(as.matrix(frame), gs,
              method= "ssgsea",
              # rnaseq=TRUE,
              # abs.ranking=FALSE,
              # min.sz=1,
              # max.sz=Inf,
              # no.bootstraps=0,
              # bootstrap.percent = .632,
              # parallel.sz=0,
              # parallel.type="SOCK",
              # mx.diff=TRUE,
              # tau=switch(method, gsva=1, ssgsea=0.25, NA),
              # kernel=TRUE,
              ssgsea.norm=TRUE,
              verbose=TRUE)

ht <- Heatmap(as.matrix(score),cluster_rows = TRUE,cluster_columns = FALSE)
pdf(file = "~/project/wangfangyu/tcga_core_ssGSEA_remove_2-1.pdf",width = 6,height = 3)
draw(ht)
dev.off()



# GO analysis -------------------------------------------------------------
rm(list = ls())
file <- read.delim(file = "~/project/wangfangyu/cluster18_diff.txt",header = T)
file1 <- read.delim(file = "~/project/wangfangyu/cluster19_diff.txt",header = T)

library(ggplot2)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)

plotbp <- function(data,tf) {
  ggplot(data,aes(x = -log10(pvalue),y = reorder(Description,-log10(pvalue)),fill=-log10(pvalue))) +
    geom_bar(stat = "identity",width = 0.8,position = position_dodge(width = 0.8)) +
    theme_classic(base_size = 14,base_family = "sans") +
    #scale_x_continuous(expand = c(0,0),limits = c(0,40)) + # cluster18
    scale_x_continuous(expand = c(0,0),limits = c(0,120)) + # cluster19
    #scale_x_continuous(expand = c(0,0),limits = c(0,max(ceiling(-log10(bp_0.05$pvalue))))) +
    scale_fill_gradient(low = "black",high = "black") +
    labs(title = paste("Biological process of cluster",tf,sep = ""),x="-log10(pValue)",y=NULL) +
    guides(fill=FALSE) +
    theme(aspect.ratio = 1/0.6,
      plot.title = element_text(size = 14,hjust = 1,vjust = 0.5),
          plot.margin = margin(t=0.1,r=0.5,b=0.1,l=0.1,unit = "cm"),
          axis.text = element_text(colour = 'black',size = 14))
}

bp <- enrichGO(gene = file1$Gene_name,
               keyType = "SYMBOL",
               OrgDb = "org.Mm.eg.db",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05
)
bp_19 <- subset(bp@result,bp@result$pvalue < 0.05,)
plotbp(data = bp_19[c(1:3,5,6,29,38,41),],tf = "19")
ggsave(filename = paste("~/project/wangfangyu/cluster19.pdf",sep = ""),height = 3,width = 6)

bp_0.05 <- bp_0.05[c(2:4,6,7,9,11,12),]
plotbp(data = bp_0.05,tf = "18")
ggsave(filename = paste("~/project/wangfangyu/cluster18.pdf",sep = ""),height = 3,width = 6)

write.table(bp_0.05,file = "~/project/wangfangyu/cluster18_BP.txt",row.names = FALSE,sep = "\t")
write.table(bp_19,file = "~/project/wangfangyu/cluster19_BP.txt",row.names = FALSE,sep = "\t")









